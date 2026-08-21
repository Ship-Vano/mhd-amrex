//
// MhdAmr.cpp — реализация AMR-решателя 2D идеальной МГД (см. MhdAmr.H).
//
// Код ориентирован на AMReX >= 24.09 (дизайн сверен с официальным туториалом
// Amr/Advection_AmrCore и документацией
// https://amrex-codes.github.io/amrex/docs_html/AmrCore.html).
//
#include "MhdAmr.H"

#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_Interpolater.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_TagBox.H>
#include <AMReX_Loop.H>
#include <AMReX_Utility.H>
#ifdef AMREX_USE_HDF5
#include <AMReX_PlotFileUtilHDF5.H>
#endif

#include <filesystem>
#include <limits>
#include <algorithm>

#include "kernels/Hlld.H"
#include "kernels/Reconstruction.H"
#include "kernels/CtUpdate.H"

using namespace amrex;

namespace mhd {

namespace {

constexpr int NGROW = 3;   // фантомных слоёв у клеточных данных: хватает для
                           // MUSCL-реконструкции на гранях, выходящих на 1
                           // ячейку за валидную область (нужно узловым ЭДС)

// Текущая задача — для функции заполнения ext_dir-границ (CpuBndryFuncFab).
const Problem* g_problem = nullptr;
double g_gamma = 5.0 / 3.0;

// Заполнение физических границ клеточных величин типа ext_dir («исторические»
// ГУ: замороженные значения начального условия, как в тесте Брио–Ву из ВКРБ).
// Сигнатура — amrex::UserFillBox (см. AMReX_PhysBCFunct.H, AMReX >= 22.xx):
// CpuBndryFuncFab сам обрабатывает foextrap/reflect_* и вызывает эту функцию
// только при наличии ext_dir-компонент; bcr указывает на BCRec первой
// заполняемой компоненты (bcomp = 0).
void ext_dir_fill(amrex::Box const& bx, amrex::Array4<amrex::Real> const& arr,
                  int dcomp, int numcomp,
                  amrex::GeometryData const& geom, amrex::Real /*time*/,
                  const amrex::BCRec* bcr, int /*bcomp*/, int /*orig_comp*/)
{
    const Box& domain = geom.Domain();
    const double dx0 = geom.CellSize(0), dx1 = geom.CellSize(1);
    const double plo0 = geom.ProbLo(0),  plo1 = geom.ProbLo(1);

    amrex::LoopOnCpu(bx, [&] (int i, int j, int k)
    {
        if (domain.contains(IntVect(AMREX_D_DECL(i, j, k)))) return;
        bool ext = false;
        for (int dim = 0; dim < AMREX_SPACEDIM; ++dim) {
            const int lo = domain.smallEnd(dim), hi = domain.bigEnd(dim);
            const int idx = (dim == 0) ? i : j;
            if ((idx < lo && bcr[0].lo(dim) == BCType::ext_dir) ||
                (idx > hi && bcr[0].hi(dim) == BCType::ext_dir)) ext = true;
        }
        if (!ext) return;
        const double x = plo0 + (i + 0.5) * dx0;
        const double y = plo1 + (j + 0.5) * dx1;
        double q[NPRIM] = {0}, uc[NCONS];
        g_problem->prim(x, y, q);
        if (g_problem->direct_b) {
            q[QBX] = 0.5 * (g_problem->Bx_face(x - 0.5*dx0, y) + g_problem->Bx_face(x + 0.5*dx0, y));
            q[QBY] = 0.5 * (g_problem->By_face(x, y - 0.5*dx1) + g_problem->By_face(x, y + 0.5*dx1));
        } else {
            q[QBX] =  (g_problem->Az(x, y + 0.5*dx1) - g_problem->Az(x, y - 0.5*dx1)) / dx1;
            q[QBY] = -(g_problem->Az(x + 0.5*dx0, y) - g_problem->Az(x - 0.5*dx0, y)) / dx0;
        }
        prim_to_cons(q, uc, g_gamma);
        // Заполняем запрошенный диапазон компонент; в этом коде FillPatch
        // всегда вызывается на всех NCONS компонентах сразу (dcomp = 0).
        for (int n = 0; n < numcomp && (dcomp + n) < NCONS; ++n)
            arr(i, j, k, dcomp + n) = uc[dcomp + n];
    });
}

int bc_code_for(BcType t, int comp, int dim)
{
    switch (t) {
    case BcType::Periodic:  return BCType::int_dir;
    case BcType::Outflow:   return BCType::foextrap;
    case BcType::Dirichlet: return BCType::ext_dir;
    case BcType::Reflect: {
        // нечётное отражение для нормальной скорости и нормальной компоненты B
        const bool odd = (dim == 0) ? (comp == UMX || comp == UBX)
                                    : (comp == UMY || comp == UBY);
        return odd ? BCType::reflect_odd : BCType::reflect_even;
    }
    }
    return BCType::foextrap;
}

} // namespace

// ===========================================================================
// Конструирование геометрии уровня 0 и AmrInfo до вызова базового
// конструктора. ВАЖНО: AmrMesh::checkInput() выполняется внутри конструктора
// базового класса, поэтому blocking_factor/max_grid_size нельзя поправить
// сеттерами в теле нашего конструктора — они задаются заранее через AmrInfo.
// ===========================================================================
namespace {

Geometry make_level0_geometry(const SimConfig& cfg)
{
    const Box dom(IntVect(AMREX_D_DECL(0, 0, 0)),
                  IntVect(AMREX_D_DECL(cfg.n_cell[0] - 1, cfg.n_cell[1] - 1, 0)));
    const RealBox rb({AMREX_D_DECL(cfg.prob_lo[0], cfg.prob_lo[1], 0.0)},
                     {AMREX_D_DECL(cfg.prob_hi[0], cfg.prob_hi[1], 1.0)});
    const Array<int, AMREX_SPACEDIM> is_per
        {AMREX_D_DECL(cfg.bc_xlo == BcType::Periodic ? 1 : 0,
                      cfg.bc_ylo == BcType::Periodic ? 1 : 0, 0)};
    return Geometry(dom, rb, 0 /*декартовы координаты*/, is_per);
}

AmrInfo make_amr_info(const SimConfig& cfg)
{
    AmrInfo info;
    info.max_level = cfg.max_level;
    const int nlev = cfg.max_level + 1;
    info.ref_ratio.assign(std::max(cfg.max_level, 1), IntVect(cfg.ref_ratio));
    info.n_error_buf.assign(nlev, IntVect(cfg.n_error_buf));

    // AMReX требует: (а) размер домена делится на blocking_factor покоординатно,
    // (б) blocking_factor — степень двойки, (в) max_grid_size кратен ему.
    // Подбираем для каждого направления наибольшую степень двойки, делящую
    // n_cell[d] и не превышающую запрошенный blocking_factor (узкие домены
    // вроде 512×4 у Брио–Ву или 64×110 у альфвеновской волны обрабатываются
    // автоматически).
    IntVect bf, mgs;
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        int p = 1;
        while (2 * p <= cfg.blocking_factor && cfg.n_cell[d] % (2 * p) == 0) p *= 2;
        bf[d]  = p;
        mgs[d] = std::max(p, (cfg.max_grid_size / p) * p);   // кратен bf
        if (p != cfg.blocking_factor) {
            amrex::Print() << "[mhd2d] blocking_factor по направлению " << d
                           << " уменьшен до " << p << " (n_cell=" << cfg.n_cell[d]
                           << " не делится на " << cfg.blocking_factor << ")\n";
        }
    }
    info.blocking_factor.assign(nlev, bf);
    info.max_grid_size.assign(nlev, mgs);
    return info;
}

} // namespace

MhdAmr::MhdAmr(const SimConfig& cfg)
    : AmrCore(make_level0_geometry(cfg), make_amr_info(cfg)),
      cfg_(cfg), prob_(make_problem(cfg))
{
    g_problem = &prob_;
    g_gamma   = cfg.gamma;

    const int nlev = max_level + 1;
    state_.resize(nlev);  state0_.resize(nlev);
    bface_.resize(nlev);  bface0_.resize(nlev);
    flux_.resize(nlev);   emf_.resize(nlev);

    // BCRec клеточных величин для FillPatch
    bcrec_.resize(NCONS);
    const BcType lo[2] = { cfg.bc_xlo, cfg.bc_ylo };
    const BcType hi[2] = { cfg.bc_xhi, cfg.bc_yhi };
    for (int n = 0; n < NCONS; ++n)
        for (int dim = 0; dim < AMREX_SPACEDIM; ++dim) {
            bcrec_[n].setLo(dim, bc_code_for(lo[dim], n, dim));
            bcrec_[n].setHi(dim, bc_code_for(hi[dim], n, dim));
        }
}

// ---------------------------------------------------------------------------
void MhdAmr::AllocLevel(int lev, const BoxArray& ba, const DistributionMapping& dm)
{
    state_[lev].define(ba, dm, NCONS, NGROW);
    state0_[lev].define(ba, dm, NCONS, NGROW);
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        const BoxArray fba = amrex::convert(ba, IntVect::TheDimensionVector(d));
        bface_[lev][d].define(fba, dm, 1, NGROW);
        bface0_[lev][d].define(fba, dm, 1, NGROW);
        flux_[lev][d].define(fba, dm, NCONS, 1);
    }
    emf_[lev].define(amrex::convert(ba, IntVect::TheNodeVector()), dm, 1, 0);
}

void MhdAmr::ClearLevel(int lev)
{
    state_[lev].clear(); state0_[lev].clear(); emf_[lev].clear();
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        bface_[lev][d].clear(); bface0_[lev][d].clear(); flux_[lev][d].clear();
    }
}

// НУ уровня: грани — через векторный потенциал Az (div B = 0 машинно),
// клеточные B — RT0-интерполяция граней (как в схеме Авдеевой–Лукина).
void MhdAmr::InitLevelData(int lev)
{
    const auto problo = Geom(lev).ProbLoArray();
    const auto dx     = Geom(lev).CellSizeArray();
    const Problem& P  = prob_;
    const double gam  = cfg_.gamma;

    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        for (MFIter mfi(bface_[lev][d]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();      // вместе с фантомами
            auto b = bface_[lev][d].array(mfi);
            amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
                const double xf = problo[0] + i * dx[0] + (d == 0 ? 0.0 : 0.5*dx[0]);
                const double yf = problo[1] + j * dx[1] + (d == 1 ? 0.0 : 0.5*dx[1]);
                if (P.direct_b) {
                    b(i,j,k) = (d == 0) ? P.Bx_face(xf, yf) : P.By_face(xf, yf);
                } else if (d == 0) {
                    b(i,j,k) =  (P.Az(xf, yf + 0.5*dx[1]) - P.Az(xf, yf - 0.5*dx[1])) / dx[1];
                } else {
                    b(i,j,k) = -(P.Az(xf + 0.5*dx[0], yf) - P.Az(xf - 0.5*dx[0], yf)) / dx[0];
                }
            });
        }
    }
    for (MFIter mfi(state_[lev]); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.fabbox();
        auto u   = state_[lev].array(mfi);
        auto bxf = bface_[lev][0].const_array(mfi);
        auto byf = bface_[lev][1].const_array(mfi);
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
            const double x = problo[0] + (i + 0.5) * dx[0];
            const double y = problo[1] + (j + 0.5) * dx[1];
            double q[NPRIM] = {0}, uc[NCONS];
            P.prim(x, y, q);
            q[QBX] = 0.5 * (bxf(i,j,k) + bxf(i+1,j,k));
            q[QBY] = 0.5 * (byf(i,j,k) + byf(i,j+1,k));
            prim_to_cons(q, uc, gam);
            for (int n = 0; n < NCONS; ++n) u(i,j,k,n) = uc[n];
        });
    }
}

void MhdAmr::MakeNewLevelFromScratch(int lev, Real, const BoxArray& ba,
                                     const DistributionMapping& dm)
{
    AllocLevel(lev, ba, dm);
    InitLevelData(lev);
}

void MhdAmr::MakeNewLevelFromCoarse(int lev, Real time, const BoxArray& ba,
                                    const DistributionMapping& dm)
{
    AllocLevel(lev, ba, dm);
    // Клеточные величины — консервативная интерполяция с грубого уровня
    {
        PhysBCFunct<CpuBndryFuncFab> cbc(Geom(lev-1), bcrec_, CpuBndryFuncFab(ext_dir_fill));
        PhysBCFunct<CpuBndryFuncFab> fbc(Geom(lev),   bcrec_, CpuBndryFuncFab(ext_dir_fill));
        amrex::InterpFromCoarseLevel(state_[lev], time, state_[lev-1], 0, 0, NCONS,
                                     Geom(lev-1), Geom(lev), cbc, 0, fbc, 0,
                                     refRatio(lev-1), &cell_cons_interp, bcrec_, 0);
    }
    // Граневые B — бездивергентная интерполяция AMReX (face_divfree_interp):
    // дискретная дивергенция мелких ячеек равна (нулевой) дивергенции грубой
    {
        Array<MultiFab*, AMREX_SPACEDIM> fmf {AMREX_D_DECL(&bface_[lev][0],   &bface_[lev][1],   nullptr)};
        Array<MultiFab*, AMREX_SPACEDIM> cmf {AMREX_D_DECL(&bface_[lev-1][0], &bface_[lev-1][1], nullptr)};
        Array<PhysBCFunctNoOp, AMREX_SPACEDIM> nbc;
        Array<Vector<BCRec>, AMREX_SPACEDIM> fbcr;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) fbcr[d].resize(1);
        amrex::InterpFromCoarseLevel(fmf, IntVect(NGROW), time, cmf, 0, 0, 1,
                                     Geom(lev-1), Geom(lev), nbc, 0, nbc, 0,
                                     refRatio(lev-1), &face_divfree_interp, fbcr, 0);
    }
    FillPhysicalFaceBoundary(lev);
    SyncCellB(lev);
}

void MhdAmr::RemakeLevel(int lev, Real time, const BoxArray& ba,
                         const DistributionMapping& dm)
{
    // Новые контейнеры; данные — FillPatch'ем из старых (свой + грубый уровень)
    MultiFab new_state(ba, dm, NCONS, NGROW);
    FillPatchCells(lev, new_state, time);

    Array<MultiFab, AMREX_SPACEDIM> new_b;
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
        new_b[d].define(amrex::convert(ba, IntVect::TheDimensionVector(d)), dm, 1, NGROW);
    Array<MultiFab*, AMREX_SPACEDIM> nbp {AMREX_D_DECL(&new_b[0], &new_b[1], nullptr)};
    FillPatchFaces(lev, nbp, time);

    ClearLevel(lev);
    AllocLevel(lev, ba, dm);
    MultiFab::Copy(state_[lev], new_state, 0, 0, NCONS, NGROW);
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
        MultiFab::Copy(bface_[lev][d], new_b[d], 0, 0, 1, NGROW);
    FillPhysicalFaceBoundary(lev);
    SyncCellB(lev);
}

// ---------------------------------------------------------------------------
// Критерий измельчения: относительный градиент плотности и/или ток jz = (∇×B)z
// ---------------------------------------------------------------------------
void MhdAmr::ErrorEst(int lev, TagBoxArray& tags, Real /*time*/, int /*ngrow*/)
{
    const Real grho = cfg_.refine_grad_rho;
    const Real gcur = cfg_.refine_current;
    const auto dx   = Geom(lev).CellSizeArray();

    // Для разностей через соседей нужны актуальные фантомы
    MultiFab tmp(grids[lev], dmap[lev], NCONS, 1);
    FillPatchCells(lev, tmp, t_);

    for (MFIter mfi(tmp, TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.tilebox();
        auto u   = tmp.const_array(mfi);
        auto tag = tags.array(mfi);
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
            const Real r  = u(i,j,k,URHO);
            const Real gx = std::abs(u(i+1,j,k,URHO) - u(i-1,j,k,URHO));
            const Real gy = std::abs(u(i,j+1,k,URHO) - u(i,j-1,k,URHO));
            bool t = (Real(0.5) * (gx + gy) / r > grho);
            if (gcur > Real(0.0)) {
                const Real jz = (u(i+1,j,k,UBY) - u(i-1,j,k,UBY)) / (2*dx[0])
                              - (u(i,j+1,k,UBX) - u(i,j-1,k,UBX)) / (2*dx[1]);
                t = t || (std::abs(jz) * std::min(dx[0], dx[1]) > gcur);
            }
            if (t) tag(i,j,k) = TagBox::SET;
        });
    }
}

// ---------------------------------------------------------------------------
// FillPatch: клеточные величины (фантомы внутри уровня, с грубого уровня и
// физические ГУ через BCRec + ext_dir-наполнитель)
// ---------------------------------------------------------------------------
void MhdAmr::FillPatchCells(int lev, MultiFab& mf, Real time)
{
    PhysBCFunct<CpuBndryFuncFab> fbc(Geom(lev), bcrec_, CpuBndryFuncFab(ext_dir_fill));
    if (lev == 0) {
        amrex::FillPatchSingleLevel(mf, time, {&state_[0]}, {time}, 0, 0, NCONS,
                                    Geom(0), fbc, 0);
    } else {
        PhysBCFunct<CpuBndryFuncFab> cbc(Geom(lev-1), bcrec_, CpuBndryFuncFab(ext_dir_fill));
        amrex::FillPatchTwoLevels(mf, time,
                                  {&state_[lev-1]}, {time}, {&state_[lev]}, {time},
                                  0, 0, NCONS, Geom(lev-1), Geom(lev),
                                  cbc, 0, fbc, 0,
                                  refRatio(lev-1), &cell_cons_interp, bcrec_, 0);
    }
}

// FillPatch: граневые компоненты B (face_divfree_interp на стыке уровней)
void MhdAmr::FillPatchFaces(int lev, Array<MultiFab*, AMREX_SPACEDIM> bf, Real time)
{
    if (lev == 0) {
        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            if (bf[d] != &bface_[0][d]) {
                MultiFab::Copy(*bf[d], bface_[0][d], 0, 0, 1, 0);
            }
            bf[d]->FillBoundary(Geom(0).periodicity());
        }
    } else {
        Array<PhysBCFunctNoOp, AMREX_SPACEDIM> nbc;
        Array<Vector<BCRec>, AMREX_SPACEDIM> fbcr;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) fbcr[d].resize(1);
        // ВНИМАНИЕ: в этой перегрузке компонентные смещения — Array<int,DIM>
        const Array<int, AMREX_SPACEDIM> zerocomp {AMREX_D_DECL(0, 0, 0)};
        Vector<Array<MultiFab*, AMREX_SPACEDIM>> cmf
            { {AMREX_D_DECL(&bface_[lev-1][0], &bface_[lev-1][1], nullptr)} };
        Vector<Array<MultiFab*, AMREX_SPACEDIM>> fmf
            { {AMREX_D_DECL(&bface_[lev][0], &bface_[lev][1], nullptr)} };
        amrex::FillPatchTwoLevels(bf, IntVect(NGROW), time, cmf, {time}, fmf, {time},
                                  0, 0, 1, Geom(lev-1), Geom(lev),
                                  nbc, zerocomp, nbc, zerocomp,
                                  refRatio(lev-1), &face_divfree_interp, fbcr, zerocomp);
    }
}

// Физические ГУ для граневых B (фантомные грани). Периодика — FillBoundary.
void MhdAmr::FillPhysicalFaceBoundary(int lev)
{
    const Box& domain = Geom(lev).Domain();
    const auto problo = Geom(lev).ProbLoArray();
    const auto dx     = Geom(lev).CellSizeArray();
    const BcType lo[2] = { cfg_.bc_xlo, cfg_.bc_ylo };
    const BcType hi[2] = { cfg_.bc_xhi, cfg_.bc_yhi };
    const Problem& P = prob_;

    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        bface_[lev][d].FillBoundary(Geom(lev).periodicity());
        const Box fdomain = amrex::convert(domain, IntVect::TheDimensionVector(d));
        for (MFIter mfi(bface_[lev][d]); mfi.isValid(); ++mfi) {
            const Box& fb = mfi.fabbox();
            auto b = bface_[lev][d].array(mfi);
            amrex::LoopOnCpu(fb, [&] (int i, int j, int k) {
                IntVect iv(AMREX_D_DECL(i, j, k));
                for (int dim = 0; dim < AMREX_SPACEDIM; ++dim) {
                    const int dlo = fdomain.smallEnd(dim), dhi = fdomain.bigEnd(dim);
                    BcType bt;
                    if      (iv[dim] < dlo) bt = lo[dim];
                    else if (iv[dim] > dhi) bt = hi[dim];
                    else continue;
                    if (bt == BcType::Periodic) continue;
                    if (bt == BcType::Outflow) {
                        IntVect src = iv;
                        src[dim] = std::clamp(iv[dim], dlo, dhi);
                        b(i,j,k) = b(src[0], src[1], 0);
                    } else if (bt == BcType::Reflect) {
                        // нормальная к границе компонента — нечётная,
                        // касательная — чётная (зеркальное отражение поля)
                        const bool normal = (dim == d);
                        IntVect src = iv;
                        if (normal) src[dim] = (iv[dim] < dlo) ? 2*dlo - iv[dim] : 2*dhi - iv[dim];
                        else        src[dim] = (iv[dim] < dlo) ? 2*dlo - iv[dim] - 1 : 2*dhi - iv[dim] + 1;
                        b(i,j,k) = (normal ? -1.0 : 1.0) * b(src[0], src[1], 0);
                    } else {  // Dirichlet: «исторические» значения из НУ
                        const double xf = problo[0] + i * dx[0] + (d == 0 ? 0.0 : 0.5*dx[0]);
                        const double yf = problo[1] + j * dx[1] + (d == 1 ? 0.0 : 0.5*dx[1]);
                        if (P.direct_b) {
                            b(i,j,k) = (d == 0) ? P.Bx_face(xf, yf) : P.By_face(xf, yf);
                        } else if (d == 0) {
                            b(i,j,k) =  (P.Az(xf, yf+0.5*dx[1]) - P.Az(xf, yf-0.5*dx[1])) / dx[1];
                        } else {
                            b(i,j,k) = -(P.Az(xf+0.5*dx[0], yf) - P.Az(xf-0.5*dx[0], yf)) / dx[0];
                        }
                    }
                    break;
                }
            });
        }
    }
}

// Интерполяция граневых B в центры ячеек (аналог базиса Равьяра–Тома):
// на прямоугольной ячейке RT0-восстановление даёт среднее двух граней.
void MhdAmr::SyncCellB(int lev)
{
    for (MFIter mfi(state_[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box bx = mfi.growntilebox(NGROW - 1);
        auto u   = state_[lev].array(mfi);
        auto bxf = bface_[lev][0].const_array(mfi);
        auto byf = bface_[lev][1].const_array(mfi);
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
            u(i,j,k,UBX) = 0.5 * (bxf(i,j,k) + bxf(i+1,j,k));
            u(i,j,k,UBY) = 0.5 * (byf(i,j,k) + byf(i,j+1,k));
        });
    }
}

// ---------------------------------------------------------------------------
// HLLD-потоки на гранях и узловые ЭДС одного уровня
// ---------------------------------------------------------------------------
void MhdAmr::ComputeFluxesAndEmf(int lev)
{
    const Limiter lim = cfg_.limiter;
    const EmfAveraging emode = cfg_.emf;
    const double gam = cfg_.gamma;

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    for (MFIter mfi(state_[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto u   = state_[lev].const_array(mfi);
        auto bxf = bface_[lev][0].const_array(mfi);
        auto byf = bface_[lev][1].const_array(mfi);
        auto fx  = flux_[lev][0].array(mfi);
        auto fy  = flux_[lev][1].array(mfi);
        auto ez  = emf_[lev].array(mfi);

        // примитивы в ячейке (i,j) по запросу
        auto qprim = [=] (int i, int j, double* q) {
            double uc[NCONS];
            for (int n = 0; n < NCONS; ++n) uc[n] = u(i, j, 0, n);
            cons_to_prim(uc, q, gam);
        };

        // --- x-потоки: грани валидной области + 1 слой (нужно узловым ЭДС) --
        {
            const Box xbx = mfi.grownnodaltilebox(0, 1);
            amrex::LoopOnCpu(xbx, [&] (int i, int j, int k) {
                double qm[NPRIM], q0[NPRIM], qp[NPRIM], qq[NPRIM];
                double qL[NPRIM], qR[NPRIM], f[NCONS];
                qprim(i-2, j, qm); qprim(i-1, j, q0); qprim(i, j, qp); qprim(i+1, j, qq);
                for (int n = 0; n < NPRIM; ++n) {
                    qL[n] = face_value_plus (qm[n], q0[n], qp[n], lim);
                    qR[n] = face_value_minus(q0[n], qp[n], qq[n], lim);
                }
                hlld_flux(qL, qR, bxf(i,j,k), f, gam);   // Bn — из staggered-массива!
                for (int n = 0; n < NCONS; ++n) fx(i,j,k,n) = f[n];
            });
        }
        // --- y-потоки: локальный поворот осей (u'=v, v'=−u, Bx'=By, By'=−Bx) -
        {
            const Box ybx = mfi.grownnodaltilebox(1, 1);
            amrex::LoopOnCpu(ybx, [&] (int i, int j, int k) {
                double qm[NPRIM], q0[NPRIM], qp[NPRIM], qq[NPRIM];
                double qL[NPRIM], qR[NPRIM], rL[NPRIM], rR[NPRIM], f[NCONS];
                qprim(i, j-2, qm); qprim(i, j-1, q0); qprim(i, j, qp); qprim(i, j+1, qq);
                for (int n = 0; n < NPRIM; ++n) {
                    qL[n] = face_value_plus (qm[n], q0[n], qp[n], lim);
                    qR[n] = face_value_minus(q0[n], qp[n], qq[n], lim);
                }
                auto rot = [] (const double* q, double* r) {
                    r[QRHO]=q[QRHO]; r[QP]=q[QP]; r[QW]=q[QW]; r[QBZ]=q[QBZ];
                    r[QU]=q[QV]; r[QV]=-q[QU]; r[QBX]=q[QBY]; r[QBY]=-q[QBX];
                };
                rot(qL, rL); rot(qR, rR);
                hlld_flux(rL, rR, byf(i,j,k), f, gam);
                fy(i,j,k,URHO)=f[URHO]; fy(i,j,k,UENE)=f[UENE];
                fy(i,j,k,UMZ)=f[UMZ];   fy(i,j,k,UBZ)=f[UBZ];
                fy(i,j,k,UMX)=-f[UMY];  fy(i,j,k,UMY)=f[UMX];
                fy(i,j,k,UBX)=-f[UBY];  fy(i,j,k,UBY)=f[UBX];   // fy[UBX] = +Ez
            });
        }
        // --- узловые ЭДС Ez(i−1/2, j−1/2): усреднение ЭДС примыкающих граней —
        // декартов аналог усреднения по рёбрам, сходящимся в вершине (ф.(8)
        // статьи Авдеевой–Лукина); GardinerStone добавляет клеточную поправку.
        // ЭДС на x-грани: Ez = −Fx[UBY]; на y-грани: Ez = +Fy[UBX].
        {
            const Box nbx = mfi.tilebox(IntVect::TheNodeVector());
            amrex::LoopOnCpu(nbx, [&] (int i, int j, int k) {
                double qmm[NPRIM], qpm[NPRIM], qmp[NPRIM], qpp[NPRIM];
                qprim(i-1, j-1, qmm); qprim(i, j-1, qpm);
                qprim(i-1, j,   qmp); qprim(i, j,   qpp);
                ez(i,j,k) = corner_emf(-fx(i, j-1, k, UBY), -fx(i, j, k, UBY),
                                        fy(i-1, j, k, UBX),  fy(i, j, k, UBX),
                                        cell_emf_z(qmm), cell_emf_z(qpm),
                                        cell_emf_z(qmp), cell_emf_z(qpp), emode);
            });
        }
    }
}

// Инжекция узловых ЭДС мелкого уровня в совпадающие узлы грубого. Поскольку
// расчёт ведётся БЕЗ подциклирования (общий Δt), после инжекции обновление
// грубой грани на границе уровней в точности равно среднему обновлений
// накрывающих её мелких граней → average_down_faces не вносит дивергенцию,
// и div B = 0 сохраняется на всей иерархии (вывод — в REPORT.md).
void MhdAmr::SyncEmfAcrossLevels()
{
    for (int lev = finest_level; lev >= 1; --lev) {
        const BoxArray cba = amrex::coarsen(emf_[lev].boxArray(), refRatio(lev-1));
        MultiFab cemf(cba, emf_[lev].DistributionMap(), 1, 0);
        amrex::average_down_nodal(emf_[lev], cemf, refRatio(lev-1));
        emf_[lev-1].ParallelCopy(cemf, 0, 0, 1);
    }
}

// ---------------------------------------------------------------------------
void MhdAmr::ApplyUpdates(int lev, Real dt)
{
    const auto dx = Geom(lev).CellSizeArray();
    const Real lx = dt / dx[0], ly = dt / dx[1];

#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
    {
        // Газовые величины: годуновский баланс потоков (ур. (1)–(3) статьи)
        for (MFIter mfi(state_[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();
            auto u  = state_[lev].array(mfi);
            auto fx = flux_[lev][0].const_array(mfi);
            auto fy = flux_[lev][1].const_array(mfi);
            amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
                for (int n = 0; n < NCONS; ++n) {
                    if (n == UBX || n == UBY) continue;   // плоскостное B — через CT
                    u(i,j,k,n) -= lx * (fx(i+1,j,k,n) - fx(i,j,k,n))
                                + ly * (fy(i,j+1,k,n) - fy(i,j,k,n));
                }
            });
        }
        // Закон Фарадея по теореме Стокса (ф.(8) Авдеевой–Лукина). Обходим все
        // грани каждого FAB целиком (validbox face-типа): дублируемые на стыках
        // боксов грани получают идентичные значения — ЭДС детерминированы.
        for (MFIter mfi(bface_[lev][0], false); mfi.isValid(); ++mfi) {
            const Box& fb = mfi.validbox();
            auto b  = bface_[lev][0].array(mfi);
            auto ez = emf_[lev].const_array(mfi);
            amrex::LoopOnCpu(fb, [&] (int i, int j, int k) {
                b(i,j,k) -= ly * (ez(i, j+1, k) - ez(i, j, k));
            });
        }
        for (MFIter mfi(bface_[lev][1], false); mfi.isValid(); ++mfi) {
            const Box& fb = mfi.validbox();
            auto b  = bface_[lev][1].array(mfi);
            auto ez = emf_[lev].const_array(mfi);
            amrex::LoopOnCpu(fb, [&] (int i, int j, int k) {
                b(i,j,k) += lx * (ez(i+1, j, k) - ez(i, j, k));
            });
        }
    }
    SyncCellB(lev);
}

// ---------------------------------------------------------------------------
Real MhdAmr::ComputeDt() const
{
    Real dt = std::numeric_limits<Real>::max();
    for (int lev = 0; lev <= finest_level; ++lev) {
        const auto dx = Geom(lev).CellSizeArray();
        for (MFIter mfi(state_[lev]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            auto u = state_[lev].const_array(mfi);
            amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
                double uc[NCONS], q[NPRIM];
                for (int n = 0; n < NCONS; ++n) uc[n] = u(i,j,k,n);
                cons_to_prim(uc, q, cfg_.gamma);
                const double B2 = q[QBX]*q[QBX] + q[QBY]*q[QBY] + q[QBZ]*q[QBZ];
                const double cfx = fast_speed(q[QRHO], q[QP], q[QBX], B2, cfg_.gamma);
                const double cfy = fast_speed(q[QRHO], q[QP], q[QBY], B2, cfg_.gamma);
                dt = std::min(dt, Real(dx[0] / (std::abs(q[QU]) + cfx)));
                dt = std::min(dt, Real(dx[1] / (std::abs(q[QV]) + cfy)));
            });
        }
    }
    ParallelDescriptor::ReduceRealMin(dt);
    return cfg_.cfl * dt;
}

Real MhdAmr::MaxDivB(int lev) const
{
    const auto dx = Geom(lev).CellSizeArray();
    Real m = 0.0;
    for (MFIter mfi(state_[lev]); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto bxf = bface_[lev][0].const_array(mfi);
        auto byf = bface_[lev][1].const_array(mfi);
        amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
            const Real d = (bxf(i+1,j,k) - bxf(i,j,k)) / dx[0]
                         + (byf(i,j+1,k) - byf(i,j,k)) / dx[1];
            m = std::max(m, std::abs(d));
        });
    }
    ParallelDescriptor::ReduceRealMax(m);
    return m;
}

// ---------------------------------------------------------------------------
// Одна стадия Эйлера на всей иерархии: фантомы → потоки/ЭДС на каждом уровне →
// синхронизация ЭДС между уровнями → обновления → average_down.
// ---------------------------------------------------------------------------
void MhdAmr::EulerStage(Real dt)
{
    for (int lev = 0; lev <= finest_level; ++lev) {
        MultiFab tmp(grids[lev], dmap[lev], NCONS, NGROW);
        FillPatchCells(lev, tmp, t_);
        MultiFab::Copy(state_[lev], tmp, 0, 0, NCONS, NGROW);

        Array<MultiFab*, AMREX_SPACEDIM> bf
            {AMREX_D_DECL(&bface_[lev][0], &bface_[lev][1], nullptr)};
        FillPatchFaces(lev, bf, t_);
        FillPhysicalFaceBoundary(lev);
        SyncCellB(lev);

        ComputeFluxesAndEmf(lev);
    }
    SyncEmfAcrossLevels();
    for (int lev = 0; lev <= finest_level; ++lev) ApplyUpdates(lev, dt);
    AverageDownAll();
}

void MhdAmr::AverageDownAll()
{
    for (int lev = finest_level; lev >= 1; --lev) {
        amrex::average_down(state_[lev], state_[lev-1], Geom(lev), Geom(lev-1),
                            0, NCONS, refRatio(lev-1));
        Array<const MultiFab*, AMREX_SPACEDIM> fb
            {AMREX_D_DECL(&bface_[lev][0], &bface_[lev][1], nullptr)};
        Array<MultiFab*, AMREX_SPACEDIM> cb
            {AMREX_D_DECL(&bface_[lev-1][0], &bface_[lev-1][1], nullptr)};
        amrex::average_down_faces(fb, cb, refRatio(lev-1), 0);
        SyncCellB(lev - 1);
    }
}

// SSP-RK2 (метод Хойна): Uⁿ⁺¹ = ½(Uⁿ + (Uⁿ + ΔtL)(+ΔtL)). Каждая стадия —
// CT-обновление; выпуклая комбинация бездивергентных полей бездивергентна,
// поэтому интегратор не нарушает div B = 0 (обоснование выбора — REPORT.md).
void MhdAmr::AdvanceHierarchy(Real dt)
{
    if (cfg_.integrator == Integrator::RK2) {
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::Copy(state0_[lev], state_[lev], 0, 0, NCONS, NGROW);
            for (int d = 0; d < AMREX_SPACEDIM; ++d)
                MultiFab::Copy(bface0_[lev][d], bface_[lev][d], 0, 0, 1, NGROW);
        }
        EulerStage(dt);
        EulerStage(dt);
        for (int lev = 0; lev <= finest_level; ++lev) {
            MultiFab::LinComb(state_[lev], 0.5, state0_[lev], 0,
                              0.5, state_[lev], 0, 0, NCONS, NGROW);
            for (int d = 0; d < AMREX_SPACEDIM; ++d)
                MultiFab::LinComb(bface_[lev][d], 0.5, bface0_[lev][d], 0,
                                  0.5, bface_[lev][d], 0, 0, 1, NGROW);
            SyncCellB(lev);
        }
    } else {
        EulerStage(dt);   // 1-й порядок — как в исходной статье
    }
}

// ---------------------------------------------------------------------------
// Запись результатов: производные величины + плотность/давление/скорость/B/divB.
// Формат: нативный plotfile AMReX (читается ParaView ≥ 5.7 напрямую) или HDF5
// (WriteMultiLevelPlotfileHDF5; нужна сборка AMReX с -DAMReX_HDF5=ON).
// ---------------------------------------------------------------------------
void MhdAmr::WritePlotFile(int step, Real time)
{
    const Vector<std::string> names
        {"rho", "u", "v", "w", "p", "Bx", "By", "Bz", "divB"};
    const int nout = names.size();

    Vector<MultiFab> out(finest_level + 1);
    for (int lev = 0; lev <= finest_level; ++lev) {
        out[lev].define(grids[lev], dmap[lev], nout, 0);
        const auto dx = Geom(lev).CellSizeArray();
        for (MFIter mfi(out[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();
            auto o = out[lev].array(mfi);
            auto u = state_[lev].const_array(mfi);
            auto bxf = bface_[lev][0].const_array(mfi);
            auto byf = bface_[lev][1].const_array(mfi);
            const double gam = cfg_.gamma;
            amrex::LoopOnCpu(bx, [&] (int i, int j, int k) {
                double uc[NCONS], q[NPRIM];
                for (int n = 0; n < NCONS; ++n) uc[n] = u(i,j,k,n);
                cons_to_prim(uc, q, gam);
                o(i,j,k,0)=q[QRHO]; o(i,j,k,1)=q[QU]; o(i,j,k,2)=q[QV];
                o(i,j,k,3)=q[QW];   o(i,j,k,4)=q[QP];
                o(i,j,k,5)=q[QBX];  o(i,j,k,6)=q[QBY]; o(i,j,k,7)=q[QBZ];
                o(i,j,k,8) = (bxf(i+1,j,k)-bxf(i,j,k))/dx[0]
                           + (byf(i,j+1,k)-byf(i,j,k))/dx[1];
            });
        }
    }

    if (amrex::ParallelDescriptor::IOProcessor())
        std::filesystem::create_directories(cfg_.output_dir);
    const std::string fname = cfg_.output_dir + "/" + amrex::Concatenate(cfg_.plot_prefix, step, 5);
    Vector<int> steps(finest_level + 1, step);
#ifdef AMREX_USE_HDF5
    if (cfg_.format == OutputFormat::Hdf5) {
        amrex::WriteMultiLevelPlotfileHDF5(fname, finest_level + 1,
                                           amrex::GetVecOfConstPtrs(out), names,
                                           Geom(), time, steps, refRatio());
        amrex::Print() << "  >> HDF5 plotfile: " << fname << ".h5\n";
        return;
    }
#else
    if (cfg_.format == OutputFormat::Hdf5) {
        amrex::Print() << "  [warn] AMReX собран без HDF5 — пишу нативный plotfile\n";
    }
#endif
    amrex::WriteMultiLevelPlotfile(fname, finest_level + 1,
                                   amrex::GetVecOfConstPtrs(out), names,
                                   Geom(), time, steps, refRatio());
    amrex::Print() << "  >> plotfile: " << fname << "\n";
}

// ---------------------------------------------------------------------------
void MhdAmr::InitData()
{
    InitFromScratch(0.0);    // строит иерархию: MakeNewLevelFromScratch + ErrorEst
    AverageDownAll();
    WritePlotFile(0, 0.0);
}

void MhdAmr::Evolve()
{
    Real next_plot_t = (cfg_.plot_dt > 0) ? cfg_.plot_dt : -1.0;

    while (t_ < cfg_.t_max && step_ < cfg_.max_steps) {
        if (max_level > 0 && cfg_.regrid_int > 0 &&
            step_ > 0 && step_ % cfg_.regrid_int == 0) {
            regrid(0, t_);
        }
        Real dt = std::min(ComputeDt(), cfg_.t_max - t_);
        AdvanceHierarchy(dt);
        t_ += dt; ++step_;

        if (step_ % cfg_.diag_int == 0) {
            Real divb = 0.0;
            for (int lev = 0; lev <= finest_level; ++lev)
                divb = std::max(divb, MaxDivB(lev));
            amrex::Print() << "step " << step_ << "  t=" << t_
                           << "  dt=" << dt << "  max|divB|=" << divb
                           << "  levels=" << finest_level + 1 << "\n";
        }
        const bool plot_now =
            (cfg_.plot_int > 0 && step_ % cfg_.plot_int == 0) ||
            (cfg_.plot_dt > 0 && t_ >= next_plot_t - 1e-14);
        if (plot_now) {
            WritePlotFile(step_, t_);
            if (cfg_.plot_dt > 0) next_plot_t += cfg_.plot_dt;
        }
    }
    WritePlotFile(step_, t_);
    amrex::Print() << "Evolve finished: " << step_ << " steps, t=" << t_ << "\n";
}

} // namespace mhd
