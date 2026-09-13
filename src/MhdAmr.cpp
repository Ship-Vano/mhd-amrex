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
#include <AMReX_GpuAtomic.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_Reduce.H>
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

// Device implementation of ext_dir.  AMReX's GpuBndryFuncFab handles normal
// outflow/reflect fill first; this functor overwrites only frozen ext_dir data.
struct ExtDirGpuFill {
    DeviceProblem problem;
    Real gamma;

    AMREX_GPU_DEVICE void operator() (const IntVect& iv, Array4<Real> const& arr,
                                      int dcomp, int numcomp, GeometryData const& geom,
                                      Real, const BCRec* bcr, int, int orig_comp) const noexcept
    {
        const Box& domain = geom.Domain();
        bool ext = false;
        for (int dim = 0; dim < AMREX_SPACEDIM; ++dim) {
            const int idx = (dim == 0) ? iv[0] : iv[1];
            if ((idx < domain.smallEnd(dim) && bcr[0].lo(dim) == BCType::ext_dir) ||
                (idx > domain.bigEnd(dim) && bcr[0].hi(dim) == BCType::ext_dir)) ext = true;
        }
        if (!ext) return;
        const Real dx0 = geom.CellSize(0), dx1 = geom.CellSize(1);
        const Real x = geom.ProbLo(0) + (Real(iv[0]) + Real(0.5))*dx0;
        const Real y = geom.ProbLo(1) + (Real(iv[1]) + Real(0.5))*dx1;
        Real q[NPRIM] = {}, uc[NCONS];
        problem.prim(x, y, q);
        if (problem.direct_b()) {
            q[QBX] = Real(0.5)*(problem.bx_face(x-Real(0.5)*dx0,y)+problem.bx_face(x+Real(0.5)*dx0,y));
            q[QBY] = Real(0.5)*(problem.by_face(x,y-Real(0.5)*dx1)+problem.by_face(x,y+Real(0.5)*dx1));
        } else {
            q[QBX] = (problem.az(x,y+Real(0.5)*dx1)-problem.az(x,y-Real(0.5)*dx1))/dx1;
            q[QBY] =-(problem.az(x+Real(0.5)*dx0,y)-problem.az(x-Real(0.5)*dx0,y))/dx0;
        }
        prim_to_cons(q, uc, gamma);
        // dcomp -- номер компоненты в приёмнике, orig_comp -- в состоянии.
        // Сейчас все вызовы заполняют весь диапазон, и оба равны нулю, но брать
        // источник по dcomp верно лишь по совпадению: заполнение подмножества
        // компонент молча прочитало бы не те.
        for (int n = 0; n < numcomp && orig_comp+n < NCONS; ++n)
            arr(iv,dcomp+n) = uc[orig_comp+n];
    }
};

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
    // Регистр потоков строится на конкретных BoxArray/DistributionMapping двух
    // соседних уровней, поэтому любая перестройка сетки его обесценивает.
    flux_reg_stale_ = true;
}

void MhdAmr::ClearLevel(int lev)
{
    state_[lev].clear(); state0_[lev].clear(); emf_[lev].clear();
    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        bface_[lev][d].clear(); bface0_[lev][d].clear(); flux_[lev][d].clear();
    }
    flux_reg_stale_ = true;
}

// НУ уровня: грани — через векторный потенциал Az (div B = 0 машинно),
// клеточные B — RT0-интерполяция граней (как в схеме Авдеевой–Лукина).
void MhdAmr::InitLevelData(int lev)
{
    BL_PROFILE("MhdAmr::InitLevelData");
    const auto problo = Geom(lev).ProbLoArray();
    const auto dx     = Geom(lev).CellSizeArray();
    const DeviceProblem P = prob_;
    const Real gam  = cfg_.gamma;

    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        for (MFIter mfi(bface_[lev][d]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.fabbox();      // вместе с фантомами
            auto b = bface_[lev][d].array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                const Real xf = problo[0] + i * dx[0] + (d == 0 ? Real(0.0) : Real(0.5)*dx[0]);
                const Real yf = problo[1] + j * dx[1] + (d == 1 ? Real(0.0) : Real(0.5)*dx[1]);
                if (P.direct_b()) {
                    b(i,j,k) = (d == 0) ? P.bx_face(xf, yf) : P.by_face(xf, yf);
                } else if (d == 0) {
                    b(i,j,k) =  (P.az(xf, yf + Real(0.5)*dx[1]) - P.az(xf, yf - Real(0.5)*dx[1])) / dx[1];
                } else {
                    b(i,j,k) = -(P.az(xf + Real(0.5)*dx[0], yf) - P.az(xf - Real(0.5)*dx[0], yf)) / dx[0];
                }
            });
        }
    }
    for (MFIter mfi(state_[lev]); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.fabbox();
        auto u   = state_[lev].array(mfi);
        auto bxf = bface_[lev][0].const_array(mfi);
        auto byf = bface_[lev][1].const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real x = problo[0] + (Real(i) + Real(0.5)) * dx[0];
            const Real y = problo[1] + (Real(j) + Real(0.5)) * dx[1];
            Real q[NPRIM] = {}, uc[NCONS];
            P.prim(x, y, q);
            q[QBX] = Real(0.5) * (bxf(i,j,k) + bxf(i+1,j,k));
            q[QBY] = Real(0.5) * (byf(i,j,k) + byf(i,j+1,k));
            prim_to_cons(q, uc, gam);
            for (int n = 0; n < NCONS; ++n) u(i,j,k,n) = uc[n];
        });
    }
}

// Пересчёт клеточного B после интерполяции при regrid.
//
// Полная энергия переносится cell_cons_interp как консервативный скаляр, а
// граневое поле -- face_divfree_interp. Их магнитные части не совпадают:
// измерено, что ½|B|² скачет на ~1e-5 относительных на каждом перестроении.
// Если оставить UENE неизменной, эту разницу целиком поглощает тепловая
// энергия, то есть давление. Здесь UENE сдвигается на изменение ½|B|², и тогда
// точно сохраняется газовая (тепловая + кинетическая) энергия, а полная --
// меняется на ошибку интерполяции поля. Что именно сохранять -- решение D-007.
void MhdAmr::SyncCellBAfterRegrid(int lev)
{
    if (!cfg_.regrid_preserve_pressure) { SyncCellB(lev); return; }

    // ВАЖНО: AmrCore выставляет grids[lev]/dmap[lev] уже ПОСЛЕ возврата из
    // MakeNewLevelFromCoarse/RemakeLevel, поэтому раскладку берём у самого
    // state_[lev], который только что определил AllocLevel.
    MultiFab me_before(state_[lev].boxArray(), state_[lev].DistributionMap(), 1, 0);
    for (MFIter mfi(state_[lev]); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto u = state_[lev].const_array(mfi);
        auto m = me_before.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            m(i,j,k) = Real(0.5) * (u(i,j,k,UBX)*u(i,j,k,UBX)
                                  + u(i,j,k,UBY)*u(i,j,k,UBY)
                                  + u(i,j,k,UBZ)*u(i,j,k,UBZ));
        });
    }
    SyncCellB(lev);
    for (MFIter mfi(state_[lev]); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto u = state_[lev].array(mfi);
        auto m = me_before.const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
            const Real me_after = Real(0.5) * (u(i,j,k,UBX)*u(i,j,k,UBX)
                                             + u(i,j,k,UBY)*u(i,j,k,UBY)
                                             + u(i,j,k,UBZ)*u(i,j,k,UBZ));
            u(i,j,k,UENE) += me_after - m(i,j,k);
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
        const ExtDirGpuFill fill {prob_, cfg_.gamma};
        PhysBCFunct<GpuBndryFuncFab<ExtDirGpuFill>> cbc(Geom(lev-1), bcrec_, GpuBndryFuncFab<ExtDirGpuFill>(fill));
        PhysBCFunct<GpuBndryFuncFab<ExtDirGpuFill>> fbc(Geom(lev),   bcrec_, GpuBndryFuncFab<ExtDirGpuFill>(fill));
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
    SyncCellBAfterRegrid(lev);
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
    SyncCellBAfterRegrid(lev);
}

// ---------------------------------------------------------------------------
// Критерий измельчения: относительный градиент плотности и/или ток jz = (∇×B)z
// ---------------------------------------------------------------------------
void MhdAmr::ErrorEst(int lev, TagBoxArray& tags, Real /*time*/, int /*ngrow*/)
{
    BL_PROFILE("MhdAmr::ErrorEst");
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
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
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
    BL_PROFILE("MhdAmr::FillPatchCells");
    const ExtDirGpuFill fill {prob_, cfg_.gamma};
    PhysBCFunct<GpuBndryFuncFab<ExtDirGpuFill>> fbc(Geom(lev), bcrec_, GpuBndryFuncFab<ExtDirGpuFill>(fill));
    if (lev == 0) {
        amrex::FillPatchSingleLevel(mf, time, {&state_[0]}, {time}, 0, 0, NCONS,
                                    Geom(0), fbc, 0);
    } else {
        PhysBCFunct<GpuBndryFuncFab<ExtDirGpuFill>> cbc(Geom(lev-1), bcrec_, GpuBndryFuncFab<ExtDirGpuFill>(fill));
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
    BL_PROFILE("MhdAmr::FillPatchFaces");
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
    const int lo0 = static_cast<int>(cfg_.bc_xlo), lo1 = static_cast<int>(cfg_.bc_ylo);
    const int hi0 = static_cast<int>(cfg_.bc_xhi), hi1 = static_cast<int>(cfg_.bc_yhi);
    const DeviceProblem P = prob_;

    for (int d = 0; d < AMREX_SPACEDIM; ++d) {
        bface_[lev][d].FillBoundary(Geom(lev).periodicity());
        const Box fdomain = amrex::convert(domain, IntVect::TheDimensionVector(d));
        for (MFIter mfi(bface_[lev][d]); mfi.isValid(); ++mfi) {
            const Box& fb = mfi.fabbox();
            auto b = bface_[lev][d].array(mfi);
            amrex::ParallelFor(fb, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                for (int dim = 0; dim < AMREX_SPACEDIM; ++dim) {
                    const int dlo = fdomain.smallEnd(dim), dhi = fdomain.bigEnd(dim);
                    BcType bt;
                    const int idx = (dim == 0) ? i : j;
                    if      (idx < dlo) bt = static_cast<BcType>((dim == 0) ? lo0 : lo1);
                    else if (idx > dhi) bt = static_cast<BcType>((dim == 0) ? hi0 : hi1);
                    else continue;
                    if (bt == BcType::Periodic) continue;
                    if (bt == BcType::Outflow) {
                        const int si = (dim == 0) ? ((i < dlo) ? dlo : ((i > dhi) ? dhi : i)) : i;
                        const int sj = (dim == 1) ? ((j < dlo) ? dlo : ((j > dhi) ? dhi : j)) : j;
                        b(i,j,k) = b(si,sj,0);
                    } else if (bt == BcType::Reflect) {
                        // нормальная к границе компонента — нечётная,
                        // касательная — чётная (зеркальное отражение поля)
                        const bool normal = (dim == d);
                        int si = i, sj = j;
                        const int src = normal ? ((idx < dlo) ? 2*dlo-idx : 2*dhi-idx)
                                               : ((idx < dlo) ? 2*dlo-idx-1 : 2*dhi-idx+1);
                        if (dim == 0) si = src; else sj = src;
                        b(i,j,k) = (normal ? Real(-1.0) : Real(1.0)) * b(si,sj,0);
                    } else {  // Dirichlet: «исторические» значения из НУ
                        const Real xf = problo[0] + i * dx[0] + (d == 0 ? Real(0) : Real(0.5)*dx[0]);
                        const Real yf = problo[1] + j * dx[1] + (d == 1 ? Real(0) : Real(0.5)*dx[1]);
                        if (P.direct_b()) {
                            b(i,j,k) = (d == 0) ? P.bx_face(xf, yf) : P.by_face(xf, yf);
                        } else if (d == 0) {
                            b(i,j,k) =  (P.az(xf, yf+Real(0.5)*dx[1]) - P.az(xf, yf-Real(0.5)*dx[1])) / dx[1];
                        } else {
                            b(i,j,k) = -(P.az(xf+Real(0.5)*dx[0], yf) - P.az(xf-Real(0.5)*dx[0], yf)) / dx[0];
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
    BL_PROFILE("MhdAmr::SyncCellB");
    for (MFIter mfi(state_[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        const Box bx = mfi.growntilebox(NGROW - 1);
        auto u   = state_[lev].array(mfi);
        auto bxf = bface_[lev][0].const_array(mfi);
        auto byf = bface_[lev][1].const_array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
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
    BL_PROFILE("MhdAmr::ComputeFluxesAndEmf");
    const Limiter lim = cfg_.limiter;
    const EmfAveraging emode = cfg_.emf;
    const Real gam = cfg_.gamma;

    // DeviceScalar is device-resident on CUDA and ordinary host storage on CPU.
    // HostDevice atomics preserve the diagnostic invariant in both execution spaces:
    // on device they are real atomics, on host `#pragma omp atomic update`. Именно
    // поэтому OpenMP-распараллеливание по тайлам ниже безопасно для счётчиков --
    // ради него прежняя версия и держала здесь reduction(+:...).
    Gpu::DeviceScalar<Long> level_fallbacks(0);
    Gpu::DeviceScalar<Long> level_floors(0);
    Long* const fallback_counter = level_fallbacks.dataPtr();
    Long* const floor_counter = level_floors.dataPtr();
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(state_[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
        auto u   = state_[lev].const_array(mfi);
        auto bxf = bface_[lev][0].const_array(mfi);
        auto byf = bface_[lev][1].const_array(mfi);
        auto fx  = flux_[lev][0].array(mfi);
        auto fy  = flux_[lev][1].array(mfi);
        auto ez  = emf_[lev].array(mfi);

        // примитивы в ячейке (i,j) по запросу
        auto qprim = [=] AMREX_GPU_DEVICE (int i, int j, Real* q) noexcept -> int {
            Real uc[NCONS];
            for (int n = 0; n < NCONS; ++n) uc[n] = u(i, j, 0, n);
            int floors = 0;
            cons_to_prim(uc, q, gam, Limits{}, &floors);
            return floors;
        };

        // --- x-потоки: грани валидной области + 1 слой (нужно узловым ЭДС) --
        {
            const Box xbx = mfi.grownnodaltilebox(0, 1);
            amrex::For(xbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                Real qm[NPRIM], q0[NPRIM], qp[NPRIM], qq[NPRIM];
                Real qL[NPRIM], qR[NPRIM], f[NCONS];
                int floors = qprim(i-2,j,qm)+qprim(i-1,j,q0)+qprim(i,j,qp)+qprim(i+1,j,qq);
                for (int n = 0; n < NPRIM; ++n) {
                    qL[n] = face_value_plus (qm[n], q0[n], qp[n], lim);
                    qR[n] = face_value_minus(q0[n], qp[n], qq[n], lim);
                }
                int fallback = 0;
                hlld_flux(qL, qR, bxf(i,j,k), f, gam, Limits{}, &fallback);
                for (int n = 0; n < NCONS; ++n) fx(i,j,k,n) = f[n];
                if (fallback) HostDevice::Atomic::Add(fallback_counter, Long(fallback));
                if (floors) HostDevice::Atomic::Add(floor_counter, Long(floors));
            });
        }
        // --- y-потоки: локальный поворот осей (u'=v, v'=−u, Bx'=By, By'=−Bx) -
        {
            const Box ybx = mfi.grownnodaltilebox(1, 1);
            amrex::For(ybx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                Real qm[NPRIM], q0[NPRIM], qp[NPRIM], qq[NPRIM];
                Real qL[NPRIM], qR[NPRIM], rL[NPRIM], rR[NPRIM], f[NCONS];
                int floors = qprim(i,j-2,qm)+qprim(i,j-1,q0)+qprim(i,j,qp)+qprim(i,j+1,qq);
                for (int n = 0; n < NPRIM; ++n) {
                    qL[n] = face_value_plus (qm[n], q0[n], qp[n], lim);
                    qR[n] = face_value_minus(q0[n], qp[n], qq[n], lim);
                }
                auto rot = [] AMREX_GPU_DEVICE (const Real* q, Real* r) noexcept {
                    r[QRHO]=q[QRHO]; r[QP]=q[QP]; r[QW]=q[QW]; r[QBZ]=q[QBZ];
                    r[QU]=q[QV]; r[QV]=-q[QU]; r[QBX]=q[QBY]; r[QBY]=-q[QBX];
                };
                rot(qL, rL); rot(qR, rR);
                int fallback = 0;
                hlld_flux(rL, rR, byf(i,j,k), f, gam, Limits{}, &fallback);
                fy(i,j,k,URHO)=f[URHO]; fy(i,j,k,UENE)=f[UENE];
                fy(i,j,k,UMZ)=f[UMZ];   fy(i,j,k,UBZ)=f[UBZ];
                fy(i,j,k,UMX)=-f[UMY];  fy(i,j,k,UMY)=f[UMX];
                fy(i,j,k,UBX)=-f[UBY];  fy(i,j,k,UBY)=f[UBX];   // fy[UBX] = +Ez
                if (fallback) HostDevice::Atomic::Add(fallback_counter, Long(fallback));
                if (floors) HostDevice::Atomic::Add(floor_counter, Long(floors));
            });
        }
        // --- узловые ЭДС Ez(i−1/2, j−1/2): усреднение ЭДС примыкающих граней —
        // декартов аналог усреднения по рёбрам, сходящимся в вершине (ф.(8)
        // статьи Авдеевой–Лукина); GardinerStone добавляет клеточную поправку.
        // ЭДС на x-грани: Ez = −Fx[UBY]; на y-грани: Ez = +Fy[UBX].
        {
            const Box nbx = mfi.tilebox(IntVect::TheNodeVector());
            amrex::For(nbx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                Real qmm[NPRIM], qpm[NPRIM], qmp[NPRIM], qpp[NPRIM];
                const int floors = qprim(i-1,j-1,qmm)+qprim(i,j-1,qpm)
                                 + qprim(i-1,j,qmp)+qprim(i,j,qpp);
                ez(i,j,k) = corner_emf(-fx(i, j-1, k, UBY), -fx(i, j, k, UBY),
                                        fy(i-1, j, k, UBX),  fy(i, j, k, UBX),
                                        cell_emf_z(qmm), cell_emf_z(qpm),
                                        cell_emf_z(qmp), cell_emf_z(qpp), emode);
                if (floors) HostDevice::Atomic::Add(floor_counter, Long(floors));
            });
        }
    }
    hlld_fallbacks_ += level_fallbacks.dataValue();
    floor_events_   += level_floors.dataValue();
}

// Число ячеек с ρ ≤ small_rho либо p ≤ small_pres на всей иерархии (после
// текущего шага). Молчаливый floor в cons_to_prim не считается доказательством
// устойчивости — этот счётчик делает его наблюдаемым (инвариант 1, AGENTS.md).
amrex::Long MhdAmr::CountNonPositiveCells() const
{
    const Limits lim;
    Gpu::DeviceScalar<Long> counter(0);
    Long* const count = counter.dataPtr();
    for (int lev = 0; lev <= finest_level; ++lev) {
        for (MFIter mfi(state_[lev]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            auto u = state_[lev].const_array(mfi);
            const Real gam = cfg_.gamma;
            amrex::For(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                Real uc[NCONS];
                for (int c = 0; c < NCONS; ++c) uc[c] = u(i,j,k,c);
                if (!(uc[URHO] > lim.small_rho) ||
                    !(pressure_from_cons(uc, gam) > lim.small_pres)) HostDevice::Atomic::Add(count, Long(1));
            });
        }
    }
    amrex::Long n = counter.dataValue();
    ParallelDescriptor::ReduceLongSum(n);
    return n;
}

// Инжекция узловых ЭДС мелкого уровня в совпадающие узлы грубого. Поскольку
// расчёт ведётся БЕЗ подциклирования (общий Δt), после инжекции обновление
// грубой грани на границе уровней в точности равно среднему обновлений
// накрывающих её мелких граней → average_down_faces не вносит дивергенцию,
// и div B = 0 сохраняется на всей иерархии (вывод — в REPORT.md).
void MhdAmr::SyncEmfAcrossLevels()
{
    BL_PROFILE("MhdAmr::SyncEmfAcrossLevels");
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
    BL_PROFILE("MhdAmr::ApplyUpdates");
    const auto dx = Geom(lev).CellSizeArray();
    const Real lx = dt / dx[0], ly = dt / dx[1];

#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    {
        // Газовые величины: годуновский баланс потоков (ур. (1)–(3) статьи)
        for (MFIter mfi(state_[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.tilebox();
            auto u  = state_[lev].array(mfi);
            auto fx = flux_[lev][0].const_array(mfi);
            auto fy = flux_[lev][1].const_array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
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
            amrex::ParallelFor(fb, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                b(i,j,k) -= ly * (ez(i, j+1, k) - ez(i, j, k));
            });
        }
        for (MFIter mfi(bface_[lev][1], false); mfi.isValid(); ++mfi) {
            const Box& fb = mfi.validbox();
            auto b  = bface_[lev][1].array(mfi);
            auto ez = emf_[lev].const_array(mfi);
            amrex::ParallelFor(fb, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                b(i,j,k) += lx * (ez(i+1, j, k) - ez(i, j, k));
            });
        }
    }
    SyncCellB(lev);
}

// ---------------------------------------------------------------------------
Real MhdAmr::ComputeDt() const
{
    BL_PROFILE("MhdAmr::ComputeDt");
    Real dt = std::numeric_limits<Real>::max();
    const Real gam = cfg_.gamma;
    for (int lev = 0; lev <= finest_level; ++lev) {
        const auto dx = Geom(lev).CellSizeArray();
        // Объекты редукции создаются один раз на уровень, а не на каждый бокс:
        // на GPU конструктор выделяет device-память, а чтение .value() -- это
        // синхронизация, и внутри цикла по боксам они шли бы на каждый бокс.
        ReduceOps<ReduceOpMin> reduce_op;
        ReduceData<Real> reduce_data(reduce_op);
        using ReduceTuple = ReduceData<Real>::Type;
        for (MFIter mfi(state_[lev]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            auto u = state_[lev].const_array(mfi);
            reduce_op.eval(bx, reduce_data, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> ReduceTuple {
                Real uc[NCONS], q[NPRIM];
                for (int n = 0; n < NCONS; ++n) uc[n] = u(i,j,k,n);
                cons_to_prim(uc, q, gam);
                const Real B2 = q[QBX]*q[QBX] + q[QBY]*q[QBY] + q[QBZ]*q[QBZ];
                const Real cfx = fast_speed(q[QRHO], q[QP], q[QBX], B2, gam);
                const Real cfy = fast_speed(q[QRHO], q[QP], q[QBY], B2, gam);
                const Real dtx = dx[0] / (std::abs(q[QU]) + cfx);
                const Real dty = dx[1] / (std::abs(q[QV]) + cfy);
                return { dtx < dty ? dtx : dty };
            });
        }
        dt = std::min(dt, amrex::get<0>(reduce_data.value(reduce_op)));
    }
    ParallelDescriptor::ReduceRealMin(dt);
    return cfg_.cfl * dt;
}

Real MhdAmr::MaxDivB(int lev) const
{
    const auto dx = Geom(lev).CellSizeArray();
    Real m = 0.0;
    ReduceOps<ReduceOpMax> reduce_op;
    ReduceData<Real> reduce_data(reduce_op);
    using ReduceTuple = ReduceData<Real>::Type;
    for (MFIter mfi(state_[lev]); mfi.isValid(); ++mfi) {
        const Box& bx = mfi.validbox();
        auto bxf = bface_[lev][0].const_array(mfi);
        auto byf = bface_[lev][1].const_array(mfi);
        reduce_op.eval(bx, reduce_data, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> ReduceTuple {
            const Real d = (bxf(i+1,j,k) - bxf(i,j,k)) / dx[0]
                         + (byf(i,j+1,k) - byf(i,j,k)) / dx[1];
            return { std::abs(d) };
        });
    }
    m = std::max(m, amrex::get<0>(reduce_data.value(reduce_op)));
    ParallelDescriptor::ReduceRealMax(m);
    return m;
}

// ---------------------------------------------------------------------------
// Одна стадия Эйлера на всей иерархии: фантомы → потоки/ЭДС на каждом уровне →
// синхронизация ЭДС между уровнями → обновления → average_down.
// ---------------------------------------------------------------------------
void MhdAmr::EulerStage(Real dt)
{
    BL_PROFILE("MhdAmr::EulerStage");
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
    const bool do_reflux = cfg_.reflux && finest_level > 0;
    if (do_reflux) {
        if (flux_reg_stale_) DefineFluxRegisters();
        AccumulateFluxRegisters(dt);
    }
    SyncEmfAcrossLevels();
    for (int lev = 0; lev <= finest_level; ++lev) ApplyUpdates(lev, dt);
    // Порядок важен: рефлюкс правит грубые ячейки У СТЫКА (не перекрытые),
    // average_down затем перезаписывает только перекрытые — они не конфликтуют.
    if (do_reflux) RefluxAll();
    AverageDownAll();
}

// ---------------------------------------------------------------------------
// Согласование газовых потоков на границах уровней (reflux).
//
// Шаг не подциклируется: Δt одинаков на всех уровнях, поэтому регистр потоков
// накапливает за одну стадию грубый вклад (CrseAdd) и мелкий (FineAdd), а
// Reflux добавляет разность к грубым ячейкам, примыкающим к мелкой сетке.
// Величины UBX, UBY исключены: плоскостное поле ведёт CT через узловые ЭДС,
// и рефлюкс газового потока разрушил бы дискретную бездивергентность.
// ---------------------------------------------------------------------------
void MhdAmr::DefineFluxRegisters()
{
    flux_reg_.resize(max_level + 1);
    for (int lev = 1; lev <= finest_level; ++lev) {
        flux_reg_[lev] = std::make_unique<amrex::YAFluxRegister>(
            grids[lev], grids[lev-1], dmap[lev], dmap[lev-1],
            Geom(lev), Geom(lev-1), refRatio(lev-1), lev, NCONS);
    }
    for (int lev = finest_level + 1; lev <= max_level; ++lev) flux_reg_[lev].reset();
    flux_reg_stale_ = false;
}

void MhdAmr::AccumulateFluxRegisters(Real dt)
{
    BL_PROFILE("MhdAmr::AccumulateFluxRegisters");
    for (int lev = 1; lev <= finest_level; ++lev) flux_reg_[lev]->reset();

#if defined(AMREX_USE_GPU)
    constexpr RunOn flux_register_run_on = RunOn::Gpu;
#else
    constexpr RunOn flux_register_run_on = RunOn::Cpu;
#endif
    for (int lev = 0; lev <= finest_level; ++lev) {
        const auto dxa = Geom(lev).CellSizeArray();
        const Real dx[AMREX_SPACEDIM] = {AMREX_D_DECL(dxa[0], dxa[1], dxa[2])};
        // MFIter обходит клеточную раскладку уровня lev: для CrseAdd это
        // грубая сторона регистра lev+1, для FineAdd — мелкая сторона lev.
        for (MFIter mfi(state_[lev]); mfi.isValid(); ++mfi) {
            const std::array<FArrayBox const*, AMREX_SPACEDIM> f
                {AMREX_D_DECL(&flux_[lev][0][mfi], &flux_[lev][1][mfi], &flux_[lev][2][mfi])};
            if (lev < finest_level)
                flux_reg_[lev+1]->CrseAdd(mfi, f, dx, dt, flux_register_run_on);
            if (lev > 0)
                flux_reg_[lev]->FineAdd(mfi, f, dx, dt, flux_register_run_on);
        }
    }
}

void MhdAmr::RefluxAll()
{
    BL_PROFILE("MhdAmr::RefluxAll");
    for (int lev = 1; lev <= finest_level; ++lev) {
        flux_reg_[lev]->Reflux(state_[lev-1], URHO, URHO, UENE - URHO + 1);  // ρ, ρv, e
        flux_reg_[lev]->Reflux(state_[lev-1], UBZ,  UBZ,  1);                // Bz
    }
}

// Диапазоны rho и p по иерархии, без перекрытых грубых ячеек. Нужны и для
// сравнения с литературными цветовыми шкалами (ОТ, цилиндр, взрыв), и для
// теста постоянного состояния, где max-min обязан быть строго нулём.
// Σ по иерархии магнитной (½|B|²) и тепловой (e − ½ρ|v|² − ½|B|²) энергий.
// Полная энергия интерполируется при regrid как консервативный скаляр, а поле B
// -- отдельным бездивергентным оператором; эти две части могут разъехаться, и
// разделение позволяет это увидеть, а не предполагать.
amrex::Real MhdAmr::TotalEnergyPart(int which) const
{
    Real total = 0.0;
    for (int lev = 0; lev <= finest_level; ++lev) {
        const auto dx = Geom(lev).CellSizeArray();
        const Real dv = dx[0] * dx[1];
        iMultiFab covered;
        const bool has_finer = (lev < finest_level);
        if (has_finer)
            covered = amrex::makeFineMask(grids[lev], dmap[lev], grids[lev+1],
                                          refRatio(lev), 0, 1);
        // Объект редукции -- на уровень, а не на бокс: на GPU его конструктор
        // выделяет device-память, а reduce_data.value() синхронизируется с
        // устройством, и внутри цикла по боксам это происходило бы на каждом.
        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<Real> reduce_data(reduce_op);
        using ReduceTuple = ReduceData<Real>::Type;
        for (MFIter mfi(state_[lev]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            auto u = state_[lev].const_array(mfi);
            Array4<const int> cov{};
            if (has_finer) cov = covered.const_array(mfi);
            reduce_op.eval(bx, reduce_data, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> ReduceTuple {
                if (has_finer && cov(i,j,k)) return { Real(0.0) };
                const Real rho = u(i,j,k,URHO);
                const Real me = Real(0.5) * (u(i,j,k,UBX)*u(i,j,k,UBX)
                                           + u(i,j,k,UBY)*u(i,j,k,UBY)
                                           + u(i,j,k,UBZ)*u(i,j,k,UBZ));
                const Real ke = Real(0.5) * (u(i,j,k,UMX)*u(i,j,k,UMX)
                                           + u(i,j,k,UMY)*u(i,j,k,UMY)
                                           + u(i,j,k,UMZ)*u(i,j,k,UMZ)) / rho;
                return { (which == 0) ? me : (u(i,j,k,UENE) - ke - me) };
            });
        }
        total += amrex::get<0>(reduce_data.value(reduce_op)) * dv;
    }
    ParallelDescriptor::ReduceRealSum(total);
    return total;
}

void MhdAmr::PrintStateRanges() const
{
    Real rho_lo =  std::numeric_limits<Real>::max();
    Real rho_hi = -std::numeric_limits<Real>::max();
    Real p_lo   =  std::numeric_limits<Real>::max();
    Real p_hi   = -std::numeric_limits<Real>::max();
    for (int lev = 0; lev <= finest_level; ++lev) {
        iMultiFab covered;
        const bool has_finer = (lev < finest_level);
        if (has_finer) {
            covered = amrex::makeFineMask(grids[lev], dmap[lev], grids[lev+1],
                                          refRatio(lev), 0, 1);
        }
        const Real gam = cfg_.gamma;
        const Real hi_identity = -std::numeric_limits<Real>::max();
        const Real lo_identity = std::numeric_limits<Real>::max();
        // См. ComputeDt: объект редукции создаётся на уровень, не на бокс.
        ReduceOps<ReduceOpMin,ReduceOpMax,ReduceOpMin,ReduceOpMax> reduce_op;
        ReduceData<Real,Real,Real,Real> reduce_data(reduce_op);
        using ReduceTuple = ReduceData<Real,Real,Real,Real>::Type;
        for (MFIter mfi(state_[lev]); mfi.isValid(); ++mfi) {
            const Box& bx = mfi.validbox();
            auto u = state_[lev].const_array(mfi);
            Array4<const int> cov{};
            if (has_finer) cov = covered.const_array(mfi);
            reduce_op.eval(bx, reduce_data, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept -> ReduceTuple {
                if (has_finer && cov(i,j,k)) return {lo_identity,hi_identity,lo_identity,hi_identity};
                Real uc[NCONS], q[NPRIM];
                for (int n = 0; n < NCONS; ++n) uc[n] = u(i,j,k,n);
                cons_to_prim(uc, q, gam);
                return {q[QRHO],q[QRHO],q[QP],q[QP]};
            });
        }
        const auto values = reduce_data.value(reduce_op);
        rho_lo = std::min(rho_lo, amrex::get<0>(values));
        rho_hi = std::max(rho_hi, amrex::get<1>(values));
        p_lo   = std::min(p_lo,   amrex::get<2>(values));
        p_hi   = std::max(p_hi,   amrex::get<3>(values));
    }
    ParallelDescriptor::ReduceRealMin(rho_lo);
    ParallelDescriptor::ReduceRealMax(rho_hi);
    ParallelDescriptor::ReduceRealMin(p_lo);
    ParallelDescriptor::ReduceRealMax(p_hi);
    amrex::Print() << "ranges: rho_min=" << rho_lo << " rho_max=" << rho_hi
                   << " p_min=" << p_lo << " p_max=" << p_hi << "\n";
}

Real MhdAmr::TotalConserved(int comp) const
{
    Real total = 0.0;
    for (int lev = 0; lev <= finest_level; ++lev) {
        const auto dx = Geom(lev).CellSizeArray();
        MultiFab tmp(grids[lev], dmap[lev], 1, 0);
        MultiFab::Copy(tmp, state_[lev], comp, 0, 1, 0);
        if (lev < finest_level) {
            // Перекрытые грубые ячейки не считаем — их представляет мелкий уровень.
            MultiFab mask = amrex::makeFineMask(grids[lev], dmap[lev], grids[lev+1],
                                                refRatio(lev), Real(1.0), Real(0.0));
            MultiFab::Multiply(tmp, mask, 0, 0, 1, 0);
        }
        total += tmp.sum(0) * dx[0] * dx[1];   // MultiFab::sum уже редуцирует по рангам
    }
    return total;
}

void MhdAmr::AverageDownAll()
{
    BL_PROFILE("MhdAmr::AverageDownAll");
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
    BL_PROFILE("MhdAmr::WritePlotFile");
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
            const Real gam = cfg_.gamma;
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {
                Real uc[NCONS], q[NPRIM];
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
    // Опорные интегралы для диагностики консервативности (см. Evolve).
    for (int n = 0; n < NCONS; ++n) conserved0_[n] = TotalConserved(n);
    conserved0_set_ = true;
    if (cfg_.write_plotfiles) WritePlotFile(0, 0.0);
}

void MhdAmr::Evolve()
{
    Real next_plot_t = (cfg_.plot_dt > 0) ? cfg_.plot_dt : -1.0;
    int last_plot_step = -1;

    while (t_ < cfg_.t_max && step_ < cfg_.max_steps) {
        if (max_level > 0 && cfg_.regrid_int > 0 &&
            step_ > 0 && step_ % cfg_.regrid_int == 0) {
            // Перестроение сетки обязано сохранять интегралы: интерполяция
            // cell_cons_interp консервативна по построению. Измеряем скачок,
            // а не предполагаем его отсутствие.
            const Real mass_before = TotalConserved(URHO);
            const Real ene_before  = TotalConserved(UENE);
            const Real mag_before  = TotalEnergyPart(0);
            const Real th_before   = TotalEnergyPart(1);
            regrid(0, t_);
            const Real mass_after = TotalConserved(URHO);
            const Real ene_after  = TotalConserved(UENE);
            const Real mag_after  = TotalEnergyPart(0);
            const Real th_after   = TotalEnergyPart(1);
            if (std::abs(mag_after - mag_before) / std::abs(mag_before) > 1.0e-14)
                amrex::Print() << "  regrid step " << step_
                               << ": d(total E)/E=" << std::abs(ene_after-ene_before)/std::abs(ene_before)
                               << "  d(magnetic)/mag=" << std::abs(mag_after-mag_before)/std::abs(mag_before)
                               << "  d(thermal)/th=" << std::abs(th_after-th_before)/std::abs(th_before)
                               << "\n";
            const Real dm = std::abs(mass_after - mass_before) / std::abs(mass_before);
            const Real de = std::abs(ene_after - ene_before) / std::abs(ene_before);
            regrid_mass_jump_ = std::max(regrid_mass_jump_, dm);
            regrid_ene_jump_  = std::max(regrid_ene_jump_, de);
            if (cfg_.diag_int > 0 && (dm > 1.0e-12 || de > 1.0e-12)) {
                amrex::Print() << "  regrid at step " << step_
                               << ": rho jump=" << dm << " ene jump=" << de << "\n";
            }
        }
        Real dt = std::min(ComputeDt(), cfg_.t_max - t_);
        AdvanceHierarchy(dt);
        t_ += dt; ++step_;

        if (step_ % cfg_.diag_int == 0) {
            Real divb = 0.0;
            for (int lev = 0; lev <= finest_level; ++lev)
                divb = std::max(divb, MaxDivB(lev));   // MaxDivB уже делает ReduceRealMax
            amrex::Long fb = hlld_fallbacks_;
            ParallelDescriptor::ReduceLongSum(fb);
            amrex::Long fe = floor_events_;
            ParallelDescriptor::ReduceLongSum(fe);
            const amrex::Long nonpos = CountNonPositiveCells();
            amrex::Print() << "step " << step_ << "  t=" << t_
                           << "  dt=" << dt << "  max|divB|=" << divb
                           << "  levels=" << finest_level + 1
                           << "  hlld_fallbacks=" << fb
                           << "  floor_events=" << fe
                           << "  nonpositive_cells=" << nonpos
                           << "  rho_drift=" << (conserved0_set_
                                 ? std::abs(TotalConserved(URHO) - conserved0_[URHO])
                                   / std::abs(conserved0_[URHO]) : Real(0.0))
                           << "\n";
        }
        const bool plot_now =
            (cfg_.plot_int > 0 && step_ % cfg_.plot_int == 0) ||
            (cfg_.plot_dt > 0 && t_ >= next_plot_t - 1e-14);
        if (cfg_.write_plotfiles && plot_now) {
            WritePlotFile(step_, t_);
            last_plot_step = step_;
            if (cfg_.plot_dt > 0) next_plot_t += cfg_.plot_dt;
        }
    }
    if (cfg_.write_plotfiles && last_plot_step != step_) WritePlotFile(step_, t_);
    amrex::Long fb_total = hlld_fallbacks_;
    ParallelDescriptor::ReduceLongSum(fb_total);
    amrex::Long floor_total = floor_events_;
    ParallelDescriptor::ReduceLongSum(floor_total);
    amrex::Print() << "Evolve finished: " << step_ << " steps, t=" << t_
                   << ", hlld_fallbacks=" << fb_total
                   << ", floor_events=" << floor_total
                   << ", nonpositive_cells=" << CountNonPositiveCells() << "\n";

    // Диагностика консервативности иерархии. На полностью периодической задаче
    // потока через внешнюю границу нет, поэтому Σρ и Σe должны сохраняться с
    // точностью округления — но только если газовые потоки согласованы на
    // стыках уровней (reflux). Без рефлюкса дрейф на несколько порядков выше:
    // именно это измеряет tests/check_amr_conservation.py.
    if (conserved0_set_) {
        const auto drift = [&](int n) {
            const Real now = TotalConserved(n);
            const Real ref = conserved0_[n];
            return (std::abs(ref) > 0.0) ? std::abs(now - ref) / std::abs(ref)
                                         : std::abs(now - ref);
        };
        // Суммарный импульс периодического вихря Орзага–Танга равен нулю, и
        // относительная невязка по нему не определена. Нормируем на полную
        // массу: получается ошибка в единицах скорости.
        const auto momentum_drift = [&](int n) {
            return std::abs(TotalConserved(n) - conserved0_[n]) / std::abs(conserved0_[URHO]);
        };
        PrintStateRanges();
        {   // финальная норма div B по всей иерархии (NUM-004)
            Real divb = 0.0, bmax = 0.0;
            for (int lev = 0; lev <= finest_level; ++lev) {
                divb = std::max(divb, MaxDivB(lev));
                for (int d = 0; d < AMREX_SPACEDIM; ++d)
                    bmax = std::max(bmax, bface_[lev][d].norm0());
            }
            ParallelDescriptor::ReduceRealMax(bmax);
            const Real dxmin = std::min(Geom(finest_level).CellSize(0),
                                        Geom(finest_level).CellSize(1));
            amrex::Print() << "divb: max_abs=" << divb
                           << " normalized=" << (bmax > 0.0 ? dxmin * divb / bmax : divb)
                           << "\n";
        }
        amrex::Print() << "regrid_jump: rho=" << regrid_mass_jump_
                       << " ene=" << regrid_ene_jump_ << "\n";
        amrex::Print() << "conservation: levels=" << finest_level + 1
                       << " reflux=" << (cfg_.reflux ? 1 : 0)
                       << " rho_rel_drift=" << drift(URHO)
                       << " mx_vel_drift=" << momentum_drift(UMX)
                       << " my_vel_drift=" << momentum_drift(UMY)
                       << " ene_rel_drift=" << drift(UENE) << "\n";
        // Доля перекрытия важна для интерпретации: если мелкий уровень покрывает
        // весь периодический домен, стыка уровней нет и рефлюкс тождественно
        // ничего не делает — такой прогон не проверяет консервативность.
        for (int lev = 0; lev <= finest_level; ++lev) {
            Long covered = 0;
            if (lev < finest_level) {
                const IntVect& rr = refRatio(lev);
                covered = grids[lev+1].numPts() / (rr[0] * rr[1]);
            }
            amrex::Print() << "  level " << lev << ": cells=" << grids[lev].numPts()
                           << " covered_by_finer=" << covered << "\n";
        }
    }
}

} // namespace mhd
