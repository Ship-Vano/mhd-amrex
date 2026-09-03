#ifndef MHD_AMR_H_
#define MHD_AMR_H_
//
// MhdAmr.H — AMR-решатель 2D идеальной МГД на базе amrex::AmrCore.
//
// Схема — декартов вариант метода Авдеевой–Лукина (2019):
//   * газовые величины (ρ, ρv, e) и Bz — в центрах ячеек (MultiFab `state_`);
//   * нормальные компоненты Bx, By — на гранях (face-centered MultiFab'ы
//     `bface_`), эволюционируют прямой аппроксимацией закона Фарадея через
//     теорему Стокса с ЭДС из HLLD-потоков → div B = 0 машинно;
//   * интегратор SSP-RK2 (Хойна), CFL ≈ 0.5;
//   * AMR без подциклирования по времени (один Δt на всех уровнях):
//     это позволяет точно согласовать ЭДС на границах уровней — узловые ЭДС
//     мелкого уровня инжектируются в совпадающие узлы грубого до применения
//     закона Фарадея, после чего average_down(_faces) сохраняет div B = 0
//     на всей иерархии (см. REPORT.md, раздел AMR).
//
// Ссылки на используемые средства AMReX:
//   AmrCore / FillPatchTwoLevels / face_divfree_interp / average_down_faces —
//   https://amrex-codes.github.io/amrex/docs_html/AmrCore.html
//
#include <AMReX_AmrCore.H>
#include <AMReX_MultiFab.H>
#include <AMReX_BCRec.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_YAFluxRegister.H>

#include <memory>

#include "Config.H"
#include "Problems.H"

namespace mhd {

class MhdAmr : public amrex::AmrCore
{
public:
    MhdAmr(const SimConfig& cfg);
    ~MhdAmr() override = default;

    void InitData();                 // НУ на всей иерархии
    void Evolve();                   // основной цикл по времени

protected:
    // --- интерфейс AmrCore -------------------------------------------------
    void MakeNewLevelFromScratch(int lev, amrex::Real time,
                                 const amrex::BoxArray& ba,
                                 const amrex::DistributionMapping& dm) override;
    void MakeNewLevelFromCoarse (int lev, amrex::Real time,
                                 const amrex::BoxArray& ba,
                                 const amrex::DistributionMapping& dm) override;
    void RemakeLevel            (int lev, amrex::Real time,
                                 const amrex::BoxArray& ba,
                                 const amrex::DistributionMapping& dm) override;
    void ClearLevel             (int lev) override;
    void ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;

private:
    // --- шаг по времени ------------------------------------------------------
    amrex::Real ComputeDt() const;
    void AdvanceHierarchy(amrex::Real dt);              // один шаг RK2 на всех уровнях
    void EulerStage(amrex::Real dt);                    // одна стадия: потоки→ЭДС→обновление
    void ComputeFluxesAndEmf(int lev);                  // HLLD-потоки и узловые Ez
    void SyncEmfAcrossLevels();                         // инжекция мелких Ez в грубые узлы
    void ApplyUpdates(int lev, amrex::Real dt);         // газ (Годунов) + Фарадей (Стокс)
    void AverageDownAll();                              // согласование уровней

    // --- консервативность на стыке уровней (reflux) -------------------------
    // Без подциклирования шаг Δt одинаков на всех уровнях, поэтому невязка
    // грубой и мелкой аппроксимаций потока на границе уровней снимается одним
    // регистром потоков: грубая ячейка у стыка получает поправку
    //   ΔU = ±Δt/Δx_c · (F_c − ⟨F_f⟩).
    // Плоскостные Bx, By исключены: их ведёт CT, и рефлюкс разрушил бы div B = 0.
    void DefineFluxRegisters();                         // (пере)определение после regrid
    void AccumulateFluxRegisters(amrex::Real dt);       // CrseAdd / FineAdd
    void RefluxAll();                                   // коррекция грубых ячеек
    // Σ по иерархии без учёта перекрытых грубых ячеек: диагностика сохранения.
    amrex::Real TotalConserved(int comp) const;
    void PrintStateRanges() const;   // min/max rho и p по иерархии

    // --- сервис ---------------------------------------------------------------
    void AllocLevel(int lev, const amrex::BoxArray& ba,
                    const amrex::DistributionMapping& dm);
    void InitLevelData(int lev);                        // НУ (через Az — бездивергентно)
    void FillPatchCells(int lev, amrex::MultiFab& mf, amrex::Real time);
    void FillPatchFaces(int lev,
                        amrex::Array<amrex::MultiFab*, AMREX_SPACEDIM> bf,
                        amrex::Real time);
    void FillPhysicalFaceBoundary(int lev);             // ГУ для граневых B
    void SyncCellB(int lev);                            // интерполяция граней в центры (RT0)
    amrex::Real MaxDivB(int lev) const;
    amrex::Long CountNonPositiveCells() const;          // ячейки, требующие floor ρ/p
    void WritePlotFile(int step, amrex::Real time);

    SimConfig cfg_;
    Problem   prob_;

    // данные уровня: [lev]
    amrex::Vector<amrex::MultiFab> state_;                                   // NCONS, 2 ghost
    amrex::Vector<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>> bface_;     // Bx, By
    // рабочие массивы стадии
    amrex::Vector<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>> flux_;      // HLLD-потоки
    amrex::Vector<amrex::MultiFab> emf_;                                     // узловые Ez
    // копии уровня t^n для RK2
    amrex::Vector<amrex::MultiFab> state0_;
    amrex::Vector<amrex::Array<amrex::MultiFab, AMREX_SPACEDIM>> bface0_;

    // flux_reg_[lev] — регистр стыка уровней lev-1 / lev (определён при lev >= 1).
    // Владение через unique_ptr: YAFluxRegister::define() не идемпотентен
    // (дописывает m_cfp_localindex/m_cfp_fab), поэтому после каждой
    // перестройки сетки регистр создаётся заново, а не переопределяется.
    amrex::Vector<std::unique_ptr<amrex::YAFluxRegister>> flux_reg_;
    bool flux_reg_stale_ = true;            // сетки изменились -> переопределить

    amrex::Vector<amrex::BCRec> bcrec_;     // ГУ клеточных величин для FillPatch
    amrex::Real t_ = 0.0;
    int step_ = 0;
    // Полные сохраняющиеся величины при t = 0 (для диагностики дрейфа).
    std::array<amrex::Real, NCONS> conserved0_ {};
    bool conserved0_set_ = false;

    // Диагностика положительности (инвариант 1 в AGENTS.md): накопительный
    // счётчик срабатываний HLLD→HLL за весь расчёт. На канонических тестах = 0.
    amrex::Long hlld_fallbacks_ = 0;
};

} // namespace mhd
#endif
