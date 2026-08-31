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

    amrex::Vector<amrex::BCRec> bcrec_;     // ГУ клеточных величин для FillPatch
    amrex::Real t_ = 0.0;
    int step_ = 0;

    // Диагностика положительности (инвариант 1 в AGENTS.md): накопительный
    // счётчик срабатываний HLLD→HLL за весь расчёт. На канонических тестах = 0.
    amrex::Long hlld_fallbacks_ = 0;
};

} // namespace mhd
#endif
