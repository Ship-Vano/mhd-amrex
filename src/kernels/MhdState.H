#ifndef MHD_STATE_H_
#define MHD_STATE_H_
//
// MhdState.H — базовые определения переменных идеальной МГД и уравнение состояния.
//
// Ядра решателя (HLLD, реконструкция, CT-обновление) написаны без зависимостей
// от AMReX: при сборке внутри AMReX макрос MHD_HD разворачивается в
// AMREX_GPU_HOST_DEVICE (ядра пригодны для GPU-запуска через amrex::ParallelFor),
// при автономной сборке (verification driver) — в пустую строку.
//
#ifdef AMREX_USE_GPU
#include <AMReX_GpuQualifiers.H>
#endif
#ifdef MHD_WITH_AMREX
#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>
#define MHD_HD AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
namespace mhd { using Real = amrex::Real; }
#else
#define MHD_HD inline
namespace mhd { using Real = double; }
#endif

#include <cmath>

namespace mhd {

// ---------------------------------------------------------------------------
// Нумерация консервативных переменных (cell-centered MultiFab).
// Магнитное поле хранится ДВАЖДЫ:
//   * нормальные компоненты Bx, By - на гранях (отдельные face-centered массивы,
//     первичные данные, эволюционируют по закону Фарадея => div B = 0);
//   * клеточные значения UBX..UBZ - интерполяция граней в центры
//     (аналог восстановления Равьяра–Тома из статьи Авдеевой–Лукина;
//     на прямоугольной ячейке RT0-интерполяция даёт среднее двух граней).
// Клеточные B нужны решателю Римана и выводу результатов.
// ---------------------------------------------------------------------------
enum ConsVar : int {
    URHO = 0,   // плотность ρ
    UMX,        // импульс ρu
    UMY,        // импульс ρv
    UMZ,        // импульс ρw   (2.5D: всё зависит только от x,y)
    UENE,       // полная энергия e = p/(γ-1) + ρ|v|²/2 + |B|²/2
    UBX,        // Bx (клеточная копия)
    UBY,        // By (клеточная копия)
    UBZ,        // Bz (чисто клеточная переменная, эволюционирует консервативно)
    NCONS
};

// Примитивные переменные для реконструкции
enum PrimVar : int {
    QRHO = 0, QU, QV, QW, QP, QBX, QBY, QBZ, NPRIM
};

// Малые ограничители для положительности
struct Limits {
    Real small_rho  = 1.0e-12;
    Real small_pres = 1.0e-12;
};

// --------- Уравнение состояния идеального газа ----------------------------
MHD_HD Real pressure_from_cons(const Real* u, Real gamma)
{
    const Real rho = u[URHO];
    const Real ke  = Real(0.5) * (u[UMX]*u[UMX] + u[UMY]*u[UMY] + u[UMZ]*u[UMZ]) / rho;
    const Real me  = Real(0.5) * (u[UBX]*u[UBX] + u[UBY]*u[UBY] + u[UBZ]*u[UBZ]);
    return (gamma - Real(1.0)) * (u[UENE] - ke - me);
}

MHD_HD void cons_to_prim(const Real* u, Real* q, Real gamma, const Limits& lim = Limits{})
{
    const Real rho = (u[URHO] > lim.small_rho) ? u[URHO] : lim.small_rho;
    q[QRHO] = rho;
    q[QU]   = u[UMX] / rho;
    q[QV]   = u[UMY] / rho;
    q[QW]   = u[UMZ] / rho;
    q[QBX]  = u[UBX];
    q[QBY]  = u[UBY];
    q[QBZ]  = u[UBZ];
    Real p  = pressure_from_cons(u, gamma);
    q[QP]   = (p > lim.small_pres) ? p : lim.small_pres;
}

MHD_HD void prim_to_cons(const Real* q, Real* u, Real gamma)
{
    u[URHO] = q[QRHO];
    u[UMX]  = q[QRHO] * q[QU];
    u[UMY]  = q[QRHO] * q[QV];
    u[UMZ]  = q[QRHO] * q[QW];
    u[UBX]  = q[QBX];
    u[UBY]  = q[QBY];
    u[UBZ]  = q[QBZ];
    u[UENE] = q[QP] / (gamma - Real(1.0))
            + Real(0.5) * q[QRHO] * (q[QU]*q[QU] + q[QV]*q[QV] + q[QW]*q[QW])
            + Real(0.5) * (q[QBX]*q[QBX] + q[QBY]*q[QBY] + q[QBZ]*q[QBZ]);
}

// Быстрая магнитозвуковая скорость в направлении n (Bn — нормальная компонента B)
MHD_HD Real fast_speed(Real rho, Real p, Real Bn, Real B2, Real gamma)
{
    const Real a2  = gamma * p / rho;          // квадрат скорости звука
    const Real ca2 = B2 / rho;                 // |B|²/ρ
    const Real cn2 = Bn * Bn / rho;            // альфвеновская по нормали
    const Real s   = a2 + ca2;
    Real disc = s * s - Real(4.0) * a2 * cn2;
    disc = (disc > Real(0.0)) ? std::sqrt(disc) : Real(0.0);
    return std::sqrt(Real(0.5) * (s + disc));
}

} // namespace mhd
#endif
