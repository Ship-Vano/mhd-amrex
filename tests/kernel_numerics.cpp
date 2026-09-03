//
// kernel_numerics.cpp — расширенные регрессии вычислительных ядер AMReX-решателя
// (T04). Проверяются инварианты, ветви HLLD, положительность-сохраняющий откат
// на HLL со счётчиком, ограничители, узловая ЭДС и порядок SSP-RK2.
//
// Сборка через CMake (target mhd2d_kernel_numerics); те же заголовки, что и в
// боевом слое.
//
#include "MhdState.H"
#include "Hlld.H"
#include "Reconstruction.H"
#include "CtUpdate.H"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>

using namespace mhd;

namespace {

int g_failures = 0;

void check(bool ok, const std::string& what)
{
    if (!ok) { std::fprintf(stderr, "FAIL: %s\n", what.c_str()); ++g_failures; }
}

void close(Real a, Real b, Real tol, const std::string& what)
{
    const Real d = std::fabs(a - b);
    const Real s = tol * (Real(1.0) + std::fabs(a) + std::fabs(b));
    if (!(d <= s)) {
        std::fprintf(stderr, "FAIL: %s : |%.17g - %.17g| = %.3g > %.3g\n",
                     what.c_str(), a, b, d, s);
        ++g_failures;
    }
}

bool finite_state(const Real* x, int n)
{
    for (int i = 0; i < n; ++i) if (!std::isfinite(x[i])) return false;
    return true;
}

// Детерминированный ГПСЧ, чтобы регрессия не зависела от платформенного rand().
struct Rng {
    unsigned long long s;
    explicit Rng(unsigned long long seed) : s(seed ? seed : 1) {}
    Real uniform(Real lo, Real hi) {
        s ^= s << 13; s ^= s >> 7; s ^= s << 17;
        const Real u = Real((s >> 11) & ((1ULL << 53) - 1)) / Real(1ULL << 53);
        return lo + (hi - lo) * u;
    }
};

// ---------------------------------------------------------------------------
void test_conversions()
{
    Rng rng(12345);
    const Real gamma = 5.0 / 3.0;
    for (int trial = 0; trial < 2000; ++trial) {
        Real q[NPRIM];
        q[QRHO] = rng.uniform(1e-4, 20.0);
        q[QU]   = rng.uniform(-8.0, 8.0);
        q[QV]   = rng.uniform(-8.0, 8.0);
        q[QW]   = rng.uniform(-8.0, 8.0);
        q[QP]   = rng.uniform(1e-4, 50.0);
        q[QBX]  = rng.uniform(-6.0, 6.0);
        q[QBY]  = rng.uniform(-6.0, 6.0);
        q[QBZ]  = rng.uniform(-6.0, 6.0);
        Real u[NCONS], back[NPRIM];
        prim_to_cons(q, u, gamma);
        cons_to_prim(u, back, gamma);
        for (int n = 0; n < NPRIM; ++n)
            close(back[n], q[n], 1e-12, "prim<->cons round trip");
        // Давление должно ровно один раз вычитать кинетическую и магнитную энергию.
        close(pressure_from_cons(u, gamma), q[QP], 1e-12, "pressure_from_cons");
    }
}

void test_fast_speed()
{
    const Real gamma = 5.0 / 3.0;
    const Real rho = 1.3, p = 0.7;
    // Гидро-предел: B = 0 -> скорость звука.
    close(fast_speed(rho, p, 0.0, 0.0, gamma), std::sqrt(gamma * p / rho),
          1e-13, "fast speed hydro limit");
    // Чисто нормальное поле: быстрая волна расщепляется на звук и альфвен ->
    // cf = max(a, ca_n).
    const Real Bn = 0.9;
    const Real a  = std::sqrt(gamma * p / rho);
    const Real ca = std::sqrt(Bn * Bn / rho);
    close(fast_speed(rho, p, Bn, Bn * Bn, gamma), std::max(a, ca), 1e-13,
          "fast speed normal-field limit = max(a, ca_n)");
    // Поперечное поле: дискриминант максимален, cf^2 = a^2 + cA^2.
    const Real Bt = 1.1;
    close(fast_speed(rho, p, 0.0, Bt * Bt, gamma),
          std::sqrt(gamma * p / rho + Bt * Bt / rho), 1e-13,
          "fast speed transverse-field limit");
    // Монотонность по |B|.
    Real prev = fast_speed(rho, p, 0.2, 0.04, gamma);
    for (Real b2 = 0.1; b2 < 5.0; b2 += 0.1) {
        const Real c = fast_speed(rho, p, 0.2, b2, gamma);
        check(c >= prev - 1e-12, "fast speed monotone in B^2");
        prev = c;
    }
    // cf >= |ca_n| и cf >= a всегда.
    Rng rng(777);
    for (int i = 0; i < 500; ++i) {
        const Real r = rng.uniform(0.1, 5.0), pp = rng.uniform(0.1, 5.0);
        const Real bn = rng.uniform(-3.0, 3.0), bt = rng.uniform(0.0, 3.0);
        const Real b2 = bn * bn + bt * bt;
        const Real cf = fast_speed(r, pp, bn, b2, gamma);
        check(cf >= std::sqrt(gamma * pp / r) - 1e-10, "cf >= a");
        check(cf >= std::fabs(bn) / std::sqrt(r) - 1e-10, "cf >= |ca_n|");
    }
}

void mhd_flux_from_prim(const Real* q, Real Bn, Real* f, Real gamma)
{
    Real qq[NPRIM];
    for (int n = 0; n < NPRIM; ++n) qq[n] = q[n];
    qq[QBX] = Bn;
    mhd_flux_x(qq, f, gamma);
}

void test_hlld_consistency()
{
    Rng rng(2024);
    const Real gamma = 5.0 / 3.0;
    for (int trial = 0; trial < 3000; ++trial) {
        Real q[NPRIM];
        q[QRHO] = rng.uniform(0.05, 8.0);
        q[QU]   = rng.uniform(-3.0, 3.0);
        q[QV]   = rng.uniform(-3.0, 3.0);
        q[QW]   = rng.uniform(-3.0, 3.0);
        q[QP]   = rng.uniform(0.05, 8.0);
        q[QBY]  = rng.uniform(-3.0, 3.0);
        q[QBZ]  = rng.uniform(-3.0, 3.0);
        const Real Bn = rng.uniform(-3.0, 3.0);
        q[QBX] = Bn;

        Real f[NCONS], fexact[NCONS];
        int fb = 0;
        hlld_flux(q, q, Bn, f, gamma, Limits{}, &fb);
        mhd_flux_from_prim(q, Bn, fexact, gamma);
        check(fb == 0, "HLLD fell back on an admissible equal state");
        for (int n = 0; n < NCONS; ++n)
            close(f[n], fexact[n], 1e-10, "HLLD consistency F(q,q)=F_phys");
        check(std::fabs(f[UBX]) < 1e-12, "HLLD transports normal B");

        // hll_flux — тот же инвариант.
        Real fh[NCONS];
        hll_flux(q, q, Bn, fh, gamma);
        for (int n = 0; n < NCONS; ++n)
            close(fh[n], fexact[n], 1e-10, "HLL consistency F(q,q)=F_phys");
    }
}

void test_hlld_supersonic_branches()
{
    const Real gamma = 5.0 / 3.0;
    // Полностью сверхзвуковой поток вправо: SL >= 0 -> f == F(qL).
    {
        Real qL[NPRIM] = {1.0, 12.0, 0.3, -0.1, 1.0, 0.4, 0.2, 0.1};
        Real qR[NPRIM] = {0.8,  11.0, -0.2, 0.0, 0.7, 0.4, -0.3, 0.05};
        Real f[NCONS], fexact[NCONS];
        hlld_flux(qL, qR, 0.4, f, gamma);
        mhd_flux_from_prim(qL, 0.4, fexact, gamma);
        for (int n = 0; n < NCONS; ++n)
            close(f[n], fexact[n], 1e-10, "HLLD right-supersonic == F(qL)");
    }
    // Полностью сверхзвуковой поток влево: SR <= 0 -> f == F(qR).
    {
        Real qL[NPRIM] = {1.0, -12.0, 0.3, -0.1, 1.0, 0.4, 0.2, 0.1};
        Real qR[NPRIM] = {0.8, -11.0, -0.2, 0.0, 0.7, 0.4, -0.3, 0.05};
        Real f[NCONS], fexact[NCONS];
        hlld_flux(qL, qR, 0.4, f, gamma);
        mhd_flux_from_prim(qR, 0.4, fexact, gamma);
        for (int n = 0; n < NCONS; ++n)
            close(f[n], fexact[n], 1e-10, "HLLD left-supersonic == F(qR)");
    }
}

void test_hlld_degeneracies()
{
    const Real gamma = 5.0 / 3.0;
    // Bn = 0 (нет альфвеновских волн) — поток конечен, отката нет.
    {
        Real qL[NPRIM] = {1.0, 0.0, 0.4, 0.0, 1.0, 0.0, 0.6, -0.1};
        Real qR[NPRIM] = {0.9, 0.1, -0.3, 0.2, 0.8, 0.0, -0.5, 0.3};
        Real f[NCONS]; int fb = 0;
        hlld_flux(qL, qR, 0.0, f, gamma, Limits{}, &fb);
        check(finite_state(f, NCONS), "HLLD Bn=0 finite");
        check(fb == 0, "HLLD Bn=0 must not fall back");
        check(std::fabs(f[UBX]) < 1e-12, "HLLD Bn=0 zero normal-B flux");
    }
    // Brio–Wu состояние на границе раздела (gamma = 2).
    {
        Real qL[NPRIM] = {1.0, 0.0, 0.0, 0.0, 1.0, 0.75, 1.0, 0.0};
        Real qR[NPRIM] = {0.125, 0.0, 0.0, 0.0, 0.1, 0.75, -1.0, 0.0};
        Real f[NCONS]; int fb = 0;
        hlld_flux(qL, qR, 0.75, f, 2.0, Limits{}, &fb);
        check(finite_state(f, NCONS), "HLLD Brio-Wu finite");
        check(fb == 0, "HLLD Brio-Wu must not fall back");
        check(std::fabs(f[UBX]) < 1e-12, "HLLD Brio-Wu zero normal-B flux");
    }
    // Вращательный разрыв: сильный сдвиг тангенциальной скорости при сильном Bn.
    {
        Real qL[NPRIM] = {1.0, 0.0, 4.0, 0.0, 1.0, 2.0, 0.0, 0.0};
        Real qR[NPRIM] = {1.0, 0.0, -4.0, 0.0, 1.0, 2.0, 0.0, 0.0};
        Real f[NCONS]; int fb = 0;
        hlld_flux(qL, qR, 2.0, f, gamma, Limits{}, &fb);
        check(finite_state(f, NCONS), "HLLD RD finite");
    }
}

void test_hlld_finiteness_sweep()
{
    Rng rng(99991);
    const Real gamma = 5.0 / 3.0;
    int fb_total = 0;
    for (int trial = 0; trial < 20000; ++trial) {
        Real qL[NPRIM], qR[NPRIM];
        for (int side = 0; side < 2; ++side) {
            Real* q = side ? qR : qL;
            q[QRHO] = std::pow(Real(10.0), rng.uniform(-3.0, 1.3));
            q[QU]   = rng.uniform(-15.0, 15.0);
            q[QV]   = rng.uniform(-15.0, 15.0);
            q[QW]   = rng.uniform(-15.0, 15.0);
            q[QP]   = std::pow(Real(10.0), rng.uniform(-4.0, 2.0));
            q[QBX]  = rng.uniform(-8.0, 8.0);
            q[QBY]  = rng.uniform(-8.0, 8.0);
            q[QBZ]  = rng.uniform(-8.0, 8.0);
        }
        const Real Bn = rng.uniform(-8.0, 8.0);
        Real f[NCONS]; int fb = 0;
        hlld_flux(qL, qR, Bn, f, gamma, Limits{}, &fb);
        check(finite_state(f, NCONS), "HLLD wrapper finite on wide parametric sweep");
        check(std::fabs(f[UBX]) < 1e-10, "HLLD wrapper zero normal-B flux");
        fb_total += fb;
    }
    std::printf("  finiteness sweep: %d/20000 states used the HLL fallback\n", fb_total);
}

void test_positivity_fallback()
{
    const Real gamma = 5.0 / 3.0;
    const Limits lim;

    // Гладкое допустимое состояние — отката быть не должно.
    {
        Real qL[NPRIM] = {1.0, 0.1, 0.0, 0.0, 1.0, 0.5, 0.3, 0.1};
        Real qR[NPRIM] = {0.9, -0.1, 0.05, 0.0, 0.8, 0.5, -0.2, 0.1};
        Real f[NCONS]; int fb = 0;
        hlld_flux(qL, qR, 0.5, f, gamma, lim, &fb);
        check(fb == 0, "smooth admissible state must not fall back");
    }

    // Давление ровно на floor (post-floor путь) — HLLD raw отвергает, откат на HLL.
    {
        Real qL[NPRIM] = {1.0, 0.0, 0.0, 0.0, lim.small_pres, 0.1, 0.1, 0.0};
        Real qR[NPRIM] = {0.5, 0.0, 0.0, 0.0, 0.5, 0.1, -0.1, 0.0};
        Real f[NCONS]; int fb = 0;
        check(!hlld_flux_raw(qL, qR, 0.1, f, gamma, lim), "raw HLLD rejects floored pressure");
        hlld_flux(qL, qR, 0.1, f, gamma, lim, &fb);
        check(fb == 1, "floored-pressure state must expose the HLL fallback");
        check(finite_state(f, NCONS), "HLL fallback flux is finite for floored state");
    }

    // Нулевое/отрицательное входное давление (симптом сбойного предыдущего шага).
    {
        Real qL[NPRIM] = {1.0, 0.0, 0.0, 0.0, 0.0, 0.1, 0.1, 0.0};
        Real qR[NPRIM] = {0.5, 0.0, 0.0, 0.0, 0.5, 0.1, -0.1, 0.0};
        Real f[NCONS]; int fb = 0;
        hlld_flux(qL, qR, 0.1, f, gamma, lim, &fb);
        check(fb == 1 && finite_state(f, NCONS), "p=0 input -> finite HLL fallback");
    }

    // Околовакуумное состояние ниже floor'ов.
    {
        Real qv[NPRIM] = {1e-13, 0.0, 0.0, 0.0, 1e-13, 0.0, 0.0, 0.0};
        Real f[NCONS]; int fb = 0;
        hlld_flux(qv, qv, 0.0, f, gamma, lim, &fb);
        check(fb == 1, "near-vacuum state must expose the HLL fallback");
        check(finite_state(f, NCONS), "near-vacuum HLL fallback flux is finite");
    }
}

void test_limiter()
{
    // Постоянные данные -> нулевой наклон.
    for (auto L : {Limiter::None, Limiter::MinMod, Limiter::MC, Limiter::VanLeer})
        close(limited_slope(0.0, 0.0, L), 0.0, 1e-15, "limiter: constant -> 0 slope");

    // Локальный экстремум (разные знаки) -> нулевой наклон (монотонность).
    for (auto L : {Limiter::MinMod, Limiter::MC, Limiter::VanLeer}) {
        close(limited_slope(1.0, -1.0, L), 0.0, 1e-15, "limiter: extremum -> 0");
        close(limited_slope(-2.0, 3.0, L), 0.0, 1e-15, "limiter: extremum -> 0");
    }

    // Нечётность: slope(-a,-b) = -slope(a,b).
    Rng rng(51);
    for (int i = 0; i < 400; ++i) {
        const Real a = rng.uniform(-4.0, 4.0), b = rng.uniform(-4.0, 4.0);
        for (auto L : {Limiter::MinMod, Limiter::MC, Limiter::VanLeer}) {
            close(limited_slope(-a, -b, L), -limited_slope(a, b, L), 1e-13,
                  "limiter odd symmetry");
        }
    }

    // TVD-границы: |slope| <= 2 min(|dl|,|dr|) для одинаковых знаков.
    for (int i = 0; i < 400; ++i) {
        const Real a = rng.uniform(0.01, 4.0), b = rng.uniform(0.01, 4.0);
        const Real m = 2.0 * std::min(a, b);
        for (auto L : {Limiter::MinMod, Limiter::MC, Limiter::VanLeer}) {
            const Real s = limited_slope(a, b, L);
            check(std::fabs(s) <= m + 1e-12, "limiter TVD bound");
            check(s >= -1e-12, "limiter same-sign slope keeps sign");
        }
        check(std::fabs(limited_slope(a, b, Limiter::None)) < 1e-15, "None limiter -> 0");
    }

    // MinMod на монотонных данных = меньший по модулю наклон.
    close(limited_slope(2.0, 1.0, Limiter::MinMod), 1.0, 1e-15, "minmod picks smaller");
    close(limited_slope(1.0, 3.0, Limiter::MinMod), 1.0, 1e-15, "minmod picks smaller");

    // Реконструкция линейного профиля точна: q(x)=a+b*x, грань i+1/2 в q0+b/2.
    const Real b = 0.7;
    const Real qm = 2.0 - b, q0 = 2.0, qp = 2.0 + b;
    for (auto L : {Limiter::MinMod, Limiter::MC, Limiter::VanLeer}) {
        close(face_value_plus(qm, q0, qp, L), q0 + 0.5 * b, 1e-13, "recon linear exact (plus)");
        close(face_value_minus(qm, q0, qp, L), q0 - 0.5 * b, 1e-13, "recon linear exact (minus)");
    }
    // Постоянный профиль -> грань = ячейка.
    close(face_value_plus(3.0, 3.0, 3.0, Limiter::MC), 3.0, 1e-15, "recon constant (plus)");
    close(face_value_minus(3.0, 3.0, 3.0, Limiter::MC), 3.0, 1e-15, "recon constant (minus)");
}

void test_corner_emf()
{
    // Постоянная ЭДС на всех гранях и ячейках -> та же ЭДС в узле (обе схемы).
    close(corner_emf(2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, EmfAveraging::BalsaraSpicer),
          2.5, 1e-15, "Balsara-Spicer constant state");
    close(corner_emf(2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, EmfAveraging::GardinerStone),
          2.5, 1e-15, "Gardiner-Stone constant state");
    // Balsara–Spicer — арифметическое среднее четырёх граней.
    close(corner_emf(1.0, 2.0, 3.0, 4.0, 9.0, 9.0, 9.0, 9.0, EmfAveraging::BalsaraSpicer),
          2.5, 1e-15, "Balsara-Spicer face average");
    // Gardiner–Stone: 2<faces> - <cells>.
    close(corner_emf(1.0, 2.0, 3.0, 4.0, 1.0, 2.0, 3.0, 4.0, EmfAveraging::GardinerStone),
          2.0 * 2.5 - 2.5, 1e-15, "Gardiner-Stone epsilon=2 form");
    // cell_emf_z = v*Bx - u*By.
    Real q[NPRIM] = {1.0, 0.3, -0.4, 0.0, 1.0, 0.5, 0.2, 0.0};
    close(cell_emf_z(q), q[QV] * q[QBX] - q[QU] * q[QBY], 1e-15, "cell_emf_z formula");

    // --- одномерный предел (Gardiner & Stone 2005, §3.1) ---------------------
    // Требование T04/NEW-006: при исчезающих поперечных градиентах угловая ЭДС
    // обязана в точности сводиться к одномерному значению с грани.
    //
    // Течение зависит только от x. Тогда обе x-грани у узла несут одно и то же
    // upwind-значение E*, а y-грани и ячейки слева/справа несут свои клеточные
    // значения E_L и E_R (по y ничего не меняется, скачка через y-грань нет).
    {
        const Real e_star = 1.7, eL = 0.4, eR = -0.9;
        const Real gs = corner_emf(e_star, e_star, eL, eR, eL, eR, eL, eR,
                                   EmfAveraging::GardinerStone);
        close(gs, e_star, 1e-15, "Gardiner-Stone reduces to the 1-D upwind EMF");

        // Balsara-Spicer в том же пределе даёт 1/2 E* + 1/4 (E_L + E_R) и
        // одномерное значение НЕ воспроизводит. Это не придирка, а причина, по
        // которой схемой по умолчанию выбрана Gardiner-Stone; фиксируем факт,
        // чтобы он не был потерян при рефакторинге.
        const Real bs = corner_emf(e_star, e_star, eL, eR, eL, eR, eL, eR,
                                   EmfAveraging::BalsaraSpicer);
        close(bs, 0.5 * e_star + 0.25 * (eL + eR), 1e-15,
              "Balsara-Spicer 1-D limit is the documented averaged value");
        if (std::fabs(bs - e_star) < 1e-12) {
            std::fprintf(stderr, "Balsara-Spicer unexpectedly reproduces the 1-D "
                                 "upwind EMF; the test no longer distinguishes "
                                 "the two schemes\n");
            std::exit(1);
        }
        // Вырожденный случай E_L = E_R = E*: обе схемы обязаны совпасть.
        close(corner_emf(e_star, e_star, e_star, e_star, e_star, e_star, e_star,
                         e_star, EmfAveraging::BalsaraSpicer),
              e_star, 1e-15, "Balsara-Spicer 1-D limit with no transverse jump");
    }
}

// Одна стадия Эйлера для линейного оператора L(y) = lambda*y.
// SSP-RK2 (Хойн), в точности как в драйверах:
//   y1   = y + dt*L(y)
//   y^n+1= 1/2 ( y + y1 + dt*L(y1) )
double ssprk2_step(double y, double lambda, double dt)
{
    const double y1 = y + dt * lambda * y;
    return 0.5 * (y + y1 + dt * lambda * y1);
}

void test_ssprk2_order()
{
    // Точный шаг: y^{n+1}/y^n = 1 + a + a^2/2  (двучлен Тейлора e^a).
    const double lambda = -1.3;
    for (double dt : {0.1, 0.05, 0.2}) {
        const double a = lambda * dt;
        close(ssprk2_step(1.0, lambda, dt), 1.0 + a + 0.5 * a * a, 1e-14,
              "SSP-RK2 amplification polynomial 1 + a + a^2/2");
    }

    // Наблюдаемый порядок сходимости на y' = lambda y, y(0)=1, T=1.
    const double T = 1.0, exact = std::exp(lambda * T);
    double prev_err = 0.0, prev_h = 0.0;
    double min_rate = 1e9;
    for (int steps : {20, 40, 80, 160, 320}) {
        const double h = T / steps;
        double y = 1.0;
        for (int k = 0; k < steps; ++k) y = ssprk2_step(y, lambda, h);
        const double err = std::fabs(y - exact);
        if (prev_err > 0.0) {
            const double rate = std::log(prev_err / err) / std::log(prev_h / h);
            min_rate = std::min(min_rate, rate);
        }
        prev_err = err; prev_h = h;
    }
    std::printf("  SSP-RK2 observed ODE order (min over refinements): %.3f\n", min_rate);
    check(min_rate > 1.95, "SSP-RK2 observed second-order convergence");
}

} // namespace

int main()
{
    test_conversions();
    test_fast_speed();
    test_hlld_consistency();
    test_hlld_supersonic_branches();
    test_hlld_degeneracies();
    test_hlld_finiteness_sweep();
    test_positivity_fallback();
    test_limiter();
    test_corner_emf();
    test_ssprk2_order();

    if (g_failures) {
        std::fprintf(stderr, "kernel numerics: %d FAILURE(S)\n", g_failures);
        return 1;
    }
    std::puts("kernel numerics: PASS");
    return 0;
}
