#ifndef MHD_HLLD_H_
#define MHD_HLLD_H_
//
// Hlld.H — приближённый решатель задачи Римана HLLD для идеальной МГД
// (Miyoshi & Kusano, JCP 208 (2005) 315–344) — тот же поток, что используется
// в схеме Авдеевой–Лукина (2019) и в ВКРБ Шаманова И.П.
//
// Решатель работает в локальной системе координат грани (x' — нормаль,
// y' — касательная): на структурированной сетке это сводится к перестановке
// компонент векторов. Помимо потока газовых величин решатель возвращает
// ЭДС  Ez' = -(v×B)·ez' на грани — она равна компоненте потока для By'
// со знаком минус и используется CT-обновлением магнитного поля.
//
// Устойчивость. Промежуточные состояния HLLD вырождаются, когда
// (S-u)(S-SM) → Bn²/ρ или когда веер Римана близок к вакууму. В таких случаях
// HLLD может дать неконечный или неположительный поток. hlld_flux() проверяет
// результат предикатом допустимости и при отказе переключается на
// положительность-сохраняющий поток HLL (HLLE) с теми же оценками крайних
// скоростей. Счётчик срабатываний (n_fallback) пишется в диагностику: на
// гладких канонических тестах он обязан оставаться нулём.
//
#include "MhdState.H"

namespace mhd {

// Поток идеальной МГД в направлении x' от примитивного состояния.
// ВАЖНО: поток для нормальной компоненты Bn тождественно равен нулю
// (Bn непрерывно на грани), поток для Bt равен -Ez = u*Bt - v*Bn.
MHD_HD void mhd_flux_x(const Real* q, Real* f, Real gamma)
{
    const Real rho = q[QRHO];
    const Real u = q[QU], v = q[QV], w = q[QW];
    const Real Bx = q[QBX], By = q[QBY], Bz = q[QBZ];
    const Real p  = q[QP];
    const Real B2 = Bx*Bx + By*By + Bz*Bz;
    const Real pT = p + Real(0.5) * B2;
    const Real e  = p / (gamma - Real(1.0)) + Real(0.5)*rho*(u*u+v*v+w*w) + Real(0.5)*B2;
    const Real vB = u*Bx + v*By + w*Bz;

    f[URHO] = rho * u;
    f[UMX]  = rho * u * u + pT - Bx * Bx;
    f[UMY]  = rho * u * v - Bx * By;
    f[UMZ]  = rho * u * w - Bx * Bz;
    f[UENE] = (e + pT) * u - Bx * vB;
    f[UBX]  = Real(0.0);
    f[UBY]  = u * By - v * Bx;   // = -Ez
    f[UBZ]  = u * Bz - w * Bx;   // = +Ey
}

// Оценка крайних скоростей веера (Davis-type, шире эйнфельдтовой →
// гарантирует SL ≤ λ_min, SR ≥ λ_max, значит HLL-поток положителен при CFL).
MHD_HD void mhd_wave_speeds(const Real* qL, const Real* qR, Real Bn, Real gamma,
                            Real& SL, Real& SR)
{
    const Real B2L = Bn*Bn + qL[QBY]*qL[QBY] + qL[QBZ]*qL[QBZ];
    const Real B2R = Bn*Bn + qR[QBY]*qR[QBY] + qR[QBZ]*qR[QBZ];
    const Real cfL = fast_speed(qL[QRHO], qL[QP], Bn, B2L, gamma);
    const Real cfR = fast_speed(qR[QRHO], qR[QP], Bn, B2R, gamma);
    const Real cmax = (cfL > cfR) ? cfL : cfR;
    SL = ((qL[QU] < qR[QU]) ? qL[QU] : qR[QU]) - cmax;
    SR = ((qL[QU] > qR[QU]) ? qL[QU] : qR[QU]) + cmax;
}

// --------------------------------------------------------------------------
// HLL (HLLE) поток — один усреднённый состояние в веере. Положительность
// плотности/давления сохраняется при SL ≤ λ_min, SR ≥ λ_max (Einfeldt et al.).
// Используется как безопасный откат HLLD и как режим аблации N0.
// --------------------------------------------------------------------------
MHD_HD void hll_flux(const Real* qL_in, const Real* qR_in, Real Bn,
                     Real* f, Real gamma)
{
    Real qL[NPRIM], qR[NPRIM];
    for (int n = 0; n < NPRIM; ++n) { qL[n] = qL_in[n]; qR[n] = qR_in[n]; }
    qL[QBX] = Bn; qR[QBX] = Bn;

    Real SL, SR;
    mhd_wave_speeds(qL, qR, Bn, gamma, SL, SR);

    Real UL[NCONS], UR[NCONS], FL[NCONS], FR[NCONS];
    prim_to_cons(qL, UL, gamma);
    prim_to_cons(qR, UR, gamma);
    mhd_flux_x(qL, FL, gamma);
    mhd_flux_x(qR, FR, gamma);

    if (SL >= Real(0.0)) { for (int n = 0; n < NCONS; ++n) f[n] = FL[n]; return; }
    if (SR <= Real(0.0)) { for (int n = 0; n < NCONS; ++n) f[n] = FR[n]; return; }

    const Real inv = Real(1.0) / (SR - SL);
    for (int n = 0; n < NCONS; ++n)
        f[n] = (SR * FL[n] - SL * FR[n] + SL * SR * (UR[n] - UL[n])) * inv;
    f[UBX] = Real(0.0);   // Bn непрерывно — поток нормальной компоненты нулевой
}

// ---------------------------------------------------------------------------
// HLLD, «сырой» вариант: вход — примитивные состояния слева/справа в локальных
// осях грани, нормальная компонента Bn (одна на грань). Возвращает true, если
// поток допустим (конечен и веер не вырожден/не вакуумный), иначе false —
// вызывающий hlld_flux() тогда откатывается на hll_flux().
// f[UBY] = -Ez на грани.
// ---------------------------------------------------------------------------
MHD_HD bool hlld_flux_raw(const Real* qL_in, const Real* qR_in, Real Bn,
                          Real* f, Real gamma, const Limits& lim = Limits{})
{
    // Локальные копии: нормальная компонента поля общая для двух состояний.
    Real qL[NPRIM], qR[NPRIM];
    for (int n = 0; n < NPRIM; ++n) { qL[n] = qL_in[n]; qR[n] = qR_in[n]; }
    qL[QBX] = Bn; qR[QBX] = Bn;

    const Real rL = qL[QRHO], rR = qR[QRHO];
    const Real uL = qL[QU],  uR = qR[QU];
    const Real pL = qL[QP],  pR = qR[QP];
    if (!(rL > lim.small_rho) || !(rR > lim.small_rho) ||
        !(pL > lim.small_pres) || !(pR > lim.small_pres)) return false;
    const Real B2L = Bn*Bn + qL[QBY]*qL[QBY] + qL[QBZ]*qL[QBZ];
    const Real B2R = Bn*Bn + qR[QBY]*qR[QBY] + qR[QBZ]*qR[QBZ];
    const Real pTL = pL + Real(0.5) * B2L;
    const Real pTR = pR + Real(0.5) * B2R;

    Real SL, SR;
    mhd_wave_speeds(qL, qR, Bn, gamma, SL, SR);

    Real UL[NCONS], UR[NCONS], FL[NCONS], FR[NCONS];
    prim_to_cons(qL, UL, gamma);
    prim_to_cons(qR, UR, gamma);
    mhd_flux_x(qL, FL, gamma);
    mhd_flux_x(qR, FR, gamma);

    if (SL >= Real(0.0)) { for (int n = 0; n < NCONS; ++n) f[n] = FL[n]; return true; }
    if (SR <= Real(0.0)) { for (int n = 0; n < NCONS; ++n) f[n] = FR[n]; return true; }

    // Положительность веера: плотность HLL-осреднённого состояния должна быть
    // строго положительна (Einfeldt et al.). Иначе звёздные состояния HLLD
    // ненадёжны — откат на HLL.
    const Real inv_s = Real(1.0) / (SR - SL);
    const Real rho_hll = (SR*UR[URHO] - SL*UL[URHO] - (FR[URHO] - FL[URHO])) * inv_s;
    if (!(rho_hll > lim.small_rho)) return false;

    // Скорость контактного разрыва SM (формула (38) M&K05)
    const Real num = (SR - uR)*rR*uR - (SL - uL)*rL*uL - pTR + pTL;
    const Real den = (SR - uR)*rR - (SL - uL)*rL;   // > 0 структурно (SL<uL, SR>uR)
    const Real SM  = num / den;

    // Полное давление в звёздных областях (формула (41))
    const Real pTs = ((SR - uR)*rR*pTL - (SL - uL)*rL*pTR
                      + rL*rR*(SR - uR)*(SL - uL)*(uR - uL)) / den;

    const Real Bn2 = Bn * Bn;

    // --- Внешние звёздные состояния U*L, U*R (формулы (43)–(48)) ----------
    auto star_state = [&] (const Real* q, const Real* U, Real S, Real pT,
                           Real* Us, Real& sqrt_rs) {
        const Real r  = q[QRHO];
        const Real u  = q[QU];
        const Real sv = S - u;     // S - u
        const Real sm = S - SM;    // S - SM
        const Real rs = r * sv / sm;
        const Real d  = r * sv * sm - Bn2;

        Real v_s, w_s, By_s, Bz_s;
        if (std::fabs(d) < Real(1.0e-12) * r * sv * sv + Real(1.0e-30)) {
            // Вырожденный случай (S-SM)(S-u) ≈ Bn²/ρ: тангенциальные величины
            // не меняются (см. обсуждение после (51) в M&K05)
            v_s = q[QV]; w_s = q[QW]; By_s = q[QBY]; Bz_s = q[QBZ];
        } else {
            const Real c1 = (SM - u) / d;
            const Real c2 = (r * sv * sv - Bn2) / d;
            v_s  = q[QV] - Bn * q[QBY] * c1;
            w_s  = q[QW] - Bn * q[QBZ] * c1;
            By_s = q[QBY] * c2;
            Bz_s = q[QBZ] * c2;
        }
        const Real vBs = SM*Bn + v_s*By_s + w_s*Bz_s;
        const Real vB  = u*Bn + q[QV]*q[QBY] + q[QW]*q[QBZ];
        const Real e   = U[UENE];
        const Real es  = (sv * e - (q[QP] + Real(0.5)*(Bn2 + q[QBY]*q[QBY] + q[QBZ]*q[QBZ])) * u
                          + pT * SM + Bn * (vB - vBs)) / sm;

        Us[URHO] = rs;
        Us[UMX]  = rs * SM;
        Us[UMY]  = rs * v_s;
        Us[UMZ]  = rs * w_s;
        Us[UENE] = es;
        Us[UBX]  = Bn;
        Us[UBY]  = By_s;
        Us[UBZ]  = Bz_s;
        sqrt_rs  = std::sqrt(rs > lim.small_rho ? rs : lim.small_rho);
    };

    Real UsL[NCONS], UsR[NCONS], srL, srR;
    // ВАЖНО: в энергию звёздных состояний входит ОБЩЕЕ звёздное давление pTs
    // (формула (48) Miyoshi & Kusano 2005), а не pTL/pTR.
    star_state(qL, UL, SL, pTs, UsL, srL);
    star_state(qR, UR, SR, pTs, UsR, srR);

    if (!(UsL[URHO] > lim.small_rho) || !(UsR[URHO] > lim.small_rho)) return false;

    // Альфвеновские скорости в звёздных областях (формула (51))
    const Real SsL = SM - std::fabs(Bn) / srL;
    const Real SsR = SM + std::fabs(Bn) / srR;

    if (SsL >= Real(0.0)) {            // область U*L
        for (int n = 0; n < NCONS; ++n) f[n] = FL[n] + SL * (UsL[n] - UL[n]);
    } else if (SsR <= Real(0.0)) {     // область U*R
        for (int n = 0; n < NCONS; ++n) f[n] = FR[n] + SR * (UsR[n] - UR[n]);
    } else if (Bn2 > Real(1.0e-30)) {
        // --- Внутренние состояния U**L, U**R (формулы (59)–(63)) ----------
        Real Uss[NCONS];
        const Real sgn = (Bn > Real(0.0)) ? Real(1.0) : Real(-1.0);
        const Real inv = Real(1.0) / (srL + srR);
        const Real vL = UsL[UMY]/UsL[URHO], wL = UsL[UMZ]/UsL[URHO];
        const Real vR = UsR[UMY]/UsR[URHO], wR = UsR[UMZ]/UsR[URHO];

        const Real vss  = (srL*vL + srR*vR + (UsR[UBY] - UsL[UBY]) * sgn) * inv;
        const Real wss  = (srL*wL + srR*wR + (UsR[UBZ] - UsL[UBZ]) * sgn) * inv;
        const Real Byss = (srL*UsR[UBY] + srR*UsL[UBY] + srL*srR*(vR - vL) * sgn) * inv;
        const Real Bzss = (srL*UsR[UBZ] + srR*UsL[UBZ] + srL*srR*(wR - wL) * sgn) * inv;

        const bool left = (SM >= Real(0.0));
        const Real* Us  = left ? UsL : UsR;
        const Real  sg2 = left ? -sgn : sgn;   // знак при скачке энергии (формула (63))
        const Real  rs  = Us[URHO];
        const Real vBs  = SM*Bn + (Us[UMY]/rs)*Us[UBY] + (Us[UMZ]/rs)*Us[UBZ];
        const Real vBss = SM*Bn + vss*Byss + wss*Bzss;
        const Real srs  = std::sqrt(rs);

        Uss[URHO] = rs;
        Uss[UMX]  = rs * SM;
        Uss[UMY]  = rs * vss;
        Uss[UMZ]  = rs * wss;
        Uss[UENE] = Us[UENE] + sg2 * srs * (vBs - vBss);
        Uss[UBX]  = Bn;
        Uss[UBY]  = Byss;
        Uss[UBZ]  = Bzss;

        if (left) {
            for (int n = 0; n < NCONS; ++n)
                f[n] = FL[n] + SsL*Uss[n] - (SsL - SL)*UsL[n] - SL*UL[n];
        } else {
            for (int n = 0; n < NCONS; ++n)
                f[n] = FR[n] + SsR*Uss[n] - (SsR - SR)*UsR[n] - SR*UR[n];
        }
    } else {
        // Bn = 0: альфвеновские волны вырождаются, HLLD → HLLC
        const Real* Us = (SM >= Real(0.0)) ? UsL : UsR;
        if (SM >= Real(0.0)) {
            for (int n = 0; n < NCONS; ++n) f[n] = FL[n] + SL * (Us[n] - UL[n]);
        } else {
            for (int n = 0; n < NCONS; ++n) f[n] = FR[n] + SR * (Us[n] - UR[n]);
        }
    }

    for (int n = 0; n < NCONS; ++n)
        if (!std::isfinite(f[n])) return false;
    return true;
}

// ---------------------------------------------------------------------------
// HLLD с положительность-сохраняющим откатом. n_fallback (если не nullptr)
// инкрементируется при каждом переключении на HLL. На GPU счётчик нужно
// заменить на atomic-редукцию.
// ---------------------------------------------------------------------------
MHD_HD void hlld_flux(const Real* qL, const Real* qR, Real Bn,
                      Real* f, Real gamma, const Limits& lim = Limits{},
                      int* n_fallback = nullptr)
{
    if (hlld_flux_raw(qL, qR, Bn, f, gamma, lim)) return;
    hll_flux(qL, qR, Bn, f, gamma);
    if (n_fallback) *n_fallback += 1;
}

} // namespace mhd
#endif
