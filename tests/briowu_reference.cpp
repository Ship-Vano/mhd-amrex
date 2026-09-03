// briowu_reference.cpp — независимый эталон для задачи Брио–Ву (фаза T06).
//
// ЗАЧЕМ ОТДЕЛЬНАЯ РЕАЛИЗАЦИЯ. Эталон, посчитанный тем же решателем, что и
// проверяемый результат, доказывает только самосогласованность. Точный
// 7-волновой римановский решатель здесь неприменим: решение Брио–Ву содержит
// СОСТАВНУЮ волну (медленная ударная + волна разрежения), а классические точные
// решатели строятся в предположении регулярных волн (D-003).
//
// Поэтому эталон строится схемой, не имеющей с проверяемой ни одной общей
// детали: центральная схема Куртганова–Тадмора (KT). Она не раскладывает
// решение по волнам, не использует HLLD и вообще не решает задачу Римана —
// нужен только физический поток f(u) и оценка максимальной скорости. Файл
// намеренно самодостаточен: он НЕ включает src/kernels, чтобы ошибка в общем
// коде не воспроизвелась одинаково в обеих ветвях.
//
//   briowu_reference <Nx> <t_end> <out.csv> [cfl]
//
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <string>
#include <vector>

namespace {

constexpr int NV = 7;                 // rho, rho u, rho v, rho w, E, By, Bz
constexpr int RHO = 0, MX = 1, MY = 2, MZ = 3, EN = 4, BY = 5, BZ = 6;

double g_gamma = 2.0;
double g_bx    = 0.75;                // нормальная компонента постоянна в 1D

struct Prim { double rho, u, v, w, p, by, bz; };

Prim to_prim(const double* q)
{
    Prim s;
    s.rho = q[RHO];
    s.u = q[MX] / s.rho;  s.v = q[MY] / s.rho;  s.w = q[MZ] / s.rho;
    s.by = q[BY];  s.bz = q[BZ];
    const double kin = 0.5 * s.rho * (s.u*s.u + s.v*s.v + s.w*s.w);
    const double mag = 0.5 * (g_bx*g_bx + s.by*s.by + s.bz*s.bz);
    s.p = (g_gamma - 1.0) * (q[EN] - kin - mag);
    return s;
}

// Физический поток идеальной МГД вдоль x при постоянном Bx.
void phys_flux(const double* q, double* f)
{
    const Prim s = to_prim(q);
    const double b2  = g_bx*g_bx + s.by*s.by + s.bz*s.bz;
    const double ptot = s.p + 0.5 * b2;
    const double bdotv = g_bx*s.u + s.by*s.v + s.bz*s.w;
    f[RHO] = q[MX];
    f[MX]  = q[MX] * s.u + ptot - g_bx * g_bx;
    f[MY]  = q[MY] * s.u - g_bx * s.by;
    f[MZ]  = q[MZ] * s.u - g_bx * s.bz;
    f[EN]  = (q[EN] + ptot) * s.u - g_bx * bdotv;
    f[BY]  = s.by * s.u - g_bx * s.v;
    f[BZ]  = s.bz * s.u - g_bx * s.w;
}

double max_speed(const double* q)
{
    const Prim s = to_prim(q);
    const double a2  = g_gamma * std::max(s.p, 1e-14) / s.rho;
    const double b2  = (g_bx*g_bx + s.by*s.by + s.bz*s.bz) / s.rho;
    const double bn2 = g_bx * g_bx / s.rho;
    const double d   = std::max(0.0, (a2 + b2) * (a2 + b2) - 4.0 * a2 * bn2);
    const double cf  = std::sqrt(0.5 * (a2 + b2 + std::sqrt(d)));
    return std::abs(s.u) + cf;
}

double minmod3(double a, double b, double c)
{
    if (a > 0 && b > 0 && c > 0) return std::min(a, std::min(b, c));
    if (a < 0 && b < 0 && c < 0) return std::max(a, std::max(b, c));
    return 0.0;
}

}  // namespace

int main(int argc, char** argv)
{
    if (argc < 4) {
        std::fprintf(stderr, "usage: briowu_reference <Nx> <t_end> <out.csv> [cfl]\n");
        return 2;
    }
    const int    nx    = std::atoi(argv[1]);
    const double t_end = std::atof(argv[2]);
    const std::string out = argv[3];
    const double cfl   = (argc > 4) ? std::atof(argv[4]) : 0.4;

    const int ng = 3;
    const int n  = nx + 2 * ng;
    const double dx = 1.0 / nx;
    std::vector<double> u(n * NV), u0(n * NV), rhs(n * NV);

    // Начальные данные Брио–Ву: gamma = 2, Bx = 0.75, разрыв в x = 0.5.
    auto set_state = [&](double* q, double rho, double p, double by) {
        q[RHO] = rho; q[MX] = q[MY] = q[MZ] = 0.0;
        q[BY] = by;   q[BZ] = 0.0;
        q[EN] = p / (g_gamma - 1.0) + 0.5 * (g_bx*g_bx + by*by);
    };
    for (int i = 0; i < n; ++i) {
        const double x = (i - ng + 0.5) * dx;
        if (x < 0.5) set_state(&u[i*NV], 1.0,   1.0,  1.0);
        else         set_state(&u[i*NV], 0.125, 0.1, -1.0);
    }

    auto fill_ghost = [&](std::vector<double>& s) {
        for (int g = 0; g < ng; ++g)
            for (int k = 0; k < NV; ++k) {
                s[g*NV + k]              = s[ng*NV + k];              // замороженные
                s[(n-1-g)*NV + k]        = s[(n-ng-1)*NV + k];
            }
    };

    // Полудискретный оператор KT: реконструкция minmod, поток
    //   H = ½(f(u⁺)+f(u⁻)) − ½a(u⁺−u⁻),  a = max(λ⁺, λ⁻).
    auto compute_rhs = [&](std::vector<double>& s, std::vector<double>& out_rhs) {
        fill_ghost(s);
        std::vector<double> uL(n * NV), uR(n * NV);
        for (int i = 1; i < n - 1; ++i)
            for (int k = 0; k < NV; ++k) {
                const double dm = s[i*NV+k]     - s[(i-1)*NV+k];
                const double dp = s[(i+1)*NV+k] - s[i*NV+k];
                const double sl = minmod3(2.0*dm, 0.5*(dm+dp), 2.0*dp);
                uL[i*NV+k] = s[i*NV+k] - 0.5 * sl;   // значение на левой грани ячейки
                uR[i*NV+k] = s[i*NV+k] + 0.5 * sl;   // на правой
            }
        std::vector<double> flux(n * NV, 0.0);
        for (int i = 1; i < n - 2; ++i) {            // грань между i и i+1
            const double* qm = &uR[i*NV];            // слева от грани
            const double* qp = &uL[(i+1)*NV];        // справа
            double fm[NV], fp[NV];
            phys_flux(qm, fm); phys_flux(qp, fp);
            const double a = std::max(max_speed(qm), max_speed(qp));
            for (int k = 0; k < NV; ++k)
                flux[i*NV+k] = 0.5 * (fm[k] + fp[k]) - 0.5 * a * (qp[k] - qm[k]);
        }
        for (int i = 2; i < n - 2; ++i)
            for (int k = 0; k < NV; ++k)
                out_rhs[i*NV+k] = -(flux[i*NV+k] - flux[(i-1)*NV+k]) / dx;
    };

    double t = 0.0;
    long steps = 0;
    while (t < t_end) {
        fill_ghost(u);
        double amax = 0.0;
        for (int i = ng; i < n - ng; ++i) amax = std::max(amax, max_speed(&u[i*NV]));
        double dt = cfl * dx / amax;
        if (t + dt > t_end) dt = t_end - t;

        u0 = u;
        compute_rhs(u, rhs);                                   // SSP-RK2, стадия 1
        for (int i = 0; i < n * NV; ++i) u[i] = u0[i] + dt * rhs[i];
        compute_rhs(u, rhs);                                   // стадия 2
        for (int i = 0; i < n * NV; ++i)
            u[i] = 0.5 * (u0[i] + u[i] + dt * rhs[i]);
        t += dt; ++steps;
    }

    double rho_min = 1e300, p_min = 1e300;
    FILE* f = std::fopen(out.c_str(), "w");
    if (!f) { std::fprintf(stderr, "cannot open %s\n", out.c_str()); return 1; }
    std::fprintf(f, "x,rho,u,v,w,p,By,Bz\n");
    for (int i = ng; i < n - ng; ++i) {
        const Prim s = to_prim(&u[i*NV]);
        const double x = (i - ng + 0.5) * dx;
        rho_min = std::min(rho_min, s.rho);
        p_min   = std::min(p_min, s.p);
        std::fprintf(f, "%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g,%.10g\n",
                     x, s.rho, s.u, s.v, s.w, s.p, s.by, s.bz);
    }
    std::fclose(f);
    std::printf("[briowu_ref KT] Nx=%d steps=%ld t=%.6f rho_min=%.6g p_min=%.6g -> %s\n",
                nx, steps, t, rho_min, p_min, out.c_str());
    return 0;
}
