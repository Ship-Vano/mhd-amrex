#ifndef MHD_PROBLEMS_H_
#define MHD_PROBLEMS_H_
//
// Problems.H — каталог начальных условий.
//
// Для гарантии div B = 0 в начальных данных плоскостные компоненты поля
// задаются через векторный потенциал Az(x,y):
//     Bx(x-грань) =  [Az(x, y+Δy/2) − Az(x, y−Δy/2)] / Δy,
//     By(y-грань) = −[Az(x+Δx/2, y) − Az(x−Δx/2, y)] / Δx.
// Такая разностная форма ∇×A даёт дискретную дивергенцию, равную нулю с
// машинной точностью на staggered-сетке (телескопическое сокращение Az).
// Для кусочно-постоянных полей (Brio–Wu) используется прямое задание Bn.
//
#include <functional>
#include <cmath>
#include <stdexcept>

#include "Config.H"
#include "kernels/MhdState.H"

namespace mhd {

struct Problem {
    // Примитивные газовые величины + Bz в точке (x, y); QBX/QBY игнорируются
    std::function<void(double x, double y, double* q)> prim;
    // Векторный потенциал Az (если direct_b == false)
    std::function<double(double x, double y)> Az;
    // Прямое задание нормальных компонент (если direct_b == true)
    bool direct_b = false;
    std::function<double(double x, double y)> Bx_face, By_face;
};

inline Problem make_problem(const SimConfig& c)
{
    const double pi = M_PI;
    Problem p;

    if (c.problem == "orszag_tang") {
        // Вихрь Орзага–Танга: ρ=25/36π, p=5/12π, γ=5/3, период [0,1]², t_end=0.5
        const double B0 = c.param("B0", 1.0/std::sqrt(4.0*M_PI));
        p.prim = [pi](double x, double y, double* q) {
            q[QRHO] = 25.0 / (36.0 * pi);
            q[QP]   = 5.0  / (12.0 * pi);
            q[QU]   = -std::sin(2*pi*y);
            q[QV]   =  std::sin(2*pi*x);
            q[QW]   = 0.0;  q[QBZ] = 0.0;
        };
        p.Az = [pi, B0](double x, double y) {
            return B0 * (std::cos(4*pi*x)/(4*pi) + std::cos(2*pi*y)/(2*pi));
        };
    }
    else if (c.problem == "brio_wu") {
        // Ударная труба Брио–Ву (2D-полоса): γ=2, Bx=0.75, t_end=0.1
        const double x0  = c.param("x_interface", 0.5);
        const double bxc = c.param("Bx", 0.75);
        p.prim = [x0](double x, double, double* q) {
            if (x < x0) { q[QRHO] = 1.0;   q[QP] = 1.0; }
            else        { q[QRHO] = 0.125; q[QP] = 0.1; }
            q[QU] = q[QV] = q[QW] = 0.0;  q[QBZ] = 0.0;
        };
        p.direct_b = true;
        p.Bx_face = [bxc](double, double) { return bxc; };
        p.By_face = [x0](double x, double) { return (x < x0) ? 1.0 : -1.0; };
    }
    else if (c.problem == "alfven_wave") {
        // Циркулярно поляризованная альфвеновская волна под углом α
        // (по умолчанию α = 30°), область [0, 1/cosα]×[0, 1/sinα], период.
        // B∥=1, B⊥ = A sin φ, Bz = A cos φ, φ = 2π(x cosα + y sinα); v⊥=B⊥, vz=Bz.
        const double alpha = c.param("alpha_deg", 30.0) * pi / 180.0;
        const double A     = c.param("amplitude", 0.1);
        const double ca = std::cos(alpha), sa = std::sin(alpha);
        p.prim = [=](double x, double y, double* q) {
            const double phi = 2*pi*(x*ca + y*sa);
            const double vp = A*std::sin(phi), vz = A*std::cos(phi);
            q[QRHO] = 1.0;  q[QP] = 0.1;
            q[QU] = -vp*sa; q[QV] = vp*ca; q[QW] = vz;
            q[QBZ] = vz;
        };
        p.Az = [=](double x, double y) {
            const double phi = 2*pi*(x*ca + y*sa);
            return ca*y - sa*x + A/(2*pi)*std::cos(phi);   // B∥-фон + волна
        };
    }
    else if (c.problem == "rotor") {
        // Вращающийся цилиндр (Tóth 2000 / ВКРБ §4.1.2)
        const double r0 = c.param("r0", 0.1), r1 = c.param("r1", 0.115);
        const double v0 = c.param("v0", 2.0);
        const double b0 = c.param("B0", 5.0/std::sqrt(4.0*pi));
        p.prim = [=](double x, double y, double* q) {
            const double rx = x-0.5, ry = y-0.5, r = std::sqrt(rx*rx+ry*ry);
            q[QP] = 1.0; q[QW] = 0.0; q[QBZ] = 0.0;
            if (r < r0)      { q[QRHO]=10.0; q[QU]=-v0*ry/r0; q[QV]=v0*rx/r0; }
            else if (r < r1) { const double f=(r1-r)/(r1-r0);
                               q[QRHO]=1.0+9.0*f; q[QU]=-f*v0*ry/r; q[QV]=f*v0*rx/r; }
            else             { q[QRHO]=1.0; q[QU]=q[QV]=0.0; }
        };
        p.Az = [b0](double, double y) { return b0 * y; };   // Bx = b0, By = 0
    }
    else if (c.problem == "magnetic_loop") {
        // Перенос петли магнитного поля (ВКРБ §4.1.4): адвекция слабого поля
        const double A0 = c.param("A0", 1.0e-3), R0 = c.param("R0", 0.3);
        p.prim = [=](double /*x*/, double /*y*/, double* q) {
            q[QRHO]=1.0; q[QP]=1.0; q[QU]=c.param("u0",2.0);
            q[QV]=c.param("v0",1.0); q[QW]=0.0; q[QBZ]=0.0;
        };
        p.Az = [=](double x, double y) {
            const double r = std::sqrt(x*x + y*y);
            return (r <= R0) ? A0*(R0 - r) : 0.0;
        };
    }
    else if (c.problem == "uniform_const") {
        // Постоянное состояние. Схема обязана воспроизводить его точно: любой
        // дрейф означает несогласованность реконструкции, потока или CT, а не
        // физику. Параметры задаются в конфиге, чтобы проверять и случай
        // B = 0, и наклонное поле с ненулевой скоростью.
        const double rho0 = c.param("rho0", 1.0), p0 = c.param("p0", 1.0);
        const double u0 = c.param("u0", 1.0), v0 = c.param("v0", -0.5);
        const double w0 = c.param("w0", 0.25);
        const double bx0 = c.param("Bx", 0.75), by0 = c.param("By", -0.4);
        const double bz0 = c.param("Bz", 0.3);
        p.prim = [=](double, double, double* q) {
            q[QRHO] = rho0; q[QP] = p0;
            q[QU] = u0; q[QV] = v0; q[QW] = w0; q[QBZ] = bz0;
        };
        // Az = Bx*y - By*x  =>  dAz/dy = Bx,  -dAz/dx = By
        p.Az = [bx0, by0](double x, double y) { return bx0 * y - by0 * x; };
    }
    else if (c.problem == "mhd_blast") {
        // Сильный МГД-взрыв в декартовой постановке (Balsara--Spicer 1999;
        // тест Athena). Окружающее β = 2p/|B|^2 = 0.2 -- жёсткая проверка
        // положительности: без отката HLLD->HLL давление уходит в минус.
        // Область [-0.5, 0.5]^2, γ=5/3, t_end = 0.2.
        const double rho0 = c.param("rho0", 1.0);
        const double pin  = c.param("p_in", 10.0);
        const double pout = c.param("p_out", 0.1);
        const double rb   = c.param("r_blast", 0.1);
        const double b0   = c.param("B0", 1.0);
        const double bx0 = b0 / std::sqrt(2.0), by0 = b0 / std::sqrt(2.0);
        p.prim = [=](double x, double y, double* q) {
            q[QRHO] = rho0;
            q[QP]   = (std::sqrt(x*x + y*y) < rb) ? pin : pout;
            q[QU] = q[QV] = q[QW] = 0.0; q[QBZ] = 0.0;
        };
        p.Az = [bx0, by0](double x, double y) { return bx0 * y - by0 * x; };
    }
    else {
        throw std::runtime_error("Unknown problem: " + c.problem);
    }
    return p;
}

} // namespace mhd
#endif
