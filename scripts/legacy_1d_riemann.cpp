// legacy_1d_riemann.cpp -- first-order 1D HLLD Riemann driver against the
// corrected legacy flux (MHDSolver1D::HLLD_flux_corrected).
//
// Purpose: a like-for-like comparison of the legacy corrected HLLD flux with
// the AMReX kernel on the same 1D Riemann problems from the VKR (Brio-Wu,
// Dai-Woodward).  Piecewise-constant, forward Euler, frozen boundaries --
// the "N0" scheme -- so the difference is the flux implementation and the
// mesh, not the order.
//
// Compiled by scripts/run_legacy_corrected.py --riemann against a fresh
// clone+overlay of the immutable legacy source.  Not a standalone build.
//
//   legacy_1d_riemann <brio_wu|dai_woodward> <Nx> <cfl> <out.csv>

#include "MHDSolver1D.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

namespace {

struct Case {
    double gamma, x_lo, x_hi, t_end, Bx;
    std::vector<double> qL, qR;   // rho,u,v,w,p,By,Bz
};

Case make_case(const std::string& name) {
    const double s = std::sqrt(4.0 * M_PI);
    if (name == "brio_wu")
        return {2.0, -0.5, 0.5, 0.1, 0.75,
                {1.0,   0.0, 0.0, 0.0, 1.0,  1.0, 0.0},
                {0.125, 0.0, 0.0, 0.0, 0.1, -1.0, 0.0}};
    if (name == "dai_woodward")
        return {5.0 / 3.0, -0.5, 0.5, 0.2, 4.0 / s,
                {1.08, 1.2, 0.01, 0.5, 0.95, 3.6 / s, 2.0 / s},
                {1.0,  0.0, 0.0,  0.0, 1.0,  4.0 / s, 2.0 / s}};
    std::fprintf(stderr, "unknown case %s\n", name.c_str());
    std::exit(2);
}

std::vector<double> cons_from(const Case& c, const double* q) {
    // state_from_primitive_vars(rho,u,v,w,p,Bx,By,Bz,gamma)
    return state_from_primitive_vars(q[0], q[1], q[2], q[3], q[4], c.Bx, q[5], q[6], c.gamma);
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 5) {
        std::fprintf(stderr, "usage: %s <brio_wu|dai_woodward> <Nx> <cfl> <out.csv>\n", argv[0]);
        return 2;
    }
    const Case c = make_case(argv[1]);
    const int Nx = std::atoi(argv[2]);
    const double cfl = std::atof(argv[3]);
    const double dx = (c.x_hi - c.x_lo) / Nx;

    std::vector<std::vector<double>> U(Nx), Un(Nx, std::vector<double>(8, 0.0));
    for (int i = 0; i < Nx; ++i) {
        const double x = c.x_lo + (i + 0.5) * dx;
        U[i] = cons_from(c, (x < 0.0 ? c.qL.data() : c.qR.data()));
    }
    const std::vector<double> bcL = U.front(), bcR = U.back();

    std::size_t fallbacks = 0;
    double t = 0.0;
    int step = 0;
    while (t < c.t_end) {
        double smax = 0.0;
        for (const auto& u : U)
            smax = std::max(smax, std::fabs(u[1] / u[0]) + cfast(u, c.gamma));
        double dt = cfl * dx / smax;
        if (t + dt > c.t_end) dt = c.t_end - t;

        std::vector<std::vector<double>> F(Nx + 1);
        for (int f = 0; f <= Nx; ++f) {
            const std::vector<double>& l = (f == 0)      ? bcL : U[f - 1];
            const std::vector<double>& r = (f == Nx)     ? bcR : U[f];
            F[f] = HLLD_flux_corrected(l, r, c.gamma, &fallbacks);
        }
        for (int i = 0; i < Nx; ++i)
            for (int k = 0; k < 8; ++k)
                Un[i][k] = U[i][k] - dt / dx * (F[i + 1][k] - F[i][k]);
        U.swap(Un);
        t += dt;
        ++step;
    }

    FILE* out = std::fopen(argv[4], "w");
    std::fprintf(out, "x,rho,u,v,w,p,Bx,By,Bz\n");
    double rho_min = 1e300, p_min = 1e300;
    for (int i = 0; i < Nx; ++i) {
        const auto& u = U[i];
        const double rho = u[0];
        const double vx = u[1] / rho, vy = u[2] / rho, vz = u[3] / rho;
        const double p = pressure(u, c.gamma);
        rho_min = std::min(rho_min, rho);
        p_min = std::min(p_min, p);
        std::fprintf(out, "%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e\n",
                     c.x_lo + (i + 0.5) * dx, rho, vx, vy, vz, p, u[5], u[6], u[7]);
    }
    std::fclose(out);
    std::printf("[legacy_1d %s] Nx=%d steps=%d t=%.4f rho_min=%.6g p_min=%.6g fallbacks=%zu\n",
                argv[1], Nx, step, t, rho_min, p_min, fallbacks);
    return 0;
}
