#ifndef MHD_PROBLEMS_H_
#define MHD_PROBLEMS_H_
//
// Canonical initial and physical boundary conditions.  DeviceProblem is a
// trivially-copyable runtime contract captured by GPU kernels by value: no
// std::function, std::map, or pointer to MhdAmr reaches device code.

#include <cmath>
#include <stdexcept>
#include <type_traits>

#include "Config.H"
#include "kernels/MhdState.H"

namespace mhd {

enum class ProblemKind : int { OrszagTang, BrioWu, AlfvenWave, Rotor,
                               MagneticLoop, UniformConst, MhdBlast };

struct DeviceProblem {
    ProblemKind kind = ProblemKind::OrszagTang;
    // Compact parameters for every current canonical case.  Field names must
    // not leak into device code; their mapping is documented in make_problem.
    Real a = Real(0), b = Real(0), c = Real(0), d = Real(0);
    Real e = Real(0), f = Real(0), g = Real(0), h = Real(0);

    MHD_HD bool direct_b() const noexcept { return kind == ProblemKind::BrioWu; }

    // Primitive gas state plus Bz; QBX/QBY come from staggered face fields.
    MHD_HD void prim(Real x, Real y, Real* q) const noexcept {
        constexpr Real pi = Real(3.141592653589793238462643383279502884);
        switch (kind) {
        case ProblemKind::OrszagTang:
            q[QRHO]=Real(25)/(Real(36)*pi); q[QP]=Real(5)/(Real(12)*pi);
            q[QU]=-std::sin(Real(2)*pi*y); q[QV]=std::sin(Real(2)*pi*x);
            q[QW]=q[QBZ]=Real(0); break;
        case ProblemKind::BrioWu:
            q[QRHO]=(x<a)?Real(1):Real(0.125); q[QP]=(x<a)?Real(1):Real(0.1);
            q[QU]=q[QV]=q[QW]=q[QBZ]=Real(0); break;
        case ProblemKind::AlfvenWave: {
            const Real phi=Real(2)*pi*(x*a+y*b), vp=c*std::sin(phi), vz=c*std::cos(phi);
            q[QRHO]=Real(1); q[QP]=Real(0.1); q[QU]=-vp*b; q[QV]=vp*a; q[QW]=q[QBZ]=vz;
            break;
        }
        case ProblemKind::Rotor: {
            const Real rx=x-Real(0.5), ry=y-Real(0.5), r=std::sqrt(rx*rx+ry*ry);
            q[QP]=Real(1); q[QW]=q[QBZ]=Real(0);
            if (r<a) { q[QRHO]=Real(10); q[QU]=-c*ry/a; q[QV]=c*rx/a; }
            else if (r<b) { const Real w=(b-r)/(b-a); q[QRHO]=Real(1)+Real(9)*w;
                             q[QU]=-w*c*ry/r; q[QV]=w*c*rx/r; }
            else { q[QRHO]=Real(1); q[QU]=q[QV]=Real(0); }
            break;
        }
        case ProblemKind::MagneticLoop:
            q[QRHO]=q[QP]=Real(1); q[QU]=c; q[QV]=d; q[QW]=q[QBZ]=Real(0); break;
        case ProblemKind::UniformConst:
            q[QRHO]=a; q[QP]=b; q[QU]=c; q[QV]=d; q[QW]=e; q[QBZ]=h; break;
        case ProblemKind::MhdBlast:
            q[QRHO]=a; q[QP]=(std::sqrt(x*x+y*y)<d)?b:c; q[QU]=q[QV]=q[QW]=q[QBZ]=Real(0); break;
        }
    }

    // Vector potential for staggered B; Brio--Wu uses direct face fields.
    MHD_HD Real az(Real x, Real y) const noexcept {
        constexpr Real pi = Real(3.141592653589793238462643383279502884);
        switch (kind) {
        case ProblemKind::OrszagTang:
            return a*(std::cos(Real(4)*pi*x)/(Real(4)*pi)+std::cos(Real(2)*pi*y)/(Real(2)*pi));
        case ProblemKind::AlfvenWave:
            return a*y-b*x+c/(Real(2)*pi)*std::cos(Real(2)*pi*(x*a+y*b));
        case ProblemKind::Rotor: return d*y;
        case ProblemKind::MagneticLoop: {
            const Real r=std::sqrt(x*x+y*y); return (r<=b)?a*(b-r):Real(0);
        }
        case ProblemKind::UniformConst: return f*y-g*x;
        case ProblemKind::MhdBlast: return e*y-f*x;
        case ProblemKind::BrioWu: return Real(0);
        }
        return Real(0);
    }

    MHD_HD Real bx_face(Real, Real) const noexcept { return b; }
    MHD_HD Real by_face(Real x, Real) const noexcept { return (x<a)?Real(1):Real(-1); }
};

static_assert(std::is_trivially_copyable<DeviceProblem>::value,
              "DeviceProblem must remain safe to capture by value in GPU lambdas");

inline DeviceProblem make_problem(const SimConfig& c) {
    constexpr Real pi = Real(3.141592653589793238462643383279502884);
    DeviceProblem p;
    if (c.problem == "orszag_tang") {
        p.kind=ProblemKind::OrszagTang; p.a=c.param("B0",Real(1)/std::sqrt(Real(4)*pi));
    } else if (c.problem == "brio_wu") {
        p.kind=ProblemKind::BrioWu; p.a=c.param("x_interface",Real(0.5)); p.b=c.param("Bx",Real(0.75));
    } else if (c.problem == "alfven_wave") {
        p.kind=ProblemKind::AlfvenWave;
        const Real alpha=c.param("alpha_deg",Real(30))*pi/Real(180);
        p.a=std::cos(alpha); p.b=std::sin(alpha); p.c=c.param("amplitude",Real(0.1));
    } else if (c.problem == "rotor") {
        p.kind=ProblemKind::Rotor; p.a=c.param("r0",Real(0.1)); p.b=c.param("r1",Real(0.115));
        p.c=c.param("v0",Real(2)); p.d=c.param("B0",Real(5)/std::sqrt(Real(4)*pi));
    } else if (c.problem == "magnetic_loop") {
        p.kind=ProblemKind::MagneticLoop; p.a=c.param("A0",Real(1.e-3)); p.b=c.param("R0",Real(0.3));
        p.c=c.param("u0",Real(2)); p.d=c.param("v0",Real(1));
    } else if (c.problem == "uniform_const") {
        p.kind=ProblemKind::UniformConst; p.a=c.param("rho0",Real(1)); p.b=c.param("p0",Real(1));
        p.c=c.param("u0",Real(1)); p.d=c.param("v0",Real(-0.5)); p.e=c.param("w0",Real(0.25));
        p.f=c.param("Bx",Real(0.75)); p.g=c.param("By",Real(-0.4)); p.h=c.param("Bz",Real(0.3));
    } else if (c.problem == "mhd_blast") {
        p.kind=ProblemKind::MhdBlast; p.a=c.param("rho0",Real(1)); p.b=c.param("p_in",Real(10));
        p.c=c.param("p_out",Real(0.1)); p.d=c.param("r_blast",Real(0.1));
        const Real b0=c.param("B0",Real(1)); p.e=b0/std::sqrt(Real(2)); p.f=p.e;
    } else {
        throw std::runtime_error("Unknown problem: " + c.problem);
    }
    return p;
}

} // namespace mhd
#endif
