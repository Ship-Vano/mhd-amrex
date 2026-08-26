#include "CtUpdate.H"
#include "MhdState.H"
#include "Reconstruction.H"

#include <cmath>
#include <cstdio>

namespace {

bool close(mhd::Real actual, mhd::Real expected, mhd::Real tolerance = 1.0e-13)
{
    return std::fabs(actual - expected) <= tolerance;
}

int fail(const char* name, mhd::Real actual, mhd::Real expected)
{
    std::fprintf(stderr, "kernel unit failure: %s: %.17g != %.17g\n", name, actual, expected);
    return 1;
}

} // namespace

int main()
{
    using namespace mhd;
    constexpr Real gamma = 5.0 / 3.0;
    const Real q[NPRIM] = {1.25, 0.3, -0.2, 0.1, 0.8, 0.7, -0.4, 0.2};
    Real u[NCONS], recovered[NPRIM];
    prim_to_cons(q, u, gamma);
    cons_to_prim(u, recovered, gamma);
    for (int n = 0; n < NPRIM; ++n)
        if (!close(recovered[n], q[n])) return fail("primitive-conservative round trip", recovered[n], q[n]);

    const Real sound_speed = std::sqrt(gamma * q[QP] / q[QRHO]);
    if (!close(fast_speed(q[QRHO], q[QP], 0.0, 0.0, gamma), sound_speed))
        return fail("fast speed B=0", fast_speed(q[QRHO], q[QP], 0.0, 0.0, gamma), sound_speed);

    if (!close(limited_slope(2.0, 1.0, Limiter::MinMod), 1.0)) return fail("minmod", limited_slope(2.0, 1.0, Limiter::MinMod), 1.0);
    if (!close(limited_slope(-1.0, 1.0, Limiter::MC), 0.0)) return fail("MC extremum", limited_slope(-1.0, 1.0, Limiter::MC), 0.0);
    if (!close(face_value_plus(1.0, 2.0, 3.0, Limiter::MC), 2.5)) return fail("MUSCL plus", face_value_plus(1.0, 2.0, 3.0, Limiter::MC), 2.5);

    if (!close(corner_emf(1.0, 2.0, 3.0, 4.0, 9.0, 9.0, 9.0, 9.0,
                          EmfAveraging::BalsaraSpicer), 2.5))
        return fail("Balsara-Spicer EMF", corner_emf(1.0, 2.0, 3.0, 4.0, 9.0, 9.0, 9.0, 9.0, EmfAveraging::BalsaraSpicer), 2.5);
    if (!close(corner_emf(1.0, 2.0, 3.0, 4.0, 1.0, 2.0, 3.0, 4.0,
                          EmfAveraging::GardinerStone), 2.5))
        return fail("Gardiner-Stone EMF", corner_emf(1.0, 2.0, 3.0, 4.0, 1.0, 2.0, 3.0, 4.0, EmfAveraging::GardinerStone), 2.5);

    std::puts("kernel unit: PASS");
    return 0;
}
