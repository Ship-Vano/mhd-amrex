#ifndef MHD_RECONSTRUCTION_H_
#define MHD_RECONSTRUCTION_H_
//
// Reconstruction.H — кусочно-линейная MUSCL-реконструкция примитивных
// переменных с TVD-ограничителями. Повышает схему до 2-го порядка по
// пространству (исходная схема Авдеевой–Лукина — 1-го порядка; реконструкция
// её естественное улучшение и отключается опцией reconstruction="constant").
//
#include "MhdState.H"

namespace mhd {

enum class Limiter : int { None = 0, MinMod = 1, MC = 2, VanLeer = 3 };

MHD_HD Real limited_slope(Real dl, Real dr, Limiter lim)
{
    switch (lim) {
    case Limiter::None:
        return Real(0.0);                       // кусочно-постоянная (1-й порядок)
    case Limiter::MinMod: {
        if (dl * dr <= Real(0.0)) return Real(0.0);
        return (std::fabs(dl) < std::fabs(dr)) ? dl : dr;
    }
    case Limiter::VanLeer: {
        if (dl * dr <= Real(0.0)) return Real(0.0);
        return Real(2.0) * dl * dr / (dl + dr);
    }
    case Limiter::MC:
    default: {
        if (dl * dr <= Real(0.0)) return Real(0.0);
        const Real dc = Real(0.5) * (dl + dr);
        const Real d2 = Real(2.0) * ((std::fabs(dl) < std::fabs(dr)) ? dl : dr);
        return (std::fabs(dc) < std::fabs(d2)) ? dc : d2;
    }
    }
}

// Реконструкция значения переменной на левой/правой стороне грани i+1/2
// из трёх соседних клеточных значений qm=q_{i-1}, q0=q_i, qp=q_{i+1}.
MHD_HD Real face_value_plus (Real qm, Real q0, Real qp, Limiter lim)  // сторона "L" грани i+1/2
{
    return q0 + Real(0.5) * limited_slope(q0 - qm, qp - q0, lim);
}
MHD_HD Real face_value_minus(Real qm, Real q0, Real qp, Limiter lim)  // сторона "R" грани i-1/2
{
    return q0 - Real(0.5) * limited_slope(q0 - qm, qp - q0, lim);
}

} // namespace mhd
#endif
