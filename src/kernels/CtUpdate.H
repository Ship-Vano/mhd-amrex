#ifndef MHD_CT_H_
#define MHD_CT_H_
//
// CtUpdate.H — бездивергентное обновление магнитного поля (constrained
// transport). Это декартов аналог схемы Авдеевой–Лукина (2019): там нормальные
// компоненты B живут на рёбрах треугольников, ЭДС усредняется в вершины из
// HLLD-потоков соседних рёбер, и закон Фарадея аппроксимируется напрямую через
// теорему Стокса. На прямоугольной сетке:
//
//   * Bx(i±1/2, j) — на x-гранях, By(i, j±1/2) — на y-гранях;
//   * HLLD-поток на x-грани содержит компоненту F[By'] = u·By − v·Bx = −Ez,
//     на y-грани — компоненту G[Bx'] = v·Bx − u·By = +Ez (в глобальных осях);
//   * ЭДС в узле (i+1/2, j+1/2) — среднее ЭДС четырёх примыкающих граней
//     (точно как усреднение по рёбрам, сходящимся в вершину, в статье);
//   * закон Фарадея по Стоксу для каждой грани:
//        Bx^{n+1}_{i+1/2,j} = Bx − (Δt/Δy)·(Ez_{i+1/2,j+1/2} − Ez_{i+1/2,j−1/2})
//        By^{n+1}_{i,j+1/2} = By + (Δt/Δx)·(Ez_{i+1/2,j+1/2} − Ez_{i−1/2,j+1/2})
//
// Дискретная дивергенция  (ΔBx/Δx + ΔBy/Δy)  при этом сохраняется точно
// (телескопическое сокращение узловых ЭДС) — div B = 0 на уровне ошибок
// округления без каких-либо чисток.
//
// Вариант усреднения ЭДС выбирается опцией:
//   * BalsaraSpicer — простое арифметическое среднее 4 граней (как в статье);
//   * GardinerStone — «ε=2»-вариант Gardiner & Stone (2005):
//        Ez_corner = 2·<Ez_faces> − <Ez_cells>,
//     где <Ez_cells> — среднее клеточных значений Ez = v·Bx − u·By по четырём
//     ячейкам вокруг узла. Восстанавливает правильный порядок аппроксимации
//     и добавляет диссипацию, полезную при CFL≈0.5.
//
#include "MhdState.H"

namespace mhd {

enum class EmfAveraging : int { BalsaraSpicer = 0, GardinerStone = 1 };

// ЭДС в узле из четырёх граневых значений (e_xlo/e_xhi — ЭДС на x-гранях ниже/выше
// узла по y; e_ylo/e_yhi — на y-гранях левее/правее узла по x) и четырёх клеточных.
MHD_HD Real corner_emf(Real ez_xface_jm, Real ez_xface_jp,
                       Real ez_yface_im, Real ez_yface_ip,
                       Real ez_cell_mm, Real ez_cell_pm,
                       Real ez_cell_mp, Real ez_cell_pp,
                       EmfAveraging mode)
{
    const Real face_avg = Real(0.25) * (ez_xface_jm + ez_xface_jp + ez_yface_im + ez_yface_ip);
    if (mode == EmfAveraging::BalsaraSpicer) {
        return face_avg;
    }
    const Real cell_avg = Real(0.25) * (ez_cell_mm + ez_cell_pm + ez_cell_mp + ez_cell_pp);
    return Real(2.0) * face_avg - cell_avg;
}

// Клеточная ЭДС Ez = -(v×B)_z = v·Bx − u·By
MHD_HD Real cell_emf_z(const Real* q) { return q[QV]*q[QBX] - q[QU]*q[QBY]; }

} // namespace mhd
#endif
