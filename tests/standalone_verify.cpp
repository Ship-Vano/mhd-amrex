//
// standalone_verify.cpp — автономный (без AMReX) верификационный драйвер.
//
// Использует РОВНО ТЕ ЖЕ вычислительные ядра (Hlld.H, Reconstruction.H,
// CtUpdate.H), что и AMReX-слой проекта, на однородной декартовой сетке.
// Назначение: численная проверка схемы (Brio–Wu, Orszag–Tang, альфвеновская
// волна), измерение div B и ошибок сходимости без необходимости собирать AMReX.
//
// Сборка:  g++ -O3 -std=c++17 -I../src/kernels standalone_verify.cpp -o verify
// Запуск:  ./verify briowu | ot <N> | alfven <N> | rotor <N>
//
#include "MhdState.H"
#include "Hlld.H"
#include "Reconstruction.H"
#include "CtUpdate.H"

#include <vector>
#include <string>
#include <cstdio>
#include <cstring>
#include <functional>
#include <algorithm>

using namespace mhd;

// ---------------------------------------------------------------------------
// Простая структура сетки с фантомными ячейками
// ---------------------------------------------------------------------------
enum class BC { Periodic, Outflow, Reflect, Dirichlet };

struct Grid {
    int nx, ny, ng = 2;
    Real x0, x1, y0, y1, dx, dy;
    BC bclo[2], bchi[2];

    std::vector<Real> U;          // (nx+2ng)*(ny+2ng)*NCONS, клеточные величины
    std::vector<Real> bx;         // (nx+1+2ng)*(ny+2ng), x-грани
    std::vector<Real> by;         // (nx+2ng)*(ny+1+2ng), y-грани
    std::vector<Real> U0, bx0, by0;   // состояние на t^n (для RK2 и Dirichlet)

    int sxc, syc, sxf;            // страйды

    Grid(int nx_, int ny_, Real x0_, Real x1_, Real y0_, Real y1_)
        : nx(nx_), ny(ny_), x0(x0_), x1(x1_), y0(y0_), y1(y1_)
    {
        dx = (x1 - x0) / nx; dy = (y1 - y0) / ny;
        sxc = nx + 2*ng; syc = ny + 2*ng; sxf = nx + 1 + 2*ng;
        U.assign(size_t(sxc)*syc*NCONS, 0.0);
        bx.assign(size_t(sxf)*syc, 0.0);
        by.assign(size_t(sxc)*(ny+1+2*ng), 0.0);
    }
    // Доступ: клеточные индексы i∈[-ng, nx+ng), face-индексы x: i∈[-ng, nx+ng]
    Real& u (int i, int j, int n) { return U [ (size_t(n)*syc + (j+ng))*sxc + (i+ng) ]; }
    Real& fx(int i, int j)        { return bx[ size_t(j+ng)*sxf + (i+ng) ]; }
    Real& fy(int i, int j)        { return by[ size_t(j+ng)*sxc + (i+ng) ]; }
    Real xc(int i) const { return x0 + (i + 0.5) * dx; }
    Real yc(int j) const { return y0 + (j + 0.5) * dy; }
};

// ---------------------------------------------------------------------------
// Граничные условия (фантомные ячейки и фантомные грани)
// ---------------------------------------------------------------------------
static void fill_ghosts(Grid& g)
{
    const int ng = g.ng, nx = g.nx, ny = g.ny;
    // --- клеточные величины, x-направление
    for (int j = -ng; j < ny+ng; ++j) {
        for (int k = 1; k <= ng; ++k) {
            for (int n = 0; n < NCONS; ++n) {
                // lo
                switch (g.bclo[0]) {
                case BC::Periodic:  g.u(-k, j, n) = g.u(nx-k, j, n); break;
                case BC::Outflow:   g.u(-k, j, n) = g.u(0, j, n);    break;
                case BC::Reflect:   g.u(-k, j, n) = (n==UMX||n==UBX) ? -g.u(k-1, j, n) : g.u(k-1, j, n); break;
                case BC::Dirichlet: /* заморожено: значения уже лежат в U из НУ */ break;
                }
                // hi
                switch (g.bchi[0]) {
                case BC::Periodic:  g.u(nx-1+k, j, n) = g.u(k-1, j, n); break;
                case BC::Outflow:   g.u(nx-1+k, j, n) = g.u(nx-1, j, n); break;
                case BC::Reflect:   g.u(nx-1+k, j, n) = (n==UMX||n==UBX) ? -g.u(nx-k, j, n) : g.u(nx-k, j, n); break;
                case BC::Dirichlet: break;
                }
            }
        }
    }
    // --- клеточные величины, y-направление
    for (int i = -ng; i < nx+ng; ++i) {
        for (int k = 1; k <= ng; ++k) {
            for (int n = 0; n < NCONS; ++n) {
                switch (g.bclo[1]) {
                case BC::Periodic:  g.u(i, -k, n) = g.u(i, ny-k, n); break;
                case BC::Outflow:   g.u(i, -k, n) = g.u(i, 0, n);    break;
                case BC::Reflect:   g.u(i, -k, n) = (n==UMY||n==UBY) ? -g.u(i, k-1, n) : g.u(i, k-1, n); break;
                case BC::Dirichlet: break;
                }
                switch (g.bchi[1]) {
                case BC::Periodic:  g.u(i, ny-1+k, n) = g.u(i, k-1, n); break;
                case BC::Outflow:   g.u(i, ny-1+k, n) = g.u(i, ny-1, n); break;
                case BC::Reflect:   g.u(i, ny-1+k, n) = (n==UMY||n==UBY) ? -g.u(i, ny-k, n) : g.u(i, ny-k, n); break;
                case BC::Dirichlet: break;
                }
            }
        }
    }
    // --- грани bx: фантомы по обоим направлениям
    for (int j = -ng; j < ny+ng; ++j) {
        for (int k = 1; k <= ng; ++k) {
            if (g.bclo[0] == BC::Periodic) g.fx(-k, j) = g.fx(nx-k, j);
            else if (g.bclo[0] == BC::Outflow) g.fx(-k, j) = g.fx(0, j);
            else if (g.bclo[0] == BC::Reflect) g.fx(-k, j) = -g.fx(k, j);
            if (g.bchi[0] == BC::Periodic) g.fx(nx+k, j) = g.fx(k, j);
            else if (g.bchi[0] == BC::Outflow) g.fx(nx+k, j) = g.fx(nx, j);
            else if (g.bchi[0] == BC::Reflect) g.fx(nx+k, j) = -g.fx(nx-k, j);
        }
    }
    for (int i = -ng; i <= nx+ng-1+1; ++i) {
        for (int k = 1; k <= ng; ++k) {
            if (i > nx) continue;
            if (g.bclo[1] == BC::Periodic) g.fx(i, -k) = g.fx(i, ny-k);
            else if (g.bclo[1] == BC::Outflow) g.fx(i, -k) = g.fx(i, 0);
            else if (g.bclo[1] == BC::Reflect) g.fx(i, -k) = g.fx(i, k-1);
            if (g.bchi[1] == BC::Periodic) g.fx(i, ny-1+k) = g.fx(i, k-1);
            else if (g.bchi[1] == BC::Outflow) g.fx(i, ny-1+k) = g.fx(i, ny-1);
            else if (g.bchi[1] == BC::Reflect) g.fx(i, ny-1+k) = g.fx(i, ny-k);
        }
    }
    // --- грани by
    for (int i = -ng; i < nx+ng; ++i) {
        for (int k = 1; k <= ng; ++k) {
            if (g.bclo[1] == BC::Periodic) g.fy(i, -k) = g.fy(i, ny-k);
            else if (g.bclo[1] == BC::Outflow) g.fy(i, -k) = g.fy(i, 0);
            else if (g.bclo[1] == BC::Reflect) g.fy(i, -k) = -g.fy(i, k);
            if (g.bchi[1] == BC::Periodic) g.fy(i, ny+k) = g.fy(i, k);
            else if (g.bchi[1] == BC::Outflow) g.fy(i, ny+k) = g.fy(i, ny);
            else if (g.bchi[1] == BC::Reflect) g.fy(i, ny+k) = -g.fy(i, ny-k);
        }
    }
    for (int j = -ng; j <= ny; ++j) {
        for (int k = 1; k <= ng; ++k) {
            if (g.bclo[0] == BC::Periodic) g.fy(-k, j) = g.fy(nx-k, j);
            else if (g.bclo[0] == BC::Outflow) g.fy(-k, j) = g.fy(0, j);
            else if (g.bclo[0] == BC::Reflect) g.fy(-k, j) = g.fy(k-1, j);
            if (g.bchi[0] == BC::Periodic) g.fy(nx-1+k, j) = g.fy(k-1, j);
            else if (g.bchi[0] == BC::Outflow) g.fy(nx-1+k, j) = g.fy(nx-1, j);
            else if (g.bchi[0] == BC::Reflect) g.fy(nx-1+k, j) = g.fy(nx-k, j);
        }
    }
}

// Интерполяция граневых B в центры ячеек (RT0 на прямоугольнике → среднее граней)
static void sync_cell_B(Grid& g)
{
    const int ng = g.ng;
    for (int j = -ng; j < g.ny+ng; ++j)
        for (int i = -ng; i < g.nx+ng; ++i) {
            if (i+1 <= g.nx+ng) g.u(i, j, UBX) = 0.5 * (g.fx(i, j) + g.fx(i+1, j));
            if (j+1 <= g.ny+ng) g.u(i, j, UBY) = 0.5 * (g.fy(i, j) + g.fy(i, j+1));
        }
}

// ---------------------------------------------------------------------------
// Один этап явного оператора L(U): возвращает новое состояние U ← U + dt·L(U)
// ---------------------------------------------------------------------------
struct Workspace {
    std::vector<Real> Q;               // примитивы во всех ячейках
    std::vector<Real> Fx, Fy;          // потоки на гранях
    std::vector<Real> Ez;              // узловые ЭДС
};

// Счётчик срабатываний положительность-сохраняющего отката HLLD -> HLL.
// На гладких канонических тестах обязан оставаться нулём (gate T04).
static int g_hlld_fallbacks = 0;

static Real max_signal_dt(Grid& g, Real gamma, Real cfl)
{
    Real dt = 1.0e30;
    for (int j = 0; j < g.ny; ++j)
        for (int i = 0; i < g.nx; ++i) {
            Real uc[NCONS], q[NPRIM];
            for (int n = 0; n < NCONS; ++n) uc[n] = g.u(i, j, n);
            cons_to_prim(uc, q, gamma);
            const Real B2 = q[QBX]*q[QBX] + q[QBY]*q[QBY] + q[QBZ]*q[QBZ];
            const Real cfx = fast_speed(q[QRHO], q[QP], q[QBX], B2, gamma);
            const Real cfy = fast_speed(q[QRHO], q[QP], q[QBY], B2, gamma);
            dt = std::min(dt, g.dx / (std::fabs(q[QU]) + cfx));
            dt = std::min(dt, g.dy / (std::fabs(q[QV]) + cfy));
        }
    return cfl * dt;
}

static void euler_stage(Grid& g, Workspace& w, Real dt, Real gamma,
                        Limiter lim, EmfAveraging emf_mode)
{
    const int ng = g.ng, nx = g.nx, ny = g.ny;
    const int sxc = g.sxc, sxf = g.sxf;

    fill_ghosts(g);
    sync_cell_B(g);
    fill_ghosts(g);   // повторно: фантомные клеточные B после sync

    // Примитивы во всех ячейках (включая фантомные)
    w.Q.assign(size_t(sxc)*g.syc*NPRIM, 0.0);
    auto q_at = [&](int i, int j) -> Real* {
        return &w.Q[ (size_t(j+ng)*sxc + (i+ng)) * NPRIM ];
    };
    for (int j = -ng; j < ny+ng; ++j)
        for (int i = -ng; i < nx+ng; ++i) {
            Real uc[NCONS];
            for (int n = 0; n < NCONS; ++n) uc[n] = g.u(i, j, n);
            cons_to_prim(uc, q_at(i, j), gamma);
        }

    // --- потоки на x-гранях: i∈[0..nx], j∈[-1..ny] -------------------------
    w.Fx.assign(size_t(sxf)*g.syc*NCONS, 0.0);
    auto fxv = [&](int i, int j, int n) -> Real& {
        return w.Fx[ (size_t(n)*g.syc + (j+ng))*sxf + (i+ng) ];
    };
    for (int j = -1; j <= ny; ++j) {
        for (int i = 0; i <= nx; ++i) {
            Real qL[NPRIM], qR[NPRIM];
            for (int n = 0; n < NPRIM; ++n) {
                qL[n] = face_value_plus (q_at(i-2, j)[n], q_at(i-1, j)[n], q_at(i, j)[n], lim);
                qR[n] = face_value_minus(q_at(i-1, j)[n], q_at(i, j)[n], q_at(i+1, j)[n], lim);
            }
            Real f[NCONS];
            hlld_flux(qL, qR, g.fx(i, j), f, gamma, Limits{}, &g_hlld_fallbacks);
            for (int n = 0; n < NCONS; ++n) fxv(i, j, n) = f[n];
        }
    }
    // --- потоки на y-гранях: перестановка осей (x'≡y, y'≡−x → компоненты) --
    w.Fy.assign(size_t(sxc)*(ny+1+2*ng)*NCONS, 0.0);
    auto fyv = [&](int i, int j, int n) -> Real& {
        return w.Fy[ (size_t(n)*(ny+1+2*ng) + (j+ng))*sxc + (i+ng) ];
    };
    for (int j = 0; j <= ny; ++j) {
        for (int i = -1; i <= nx; ++i) {
            Real qL[NPRIM], qR[NPRIM], rL[NPRIM], rR[NPRIM];
            for (int n = 0; n < NPRIM; ++n) {
                qL[n] = face_value_plus (q_at(i, j-2)[n], q_at(i, j-1)[n], q_at(i, j)[n], lim);
                qR[n] = face_value_minus(q_at(i, j-1)[n], q_at(i, j)[n], q_at(i, j+1)[n], lim);
            }
            // Поворот в локальные оси грани: u'=v, v'=−u, Bx'=By, By'=−Bx
            auto rot = [](const Real* q, Real* r) {
                r[QRHO]=q[QRHO]; r[QP]=q[QP]; r[QW]=q[QW]; r[QBZ]=q[QBZ];
                r[QU]= q[QV]; r[QV]=-q[QU];
                r[QBX]= q[QBY]; r[QBY]=-q[QBX];
            };
            rot(qL, rL); rot(qR, rR);
            Real f[NCONS];
            hlld_flux(rL, rR, g.fy(i, j), f, gamma, Limits{}, &g_hlld_fallbacks);
            // Обратный поворот потока импульса/поля
            Real fg[NCONS];
            fg[URHO]=f[URHO]; fg[UENE]=f[UENE]; fg[UMZ]=f[UMZ]; fg[UBZ]=f[UBZ];
            fg[UMX] = -f[UMY];  fg[UMY] = f[UMX];
            fg[UBX] = -f[UBY];  fg[UBY] = f[UBX];   // fg[UBX] = +Ez (поток Bx по y)
            for (int n = 0; n < NCONS; ++n) fyv(i, j, n) = fg[n];
        }
    }

    // --- узловые ЭДС Ez(i+1/2, j+1/2), узлы: i∈[0..nx], j∈[0..ny] ----------
    // ЭДС на x-грани:  Ez = −F[UBY];  ЭДС на y-грани:  Ez = +G[UBX].
    w.Ez.assign(size_t(nx+1)*(ny+1), 0.0);
    auto ez = [&](int i, int j) -> Real& { return w.Ez[ size_t(j)*(nx+1) + i ]; };
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i <= nx; ++i) {
            const Real exm = -fxv(i, j-1, UBY);   // x-грань ниже узла
            const Real exp_= -fxv(i, j,   UBY);   // x-грань выше узла
            const Real eym =  fyv(i-1, j, UBX);   // y-грань левее узла
            const Real eyp =  fyv(i,   j, UBX);   // y-грань правее узла
            ez(i, j) = corner_emf(exm, exp_, eym, eyp,
                                  cell_emf_z(q_at(i-1, j-1)), cell_emf_z(q_at(i, j-1)),
                                  cell_emf_z(q_at(i-1, j)),   cell_emf_z(q_at(i, j)),
                                  emf_mode);
        }

    // --- обновление газовых величин (Годунов: баланс потоков по ячейке) ----
    const Real lx = dt / g.dx, ly = dt / g.dy;
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i < nx; ++i)
            for (int n = 0; n < NCONS; ++n) {
                if (n == UBX || n == UBY) continue;   // эволюционируют через CT
                g.u(i, j, n) -= lx * (fxv(i+1, j, n) - fxv(i, j, n))
                              + ly * (fyv(i, j+1, n) - fyv(i, j, n));
            }

    // --- закон Фарадея по теореме Стокса (бездивергентное обновление) ------
    for (int j = 0; j < ny; ++j)
        for (int i = 0; i <= nx; ++i)
            g.fx(i, j) -= ly * (ez(i, j+1) - ez(i, j));
    for (int j = 0; j <= ny; ++j)
        for (int i = 0; i < nx; ++i)
            g.fy(i, j) += lx * (ez(i+1, j) - ez(i, j));

    sync_cell_B(g);
}

// SSP-RK2 (метод Хойна): U¹ = U + dt·L(U);  Uⁿ⁺¹ = ½(Uⁿ + U¹ + dt·L(U¹)).
// Каждая стадия — CT-обновление, выпуклая комбинация бездивергентных полей
// бездивергентна, поэтому div B = 0 сохраняется точно.
static void rk2_step(Grid& g, Workspace& w, Real dt, Real gamma,
                     Limiter lim, EmfAveraging emf_mode)
{
    g.U0 = g.U; g.bx0 = g.bx; g.by0 = g.by;
    euler_stage(g, w, dt, gamma, lim, emf_mode);
    euler_stage(g, w, dt, gamma, lim, emf_mode);
    for (size_t k = 0; k < g.U.size();  ++k) g.U[k]  = 0.5 * (g.U0[k]  + g.U[k]);
    for (size_t k = 0; k < g.bx.size(); ++k) g.bx[k] = 0.5 * (g.bx0[k] + g.bx[k]);
    for (size_t k = 0; k < g.by.size(); ++k) g.by[k] = 0.5 * (g.by0[k] + g.by[k]);
    sync_cell_B(g);
}

static Real max_divB(Grid& g)
{
    Real m = 0.0;
    for (int j = 0; j < g.ny; ++j)
        for (int i = 0; i < g.nx; ++i) {
            Real d = (g.fx(i+1, j) - g.fx(i, j)) / g.dx
                   + (g.fy(i, j+1) - g.fy(i, j)) / g.dy;
            m = std::max(m, std::fabs(d));
        }
    return m;
}

// ---------------------------------------------------------------------------
// Начальные условия
// ---------------------------------------------------------------------------
static void init_from_prim(Grid& g, Real gamma,
    const std::function<void(Real,Real,Real*)>& prim,        // (x,y) → q[NPRIM]
    const std::function<Real(Real,Real)>& bx_face,           // Bx на x-грани
    const std::function<Real(Real,Real)>& by_face)           // By на y-грани
{
    const int ng = g.ng;
    for (int j = -ng; j < g.ny+ng; ++j)
        for (int i = -ng; i <= g.nx+ng; ++i)
            g.fx(i, j) = bx_face(g.x0 + i*g.dx, g.yc(j));
    for (int j = -ng; j <= g.ny+ng; ++j)
        for (int i = -ng; i < g.nx+ng; ++i)
            g.fy(i, j) = by_face(g.xc(i), g.y0 + j*g.dy);
    for (int j = -ng; j < g.ny+ng; ++j)
        for (int i = -ng; i < g.nx+ng; ++i) {
            Real q[NPRIM], uc[NCONS];
            prim(g.xc(i), g.yc(j), q);
            // клеточные B — из граней (бездивергентная интерполяция)
            q[QBX] = 0.5 * (g.fx(i, j) + g.fx(i+1, j));
            q[QBY] = 0.5 * (g.fy(i, j) + g.fy(i, j+1));
            prim_to_cons(q, uc, gamma);
            for (int n = 0; n < NCONS; ++n) g.u(i, j, n) = uc[n];
        }
}

enum class TimeInt { Euler, RK2 };

static void run(Grid& g, Real gamma, Real cfl, Real tmax,
                Limiter lim, EmfAveraging emf, const char* tag,
                TimeInt ti = TimeInt::RK2)
{
    Workspace w;
    Real t = 0.0; int step = 0; Real divmax_hist = 0.0;
    while (t < tmax) {
        Real dt = std::min(max_signal_dt(g, gamma, cfl), tmax - t);
        if (ti == TimeInt::RK2) rk2_step(g, w, dt, gamma, lim, emf);
        else                    euler_stage(g, w, dt, gamma, lim, emf);
        t += dt; ++step;
        divmax_hist = std::max(divmax_hist, max_divB(g));
        if (step % 100 == 0)
            std::printf("[%s] step %5d  t=%.4f  dt=%.3e  max|divB|=%.3e\n",
                        tag, step, t, dt, max_divB(g));
    }
    std::printf("[%s] DONE: %d steps, t=%.4f, max|divB| over run = %.3e, hlld_fallbacks=%d\n",
                tag, step, t, divmax_hist, g_hlld_fallbacks);
}

static void dump_field(Grid& g, Real gamma, const char* fname)
{
    FILE* f = std::fopen(fname, "w");
    std::fprintf(f, "x,y,rho,u,v,w,p,Bx,By,Bz,divB\n");
    for (int j = 0; j < g.ny; ++j)
        for (int i = 0; i < g.nx; ++i) {
            Real uc[NCONS], q[NPRIM];
            for (int n = 0; n < NCONS; ++n) uc[n] = g.u(i, j, n);
            cons_to_prim(uc, q, gamma);
            Real d = (g.fx(i+1,j)-g.fx(i,j))/g.dx + (g.fy(i,j+1)-g.fy(i,j))/g.dy;
            std::fprintf(f, "%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.8e,%.3e\n",
                         g.xc(i), g.yc(j), q[QRHO], q[QU], q[QV], q[QW], q[QP],
                         q[QBX], q[QBY], q[QBZ], d);
        }
    std::fclose(f);
}

// ===========================================================================
int main(int argc, char** argv)
{
    const std::string test = (argc > 1) ? argv[1] : "briowu";

    if (test == "briowu") {
        // Брио–Ву: γ=2, Bx=0.75, t=0.1; «исторические» (Dirichlet) ГУ по x
        const Real gamma = 2.0;
        Grid g(512, 4, 0.0, 1.0, 0.0, 0.01);
        g.bclo[0] = g.bchi[0] = BC::Dirichlet;
        g.bclo[1] = g.bchi[1] = BC::Periodic;
        init_from_prim(g, gamma,
            [](Real x, Real, Real* q) {
                if (x < 0.5) { q[QRHO]=1.0;   q[QP]=1.0; }
                else         { q[QRHO]=0.125; q[QP]=0.1; }
                q[QU]=q[QV]=q[QW]=0.0; q[QBZ]=0.0;
            },
            [](Real, Real){ return 0.75; },
            [](Real x, Real){ return (x < 0.5) ? 1.0 : -1.0; });
        run(g, gamma, 0.4, 0.1, Limiter::MC, EmfAveraging::GardinerStone, "briowu");
        dump_field(g, gamma, "out_briowu.csv");
    }
    else if (test == "briowu1d") {
        // Параметризованная Брио–Ву-полоса для ablation-сравнения со схемой
        // legacy (N0..N3 из ТЗ). Аргументы:
        //   briowu1d <Nx> <none|minmod|mc|vanleer> <euler|rk2> <bs|gs> <cfl> [out.csv]
        const Real gamma = 2.0;
        const int  Nx    = (argc > 2) ? std::atoi(argv[2]) : 400;
        const std::string slim = (argc > 3) ? argv[3] : "mc";
        const std::string sti  = (argc > 4) ? argv[4] : "rk2";
        const std::string semf = (argc > 5) ? argv[5] : "gs";
        const Real cfl   = (argc > 6) ? std::atof(argv[6]) : 0.1;
        const std::string out = (argc > 7) ? argv[7] : "out_briowu1d.csv";
        Limiter lim = Limiter::MC;
        if      (slim == "none")    lim = Limiter::None;
        else if (slim == "minmod")  lim = Limiter::MinMod;
        else if (slim == "vanleer") lim = Limiter::VanLeer;
        const TimeInt ti = (sti == "euler") ? TimeInt::Euler : TimeInt::RK2;
        const EmfAveraging emf = (semf == "bs") ? EmfAveraging::BalsaraSpicer
                                                : EmfAveraging::GardinerStone;
        Grid g(Nx, 4, 0.0, 1.0, 0.0, 4.0 / Nx);   // 4 ячейки по y, dy = dx
        g.bclo[0] = g.bchi[0] = BC::Dirichlet;
        g.bclo[1] = g.bchi[1] = BC::Periodic;
        init_from_prim(g, gamma,
            [](Real x, Real, Real* q) {
                if (x < 0.5) { q[QRHO]=1.0;   q[QP]=1.0; }
                else         { q[QRHO]=0.125; q[QP]=0.1; }
                q[QU]=q[QV]=q[QW]=0.0; q[QBZ]=0.0;
            },
            [](Real, Real){ return 0.75; },
            [](Real x, Real){ return (x < 0.5) ? 1.0 : -1.0; });
        run(g, gamma, cfl, 0.1, lim, emf, "briowu1d", ti);
        dump_field(g, gamma, out.c_str());
    }
    else if (test == "ot") {
        // Вихрь Орзага–Танга: периодический квадрат, γ=5/3, t=0.5
        const int N = (argc > 2) ? std::atoi(argv[2]) : 192;
        const Real gamma = 5.0/3.0;
        const Real pi = M_PI;
        Grid g(N, N, 0.0, 1.0, 0.0, 1.0);
        g.bclo[0]=g.bchi[0]=g.bclo[1]=g.bchi[1]=BC::Periodic;
        // B = ∇×A, Az = B0( cos(4πx)/(4π) + cos(2πy)/(2π) ) → Bx=−B0 sin2πy, By=B0 sin4πx
        const Real B0 = 1.0/std::sqrt(4.0*M_PI);  // нормировка без 4π в уравнениях
        init_from_prim(g, gamma,
            [pi](Real x, Real y, Real* q) {
                q[QRHO]=25.0/(36.0*pi); q[QP]=5.0/(12.0*pi);
                q[QU]=-std::sin(2*pi*y); q[QV]=std::sin(2*pi*x); q[QW]=0.0; q[QBZ]=0.0;
            },
            [pi,B0](Real, Real y){ return -B0*std::sin(2*pi*y); },
            [pi,B0](Real x, Real){ return  B0*std::sin(4*pi*x); });
        run(g, gamma, 0.4, 0.5, Limiter::MC, EmfAveraging::GardinerStone, "ot");
        dump_field(g, gamma, "out_ot.csv");
    }
    else if (test == "alfven") {
        // Циркулярно поляризованная альфвеновская волна под углом α=30°
        const int N = (argc > 2) ? std::atoi(argv[2]) : 32;
        const Real gamma = 5.0/3.0;
        const Real pi = M_PI;
        const Real alpha = pi/6.0, ca = std::cos(alpha), sa = std::sin(alpha);
        const Real Lx = 1.0/ca, Ly = 1.0/sa;
        const int ny = int(std::lround(N * Ly / Lx));
        Grid g(N, ny, 0.0, Lx, 0.0, Ly);
        g.bclo[0]=g.bchi[0]=g.bclo[1]=g.bchi[1]=BC::Periodic;
        // B∥=1, B⊥=0.1 sinφ, Bz=0.1 cosφ, φ=2π(x cosα + y sinα);  v⊥=B⊥, vz=Bz, v∥=0
        // Векторный потенциал: Az = B∥(y cosα − x sinα) − (0.1/2π) cos φ
        //   → Bx = ∂Az/∂y = B∥cosα − 0.1 sinα·(−sinφ)?  Проверка: Bx = B∥cosα − B⊥ sinα ✓
        auto Az = [=](Real x, Real y) {
            const Real phi = 2*pi*(x*ca + y*sa);
            return ca*y - sa*x + 0.1/(2*pi)*std::cos(phi);
        };
        // Граневые значения как точные средние по грани через разности Az:
        // Bx(face) = (Az(x, y+dy/2) − Az(x, y−dy/2))/dy — бездивергентно машинно.
        const Real hdx = 0.5*g.dx, hdy = 0.5*g.dy;
        init_from_prim(g, gamma,
            [=](Real x, Real y, Real* q) {
                const Real phi = 2*pi*(x*ca + y*sa);
                const Real vperp = 0.1*std::sin(phi), vz = 0.1*std::cos(phi);
                q[QRHO]=1.0; q[QP]=0.1;
                q[QU] = -vperp*sa; q[QV] = vperp*ca; q[QW] = vz;
                q[QBZ] = vz;   // Bz = 0.1 cosφ
            },
            [=](Real x, Real y){ return (Az(x, y+hdy) - Az(x, y-hdy)) / g.dy; },
            [=](Real x, Real y){ return -(Az(x+hdx, y) - Az(x-hdx, y)) / g.dx; });

        // Сохраним точное (начальное) решение для оценки ошибки при t=1
        std::vector<Real> exact = g.U;
        run(g, gamma, 0.4, 1.0, Limiter::MC, EmfAveraging::GardinerStone, "alfven");

        // L1/L2-ошибки по компонентам B⊥, Bz, v⊥, vz (волна вернулась в НУ)
        Real l1 = 0, l2 = 0; int cnt = 0;
        Real l1_bperp = 0, l2_bperp = 0;
        for (int j = 0; j < g.ny; ++j)
            for (int i = 0; i < g.nx; ++i) {
                Real uc[NCONS], q[NPRIM], u0[NCONS], q0[NPRIM];
                for (int n = 0; n < NCONS; ++n) {
                    uc[n] = g.u(i, j, n);
                    u0[n] = exact[(size_t(n)*g.syc + (j+g.ng))*g.sxc + (i+g.ng)];
                }
                cons_to_prim(uc, q, gamma); cons_to_prim(u0, q0, gamma);
                const Real bperp  = q[QBY]*ca - q[QBX]*sa;
                const Real bperp0 = q0[QBY]*ca - q0[QBX]*sa;
                Real esum = std::fabs(bperp - bperp0) + std::fabs(q[QBZ] - q0[QBZ])
                          + std::fabs((q[QV]*ca - q[QU]*sa) - (q0[QV]*ca - q0[QU]*sa))
                          + std::fabs(q[QW] - q0[QW]);
                l1 += esum; l2 += esum*esum; ++cnt;
                l1_bperp += std::fabs(bperp - bperp0);
                l2_bperp += (bperp - bperp0)*(bperp - bperp0);
            }
        std::printf("[alfven N=%d] L1(sum4)=%.6e  L2(sum4)=%.6e  L1(Bperp)=%.6e  L2(Bperp)=%.6e\n",
                    N, l1/cnt, std::sqrt(l2/cnt), l1_bperp/cnt, std::sqrt(l2_bperp/cnt));
        char fn[64]; std::snprintf(fn, 64, "out_alfven_%d.csv", N);
        dump_field(g, gamma, fn);
    }
    else if (test == "rotor") {
        // Вращающийся цилиндр (rotor, Tóth 2000): проверка устойчивости и div B
        const int N = (argc > 2) ? std::atoi(argv[2]) : 128;
        const Real gamma = 1.4;
        Grid g(N, N, 0.0, 1.0, 0.0, 1.0);
        g.bclo[0]=g.bchi[0]=g.bclo[1]=g.bchi[1]=BC::Outflow;
        const Real r0 = 0.1, r1 = 0.115, v0 = 2.0;
        init_from_prim(g, gamma,
            [=](Real x, Real y, Real* q) {
                const Real r = std::sqrt((x-0.5)*(x-0.5) + (y-0.5)*(y-0.5));
                q[QP] = 1.0; q[QW]=0.0; q[QBZ]=0.0;
                if (r < r0) {
                    q[QRHO]=10.0; q[QU]=-v0*(y-0.5)/r0; q[QV]=v0*(x-0.5)/r0;
                } else if (r < r1) {
                    const Real f = (r1 - r)/(r1 - r0);
                    q[QRHO]=1.0+9.0*f; q[QU]=-f*v0*(y-0.5)/r; q[QV]=f*v0*(x-0.5)/r;
                } else { q[QRHO]=1.0; q[QU]=q[QV]=0.0; }
            },
            [](Real, Real){ return 5.0/std::sqrt(4.0*M_PI); },
            [](Real, Real){ return 0.0; });
        run(g, gamma, 0.4, 0.15, Limiter::MC, EmfAveraging::GardinerStone, "rotor");
        dump_field(g, gamma, "out_rotor.csv");
    }
    else {
        std::fprintf(stderr, "unknown test '%s'\n", test.c_str());
        return 1;
    }
    return 0;
}
