#ifndef MHD_CONFIG_H_
#define MHD_CONFIG_H_
//
// Config.H — конфигурация расчёта, читаемая из JSON-файла (nlohmann/json).
// Все параметры задачи (НУ, ГУ, сетка, AMR, схема, вывод) задаются без
// перекомпиляции:   ./mhd2d inputs/orszag_tang.json
//
#include <array>
#include <string>
#include <map>

#include "kernels/MhdState.H"
#include "kernels/Reconstruction.H"
#include "kernels/CtUpdate.H"

namespace mhd {

enum class BcType : int { Periodic, Outflow, Reflect, Dirichlet };
enum class Integrator : int { Euler, RK2 };           // см. обсуждение в REPORT.md
enum class OutputFormat : int { Native, Hdf5 };

struct SimConfig {
    // --- задача ---------------------------------------------------------
    std::string problem = "orszag_tang";
    std::map<std::string, double> pp;     // problem_params: произвольные числа

    // --- геометрия и базовая сетка ---------------------------------------
    std::array<double, 2> prob_lo {0.0, 0.0};
    std::array<double, 2> prob_hi {1.0, 1.0};
    std::array<int, 2>    n_cell  {256, 256};

    // --- AMR --------------------------------------------------------------
    int max_level       = 0;       // 0 — однородная сетка
    int ref_ratio       = 2;
    int regrid_int      = 4;       // период перестроения сетки (в шагах)
    int blocking_factor = 8;
    int max_grid_size   = 64;
    int n_error_buf     = 2;
    double refine_grad_rho = 0.25; // порог |Δρ|/ρ для тегирования
    double refine_current  = -1.0; // порог |jz|·h (выкл. если < 0)

    // --- граничные условия (по сторонам) ----------------------------------
    BcType bc_xlo = BcType::Periodic, bc_xhi = BcType::Periodic;
    BcType bc_ylo = BcType::Periodic, bc_yhi = BcType::Periodic;

    // --- время -------------------------------------------------------------
    double cfl       = 0.5;
    double t_max     = 0.5;
    int    max_steps = 1000000;
    Integrator integrator = Integrator::RK2;

    // --- схема --------------------------------------------------------------
    double gamma = 5.0 / 3.0;
    Limiter limiter = Limiter::MC;
    EmfAveraging emf = EmfAveraging::GardinerStone;

    // --- вывод ---------------------------------------------------------------
    int    plot_int  = -1;        // каждые N шагов (если > 0)
    double plot_dt   = -1.0;      // либо каждые plot_dt по времени
    std::string plot_prefix = "plt";
    std::string output_dir  = ".";  // директория для plotfile'ов
    OutputFormat format = OutputFormat::Native;
    int diag_int = 10;            // период печати div B / минимумов
    bool write_plotfiles = true;  // false: compute-only / benchmark mode

    static SimConfig from_json_file(const std::string& path);
    double param(const std::string& key, double def) const {
        auto it = pp.find(key);
        return (it == pp.end()) ? def : it->second;
    }
};

} // namespace mhd
#endif
