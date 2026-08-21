//
// main.cpp — точка входа.  Запуск:  ./mhd2d <config.json>
//
// Вся постановка задачи (НУ, ГУ, сетка, AMR, схема, вывод) задаётся JSON-файлом
// без перекомпиляции; примеры — в каталоге inputs/.
//
#include <vector>
#include <string>
#include <cstdio>

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#include "Config.H"
#include "MhdAmr.H"

int main(int argc, char* argv[])
{
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s <config.json> [amrex.opt=value ...]\n", argv[0]);
        return 1;
    }
    const std::string cfg_path = argv[1];

    // ВАЖНО: AMReX трактует первый позиционный аргумент командной строки как
    // СВОЙ inputs-файл и парсит его в формате ParmParse — JSON он не понимает
    // (abort "value with no defn: {"). Поэтому путь к JSON исключается из
    // argv, а остальные аргументы (вида amrex.v=1, key=value) передаются
    // дальше — ими можно переопределять служебные параметры AMReX.
    std::vector<char*> args;
    args.push_back(argv[0]);
    for (int i = 2; i < argc; ++i) args.push_back(argv[i]);
    int    amrex_argc = static_cast<int>(args.size());
    char** amrex_argv = args.data();

    // AMReX-инициализация (MPI стартует внутри, если включён).
    amrex::Initialize(amrex_argc, amrex_argv);
    {
        mhd::SimConfig cfg = mhd::SimConfig::from_json_file(cfg_path);

        amrex::Print() << "=== 2D ideal MHD (Avdeeva–Lukin staggered CT + HLLD) ===\n"
                       << "problem      : " << cfg.problem << "\n"
                       << "base grid    : " << cfg.n_cell[0] << " x " << cfg.n_cell[1] << "\n"
                       << "max AMR level: " << cfg.max_level << "\n"
                       << "CFL          : " << cfg.cfl << ",  t_max: " << cfg.t_max << "\n";

        mhd::MhdAmr sim(cfg);
        sim.InitData();
        sim.Evolve();
    }
    amrex::Finalize();
    return 0;
}
