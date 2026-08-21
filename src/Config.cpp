//
// Config.cpp — разбор JSON-конфигурации (nlohmann/json, подтягивается CMake'ом).
//
#include "Config.H"

#include <nlohmann/json.hpp>
#include <fstream>
#include <stdexcept>

namespace mhd {

namespace {

BcType parse_bc(const std::string& s)
{
    if (s == "periodic")  return BcType::Periodic;
    if (s == "outflow")   return BcType::Outflow;
    if (s == "reflect")   return BcType::Reflect;
    if (s == "dirichlet") return BcType::Dirichlet;   // «исторические» ГУ: значения НУ заморожены
    throw std::runtime_error("Unknown bc type: " + s);
}

Limiter parse_limiter(const std::string& s)
{
    if (s == "none")    return Limiter::None;     // 1-й порядок, как в исходной статье
    if (s == "minmod")  return Limiter::MinMod;
    if (s == "mc")      return Limiter::MC;
    if (s == "vanleer") return Limiter::VanLeer;
    throw std::runtime_error("Unknown limiter: " + s);
}

} // namespace

SimConfig SimConfig::from_json_file(const std::string& path)
{
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open config file: " + path);
    nlohmann::json j; in >> j;

    SimConfig c;
    c.problem = j.value("problem", c.problem);

    if (j.contains("problem_params"))
        for (auto& [k, v] : j["problem_params"].items())
            c.pp[k] = v.get<double>();

    if (j.contains("geometry")) {
        auto& g = j["geometry"];
        if (g.contains("prob_lo")) { c.prob_lo = { g["prob_lo"][0], g["prob_lo"][1] }; }
        if (g.contains("prob_hi")) { c.prob_hi = { g["prob_hi"][0], g["prob_hi"][1] }; }
        if (g.contains("n_cell"))  { c.n_cell  = { g["n_cell"][0],  g["n_cell"][1]  }; }
    }
    if (j.contains("amr")) {
        auto& a = j["amr"];
        c.max_level       = a.value("max_level",       c.max_level);
        c.ref_ratio       = a.value("ref_ratio",       c.ref_ratio);
        c.regrid_int      = a.value("regrid_int",      c.regrid_int);
        c.blocking_factor = a.value("blocking_factor", c.blocking_factor);
        c.max_grid_size   = a.value("max_grid_size",   c.max_grid_size);
        c.n_error_buf     = a.value("n_error_buf",     c.n_error_buf);
        c.refine_grad_rho = a.value("refine_grad_rho", c.refine_grad_rho);
        c.refine_current  = a.value("refine_current",  c.refine_current);
    }
    if (j.contains("bc")) {
        auto& b = j["bc"];
        c.bc_xlo = parse_bc(b.value("x_lo", "periodic"));
        c.bc_xhi = parse_bc(b.value("x_hi", "periodic"));
        c.bc_ylo = parse_bc(b.value("y_lo", "periodic"));
        c.bc_yhi = parse_bc(b.value("y_hi", "periodic"));
        // Периодичность обязана быть парной — иначе геометрия AMReX некорректна
        if ((c.bc_xlo == BcType::Periodic) != (c.bc_xhi == BcType::Periodic) ||
            (c.bc_ylo == BcType::Periodic) != (c.bc_yhi == BcType::Periodic))
            throw std::runtime_error("Periodic BC must be set on both sides of a dimension");
    }
    if (j.contains("time")) {
        auto& t = j["time"];
        c.cfl       = t.value("cfl",       c.cfl);
        c.t_max     = t.value("t_max",     c.t_max);
        c.max_steps = t.value("max_steps", c.max_steps);
        const std::string integ = t.value("integrator", "rk2");
        if      (integ == "rk2")   c.integrator = Integrator::RK2;
        else if (integ == "euler") c.integrator = Integrator::Euler;
        else throw std::runtime_error("Unknown integrator: " + integ);
    }
    if (j.contains("scheme")) {
        auto& s = j["scheme"];
        c.gamma   = s.value("gamma", c.gamma);
        c.limiter = parse_limiter(s.value("limiter", "mc"));
        const std::string e = s.value("emf_averaging", "gardiner_stone");
        c.emf = (e == "balsara_spicer") ? EmfAveraging::BalsaraSpicer
                                        : EmfAveraging::GardinerStone;
    }
    if (j.contains("output")) {
        auto& o = j["output"];
        c.plot_int    = o.value("plot_int",    c.plot_int);
        c.plot_dt     = o.value("plot_dt",     c.plot_dt);
        c.plot_prefix = o.value("prefix",      c.plot_prefix);
        c.output_dir  = o.value("output_dir",  c.output_dir);
        c.diag_int    = o.value("diag_int",    c.diag_int);
        const std::string f = o.value("format", "native");
        c.format = (f == "hdf5") ? OutputFormat::Hdf5 : OutputFormat::Native;
    }
    return c;
}

} // namespace mhd
