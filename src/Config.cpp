// Config.cpp -- строгий разбор JSON-конфигурации.
#include "Config.H"

#include <nlohmann/json.hpp>

#include <cmath>
#include <fstream>
#include <initializer_list>
#include <stdexcept>

namespace mhd {
namespace {

using json = nlohmann::json;

[[noreturn]] void fail(const std::string& path, const std::string& message)
{
    throw std::runtime_error("Invalid config " + path + ": " + message);
}

void require_object(const json& value, const std::string& path)
{
    if (!value.is_object()) fail(path, "must be an object");
}

void reject_unknown_keys(const json& object, std::initializer_list<const char*> allowed,
                         const std::string& path)
{
    require_object(object, path);
    for (const auto& [key, ignored] : object.items()) {
        bool known = false;
        for (const char* candidate : allowed) known = known || key == candidate;
        if (!known) fail(path, "unknown key '" + key + "'");
    }
}

template <class T>
T required_number(const json& value, const std::string& path)
{
    if (!value.is_number()) fail(path, "must be a number");
    return value.get<T>();
}

void read_real_pair(const json& value, std::array<double, 2>& target,
                    const std::string& path)
{
    if (!value.is_array() || value.size() != 2) fail(path, "must contain exactly two numbers");
    target = {required_number<double>(value[0], path + "[0]"),
              required_number<double>(value[1], path + "[1]")};
}

void read_int_pair(const json& value, std::array<int, 2>& target,
                   const std::string& path)
{
    if (!value.is_array() || value.size() != 2) fail(path, "must contain exactly two integers");
    if (!value[0].is_number_integer() || !value[1].is_number_integer())
        fail(path, "must contain exactly two integers");
    target = {value[0].get<int>(), value[1].get<int>()};
}

BcType parse_bc(const std::string& s)
{
    if (s == "periodic") return BcType::Periodic;
    if (s == "outflow") return BcType::Outflow;
    if (s == "reflect") return BcType::Reflect;
    if (s == "dirichlet") return BcType::Dirichlet;
    fail("bc", "unknown boundary type '" + s + "'");
}

Limiter parse_limiter(const std::string& s)
{
    if (s == "none") return Limiter::None;
    if (s == "minmod") return Limiter::MinMod;
    if (s == "mc") return Limiter::MC;
    if (s == "vanleer") return Limiter::VanLeer;
    fail("scheme.limiter", "unknown limiter '" + s + "'");
}

void validate(const SimConfig& c)
{
    if (c.problem.empty()) fail("problem", "must not be empty");
    for (int d = 0; d < 2; ++d) {
        if (!(c.prob_hi[d] > c.prob_lo[d])) fail("geometry", "prob_hi must exceed prob_lo");
        if (c.n_cell[d] <= 0) fail("geometry.n_cell", "entries must be positive");
    }
    if (c.max_level < 0 || c.ref_ratio < 2 || c.regrid_int < 0 ||
        c.blocking_factor <= 0 || c.max_grid_size <= 0 || c.n_error_buf < 0)
        fail("amr", "invalid level, ratio, grid or buffer value");
    if (!(c.cfl > 0.0 && c.cfl <= 1.0) || !(c.t_max >= 0.0) || c.max_steps <= 0)
        fail("time", "require 0 < cfl <= 1, t_max >= 0 and max_steps > 0");
    if (!(c.gamma > 1.0) || !std::isfinite(c.gamma)) fail("scheme.gamma", "must be finite and > 1");
    if (c.plot_int == 0 || c.plot_dt == 0.0 || c.diag_int <= 0)
        fail("output", "plot_int must not be 0, plot_dt must not be 0, diag_int must be positive");
    if (c.output_dir.empty() || c.plot_prefix.empty()) fail("output", "prefix and output_dir must not be empty");
}

} // namespace

SimConfig SimConfig::from_json_file(const std::string& path)
{
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open config file: " + path);
    json j;
    try {
        in >> j;
    } catch (const json::exception& error) {
        throw std::runtime_error("Cannot parse config " + path + ": " + error.what());
    }

    reject_unknown_keys(j, {"problem", "problem_params", "geometry", "amr", "bc", "time", "scheme", "output"}, "root");
    SimConfig c;
    try {
        if (j.contains("problem")) {
            if (!j["problem"].is_string()) fail("problem", "must be a string");
            c.problem = j["problem"].get<std::string>();
        }
        if (j.contains("problem_params")) {
            require_object(j["problem_params"], "problem_params");
            for (const auto& [key, value] : j["problem_params"].items())
                c.pp[key] = required_number<double>(value, "problem_params." + key);
        }
        if (j.contains("geometry")) {
            const auto& g = j["geometry"];
            reject_unknown_keys(g, {"prob_lo", "prob_hi", "n_cell"}, "geometry");
            if (g.contains("prob_lo")) read_real_pair(g["prob_lo"], c.prob_lo, "geometry.prob_lo");
            if (g.contains("prob_hi")) read_real_pair(g["prob_hi"], c.prob_hi, "geometry.prob_hi");
            if (g.contains("n_cell")) read_int_pair(g["n_cell"], c.n_cell, "geometry.n_cell");
        }
        if (j.contains("amr")) {
            const auto& a = j["amr"];
            reject_unknown_keys(a, {"max_level", "ref_ratio", "regrid_int", "blocking_factor", "max_grid_size", "n_error_buf", "refine_grad_rho", "refine_current", "reflux", "regrid_preserve_pressure"}, "amr");
            c.max_level = a.value("max_level", c.max_level);
            c.ref_ratio = a.value("ref_ratio", c.ref_ratio);
            c.regrid_int = a.value("regrid_int", c.regrid_int);
            c.blocking_factor = a.value("blocking_factor", c.blocking_factor);
            c.max_grid_size = a.value("max_grid_size", c.max_grid_size);
            c.n_error_buf = a.value("n_error_buf", c.n_error_buf);
            c.refine_grad_rho = a.value("refine_grad_rho", c.refine_grad_rho);
            c.refine_current = a.value("refine_current", c.refine_current);
            c.reflux = a.value("reflux", c.reflux);
            c.regrid_preserve_pressure = a.value("regrid_preserve_pressure", c.regrid_preserve_pressure);
        }
        if (j.contains("bc")) {
            const auto& b = j["bc"];
            reject_unknown_keys(b, {"x_lo", "x_hi", "y_lo", "y_hi"}, "bc");
            c.bc_xlo = parse_bc(b.value("x_lo", "periodic"));
            c.bc_xhi = parse_bc(b.value("x_hi", "periodic"));
            c.bc_ylo = parse_bc(b.value("y_lo", "periodic"));
            c.bc_yhi = parse_bc(b.value("y_hi", "periodic"));
            if ((c.bc_xlo == BcType::Periodic) != (c.bc_xhi == BcType::Periodic) ||
                (c.bc_ylo == BcType::Periodic) != (c.bc_yhi == BcType::Periodic))
                fail("bc", "periodic boundary must be set on both sides of a dimension");
        }
        if (j.contains("time")) {
            const auto& t = j["time"];
            reject_unknown_keys(t, {"cfl", "t_max", "max_steps", "integrator"}, "time");
            c.cfl = t.value("cfl", c.cfl);
            c.t_max = t.value("t_max", c.t_max);
            c.max_steps = t.value("max_steps", c.max_steps);
            const std::string integrator = t.value("integrator", "rk2");
            if (integrator == "rk2") c.integrator = Integrator::RK2;
            else if (integrator == "euler") c.integrator = Integrator::Euler;
            else fail("time.integrator", "unknown integrator '" + integrator + "'");
        }
        if (j.contains("scheme")) {
            const auto& s = j["scheme"];
            reject_unknown_keys(s, {"gamma", "limiter", "emf_averaging"}, "scheme");
            c.gamma = s.value("gamma", c.gamma);
            c.limiter = parse_limiter(s.value("limiter", "mc"));
            const std::string emf = s.value("emf_averaging", "gardiner_stone");
            if (emf == "balsara_spicer") c.emf = EmfAveraging::BalsaraSpicer;
            else if (emf == "gardiner_stone") c.emf = EmfAveraging::GardinerStone;
            else fail("scheme.emf_averaging", "unknown mode '" + emf + "'");
        }
        if (j.contains("output")) {
            const auto& o = j["output"];
            reject_unknown_keys(o, {"plot_int", "plot_dt", "prefix", "output_dir", "format", "diag_int", "write_plotfiles"}, "output");
            c.plot_int = o.value("plot_int", c.plot_int);
            c.plot_dt = o.value("plot_dt", c.plot_dt);
            c.plot_prefix = o.value("prefix", c.plot_prefix);
            c.output_dir = o.value("output_dir", c.output_dir);
            c.diag_int = o.value("diag_int", c.diag_int);
            c.write_plotfiles = o.value("write_plotfiles", c.write_plotfiles);
            const std::string format = o.value("format", "native");
            if (format == "native") c.format = OutputFormat::Native;
            else if (format == "hdf5") c.format = OutputFormat::Hdf5;
            else fail("output.format", "unknown format '" + format + "'");
        }
        validate(c);
    } catch (const json::exception& error) {
        throw std::runtime_error("Invalid config " + path + ": " + error.what());
    }
    return c;
}

} // namespace mhd
