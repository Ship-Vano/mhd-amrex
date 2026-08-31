#!/usr/bin/env python3
"""Compute compact, area-weighted diagnostics from a legacy VTU result."""

from __future__ import annotations

import argparse
import json
import math
import xml.etree.ElementTree as ET
from pathlib import Path


def numbers(text: str | None, cast=float):
    return [cast(word) for word in (text or "").split()]


def read_mesh(mesh: Path):
    nodes: dict[int, tuple[float, float]] = {}
    elements: list[tuple[int, int, int]] = []
    section = ""
    for line in mesh.read_text(encoding="utf-8").splitlines():
        if line in {"$Nodes", "$Elements"}:
            section = line
            continue
        if line in {"$EndNodes", "$EndElements"}:
            section = ""
            continue
        columns = line.split()
        if section == "$Nodes" and len(columns) == 4:
            nodes[int(columns[0])] = (float(columns[1]), float(columns[2]))
        elif section == "$Elements" and len(columns) == 5:
            elements.append(tuple(map(int, columns[2:5])))
    return nodes, elements


def read_vtu(vtu: Path):
    root = ET.parse(vtu).getroot()
    arrays = root.findall(".//DataArray")
    state = next((a for a in arrays if a.attrib.get("Name") == "elemUs"), None)
    if state is None:
        raise ValueError(f"{vtu}: elemUs cell data is absent")
    values = numbers(state.text)
    if len(values) % 9:
        raise ValueError(f"{vtu}: elemUs component count is not a multiple of 9")
    return [values[i : i + 9] for i in range(0, len(values), 9)]


def primitive(u: list[float], gamma: float):
    rho = u[0]
    if rho <= 0 or not math.isfinite(rho):
        return rho, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan, math.nan
    vx, vy, vz = u[1] / rho, u[2] / rho, u[3] / rho
    bx, by, bz = u[5], u[6], u[7]
    pressure = (gamma - 1.0) * (u[4] - 0.5 * rho * (vx * vx + vy * vy + vz * vz) - 0.5 * (bx * bx + by * by + bz * bz))
    return rho, pressure, vx, vy, vz, bx, by, bz


def diagnostics(mesh: Path, vtu: Path, case: str, gamma: float,
                reference_vtu: Path | None = None):
    nodes, elements = read_mesh(mesh)
    states = read_vtu(vtu)
    if len(elements) != len(states):
        raise ValueError(f"mesh has {len(elements)} triangles, VTU has {len(states)} states")
    reference_states = read_vtu(reference_vtu) if reference_vtu else None
    if reference_states is not None and len(reference_states) != len(states):
        raise ValueError(f"{reference_vtu}: cell count differs from {vtu}")
    area_total = 0.0
    mag_energy = 0.0
    mass = momentum_x = momentum_y = momentum_z = total_energy = 0.0
    bz_abs_max = 0.0
    rho_values, p_values = [], []
    cp_l1_numer, cp_l2_numer, cp_linf = 0.0, 0.0, 0.0
    cp_component_l1 = {
        name: 0.0 for name in (
            "rho", "pressure", "v_parallel", "v_perp", "vz",
            "b_parallel", "b_perp", "bz",
        )
    }
    cp_component_l2 = {name: 0.0 for name in cp_component_l1}
    cp_component_linf = {name: 0.0 for name in cp_component_l1}
    cp_component_reference_l1 = {name: 0.0 for name in cp_component_l1}
    loop_return_l1 = loop_return_l2 = loop_return_linf = loop_reference_l1 = 0.0
    alpha = math.pi / 6.0
    for index, (element, state) in enumerate(zip(elements, states)):
        (x0, y0), (x1, y1), (x2, y2) = (nodes[index] for index in element)
        area = 0.5 * abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0))
        cx, cy = (x0 + x1 + x2) / 3.0, (y0 + y1 + y2) / 3.0
        rho, pressure, vx, vy, vz, bx, by, bz = primitive(state, gamma)
        area_total += area
        mag_energy += area * 0.5 * (bx * bx + by * by + bz * bz)
        mass += area * state[0]
        momentum_x += area * state[1]
        momentum_y += area * state[2]
        momentum_z += area * state[3]
        total_energy += area * state[4]
        bz_abs_max = max(bz_abs_max, abs(bz))
        rho_values.append(rho)
        p_values.append(pressure)
        if case == "cp_alfven":
            phase = 2.0 * math.pi * (cx * math.cos(alpha) + cy * math.sin(alpha))
            sine = math.sin(phase)
            cosine = math.cos(phase)
            # The legacy task uses Toth's 30-degree travelling wave.  Work
            # in (parallel, perpendicular, z) coordinates so the metric is
            # comparable to Table II of Toth (2000), rather than conflating
            # a physical rotation with numerical error in Cartesian Bx/By.
            observed = {
                "rho": rho,
                "pressure": pressure,
                "v_parallel": vx * math.cos(alpha) + vy * math.sin(alpha),
                "v_perp": -vx * math.sin(alpha) + vy * math.cos(alpha),
                "vz": vz,
                "b_parallel": bx * math.cos(alpha) + by * math.sin(alpha),
                "b_perp": -bx * math.sin(alpha) + by * math.cos(alpha),
                "bz": bz,
            }
            reference = {
                "rho": 1.0,
                "pressure": 0.1,
                "v_parallel": 0.0,
                "v_perp": 0.1 * sine,
                "vz": 0.1 * cosine,
                "b_parallel": 1.0,
                "b_perp": 0.1 * sine,
                "bz": 0.1 * cosine,
            }
            component_error = {
                name: abs(observed[name] - reference[name]) for name in observed
            }
            err = max(component_error.values())
            cp_l1_numer += area * err
            cp_l2_numer += area * err * err
            cp_linf = max(cp_linf, err)
            for name, error in component_error.items():
                cp_component_l1[name] += area * error
                cp_component_l2[name] += area * error * error
                cp_component_linf[name] = max(cp_component_linf[name], error)
                cp_component_reference_l1[name] += area * abs(reference[name])
        if case == "magnetic_loop" and reference_states is not None:
            _, _, _, _, _, reference_bx, reference_by, _ = primitive(reference_states[index], gamma)
            magnetic_difference = math.hypot(bx - reference_bx, by - reference_by)
            loop_return_l1 += area * magnetic_difference
            loop_return_l2 += area * magnetic_difference * magnetic_difference
            loop_return_linf = max(loop_return_linf, magnetic_difference)
            loop_reference_l1 += area * math.hypot(reference_bx, reference_by)
    finite = all(math.isfinite(value) for state in states for value in state)
    result = {
        "case": case,
        "cells": len(states),
        "area": area_total,
        "finite": finite,
        "rho_min": min(rho_values),
        "rho_max": max(rho_values),
        "pressure_min": min(p_values),
        "pressure_max": max(p_values),
        "magnetic_energy": mag_energy if finite else None,
        "mass": mass if finite else None,
        "momentum_x": momentum_x if finite else None,
        "momentum_y": momentum_y if finite else None,
        "momentum_z": momentum_z if finite else None,
        "total_energy": total_energy if finite else None,
        "bz_abs_max": bz_abs_max if finite else None,
    }
    if case == "cp_alfven":
        result["cp_return_error_l1"] = cp_l1_numer / area_total if finite else None
        result["cp_return_error_l2"] = math.sqrt(cp_l2_numer / area_total) if finite else None
        result["cp_return_error_linf"] = cp_linf if finite else None
        result["cp_component_errors"] = {
            name: {
                "l1": cp_component_l1[name] / area_total if finite else None,
                "l2": math.sqrt(cp_component_l2[name] / area_total) if finite else None,
                "linf": cp_component_linf[name] if finite else None,
                "relative_l1": (
                    cp_component_l1[name] / cp_component_reference_l1[name]
                    if finite and cp_component_reference_l1[name] > 0.0 else None
                ),
            }
            for name in cp_component_l1
        }
        # Equation (45) and Table II in Toth (2000) average relative L1
        # errors over the nonzero wave variables v_perp, vz, B_perp and Bz.
        toth_components = ("v_perp", "vz", "b_perp", "bz")
        result["cp_toth_average_relative_l1"] = (
            sum(cp_component_l1[name] / cp_component_reference_l1[name]
                for name in toth_components) / len(toth_components)
            if finite and all(cp_component_reference_l1[name] > 0.0 for name in toth_components)
            else None
        )
    if case == "magnetic_loop" and reference_states is not None:
        result["magnetic_loop_return_b_error"] = {
            "l1": loop_return_l1 / area_total if finite else None,
            "l2": math.sqrt(loop_return_l2 / area_total) if finite else None,
            "linf": loop_return_linf if finite else None,
            "relative_l1": (
                loop_return_l1 / loop_reference_l1
                if finite and loop_reference_l1 > 0.0 else None
            ),
            "reference": str(reference_vtu),
        }
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--vtu", type=Path, required=True)
    parser.add_argument("--case", choices=("brio_wu", "cp_alfven", "magnetic_loop"), required=True)
    parser.add_argument("--gamma", type=float, required=True)
    parser.add_argument("--reference-vtu", type=Path,
                        help="same-mesh reference state for a return-error diagnostic")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(
        json.dumps(diagnostics(args.mesh, args.vtu, args.case, args.gamma, args.reference_vtu),
                   indent=2, sort_keys=True, allow_nan=False) + "\n"
    )


if __name__ == "__main__":
    main()
