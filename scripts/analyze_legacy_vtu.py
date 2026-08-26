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


def diagnostics(mesh: Path, vtu: Path, case: str, gamma: float):
    nodes, elements = read_mesh(mesh)
    states = read_vtu(vtu)
    if len(elements) != len(states):
        raise ValueError(f"mesh has {len(elements)} triangles, VTU has {len(states)} states")
    area_total = 0.0
    mag_energy = 0.0
    rho_values, p_values = [], []
    cp_l1_numer, cp_l2_numer, cp_linf = 0.0, 0.0, 0.0
    alpha = math.pi / 6.0
    for element, state in zip(elements, states):
        (x0, y0), (x1, y1), (x2, y2) = (nodes[index] for index in element)
        area = 0.5 * abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0))
        cx, cy = (x0 + x1 + x2) / 3.0, (y0 + y1 + y2) / 3.0
        rho, pressure, vx, vy, vz, bx, by, bz = primitive(state, gamma)
        area_total += area
        mag_energy += area * 0.5 * (bx * bx + by * by + bz * bz)
        rho_values.append(rho)
        p_values.append(pressure)
        if case == "cp_alfven":
            phase = 2.0 * math.pi * (cx * math.cos(alpha) + cy * math.sin(alpha))
            reference = [
                1.0,
                -0.1 * math.sin(phase) * math.sin(alpha),
                0.1 * math.sin(phase) * math.cos(alpha),
                0.1 * math.cos(phase),
                math.cos(alpha) - 0.1 * math.sin(phase) * math.sin(alpha),
                math.sin(alpha) + 0.1 * math.sin(phase) * math.cos(alpha),
                0.1 * math.cos(phase),
            ]
            observed = [rho, vx, vy, vz, bx, by, bz]
            err = max(abs(a - b) for a, b in zip(observed, reference))
            cp_l1_numer += area * err
            cp_l2_numer += area * err * err
            cp_linf = max(cp_linf, err)
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
    }
    if case == "cp_alfven":
        result["cp_return_error_l1"] = cp_l1_numer / area_total if finite else None
        result["cp_return_error_l2"] = math.sqrt(cp_l2_numer / area_total) if finite else None
        result["cp_return_error_linf"] = cp_linf if finite else None
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--vtu", type=Path, required=True)
    parser.add_argument("--case", choices=("brio_wu", "cp_alfven", "magnetic_loop"), required=True)
    parser.add_argument("--gamma", type=float, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(json.dumps(diagnostics(args.mesh, args.vtu, args.case, args.gamma), indent=2, sort_keys=True, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
