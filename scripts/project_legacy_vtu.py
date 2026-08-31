#!/usr/bin/env python3
"""Area-weight a legacy 2-D VTU solution into a reproducible 1-D profile.

This is primarily for Brio--Wu on an irregular triangular strip: values are
averaged by triangle area in x-bins, so the resulting profile is not a slice
chosen by eye.  The output is a raw CSV plus a compact JSON summary; neither
is an assertion that a qualitative literature plot is a pointwise reference.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

from analyze_legacy_vtu import primitive, read_mesh, read_vtu


FIELDS = ("rho", "pressure", "vx", "vy", "vz", "bx", "by", "bz")


def triangle_centroid_and_area(nodes: dict[int, tuple[float, float]], triangle: tuple[int, int, int]):
    (x0, y0), (x1, y1), (x2, y2) = (nodes[index] for index in triangle)
    area = 0.5 * abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0))
    return (x0 + x1 + x2) / 3.0, (y0 + y1 + y2) / 3.0, area


def project(mesh: Path, vtu: Path, gamma: float, bins: int):
    nodes, elements = read_mesh(mesh)
    states = read_vtu(vtu)
    if len(elements) != len(states):
        raise ValueError(f"mesh has {len(elements)} triangles but VTU has {len(states)} cell states")
    if bins < 4:
        raise ValueError("at least four x-bins are required")
    xlo = min(x for x, _ in nodes.values())
    xhi = max(x for x, _ in nodes.values())
    if not xhi > xlo:
        raise ValueError("mesh has zero x extent")

    accum = [{"area": 0.0, "x_area": 0.0, "cells": 0,
              **{field: 0.0 for field in FIELDS},
              **{f"{field}_sq": 0.0 for field in FIELDS}}
             for _ in range(bins)]
    for element, state in zip(elements, states):
        x, _, area = triangle_centroid_and_area(nodes, element)
        values = primitive(state, gamma)
        if not all(math.isfinite(value) for value in values):
            raise ValueError("cannot project a non-finite state")
        index = min(bins - 1, int((x - xlo) / (xhi - xlo) * bins))
        bucket = accum[index]
        bucket["area"] += area
        bucket["x_area"] += area * x
        bucket["cells"] += 1
        for field, value in zip(FIELDS, values):
            bucket[field] += area * value
            bucket[f"{field}_sq"] += area * value * value

    rows = []
    for index, bucket in enumerate(accum):
        if bucket["area"] == 0.0:
            continue
        row = {
            "bin": index,
            "x": bucket["x_area"] / bucket["area"],
            "area": bucket["area"],
            "cells": bucket["cells"],
        }
        for field in FIELDS:
            mean = bucket[field] / bucket["area"]
            variance = max(0.0, bucket[f"{field}_sq"] / bucket["area"] - mean * mean)
            row[field] = mean
            row[f"{field}_std"] = math.sqrt(variance)
        rows.append(row)
    return rows, {"cells": len(states), "bins_requested": bins, "bins_nonempty": len(rows), "xlo": xlo, "xhi": xhi}


def front_candidates(rows: list[dict[str, float]], field: str, maximum: int = 8):
    candidates = []
    for index in range(1, len(rows) - 1):
        dx = rows[index + 1]["x"] - rows[index - 1]["x"]
        if dx <= 0.0:
            continue
        gradient = abs((rows[index + 1][field] - rows[index - 1][field]) / dx)
        candidates.append((gradient, index))
    selected = []
    for gradient, index in sorted(candidates, reverse=True):
        if all(abs(index - old_index) > 3 for _, old_index in selected):
            selected.append((gradient, index))
        if len(selected) == maximum:
            break
    return [{"x": rows[index]["x"], "abs_central_gradient": gradient} for gradient, index in selected]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", type=Path, required=True)
    parser.add_argument("--vtu", type=Path, required=True)
    parser.add_argument("--gamma", type=float, required=True)
    parser.add_argument("--bins", type=int, default=200)
    parser.add_argument("--csv", type=Path, required=True)
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()
    rows, summary = project(args.mesh, args.vtu, args.gamma, args.bins)
    args.csv.parent.mkdir(parents=True, exist_ok=True)
    args.summary.parent.mkdir(parents=True, exist_ok=True)
    columns = ("bin", "x", "area", "cells", *FIELDS, *(f"{field}_std" for field in FIELDS))
    with args.csv.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    summary["front_candidates"] = {
        field: front_candidates(rows, field) for field in ("rho", "pressure", "vy", "by")
    }
    summary["max_within_bin_std"] = {
        field: max(row[f"{field}_std"] for row in rows) for field in FIELDS
    }
    args.summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
