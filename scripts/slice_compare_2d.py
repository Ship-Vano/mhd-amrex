#!/usr/bin/env python3
"""Extract a 1-D slice from a legacy VTU and/or an AMReX standalone CSV.

Purpose: put the canonical literature slices side by side --
  * Orszag-Tang: pressure / density along y = 0.3125
    (Avdeeva & Lukin 2019, Fig. 8; VKRB Shamanov, Figs. 24, 26)
  * rotating cylinder: density / pressure along the diagonal x = y
    (Avdeeva & Lukin 2019, Fig. 6; VKRB Shamanov, Figs. 20, 22)

The literature plots are qualitative (no digitised data), so this script
reports the profile plus feature statistics (extrema, spurious-oscillation
count, front sharpness) rather than an L-error against a picture.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def _read_legacy(mesh: Path, vtu: Path, gamma: float):
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from analyze_legacy_vtu import primitive, read_mesh, read_vtu
    nodes, elements = read_mesh(mesh)
    states = read_vtu(vtu)
    pts = []
    for element, state in zip(elements, states):
        (x0, y0), (x1, y1), (x2, y2) = (nodes[i] for i in element)
        cx, cy = (x0 + x1 + x2) / 3.0, (y0 + y1 + y2) / 3.0
        rho, p, vx, vy, vz, bx, by, bz = primitive(state, gamma)
        pts.append((cx, cy, rho, p, vx, vy, bx, by))
    return pts


def _read_amrex(csv_path: Path):
    pts = []
    for r in csv.DictReader(csv_path.open()):
        pts.append((float(r["x"]), float(r["y"]), float(r["rho"]), float(r["p"]),
                    float(r["u"]), float(r["v"]), float(r["Bx"]), float(r["By"])))
    return pts


def _slice(pts, mode: str, coord: float, half_band: float):
    picked = []
    for (x, y, rho, p, vx, vy, bx, by) in pts:
        if mode == "horizontal":
            if abs(y - coord) <= half_band:
                picked.append((x, rho, p, vx, vy, bx, by))
        else:  # diagonal x = y
            if abs(x - y) <= half_band * math.sqrt(2.0):
                picked.append((x, rho, p, vx, vy, bx, by))
    picked.sort()
    return picked


def _features(picked, key_index: int):
    """Feature stats for one field along the slice."""
    xs = [row[0] for row in picked]
    ys = [row[key_index] for row in picked]
    if len(ys) < 5:
        return {"points": len(ys)}
    # bin to a fixed regular grid to suppress the triangular-mesh scatter
    n = 64
    lo, hi = min(xs), max(xs)
    grid = [[] for _ in range(n)]
    for x, v in zip(xs, ys):
        grid[min(n - 1, int((x - lo) / (hi - lo) * n))].append(v)
    prof = [sum(b) / len(b) if b else math.nan for b in grid]
    prof = [v for v in prof if not math.isnan(v)]
    diffs = [b - a for a, b in zip(prof, prof[1:])]
    sign_changes = sum(1 for a, b in zip(diffs, diffs[1:]) if a * b < 0)
    tv = sum(abs(d) for d in diffs)
    return {
        "points": len(ys),
        "min": min(ys), "max": max(ys),
        "binned_min": min(prof), "binned_max": max(prof),
        "total_variation": tv,
        "interior_extrema": sign_changes,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--case", choices=("orszag_tang", "rotor"), required=True)
    ap.add_argument("--legacy-mesh", type=Path)
    ap.add_argument("--legacy-vtu", type=Path)
    ap.add_argument("--amrex-csv", type=Path)
    ap.add_argument("--gamma", type=float, default=5.0 / 3.0)
    ap.add_argument("--half-band", type=float, default=0.02)
    ap.add_argument("--out", type=Path, required=True)
    args = ap.parse_args()

    mode, coord = ("horizontal", 0.3125) if args.case == "orszag_tang" else ("diagonal", 0.0)
    fields = {"rho": 1, "pressure": 2}

    result = {"case": args.case, "slice": mode, "coord": coord,
              "reference_note": (
                  "Orszag-Tang y=0.3125 slice vs Avdeeva&Lukin 2019 Fig.8 / "
                  "VKRB Figs.24,26; rotor x=y diagonal vs Avdeeva Fig.6 / VKRB "
                  "Figs.20,22. Literature plots are qualitative."),
              "solvers": {}}

    if args.legacy_mesh and args.legacy_vtu:
        picked = _slice(_read_legacy(args.legacy_mesh, args.legacy_vtu, args.gamma),
                        mode, coord, args.half_band)
        result["solvers"]["legacy_corrected"] = {
            "slice_points": len(picked),
            **{f: _features(picked, i) for f, i in fields.items()},
            "profile": [{"x": round(r[0], 5), "rho": r[1], "p": r[2]} for r in picked],
        }
    if args.amrex_csv:
        picked = _slice(_read_amrex(args.amrex_csv), mode, coord, args.half_band)
        result["solvers"]["amrex"] = {
            "slice_points": len(picked),
            **{f: _features(picked, i) for f, i in fields.items()},
            "profile": [{"x": round(r[0], 5), "rho": r[1], "p": r[2]} for r in picked],
        }

    args.out.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    for solver, data in result["solvers"].items():
        print(f"{solver}: {data['slice_points']} slice points")
        for f in fields:
            s = data[f]
            print(f"  {f:9s} range [{s['min']:.4f}, {s['max']:.4f}] "
                  f"binned [{s['binned_min']:.4f}, {s['binned_max']:.4f}] "
                  f"TV={s['total_variation']:.3f} interior_extrema={s['interior_extrema']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
