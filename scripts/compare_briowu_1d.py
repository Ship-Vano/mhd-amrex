#!/usr/bin/env python3
"""Controlled 1-D Brio-Wu comparison: legacy_corrected vs AMReX scheme ablation.

Level A of the LEGACY_BASELINE comparison protocol. Every profile is on the
same physical strip (gamma=2, Bx=0.75, frozen x boundaries, t=0.1) and is
collapsed to a 1-D x-profile; metrics are taken against a shared provisional
reference (a grid-converged AMReX N3 run) interpolated onto each profile's
x-grid.

The reference is PROVISIONAL: the branch is fixed (compound / non-regular,
owner decision D-003) but an independent exact MHD Riemann solution is a T06
deliverable. Front positions and errors here are self-consistent scheme
comparisons, not literature-validated absolutes.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


FIELDS = ("rho", "u", "p", "by")


def load_amrex_strip(path: Path) -> dict:
    """standalone_verify CSV -> 1-D profile, averaging the y-strip per column."""
    xs: dict[float, dict[str, list[float]]] = {}
    with path.open() as handle:
        for row in csv.DictReader(handle):
            x = float(row["x"])
            bucket = xs.setdefault(round(x, 12), {k: [] for k in FIELDS})
            bucket["rho"].append(float(row["rho"]))
            bucket["u"].append(float(row["u"]))
            bucket["p"].append(float(row["p"]))
            bucket["by"].append(float(row["By"]))
    x_sorted = sorted(xs)
    prof = {"x": x_sorted}
    for key in FIELDS:
        prof[key] = [sum(xs[x][key]) / len(xs[x][key]) for x in x_sorted]
    prof["strip_std_max"] = {
        key: max(_std(xs[x][key]) for x in x_sorted) for key in FIELDS
    }
    return prof


def load_legacy_profile(path: Path) -> dict:
    """project_legacy_vtu.py area-weighted bin CSV -> the same 1-D profile shape."""
    rows = list(csv.DictReader(path.open()))
    prof = {
        "x": [float(r["x"]) for r in rows],
        "rho": [float(r["rho"]) for r in rows],
        "u": [float(r["vx"]) for r in rows],
        "p": [float(r["pressure"]) for r in rows],
        "by": [float(r["by"]) for r in rows],
    }
    prof["strip_std_max"] = {
        "rho": max(float(r["rho_std"]) for r in rows),
        "u": max(float(r["vx_std"]) for r in rows),
        "p": max(float(r["pressure_std"]) for r in rows),
        "by": max(float(r["by_std"]) for r in rows),
    }
    return prof


def _std(values: list[float]) -> float:
    n = len(values)
    if n < 2:
        return 0.0
    mean = sum(values) / n
    return math.sqrt(sum((v - mean) ** 2 for v in values) / n)


def interp(x_ref: list[float], y_ref: list[float], x: float) -> float:
    if x <= x_ref[0]:
        return y_ref[0]
    if x >= x_ref[-1]:
        return y_ref[-1]
    lo, hi = 0, len(x_ref) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if x_ref[mid] <= x:
            lo = mid
        else:
            hi = mid
    t = (x - x_ref[lo]) / (x_ref[hi] - x_ref[lo])
    return y_ref[lo] * (1.0 - t) + y_ref[hi] * t


def total_variation(values: list[float]) -> float:
    return sum(abs(b - a) for a, b in zip(values, values[1:]))


def metrics(profile: dict, reference: dict) -> dict:
    out: dict[str, dict[str, float]] = {}
    for key in FIELDS:
        xs = profile["x"]
        obs = profile[key]
        ref_on_x = [interp(reference["x"], reference[key], x) for x in xs]
        n = len(xs)
        l1 = sum(abs(o - r) for o, r in zip(obs, ref_on_x)) / n
        l2 = math.sqrt(sum((o - r) ** 2 for o, r in zip(obs, ref_on_x)) / n)
        linf = max(abs(o - r) for o, r in zip(obs, ref_on_x))
        tv_obs = total_variation(obs)
        tv_ref = total_variation(reference[key])
        ref_lo, ref_hi = min(reference[key]), max(reference[key])
        overshoot = max(0.0, max(obs) - ref_hi)
        undershoot = max(0.0, ref_lo - min(obs))
        out[key] = {
            "l1": l1, "l2": l2, "linf": linf,
            "tv": tv_obs, "tv_ref": tv_ref,
            "tv_excess": max(0.0, tv_obs - tv_ref),
            "overshoot": overshoot, "undershoot": undershoot,
            "strip_std_max": profile["strip_std_max"][key],
        }
    return out


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, required=True,
                        help="AMReX standalone CSV used as the provisional reference")
    parser.add_argument("--amrex", action="append", default=[],
                        metavar="TAG:CSV", help="an AMReX ablation profile")
    parser.add_argument("--legacy", action="append", default=[],
                        metavar="TAG:CSV", help="a legacy area-weighted bin profile")
    parser.add_argument("--summary", type=Path, required=True)
    args = parser.parse_args()

    reference = load_amrex_strip(args.reference)
    results: dict[str, dict] = {}

    for spec in args.amrex:
        tag, path = spec.split(":", 1)
        results[tag] = {"solver": "amrex", "source": path,
                        "metrics": metrics(load_amrex_strip(Path(path)), reference)}
    for spec in args.legacy:
        tag, path = spec.split(":", 1)
        results[tag] = {"solver": "legacy_corrected", "source": path,
                        "metrics": metrics(load_legacy_profile(Path(path)), reference)}

    summary = {
        "case": "brio_wu_1d",
        "gamma": 2.0, "Bx": 0.75, "t_end": 0.1, "cfl": 0.1,
        "reference": {
            "source": str(args.reference),
            "kind": "provisional grid-converged AMReX N3 (MUSCL-MC + SSP-RK2 + "
                    "Gardiner-Stone CT); compound/non-regular branch (D-003); "
                    "NOT an independent exact Riemann solution (T06)",
            "points": len(reference["x"]),
        },
        "profiles": results,
    }
    args.summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    # Human-readable table on stdout.
    hdr = f"{'profile':<10} {'solver':<16} " + " ".join(f"{f+'.L1':>10}" for f in FIELDS) + \
          f" {'by.TVxs':>9} {'by.under':>9} {'p.under':>9}"
    print(hdr)
    print("-" * len(hdr))
    for tag, row in results.items():
        m = row["metrics"]
        print(f"{tag:<10} {row['solver']:<16} " +
              " ".join(f"{m[f]['l1']:>10.3e}" for f in FIELDS) +
              f" {m['by']['tv_excess']:>9.3e} {m['by']['undershoot']:>9.3e}"
              f" {m['p']['undershoot']:>9.3e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
