#!/usr/bin/env python3
"""Like-for-like 1-D Riemann comparison: legacy_corrected HLLD vs AMReX kernel.

Both the legacy driver (scripts/legacy_1d_riemann.cpp) and the AMReX N0 mode
are piecewise-constant + forward Euler + HLLD on the same grid, so N0 vs
legacy isolates the flux implementation and the mesh. N3 (MUSCL + SSP-RK2 +
Gardiner-Stone) shows what the second-order scheme adds.

Cases: brio_wu (gamma=2, Bx=0.75, t=0.1) and dai_woodward (gamma=5/3,
Bx=4/sqrt(4pi), t=0.2) from the VKR (Shamanov 2025, ch. 2). The reference is a
grid-converged AMReX N3 run -- provisional, not an independent exact solution.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

FIELDS = ("rho", "u", "p", "By")


def _load_amrex(path: Path, x_shift: float) -> dict:
    """standalone briowu1d/dw1d CSV -> 1-D profile (y-strip mean, per x)."""
    cols: dict[float, dict] = {}
    for r in csv.DictReader(path.open()):
        x = round(float(r["x"]) + x_shift, 10)
        b = cols.setdefault(x, {k: 0.0 for k in FIELDS} | {"n": 0})
        b["rho"] += float(r["rho"]); b["u"] += float(r["u"])
        b["p"] += float(r["p"]);     b["By"] += float(r["By"])
        b["n"] += 1
    xs = sorted(cols)
    return {"x": xs, **{k: [cols[x][k] / cols[x]["n"] for x in xs] for k in FIELDS}}


def _load_legacy(path: Path) -> dict:
    rows = list(csv.DictReader(path.open()))
    return {"x": [float(r["x"]) for r in rows],
            "rho": [float(r["rho"]) for r in rows],
            "u": [float(r["u"]) for r in rows],
            "p": [float(r["p"]) for r in rows],
            "By": [float(r["By"]) for r in rows]}


def _interp(xr, yr, x):
    if x <= xr[0]:
        return yr[0]
    if x >= xr[-1]:
        return yr[-1]
    lo, hi = 0, len(xr) - 1
    while hi - lo > 1:
        mid = (lo + hi) // 2
        (lo, hi) = (mid, hi) if xr[mid] <= x else (lo, mid)
    t = (x - xr[lo]) / (xr[hi] - xr[lo])
    return yr[lo] * (1 - t) + yr[hi] * t


def _tv(v):
    return sum(abs(b - a) for a, b in zip(v, v[1:]))


def _metrics(prof, ref):
    out = {}
    for k in FIELDS:
        xs, obs = prof["x"], prof[k]
        ro = [_interp(ref["x"], ref[k], x) for x in xs]
        n = len(xs)
        out[k] = {
            "l1": sum(abs(a - b) for a, b in zip(obs, ro)) / n,
            "l2": math.sqrt(sum((a - b) ** 2 for a, b in zip(obs, ro)) / n),
            "linf": max(abs(a - b) for a, b in zip(obs, ro)),
            "tv": _tv(obs), "tv_ref": _tv(ref[k]),
            "tv_excess": max(0.0, _tv(obs) - _tv(ref[k])),
            "overshoot": max(0.0, max(obs) - max(ref[k])),
            "undershoot": max(0.0, min(ref[k]) - min(obs)),
        }
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--case", choices=("brio_wu", "dai_woodward"), required=True)
    ap.add_argument("--reference", type=Path, required=True, help="AMReX N3 fine CSV")
    ap.add_argument("--legacy", type=Path, required=True, help="legacy_1d_riemann CSV")
    ap.add_argument("--amrex-n0", type=Path, required=True)
    ap.add_argument("--amrex-n3", type=Path, required=True)
    ap.add_argument("--amrex-x-shift", type=float, default=0.0,
                    help="add to AMReX x (briowu1d uses [0,1], legacy uses [-0.5,0.5] -> -0.5)")
    ap.add_argument("--summary", type=Path, required=True)
    args = ap.parse_args()

    ref = _load_amrex(args.reference, args.amrex_x_shift)
    profiles = {
        "legacy_corrected": _load_legacy(args.legacy),
        "amrex_N0": _load_amrex(args.amrex_n0, args.amrex_x_shift),
        "amrex_N3": _load_amrex(args.amrex_n3, args.amrex_x_shift),
    }
    summary = {
        "case": args.case,
        "reference": {"source": str(args.reference),
                      "kind": "provisional grid-converged AMReX N3 "
                              "(MUSCL-MC + SSP-RK2 + Gardiner-Stone); not an "
                              "independent exact Riemann solution",
                      "points": len(ref["x"])},
        "profiles": {name: {"points": len(p["x"]), "metrics": _metrics(p, ref)}
                     for name, p in profiles.items()},
    }
    args.summary.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    hdr = f"{'profile':<18} " + " ".join(f"{f+'.L1':>10}" for f in FIELDS) + \
          f" {'By.TVxs':>9} {'p.under':>9}"
    print(f"# {args.case}")
    print(hdr)
    print("-" * len(hdr))
    for name, row in summary["profiles"].items():
        m = row["metrics"]
        print(f"{name:<18} " + " ".join(f"{m[f]['l1']:>10.3e}" for f in FIELDS) +
              f" {m['By']['tv_excess']:>9.3e} {m['p']['undershoot']:>9.3e}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
