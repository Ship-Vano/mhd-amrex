#!/usr/bin/env python3
"""Phase T06: annotate Brio-Wu fronts and score a solution against an
independent reference.

The reference must NOT be the scheme under test. `tests/briowu_reference.cpp`
provides one: a Kurganov-Tadmor central scheme that shares no code with the
HLLD path and does not decompose the solution into waves at all. That matters
here because the Brio-Wu solution contains a *compound* wave (slow shock with an
attached rarefaction, decision D-003), which classical exact Riemann solvers
cannot represent -- so "exact solver" is not an option and an independent
*scheme* is.

Metrics per variable: L1 / L2 / Linf against the reference, total variation and
its excess over the reference, overshoot/undershoot beyond the reference range,
and the number of sign changes of the first difference (a proxy for spurious
oscillation). Per front: position error in cells and the 10-90 % width, so
diffusion and misplacement are reported separately rather than mixed into one
error number.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

VARIABLES = ("rho", "u", "p", "By")


def load(path: Path, nx: int | None = None) -> list[dict]:
    rows = list(csv.DictReader(path.open()))
    if not rows:
        raise SystemExit(f"{path}: empty")
    # A 2-D strip dump repeats each x for every y; collapse to a 1-D profile.
    if nx:
        acc: dict[int, dict] = {}
        for r in rows:
            k = int(round(float(r["x"]) * nx - 0.5))
            b = acc.setdefault(k, {v: 0.0 for v in ("x",) + VARIABLES} | {"n": 0})
            for key in ("x",) + VARIABLES:
                b[key] += float(r[key])
            b["n"] += 1
        return [{k: b[k] / b["n"] for k in ("x",) + VARIABLES}
                for b in (acc[i] for i in sorted(acc))]
    return [{k: float(r[k]) for k in ("x",) + VARIABLES} for r in rows]


def interp(ref: list[dict], xs: list[float], var: str) -> list[float]:
    rx = [r["x"] for r in ref]
    rv = [r[var] for r in ref]
    out = []
    for x in xs:
        if x <= rx[0]:
            out.append(rv[0]); continue
        if x >= rx[-1]:
            out.append(rv[-1]); continue
        lo, hi = 0, len(rx) - 1
        while hi - lo > 1:
            mid = (lo + hi) // 2
            if rx[mid] <= x:
                lo = mid
            else:
                hi = mid
        t = (x - rx[lo]) / (rx[hi] - rx[lo])
        out.append(rv[lo] * (1.0 - t) + rv[hi] * t)
    return out


def total_variation(values: list[float]) -> float:
    return sum(abs(values[i + 1] - values[i]) for i in range(len(values) - 1))


def sign_changes(values: list[float], eps: float) -> int:
    diffs = [values[i + 1] - values[i] for i in range(len(values) - 1)]
    diffs = [d for d in diffs if abs(d) > eps]
    return sum(1 for i in range(len(diffs) - 1) if diffs[i] * diffs[i + 1] < 0.0)


def find_fronts(ref: list[dict], var: str, count: int) -> list[dict]:
    """Locate the `count` sharpest features of `var` in the reference.

    Fronts are taken where |d var / dx| is locally largest, with a minimum
    separation so a single front is not reported twice.
    """
    xs = [r["x"] for r in ref]
    vs = [r[var] for r in ref]
    grad = [abs(vs[i + 1] - vs[i]) for i in range(len(vs) - 1)]
    order = sorted(range(len(grad)), key=lambda i: grad[i], reverse=True)
    picked: list[int] = []
    min_sep = max(4, len(grad) // 60)
    for i in order:
        if all(abs(i - j) > min_sep for j in picked):
            picked.append(i)
        if len(picked) == count:
            break
    picked.sort()
    return [{"variable": var, "x": 0.5 * (xs[i] + xs[i + 1]),
             "jump": vs[i + 1] - vs[i]} for i in picked]


def front_width_cells(rows: list[dict], var: str, x0: float, jump: float,
                      window: float) -> float | None:
    """Cells spanned by the 10-90 % transition of a front near x0."""
    sel = [r for r in rows if abs(r["x"] - x0) <= window]
    if len(sel) < 3:
        return None
    lo_val, hi_val = sel[0][var], sel[-1][var]
    if abs(hi_val - lo_val) < 0.1 * abs(jump):
        return None
    lo10 = lo_val + 0.1 * (hi_val - lo_val)
    lo90 = lo_val + 0.9 * (hi_val - lo_val)
    xs = [r["x"] for r in sel]
    vals = [r[var] for r in sel]

    def cross(target: float) -> float | None:
        for i in range(len(vals) - 1):
            a, b = vals[i], vals[i + 1]
            if (a - target) * (b - target) <= 0.0 and a != b:
                t = (target - a) / (b - a)
                return xs[i] + t * (xs[i + 1] - xs[i])
        return None

    x10, x90 = cross(lo10), cross(lo90)
    if x10 is None or x90 is None:
        return None
    dx = (rows[1]["x"] - rows[0]["x"]) if len(rows) > 1 else 1.0
    return abs(x90 - x10) / dx


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--reference", required=True, type=Path)
    ap.add_argument("--candidate", required=True, type=Path)
    ap.add_argument("--candidate-nx", type=int,
                    help="collapse a 2-D strip dump with this many x cells")
    ap.add_argument("--reference-nx", type=int)
    ap.add_argument("--label", default="candidate")
    ap.add_argument("--candidate-x-shift", type=float, default=0.0,
                    help="added to the candidate x column; the legacy driver uses "
                         "[-0.5, 0.5] while the reference uses [0, 1]")
    ap.add_argument("--annotations", type=Path,
                    help="write/read the versioned front annotation JSON")
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()

    ref = load(args.reference, args.reference_nx)
    cand = load(args.candidate, args.candidate_nx)
    if args.candidate_x_shift:
        for r in cand:
            r["x"] += args.candidate_x_shift
    xs = [r["x"] for r in cand]
    dx = xs[1] - xs[0] if len(xs) > 1 else 1.0

    # Brio-Wu at t=0.1 has five features left to right: fast rarefaction,
    # compound wave, contact discontinuity, slow shock, fast rarefaction.
    # rho shows the contact; By shows the compound wave and the slow shock.
    fronts = find_fronts(ref, "rho", 3) + find_fronts(ref, "By", 3)
    if args.annotations:
        if args.annotations.is_file():
            fronts = json.loads(args.annotations.read_text())["fronts"]
        else:
            args.annotations.write_text(json.dumps(
                {"schema": 1, "case": "brio_wu", "t": 0.1,
                 "reference": str(args.reference), "fronts": fronts},
                indent=2) + "\n")

    result = {"label": args.label, "reference": str(args.reference),
              "candidate": str(args.candidate), "cells": len(xs),
              "variables": {}, "fronts": []}

    for var in VARIABLES:
        rv = interp(ref, xs, var)
        cv = [r[var] for r in cand]
        n = len(xs)
        diff = [cv[i] - rv[i] for i in range(n)]
        scale = max(abs(v) for v in rv) or 1.0
        tv_c, tv_r = total_variation(cv), total_variation(rv)
        result["variables"][var] = {
            "l1": sum(abs(d) for d in diff) / n,
            "l2": math.sqrt(sum(d * d for d in diff) / n),
            "linf": max(abs(d) for d in diff),
            "relative_l1": sum(abs(d) for d in diff) / n / scale,
            "tv": tv_c,
            "tv_reference": tv_r,
            "tv_excess": tv_c - tv_r,
            "overshoot": max(0.0, max(cv) - max(rv)),
            "undershoot": max(0.0, min(rv) - min(cv)),
            "first_difference_sign_changes": sign_changes(cv, 1.0e-6 * scale),
        }

    for fr in fronts:
        var = fr["variable"]
        window = 25.0 * dx
        near = [r for r in cand if abs(r["x"] - fr["x"]) <= window]
        pos = None
        if len(near) > 2:
            grads = [(abs(near[i + 1][var] - near[i][var]),
                      0.5 * (near[i + 1]["x"] + near[i]["x"]))
                     for i in range(len(near) - 1)]
            pos = max(grads)[1]
        result["fronts"].append({
            "variable": var,
            "reference_x": fr["x"],
            "candidate_x": pos,
            "position_error_cells": (abs(pos - fr["x"]) / dx) if pos is not None else None,
            "width_cells_10_90": front_width_cells(cand, var, fr["x"], fr["jump"], window),
            "reference_jump": fr["jump"],
        })

    text = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(text)
    print(f"=== {args.label} ({len(xs)} cells) vs {args.reference.name} ===")
    for var, m in result["variables"].items():
        print(f"  {var:4s} relL1={m['relative_l1']:.4e}  Linf={m['linf']:.4e}  "
              f"TVexc={m['tv_excess']:+.4e}  over={m['overshoot']:.3e}  "
              f"under={m['undershoot']:.3e}  signchg={m['first_difference_sign_changes']}")
    for f in result["fronts"]:
        pe = f["position_error_cells"]
        w = f["width_cells_10_90"]
        print(f"  front {f['variable']:3s} x_ref={f['reference_x']:.4f}  "
              f"pos_err={'n/a' if pe is None else f'{pe:.2f}'} cells  "
              f"width={'n/a' if w is None else f'{w:.2f}'} cells")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
