#!/usr/bin/env python3
"""Gate T06: cross-validate two independent schemes on Brio-Wu.

The reference (Kurganov-Tadmor central, tests/briowu_reference.cpp) shares no
code with the solver under test and performs no wave decomposition, so
agreement between them is real evidence rather than self-consistency. The
Brio-Wu solution contains a compound wave (D-003), which is exactly why an
exact Riemann solver cannot serve as the reference.

Two things are asserted:
  1. the reference converges under refinement (an unconverged reference is
     worthless, and this catches a broken build of it);
  2. the HLLD solver agrees with it to a stated relative L1, and its state
     extrema match -- so a scheme regression that shifts the compound wave or
     clips a plateau fails here even if it stays smooth and positive.
"""
from __future__ import annotations

import argparse
import csv
import subprocess
import sys
import tempfile
from pathlib import Path

VARS = ("rho", "u", "p", "By")


def load(path: Path, nx: int | None = None) -> list[dict]:
    rows = list(csv.DictReader(path.open()))
    if nx:
        acc: dict[int, dict] = {}
        for r in rows:
            k = int(round(float(r["x"]) * nx - 0.5))
            b = acc.setdefault(k, {v: 0.0 for v in ("x",) + VARS} | {"n": 0})
            for key in ("x",) + VARS:
                b[key] += float(r[key])
            b["n"] += 1
        return [{k: b[k] / b["n"] for k in ("x",) + VARS}
                for b in (acc[i] for i in sorted(acc))]
    return [{k: float(r[k]) for k in ("x",) + VARS} for r in rows]


def interp(ref, xs, var):
    rx = [r["x"] for r in ref]; rv = [r[var] for r in ref]; out = []
    for x in xs:
        if x <= rx[0]: out.append(rv[0]); continue
        if x >= rx[-1]: out.append(rv[-1]); continue
        lo, hi = 0, len(rx) - 1
        while hi - lo > 1:
            m = (lo + hi) // 2
            if rx[m] <= x: lo = m
            else: hi = m
        t = (x - rx[lo]) / (rx[hi] - rx[lo])
        out.append(rv[lo] * (1 - t) + rv[hi] * t)
    return out


def rel_l1(cand, ref_rows, var):
    xs = [r["x"] for r in cand]
    rv = interp(ref_rows, xs, var)
    cv = [r[var] for r in cand]
    scale = max(abs(v) for v in rv) or 1.0
    return sum(abs(cv[i] - rv[i]) for i in range(len(xs))) / len(xs) / scale


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference", required=True)
    ap.add_argument("--verify", required=True)
    ap.add_argument("--fronts", required=True)
    ap.add_argument("--coarse", type=int, default=400)
    ap.add_argument("--fine", type=int, default=1600)
    ap.add_argument("--candidate-nx", type=int, default=800)
    ap.add_argument("--max-rel-l1", type=float, default=6.0e-3)
    args = ap.parse_args()

    with tempfile.TemporaryDirectory() as tmp:
        d = Path(tmp)
        for n in (args.coarse, args.fine):
            r = subprocess.run([args.reference, str(n), "0.1", str(d / f"ref{n}.csv"), "0.4"],
                               text=True, capture_output=True)
            print(r.stdout.strip())
            if r.returncode:
                print(r.stderr, file=sys.stderr); return r.returncode
        r = subprocess.run([args.verify, "briowu1d", str(args.candidate_nx),
                            "mc", "rk2", "gs", "0.1", str(d / "hlld.csv")],
                           text=True, capture_output=True)
        if r.returncode:
            print(r.stdout + r.stderr, file=sys.stderr); return r.returncode

        coarse = load(d / f"ref{args.coarse}.csv")
        fine = load(d / f"ref{args.fine}.csv")
        hlld = load(d / "hlld.csv", args.candidate_nx)

        # 1. the reference must converge under refinement
        self_err = sum(rel_l1(coarse, fine, v) for v in VARS) / len(VARS)
        if not (self_err < 2.0e-2):
            print(f"regression: reference not converging, self rel L1 = {self_err:.4e}",
                  file=sys.stderr)
            return 1

        # 2. the independent schemes must agree
        failures = []
        for v in VARS:
            e = rel_l1(hlld, fine, v)
            print(f"  {v:4s} rel_L1(HLLD vs KT) = {e:.4e}")
            if not (e <= args.max_rel_l1):
                failures.append(f"{v}: rel L1 {e:.4e} exceeds {args.max_rel_l1:.4e}")

        for v, tol in (("rho", 0.02), ("p", 0.02)):
            a = min(r[v] for r in hlld); b = min(r[v] for r in fine)
            if abs(a - b) > tol * max(abs(b), 1e-30):
                failures.append(f"min {v}: HLLD {a:.6g} vs KT {b:.6g}")

        if failures:
            for f in failures:
                print(f"regression: {f}", file=sys.stderr)
            return 1

    print(f"reference self-convergence {self_err:.3e}; two independent schemes agree")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
