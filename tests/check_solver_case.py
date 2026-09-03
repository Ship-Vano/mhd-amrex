#!/usr/bin/env python3
"""Run one solver configuration and assert its physical invariants.

Covers the canonical-case matrix required by PROGRAM_PLAN (constant state,
Orszag-Tang, rotor, MHD blast, magnetic loop). The invariants asserted here are
the ones that hold regardless of resolution:

  * the run reaches its configured end time;
  * no cell needs a rho/p floor at any point (`nonpositive_cells`);
  * the normalized div B stays at roundoff (NUM-004);
  * rho and p stay inside a stated band -- for the canonical problems these are
    the literature colour-scale ranges, so a scheme regression that shifts a
    shock or damps a feature is caught, not just a crash.

`--exact-constant` additionally demands rho_min == rho_max and p_min == p_max
bitwise: a constant state must be preserved exactly, and any drift means the
reconstruction, the Riemann solve or the CT update is inconsistent.
"""
from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path


def parse(pattern: str, text: str, name: str) -> float:
    match = re.search(pattern, text)
    if match is None:
        raise SystemExit(f"solver printed no {name} diagnostic")
    return float(match.group(1))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--executable", required=True)
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--rho-min", type=float)
    ap.add_argument("--rho-max", type=float)
    ap.add_argument("--p-min", type=float)
    ap.add_argument("--p-max", type=float)
    ap.add_argument("--tol", type=float, default=0.05,
                    help="relative slack on the rho/p band")
    ap.add_argument("--max-divb", type=float, default=1.0e-12)
    ap.add_argument("--max-fallbacks", type=int, default=0)
    ap.add_argument("--max-floors", type=int, default=0,
                    help="NEW-003: rho/p floor events inside cons_to_prim, counted "
                         "separately from HLLD fallbacks and from post-step "
                         "nonpositive cells")
    ap.add_argument("--exact-constant", action="store_true")
    args = ap.parse_args()

    run = subprocess.run([args.executable, str(args.config)], text=True,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print(run.stdout[-2500:])
    if run.returncode != 0:
        return run.returncode

    if "Evolve finished" not in run.stdout:
        print("regression: run did not finish", file=sys.stderr)
        return 1

    fallbacks = int(parse(r"hlld_fallbacks=(\d+)", run.stdout, "hlld_fallbacks"))
    floors = int(parse(r"floor_events=(\d+)", run.stdout, "floor_events"))
    nonpos = int(parse(r"nonpositive_cells=(\d+)", run.stdout, "nonpositive_cells"))
    divb = parse(r"normalized=([0-9.eE+-]+)", run.stdout, "divb")
    rho_lo = parse(r"rho_min=([0-9.eE+-]+)", run.stdout, "rho_min")
    rho_hi = parse(r"rho_max=([0-9.eE+-]+)", run.stdout, "rho_max")
    p_lo = parse(r"p_min=([0-9.eE+-]+)", run.stdout, "p_min")
    p_hi = parse(r"p_max=([0-9.eE+-]+)", run.stdout, "p_max")

    failures = []
    if nonpos != 0:
        failures.append(f"{nonpos} cells at or below the rho/p floor")
    if fallbacks > args.max_fallbacks:
        failures.append(f"{fallbacks} HLLD->HLL fallbacks exceeds {args.max_fallbacks}")
    if floors > args.max_floors:
        failures.append(f"{floors} rho/p floor events exceeds {args.max_floors}")
    if not (divb <= args.max_divb):
        failures.append(f"normalized div B {divb:.3e} exceeds {args.max_divb:.3e}")
    if not (rho_lo > 0.0 and p_lo > 0.0):
        failures.append(f"nonpositive state: rho_min={rho_lo:.6g} p_min={p_lo:.6g}")

    if args.exact_constant:
        # Bitwise, not "within a tolerance": a uniform state is an exact
        # solution of the discrete scheme, so any spread is a real defect.
        if rho_lo != rho_hi or p_lo != p_hi:
            failures.append(f"constant state not preserved exactly: "
                            f"rho [{rho_lo!r}, {rho_hi!r}], p [{p_lo!r}, {p_hi!r}]")

    for value, bound, name, side in ((rho_lo, args.rho_min, "rho_min", "lo"),
                                     (rho_hi, args.rho_max, "rho_max", "hi"),
                                     (p_lo, args.p_min, "p_min", "lo"),
                                     (p_hi, args.p_max, "p_max", "hi")):
        if bound is None:
            continue
        slack = args.tol * max(abs(bound), 1.0e-30)
        if abs(value - bound) > slack:
            failures.append(f"{name}={value:.6g} outside {bound:.6g} +- {slack:.3g}")

    if failures:
        for f in failures:
            print(f"regression: {f}", file=sys.stderr)
        return 1

    print(f"OK  rho [{rho_lo:.6g}, {rho_hi:.6g}]  p [{p_lo:.6g}, {p_hi:.6g}]  "
          f"divB(norm)={divb:.3e}  fallbacks={fallbacks}  floors={floors}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
