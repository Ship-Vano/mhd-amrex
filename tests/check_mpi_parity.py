#!/usr/bin/env python3
"""Mandatory matrix: MPI decomposition parity.

The scheme must not depend on how the domain is split across ranks. Ghost
exchange, the nodal EMF sync and the flux registers all touch rank boundaries,
so a decomposition-dependent bug there is both easy to introduce and invisible
on a single rank.

The check is deliberately strict: the hierarchy-wide rho/p extrema and the
conservation drift must match the serial run to 1e-12, not merely "look
similar". Skips cleanly (exit 0) if mpirun is unavailable, so the suite still
runs on a build without MPI.
"""
from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys


def extract(text: str) -> dict:
    out = {}
    m = re.search(r"rho_min=([0-9.eE+-]+) rho_max=([0-9.eE+-]+) "
                  r"p_min=([0-9.eE+-]+) p_max=([0-9.eE+-]+)", text)
    if m:
        out.update({k: float(m.group(i)) for i, k in
                    enumerate(("rho_min", "rho_max", "p_min", "p_max"), start=1)})
    m = re.search(r"rho_rel_drift=([0-9.eE+-]+)", text)
    if m:
        out["rho_rel_drift"] = float(m.group(1))
    m = re.search(r"normalized=([0-9.eE+-]+)", text)
    if m:
        out["divb"] = float(m.group(1))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--executable", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--ranks", default="2,4")
    ap.add_argument("--tol", type=float, default=1.0e-12)
    args = ap.parse_args()

    if shutil.which("mpirun") is None:
        print("mpirun not available; skipping MPI decomposition parity")
        return 0

    def run(cmd):
        r = subprocess.run(cmd, text=True, stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT)
        if r.returncode != 0:
            print(r.stdout[-2000:], file=sys.stderr)
            raise SystemExit(f"run failed: {' '.join(cmd)}")
        return extract(r.stdout)

    serial = run([args.executable, args.config])
    if not serial:
        print("regression: serial run printed no diagnostics", file=sys.stderr)
        return 1
    print(f"  serial: {serial}")

    failures = []
    for n in (int(x) for x in args.ranks.split(",")):
        got = run(["mpirun", "-np", str(n), "--oversubscribe",
                   args.executable, args.config])
        for key, ref in serial.items():
            cur = got.get(key)
            if cur is None:
                failures.append(f"np={n}: missing {key}")
            elif abs(cur - ref) > args.tol * max(abs(ref), 1.0):
                failures.append(f"np={n}: {key} {cur:.12g} != serial {ref:.12g}")
        print(f"  np={n}: {got}")

    if failures:
        for f in failures:
            print(f"regression: {f}", file=sys.stderr)
        return 1
    print("MPI decomposition parity holds")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
