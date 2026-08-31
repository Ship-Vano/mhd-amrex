#!/usr/bin/env python3
"""Regression: circularly polarised Alfvén wave converges at ~2nd order.

The standalone driver evolves the Tóth 30-degree CP Alfvén wave for one full
return (t=1) and prints its area-averaged L1 error over (B_perp, Bz, v_perp,
vz). Refining the grid must reduce that error at the accepted rate (T04 gate:
observed order p >= 1.8) — this is what separates the MUSCL + SSP-RK2 scheme
from the first-order legacy baseline.
"""
import argparse
import math
import re
import subprocess
import sys


def l1_sum4(executable: str, n: int) -> float:
    out = subprocess.run([executable, "alfven", str(n)], text=True,
                         stdout=subprocess.PIPE, stderr=subprocess.STDOUT, check=True).stdout
    match = re.search(r"L1\(sum4\)=([0-9.eE+-]+)", out)
    if match is None:
        raise SystemExit(f"alfven {n}: L1(sum4) not found in:\n{out}")
    fallbacks = re.search(r"hlld_fallbacks=(\d+)", out)
    if fallbacks is None or int(fallbacks.group(1)) != 0:
        raise SystemExit(f"alfven {n}: unexpected HLLD fallback(s)\n{out}")
    return float(match.group(1))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", required=True)
    parser.add_argument("--coarse", type=int, default=64)
    parser.add_argument("--fine", type=int, default=128)
    parser.add_argument("--min-order", type=float, default=1.8)
    args = parser.parse_args()

    e_coarse = l1_sum4(args.executable, args.coarse)
    e_fine = l1_sum4(args.executable, args.fine)
    order = math.log(e_coarse / e_fine) / math.log(args.fine / args.coarse)
    print(f"alfven L1(sum4): N={args.coarse} {e_coarse:.6e}, "
          f"N={args.fine} {e_fine:.6e}, observed order = {order:.3f}")
    if order < args.min_order:
        print(f"regression: observed order {order:.3f} < required {args.min_order}",
              file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
