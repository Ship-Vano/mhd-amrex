#!/usr/bin/env python3
"""NEW-004: the ablation axes must be reachable from the JSON config *and* must
actually change the result.

CMP-006 wants first order / MUSCL / MUSCL+SSP-RK2 / full CT selectable without
recompiling. All four axes already exist as config keys -- scheme.limiter
(none => piecewise constant, i.e. first order in space), time.integrator, and
scheme.emf_averaging. What was missing is proof that they take effect: a knob
that parses and is then ignored is worse than no knob, because it makes an
ablation table meaningless.

So this test runs the four canonical profiles and requires that consecutive
ones differ. It deliberately compares the resulting state extrema rather than
just checking the run succeeds.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

MODES = [
    ("N0 constant + Euler + Balsara-Spicer",    "none",   "euler", "balsara_spicer"),
    ("N1 MUSCL(minmod) + Euler + BS",           "minmod", "euler", "balsara_spicer"),
    ("N2 MUSCL(mc) + SSP-RK2 + BS",             "mc",     "rk2",   "balsara_spicer"),
    ("N3 MUSCL(mc) + SSP-RK2 + Gardiner-Stone", "mc",     "rk2",   "gardiner_stone"),
]

BASE = {
    "problem": "orszag_tang",
    "geometry": {"prob_lo": [0.0, 0.0], "prob_hi": [1.0, 1.0], "n_cell": [64, 64]},
    "amr": {"max_level": 0},
    "bc": {"x_lo": "periodic", "x_hi": "periodic",
           "y_lo": "periodic", "y_hi": "periodic"},
    "time": {"cfl": 0.4, "t_max": 0.1, "integrator": "rk2"},
    "scheme": {"gamma": 1.6666666666666667, "limiter": "mc",
               "emf_averaging": "gardiner_stone"},
    "output": {"prefix": "plt_abl", "output_dir": ".", "format": "native",
               "diag_int": 100000, "write_plotfiles": False},
}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--executable", required=True)
    args = ap.parse_args()

    results = []
    with tempfile.TemporaryDirectory() as tmp:
        cfg = Path(tmp) / "ablation.json"
        for label, limiter, integrator, emf in MODES:
            c = json.loads(json.dumps(BASE))
            c["scheme"]["limiter"] = limiter
            c["scheme"]["emf_averaging"] = emf
            c["time"]["integrator"] = integrator
            c["output"]["output_dir"] = tmp
            cfg.write_text(json.dumps(c, indent=2))
            run = subprocess.run([args.executable, str(cfg)], text=True,
                                 stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            if run.returncode != 0:
                print(run.stdout[-2000:], file=sys.stderr)
                print(f"regression: mode '{label}' failed to run", file=sys.stderr)
                return 1
            m = re.search(r"rho_min=([0-9.eE+-]+) rho_max=([0-9.eE+-]+) "
                          r"p_min=([0-9.eE+-]+) p_max=([0-9.eE+-]+)", run.stdout)
            if m is None:
                print(f"regression: mode '{label}' printed no ranges", file=sys.stderr)
                return 1
            vals = tuple(float(m.group(i)) for i in range(1, 5))
            results.append((label, vals))
            print(f"  {label:<42} rho [{vals[0]:.6f}, {vals[1]:.6f}] "
                  f"p [{vals[2]:.6f}, {vals[3]:.6f}]")

    failures = []
    for (la, va), (lb, vb) in zip(results, results[1:]):
        if max(abs(a - b) for a, b in zip(va, vb)) < 1.0e-10:
            failures.append(f"'{la}' and '{lb}' give identical results: "
                            f"an ablation knob is being ignored")
    if failures:
        for f in failures:
            print(f"regression: {f}", file=sys.stderr)
        return 1

    print("all four ablation profiles are selectable from JSON and distinct")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
