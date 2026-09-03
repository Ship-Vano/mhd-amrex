#!/usr/bin/env python3
"""Gate T07: consistency of gas fluxes across coarse-fine boundaries (reflux).

The AMR hierarchy is advanced without subcycling, so a coarse cell next to a
refined patch is updated with the *coarse* flux while the fine cells behind
that same face are updated with the *fine* fluxes.  Unless the mismatch is
refluxed, the hierarchy is not conservative.

On a fully periodic problem nothing leaves the domain, so the hierarchy sums
of rho and of total energy must be constant to roundoff.  This test asserts
two things, and both matter:

  1. with reflux enabled the drift is at roundoff;
  2. with reflux disabled the drift is orders of magnitude larger.

Point 2 is what keeps the test honest: without it, a reflux implementation
that silently did nothing would still pass.
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path


FIELDS = ("rho_rel_drift", "mx_vel_drift", "my_vel_drift", "ene_rel_drift")


def run(executable: str, config: Path) -> dict:
    result = subprocess.run([executable, str(config)], text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if result.returncode != 0:
        print(result.stdout)
        raise SystemExit(f"solver exited {result.returncode} for {config}")
    line = None
    for candidate in result.stdout.splitlines():
        if candidate.startswith("conservation:"):
            line = candidate
    if line is None:
        print(result.stdout)
        raise SystemExit("solver printed no conservation diagnostic")
    out = {}
    for key in FIELDS + ("levels", "reflux"):
        match = re.search(rf"{key}=([0-9.eE+-]+)", line)
        if match is None:
            raise SystemExit(f"conservation diagnostic missing {key}: {line}")
        out[key] = float(match.group(1))
    coverage = [int(m) for m in re.findall(r"covered_by_finer=(\d+)", result.stdout)]
    base = re.search(r"level 0: cells=(\d+)", result.stdout)
    out["covered"] = coverage[0] if coverage else 0
    out["base_cells"] = int(base.group(1)) if base else 0
    return out


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", required=True)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--tol", type=float, default=1.0e-12,
                        help="max relative drift with reflux enabled")
    parser.add_argument("--min-defect", type=float, default=1.0e-6,
                        help="min drift that must appear when reflux is disabled")
    args = parser.parse_args()

    on = run(args.executable, args.config)

    # The comparison is only meaningful if a coarse-fine boundary exists at all:
    # a fine level covering the whole periodic domain has no interface, and then
    # reflux is a no-op that would pass a naive test vacuously.
    if on["levels"] < 2:
        print(f"regression: expected a refined hierarchy, got levels={on['levels']}",
              file=sys.stderr)
        return 1
    if not (0 < on["covered"] < on["base_cells"]):
        print(f"regression: fine level must partially cover the base grid "
              f"(covered={on['covered']} of {on['base_cells']}); "
              f"no coarse-fine interface means reflux is untested", file=sys.stderr)
        return 1

    for key in ("rho_rel_drift", "ene_rel_drift"):
        if not (on[key] <= args.tol):
            print(f"regression: {key}={on[key]:.6g} exceeds {args.tol:.6g} with reflux on",
                  file=sys.stderr)
            return 1

    # Same configuration, reflux disabled: the defect must be visible.
    config = json.loads(args.config.read_text())
    config["amr"]["reflux"] = False
    with tempfile.TemporaryDirectory() as tmp:
        disabled = Path(tmp) / "no_reflux.json"
        disabled.write_text(json.dumps(config, indent=2) + "\n")
        off = run(args.executable, disabled)

    defect = max(off["rho_rel_drift"], off["ene_rel_drift"])
    if not (defect >= args.min_defect):
        print(f"regression: disabling reflux changed nothing (drift={defect:.6g}); "
              f"the test is not exercising the coarse-fine correction", file=sys.stderr)
        return 1

    print(f"reflux on : " + "  ".join(f"{k}={on[k]:.4g}" for k in FIELDS))
    print(f"reflux off: " + "  ".join(f"{k}={off[k]:.4g}" for k in FIELDS))
    print(f"coarse-fine interface present: {on['covered']} of {on['base_cells']} "
          f"base cells covered")
    print(f"mass drift improved {off['rho_rel_drift'] / max(on['rho_rel_drift'], 1e-300):.3g}x, "
          f"energy drift improved {off['ene_rel_drift'] / max(on['ene_rel_drift'], 1e-300):.3g}x")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
