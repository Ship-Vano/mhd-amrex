#!/usr/bin/env python3
"""Run the corrected legacy HLLD flux on 1-D Riemann problems from the VKR.

Clones the immutable legacy source at the official commit, applies the
corrected-physics overlay with a check, compiles the bundled first-order 1-D
driver `scripts/legacy_1d_riemann.cpp` against the patched `MHDSolver1D.cpp`,
and runs Brio-Wu and Dai-Woodward. Emits a manifest with source / overlay /
output SHA-256. The official source directory is never modified.

This is the "N0" scheme (piecewise constant + forward Euler + HLLD), so a
comparison with the AMReX `hlld_flux` on the same grid isolates the flux
implementation and the mesh, not the order.
"""
from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import platform
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OFFICIAL_COMMIT = "9d0f60ea8576fac5d6f28c4dec142236d76131d6"
CASES = {"brio_wu": (0.1, 0.1), "dai_woodward": (0.2, 0.2)}   # (t_end, cfl)


def run(cmd, **kw):
    return subprocess.run(cmd, text=True, capture_output=True, check=False, **kw)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    h.update(path.read_bytes())
    return h.hexdigest()


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source", type=Path, required=True)
    ap.add_argument("--artifact-dir", type=Path, required=True)
    ap.add_argument("--overlay", type=Path,
                    default=ROOT / "legacy/patches/0001-legacy-corrected-physics.patch")
    ap.add_argument("--compiler", default="/opt/homebrew/opt/gcc/bin/g++-15")
    ap.add_argument("--nx", type=int, default=400)
    args = ap.parse_args()

    source = args.source.resolve()
    art = args.artifact_dir.resolve()
    if art.exists():
        raise SystemExit(f"refusing to overwrite {art}")
    head = run(["git", "-C", str(source), "rev-parse", "HEAD"]).stdout.strip()
    if head != OFFICIAL_COMMIT:
        raise SystemExit(f"source HEAD {head} != official {OFFICIAL_COMMIT}")
    if run(["git", "-C", str(source), "status", "--porcelain"]).stdout.strip():
        raise SystemExit("source worktree is dirty")

    art.mkdir(parents=True)
    wt = art / "source"
    errors: list[str] = []
    if run(["git", "clone", "-q", "--no-hardlinks", str(source), str(wt)]).returncode:
        raise SystemExit("clone failed")
    run(["git", "checkout", "-q", "--detach", OFFICIAL_COMMIT], cwd=wt)
    if run(["git", "apply", "--check", "--whitespace=error", str(args.overlay)], cwd=wt).returncode:
        raise SystemExit("overlay --check failed")
    run(["git", "apply", "--whitespace=error", str(args.overlay)], cwd=wt)

    driver_src = ROOT / "scripts/legacy_1d_riemann.cpp"
    binary = art / "legacy_1d_riemann"
    compile_cmd = [
        args.compiler, "-O2", "-std=c++23",
        "-I", str(wt / "src/MHDsolver"), "-I", str(wt / "src/geometry"),
        "-I", str(wt / "src/include"),
        str(driver_src),
        str(wt / "src/MHDsolver/MHDSolver1D.cpp"),
        str(wt / "src/include/utility/utility.cpp"),
        "-o", str(binary),
    ]
    cc = run(compile_cmd)
    (art / "compile.log").write_text(cc.stdout + cc.stderr)
    if cc.returncode:
        errors.append("compile_failed")

    results: dict[str, dict] = {}
    if not errors:
        for name, (t_end, cfl) in CASES.items():
            csv_path = art / f"{name}_1d.csv"
            r = run([str(binary), name, str(args.nx), str(cfl), str(csv_path)])
            (art / f"{name}.log").write_text(r.stdout + r.stderr)
            m = re.search(r"steps=(\d+) t=([\d.]+) rho_min=([\d.eE+-]+) "
                          r"p_min=([\d.eE+-]+) fallbacks=(\d+)", r.stdout)
            row = {"exit_code": r.returncode, "csv_sha256": sha256(csv_path)
                   if csv_path.is_file() else None}
            if m:
                row.update(steps=int(m.group(1)), reported_t=float(m.group(2)),
                           rho_min=float(m.group(3)), p_min=float(m.group(4)),
                           fallbacks=int(m.group(5)))
                row["finite_positive"] = row["rho_min"] > 0.0 and row["p_min"] > 0.0
            else:
                errors.append(f"{name}_parse_failed")
            results[name] = row

    manifest = {
        "schema_version": 1,
        "profile": "legacy_corrected",
        "kind": "1D Riemann, first-order HLLD (N0 scheme), frozen boundaries",
        "official_commit": OFFICIAL_COMMIT,
        "official_source_sha256": hashlib.sha256(
            subprocess.run(["git", "-C", str(source), "archive", OFFICIAL_COMMIT],
                           capture_output=True).stdout).hexdigest(),
        "overlay_sha256": sha256(args.overlay),
        "driver_sha256": sha256(driver_src),
        "nx": args.nx,
        "compiler": run([args.compiler, "--version"]).stdout.splitlines()[0],
        "platform": platform.platform(),
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "results": results,
        "quality_gate": {
            "status": "pass" if not errors and all(
                r.get("finite_positive") for r in results.values()) else "fail",
            "errors": errors,
        },
    }
    (art / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"artifact_dir": str(art),
                      "status": manifest["quality_gate"]["status"],
                      "results": results}, indent=2, sort_keys=True))
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
