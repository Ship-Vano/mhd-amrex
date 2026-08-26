#!/usr/bin/env python3
"""Create and run an isolated, immutable ``legacy_vkr`` controlled case.

The source directory is checked before use and never modified.  A detached
clone is built in the chosen raw-artifact directory; that directory contains
the mesh, JSON config, console log, source/archive checksums, outputs and a
machine-readable manifest.  Existing artifact directories are refused.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import platform
import re
import shutil
import subprocess
import sys
from pathlib import Path

from legacy_vkr_mesh import write_rectangular_tri_mesh


OFFICIAL_COMMIT = "9d0f60ea8576fac5d6f28c4dec142236d76131d6"
CASES = {
    "brio_wu": {"taskType": 1, "finalTime": 0.1, "gamma": 2.0, "domain": (0.0, 1.0, -0.01, 0.01), "resolution": (128, 4)},
    "cp_alfven": {"taskType": 8, "finalTime": 1.0, "gamma": 5.0 / 3.0, "domain": (0.0, 2.0 / 3.0**0.5, 0.0, 2.0), "resolution": (32, 56)},
    "magnetic_loop": {"taskType": 9, "finalTime": 2.0, "gamma": 5.0 / 3.0, "domain": (-1.0, 1.0, -0.5, 0.5), "resolution": (64, 32)},
}


def command(args: list[str], cwd: Path | None = None, capture: bool = False):
    return subprocess.run(args, cwd=cwd, check=True, text=True, capture_output=capture)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def tool_version(args: list[str]) -> str:
    try:
        return command(args, capture=True).stdout.strip() or command(args, capture=True).stderr.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True, help="clean MHD2D source at the official commit")
    parser.add_argument("--case", choices=CASES, required=True)
    parser.add_argument("--artifact-dir", type=Path, required=True)
    parser.add_argument("--compiler", default=shutil.which("g++-15") or "/opt/homebrew/opt/gcc/bin/g++-15")
    parser.add_argument("--jobs", type=int, default=4)
    args = parser.parse_args()
    source, artifact = args.source.resolve(), args.artifact_dir.resolve()
    if artifact.exists():
        raise SystemExit(f"refusing to overwrite existing artifact directory: {artifact}")
    head = command(["git", "-C", str(source), "rev-parse", "HEAD"], capture=True).stdout.strip()
    if head != OFFICIAL_COMMIT:
        raise SystemExit(f"source HEAD {head} is not official {OFFICIAL_COMMIT}")
    if command(["git", "-C", str(source), "status", "--porcelain"], capture=True).stdout:
        raise SystemExit("source worktree is dirty; refusing to use it")
    if not Path(args.compiler).exists():
        raise SystemExit(f"C++ compiler not found: {args.compiler}")

    artifact.mkdir(parents=True)
    worktree, build = artifact / "source", artifact / "build"
    command(["git", "clone", "--quiet", "--no-hardlinks", str(source), str(worktree)])
    command(["git", "checkout", "--quiet", "--detach", OFFICIAL_COMMIT], cwd=worktree)
    case = CASES[args.case]
    input_dir, output_dir = worktree / "InputData", worktree / "OutputData"
    input_dir.mkdir()
    output_dir.mkdir()
    xlo, xhi, ylo, yhi = case["domain"]
    nx, ny = case["resolution"]
    mesh = input_dir / "mesh.txt"
    write_rectangular_tri_mesh(mesh, xlo, xhi, ylo, yhi, nx, ny)
    config = {
        "taskType": case["taskType"],
        "finalTime": case["finalTime"],
        "debugDivergence": False,
        "ghostOutput": False,
        "iterationsPerFrame": 1000000000,
        "cylindrical": False,
        "gpu": False,
        "exportFileName": "OutputData/final.vtu",
    }
    config_path = input_dir / "solverConfig.json"
    config_path.write_text(json.dumps(config, indent=2, sort_keys=True) + "\n")
    command(["cmake", "-S", str(worktree), "-B", str(build), "-DCMAKE_BUILD_TYPE=Release", f"-DCMAKE_CXX_COMPILER={args.compiler}"])
    command(["cmake", "--build", str(build), "--parallel", str(args.jobs)])
    executable = build / "MHD2D"
    log_path = artifact / "console.log"
    with log_path.open("w", encoding="utf-8") as log:
        result = subprocess.run([str(executable)], cwd=worktree, text=True, stdout=log, stderr=subprocess.STDOUT)
    if result.returncode:
        raise SystemExit(f"solver exited {result.returncode}; inspect {log_path}")
    for filename in ("final.vtu", "tmpres_0.vtu"):
        shutil.copy2(output_dir / filename, artifact / filename)
    final_summary = artifact / "final_summary.json"
    command([sys.executable, str(Path(__file__).with_name("analyze_legacy_vtu.py")), "--mesh", str(mesh), "--vtu", str(artifact / "final.vtu"), "--case", args.case, "--gamma", str(case["gamma"]), "--output", str(final_summary)])
    final_diagnostics = json.loads(final_summary.read_text(encoding="utf-8"))
    log_text = log_path.read_text(encoding="utf-8")
    match = re.search(r"Final time = ([^;]+); iterations = (\d+)", log_text)
    actual_time = float(match.group(1)) if match else None
    iterations = int(match.group(2)) if match else None
    reached_target = actual_time is not None and abs(actual_time - case["finalTime"]) <= 1.0e-12
    final_state_finite = bool(final_diagnostics["finite"])
    archive = artifact / "legacy_source.tar"
    with archive.open("wb") as out:
        subprocess.run(["git", "-C", str(source), "archive", "--format=tar", OFFICIAL_COMMIT], check=True, stdout=out)
    manifest = {
        "schema_version": 1,
        "profile": "legacy_vkr",
        "source": {"official_commit": OFFICIAL_COMMIT, "archive_sha256": sha256(archive), "source_dirty": False},
        "case": {"name": args.case, "task_type": case["taskType"], "configured_final_time": case["finalTime"], "gamma": case["gamma"], "domain": {"x": [xlo, xhi], "y": [ylo, yhi]}, "mesh": {"generator": "scripts/legacy_vkr_mesh.py", "nx": nx, "ny": ny, "sha256": sha256(mesh)}},
        "build": {"compiler": str(Path(args.compiler).resolve()), "compiler_version": tool_version([args.compiler, "--version"]), "cmake_version": tool_version(["cmake", "--version"]), "build_type": "Release", "openmp": True},
        "machine": {"platform": platform.platform(), "python": sys.version, "omp_num_threads_environment": os.environ.get("OMP_NUM_THREADS", "unset")},
        "solver": {"process_exit_code": result.returncode, "reported_final_time": actual_time, "iterations": iterations, "reached_configured_final_time": reached_target, "final_state_finite": final_state_finite},
        "files": {name: sha256(artifact / name) for name in ("console.log", "final.vtu", "tmpres_0.vtu", "final_summary.json", "legacy_source.tar")},
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
    }
    (artifact / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"artifact_dir": str(artifact), "manifest": str(artifact / "manifest.json"), "final": final_diagnostics, "reached_target": reached_target}, indent=2, sort_keys=True))
    if result.returncode:
        raise SystemExit(f"solver exited {result.returncode}; inspect {log_path}")
    if not reached_target:
        raise SystemExit(f"solver did not reach configured final time {case['finalTime']}; inspect {log_path}")
    if not final_state_finite:
        raise SystemExit(f"solver reached final time with a non-finite state; inspect {log_path}")


if __name__ == "__main__":
    main()
