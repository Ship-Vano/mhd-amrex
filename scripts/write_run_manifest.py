#!/usr/bin/env python3
"""Write a deterministic provenance manifest for one solver configuration.

The manifest deliberately has no wall-clock field: equal inputs on the same
source tree and machine must serialize identically. Timing records belong in
benchmarks/raw/ and are introduced by the benchmark phase.
"""
import argparse
import hashlib
import json
import os
import pathlib
import platform
import subprocess
import sys


ROOT = pathlib.Path(__file__).resolve().parents[1]


def command_output(command: list[str]) -> str:
    try:
        return subprocess.check_output(command, cwd=ROOT, text=True, stderr=subprocess.DEVNULL).strip()
    except (OSError, subprocess.CalledProcessError):
        return "unavailable"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, help="input JSON configuration")
    parser.add_argument("--output", required=True, help="target manifest JSON")
    args = parser.parse_args()
    config_path = pathlib.Path(args.config).resolve()
    output_path = pathlib.Path(args.output).resolve()
    config_bytes = config_path.read_bytes()
    manifest = {
        "schema_version": 1,
        "solver": "mhd2d_amrex",
        "source": {
            "commit": command_output(["git", "rev-parse", "HEAD"]),
            "dirty": command_output(["git", "status", "--porcelain=v1"]) != "",
        },
        "configuration": {
            "path": os.path.relpath(config_path, ROOT),
            "sha256": hashlib.sha256(config_bytes).hexdigest(),
            "json": json.loads(config_bytes),
        },
        "environment": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "processor": platform.processor() or "unavailable",
            "python": platform.python_version(),
            "omp_num_threads": os.environ.get("OMP_NUM_THREADS", "unset"),
            "mpi_launcher": command_output(["mpirun", "--version"]).split("\n")[0],
        },
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
