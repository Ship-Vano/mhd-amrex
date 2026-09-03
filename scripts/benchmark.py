#!/usr/bin/env python3
"""Phase T09: a timing record that can be defended, not a single stopwatch read.

A benchmark number is only meaningful with its provenance and its spread, so
every run here carries:

  * the machine, compiler, build type, git commit and whether the tree was dirty;
  * >= 1 discarded warm-up plus >= 5 timed repeats;
  * median and MAD (and min/max), not a mean -- one descheduled run must not
    move the reported number;
  * compute-only and end-to-end separated, because plotfile I/O is not part of
    the scheme's cost;
  * the case configuration and the step count, so cost per cell-step is
    recomputable.

Timings from a single laptop are explicitly diagnostic: they are recorded as
such and must not be quoted as scaling results (see docs/REPRODUCIBILITY.md).
"""
from __future__ import annotations

import argparse
import json
import platform
import re
import statistics
import subprocess
import sys
import time
from pathlib import Path


def command_text(cmd: list[str]) -> str:
    try:
        r = subprocess.run(cmd, text=True, capture_output=True)
        return r.stdout.strip() or r.stderr.strip()
    except OSError as exc:
        return f"<unavailable: {exc}>"


def machine_manifest(executable: Path) -> dict:
    root = Path(__file__).resolve().parents[1]
    return {
        "platform": platform.platform(),
        "machine": platform.machine(),
        "processor": platform.processor() or command_text(
            ["sysctl", "-n", "machdep.cpu.brand_string"]),
        "python": platform.python_version(),
        "cpu_count": __import__("os").cpu_count(),
        "executable": str(executable),
        "executable_mtime": time.strftime(
            "%Y-%m-%dT%H:%M:%SZ", time.gmtime(executable.stat().st_mtime)),
        "git_commit": command_text(["git", "-C", str(root), "rev-parse", "HEAD"]),
        "git_dirty": bool(command_text(["git", "-C", str(root), "status", "--porcelain"])),
        "omp_num_threads": __import__("os").environ.get("OMP_NUM_THREADS", "<unset>"),
    }


def run_once(executable: Path, config: Path) -> tuple[float, dict]:
    start = time.perf_counter()
    r = subprocess.run([str(executable), str(config)], text=True,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    wall = time.perf_counter() - start
    if r.returncode != 0:
        print(r.stdout[-3000:], file=sys.stderr)
        raise SystemExit(f"run failed ({r.returncode}) for {config}")
    info: dict = {}
    m = re.search(r"Evolve finished: (\d+) steps, t=([0-9.eE+-]+)", r.stdout)
    if m:
        info["steps"] = int(m.group(1))
        info["t_end"] = float(m.group(2))
    for key in ("hlld_fallbacks", "floor_events", "nonpositive_cells"):
        mm = re.search(rf"{key}=(\d+)", r.stdout)
        if mm:
            info[key] = int(mm.group(1))
    return wall, info


def summarize(samples: list[float]) -> dict:
    med = statistics.median(samples)
    mad = statistics.median([abs(s - med) for s in samples])
    return {
        "median_s": med,
        "mad_s": mad,
        "min_s": min(samples),
        "max_s": max(samples),
        "relative_spread": (max(samples) - min(samples)) / med if med > 0 else None,
        "samples_s": samples,
    }


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--executable", required=True, type=Path)
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--label", required=True)
    ap.add_argument("--repeats", type=int, default=5)
    ap.add_argument("--warmup", type=int, default=1)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()

    if args.repeats < 5:
        print("refusing to report fewer than 5 timed repeats: the spread would "
              "not be estimable", file=sys.stderr)
        return 2

    config = json.loads(args.config.read_text())
    cells = config["geometry"]["n_cell"]

    for _ in range(args.warmup):
        run_once(args.executable, args.config)      # discarded: cache/page warm-up

    samples, info = [], {}
    for _ in range(args.repeats):
        wall, info = run_once(args.executable, args.config)
        samples.append(wall)

    timing = summarize(samples)
    steps = info.get("steps")
    base_cells = cells[0] * cells[1]
    record = {
        "schema_version": 1,
        "label": args.label,
        "kind": "compute_only" if not config.get("output", {}).get(
            "write_plotfiles", True) else "end_to_end",
        "config": {
            "path": str(args.config),
            "problem": config.get("problem"),
            "n_cell": cells,
            "max_level": config.get("amr", {}).get("max_level", 0),
            "reflux": config.get("amr", {}).get("reflux", True),
            "cfl": config.get("time", {}).get("cfl"),
            "t_max": config.get("time", {}).get("t_max"),
            "integrator": config.get("time", {}).get("integrator"),
        },
        "run": info,
        "timing": timing,
        "derived": {
            "base_cells": base_cells,
            # Normalized by the BASE grid only. For an AMR run the refined
            # levels carry additional cells, so this is NOT a work-normalized
            # cost and must not be compared across different max_level values;
            # it is comparable between runs on the same hierarchy.
            "us_per_base_cell_step": (timing["median_s"] * 1e6 / (base_cells * steps))
                                     if steps else None,
            "base_cell_normalization_note":
                ("base grid only; refined levels add cells, so this is not "
                 "work-normalized across different max_level"
                 if config.get("amr", {}).get("max_level", 0) > 0
                 else "uniform grid: this is the work-normalized cost"),
        },
        "warmup_discarded": args.warmup,
        "repeats": args.repeats,
        "machine": machine_manifest(args.executable),
        "caveat": "single workstation, no pinning, no exclusive access: "
                  "diagnostic timing, not a scaling result",
    }

    text = json.dumps(record, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text)
    print(f"{args.label}: median {timing['median_s']:.4f} s  "
          f"MAD {timing['mad_s']:.4f} s  spread {timing['relative_spread']:.1%}  "
          f"steps={steps}  "
          f"{record['derived']['us_per_base_cell_step']:.4f} us/cell-step"
          if steps else f"{args.label}: median {timing['median_s']:.4f} s")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
