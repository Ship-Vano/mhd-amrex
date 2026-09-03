#!/usr/bin/env python3
"""Phase T11 (single node): OpenMP and MPI scaling with honest normalisation.

Reports S_p = T_1 / T_p and E_p = S_p / p, each from median-of-repeats timings
with the spread retained, plus cells per rank. It records the CPU topology,
because on a heterogeneous machine (performance + efficiency cores) speedup is
expected to stop tracking thread count once the performance cores are full, and
a table without that context invites the wrong conclusion.

Accuracy is checked, not assumed: every configuration must reproduce the same
final rho/p ranges, so a "speedup" that silently changed the answer fails.
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import re
import statistics
import subprocess
import sys
import time
from pathlib import Path


def topology() -> dict:
    def sysctl(key):
        try:
            return int(subprocess.run(["sysctl", "-n", key], text=True,
                                      capture_output=True).stdout.strip())
        except (ValueError, OSError):
            return None
    return {
        "logical_cpus": os.cpu_count(),
        "physical_cpus": sysctl("hw.physicalcpu"),
        "performance_cores": sysctl("hw.perflevel0.physicalcpu"),
        "efficiency_cores": sysctl("hw.perflevel1.physicalcpu"),
        "platform": platform.platform(),
    }


def run(cmd: list[str], env: dict) -> tuple[float, dict]:
    start = time.perf_counter()
    r = subprocess.run(cmd, env=env, text=True,
                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    wall = time.perf_counter() - start
    if r.returncode != 0:
        print(r.stdout[-3000:], file=sys.stderr)
        raise SystemExit(f"run failed: {' '.join(cmd)}")
    m = re.search(r"rho_min=([0-9.eE+-]+) rho_max=([0-9.eE+-]+) "
                  r"p_min=([0-9.eE+-]+) p_max=([0-9.eE+-]+)", r.stdout)
    st = re.search(r"Evolve finished: (\d+) steps", r.stdout)
    return wall, {"ranges": tuple(float(m.group(i)) for i in range(1, 5)) if m else None,
                  "steps": int(st.group(1)) if st else None}


def measure(cmd: list[str], env: dict, repeats: int, warmup: int) -> tuple[dict, dict]:
    for _ in range(warmup):
        run(cmd, env)
    samples, info = [], {}
    for _ in range(repeats):
        w, info = run(cmd, env)
        samples.append(w)
    med = statistics.median(samples)
    mad = statistics.median([abs(s - med) for s in samples])
    return {"median_s": med, "mad_s": mad, "min_s": min(samples),
            "max_s": max(samples), "samples_s": samples}, info


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--executable", required=True)
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--mode", choices=("omp", "mpi"), required=True)
    ap.add_argument("--counts", default="1,2,4,8")
    ap.add_argument("--repeats", type=int, default=5)
    ap.add_argument("--warmup", type=int, default=1)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()

    counts = [int(c) for c in args.counts.split(",")]
    cfg = json.loads(args.config.read_text())
    cells = cfg["geometry"]["n_cell"]
    base_cells = cells[0] * cells[1]

    rows, baseline, reference_ranges = [], None, None
    for p in counts:
        env = dict(os.environ)
        if args.mode == "omp":
            env["OMP_NUM_THREADS"] = str(p)
            cmd = [args.executable, str(args.config)]
        else:
            env["OMP_NUM_THREADS"] = "1"
            cmd = ["mpirun", "-np", str(p), "--oversubscribe",
                   args.executable, str(args.config)]
        timing, info = measure(cmd, env, args.repeats, args.warmup)
        if baseline is None:
            baseline = timing["median_s"]
            reference_ranges = info["ranges"]
        speedup = baseline / timing["median_s"]
        # A speedup that changed the answer is not a speedup.
        same = (info["ranges"] is not None and reference_ranges is not None and
                max(abs(a - b) for a, b in zip(info["ranges"], reference_ranges)) < 1.0e-10)
        rows.append({
            "p": p, "timing": timing, "speedup": speedup, "efficiency": speedup / p,
            "cells_per_unit": base_cells / p,
            "result_identical_to_serial": same,
            "ranges": info["ranges"],
        })
        flag = "" if same else "   <-- RESULT DIFFERS FROM SERIAL"
        print(f"  {args.mode} p={p:<3d} median {timing['median_s']:.4f} s  "
              f"MAD {timing['mad_s']:.4f}  S_p={speedup:.2f}  E_p={speedup/p:.0%}{flag}")

    record = {
        "schema_version": 1,
        "mode": args.mode,
        "config": str(args.config),
        "n_cell": cells,
        "repeats": args.repeats,
        "warmup_discarded": args.warmup,
        "topology": topology(),
        "rows": rows,
        "caveat": "single workstation, no pinning, no exclusive access; on a "
                  "heterogeneous CPU speedup saturates once the performance "
                  "cores are full -- diagnostic, not a cluster scaling result",
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")

    bad = [r for r in rows if not r["result_identical_to_serial"]]
    if bad:
        print(f"regression: {len(bad)} configuration(s) changed the result",
              file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
