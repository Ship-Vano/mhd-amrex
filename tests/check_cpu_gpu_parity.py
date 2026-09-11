#!/usr/bin/env python3
"""Run the same configuration on CPU and CUDA binaries and compare diagnostics."""
from __future__ import annotations

import argparse
import math
import re
import subprocess
from pathlib import Path


PATTERNS = {
    "rho_min": r"rho_min=([0-9.eE+-]+)",
    "rho_max": r"rho_max=([0-9.eE+-]+)",
    "p_min": r"p_min=([0-9.eE+-]+)",
    "p_max": r"p_max=([0-9.eE+-]+)",
    "divb_max_abs": r"max_abs=([0-9.eE+-]+)",
    "divb_normalized": r"normalized=([0-9.eE+-]+)",
    "fallbacks": r"hlld_fallbacks=(\d+)",
    "floors": r"floor_events=(\d+)",
    "nonpositive": r"nonpositive_cells=(\d+)",
}


def run(executable: str, config: Path, label: str, output_dir: Path | None) -> str:
    completed = subprocess.run([executable, str(config)], text=True,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if output_dir is not None:
        output_dir.mkdir(parents=True, exist_ok=True)
        (output_dir / f"{label}.log").write_text(completed.stdout, encoding="utf-8")
    if completed.returncode != 0:
        raise SystemExit(f"{label} solver failed with exit code {completed.returncode}\n{completed.stdout[-4000:]}")
    if "Evolve finished" not in completed.stdout:
        raise SystemExit(f"{label} solver did not reach Evolve finished")
    return completed.stdout


def diagnostics(text: str, label: str) -> dict[str, float | int]:
    result: dict[str, float | int] = {}
    for name, pattern in PATTERNS.items():
        matches = re.findall(pattern, text)
        if not matches:
            raise SystemExit(f"{label} solver printed no {name} diagnostic")
        value = matches[-1]
        result[name] = int(value) if name in {"fallbacks", "floors", "nonpositive"} else float(value)
    return result


def close(cpu: float, gpu: float, atol: float, rtol: float) -> bool:
    return math.isfinite(cpu) and math.isfinite(gpu) and abs(cpu - gpu) <= atol + rtol * max(abs(cpu), abs(gpu))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--cpu", required=True, help="CPU-only mhd2d executable")
    ap.add_argument("--gpu", required=True, help="CUDA mhd2d executable")
    ap.add_argument("--config", required=True, type=Path)
    ap.add_argument("--atol", type=float, default=5e-11)
    ap.add_argument("--rtol", type=float, default=5e-10)
    ap.add_argument("--output-dir", type=Path,
                    help="durable directory for cpu.log and gpu.log")
    args = ap.parse_args()

    cpu = diagnostics(run(args.cpu, args.config, "cpu", args.output_dir), "cpu")
    gpu = diagnostics(run(args.gpu, args.config, "gpu", args.output_dir), "gpu")
    failures: list[str] = []
    for name in ("fallbacks", "floors", "nonpositive"):
        if cpu[name] != gpu[name]:
            failures.append(f"{name}: CPU={cpu[name]} GPU={gpu[name]}")
    for name in ("rho_min", "rho_max", "p_min", "p_max", "divb_max_abs", "divb_normalized"):
        if not close(float(cpu[name]), float(gpu[name]), args.atol, args.rtol):
            failures.append(f"{name}: CPU={cpu[name]:.17g} GPU={gpu[name]:.17g}")
    if failures:
        print("CPU/GPU parity failed (this is a numerical correctness gate):")
        print("\n".join(f"  {f}" for f in failures))
        return 1
    print("OK CPU/GPU parity: diagnostics agree within "
          f"atol={args.atol:g}, rtol={args.rtol:g}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
