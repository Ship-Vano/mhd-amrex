#!/usr/bin/env python3
"""Run the same configuration on CPU and CUDA binaries and compare diagnostics.

Not every diagnostic may be compared the same way, and getting this wrong makes
the gate look strict while checking nothing.

Counters (fallbacks, floors, nonpositive cells) are integers describing which
branches the scheme took. They must match exactly: a GPU run that takes a
different number of positivity fallbacks is running a different scheme.

State ranges (rho, p) are compared with a tolerance. Reduction order, FMA
contraction and thread count all differ between the two builds, so the last
bits will differ and then amplify over the run; the tolerance has to admit that
while still catching a real divergence.

max|divB| must NOT be compared CPU-to-GPU with a tight tolerance, and the
default tolerances would make that comparison vacuous anyway. It is a maximum
over differences of nearly equal face values: both operands are order one, the
difference is order 1e-16, so a one-ulp change moves it by tens of percent.
Measured on this project between two CPU builds differing only in operation
order, it moved 8.2e-13 -> 4.2e-12 with the solution otherwise identical. The
meaningful statement is that each build keeps divergence at roundoff, so both
are checked against an absolute threshold instead of against each other.
"""
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
    ap.add_argument("--rtol", type=float, default=5e-10,
                    help="допуск для rho/p. Для длинного нелинейного прогона его "
                         "может потребоваться ослабить -- но тогда ослабление "
                         "должно быть записано в манифест кампании, а не "
                         "подобрано молча до прохождения gate")
    ap.add_argument("--divb-max", type=float, default=1.0e-12,
                    help="порог нормированной div B, тот же, что у CPU-тестов; "
                         "проверяется у каждой сборки отдельно, а не CPU против GPU")
    ap.add_argument("--output-dir", type=Path,
                    help="durable directory for cpu.log and gpu.log")
    args = ap.parse_args()

    cpu = diagnostics(run(args.cpu, args.config, "cpu", args.output_dir), "cpu")
    gpu = diagnostics(run(args.gpu, args.config, "gpu", args.output_dir), "gpu")
    failures: list[str] = []
    for name in ("fallbacks", "floors", "nonpositive"):
        if cpu[name] != gpu[name]:
            failures.append(f"{name}: CPU={cpu[name]} GPU={gpu[name]}")
    for name in ("rho_min", "rho_max", "p_min", "p_max"):
        if not close(float(cpu[name]), float(gpu[name]), args.atol, args.rtol):
            failures.append(f"{name}: CPU={cpu[name]:.17g} GPU={gpu[name]:.17g}")
    # div B -- по абсолютному порогу у каждой сборки (см. docstring).
    for label, values in (("CPU", cpu), ("GPU", gpu)):
        if not (float(values["divb_normalized"]) <= args.divb_max):
            failures.append(f"{label} divb_normalized={values['divb_normalized']:.3e} "
                            f"exceeds {args.divb_max:.1e}")
    if failures:
        print("CPU/GPU parity failed (this is a numerical correctness gate):")
        print("\n".join(f"  {f}" for f in failures))
        return 1
    print(f"OK CPU/GPU parity: counters identical; rho/p agree within "
          f"atol={args.atol:g}, rtol={args.rtol:g}; "
          f"divb_normalized CPU={cpu['divb_normalized']:.3e} "
          f"GPU={gpu['divb_normalized']:.3e}, both <= {args.divb_max:.1e}")
    print(f"   (max|divB| CPU={cpu['divb_max_abs']:.3e} GPU={gpu['divb_max_abs']:.3e} "
          f"-- reported, not gated: see module docstring)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
