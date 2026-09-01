#!/usr/bin/env python3
"""Turn standalone diagnostic output into bounded CTest regressions."""
import argparse
import re
import subprocess
import sys


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", required=True)
    parser.add_argument("--case", choices=("briowu", "alfven", "loop", "dw1d"), required=True)
    parser.add_argument("--resolution", type=int)
    args = parser.parse_args()
    command = [args.executable, args.case]
    if args.case == "dw1d":
        command += [str(args.resolution or 128), "none", "euler", "bs", "0.2"]
    elif args.case == "loop":
        command += [str(args.resolution or 64), "0.5", "0.2"]   # half a wrap — fast
    elif args.resolution is not None:
        command.append(str(args.resolution))
    result = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print(result.stdout, end="")
    if result.returncode:
        return result.returncode
    div_values = [float(value) for value in re.findall(r"max\|divB\|(?: over run)?\s*=\s*([0-9.eE+-]+)", result.stdout)]
    if not div_values or max(div_values) > 1.0e-10:
        print("verification regression: divergence diagnostic missing or too large", file=sys.stderr)
        return 1

    # The positivity-preserving HLLD->HLL fallback must never fire on a smooth
    # canonical case (T04 gate). A non-zero count means a floored / degenerate
    # intermediate state occurred where the scheme should have stayed regular.
    fallbacks = re.search(r"hlld_fallbacks=(\d+)", result.stdout)
    if fallbacks is None:
        print("verification regression: hlld_fallbacks diagnostic missing", file=sys.stderr)
        return 1
    if int(fallbacks.group(1)) != 0:
        print(f"verification regression: {fallbacks.group(1)} HLLD fallbacks on a smooth case", file=sys.stderr)
        return 1
    if args.case == "alfven":
        match = re.search(r"L1\(Bperp\)=([0-9.eE+-]+)", result.stdout)
        if match is None or float(match.group(1)) >= 1.0e-2:
            print("verification regression: Alfvén L1(Bperp) missing or too large", file=sys.stderr)
            return 1
    if args.case == "loop":
        # Second-order field-loop transport must keep most of the magnetic
        # energy (the VKR first-order solver keeps ~18 % over two wraps).
        match = re.search(r"ratio=([0-9.eE+-]+)", result.stdout)
        if match is None or not (0.5 < float(match.group(1)) <= 1.0 + 1e-9):
            print("verification regression: field-loop E_B ratio out of range", file=sys.stderr)
            return 1
    if args.case == "dw1d":
        # First-order HLLD on the Dai-Woodward Riemann problem: finite and
        # positive. div B == 0 and fallbacks == 0 are checked above.
        m = re.search(r"rho_min=([0-9.eE+-]+) p_min=([0-9.eE+-]+)", result.stdout)
        if m is None:
            print("verification regression: dw1d summary missing", file=sys.stderr)
            return 1
        rho_min, p_min = float(m.group(1)), float(m.group(2))
        if not (rho_min > 0.0 and p_min > 0.0):
            print(f"verification regression: dw1d nonpositive state "
                  f"rho_min={rho_min} p_min={p_min}", file=sys.stderr)
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
