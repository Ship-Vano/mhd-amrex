#!/usr/bin/env python3
"""Turn standalone diagnostic output into bounded CTest regressions."""
import argparse
import re
import subprocess
import sys


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", required=True)
    parser.add_argument("--case", choices=("briowu", "alfven"), required=True)
    parser.add_argument("--resolution", type=int)
    args = parser.parse_args()
    command = [args.executable, args.case]
    if args.resolution is not None:
        command.append(str(args.resolution))
    result = subprocess.run(command, text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print(result.stdout, end="")
    if result.returncode:
        return result.returncode
    div_values = [float(value) for value in re.findall(r"max\|divB\|(?: over run)?\s*=\s*([0-9.eE+-]+)", result.stdout)]
    if not div_values or max(div_values) > 1.0e-10:
        print("verification regression: divergence diagnostic missing or too large", file=sys.stderr)
        return 1
    if args.case == "alfven":
        match = re.search(r"L1\(Bperp\)=([0-9.eE+-]+)", result.stdout)
        if match is None or float(match.group(1)) >= 1.0e-2:
            print("verification regression: Alfvén L1(Bperp) missing or too large", file=sys.stderr)
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
