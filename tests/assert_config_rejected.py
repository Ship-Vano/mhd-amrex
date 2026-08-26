#!/usr/bin/env python3
"""Assert that a deliberately malformed config is rejected before a run."""
import argparse
import subprocess
import sys


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--executable", required=True)
    parser.add_argument("--config", required=True)
    args = parser.parse_args()
    result = subprocess.run([args.executable, args.config], text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    print(result.stdout, end="")
    if result.returncode == 0 or "unknown key" not in result.stdout:
        print("config regression: malformed schema was not rejected as expected", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
