#!/usr/bin/env python3
"""Require manifests from identical inputs to be byte-equivalent."""
import argparse
import pathlib
import subprocess
import sys


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runner", required=True)
    parser.add_argument("--config", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()
    out_dir = pathlib.Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    first, second = out_dir / "first.json", out_dir / "second.json"
    for target in (first, second):
        result = subprocess.run([sys.executable, args.runner, "--config", args.config, "--output", str(target)], text=True)
        if result.returncode:
            return result.returncode
    if first.read_bytes() != second.read_bytes():
        print("manifest regression: identical inputs produced different manifests", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
