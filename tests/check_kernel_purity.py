#!/usr/bin/env python3
"""ADR 0001: src/kernels/ must not depend structurally on AMReX.

The kernel layer is the one part a future DG operator reuses unchanged, so it
may depend on AMReX only for the scalar type and the device qualifier macro.
Anything that pulls in AMReX's data model (MultiFab, Box, Array4, MFIter,
Geometry, ...) would tie the physics to the finite-volume mesh backend and
silently break the boundary the ADR establishes.

This is checkable rather than aspirational: the kernel unit tests already link
without AMReX, but a header could acquire a container include and still compile
inside the solver. Here we check the includes directly.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

ALLOWED_AMREX = {"AMReX_REAL.H", "AMReX_GpuQualifiers.H"}
FORBIDDEN_SYMBOLS = ("MultiFab", "amrex::Box", "Array4", "MFIter", "Geometry",
                     "BoxArray", "DistributionMapping", "FabArray")


def strip_comments(text: str) -> str:
    """Remove // and /* */ comments.

    Prose mentioning a type is not a dependency: MhdState.H documents that the
    conserved variables live in a cell-centred MultiFab, which is true and worth
    saying. Only code counts.
    """
    text = re.sub(r"/\*.*?\*/", " ", text, flags=re.S)
    return re.sub(r"//[^\n]*", " ", text)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--kernel-dir", required=True, type=Path)
    args = ap.parse_args()

    headers = sorted(args.kernel_dir.glob("*.H"))
    if not headers:
        print(f"no kernel headers found under {args.kernel_dir}", file=sys.stderr)
        return 1

    failures = []
    for header in headers:
        text = strip_comments(header.read_text())
        for include in re.findall(r'#\s*include\s*[<"]([^>"]+)[>"]', text):
            if include.startswith("AMReX") and include not in ALLOWED_AMREX:
                failures.append(f"{header.name}: includes {include}; ADR 0001 allows "
                                f"only {sorted(ALLOWED_AMREX)}")
        for symbol in FORBIDDEN_SYMBOLS:
            if symbol in text:
                failures.append(f"{header.name}: mentions {symbol}; the kernel layer "
                                f"must not use AMReX containers")

    if failures:
        for f in failures:
            print(f"ADR 0001 violation: {f}", file=sys.stderr)
        return 1

    print(f"kernel layer clean: {len(headers)} headers, "
          f"AMReX includes limited to {sorted(ALLOWED_AMREX)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
