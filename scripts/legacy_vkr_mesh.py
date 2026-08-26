#!/usr/bin/env python3
"""Emit the minimal triangular-mesh text format accepted by legacy MHD2D.

This is a deterministic controlled-input generator, not a reconstruction of
the unproven historical Netgen mesh.  It emits two CCW triangles per rectangle
and is intentionally dependency-free so its output can be hashed and rerun.
"""

from __future__ import annotations

import argparse
from pathlib import Path


def write_rectangular_tri_mesh(
    output: Path, xlo: float, xhi: float, ylo: float, yhi: float, nx: int, ny: int
) -> None:
    if nx < 1 or ny < 1 or not xhi > xlo or not yhi > ylo:
        raise ValueError("require xhi > xlo, yhi > ylo, nx >= 1, ny >= 1")

    dx = (xhi - xlo) / nx
    dy = (yhi - ylo) / ny
    node_id = lambda i, j: j * (nx + 1) + i + 1
    lines = ["$Nodes", str((nx + 1) * (ny + 1))]
    for j in range(ny + 1):
        for i in range(nx + 1):
            lines.append(f"{node_id(i, j)} {xlo + i * dx:.17g} {ylo + j * dy:.17g} 0")
    lines.extend(["$EndNodes", "$Elements", str(2 * nx * ny)])
    element_id = 1
    for j in range(ny):
        for i in range(nx):
            n00, n10 = node_id(i, j), node_id(i + 1, j)
            n01, n11 = node_id(i, j + 1), node_id(i + 1, j + 1)
            # Both triples are counter-clockwise in the xy plane.
            lines.append(f"{element_id} 3 {n00} {n10} {n11}")
            element_id += 1
            lines.append(f"{element_id} 3 {n00} {n11} {n01}")
            element_id += 1
    lines.append("$EndElements")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--xlo", type=float, required=True)
    parser.add_argument("--xhi", type=float, required=True)
    parser.add_argument("--ylo", type=float, required=True)
    parser.add_argument("--yhi", type=float, required=True)
    parser.add_argument("--nx", type=int, required=True)
    parser.add_argument("--ny", type=int, required=True)
    args = parser.parse_args()
    write_rectangular_tri_mesh(**vars(args))


if __name__ == "__main__":
    main()
