#!/usr/bin/env python3
"""Generate an irregular Netgen triangular mesh for the legacy text format.

The legacy solver reads only nodes and triangle connectivity.  This adapter
uses Netgen for the actual unstructured triangulation, validates orientation
and periodic boundary-node pairing, then emits that minimal format without
claiming that it is an archived historical mesh.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def _quality(points: list[tuple[float, float]], triangles: list[tuple[int, int, int]]) -> dict[str, float]:
    min_angle = 180.0
    max_aspect = 0.0
    for triangle in triangles:
        coords = [points[index - 1] for index in triangle]
        lengths = [math.dist(coords[i], coords[(i + 1) % 3]) for i in range(3)]
        semiperimeter = sum(lengths) / 2.0
        area2 = abs((coords[1][0] - coords[0][0]) * (coords[2][1] - coords[0][1]) -
                    (coords[2][0] - coords[0][0]) * (coords[1][1] - coords[0][1]))
        area = area2 / 2.0
        if area <= 0.0:
            raise ValueError("Netgen returned a degenerate triangle")
        inradius = area / semiperimeter
        circumradius = lengths[0] * lengths[1] * lengths[2] / (4.0 * area)
        max_aspect = max(max_aspect, circumradius / (2.0 * inradius))
        for i, length in enumerate(lengths):
            other_a, other_b = lengths[(i + 1) % 3], lengths[(i + 2) % 3]
            cosine = max(-1.0, min(1.0, (other_a * other_a + other_b * other_b - length * length) / (2.0 * other_a * other_b)))
            min_angle = min(min_angle, math.degrees(math.acos(cosine)))
    return {"min_angle_degrees": min_angle, "max_radius_ratio": max_aspect}


def write_netgen_rectangular_tri_mesh(
    output: Path,
    xlo: float,
    xhi: float,
    ylo: float,
    yhi: float,
    maxh: float,
    grading: float = 0.3,
    boundary_jitter: float = 0.0,
) -> dict[str, object]:
    if not (xhi > xlo and yhi > ylo and maxh > 0.0 and 0.0 < grading <= 1.0 and 0.0 <= boundary_jitter < 0.25):
        raise ValueError("require positive extents, maxh > 0, 0 < grading <= 1, and 0 <= boundary_jitter < 0.25")
    try:
        from netgen import config
        from netgen.geom2d import SplineGeometry
    except ImportError as error:
        raise RuntimeError("Netgen Python bindings are required; use the dedicated Netgen environment") from error

    # Pre-segment every side so opposite periodic sides have the same nodes.
    # Netgen then creates an irregular interior triangulation without breaking
    # legacy's coordinate-based periodic pairing.
    # Boundary perturbations can lengthen a segment.  Over-resolve the side
    # slightly so Netgen never adaptively splits only one member of a periodic
    # pair and destroys the one-to-one boundary correspondence.
    jitter_safety = 1.0 + 2.0 * boundary_jitter
    nx = math.ceil((xhi - xlo) / maxh * jitter_safety)
    ny = math.ceil((yhi - ylo) / maxh * jitter_safety)

    def subdivide(lo: float, hi: float, count: int, phase: float) -> list[float]:
        spacing = (hi - lo) / count
        values = [lo]
        for index in range(1, count):
            # Deterministic perturbation, shared by opposite sides.  This
            # prevents a structured right-triangle control mesh from being
            # mistaken for the primary unstructured legacy verification mesh.
            offset = boundary_jitter * spacing * math.sin(2.0 * math.pi * 7.0 * index / count + phase)
            values.append(lo + index * spacing + offset)
        values.append(hi)
        if not all(second > first for first, second in zip(values, values[1:])):
            raise ValueError("jitter produced non-monotone boundary nodes")
        return values

    x_values = subdivide(xlo, xhi, nx, 0.0)
    y_values = subdivide(ylo, yhi, ny, math.pi / 5.0)
    boundary = [(x, ylo) for x in x_values]
    boundary += [(xhi, y) for y in y_values[1:]]
    boundary += [(x, yhi) for x in reversed(x_values[:-1])]
    boundary += [(xlo, y) for y in reversed(y_values[1:-1])]
    geometry = SplineGeometry()
    boundary_ids = [geometry.AppendPoint(*point) for point in boundary]
    for start, end in zip(range(len(boundary_ids)), list(range(1, len(boundary_ids))) + [0]):
        geometry.Append(["line", boundary_ids[start], boundary_ids[end]], bc="outer")
    mesh = geometry.GenerateMesh(maxh=maxh, grading=grading)
    points = [(float(point[0]), float(point[1])) for point in mesh.Points()]
    triangles: list[tuple[int, int, int]] = []
    for element in mesh.Elements2D():
        vertices = tuple(int(vertex.nr) for vertex in element.vertices)
        if len(vertices) != 3:
            raise ValueError(f"Netgen returned a non-triangle element: {vertices}")
        a, b, c = (points[index - 1] for index in vertices)
        signed_area2 = (b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])
        if signed_area2 < 0.0:
            vertices = (vertices[0], vertices[2], vertices[1])
        triangles.append(vertices)

    tolerance = 1.0e-12 * max(xhi - xlo, yhi - ylo)
    left = sorted(y for x, y in points if abs(x - xlo) <= tolerance)
    right = sorted(y for x, y in points if abs(x - xhi) <= tolerance)
    bottom = sorted(x for x, y in points if abs(y - ylo) <= tolerance)
    top = sorted(x for x, y in points if abs(y - yhi) <= tolerance)
    pairs_match = lambda first, second: len(first) == len(second) and all(abs(a - b) <= tolerance for a, b in zip(first, second))
    if not pairs_match(left, right) or not pairs_match(bottom, top):
        raise ValueError("Netgen boundary nodes are not pairwise periodic-compatible")

    lines = ["$Nodes", str(len(points))]
    lines.extend(f"{index} {x:.17g} {y:.17g} 0" for index, (x, y) in enumerate(points, start=1))
    lines.extend(["$EndNodes", "$Elements", str(len(triangles))])
    lines.extend(f"{index} 3 {a} {b} {c}" for index, (a, b, c) in enumerate(triangles, start=1))
    lines.append("$EndElements")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return {
        "generator": "Netgen",
        "netgen_version": config.NETGEN_VERSION,
        "nodes": len(points),
        "triangles": len(triangles),
        "maxh": maxh,
        "grading": grading,
        "boundary_jitter": boundary_jitter,
        "boundary_subdivisions": {"x": nx, "y": ny},
        "quality": _quality(points, triangles),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--metadata", type=Path)
    parser.add_argument("--xlo", type=float, required=True)
    parser.add_argument("--xhi", type=float, required=True)
    parser.add_argument("--ylo", type=float, required=True)
    parser.add_argument("--yhi", type=float, required=True)
    parser.add_argument("--maxh", type=float, required=True)
    parser.add_argument("--grading", type=float, default=0.3)
    parser.add_argument("--boundary-jitter", type=float, default=0.0,
                        help="deterministic fraction of local boundary spacing; use <= 0.24")
    args = parser.parse_args()
    metadata = write_netgen_rectangular_tri_mesh(
        args.output, args.xlo, args.xhi, args.ylo, args.yhi, args.maxh, args.grading, args.boundary_jitter
    )
    if args.metadata:
        args.metadata.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(metadata, sort_keys=True))


if __name__ == "__main__":
    main()
