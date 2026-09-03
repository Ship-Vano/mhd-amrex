#!/usr/bin/env python3
"""Render a 2-D pseudocolour map (ParaView-style) from a standalone CSV.

Pure standard library (no numpy / matplotlib / ParaView). Reads the
`x,y,<field>,...` grid written by `mhd2d_verify` and writes a PNG with a
`jet` colour map, in the layout of the VKR thesis field maps.

  field_map.py <csv> <field> <out.png> [--gamma G] [--upscale K] [--vmin a --vmax b]

Prints the actual data range so the report caption can cite it.
"""
from __future__ import annotations

import argparse
import csv
import math
import struct
import sys
import zlib
from pathlib import Path


def _jet(t: float) -> tuple[int, int, int]:
    t = 0.0 if t < 0 else (1.0 if t > 1 else t)
    def cl(v): return int(255 * (0.0 if v < 0 else (1.0 if v > 1 else v)))
    return cl(1.5 - abs(4 * t - 3)), cl(1.5 - abs(4 * t - 2)), cl(1.5 - abs(4 * t - 1))


def _write_png(path: Path, w: int, h: int, rows: list[bytes]) -> None:
    def chunk(typ: bytes, data: bytes) -> bytes:
        return (struct.pack(">I", len(data)) + typ + data +
                struct.pack(">I", zlib.crc32(typ + data) & 0xffffffff))
    ihdr = struct.pack(">IIBBBBB", w, h, 8, 2, 0, 0, 0)   # 8-bit truecolour
    raw = b"".join(b"\x00" + r for r in rows)
    path.write_bytes(b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", ihdr) +
                     chunk(b"IDAT", zlib.compress(raw, 9)) + chunk(b"IEND", b""))


def _derive(rows: list[dict], field: str, gamma: float):
    if field in ("rho", "p", "u", "v", "w", "Bx", "By", "Bz", "divB"):
        return [float(r[field]) for r in rows]
    if field == "Bmag":
        return [math.hypot(float(r["Bx"]), float(r["By"])) for r in rows]
    if field == "pmag":
        return [0.5 * (float(r["Bx"]) ** 2 + float(r["By"]) ** 2 + float(r["Bz"]) ** 2)
                for r in rows]
    if field == "vmag":
        return [math.sqrt(float(r["u"]) ** 2 + float(r["v"]) ** 2 + float(r["w"]) ** 2)
                for r in rows]
    raise SystemExit(f"unknown field {field}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("csv", type=Path)
    ap.add_argument("field")
    ap.add_argument("out", type=Path)
    ap.add_argument("--gamma", type=float, default=5.0 / 3.0)
    ap.add_argument("--upscale", type=int, default=4)
    ap.add_argument("--vmin", type=float)
    ap.add_argument("--vmax", type=float)
    args = ap.parse_args()

    rows = list(csv.DictReader(args.csv.open()))
    xs = sorted({float(r["x"]) for r in rows})
    ys = sorted({float(r["y"]) for r in rows})
    nx, ny = len(xs), len(ys)
    if nx * ny != len(rows):
        raise SystemExit(f"{args.csv}: {len(rows)} rows != {nx}x{ny} grid")
    xi = {v: i for i, v in enumerate(xs)}
    yi = {v: i for i, v in enumerate(ys)}

    vals = _derive(rows, args.field, args.gamma)
    grid = [[0.0] * nx for _ in range(ny)]
    for r, v in zip(rows, vals):
        grid[yi[float(r["y"])]][xi[float(r["x"])]] = v

    vmin = args.vmin if args.vmin is not None else min(vals)
    vmax = args.vmax if args.vmax is not None else max(vals)
    span = (vmax - vmin) or 1.0
    print(f"{args.csv.name} {args.field}: range [{min(vals):.6g}, {max(vals):.6g}] "
          f"grid {nx}x{ny}")

    k = max(1, args.upscale)
    W, H = nx * k, ny * k
    png_rows: list[bytes] = []
    for j in range(H):                        # PNG top row first -> flip y
        gj = ny - 1 - (j // k)
        line = bytearray(W * 3)
        for i in range(W):
            r, g, b = _jet((grid[gj][i // k] - vmin) / span)
            line[3 * i:3 * i + 3] = bytes((r, g, b))
        png_rows.append(bytes(line))
    _write_png(args.out, W, H, png_rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
