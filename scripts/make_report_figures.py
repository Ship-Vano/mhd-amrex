#!/usr/bin/env python3
"""Produce compact pgfplots data files for docs/report.tex from verified runs.

Pure standard library. Every output file carries a header comment naming the
source run so the report figures stay traceable (RPT-004).

Inputs (regenerate first if absent):
  ./build/release/mhd2d_verify briowu1d 400  none euler bs 0.1 /tmp/bw_N0.csv
  ./build/release/mhd2d_verify briowu1d 400  mc   rk2   gs 0.1 /tmp/bw_N3.csv
  ./build/release/mhd2d_verify briowu1d 2048 mc   rk2   gs 0.1 /tmp/bw_ref.csv
  ./build/release/mhd2d_verify ot 128        # -> out_ot.csv
  ./build/release/mhd2d_verify rotor 128     # -> out_rotor.csv
  ./build/release/mhd2d_verify alfven {16,32,64,128}  # -> out_alfven_*.csv
  benchmarks/raw/legacy_corrected/brio_wu_struct_400x8_cmp/brio_profile.csv
  a legacy rotor VTU (task 4) + its structured mesh
"""
from __future__ import annotations

import argparse
import csv
import math
import re
from pathlib import Path

OUT = Path("docs/figures/data")


def _write(name: str, header: str, cols: list[str], rows: list[tuple]):
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / name).open("w") as f:
        f.write(f"# {header}\n")
        f.write(" ".join(cols) + "\n")
        for r in rows:
            f.write(" ".join(f"{v:.6g}" for v in r) + "\n")
    print(f"wrote {OUT / name}  ({len(rows)} rows)")


def _collapse_strip(csv_path: Path, nx: int):
    """standalone briowu1d CSV -> 1-D (x, rho, u, By, p), y-column mean, nx bins."""
    xs: dict[int, dict] = {}
    with csv_path.open() as fh:
        for row in csv.DictReader(fh):
            k = int(round(float(row["x"]) * nx - 0.5))
            b = xs.setdefault(k, {"x": 0.0, "rho": 0.0, "u": 0.0, "By": 0.0, "p": 0.0, "n": 0})
            for key in ("x", "rho", "u", "By", "p"):
                b[key] += float(row[key])
            b["n"] += 1
    out = []
    for k in sorted(xs):
        b = xs[k]
        out.append((b["x"] / b["n"], b["rho"] / b["n"], b["u"] / b["n"],
                    b["By"] / b["n"], b["p"] / b["n"]))
    return out


def _downsample(rows, target: int):
    if len(rows) <= target:
        return rows
    step = len(rows) / target
    return [rows[min(len(rows) - 1, int(i * step))] for i in range(target)]


def briowu(args):
    src = {
        "legacy": Path("benchmarks/raw/legacy_corrected/brio_wu_struct_400x8_cmp/brio_profile.csv"),
        "n0": Path(args.bw_n0), "n3": Path(args.bw_n3), "ref": Path(args.bw_ref),
    }
    # legacy area-weighted bin profile
    leg = []
    for row in csv.DictReader(src["legacy"].open()):
        leg.append((float(row["x"]), float(row["rho"]), float(row["vx"]),
                    float(row["by"]), float(row["pressure"])))
    _write("briowu1d_legacy.dat",
           f"legacy_corrected Brio-Wu 400x8 structured, t=0.1, CFL 0.1; source {src['legacy']}",
           ["x", "rho", "u", "By", "p"], leg)
    for tag, nx, ndown in (("n0", 400, 200), ("n3", 400, 200), ("ref", 2048, 400)):
        rows = _downsample(_collapse_strip(src[tag], nx), ndown)
        label = {"n0": "AMReX N0 (const+Euler+BS)", "n3": "AMReX N3 (MUSCL+SSPRK2+GS)",
                 "ref": "AMReX N3 N=2048 provisional reference"}[tag]
        _write(f"briowu1d_{tag}.dat",
               f"{label}, Nx={nx}, t=0.1, CFL 0.1; source {src[tag]}",
               ["x", "rho", "u", "By", "p"], rows)


def alfven(args):
    rows = []
    for n in (16, 32, 64, 128):
        p = Path(f"out_alfven_{n}.csv")
        # the L1 is printed by the driver, not in the CSV; recompute from CSV vs exact
        # simpler: read it from a stored value the caller passes, else recompute
        d = list(csv.DictReader(p.open()))
        # exact CP-Alfven at t=1 == initial; reconstruct B_perp error
        num = den = 0.0
        ca, sa = math.cos(math.pi / 6), math.sin(math.pi / 6)
        for r in d:
            x, y = float(r["x"]), float(r["y"])
            phi = 2 * math.pi * (x * ca + y * sa)
            bx, by, bz = float(r["Bx"]), float(r["By"]), float(r["Bz"])
            bperp = -bx * sa + by * ca
            ref = 0.1 * math.sin(phi)
            num += abs(bperp - ref)
            den += abs(ref)
        rows.append((n, num / den))
    _write("alfven_conv.dat",
           "CP-Alfven B_perp relative L1 vs N; standalone mhd2d_verify, t=1, CFL 0.4",
           ["N", "L1_Bperp"], rows)


def ot(args):
    rows_all = [r for r in csv.DictReader(Path("out_ot.csv").open())
                if abs(float(r["y"]) - 0.3125) < 0.012]
    nb = 96
    b = [[] for _ in range(nb)]
    for r in rows_all:
        b[min(nb - 1, int(float(r["x"]) * nb))].append((float(r["rho"]), float(r["p"])))
    rows = []
    for i, c in enumerate(b):
        if c:
            rows.append(((i + 0.5) / nb,
                         sum(v[0] for v in c) / len(c),
                         sum(v[1] for v in c) / len(c)))
    _write("ot_slice.dat",
           "AMReX Orszag-Tang N=128, t=0.5, y=0.3125 slice (96 x-bins); source out_ot.csv",
           ["x", "rho", "p"], rows)


def rotor(args):
    def diag(rows_iter, key):
        pts = []
        for r in rows_iter:
            x, y = r[0], r[1]
            if abs(x - y) < 0.02 * math.sqrt(2):
                pts.append((x, r[2], r[3]))
        pts.sort()
        nb = 48
        b = [[] for _ in range(nb)]
        for x, rho, p in pts:
            b[min(nb - 1, int(x * nb))].append((rho, p))
        return [((i + 0.5) / nb, sum(v[0] for v in c) / len(c), sum(v[1] for v in c) / len(c))
                for i, c in enumerate(b) if c]

    amrex = ((float(r["x"]), float(r["y"]), float(r["rho"]), float(r["p"]))
             for r in csv.DictReader(Path("out_rotor.csv").open()))
    _write("rotor_diag_amrex.dat",
           "AMReX rotor N=128, t=0.15, x=y diagonal (48 bins); source out_rotor.csv",
           ["x", "rho", "p"], diag(amrex, None))

    if args.rotor_vtu and args.rotor_mesh:
        import sys
        sys.path.insert(0, "scripts")
        from analyze_legacy_vtu import primitive, read_mesh, read_vtu
        n, e = read_mesh(Path(args.rotor_mesh))
        s = read_vtu(Path(args.rotor_vtu))
        leg = []
        for el, st in zip(e, s):
            (x0, y0), (x1, y1), (x2, y2) = (n[i] for i in el)
            cx, cy = (x0 + x1 + x2) / 3, (y0 + y1 + y2) / 3
            rho, p, *_ = primitive(st, 1.4)
            leg.append((cx, cy, rho, p))
        _write("rotor_diag_legacy.dat",
               f"legacy_corrected rotor 128x128 structured, t=0.15, CFL 0.5, x=y diagonal; "
               f"source {args.rotor_vtu}",
               ["x", "rho", "p"], diag(iter(leg), None))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--bw-n0", default="/tmp/bw_N0.csv")
    ap.add_argument("--bw-n3", default="/tmp/bw_N3.csv")
    ap.add_argument("--bw-ref", default="/tmp/bw_ref.csv")
    ap.add_argument("--rotor-vtu")
    ap.add_argument("--rotor-mesh")
    ap.add_argument("--only", choices=("briowu", "alfven", "ot", "rotor"))
    args = ap.parse_args()
    todo = [args.only] if args.only else ["briowu", "alfven", "ot", "rotor"]
    for name in todo:
        globals()[name](args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
