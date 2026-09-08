#!/usr/bin/env python3
"""РП3: прямое сравнение схемы ВКР и новой схемы на общей метрике.

ВКР [2] приводит для задачи Брио–Ву таблицу «усреднённой ошибки» на сетках
N = 64…512 с self-reference на N = 1024. Само определение усреднённой ошибки в
ВКР не выписано, поэтому здесь оно фиксируется явно и считается **одинаково для
обеих схем** — иначе сравнение бессмысленно:

    delta_N = (1/N) * sum_i sum_{v in (rho, u, p, By)} |q_v(x_i) - q_v^ref(x_i)|

то есть среднее по ячейкам от суммы модулей ошибки по четырём основным
переменным. Наблюдаемый порядок: p = log2(delta_N / delta_2N).

Считаются две версии одной и той же метрики:

  independent   эталон — центральная схема Куртганова–Тадмора при N = 6400
                (tests/briowu_reference.cpp), не разделяющая с проверяемыми
                схемами ни строки кода. Это основная таблица.
  self          эталон — сама же схема на самой мелкой сетке (протокол ВКР).
                Приводится, чтобы показать, во что обходится self-reference:
                схема сравнивается сама с собой, и порядок систематически
                завышается.

Дополнительно: ширина каждого фронта по уровню 10–90 % в ячейках. Именно она
отвечает на вопрос ТЗ «во сколько раз резче», и её нельзя подменять L1: схема
может выиграть по L1 за счёт гладких участков, не тронув разрывы.
"""
from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import briowu_fronts as bf  # noqa: E402

VARS = ("rho", "u", "p", "By")


def load_legacy_profile(path: Path) -> list[dict]:
    """Профиль legacy_corrected после проекции с треугольной сетки.

    scripts/project_legacy_vtu.py даёт другие имена столбцов (pressure, vx, by)
    и уже приводит x к [0, 1], поэтому сдвиг координаты здесь не нужен --
    в отличие от одномерного legacy-драйвера, который считает на [-0.5, 0.5].
    """
    import csv
    rename = {"rho": "rho", "vx": "u", "pressure": "p", "by": "By"}
    rows = []
    for r in csv.DictReader(path.open()):
        row = {"x": float(r["x"])}
        for src, dst in rename.items():
            row[dst] = float(r[src])
        rows.append(row)
    if not rows:
        raise SystemExit(f"{path}: пустой профиль")
    return rows


def delta_n(cand: list[dict], ref: list[dict]) -> dict:
    """Метрика ВКР и обычные L1/L2, все на сетке кандидата."""
    xs = [r["x"] for r in cand]
    n = len(xs)
    per_var, total = {}, 0.0
    l2_total = 0.0
    for var in VARS:
        rv = bf.interp(ref, xs, var)
        cv = [r[var] for r in cand]
        diff = [cv[i] - rv[i] for i in range(n)]
        scale = max(abs(v) for v in rv) or 1.0
        l1 = sum(abs(d) for d in diff) / n
        l2 = math.sqrt(sum(d * d for d in diff) / n)
        per_var[var] = {"l1": l1, "l2": l2, "relative_l1": l1 / scale,
                        "linf": max(abs(d) for d in diff)}
        total += l1
        l2_total += l2
    return {"delta_n": total, "l1_sum": total, "l2_sum": l2_total,
            "per_variable": per_var}


def widths(cand: list[dict], ref: list[dict]) -> dict:
    xs = [r["x"] for r in cand]
    dx = xs[1] - xs[0]
    out = {}
    for fr in bf.find_fronts(ref, "rho", 3) + bf.find_fronts(ref, "By", 3):
        w = bf.front_width_cells(cand, fr["variable"], fr["x"], fr["jump"], 25.0 * dx)
        out[f'{fr["variable"]}@{fr["x"]:.4f}'] = w
    return out


def run_amrex(verify: str, nx: int, out: Path, limiter: str, integrator: str,
              emf: str, cfl: float) -> float:
    t0 = time.perf_counter()
    r = subprocess.run([verify, "briowu1d", str(nx), limiter, integrator, emf,
                        str(cfl), str(out)],
                       text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    if r.returncode:
        print(r.stdout, file=sys.stderr)
        raise SystemExit(f"AMReX run failed at N={nx}")
    return time.perf_counter() - t0


def observed_order(rows: list[dict], key: str) -> None:
    """Проставляет наблюдаемый порядок между соседними сетками."""
    for i in range(1, len(rows)):
        a, b = rows[i - 1], rows[i]
        if a.get(key) and b.get(key) and b[key] > 0.0:
            ratio = b["nx"] / a["nx"]
            rows[i][key + "_order"] = math.log(a[key] / b[key]) / math.log(ratio)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--verify", required=True)
    ap.add_argument("--reference-csv", required=True, type=Path,
                    help="независимый эталон KT (mhd2d_briowu_reference)")
    ap.add_argument("--resolutions", default="64,128,256,400,512")
    ap.add_argument("--self-reference", type=int, default=1024,
                    help="сетка для self-reference (протокол ВКР)")
    ap.add_argument("--cfl", type=float, default=0.1)
    ap.add_argument("--legacy-profiles", type=Path,
                    default=Path("benchmarks/raw/rp3_convergence"),
                    help="каталог с прогонами legacy_corrected: "
                         "<dir>/legacy_bw_<N>/brio_profile.csv")
    ap.add_argument("--raw-dir", type=Path,
                    default=Path("benchmarks/raw/rp3_convergence"))
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()

    resolutions = [int(v) for v in args.resolutions.split(",")]
    args.raw_dir.mkdir(parents=True, exist_ok=True)
    ref = bf.load(args.reference_csv)

    # --- новая схема: рабочая конфигурация N3 и её же self-reference --------
    self_ref_path = args.raw_dir / f"amrex_n3_{args.self_reference}.csv"
    if not self_ref_path.is_file():
        run_amrex(args.verify, args.self_reference, self_ref_path,
                  "mc", "rk2", "gs", args.cfl)
    amrex_self_ref = bf.load(self_ref_path, args.self_reference)

    schemes: dict[str, list[dict]] = {"amrex_n3": [], "amrex_n0": [], "legacy_corrected": []}

    for nx in resolutions:
        for label, (lim, ti, emf) in (("amrex_n3", ("mc", "rk2", "gs")),
                                      ("amrex_n0", ("none", "euler", "bs"))):
            csv_path = args.raw_dir / f"{label}_{nx}.csv"
            wall = run_amrex(args.verify, nx, csv_path, lim, ti, emf, args.cfl)
            cand = bf.load(csv_path, nx)
            row = {"nx": nx, "wall_s": wall, "csv": str(csv_path)}
            row["independent"] = delta_n(cand, ref)
            row["self"] = delta_n(cand, amrex_self_ref)
            row["front_widths_cells"] = widths(cand, ref)
            schemes[label].append(row)
            print(f"  {label:16s} N={nx:<5d} delta_N(незав.)={row['independent']['delta_n']:.6f}  "
                  f"delta_N(self)={row['self']['delta_n']:.6f}  {wall:.2f} s")

        # --- историческая схема: треугольная сетка, первый порядок ----------
        legacy_csv = args.legacy_profiles / f"legacy_bw_{nx}" / "brio_profile.csv"
        if not legacy_csv.is_file():
            print(f"  legacy N={nx}: профиль отсутствует ({legacy_csv}) — пропуск")
            continue
        cand = load_legacy_profile(legacy_csv)
        row = {"nx": nx, "csv": str(legacy_csv)}
        row["independent"] = delta_n(cand, ref)
        row["front_widths_cells"] = widths(cand, ref)
        schemes["legacy_corrected"].append(row)
        print(f"  {'legacy_corrected':16s} N={nx:<5d} "
              f"delta_N(незав.)={row['independent']['delta_n']:.6f}")

    for label, rows in schemes.items():
        rows.sort(key=lambda r: r["nx"])
        for kind in ("independent", "self"):
            flat = [{"nx": r["nx"], "d": r[kind]["delta_n"]} for r in rows if kind in r]
            observed_order(flat, "d")
            for r, f in zip([r for r in rows if kind in r], flat):
                if "d_order" in f:
                    r[kind]["observed_order"] = f["d_order"]

    record = {
        "schema_version": 1,
        "case": "brio_wu",
        "t": 0.1,
        "cfl": args.cfl,
        "metric": "delta_N = mean_i sum_v |q_v - q_v_ref|, v in (rho, u, p, By)",
        "metric_note": "определение ВКР не выписано в источнике; здесь оно "
                       "зафиксировано явно и применено к обеим схемам одинаково",
        "independent_reference": str(args.reference_csv),
        "self_reference_nx": args.self_reference,
        "schemes": schemes,
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
