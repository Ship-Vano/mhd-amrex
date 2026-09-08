#!/usr/bin/env python3
"""РП1/РП2: откуда берутся выбросы у фронтов Брио–Ву и что их убирает.

Замечание к отчёту было конкретным: у схемы второго порядка на Брио–Ву видны
перелёты, а отчёт не говорил, что именно их порождает — реконструкция,
интегратор или усреднение ЭДС. Здесь это разделяется измерением, а не
рассуждением: матрица прогонов, в которой за раз меняется ровно одна деталь
схемы, и общая метрика против эталона, не имеющего с проверяемой схемой ни
одной общей строки кода (центральная схема Куртганова–Тадмора,
tests/briowu_reference.cpp).

Три среза матрицы:

  ladder      кусочно-постоянная -> +реконструкция -> +SSP-RK2 -> +Gardiner-Stone.
              Отвечает на вопрос «какая деталь схемы вносит выброс».
  limiter     minmod / van Leer / MC при прочих равных (РП1).
  integrator  Эйлер / средняя точка (не-SSP) / SSP-RK2 / TVD-RK3 (РП2).

Для каждого прогона: относительная L1 по (rho, u, p, By), избыток полной
вариации над эталоном, амплитуда перелёта/недолёта за диапазон эталона, число
смен знака первой разности и ширина фронтов 10-90 % в ячейках. Ширина и
амплитуда выброса разнесены намеренно: лимитер, который «убирает осцилляции»
размазав контактный разрыв вдвое, — не улучшение, и таблица обязана это
показывать.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import briowu_fronts as bf  # noqa: E402

VARS = ("rho", "u", "p", "By")

# (метка, лимитер, интегратор, ЭДС, переменные реконструкции, CFL или None, срез, пояснение)
MATRIX = [
    ("N0  const+Euler+BS",      "none",    "euler",    "bs", "prim", None, "ladder",
     "кусочно-постоянная реконструкция: схема ВКР в декартовом варианте"),
    ("N1  MUSCL-MC+Euler+BS",   "mc",      "euler",    "bs", "prim", None, "ladder",
     "добавлена только реконструкция"),
    ("N2  MUSCL-MC+SSP-RK2+BS", "mc",      "rk2",      "bs", "prim", None, "ladder",
     "добавлен SSP-RK2"),
    ("N3  MUSCL-MC+SSP-RK2+GS", "mc",      "rk2",      "gs", "prim", None, "ladder",
     "добавлено усреднение ЭДС Gardiner-Stone"),

    ("minmod",                  "minmod",  "rk2",      "gs", "prim", None, "limiter",
     "самый диссипативный из трёх TVD-лимитеров"),
    ("van Leer",                "vanleer", "rk2",      "gs", "prim", None, "limiter",
     "гладкий, промежуточный по компрессии"),
    ("MC",                      "mc",      "rk2",      "gs", "prim", None, "limiter",
     "monotonized central, наиболее компрессивный"),

    ("MC, примитивные",         "mc",      "rk2",      "gs", "prim", None, "recon",
     "рабочий вариант: ограничиваются ρ, u, p, B"),
    ("MC, консервативные",      "mc",      "rk2",      "gs", "cons", None, "recon",
     "ограничиваются ρ, ρv, e, B; на грань выдаются примитивы"),
    ("minmod, примитивные",     "minmod",  "rk2",      "gs", "prim", None, "recon",
     "то же для наиболее диссипативного лимитера"),
    ("minmod, консервативные",  "minmod",  "rk2",      "gs", "cons", None, "recon",
     "то же для наиболее диссипативного лимитера"),

    ("Эйлер",                   "mc",      "euler",    "gs", "prim", None, "integrator",
     "первый порядок по времени"),
    ("RK2 средняя точка",       "mc",      "midpoint", "gs", "prim", None, "integrator",
     "второй порядок, НЕ SSP: контрольный вариант"),
    ("SSP-RK2 (Хойн)",          "mc",      "rk2",      "gs", "prim", None, "integrator",
     "рабочий интегратор"),
    ("TVD-RK3",                 "mc",      "rk3",      "gs", "prim", None, "integrator",
     "три стадии, SSP-константа 1"),
]

# Различие SSP и не-SSP интегратора проявляется не при любом шаге, а при шаге,
# близком к границе устойчивости: SSP-свойство — это утверждение о том, до
# какого числа Куранта TVD-оценка шага Эйлера переносится на весь шаг. Поэтому
# сравнение ведётся развёрткой по CFL, а не в одной точке.
CFL_SWEEP_INTEGRATORS = ("rk2", "midpoint", "rk3")
CFL_SWEEP_VALUES = (0.1, 0.4, 0.6, 0.8, 0.9)


def run_case(verify: str, nx: int, lim: str, ti: str, emf: str,
             cfl: float, out: Path, rvars: str = "prim") -> tuple[float, int, int]:
    t0 = time.perf_counter()
    r = subprocess.run([verify, "briowu1d", str(nx), lim, ti, emf, str(cfl),
                        str(out), rvars],
                       text=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    wall = time.perf_counter() - t0
    if r.returncode:
        print(r.stdout, file=sys.stderr)
        raise SystemExit(f"run failed: {lim}/{ti}/{emf}/{rvars}")
    tail = r.stdout.strip().splitlines()[-1]
    def field(name: str) -> int:
        for tok in tail.split():
            if tok.startswith(name + "="):
                return int(tok.split("=")[1])
        return -1
    return wall, field("hlld_fallbacks"), field("recon_fallbacks")


def score(reference: Path, candidate: Path, nx: int, label: str) -> dict:
    ref = bf.load(reference)
    cand = bf.load(candidate, nx)
    xs = [r["x"] for r in cand]
    dx = xs[1] - xs[0]
    out = {"label": label, "cells": len(xs), "variables": {}, "fronts": []}
    for var in VARS:
        rv = bf.interp(ref, xs, var)
        cv = [r[var] for r in cand]
        n = len(xs)
        diff = [cv[i] - rv[i] for i in range(n)]
        scale = max(abs(v) for v in rv) or 1.0
        out["variables"][var] = {
            "relative_l1": sum(abs(d) for d in diff) / n / scale,
            "relative_l2": (sum(d * d for d in diff) / n) ** 0.5 / scale,
            "linf": max(abs(d) for d in diff),
            "tv_excess": bf.total_variation(cv) - bf.total_variation(rv),
            "overshoot": max(0.0, max(cv) - max(rv)),
            "undershoot": max(0.0, min(rv) - min(cv)),
            "first_difference_sign_changes": bf.sign_changes(cv, 1.0e-6 * scale),
        }
    # Фронты те же, что и в T06: три по rho (включая контактный разрыв) и три
    # по By (составная волна и медленная ударная).
    for fr in bf.find_fronts(ref, "rho", 3) + bf.find_fronts(ref, "By", 3):
        var = fr["variable"]
        out["fronts"].append({
            "variable": var,
            "reference_x": fr["x"],
            "width_cells_10_90": bf.front_width_cells(cand, var, fr["x"],
                                                      fr["jump"], 25.0 * dx),
        })
    # Сводные числа для таблицы отчёта: наибольший перелёт по всем переменным в
    # долях диапазона эталона и ширина самого резкого фронта rho (контакт).
    ref_range = {v: (max(r[v] for r in ref) - min(r[v] for r in ref)) or 1.0
                 for v in VARS}
    out["max_relative_overshoot"] = max(
        max(out["variables"][v]["overshoot"], out["variables"][v]["undershoot"])
        / ref_range[v] for v in VARS)
    rho_widths = [f["width_cells_10_90"] for f in out["fronts"]
                  if f["variable"] == "rho" and f["width_cells_10_90"] is not None]
    out["contact_width_cells"] = min(rho_widths) if rho_widths else None
    out["mean_relative_l1"] = sum(out["variables"][v]["relative_l1"]
                                  for v in VARS) / len(VARS)
    out["total_sign_changes"] = sum(
        out["variables"][v]["first_difference_sign_changes"] for v in VARS)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--verify", required=True)
    ap.add_argument("--reference-csv", required=True, type=Path,
                    help="готовый эталон KT (см. mhd2d_briowu_reference)")
    ap.add_argument("--nx", type=int, default=400)
    ap.add_argument("--cfl", type=float, default=0.1)
    ap.add_argument("--raw-dir", type=Path,
                    default=Path("benchmarks/raw/rp1_limiters"))
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()

    args.raw_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for label, lim, ti, emf, rvars, cfl_override, group, note in MATRIX:
        cfl = args.cfl if cfl_override is None else cfl_override
        # rvars входит в имя файла: иначе примитивный и консервативный прогоны
        # с одинаковыми лимитером/интегратором/ЭДС затирали бы друг друга.
        csv_path = (args.raw_dir /
                    f"briowu_{lim}_{ti}_{emf}_{rvars}_{args.nx}.csv")
        wall, hlld_fb, recon_fb = run_case(args.verify, args.nx, lim, ti, emf,
                                           cfl, csv_path, rvars)
        rec = score(args.reference_csv, csv_path, args.nx, label)
        rec.update({"group": group, "limiter": lim, "integrator": ti,
                    "emf": emf, "recon_vars": rvars, "cfl": cfl,
                    "note": note, "wall_s": wall,
                    "hlld_fallbacks": hlld_fb,
                    "recon_fallbacks": recon_fb,
                    "csv": str(csv_path)})
        rows.append(rec)
        print(f"  {label:26s} L1={rec['mean_relative_l1']:.3e}  "
              f"overshoot={rec['max_relative_overshoot']:.3e}  "
              f"contact={rec['contact_width_cells']}  "
              f"signchg={rec['total_sign_changes']}  "
              f"recon_fb={recon_fb}  {wall:.2f} s")

    record = {
        "schema_version": 1,
        "case": "brio_wu",
        "t": 0.1,
        "nx": args.nx,
        "cfl": args.cfl,
        "reference": str(args.reference_csv),
        "reference_scheme": "Kurganov-Tadmor central, no Riemann solver, "
                            "no shared code with the scheme under test",
        "rows": rows,
    }
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(record, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
