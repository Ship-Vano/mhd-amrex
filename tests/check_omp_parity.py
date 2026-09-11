#!/usr/bin/env python3
"""Обязательная матрица: независимость результата от числа потоков OpenMP.

Парный тест к `check_mpi_parity.py`. MPI-версия ловит зависимость от разбиения
по рангам, эта — от разбиения по тайлам внутри ранга: распараллеливание идёт
по MFIter, и общее состояние трогают и обновление ячеек, и накопление
диагностических счётчиков.

Что тест гарантирует: физические величины, дрейф сохраняющихся, норма div B и
все три счётчика совпадают при 1, 2 и 4 потоках. Систематическая ошибка
распараллеливания — перепутанный тайл, потерянное обновление, счётчик, который
сбрасывается на каждом потоке, — даёт расхождение сразу и здесь ловится.

Чего тест НЕ гарантирует, и это проверено экспериментом: он не детектор гонок.
Сломанная атомарность счётчиков (замена `HostDevice::Atomic::Add` на обычный
`+=`) этот тест проходит — на `mhd_blast` счётчик трогают всего 48 раз за 241
шаг, и окно гонки практически не воспроизводится. Для доказательства
отсутствия гонок нужен сборка с ThreadSanitizer, а не сравнение результатов.
Поэтому формулировка в отчёте должна быть «результат не зависит от числа
потоков», а не «гонок нет».

Задача выбирается с ненулевыми счётчиками: на задаче с нулями сравнивать было
бы нечего.
"""
from __future__ import annotations

import argparse
import os
import re
import subprocess
import sys

# Величины, которые обязаны совпасть: физика, дрейф сохраняющихся, div B и
# все три счётчика.
FLOAT_KEYS = ("rho_min", "rho_max", "p_min", "p_max", "rho_rel_drift", "divb")
INT_KEYS = ("steps", "hlld_fallbacks", "floor_events", "nonpositive_cells")


def extract(text: str) -> dict:
    out: dict = {}
    m = re.search(r"rho_min=([0-9.eE+-]+) rho_max=([0-9.eE+-]+) "
                  r"p_min=([0-9.eE+-]+) p_max=([0-9.eE+-]+)", text)
    if m:
        out.update({k: float(m.group(i)) for i, k in
                    enumerate(("rho_min", "rho_max", "p_min", "p_max"), start=1)})
    m = re.search(r"rho_rel_drift=([0-9.eE+-]+)", text)
    if m:
        out["rho_rel_drift"] = float(m.group(1))
    m = re.search(r"normalized=([0-9.eE+-]+)", text)
    if m:
        out["divb"] = float(m.group(1))
    m = re.search(r"Evolve finished: (\d+) steps, t=[0-9.eE+-]+, "
                  r"hlld_fallbacks=(\d+), floor_events=(\d+), "
                  r"nonpositive_cells=(\d+)", text)
    if m:
        out.update({k: int(m.group(i)) for i, k in enumerate(INT_KEYS, start=1)})
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--executable", required=True)
    ap.add_argument("--config", required=True)
    ap.add_argument("--threads", default="2,4")
    ap.add_argument("--tol", type=float, default=1.0e-12)
    ap.add_argument("--require-nonzero-counters", action="store_true",
                    help="падать, если на выбранной задаче счётчики нулевые: "
                         "тогда сравнивать по ним нечего")
    args = ap.parse_args()

    def run(threads: int) -> tuple[dict, str]:
        env = dict(os.environ, OMP_NUM_THREADS=str(threads))
        r = subprocess.run([args.executable, args.config], env=env, text=True,
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        if r.returncode != 0:
            print(r.stdout[-2000:], file=sys.stderr)
            raise SystemExit(f"run failed with OMP_NUM_THREADS={threads}")
        return extract(r.stdout), r.stdout

    serial, serial_log = run(1)
    missing = [k for k in FLOAT_KEYS + INT_KEYS if k not in serial]
    if missing:
        print(f"regression: serial run printed no {missing}", file=sys.stderr)
        return 1

    # Сборка без OpenMP делает тест бессодержательным: он сравнил бы прогон сам
    # с собой. Это не провал, но и не проверка -- так и сообщаем.
    if "OMP initialized" not in serial_log:
        print("built without OpenMP; thread parity is vacuous here")
        return 0

    if args.require_nonzero_counters and not (serial["hlld_fallbacks"] or
                                              serial["floor_events"]):
        print("regression: counters are zero on this case, so there is nothing "
              "to compare -- pick a case that exercises them", file=sys.stderr)
        return 1

    print(f"  threads=1: {serial}")
    failures = []
    for n in (int(x) for x in args.threads.split(",")):
        got, _ = run(n)
        for key in INT_KEYS:
            if got.get(key) != serial[key]:
                failures.append(f"threads={n}: {key} {got.get(key)} != {serial[key]}")
        for key in FLOAT_KEYS:
            cur, ref = got.get(key), serial[key]
            if cur is None:
                failures.append(f"threads={n}: missing {key}")
            elif abs(cur - ref) > args.tol * max(abs(ref), 1.0):
                failures.append(f"threads={n}: {key} {cur:.12g} != {ref:.12g}")
        print(f"  threads={n}: {got}")

    if failures:
        for f in failures:
            print(f"regression: {f}", file=sys.stderr)
        return 1
    print("OpenMP thread parity holds (physics, drift, div B, counters)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
