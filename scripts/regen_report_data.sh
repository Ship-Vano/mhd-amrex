#!/bin/sh
# Deterministically regenerate every input the report figures are built from.
#
# The figure data files record where they came from. Recording a /tmp path makes
# that record useless: the file is gone by the time anyone checks. This script
# writes all intermediates to one fixed, repo-relative directory instead, so a
# .dat header names something that can actually be reproduced.
#
#   sh scripts/regen_report_data.sh && python3 scripts/make_report_figures.py
#
set -eu

ROOT=$(cd "$(dirname "$0")/.." && pwd)
OUT="$ROOT/benchmarks/raw/report_inputs"
BUILD="${BUILD_DIR:-$ROOT/build/release}"
VERIFY="$BUILD/mhd2d_verify"
REFERENCE="$BUILD/briowu_reference"

for exe in "$VERIFY" "$REFERENCE"; do
    [ -x "$exe" ] || { echo "missing $exe; build first: cmake --build $BUILD" >&2; exit 1; }
done
mkdir -p "$OUT"
cd "$ROOT"

echo "== Brio-Wu 1-D ablation (N=400) and provisional AMReX reference (N=2048)"
"$VERIFY" briowu1d 400  none euler bs 0.1 "$OUT/bw_n0_400.csv"  > /dev/null
"$VERIFY" briowu1d 400  mc   rk2   gs 0.1 "$OUT/bw_n3_400.csv"  > /dev/null
"$VERIFY" briowu1d 2048 mc   rk2   gs 0.1 "$OUT/bw_amrex_2048.csv" > /dev/null

echo "== Brio-Wu independent reference (Kurganov-Tadmor, N=6400)"
"$REFERENCE" 6400 0.1 "$OUT/bw_kt_6400.csv" 0.4

echo "== Dai-Woodward 1-D"
"$VERIFY" dw1d 400  none euler bs 0.2 "$OUT/dw_n0_400.csv"  > /dev/null
"$VERIFY" dw1d 400  mc   rk2   gs 0.2 "$OUT/dw_n3_400.csv"  > /dev/null
"$VERIFY" dw1d 2048 mc   rk2   gs 0.2 "$OUT/dw_ref_2048.csv" > /dev/null

echo "== 2-D canonical dumps (Orszag-Tang, rotor, field loop, CP-Alfven)"
"$VERIFY" ot    128 > /dev/null && mv -f out_ot.csv    "$OUT/ot_128.csv"
"$VERIFY" rotor 128 > /dev/null && mv -f out_rotor.csv "$OUT/rotor_128.csv"
"$VERIFY" loop  128 2.0 0.1 > /dev/null && mv -f out_loop_128.csv "$OUT/loop_128.csv"
for n in 16 32 64 128; do
    "$VERIFY" alfven "$n" > /dev/null && mv -f "out_alfven_$n.csv" "$OUT/alfven_$n.csv"
done

echo
echo "regenerated into $OUT"
ls -1 "$OUT"
