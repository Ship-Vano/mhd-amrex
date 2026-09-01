# Cross-solver comparison — Brio–Wu, 1-D, uniform grid (T08 Level A, provisional)

> Broader suite-wide comparison (Brio–Wu, Dai–Woodward, CP-Alfvén, rotor,
> Orszag–Tang, magnetic loop) is in
> [`CROSS_SOLVER_COMPARISON.md`](CROSS_SOLVER_COMPARISON.md). This note keeps the
> detailed Brio–Wu N0–N3 ablation and the equal-cost timing; note its
> `legacy_corrected` profile was the 2-D triangular strip, whereas
> `CROSS_SOLVER_COMPARISON.md §1` uses a true 1-D legacy driver (which agrees
> with the AMReX kernel to ~6–10 %).

Date: 2026-08-31. Scope: **equal-spacing and equal-cost** views of the
`legacy_corrected` unstructured solver against the AMReX scheme and its N0–N3
ablation, on the classic Brio–Wu shock tube. This is Level A of the
`LEGACY_BASELINE.md` §6 protocol (pure scheme comparison on a shared 1-D
profile). It is **provisional**: the reference is a grid-converged AMReX N3 run,
not an independent exact Riemann solution — that is a T06 deliverable. Absolute
errors here are therefore self-consistent scheme comparisons, not
literature-validated.

## Setup (identical physics for every run)

| | value |
|---|---|
| states | left `(ρ,p,u,v,w,By,Bz)=(1,1,0,0,0,1,0)`, right `(0.125,0.1,0,0,0,-1,0)` |
| `Bx` | 0.75 (constant) |
| `γ` | 2 |
| `t_end` | 0.1 (last step clipped exactly) |
| CFL | 0.1 |
| x boundaries | frozen (Dirichlet = initial state) |
| y | thin strip, periodic; profile is the y-collapsed 1-D column mean, with the max within-column σ reported |
| branch | compound / non-regular (owner decision D-003) |

- **Legacy** (`legacy_corrected`, `--mesh-backend structured`): piecewise-constant
  FV, forward Euler, corrected HLLD, nodal-averaged EMF, RT0 cell-`B`, on a
  structured `nx×ny` right-triangle mesh (2 triangles / rectangle). Runs:
  `128×4` and `400×8`.
- **AMReX** (`mhd2d_verify briowu1d`): Cartesian FV, HLLD from the staggered
  `Bn`, staggered CT. Ablation per `TECHNICAL_SPECIFICATION.md` §6:

  | id | space | time | corner EMF |
  |---|---|---|---|
  | N0 | piecewise constant | forward Euler | Balsara–Spicer |
  | N1 | MUSCL (MC) | forward Euler | Balsara–Spicer |
  | N2 | MUSCL (MC) | SSP-RK2 (Heun) | Balsara–Spicer |
  | N3 | MUSCL (MC) | SSP-RK2 (Heun) | Gardiner–Stone (ε=2) |

- **Reference**: AMReX N3 at N=2048 (8840 steps), `max|divB| = 0`, 0 HLLD
  fallbacks. Provisional (see above).

Every run: `max|divB|` at roundoff, `hlld_fallbacks = 0`, `nonpositive_cells = 0`.
Metrics computed by `scripts/compare_briowu_1d.py`; machine-readable output
`benchmarks/summary/briowu_1d_comparison.json`.

## Equal grid spacing (N = 400)

L1 error vs the provisional reference, y-collapsed profile:

| run | ρ L1 | u L1 | p L1 | By L1 | By TV-excess | u over/under | max within-column σ |
|---|---:|---:|---:|---:|---:|---:|---:|
| legacy_corrected 400×8 | 4.26e-2 | 7.60e-2 | 3.72e-2 | 5.63e-2 | 0.68 | +0.10 / 0 | 0.20 (By) |
| AMReX **N0** 400 | 1.08e-2 | 2.07e-2 | 1.09e-2 | 1.33e-2 | 0 (TV ≤ ref) | 0 / 0 | 0 |
| AMReX **N1** 400 | 2.44e-3 | 7.63e-3 | 2.65e-3 | 4.33e-3 | 0.30 | — | 0 |
| AMReX **N2** 400 | 2.55e-3 | 6.64e-3 | 2.63e-3 | 3.97e-3 | 0.18 | — | 0 |
| AMReX **N3** 400 | 2.52e-3 | 6.25e-3 | 2.50e-3 | 3.96e-3 | 0.11 | +0.005 / −0.013 | 0 |

Reading:

1. **Cartesian first-order (N0) already beats legacy first-order by ~4× in L1**
   at the same nominal resolution, and does so while staying strictly monotone
   (TV ≤ reference, zero overshoot) and exactly 1-D (σ = 0). The legacy
   unstructured strip develops a 15–20 % transverse spread (σ ≈ 0.2 in `By`,
   0.16 in `p`) and a spurious `By` TV-excess of 0.68 — the two triangle
   orientations per rectangle see asymmetric stencils, and RT0 + nodal EMF add
   grid noise the Cartesian corner-EMF does not.
2. **MUSCL (N1) buys ~4× more accuracy** but introduces the expected 2nd-order
   Gibbs signature: `By` TV-excess 0.30, `p` undershoot ~2e-3.
3. **SSP-RK2 (N2) and then Gardiner–Stone CT (N3) suppress the excess variation**
   (0.30 → 0.18 → 0.11) and the undershoots (`By` undershoot → 0) with no L1
   penalty — GS adds dissipation exactly where the parasitic oscillation lives.
4. **Net at N = 400: AMReX N3 is ≈ 15× lower L1 than legacy_corrected**, better
   bounded, and exactly 1-D.

The legacy error is partly a systematic offset from the AMReX-family branch (the
reference is not independent), so 15× is a floor on the accuracy gap, not a
calibrated ratio. `128×4` legacy is worse still (ρ L1 6.8e-2, `By` TV-excess 1.18).

## Equal compute cost (diagnostic timing, single serial run, no pinning)

Median of 3 repeats after 1 warm-up, `OMP_NUM_THREADS=1`, this workstation:

| run | steps | wall | vs legacy at same dx |
|---|---:|---:|---:|
| AMReX N0 @ 128 | 852 | 0.065 s | — |
| AMReX N3 @ 128 | 858 | 0.119 s | — |
| **legacy_corrected 128×4** | 4307 | **13.39 s** | 1× |
| AMReX N0 @ 400 | 1704 | 0.43 s | — |
| AMReX N3 @ 200 | 862 | 0.26 s | — |
| AMReX N3 @ 400 | 1728 | 0.94 s | — |
| AMReX N3 @ 800 | 3457 | 3.69 s | — |
| legacy_corrected 400×8 | 10618 | 209 s (1 run) | — |

At **matched grid spacing dx = 1/128**: AMReX first-order (N0) is **≈ 206×**
faster than legacy first-order, and the full second-order scheme (N3) is
**≈ 113×** faster *and* ~10× more accurate. The wall-time gap is structural —
unstructured connectivity, `std::vector<vector<double>>` state, Eigen, a
per-edge state rotation, `omp_set_num_threads` churn, no vectorisation — not a
tuning artefact. legacy also takes ~4300 steps where AMReX takes ~850 because
legacy's arbitrary-triangle CFL bound (`A / Σλℓ`) is stricter than the Cartesian
directional bound at CFL 0.1.

So the equal-cost view is lopsided: AMReX N3 @ 400 (ρ L1 2.5e-3, 0.94 s) is
**222×** faster than legacy 400×8 (209 s) at the same grid spacing, while also
being ~15× more accurate.

**This is a diagnostic, not a benchmark** — one machine, one serial process, no
warm-up statistics beyond 2–3 repeats, no affinity. A proper equal-cost /
equal-error campaign with `S_p`, `E_p`, breakdowns and repeats is T09/T11.

## H1 status for this case

For uniform-grid Brio–Wu, on the accuracy/time axes and against a provisional
reference, **H1 is supported**: AMReX gives ~15× lower L1 at equal spacing and a
~113× (2nd-order) to ~206× (1st-order) wall-time advantage at matched dx.
Caveats that keep this from being a final result:

- reference is not independent (T06);
- one problem, one dimension, no AMR (the AMR error/overhead trade is T07/T08);
- timing is diagnostic;
- part of the legacy error is a discretization offset, not diffusion.

The ablation cleanly separates the contributions: ~4× from Cartesian vs
unstructured (N0 vs legacy), ~4× from MUSCL (N1 vs N0), and the CT/RK choices
trade a little L1-neutral dissipation for much lower excess variation.

## Reproduce

```sh
V=./build/release/mhd2d_verify
for N in 200 400 800 2048; do
  $V briowu1d $N mc rk2 gs 0.1 /tmp/bw_N3_$N.csv          # N3 + reference
done
$V briowu1d 400 none euler bs 0.1 /tmp/bw_N0_400.csv       # ablation
$V briowu1d 400 mc   euler bs 0.1 /tmp/bw_N1_400.csv
$V briowu1d 400 mc   rk2   bs 0.1 /tmp/bw_N2_400.csv

python3 scripts/run_legacy_corrected.py --source /Users/ivansamanov/Documents/MHD2D \
  --case brio_wu --mesh-backend structured --structured-nx 400 --structured-ny 8 \
  --artifact-dir benchmarks/raw/legacy_corrected/brio_wu_struct_400x8_cmp --omp-threads 1

python3 scripts/compare_briowu_1d.py --reference /tmp/bw_N3_2048.csv \
  --amrex N0_400:/tmp/bw_N0_400.csv --amrex N3_400:/tmp/bw_N3_400.csv \
  --legacy legacy_400:benchmarks/raw/legacy_corrected/brio_wu_struct_400x8_cmp/brio_profile.csv \
  --summary benchmarks/summary/briowu_1d_comparison.json
```
