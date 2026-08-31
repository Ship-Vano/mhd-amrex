# legacy_corrected: end-to-end CPU physical-correctness evidence

Date: 2026-08-31. Status of this slice: **PASS** for the CPU physical-correctness
question. This closes the "no new canonical raw runs" gap left open in
[`LEGACY_CORRECTED_T05_STATE.md`](LEGACY_CORRECTED_T05_STATE.md). The remaining
T05 open items (CUDA parity, Brio–Wu branch classification) are unchanged and
still belong to later phases.

## Question

Is the corrected legacy solver *physically correct in its results*, end to end,
on the canonical cases — not just in isolated kernel unit tests?

Historical `legacy_vkr` at commit `9d0f60ea…` is **not**: the controlled T03
runs produce NaN with negative density and pressure on every canonical case
(Brio–Wu at iteration 68 with `rho_min=-0.177`; CP-Alfvén and field-loop at
iteration 0). See [`LEGACY_VKR_T03.md`](LEGACY_VKR_T03.md).

## What was run

The sanctioned runner `scripts/run_legacy_corrected.py` was given a new
dependency-free mesh backend (`--mesh-backend structured`) because the Netgen
Python environment used by the original T05 runs
(`/private/tmp/mhd-netgen-venv`) no longer exists on this host and the Netgen
wheel fails to initialise here. The structured backend emits the exact same
minimal text format via `scripts/legacy_vkr_mesh.py` — two counter-clockwise
right triangles per rectangle — so:

- the corrected and the historical `legacy_vkr` runs now share an identical mesh
  hash for a given resolution;
- the run is reproducible on any host with GCC + CMake + Python, no Netgen;
- the effective `dx` is well defined, which the later cross-solver comparison
  (T08) needs against the AMReX Cartesian grid.

The Netgen backend is retained unchanged for the irregular-triangle robustness
stress runs. Choosing a structured mesh is a **harness** choice, not a physics
change: the overlay `legacy/patches/0001-legacy-corrected-physics.patch`
(SHA-256 `7d7314c6…`) is byte-identical to the one behind
`LEGACY_CORRECTED_T05_STATE.md`, applied with `git apply --check
--whitespace=error` from a clean clone of `9d0f60ea…` on every run.

Each run: clean clone → overlay → structured mesh → Release build (GCC 15.2) →
`ctest` (1/1 pass) → `MHD2D` → area-weighted VTU diagnostics → manifest with
source/overlay/mesh SHA-256 and a quality gate. Raw artifacts (git-ignored) are
in `benchmarks/raw/legacy_corrected/*_struct_*_t05/`; the compact aggregate is
[`benchmarks/summary/legacy_corrected_endtoend.json`](../benchmarks/summary/legacy_corrected_endtoend.json).

## Result — every canonical case is physically admissible

| Case | mesh (rect) | `t_end` | CFL | steps | reached `t_end` | finite | `rho_min` | `p_min` | `max|div B|` (abs / scaled) | HLLD→HLLE fallbacks | max conservation residual | quality gate |
|---|---|---:|---:|---:|:--:|:--:|---:|---:|---|---:|---:|:--:|
| Brio–Wu (task 1) | 128×4 | 0.1 | 0.1 | 4307 | yes | yes | 0.11196 | 0.07985 | 8.8e-12 / 6.2e-15 | 0 | 4.0e-16 | pass |
| CP-Alfvén 30° (task 8) | 32×56 | 1.0 | 0.1 | 2211 | yes | yes | 0.89327 | 0.09226 | 4.4e-13 / 2.5e-15 | 0 | 1.8e-13 | pass |
| Field loop, legacy-scaled `v=(2,1)` (task 9) | 64×32 | 2.0 | 0.1 | 10762 | yes | yes | 1.00000 | 0.99999993 | 4.8e-16 / 2.2e-15 | 0 | 1.9e-14 | pass |
| Field loop, Athena geometry (task 9) | 64×37 | 4.6188 | 0.1 | 18159 | yes | yes | 1.00000 | 0.99999959 | 6.9e-16 / 3.2e-15 | 0 | 1.3e-13 | pass |

- **Finite and positive.** `ensure_admissible_states()` runs pre-flux, post-gas
  and post-CT every step and throws on a non-finite or non-positive `rho`/`p`;
  it never fired. No silent floor is used on the CPU path.
- **Discrete solenoidality.** The scaled magnetic-flux imbalance
  `|Phi_K| / (B_ref · perimeter_K)` stays ≤ 6.2e-15 (double-precision roundoff)
  through every run, including the ~18k-step Athena field loop.
- **Conservation.** For each of mass, the three momenta and total energy, the
  integrated state change plus the integrated boundary flux closes to ≤ 1.8e-13
  (roundoff-scaled; the periodic cases close to ~1e-14).
- **HLLD robustness.** Zero HLLD→HLLE fallbacks on all canonical data: the
  corrected HLLD intermediate states stay admissible without help from the
  fallback path.
- **Brio–Wu structure.** The area-weighted `x`-profile
  (`benchmarks/raw/.../brio_wu_struct_128x4_t05/brio_profile.csv`) shows the
  canonical **compound / non-regular** wave sequence (owner decision D-003):
  left state → fast rarefaction (`rho` 1.0→0.74, `By` 1.0→0.66) → compound wave
  with `By` sign reversal near `x≈0.5` → post-compound plateau `rho≈0.5` →
  contact near `x≈0.68` → slow shock → right state `rho≈0.124`, `p≈0.10`,
  `By≈-1.0`. This matches Brio & Wu (1988) and the Athena Brio–Wu test.

## Accuracy is first order — by design, not a defect

The corrected solver keeps the historical scheme: piecewise-constant space,
forward-Euler time, nodal-averaged EMF, RT0 cell-`B` reconstruction. `minmod` /
`applyLimiter` exists in the source but is **still never called**; this is
intentional so that H1 compares a genuine first-order legacy baseline against the
second-order AMReX scheme.

Consequences, all consistent with a first-order Godunov + first-order CT method:

- CP-Alfvén: Tóth-style average relative `L1` over `(v⊥, vz, B⊥, Bz)` = **0.370**
  at 32×56; the transverse wave loses ~39 % amplitude over one return. Total
  magnetic energy is essentially conserved (ratio 0.994) because the constant
  background `B∥` dominates it.
- Field loop: magnetic energy retained is **≈ 0.17** of the initial after two
  domain wraps (legacy-scaled) / the Athena return time; the return-`B` relative
  `L1` is ≈ 1.0, i.e. the weak loop is almost fully diffused. Same qualitative
  outcome as the earlier Netgen `h=0.04` run.
- Brio–Wu on a triangular strip: the maximum within-`x`-bin standard deviation
  of the projected profile reaches ~0.20 (pressure) / ~0.10 (`rho`, `By`) near
  the strongest fronts. A strictly 1-D problem acquires this transverse scatter
  because the two triangle orientations per rectangle see asymmetric stencils.
  The AMReX Cartesian strip will not have it; the T08 comparison must therefore
  project area-weighted onto shared `x`-bins and report this spread.

## Numerical-delta ledger — addition

| Delta | Kind | Reason | Regression / evidence | Claim boundary |
|---|---|---|---|---|
| `--mesh-backend structured` in `run_legacy_corrected.py` | harness, not physics | Netgen unavailable on this host; a hashable dependency-free mesh makes the corrected runs reproducible and `dx`-comparable to AMReX | 4 manifests with `quality_gate = pass`; overlay SHA-256 unchanged; `git apply --check` clean each run | applies to the reproducibility of the run, not to the solver mathematics, which is untouched |

No change was made to the overlay or to any solver / geometry source file for
this slice.

## What is now closed vs. still open

Closed by this slice:

- `legacy_corrected` runs every canonical CPU case to completion, admissible,
  conservative and discretely divergence-free — the "no new canonical raw runs
  for the current overlay" gap from `LEGACY_CORRECTED_T05_STATE.md`.
- The corrected solver produces the correct Brio–Wu compound-wave structure.

Still open (unchanged, later phases):

- CUDA path is untouched and unverified; the runner still forces `gpu=false`.
  No CPU/GPU parity claim (T10/T12 territory).
- Brio–Wu quantitative branch-specific error table vs an independent
  compound-branch reference (T06). The branch is now chosen (D-003 =
  compound/non-regular); the reference solution and front-metric annotations are
  the T06 deliverable.
- The non-canonical convenience git worktree `legacy/mhd2d-corrected/` (branch
  `legacy_corrected_cpu_t05`, commit `516935f`) is a leftover from an earlier
  session and is untracked in `mhd-amrex`. The canonical reproducible path is
  `run_legacy_corrected.py`, which clones and patches from the immutable source.

## Reproduce

```sh
for case in brio_wu cp_alfven magnetic_loop_legacy_scaled magnetic_loop_athena; do
  python3 scripts/run_legacy_corrected.py \
    --source /Users/ivansamanov/Documents/MHD2D \
    --case "$case" --mesh-backend structured \
    --artifact-dir "benchmarks/raw/legacy_corrected/${case}_struct_t05_$(date +%s)" \
    --omp-threads 1
done
```
