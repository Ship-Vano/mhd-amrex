# T04 slice: AMReX kernel correctness + positivity-preserving HLLD fallback

Date: 2026-08-31. Status: **PARTIAL → advanced** on the T04 gate. This slice
adds the kernel unit/regression coverage and the positivity safety net that were
missing; the CT `div B` norm artifact (NUM-004) and the full canonical AMR suite
remain for later T04/T07 work.

## Motivation

`legacy_vkr` blows up on every canonical Riemann problem because a degenerate or
floored intermediate state turns into NaN and propagates (T03). The AMReX
`hlld_flux` had the same structural exposure: no validity predicate, no
positivity-preserving fallback, and `cons_to_prim` floored `rho`/`p` silently
with no counter — exactly the pattern `AGENTS.md` invariant 1 forbids
("Нельзя считать молчаливый floor доказательством корректности").

## Changes

### `src/kernels/Hlld.H`

- `mhd_wave_speeds()` — the Davis-type extremal speed estimate factored out
  (`SL = min(uL,uR) − max(cf)`, `SR = max(uL,uR) + max(cf)`); wider than
  Einfeldt, so an HLL average built on it is positivity-preserving under CFL.
  Numerically identical to the previous inline estimate.
- `hll_flux()` — the single-state HLL (HLLE) flux. Reused as the HLLD fallback
  and available as the N0 ablation (piecewise-constant + Euler + HLLE).
- `hlld_flux_raw()` — the previous HLLD body, now returning `bool`. Returns
  `false` (→ caller falls back) when:
  - an input state has `rho ≤ small_rho` or `p ≤ small_pres` (post-floor path);
  - the HLL-averaged fan density `rho_hll ≤ small_rho` (near-vacuum fan,
    Einfeldt positivity check);
  - either outer star density `rho*_{L,R} ≤ small_rho`;
  - any returned flux component is non-finite.
- `hlld_flux()` — thin wrapper: run `hlld_flux_raw`; on `false` overwrite with
  `hll_flux` and, if a non-null `int* n_fallback` was passed, increment it.
  **For every admissible state the returned flux is bit-identical to before** —
  verified: the standalone Alfvén `L1` errors and Brio–Wu `rho_min`/`p_min` are
  unchanged to the last digit, step counts identical.

The GPU port (T12) must replace the `int*` counter with an atomic / reduction.

### Diagnostics

- `tests/standalone_verify.cpp`: file-scope `g_hlld_fallbacks`, threaded into
  both HLLD calls, printed on the `DONE` line.
- `src/MhdAmr.{H,cpp}`: `hlld_fallbacks_` member accumulated (OMP reduction) in
  `ComputeFluxesAndEmf`; new `CountNonPositiveCells()` scans the hierarchy for
  `rho ≤ small_rho` or `p ≤ small_pres`. Both printed every `diag_int` steps and
  at end of run (`Evolve finished: … hlld_fallbacks=N, nonpositive_cells=M`),
  reduced across ranks.

## New CTest coverage

| test | contents |
|---|---|
| `kernel.numerics` (`tests/kernel_numerics.cpp`) | 2000-sample prim↔cons round trip and `pressure_from_cons`; `fast_speed` hydro / normal-field `max(a,ca)` / transverse `√(a²+cA²)` limits + monotonicity + `cf ≥ a, cf ≥ |ca_n|`; HLLD consistency `F(q,q)=F_phys` over 3000 states (fallback must stay 0); HLL consistency; HLLD supersonic-left/right branches; `Bn=0`, Brio–Wu, rotational-discontinuity degeneracies; 20000-sample finiteness sweep over `ρ,p ∈ [10⁻⁴,…]`, `|v|≤15`, `|B|≤8`; **positivity fallback**: smooth→0 fallbacks, `p` at floor / `p≤0` / near-vacuum → `raw` rejects and wrapper exposes exactly one fallback with finite flux; limiter constant-preservation / extremum→0 / odd symmetry / TVD bound / linear-profile exactness; `corner_emf` Balsara–Spicer & Gardiner–Stone ε=2 forms + constant state; SSP-RK2 amplification polynomial `1+a+a²/2` and observed ODE order 2.004 |
| `standalone.alfven_order` (`tests/check_alfven_order.py`) | CP-Alfvén `L1(sum4)`: N=64 `2.365e-3`, N=128 `6.255e-4`, **observed order 1.919 ≥ 1.8** |
| `standalone.briowu`, `standalone.alfven32` | now also assert `hlld_fallbacks=0` |

`ctest --preset release`: **9/9 pass**.

## Canonical-case fallback / positivity audit

All runs on the default inputs, current commit:

| case | driver | steps | `hlld_fallbacks` | `nonpositive_cells` | max\|divB\| |
|---|---|---:|---:|---:|---|
| Brio–Wu | AMReX `mhd2d` | 391 | 0 | 0 | 0 |
| Brio–Wu | standalone | 488 | 0 | — | 0 |
| CP-Alfvén N=32…128 | standalone | 41…322 | 0 | — | ≤ 6e-12 |
| Orszag–Tang N=96 | standalone | 302 | 0 | — | 3.0e-13 |
| rotor N=96 | standalone | 126 | 0 | — | 7.7e-13 |
| rotor (AMR, 2 levels) | AMReX `mhd2d` | 470 | 0 | 0 | 1.3e-11 |

The fallback never fires on canonical data — it is strictly the safety net for
the floored / degenerate path, which `kernel.numerics` exercises deterministically.

## Still open in T04 / adjacent

- NUM-004: a CTest artifact carrying the normalized uniform-CT `div B` norm
  `dx·‖div B‖∞ / max(‖B‖∞, Bref)` (currently only printed).
- Characteristic vs primitive limiting, corner-EMF 1-D consistency proof.
- Gas reflux for AMR conservation (P0) — T07, not touched here.
- Threading a per-cell floor *event* counter through `cons_to_prim` itself
  (current `CountNonPositiveCells` counts post-step cells at/below the floor,
  which is the physically meaningful "admissible run ⇒ 0" quantity).
