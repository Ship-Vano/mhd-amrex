# Are the corrected solvers physically adequate? — comparison with the literature

Date: 2026-08-31. Question from the owner: do the corrected `legacy_corrected`
and AMReX solvers give results consistent with the graphs and test tables in
(a) I. P. Shamanov's bachelor thesis (VKRB, `VKRB_Shamanov_FINITE.pdf`),
(b) E. N. Avdeeva & V. V. Lukin, *Divergence-free finite-difference method for
2D ideal MHD*, JPCS 1336 (2019) 012026, and
(c) the Athena test catalogue / primary literature (Brio–Wu 1988, Tóth 2000,
Miyoshi–Kusano 2005)?

**Short answer: yes, they are physically adequate.** Every canonical case that
runs to completion reproduces the qualitative wave structure and the
quantitative ranges reported in those sources, with the differences that are
*expected* from the scheme choices (see below). Two caveats and one real gap are
listed at the end.

Method notes:

- The literature figures are **qualitative** (no digitised data was published),
  so this is a feature-and-range comparison — extrema, front locations,
  colour-scale bounds, divergence growth, oscillation presence — not an
  L-error against a picture.
- `legacy_corrected` and the two literature codes are the **same scheme class**:
  first-order forward Euler, HLLD, Stokes-theorem Faraday update, RT0 cell-`B`.
  Avdeeva–Lukin's extremal-speed estimate `SL,R = min(uL,uR) ∓ max(cfL,cfR)`
  is *identical* to the one in `src/kernels/Hlld.H::mhd_wave_speeds`.
- The AMReX solver is **second order** (MUSCL-MC + SSP-RK2 + constrained
  transport). It should — and does — beat the first-order literature results at
  equal resolution; that gap is the subject of hypothesis H1, not a
  discrepancy.
- Meshes differ: `legacy_corrected` ran on a **structured right-triangle** mesh
  (reproducible, no Netgen), Avdeeva–Lukin used **equilateral** triangles, the
  VKRB used Netgen. So no claim of a pixel match, and the VKRB's own error
  tables (its Tables 1–4) stay `historical_unreproduced` — consistency of order
  and magnitude is shown, not bit reproduction.

## `legacy_corrected` (first order) vs the literature

| Case | This work | VKRB / Avdeeva–Lukin | Verdict |
|---|---|---|---|
| **Brio–Wu** (128×4, 400×8, γ=2, `Bx`=0.75, `t`=0.1, CFL 0.1) | correct **compound / non-regular** sequence: fast rarefaction → `By` sign-flip compound wave near x≈0.5 → post-compound plateau ρ≈0.5 → contact → slow shock → right state; finite, ρ_min≈0.11, scaled `‖div B‖` ≤ 6·10⁻¹⁵ | VKRB Fig. 6 & Table 2 ("составная волна в форме пика", order ≈ 1st); Avdeeva Fig. 4 ("shock front 2 cells, contact 3–4 cells") | ✅ structure matches |
| **CP-Alfvén 30°** (32×56, `t`=1, CFL 0.1) | Tóth-style average relative `L1` = **0.370**; isolated-1D HLLD `L1` = 0.43 / 0.24 / 0.13 at N=32/64/128 → observed order ≈ 1.2–1.3 | VKRB Table 3 (2D): `L2` = 0.202 / 0.087 / 0.035, order ≈ 1.22 / 1.29 | ✅ same first-order class, same ~30–40 % error at N≈32 |
| **Field loop** (64×32 & 64×37, `t`=2 / 4.62, CFL 0.1) | heavy diffusion, `E_B(t)/E_B(0)` ≈ 0.17; `div B`: **1.4·10⁻¹⁷ → 4.8·10⁻¹⁶** (scaled 2.2·10⁻¹⁵) | VKRB §4.1.4: "**свойство высокой диссипации**"; `divₕ(0)` = 5.92·10⁻¹⁷ → `divₕ(2)` = 2.02·10⁻¹⁵ | ✅✅ divergence numbers match to a factor of ~2; "high dissipation" reproduced |
| **Rotor** (128×128, γ=1.4, `Bx`=5/√4π, `t`=0.15, CFL 0.5) | completes; x=y diagonal is **symmetric**: ambient ρ=p=1 → compression hump p≈1.46 → **ρ spike** (≈4.4 binned / 8.3 raw) → central **p evacuation to 0.06**; `div B`: 5.3·10⁻¹⁶ → 4.1·10⁻¹² (scaled 3.3·10⁻¹⁵); 0 HLLE fallbacks | VKRB §4.1.2: `divₕ` 5.4·10⁻¹³ → 5.0·10⁻¹¹ ("возросло на 2 порядка … ошибки округления"); Figs. 19–22 ρ scale ≈ [0.55, 11], two spikes ≈ 10, p dip ≈ 0.05–0.1, humps ≈ 1.55; Avdeeva Figs. 5–6 ρ ∈ [0.571, 10.60] | ✅ structure + "div B grows ≈2 orders, stays roundoff" reproduced; the ρ peak is more diffused than the literature because N=128 ≪ their N≈400 |
| **Orszag–Tang** (128×128, `t`=0.5) | **positivity guard fires at `t≈0.4314`, element 1371** (centroid ≈ (0.354, 0.044), near the periodic y-boundary) — the corrected code refuses to emit a non-physical state. **Identical failure point at CFL 0.5 and at CFL 0.2** (iterations 1730 vs 4324, same `t`, same element) → it is *not* a timestep-stability issue but first-order under-resolution of an OT low-β region on this coarse structured mesh. | VKRB §4.1.3 & Avdeeva §4.3 ran OT to `t`=0.5 on **N≈400 / 800** (≈3–6× finer), CFL 0.2 | ⚠️ `legacy_corrected` does not complete OT at N=128; needs a finer mesh. The historical code has no positivity check (only a NaN test), so its published OT figures may contain silently-inadmissible cells. Also see the CFL-config gap below. |

Historical `legacy_vkr` (uncorrected, commit `9d0f60e`) NaNs on Brio–Wu,
CP-Alfvén and the field loop (T03). The corrected overlay makes all of these
physically admissible **and** brings them into agreement with the literature —
i.e. the corrections (HLLD flux, HLLE fallback, conversions, CFL bound,
edge↔ghost maps, RT0 sign, CT orientation, post-CT energy) *restored* the
physics rather than distorting it.

## AMReX (second order) vs the literature

| Case | This work | Literature | Verdict |
|---|---|---|---|
| **Brio–Wu** (N3, N up to 2048) | canonical compound-branch profile; ρ_min ≈ 0.115, `max\|div B\|` = 0, 0 fallbacks; ≈ 15× lower `L1` than `legacy_corrected` at equal `dx` | Avdeeva Fig. 4; Brio & Wu 1988; Athena Brio–Wu | ✅ + expected 2nd-order gain over the 1st-order references |
| **CP-Alfvén 30°** | observed order **1.92** (N=64→128); `L1(sum4)` = 6.25·10⁻⁴ at N=128 | VKRB Table 3 (1st order): 0.035 at N=128 → AMReX ≈ **56×** more accurate | ✅ second order confirmed; the large margin over the 1st-order refs is exactly what H1 predicts |
| **Orszag–Tang** (N=128, `t`=0.5, CFL 0.4) | ρ ∈ **[0.091, 0.492]**, p ∈ **[0.028, 0.505]**; `max\|div B\|` = 5·10⁻¹³, 0 fallbacks. y=0.3125 pressure slice: left plateau ≈ 0.21 → dip at x≈0.30 to ≈ 0.07 → plateau ≈ 0.10 → **peak ≈ 0.26 at x≈0.58** → deep valley x≈0.73–0.92 at ≈ 0.04 → recovery to ≈ 0.20 | Avdeeva Fig. 7 ρ ∈ [0.087, 0.489], Fig. 25 p ∈ [0.028, 0.49], **Fig. 8** p(y=0.3125); VKRB Figs. 23–26 | ✅✅ colour-scale ranges match to ~2 %; the slice reproduces Avdeeva Fig. 8 **feature-by-feature** (same count, same x-locations, same magnitudes), no spurious oscillations |
| **Rotor** (N=128, `t`=0.15) | ρ ∈ **[0.66, 11.55]**, p ∈ **[0.039, 1.73]**; `max\|div B\|` = 1.2·10⁻¹², 0 fallbacks. x=y diagonal: perfectly **symmetric**, ambient ρ=p=1, compression humps p≈1.5, two ρ spikes (11.5 raw), central p evacuation to **0.04** | Avdeeva Fig. 5 ρ ∈ [0.571, 10.60], Fig. 6; VKRB Fig. 19 ρ ∈ [0.55, 11], Figs. 20/22 | ✅✅ ranges and diagonal structure match; the ρ spike slightly overshoots 10 (2nd-order sharpening of the ρ=10 disk edge) |
| **`div B` on every case** | scaled norm at roundoff, 0 HLLE fallbacks, 0 non-positive cells | both papers' central claim: the scheme keeps `div B` at roundoff | ✅✅ |

The AMReX `hlld_flux` HLLE fallback (added this session) never fired on any of
these — consistent with Avdeeva–Lukin and the VKRB reporting no positivity
trouble on the standard suite.

## Caveats and one real gap

1. **Brio–Wu quantitative branch-specific error table** vs an *independent*
   exact/near-exact compound-branch Riemann solution is still **T06**. The
   branch is chosen (D-003 = compound/non-regular); the current cross-solver
   numbers use a provisional grid-converged AMReX reference
   (`docs/CROSS_SOLVER_BRIOWU_1D.md`).
2. **VKRB CPU/GPU timing table** (Intel i5-4440 + GTX 960, Table 4:
   765 s → 539 s etc., ≈1.4× "speed-up") is **not reproduced** — different
   hardware, and the legacy GPU path is disabled (`gpu=false`). Our diagnostic
   timings are on a different machine and are not a benchmark.
3. **Gap — `legacy_corrected` task types 4 (rotor) and 5 (Orszag–Tang) still
   hard-code `cflNum = 0.5` and ignore the config `cfl`.** The T05 "make CFL
   explicit" correction was applied only to Brio–Wu / CP-Alfvén / field loop
   (task types 1 / 8 / 9). This is cosmetic for the OT case above — a scratch
   patch honouring `configuredCfl` gave the *same* failure at CFL 0.2 — but it
   still belongs in the overlay with a regression as a follow-up T05 item, so
   that rotor/OT CFL is a documented config value rather than a hidden constant.

4. **`legacy_corrected` Orszag–Tang has not been run to completion at any
   resolution here.** N=128 structured fails at `t≈0.43`; a run at N≈400 (to
   match the literature) was not attempted this session — it is the natural
   next legacy V&V step.

## Reproduce

```sh
# AMReX standalone (2nd order)
./build/release/mhd2d_verify ot 128        # -> out_ot.csv
./build/release/mhd2d_verify rotor 128     # -> out_rotor.csv
python3 scripts/slice_compare_2d.py --case orszag_tang --amrex-csv out_ot.csv --out /tmp/ot.json
python3 scripts/slice_compare_2d.py --case rotor --amrex-csv out_rotor.csv --gamma 1.4 --out /tmp/rotor.json

# legacy_corrected (1st order) — rotor via the worktree binary on a structured mesh
python3 -c "import sys;sys.path.insert(0,'scripts');from legacy_vkr_mesh import write_rectangular_tri_mesh;from pathlib import Path;write_rectangular_tri_mesh(Path('mesh.txt'),0,1,0,1,128,128)"
# taskType 4 (rotor) / 5 (OT), see legacy/mhd2d-corrected
```
