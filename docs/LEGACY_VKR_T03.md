# T03: controlled reproduction of `legacy_vkr`

**Hypothesis.** The clean source commit `9d0f60ea8576fac5d6f28c4dec142236d76131d6`
can be built and run reproducibly on this macOS CPU host without changing its
mathematics, provided each new input and output is placed in an isolated run
directory.

## Scope and provenance

The owner has designated the latest commit above as the official `legacy_vkr`
source and authorised new controlled runs.  `scripts/run_legacy_vkr.py` refuses
a different or dirty source, makes a detached clone, archives that source and
records SHA-256 values for the archive, mesh, JSON, console log and VTU states.
It does not modify the source directory.  The generated raw directories are
under `benchmarks/raw/legacy_vkr/` and deliberately ignored by Git because VTU
files are reproducible heavy artifacts.

The historical Netgen mesh/config/result bundle on `Elements` remains
`legacy_unattributed`: its provenance to this commit is not proven.  The new
generator is intentionally a dependency-free structured triangular mesh
generator, **not** a claim that it reproduces that historical Netgen grid.  It
uses the exact simple text format consumed by `NetGeometry.cpp`; its two
triangles per rectangle are counter-clockwise and its hash is in every run
manifest.  This separates a scientifically reproducible new run from an
unsupported retrospective claim about the VCR results.

## Canonical controlled cases

| Run id | Code task | Domain / mesh | Physical time | Scientific comparator and acceptance observation |
|---|---:|---|---:|---|
| `brio_wu` | 1 | `[0,1] x [-0.01,0.01]`, 128x4 rectangles | 0.1 | The standard states, `Bx=0.75`, `gamma=2` match [Brio–Wu as specified by Athena](https://www.astro.princeton.edu/~jstone/Athena/tests/brio-wu/Brio-Wu.html). Confirm finiteness and retain the profile; a branch-specific error table remains T06 because the owner has not selected regular vs. compound/non-regular reference. |
| `cp_alfven` | 8 | `[0,2/sqrt(3)] x [0,2]`, 32x56 rectangles | 1.0 | [Athena CP Alfvén](https://www.astro.princeton.edu/~jstone/Athena/tests/cp-alfven-wave/cp-alfven.html) identifies it as an exact nonlinear smooth solution. The domain makes `x cos(pi/6)+y sin(pi/6)` periodic in both directions; compare returned state to the analytic state with area-weighted L1/L2/Linf errors. |
| `magnetic_loop` | 9 | `[-1,1] x [-0.5,0.5]`, 64x32 rectangles | 2.0 | [Athena field-loop test](https://www.astro.princeton.edu/~jstone/Athena/tests/field-loop/Field-loop.html) specifies `rho=p=1`, `gamma=5/3`, `A=1e-3`, `R0=0.3` and magnetic-energy decay as the quantitative diagnostic. This legacy run has velocity `(2,1)`, so one wrap is `t=1` and code-imposed `t=2` is two wraps. |

The test selection is also supported by the [Athena test catalogue](https://www.astro.princeton.edu/~jstone/Athena/tests/), which places Brio–Wu, nonlinear CP Alfvén and field-loop advection in its MHD verification suite.  The HLLD solver itself is correctly attributed to [Miyoshi & Kusano (2005)](https://doi.org/10.1016/j.jcp.2005.02.017), rather than using a visual match as a correctness claim.

## Scientific corrections and non-claims

- `taskType=2` is not the canonical CP setup: its code uses a different,
  effectively y-directed wave with an explicit `4 pi` convention.  T03 uses
  `taskType=8`.
- `taskType=8` computes `alpha=pi/6`, i.e. **30 degrees**.  Its comment says
  `atan(0.5)`, approximately 26.6 degrees; the calculation, not the comment,
  defines the run.
- The legacy code itself forces Brio–Wu CFL to `0.9`, while the VCR describes
  `0.1`; the controlled run therefore cannot be called a reproduction of the
  VCR convergence table.  This is an immutable numerical delta for T05, not a
  value adjusted in `legacy_vkr`.
- Printed `computeDivergence()` is not accepted as `max(abs(div B))`: the code
  takes a signed maximum.  The field-loop metric here is magnetic energy and
  state finiteness only; a valid divergence norm belongs to T05.
- Athena's illustrated field-loop geometry has a different inclined unit-speed
  flow.  The selected domain is deliberately compatible with the legacy hard
  coded velocity and final time, so it verifies the same physical advection
  mechanism without falsely claiming bitwise identity with Athena's plot.

## Reproduction command

```sh
python3 scripts/run_legacy_vkr.py \
  --source /Users/ivansamanov/Documents/MHD2D \
  --case cp_alfven \
  --artifact-dir benchmarks/raw/legacy_vkr/cp_alfven_n32x56
```

Run it once for each case in the table.  The runner accepts success only when
the process exits normally, the console's `Final time` equals the configured
physical time, and the final VTU has finite conservative state; a detected NaN
can otherwise cause this legacy code to exit with status zero.  `tmpres_0.vtu`
is retained as a raw file but is not
labelled an initial state, because an early failure at iteration zero overwrites
that filename with the final dump.

It requires CMake, Python 3, Git and
GCC with OpenMP; on this host the compiler is Homebrew `g++-15`.  The script
does not treat wall time as a benchmark and does not create a performance claim.

## T03 gate

`PASS` requires all three commands to complete and the three generated
`manifest.json` files to identify the official source/archive, exact mesh and
config, and final-state diagnostics.  The decision on the Brio–Wu Riemann
branch is explicitly outside T03, so T03 does not publish a branch-dependent
accuracy claim.

## Controlled-run result (2026-08-26)

**Status: FAILED.** The build and provenance portions succeeded, but the
scientific reproduction gate did not.  All runs used GCC 15.2, OpenMP 4.5 and
the same source archive SHA-256
`afb51261d3d3b36760573afa4a12e7567d49bcfdc3b15c504ec9604df15ff440`.

| Case | Target / reached time | Iterations | Result |
|---|---:|---:|---|
| Brio–Wu | `0.1 / 0.007406013972` | 68 | NaN; `rho_min=-0.176755`, `p_min=-0.073292362282`. |
| CP Alfvén | `1.0 / 0.003143912188` | 0 | NaN on the first update; no analytic return error can be defined. |
| Magnetic loop | `2.0 / 0.0008860063359` | 0 | NaN on the first update; no magnetic-energy decay can be measured. |

The versioned compact record is
[`benchmarks/summary/legacy_vkr_t03_controlled_runs.json`](../benchmarks/summary/legacy_vkr_t03_controlled_runs.json);
the corresponding ignored raw directories contain each full manifest, mesh,
JSON, console log and VTU hashes.  The exit status of the executable is zero in
all three cases, which is why the runner's physical-time and finiteness checks
are essential.

The evidence warrants neither comparison with the VCR convergence numbers nor
a claim that the three legacy cases reproduce the literature.  The next
allowed numerical work is T05: preserve `legacy_vkr`, create a separate
`legacy_corrected` profile, and add regressions before diagnosing and fixing
the edge-state/NaN path, positivity failure and signed-divergence diagnostic.
