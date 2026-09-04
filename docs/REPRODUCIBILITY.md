# T02 reproducibility foundation

## Build and test

```sh
cmake --preset release
cmake --build --preset release
ctest --preset release
```

The `release` preset records the CPU Release configuration with MPI and
OpenMP. Other presets: `cpu-release` (serial), `cpu-debug` (AMReX assertions
and bound checking), `mpi-release`, `profile` (`AMReX_TINY_PROFILE`),
`hdf5-release`, and `cuda-release` — the last is declared but is **not**
buildable on this workstation, which has no NVIDIA GPU.

`ctest` runs 21 tests: kernel invariants; bounded standalone regressions;
strict JSON validation; the canonical-case matrix (constant state, Orszag--Tang,
rotor, MHD blast) checked against literature `rho`/`p` ranges; AMR conservation
with and without reflux, and across regrids; MPI decomposition parity; the
independent Brio--Wu cross-validation; the ADR 0001 kernel-purity check; and a
deterministic-manifest repeat check.

Note what the standalone regressions are and are not: they pin behaviour
against previously recorded values. Independent numerical evidence comes from
`briowu.independent_reference`, which compares against a scheme sharing no code
with the solver (see `docs/BRIOWU_T06.md`).

## Strict configurations

Only documented top-level sections and keys are accepted. Unknown keys,
wrong types, invalid paired periodic boundaries and physically invalid basic
parameters stop before AMReX allocates the mesh. The intentionally malformed
`config.rejects_unknown_key` CTest verifies rejection with `unknown key`; its
test wrapper passes only when the validator itself rejects the input.

## Compute-only runs

Set `output.write_plotfiles` to `false` to disable initial, periodic and final
plotfiles. `inputs/benchmark_smoke.json` is a small example. This only removes
non-critical plotfile I/O.

Timing is measured with `scripts/benchmark.py` (machine manifest, at least one
discarded warm-up and at least five timed repeats, median and MAD) and
`scripts/scaling.py` for OpenMP/MPI speedup; see `docs/T09_TIMING.md`. All such
numbers come from a single unpinned workstation and are recorded as
**diagnostic**: they may not be quoted as scaling results. The multi-node and
cluster campaigns remain blocked on access (decision D-005).

The final plotfile is now skipped when the scheduled output already wrote the
same final state, preventing the previous `.old.*` directory side effect.

## Run manifests

```sh
python3 scripts/write_run_manifest.py \
  --config inputs/benchmark_smoke.json \
  --output benchmarks/raw/example-manifest.json
```

The manifest captures commit, dirty state, canonical JSON and its SHA-256,
plus stable local environment fields. It intentionally excludes wall-clock
time so that identical invocations are byte-equivalent; timings and repetitions
will be added only by the benchmark phase. Store raw results in
`benchmarks/raw/`, aggregates in `benchmarks/summary/`, and figures in
`benchmarks/figures/`.

## Regenerating the report

```sh
cmake --build --preset release
sh scripts/regen_report_data.sh          # all figure inputs, deterministic
python3 scripts/make_report_figures.py   # -> docs/figures/data/*.dat
cd docs && latexmk -pdf report.tex       # -> 19 pages, 0 overfull boxes
```

`scripts/make_report_figures.py` refuses to record an ephemeral (`/tmp`) source
path in a figure header, so a figure whose provenance cannot be reproduced
cannot silently enter the report.

## Legacy solver runs

The historical tree at `/Users/ivansamanov/Documents/MHD2D` is never modified.
It is pinned by the annotated tag `legacy_vkr/9d0f60e` and by the archive
checksum recorded in `docs/LEGACY_BASELINE.md`, both verified by the runner on
every invocation.

```sh
python3 scripts/run_legacy_corrected.py --source <MHD2D> --case rotor \
    --mesh-backend structured --artifact-dir benchmarks/raw/legacy_corrected/rotor_128
```

For irregular meshes use `--mesh-backend netgen --netgen-python <venv>/bin/python`;
the pinned environment is described in `docs/MESH_PIPELINE.md`.
