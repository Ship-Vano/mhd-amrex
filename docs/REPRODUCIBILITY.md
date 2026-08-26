# T02 reproducibility foundation

## Build and test

```sh
cmake --preset release
cmake --build --preset release
ctest --preset release
```

The `release` preset records the CPU Release configuration with MPI and
OpenMP. `ctest` runs direct kernel invariants, bounded standalone Brio--Wu and
Alfvén regressions, strict JSON validation, and a deterministic-manifest
repeat check. The standalone tests are regressions, not independent numerical
verification; their required evidence is still tracked in T04.

## Strict configurations

Only documented top-level sections and keys are accepted. Unknown keys,
wrong types, invalid paired periodic boundaries and physically invalid basic
parameters stop before AMReX allocates the mesh. The intentionally malformed
`config.rejects_unknown_key` CTest verifies rejection with `unknown key`; its
test wrapper passes only when the validator itself rejects the input.

## Compute-only runs

Set `output.write_plotfiles` to `false` to disable initial, periodic and final
plotfiles. `inputs/benchmark_smoke.json` is a small example. This only removes
non-critical plotfile I/O; it is not yet a benchmark harness and must not be
used to claim performance before T09.

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
