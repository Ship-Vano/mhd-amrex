#!/usr/bin/env python3
"""Cheap structural gate for the AMReX CUDA execution path.

This test deliberately does *not* claim to run CUDA.  It catches regressions
that a CPU-only CI host can see: host-only launch loops, std::function problem
data, missing GPU boundary functor, or loss of CUDA target setup.

It also guards the other direction, which the CUDA port actually broke once:
porting the hot loop to `ParallelFor` removed the `#pragma omp parallel` that
distributed MFIter tiles across threads, and single-node OpenMP speedup
silently went from 2.7x to nothing.  No numerical test can see that -- results
stay bit-identical, only the wall clock changes -- so the pragma is asserted
structurally here, together with the `Gpu::notInLaunchRegion()` guard that
keeps it from spawning host threads around GPU kernel launches.
"""
from __future__ import annotations

import argparse
from pathlib import Path


def _function_body(text: str, qualified_name: str) -> str | None:
    """Тело функции от её сигнатуры до начала следующей функции верхнего уровня.

    Грубо, но достаточно: нужен ответ на вопрос «есть ли прагма именно в этой
    функции», а не разбор C++.
    """
    start = text.find(qualified_name)
    if start < 0:
        return None
    nxt = text.find("\nvoid MhdAmr::", start + len(qualified_name))
    nxt2 = text.find("\namrex::", start + len(qualified_name))
    ends = [e for e in (nxt, nxt2) if e > 0]
    return text[start:min(ends)] if ends else text[start:]


def require(text: str, needle: str, source: Path, failures: list[str]) -> None:
    if needle not in text:
        failures.append(f"{source.relative_to(source.parents[1])}: missing {needle!r}")



def device_lambda_methods_must_be_public(solver_text: str, header_text: str,
                                         failures: list[str]) -> None:
    """nvcc: extended __device__ lambda не может жить в private/protected методе.

    Ограничение CUDA 12.4 формулируется как "The enclosing parent function ...
    cannot have private or protected access within its class". На CPU-хосте оно
    невидимо: код компилируется обычным C++ и проходит все тесты. Проявляется
    только под nvcc, и к этому моменту позади уже полная сборка AMReX, поэтому
    инвариант проверяется здесь.
    """
    # Методы, в телах которых есть device-лямбда.
    needs_public: list[str] = []
    marks = [(m, solver_text.find(m)) for m in ("\nvoid MhdAmr::", "\namrex::Real MhdAmr::",
                                               "\nReal MhdAmr::", "\namrex::Long MhdAmr::")]
    starts = []
    for prefix in ("\nvoid MhdAmr::", "\namrex::Real MhdAmr::", "\nReal MhdAmr::",
                   "\namrex::Long MhdAmr::"):
        i = 0
        while True:
            i = solver_text.find(prefix, i)
            if i < 0:
                break
            name_start = i + len(prefix)
            name_end = solver_text.find("(", name_start)
            starts.append((i, solver_text[name_start:name_end].strip()))
            i = name_end
    starts.sort()
    for idx, (pos, name) in enumerate(starts):
        end = starts[idx + 1][0] if idx + 1 < len(starts) else len(solver_text)
        if "AMREX_GPU_DEVICE" in solver_text[pos:end]:
            needs_public.append(name)

    # Имена, объявленные в public-секции заголовка.
    public_names: set[str] = set()
    section = None
    for line in header_text.splitlines():
        stripped = line.strip()
        if stripped in ("public:", "protected:", "private:"):
            section = stripped[:-1]
            continue
        if section != "public" or "(" not in line:
            continue
        head = line.split("(")[0].strip()
        if head:
            public_names.add(head.split()[-1].lstrip("*&"))

    for name in sorted(set(needs_public)):
        if name not in public_names:
            failures.append(
                f"MhdAmr::{name} contains a __device__ lambda but is not public; "
                "nvcc rejects extended device lambdas in private/protected members")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-dir", required=True, type=Path)
    args = ap.parse_args()
    root = args.source_dir.resolve()
    solver = root / "src" / "MhdAmr.cpp"
    problems = root / "src" / "Problems.H"
    cmake = root / "CMakeLists.txt"
    failures: list[str] = []

    solver_text = solver.read_text(encoding="utf-8")
    problem_text = problems.read_text(encoding="utf-8")
    cmake_text = cmake.read_text(encoding="utf-8")
    if "LoopOnCpu" in solver_text:
        failures.append("src/MhdAmr.cpp still contains a LoopOnCpu launch")
    if "std::function<" in problem_text:
        failures.append("src/Problems.H uses std::function, which cannot be captured by CUDA lambdas")
    header_text = (root / "src" / "MhdAmr.H").read_text(encoding="utf-8")
    device_lambda_methods_must_be_public(solver_text, header_text, failures)
    require(solver_text, "GpuBndryFuncFab<ExtDirGpuFill>", solver, failures)
    require(solver_text, "Gpu::DeviceScalar<Long>", solver, failures)
    require(solver_text, "ReduceOps<ReduceOpMin>", solver, failures)
    require(solver_text, "ReduceOps<ReduceOpMax>", solver, failures)
    require(solver_text, "RunOn::Gpu", solver, failures)
    require(problem_text, "struct DeviceProblem", problems, failures)
    require(problem_text, "MHD_HD", problems, failures)
    require(cmake_text, "setup_target_for_cuda_compilation(mhd2d)", cmake, failures)

    # Горячие циклы обязаны сохранять распараллеливание по тайлам, а прагма --
    # GPU-защиту. Проверяется на уровне функции, иначе тест зелёный, когда
    # прагма есть где-то ещё в файле.
    for func in ("MhdAmr::ComputeFluxesAndEmf", "MhdAmr::ApplyUpdates"):
        body = _function_body(solver_text, func)
        if body is None:
            failures.append(f"src/MhdAmr.cpp: cannot locate {func}")
            continue
        if "#pragma omp parallel" not in body:
            failures.append(
                f"{func} lost '#pragma omp parallel': MFIter tiles no longer "
                f"spread across threads and single-node OpenMP scaling dies "
                f"without any numerical test noticing")
        elif "Gpu::notInLaunchRegion()" not in body:
            failures.append(
                f"{func}: '#pragma omp parallel' must be guarded by "
                f"'if (Gpu::notInLaunchRegion())', otherwise a GPU build with "
                f"OpenMP spawns host threads that each launch kernels")

    if failures:
        for failure in failures:
            print(f"GPU portability regression: {failure}")
        return 1
    print("OK: static CUDA-portability constraints hold; GPU execution still needs a CUDA host.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
