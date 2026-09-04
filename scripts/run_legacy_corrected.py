#!/usr/bin/env python3
"""Run one isolated legacy_corrected case with a versioned patch overlay.

The official legacy tree remains read-only. This runner clones its exact
commit, applies the corrected-physics overlay, creates an irregular Netgen
mesh, builds, executes CTest and then executes the solver. It retains a
manifest even if a numerical quality gate fails, so a failure is evidence
rather than a silently discarded run.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import math
import os
import platform
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
OFFICIAL_COMMIT = "9d0f60ea8576fac5d6f28c4dec142236d76131d6"

# The two field-loop cards are intentionally separate. The first is a
# canonical Athena-compatible return geometry; the second retains the scaled
# historical velocity only as an internal legacy regression, not as Athena.
CASES: dict[str, dict[str, Any]] = {
    "brio_wu": {
        "task_type": 1,
        "final_time": 0.1,
        "gamma": 2.0,
        "domain": (0.0, 1.0, -0.01, 0.01),
        "maxh": 0.0025,
        # Dependency-free structured backend: two CCW right triangles per
        # rectangle.  These are the T03 legacy_vkr resolutions so that the
        # corrected and historical runs share an identical mesh hash.
        "structured_resolution": (128, 4),
        "grading": 0.3,
        # Netgen's interior remains irregular.  Matching unperturbed boundary
        # nodes are required by legacy's coordinate-based periodic pairing.
        "boundary_jitter": 0.0,
        "cfl": 0.1,
        "diagnostic_case": "brio_wu",
        "profile_projection": True,
    },
    "cp_alfven": {
        "task_type": 8,
        "final_time": 1.0,
        "gamma": 5.0 / 3.0,
        "domain": (0.0, 2.0 / math.sqrt(3.0), 0.0, 2.0),
        "maxh": 0.04,
        "structured_resolution": (32, 56),
        "grading": 0.3,
        "boundary_jitter": 0.0,
        "cfl": 0.1,
        "diagnostic_case": "cp_alfven",
        "profile_projection": False,
    },
    "magnetic_loop_athena": {
        "task_type": 9,
        "final_time": 8.0 / math.sqrt(3.0),
        "gamma": 5.0 / 3.0,
        "domain": (-1.0, 1.0, -1.0 / (2.0 * math.cos(math.pi / 6.0)),
                   1.0 / (2.0 * math.cos(math.pi / 6.0))),
        "maxh": 0.04,
        "structured_resolution": (64, 37),
        "grading": 0.3,
        "boundary_jitter": 0.0,
        "cfl": 0.1,
        "diagnostic_case": "magnetic_loop",
        "profile_projection": False,
        "field_loop_u": math.sin(math.pi / 3.0),
        "field_loop_v": math.cos(math.pi / 3.0),
        "field_loop_radius": 0.3,
        "field_loop_amplitude": 1.0e-3,
    },
    # Вращающийся цилиндр (Tóth) и вихрь Орзага-Танга. До этого обе карты не
    # имели воспроизводимого раннера: CFL и finalTime были зашиты в исходник,
    # поэтому конфиг игнорировался. Обе теперь управляются манифестом.
    "rotor": {
        "task_type": 4,
        "final_time": 0.15,
        "gamma": 1.4,
        "domain": (0.0, 1.0, 0.0, 1.0),
        "maxh": 0.01,
        "structured_resolution": (128, 128),
        "grading": 0.3,
        "boundary_jitter": 0.0,
        "cfl": 0.5,
        "diagnostic_case": "rotor",
        "profile_projection": False,
    },
    "orszag_tang": {
        "task_type": 5,
        "final_time": 0.5,
        "gamma": 5.0 / 3.0,
        "domain": (0.0, 1.0, 0.0, 1.0),
        "maxh": 0.01,
        "structured_resolution": (128, 128),
        "grading": 0.3,
        "boundary_jitter": 0.0,
        "cfl": 0.5,
        "diagnostic_case": "orszag_tang",
        "profile_projection": False,
    },
    "magnetic_loop_legacy_scaled": {
        "task_type": 9,
        "final_time": 2.0,
        "gamma": 5.0 / 3.0,
        "domain": (-1.0, 1.0, -0.5, 0.5),
        "maxh": 0.04,
        "structured_resolution": (64, 32),
        "grading": 0.3,
        "boundary_jitter": 0.0,
        "cfl": 0.1,
        "diagnostic_case": "magnetic_loop",
        "profile_projection": False,
        "field_loop_u": 2.0,
        "field_loop_v": 1.0,
        "field_loop_radius": 0.3,
        "field_loop_amplitude": 1.0e-3,
    },
}


def invoke(command: list[str], *, cwd: Path | None = None,
           environment: dict[str, str] | None = None,
           stdout: Any | None = None) -> subprocess.CompletedProcess[str]:
    if stdout is None:
        captured_stdout: Any = subprocess.PIPE
        captured_stderr: Any = subprocess.PIPE
    else:
        captured_stdout = stdout
        captured_stderr = subprocess.STDOUT
    return subprocess.run(
        command,
        cwd=cwd,
        env=environment,
        text=True,
        stdout=captured_stdout,
        stderr=captured_stderr,
        check=False,
    )


def command_text(command: list[str]) -> str:
    result = invoke(command)
    if result.returncode:
        return "unavailable"
    return (result.stdout or result.stderr or "").strip()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def copy_if_present(source: Path, destination: Path) -> bool:
    if not source.is_file():
        return False
    shutil.copy2(source, destination)
    return True


def read_json_if_present(path: Path) -> dict[str, Any] | None:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (FileNotFoundError, json.JSONDecodeError):
        return None


def parse_solver_log(text: str) -> dict[str, Any]:
    result: dict[str, Any] = {}
    final = re.search(r"Final time = ([^;]+); iterations = (\d+)", text)
    if final:
        result["reported_final_time"] = float(final.group(1))
        result["iterations"] = int(final.group(2))
    initial_div = re.search(
        r"Init max\|magnetic flux residual\| = ([^;]+); max\|div B\| = ([^;]+); "
        r"global scaled magnetic-flux imbalance \(B_ref = ([^)]+)\) = ([^\s]+)",
        text,
    )
    if initial_div:
        result["initial_divergence"] = {
            "max_flux": float(initial_div.group(1)),
            "max_abs": float(initial_div.group(2)),
            "reference_field": float(initial_div.group(3)),
            "max_scaled": float(initial_div.group(4)),
        }
    final_div = re.search(
        r"Final max\|magnetic flux residual\| = ([^;]+); max\|div B\| = ([^;]+); "
        r"global scaled magnetic-flux imbalance = ([^\s]+)",
        text,
    )
    if final_div:
        result["final_divergence"] = {
            "max_flux": float(final_div.group(1)),
            "max_abs": float(final_div.group(2)),
            "max_scaled": float(final_div.group(3)),
        }
    fallback = re.search(r"HLLD-to-HLLE fallbacks = (\d+); CFL candidate range = \[([^,]+), ([^\]]+)\]", text)
    if fallback:
        result["hlld_to_hlle_fallbacks"] = int(fallback.group(1))
        result["cfl_candidate_range"] = [float(fallback.group(2)), float(fallback.group(3))]
    floor = re.search(
        r"Pressure floor events = (\d+); energy added by floor = ([^\s]+)", text)
    if floor:
        # Ненулевое число -- часть результата, а не скрытая правка: пол
        # означает, что схема в этих ячейках потеряла знак внутренней энергии.
        result["pressure_floor_events"] = int(floor.group(1))
        result["pressure_floor_energy_added"] = float(floor.group(2))
    correction = re.search(
        r"CT reconstruction magnetic-energy change: signed = ([^;]+); L1 = ([^\s]+)", text
    )
    if correction:
        result["ct_reconstruction_magnetic_energy_change"] = {
            "signed": float(correction.group(1)),
            "l1": float(correction.group(2)),
        }
    balances: dict[str, dict[str, float]] = {}
    for component, delta, boundary, residual in re.findall(
        r"Conservation component (\d+): delta=([^;]+); boundary_flux_integral=([^;]+); residual=([^\s]+)",
        text,
    ):
        balances[component] = {
            "delta": float(delta),
            "boundary_flux_integral": float(boundary),
            "residual": float(residual),
        }
    if balances:
        result["conservation"] = balances
    return result


def file_hashes(directory: Path, names: list[str]) -> dict[str, str]:
    return {name: sha256(directory / name) for name in names if (directory / name).is_file()}


def json_safe(value: Any) -> Any:
    """Replace non-finite measurements by JSON null without hiding failures."""
    if isinstance(value, float):
        return value if math.isfinite(value) else None
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, list):
        return [json_safe(item) for item in value]
    return value


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True,
                        help="clean official MHD2D source tree")
    parser.add_argument("--case", choices=sorted(CASES), required=True)
    parser.add_argument("--artifact-dir", type=Path, required=True)
    parser.add_argument("--overlay", type=Path,
                        default=ROOT / "legacy/patches/0001-legacy-corrected-physics.patch")
    parser.add_argument("--netgen-python", type=Path,
                        default=Path("/private/tmp/mhd-netgen-venv/bin/python3"))
    parser.add_argument("--mesh-backend", choices=("netgen", "structured"),
                        default="netgen",
                        help="netgen: irregular triangulation (robustness stress); "
                             "structured: dependency-free CCW right-triangle mesh "
                             "(reproducible, no Netgen, shares its hash with legacy_vkr)")
    parser.add_argument("--structured-nx", type=int,
                        help="override the case structured x resolution (rectangles)")
    parser.add_argument("--structured-ny", type=int,
                        help="override the case structured y resolution (rectangles)")
    parser.add_argument("--compiler", default=shutil.which("g++-15") or
                        "/opt/homebrew/opt/gcc/bin/g++-15")
    parser.add_argument("--jobs", type=int, default=4)
    parser.add_argument("--omp-threads", type=int, default=1)
    parser.add_argument("--maxh", type=float,
                        help="override the case mesh maximum edge size")
    args = parser.parse_args()

    source = args.source.resolve()
    overlay = args.overlay.resolve()
    artifact = args.artifact_dir.resolve()
    if artifact.exists():
        raise SystemExit(f"refusing to overwrite existing artifact directory: {artifact}")
    if not overlay.is_file():
        raise SystemExit(f"overlay patch is missing: {overlay}")
    if args.mesh_backend == "netgen" and not args.netgen_python.is_file():
        raise SystemExit(f"Netgen Python interpreter is missing: {args.netgen_python}")
    # Accept both an absolute path and a bare command name resolved on PATH;
    # the resolved path is what gets recorded in the manifest.
    resolved_compiler = args.compiler if Path(args.compiler).is_file() else shutil.which(args.compiler)
    if not resolved_compiler or not Path(resolved_compiler).is_file():
        raise SystemExit(f"C++ compiler is missing: {args.compiler}")
    args.compiler = resolved_compiler
    head = command_text(["git", "-C", str(source), "rev-parse", "HEAD"])
    if head != OFFICIAL_COMMIT:
        raise SystemExit(f"source HEAD {head} is not official {OFFICIAL_COMMIT}")
    if command_text(["git", "-C", str(source), "status", "--porcelain"]):
        raise SystemExit("source worktree is dirty; refusing to use it")

    case = dict(CASES[args.case])
    if args.maxh is not None:
        if not args.maxh > 0.0:
            raise SystemExit("--maxh must be positive")
        case["maxh"] = args.maxh
    structured_nx, structured_ny = case.get("structured_resolution", (0, 0))
    if args.structured_nx is not None:
        structured_nx = args.structured_nx
    if args.structured_ny is not None:
        structured_ny = args.structured_ny
    if args.mesh_backend == "structured" and not (structured_nx > 0 and structured_ny > 0):
        raise SystemExit("structured backend needs a positive nx/ny (case default or override)")
    artifact.mkdir(parents=True)
    worktree = artifact / "source"
    build = artifact / "build"
    errors: list[str] = []
    created_files: list[str] = []

    clone = invoke(["git", "clone", "--quiet", "--no-hardlinks", str(source), str(worktree)])
    if clone.returncode:
        raise SystemExit(f"failed to clone source: {clone.stderr}")
    checkout = invoke(["git", "checkout", "--quiet", "--detach", OFFICIAL_COMMIT], cwd=worktree)
    if checkout.returncode:
        raise SystemExit(f"failed to checkout source: {checkout.stderr}")
    apply_check = invoke(["git", "apply", "--check", "--whitespace=error", str(overlay)], cwd=worktree)
    if apply_check.returncode:
        errors.append("overlay_apply_check_failed")
    else:
        apply = invoke(["git", "apply", "--whitespace=error", str(overlay)], cwd=worktree)
        if apply.returncode:
            errors.append("overlay_apply_failed")
    if not errors and invoke(["git", "diff", "--check"], cwd=worktree).returncode:
        errors.append("overlay_whitespace_check_failed")

    input_dir, output_dir = worktree / "InputData", worktree / "OutputData"
    input_dir.mkdir(exist_ok=True)
    output_dir.mkdir(exist_ok=True)
    mesh = input_dir / "mesh.txt"
    mesh_metadata = input_dir / "mesh_metadata.json"
    xlo, xhi, ylo, yhi = case["domain"]
    if args.mesh_backend == "structured":
        # Dependency-free: emit the same minimal text format Netgen would,
        # but with a deterministic CCW right-triangle pair per rectangle.
        # Reproducible on any host and directly comparable to the AMReX
        # Cartesian solver on a matched effective dx.
        sys.path.insert(0, str(ROOT / "scripts"))
        from legacy_vkr_mesh import write_rectangular_tri_mesh  # noqa: E402
        try:
            write_rectangular_tri_mesh(mesh, xlo, xhi, ylo, yhi,
                                       structured_nx, structured_ny)
            structured_meta = {
                "generator": "scripts/legacy_vkr_mesh.py",
                "backend": "structured",
                "nx": structured_nx, "ny": structured_ny,
                "nodes": (structured_nx + 1) * (structured_ny + 1),
                "triangles": 2 * structured_nx * structured_ny,
                "dx": (xhi - xlo) / structured_nx,
                "dy": (yhi - ylo) / structured_ny,
            }
            mesh_metadata.write_text(
                json.dumps(structured_meta, indent=2, sort_keys=True) + "\n",
                encoding="utf-8")
            (artifact / "mesh_generation.log").write_text(
                json.dumps(structured_meta, sort_keys=True) + "\n", encoding="utf-8")
            generated_mesh = subprocess.CompletedProcess([], 0, "", "")
        except Exception as exc:  # noqa: BLE001 - record the failure as evidence
            (artifact / "mesh_generation.log").write_text(
                f"structured mesh generation failed: {exc}\n", encoding="utf-8")
            generated_mesh = subprocess.CompletedProcess([], 1, "", str(exc))
    else:
        mesh_command = [
            str(args.netgen_python), str(ROOT / "scripts/legacy_netgen_mesh.py"),
            "--output", str(mesh), "--metadata", str(mesh_metadata),
            "--xlo", str(xlo), "--xhi", str(xhi),
            "--ylo", str(ylo), "--yhi", str(yhi),
            "--maxh", str(case["maxh"]), "--grading", str(case["grading"]),
            "--boundary-jitter", str(case["boundary_jitter"]),
        ]
        generated_mesh = invoke(mesh_command, cwd=worktree)
        (artifact / "mesh_generation.log").write_text(
            (generated_mesh.stdout or "") + (generated_mesh.stderr or ""), encoding="utf-8"
        )
    created_files.append("mesh_generation.log")
    if generated_mesh.returncode:
        errors.append("mesh_generation_failed")

    config = {
        "taskType": case["task_type"],
        "finalTime": case["final_time"],
        "cfl": case["cfl"],
        "debugDivergence": False,
        "ghostOutput": False,
        "iterationsPerFrame": 1_000_000_000,
        "fieldLoopU": case.get("field_loop_u", 2.0),
        "fieldLoopV": case.get("field_loop_v", 1.0),
        "fieldLoopRadius": case.get("field_loop_radius", 0.3),
        "fieldLoopAmplitude": case.get("field_loop_amplitude", 1.0e-3),
        "cylindrical": False,
        "gpu": False,
        # D-009: реконструкция RT0 меняет представление поля, а не греет газ,
        # поэтому по умолчанию сохраняется внутренняя энергия. Пол по давлению
        # держим включённым как страховку -- его срабатывания попадают в лог и
        # в манифест, то есть остаются частью результата.
        "ctEnergyMode": case.get("ct_energy_mode", "preserve_internal"),
        "pressureFloor": case.get("pressure_floor", 1.0e-10),
        "exportFileName": "OutputData/final.vtu",
    }
    config_path = input_dir / "solverConfig.json"
    config_path.write_text(json.dumps(config, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    build_log = artifact / "build.log"
    configure = build_result = None
    if not errors:
        with build_log.open("w", encoding="utf-8") as log:
            configure = invoke(
                ["cmake", "-S", str(worktree), "-B", str(build), "-DBUILD_TESTING=ON",
                 "-DCMAKE_BUILD_TYPE=Release", f"-DCMAKE_CXX_COMPILER={args.compiler}"],
                stdout=log,
            )
            build_result = invoke(["cmake", "--build", str(build), "--parallel", str(args.jobs)],
                                  stdout=log) if not configure.returncode else configure
    else:
        build_log.write_text("Build skipped because overlay application or mesh generation failed.\n",
                             encoding="utf-8")
    created_files.append("build.log")
    if configure is not None and (configure.returncode or build_result.returncode):
        errors.append("build_failed")

    ctest_log = artifact / "ctest.log"
    if not errors:
        with ctest_log.open("w", encoding="utf-8") as log:
            ctest_result = invoke(["ctest", "--test-dir", str(build), "--output-on-failure"], stdout=log)
        if ctest_result.returncode:
            errors.append("ctest_failed")
    else:
        ctest_log.write_text("CTest skipped because build setup failed.\n", encoding="utf-8")
        ctest_result = None
    created_files.append("ctest.log")

    solver_log = artifact / "solver.log"
    solver_result: subprocess.CompletedProcess[str] | None = None
    environment = os.environ.copy()
    environment["OMP_NUM_THREADS"] = str(args.omp_threads)
    if not errors:
        with solver_log.open("w", encoding="utf-8") as log:
            solver_result = invoke([str(build / "MHD2D")], cwd=worktree,
                                   environment=environment, stdout=log)
        if solver_result.returncode:
            errors.append(f"solver_exit_{solver_result.returncode}")
    else:
        solver_log.write_text("Solver skipped because build or CTest failed.\n", encoding="utf-8")
    created_files.append("solver.log")

    # Preserve inputs/outputs at artifact root as stable raw evidence. The
    # clone remains for forensic inspection but paths in the manifest point to
    # these compact copies.
    for source_file, target_name in (
        (mesh, "mesh.txt"),
        (mesh_metadata, "mesh_metadata.json"),
        (config_path, "solver_config.json"),
        (output_dir / "tmpres_0.vtu", "initial_pre_update.vtu"),
        (output_dir / "final.vtu", "final.vtu"),
        (output_dir / "physical_failure.json", "physical_failure.json"),
    ):
        if copy_if_present(source_file, artifact / target_name):
            created_files.append(target_name)
    copy_if_present(overlay, artifact / "applied_overlay.patch")
    created_files.append("applied_overlay.patch")

    final_summary_path = artifact / "final_summary.json"
    initial_summary_path = artifact / "initial_summary.json"
    final_diagnostics = initial_diagnostics = None
    if (artifact / "final.vtu").is_file():
        analysis_command = [
            sys.executable, str(ROOT / "scripts/analyze_legacy_vtu.py"),
            "--mesh", str(artifact / "mesh.txt"),
            "--vtu", str(artifact / "final.vtu"),
            "--case", str(case["diagnostic_case"]),
            "--gamma", str(case["gamma"]),
            "--output", str(final_summary_path),
        ]
        if case["diagnostic_case"] == "magnetic_loop" and \
                (artifact / "initial_pre_update.vtu").is_file():
            analysis_command.extend([
                "--reference-vtu", str(artifact / "initial_pre_update.vtu"),
            ])
        analysis = invoke(analysis_command)
        if analysis.returncode:
            errors.append("final_analysis_failed")
        else:
            final_diagnostics = read_json_if_present(final_summary_path)
            created_files.append("final_summary.json")
    else:
        errors.append("final_vtu_missing")
    if (artifact / "initial_pre_update.vtu").is_file():
        initial_command = [
            sys.executable, str(ROOT / "scripts/analyze_legacy_vtu.py"),
            "--mesh", str(artifact / "mesh.txt"),
            "--vtu", str(artifact / "initial_pre_update.vtu"),
            "--case", str(case["diagnostic_case"]),
            "--gamma", str(case["gamma"]),
            "--output", str(initial_summary_path),
        ]
        initial_analysis = invoke(initial_command)
        if initial_analysis.returncode:
            errors.append("initial_analysis_failed")
        else:
            initial_diagnostics = read_json_if_present(initial_summary_path)
            created_files.append("initial_summary.json")

    profile_summary: dict[str, Any] | None = None
    if case["profile_projection"] and (artifact / "final.vtu").is_file():
        profile_csv, profile_json = artifact / "brio_profile.csv", artifact / "brio_profile_summary.json"
        projection = invoke([
            sys.executable, str(ROOT / "scripts/project_legacy_vtu.py"),
            "--mesh", str(artifact / "mesh.txt"), "--vtu", str(artifact / "final.vtu"),
            "--gamma", str(case["gamma"]), "--bins", "256",
            "--csv", str(profile_csv), "--summary", str(profile_json),
        ])
        if projection.returncode:
            errors.append("brio_projection_failed")
        else:
            profile_summary = read_json_if_present(profile_json)
            created_files.extend(["brio_profile.csv", "brio_profile_summary.json"])

    archive = artifact / "official_source.tar"
    with archive.open("wb") as stream:
        archive_result = subprocess.run(
            ["git", "-C", str(source), "archive", "--format=tar", OFFICIAL_COMMIT],
            stdout=stream, stderr=subprocess.PIPE, check=False,
        )
    if archive_result.returncode:
        errors.append("source_archive_failed")
    else:
        created_files.append("official_source.tar")

    log_text = solver_log.read_text(encoding="utf-8")
    parsed_log = parse_solver_log(log_text)
    reported_time = parsed_log.get("reported_final_time")
    reached_target = (reported_time is not None and
                      # The legacy executable prints its final time with the
                      # stream's default precision.  This is a parsing
                      # tolerance only; the solver itself clips its last dt
                      # to finalTime.
                      abs(reported_time - case["final_time"]) <=
                      1.0e-8 * max(1.0, abs(case["final_time"])))
    if solver_result is not None and solver_result.returncode == 0 and not reached_target:
        errors.append("target_time_not_reached")
    finite_and_positive = bool(final_diagnostics and final_diagnostics.get("finite") and
                               final_diagnostics.get("rho_min", 0.0) > 0.0 and
                               final_diagnostics.get("pressure_min", 0.0) > 0.0)
    if final_diagnostics is not None and not finite_and_positive:
        errors.append("final_state_nonphysical")

    periodic_energy_balance_ok: bool | None = None
    if case["task_type"] in (8, 9) and solver_result is not None:
        energy_residual = (parsed_log.get("conservation", {}).get("4", {}).get("residual"))
        initial_total_energy = (initial_diagnostics or {}).get("total_energy")
        if isinstance(energy_residual, (int, float)) and \
                isinstance(initial_total_energy, (int, float)):
            periodic_energy_balance_ok = (
                abs(energy_residual) / max(abs(initial_total_energy), 1.0) <= 1.0e-10
            )
        if periodic_energy_balance_ok is not True:
            errors.append("periodic_total_energy_balance_failed")

    energy_ratio = None
    if initial_diagnostics and final_diagnostics:
        initial_energy = initial_diagnostics.get("magnetic_energy")
        final_energy = final_diagnostics.get("magnetic_energy")
        if isinstance(initial_energy, (int, float)) and initial_energy > 0.0 and \
                isinstance(final_energy, (int, float)):
            energy_ratio = final_energy / initial_energy

    manifest = {
        "schema_version": 3,
        "profile": "legacy_corrected",
        "source": {
            "official_commit": OFFICIAL_COMMIT,
            "official_source_sha256": sha256(archive) if archive.is_file() else None,
            "source_dirty": False,
            "overlay_file": "applied_overlay.patch",
            "overlay_sha256": sha256(artifact / "applied_overlay.patch"),
            "applied_diff_sha256": hashlib.sha256(
                command_text(["git", "-C", str(worktree), "diff", "--binary"]).encode()
            ).hexdigest(),
        },
        "case": {
            "id": args.case,
            "task_type": case["task_type"],
            "gamma": case["gamma"],
            "final_time": case["final_time"],
            "cfl": case["cfl"],
            "domain": {"x": [xlo, xhi], "y": [ylo, yhi]},
            "mesh": {
                "backend": args.mesh_backend,
                "generator": ("scripts/legacy_vkr_mesh.py"
                              if args.mesh_backend == "structured"
                              else "scripts/legacy_netgen_mesh.py"),
                "maxh": None if args.mesh_backend == "structured" else case["maxh"],
                "grading": None if args.mesh_backend == "structured" else case["grading"],
                "boundary_jitter": (None if args.mesh_backend == "structured"
                                    else case["boundary_jitter"]),
                "structured_resolution": ([structured_nx, structured_ny]
                                          if args.mesh_backend == "structured" else None),
                "mesh_sha256": sha256(artifact / "mesh.txt") if (artifact / "mesh.txt").is_file() else None,
                "metadata": read_json_if_present(artifact / "mesh_metadata.json"),
            },
            "field_loop": {
                "u": config["fieldLoopU"], "v": config["fieldLoopV"],
                "radius": config["fieldLoopRadius"], "amplitude": config["fieldLoopAmplitude"],
            } if case["task_type"] == 9 else None,
        },
        "build": {
            "compiler": str(Path(args.compiler).resolve()),
            "compiler_version": command_text([args.compiler, "--version"]),
            "cmake_version": command_text(["cmake", "--version"]),
            "build_type": "Release",
            "build_exit_code": None if build_result is None else build_result.returncode,
            "ctest_exit_code": None if ctest_result is None else ctest_result.returncode,
        },
        "environment": {
            "platform": platform.platform(),
            "machine": platform.machine(),
            "processor": platform.processor() or "unavailable",
            "python": sys.version,
            "netgen_python": str(args.netgen_python),
            "omp_num_threads": args.omp_threads,
        },
        "solver": {
            "process_exit_code": None if solver_result is None else solver_result.returncode,
            "reported": parsed_log,
            "reached_configured_final_time": reached_target,
            "initial_snapshot_semantics": (
                "initial_pre_update.vtu contains the pre-update state written at iteration 0; "
                "the legacy console timestamp beside it is the first candidate dt, not its state time"
            ),
            "failure_dump": read_json_if_present(artifact / "physical_failure.json"),
        },
        "diagnostics": {
            "analyzer": {
                "path": "scripts/analyze_legacy_vtu.py",
                "sha256": sha256(ROOT / "scripts/analyze_legacy_vtu.py"),
            },
            "initial": initial_diagnostics,
            "final": final_diagnostics,
            "magnetic_energy_ratio_final_over_initial": energy_ratio,
            "brio_area_weighted_profile": profile_summary,
        },
        "quality_gate": {
            "ctest_passed": ctest_result is not None and ctest_result.returncode == 0,
            "solver_exit_zero": solver_result is not None and solver_result.returncode == 0,
            "reached_target_time": reached_target,
            "finite_positive_final_state": finite_and_positive,
            "periodic_total_energy_balance_passed": periodic_energy_balance_ok,
            "status": "pass" if not errors else "fail",
            "failure_reasons": errors,
        },
        "files": file_hashes(artifact, sorted(set(created_files))),
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
    }
    manifest = json_safe(manifest)
    manifest_path = artifact / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n",
                             encoding="utf-8")
    print(json.dumps({
        "artifact_dir": str(artifact),
        "quality_status": manifest["quality_gate"]["status"],
        "failure_reasons": errors,
        "final": final_diagnostics,
        "reported": parsed_log,
    }, indent=2, sort_keys=True, allow_nan=False))
    return 1 if errors else 0


if __name__ == "__main__":
    raise SystemExit(main())
