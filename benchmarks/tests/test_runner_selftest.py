"""H. Self-test for the WP-03 production timing harness.

Exercises the real UppASD CPU executable against the INFRA_test_only case:
different Nstep values are generated correctly, slope/intercept extraction
works from real production runs, a genuinely failed run is excluded from the
fit, and the source template is left untouched throughout.

Skipped (not failed) wherever no runnable UppASD executable can be found in
this environment -- see blueprint section 16: authoritative performance
testing needs controlled hardware, but even a harness self-test must not
hard-fail on a machine that simply has no build.
"""

from __future__ import annotations

import os
import pathlib
import shutil
import subprocess

import pytest

from harness import cases
from harness import runner
from harness import steady_state

BENCHMARKS_DIR = pathlib.Path(__file__).resolve().parent.parent
REPO_ROOT = BENCHMARKS_DIR.parent
INFRA_CASE_PATH = BENCHMARKS_DIR / "cases" / "INFRA_test_only" / "case.yaml"

_CANDIDATE_BINARIES = ["bin/sd.gfortran", "bin/sd.f95"]
_CANDIDATE_LIB_DIRS = [
    pathlib.Path.home() / ".local" / "lib",
    pathlib.Path.home() / "miniconda3" / "lib",
    pathlib.Path.home() / "anaconda3" / "lib",
]


def _resolve_env(binary_path):
    """Best-effort runtime environment so a real MKL-linked binary can start.

    Purely a self-test/dev convenience for locating a working runtime on
    *this* machine; production campaigns are expected to supply a correct
    environment themselves (see ``runner.run_sample``'s ``env`` parameter).
    Forces single-threaded execution so the self-test's timing signal is not
    at the mercy of whatever OMP_NUM_THREADS happens to be set in the shell.
    """
    env = os.environ.copy()
    ldd = shutil.which("ldd")
    if ldd is not None:
        probe = subprocess.run([ldd, str(binary_path)], text=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        missing = [line.split()[0] for line in probe.stdout.splitlines() if "not found" in line]
        if missing:
            extra_dirs = [
                str(candidate) for candidate in _CANDIDATE_LIB_DIRS
                if candidate.is_dir() and any((candidate / name).exists() for name in missing)
            ]
            if not extra_dirs:
                return None
            env["LD_LIBRARY_PATH"] = os.pathsep.join(extra_dirs + [env.get("LD_LIBRARY_PATH", "")])
            check = subprocess.run(
                [ldd, str(binary_path)], text=True, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, env=env
            )
            if "not found" in check.stdout:
                return None
    env["OMP_NUM_THREADS"] = "1"
    return env


def _find_binary():
    for relative in _CANDIDATE_BINARIES:
        candidate = REPO_ROOT / relative
        if not candidate.is_file():
            continue
        env = _resolve_env(candidate)
        if env is not None:
            return candidate, env
    return None, None


@pytest.fixture(scope="module")
def binary_and_env():
    binary_path, env = _find_binary()
    if binary_path is None:
        pytest.skip("no runnable UppASD CPU executable found in this environment")
    return binary_path, env


@pytest.fixture
def infra_case():
    return cases.load_case_manifest(INFRA_CASE_PATH)


def _context(binary_path, env):
    return runner.developer_context(binary_path, omp_threads=int(env["OMP_NUM_THREADS"]))


def test_different_nstep_values_are_generated_correctly(tmp_path, infra_case, binary_and_env):
    _binary_path, _env = binary_and_env
    work_root = tmp_path / "work"
    for nstep in (25, 40):
        run_dir = cases.generate_run_directory(
            infra_case, "T0", "2x2x2", work_root, run_id=f"nstep-check-{nstep}", extra_overrides={"Nstep": nstep}
        )
        assert cases.read_keyword(run_dir.path, "Nstep") == str(nstep)


def test_slope_and_intercept_extraction_from_real_production_runs(tmp_path, infra_case, binary_and_env):
    binary_path, env = binary_and_env
    work_root = tmp_path / "work"
    results_dir = tmp_path / "results"
    context = _context(binary_path, env)

    workload_metadata = runner.resolve_workload_metadata(infra_case, "T0", "2x2x2", work_root, binary_path, env=env)
    assert workload_metadata["natom"] == 16

    nsteps = (50_000, 150_000, 300_000)
    samples = []
    for index, nstep in enumerate(nsteps):
        result = runner.run_sample(
            infra_case, "T0", "2x2x2",
            nstep=nstep, run_id=f"fit-{nstep}", campaign_id="selftest",
            sample_index=index, work_root=work_root, results_dir=results_dir,
            binary_path=binary_path, measurement_profile="DYNAMICS_ONLY",
            workload_metadata=workload_metadata, context=context, env=env,
        )
        assert result.record["run_status"] == "COMPLETED"
        assert result.record["numerical_valid"] is True
        assert result.record["nstep"] == nstep
        samples.append((nstep, result.record["process_wall_seconds"]))

    fit = steady_state.fit_multi_nstep(samples)
    assert fit.fit_method == steady_state.MULTI_NSTEP_LEAST_SQUARES
    assert fit.nstep_fit_points == sorted(nsteps)
    assert fit.steady_step_seconds > 0.0

    two_point = steady_state.two_point_estimate(nsteps[0], samples[0][1], nsteps[-1], samples[-1][1])
    # Loose cross-check only -- real wall-clock noise, not an exact model.
    assert fit.steady_step_seconds == pytest.approx(two_point.steady_step_seconds, rel=0.5)


def test_a_genuinely_failed_run_is_excluded_from_the_fit(tmp_path, infra_case, binary_and_env):
    binary_path, env = binary_and_env
    work_root = tmp_path / "work"
    results_dir = tmp_path / "results"
    context = _context(binary_path, env)

    workload_metadata = runner.resolve_workload_metadata(infra_case, "T0", "2x2x2", work_root, binary_path, env=env)

    good = runner.run_sample(
        infra_case, "T0", "2x2x2",
        nstep=10_000, run_id="ok", campaign_id="selftest", sample_index=0,
        work_root=work_root, results_dir=results_dir, binary_path=binary_path,
        measurement_profile="DYNAMICS_ONLY", workload_metadata=workload_metadata,
        context=context, env=env,
    )
    # A genuinely broken production invocation -- a real nonzero-exit process,
    # not a fabricated failure -- must be marked FAILED and excluded below.
    failed = runner.run_sample(
        infra_case, "T0", "2x2x2",
        nstep=20_000, run_id="broken", campaign_id="selftest", sample_index=1,
        work_root=work_root, results_dir=results_dir, binary_path="/bin/false",
        measurement_profile="DYNAMICS_ONLY", workload_metadata=workload_metadata,
        context=context, env=env,
    )

    assert good.record["run_status"] == "COMPLETED"
    assert failed.record["run_status"] == "FAILED"
    assert failed.record["numerical_valid"] is False
    assert failed.record["environment_valid"] is False
    assert failed.record["process_wall_seconds"] is None
    # Workload fields describe the configured cell, not the outcome -- still real.
    assert failed.record["natom"] == workload_metadata["natom"]

    fit_samples = [
        (record["nstep"], record["process_wall_seconds"])
        for record in (good.record, failed.record)
        if record["run_status"] == "COMPLETED"
    ]
    assert fit_samples == [(10_000, good.record["process_wall_seconds"])]
    with pytest.raises(steady_state.SteadyStateError):
        # Only one distinct nstep survives exclusion -- correctly unfittable.
        steady_state.fit_multi_nstep(fit_samples)


def test_source_template_remains_untouched(tmp_path, infra_case, binary_and_env):
    binary_path, env = binary_and_env
    work_root = tmp_path / "work"
    results_dir = tmp_path / "results"
    context = _context(binary_path, env)
    before = cases.compute_template_fingerprint(infra_case.template_dir)

    workload_metadata = runner.resolve_workload_metadata(infra_case, "T0", "2x2x2", work_root, binary_path, env=env)
    runner.run_sample(
        infra_case, "T0", "2x2x2",
        nstep=5_000, run_id="template-integrity", campaign_id="selftest", sample_index=0,
        work_root=work_root, results_dir=results_dir, binary_path=binary_path,
        measurement_profile="PRODUCTION_LIGHT", workload_metadata=workload_metadata,
        context=context, env=env,
    )

    after = cases.compute_template_fingerprint(infra_case.template_dir)
    assert before == after
