"""Run the LSWT-to-trajectory OAM bridge check (Blueprint B5.4).

The harness reports two related but distinct checks:

* the precession frequency of a narrow honeycomb magnon packet is compared
  with the energy from ``do_oam_lswt``;
* the settled particle-field bridge is checked with the independent
  reciprocal-supercell packet, whose analytic gradient contains the angular
  variation of the band spinor.

The production trajectory packet is a fixed sublattice-spinor control. Its
targets are ``lambda(l=0) ~= 0``, ``lambda(l=1) ~= 1``, and a shift of one;
it contains no angular variation of the band spinor and therefore no Berry
term. Its per-sublattice value is not used as a quantitative ``F_n`` check at
the tested short wavelength: the intrinsic connection is a sum of
per-sublattice terms, but linear-FEM gradients become biased as ``k.a`` grows.
B5.4 validates the independent oracle against LSWT, not the production
trajectory kernel. The accepted bridge target is ``lambda = l - 2 F_n`` for
the C1/ ``exp(+i k.r)`` particle convention.

Example::

    UPPASD=./build-oam/uppasd python3 run_bridge_checks.py \\
        --workdir /tmp/oam-bridge

Use ``--force`` only for a generated work directory that may be replaced.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np

import mkhoney
import oracle_honey
import oracle_bridge
import oracle_traj


S3 = np.sqrt(3.0)
N = 24
CLEAN_N = 48
C1 = np.array([1.0, 0.0, 0.0])
C2 = np.array([0.5, S3 / 2.0, 0.0])
BASIS = np.array([[0.0, 0.0, 0.0], [0.5, 0.5 / S3, 0.0]])
J = 1.0
D = 0.30
AMP = 0.025
SIGMA = 4.0
TIMESTEP_S = 1.0e-16
HBAR_MEV_PS = 0.6582119569

# U5a/B5.5 angularly varying honeycomb packet.  Keep this implementation in
# the harness: the archived reference under docs/dev/oam/prompts/rounds/r4 is
# an audit record, not a
# runtime dependency.
BAND_N = 90
BAND_DK = 0.15
BAND_NPHI = 192
BAND_NR = 13
BAND_D = 1.0
BAND = 0


def coords_honeycomb(n: int = N) -> np.ndarray:
    """Return the atom order used by UppASD for the generated posfile."""

    return oracle_traj.build_coords(n, n, C1, C2, BASIS)


def write_restart(directory: Path, moments: np.ndarray) -> None:
    natom = moments.shape[1]
    with (directory / "restart.in").open("w") as handle:
        handle.write("#" * 80 + "\n")
        handle.write("# File type: R\n# Simulation type: S\n")
        handle.write(f"# Number of atoms: {natom:9d}\n")
        handle.write("# Number of ensembles:         1\n")
        handle.write("#" * 80 + "\n")
        handle.write("  # iter     ens   iatom           |Mom|             M_x             M_y             M_z\n")
        for atom in range(natom):
            mx, my, mz = moments[:, atom]
            handle.write(
                f"{0:8d}{1:8d}{atom + 1:8d}  {1.0:16.8E}"
                f"{mx:24.16E}{my:24.16E}{mz:24.16E}\n"
            )


def packet(coords: np.ndarray, mode: np.ndarray, k0: float, ell: int,
           *, clean: bool = False, n: int = N) -> np.ndarray:
    """Construct a small-amplitude optical-band packet in the global frame."""

    xy = coords[:2]
    if clean:
        center = 0.5 * (n * C1[:2] + n * C2[:2])
        cell = np.column_stack((n * C1[:2], n * C2[:2]))
        dxy = xy - center[:, None]
        reduced = np.linalg.solve(cell, dxy)
        reduced -= np.round(reduced)
        dxy = cell @ reduced
    else:
        center = xy.mean(axis=1)
        dxy = xy - center[:, None]
    radius = np.linalg.norm(dxy, axis=0)
    theta = np.arctan2(dxy[1], dxy[0])
    envelope = np.exp(-0.5 * (radius / SIGMA) ** 2)
    if clean:
        envelope *= (radius / SIGMA) ** abs(ell)
    phase = np.exp(1j * (k0 * xy[0] + ell * theta))
    sublattice_mode = np.asarray([mode[0], mode[1]], complex)
    psi = AMP * sublattice_mode[np.arange(coords.shape[1]) % 2] * envelope * phase
    mz = np.sqrt(np.clip(1.0 - np.abs(psi) ** 2, 0.0, None))
    return np.vstack((psi.real, psi.imag, mz))


def write_trajectory_case(directory: Path, *, simid: str, k0: float,
                           mode: np.ndarray, ell: int, nstep: int,
                           gradient: str = "fem", n: int = N,
                           clean: bool = False, damping: float = 0.0,
                           temperature: float = 0.0) -> np.ndarray:
    """Write one trajectory input and return its prescribed complex packet."""

    directory.mkdir(parents=True, exist_ok=True)
    # This supplies the common honeycomb geometry, exchange and DMI files.
    mkhoney.write(str(directory), C2, J=J, D=D, kgrid=(30, 30), nphi=64, nr=16)
    coords = coords_honeycomb(n)
    moments = packet(coords, mode, k0, ell, clean=clean, n=n)
    write_restart(directory, moments)
    origin = coords[:2].mean(axis=1)
    with (directory / "inpsd.dat").open("w") as handle:
        handle.write(f"""simid {simid}
ncell {n} {n} 1
BC P P 0
cell {C1[0]:.10f} {C1[1]:.10f} 0.0
     {C2[0]:.10f} {C2[1]:.10f} 0.0
     0.0 0.0 1.0
Sym 0
posfile ./posfile
posfiletype C
momfile ./momfile
exchange ./jfile
dm ./dmfile
maptype 1
initmag 4
restartfile ./restart.in
ip_mode N
mode S
temp {temperature}
damping {damping}
Nstep {nstep}
timestep {TIMESTEP_S:.16e}
do_avrg N
do_prnstruct 0
do_tottraj Y
tottraj_step 1
tottraj_buff 32
do_oam_traj Y
oam_gradient {gradient}
oam_step 10
oam_buff 10
oam_origin {origin[0]:.16e} {origin[1]:.16e} 0.0
oam_weight site
oam_sigma_max 0.6
{"oam_axis 0.0 0.0 1.0" if clean else ""}
do_chern N
do_oam_lswt N
""")
    return moments[0] + 1j * moments[1]


def run_binary(binary: str, directory: Path) -> None:
    env = os.environ.copy()
    env.setdefault("OMP_NUM_THREADS", "1")
    executable = shutil.which(binary) if os.path.sep not in binary else str(Path(binary).expanduser().resolve())
    if executable is None:
        raise FileNotFoundError(f"UppASD executable not found: {binary}")
    result = subprocess.run([executable], cwd=directory, text=True,
                            capture_output=True, env=env)
    (directory / "uppasd.log").write_text(result.stdout + result.stderr)
    if result.returncode:
        tail = (result.stdout + result.stderr).splitlines()[-30:]
        raise RuntimeError(f"UppASD failed in {directory}:\n" + "\n".join(tail))


def load_lswt_reference(directory: Path, target_k: float = 2.9):
    output = np.loadtxt(directory / "oam_lswt.honey.out", comments="#")
    rows = output[np.isclose(output[:, 1], 1.0)]
    row = rows[np.argmin(np.abs(rows[:, 0] - target_k))]
    k0, energy, fishman = float(row[0]), float(row[2]), float(row[3])
    _, vectors = oracle_honey.evecs(np.array([k0, 0.0]), J=J, D=D, S=1.0)
    return k0, energy, fishman, vectors[:, 0]


def load_oam(directory: Path, simid: str) -> np.ndarray:
    rows = []
    for line in (directory / f"oam_traj.{simid}.out").read_text().splitlines():
        if line.strip() and not line.lstrip().startswith("#"):
            rows.append([float(value) for value in line.split()])
    if not rows:
        raise RuntimeError(f"No OAM rows in {directory}/oam_traj.{simid}.out")
    return np.asarray(rows)


def load_moments(directory: Path, simid: str, natom: int):
    """Read standard moment output as (steps, mx+imy arrays)."""

    grouped = {}
    filename = directory / f"moment.{simid}.out"
    for line in filename.read_text().splitlines():
        if not line.strip() or line.lstrip().startswith("#"):
            continue
        values = line.split()
        if len(values) < 7:
            continue
        step = int(float(values[0]))
        atom = int(float(values[2])) - 1
        mx, my = float(values[4]), float(values[5])
        grouped.setdefault(step, np.zeros((2, natom), float))[:, atom] = (mx, my)
    if len(grouped) < 8:
        raise RuntimeError(f"Too few moment snapshots in {filename}: {len(grouped)}")
    steps = np.array(sorted(grouped), dtype=float)
    transverse = np.stack([grouped[int(step)][0] + 1j * grouped[int(step)][1]
                            for step in steps])
    return steps, transverse


def frequency_from_projection(directory: Path, simid: str, psi0: np.ndarray):
    steps, transverse = load_moments(directory, simid, psi0.size)
    projection = transverse @ np.conjugate(psi0)
    phase = np.unwrap(np.angle(projection))
    first = max(1, len(steps) // 10)
    slope, intercept = np.polyfit(steps[first:], phase[first:], 1)
    frequency = abs(slope) / (TIMESTEP_S * 1.0e12)
    quality = np.abs(projection[first:])
    return frequency, float(np.std(phase[first:] - (slope * steps[first:] + intercept))), steps, quality


def result(label: str, passed: bool, detail: str) -> bool:
    print(f"[{'PASS' if passed else 'FAIL'}] {label}: {detail}")
    return passed


def band_ring_gauge(k: float):
    phis = 2.0 * np.pi * np.arange(BAND_NPHI) / BAND_NPHI
    us = [oracle_honey.evecs(np.array([k * np.cos(phi), k * np.sin(phi)]),
                             D=BAND_D)[1][:, BAND] for phi in phis]
    for index in range(1, BAND_NPHI):
        overlap = np.vdot(us[index - 1], us[index])
        us[index] = us[index] * np.conj(overlap) / abs(overlap)
    theta = np.angle(np.vdot(us[-1], us[0]))
    return phis, [u * np.exp(1j * index * theta / BAND_NPHI)
                  for index, u in enumerate(us)]


def band_packet_with_gradient(xy: np.ndarray, center: np.ndarray, k0: float,
                              ell: int):
    krs = np.linspace(k0 - 2.5 * BAND_DK, k0 + 2.5 * BAND_DK, BAND_NR)
    positions = xy[:2] - center[:, None]
    sublattice = np.arange(xy.shape[1]) % 2
    field = np.zeros(xy.shape[1], complex)
    grad_x = np.zeros(xy.shape[1], complex)
    grad_y = np.zeros(xy.shape[1], complex)
    previous = None
    for kr in krs:
        phis, spinors = band_ring_gauge(kr)
        if previous is not None:
            overlap = np.vdot(previous, spinors[0])
            spinors = [u * np.conj(overlap) / abs(overlap) for u in spinors]
        previous = spinors[0]
        envelope = np.exp(-(kr - k0) ** 2 / (2.0 * BAND_DK ** 2))
        for phi, spinor in zip(phis, spinors):
            wavevector = kr * np.array([np.cos(phi), np.sin(phi)])
            term = (envelope * np.exp(1j * ell * phi) *
                    spinor[sublattice] * np.exp(1j * (wavevector @ positions)))
            field += term
            grad_x += 1j * wavevector[0] * term
            grad_y += 1j * wavevector[1] * term
    scale = 0.05 / np.abs(field).max()
    return field * scale, grad_x * scale, grad_y * scale


def band_exact_lambda(xy: np.ndarray, cell: np.ndarray, inv_cell: np.ndarray,
                      field: np.ndarray, grad_x: np.ndarray,
                      grad_y: np.ndarray) -> float:
    weight = np.abs(field) ** 2
    theta = 2.0 * np.pi * ((inv_cell @ xy[:2]) % 1.0)
    reduced_center = np.mod(np.arctan2(
        (weight * np.sin(theta)).sum(axis=1),
        (weight * np.cos(theta)).sum(axis=1)) / (2.0 * np.pi), 1.0)
    center = cell @ reduced_center
    reduced = inv_cell @ (xy[:2] - center[:, None])
    reduced -= np.round(reduced)
    lever = cell @ reduced
    return float(np.sum(np.imag(np.conj(field) *
                                (lever[0] * grad_y - lever[1] * grad_x))) /
                 weight.sum())


def band_lswt_prediction(k0: float, ell: int) -> float:
    krs = np.linspace(k0 - 2.5 * BAND_DK, k0 + 2.5 * BAND_DK, BAND_NR)
    kk = np.linspace(0.0, krs[-1], 400)
    berry = oracle_honey.F_of_k(kk, BAND, D=BAND_D, nphi=384)
    berry_at_packet = np.interp(krs, kk, berry)
    weight = np.exp(-(krs - k0) ** 2 / BAND_DK ** 2) * krs
    return float(ell - 2.0 * np.sum(weight * berry_at_packet) / weight.sum())


def write_band_case(directory: Path, field: np.ndarray, simid: str,
                    gradient: str) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    tau = [np.array(BASIS[0][:2]), np.array(BASIS[1][:2])]
    (directory / "posfile").write_text("".join(
        f"{index + 1} {index + 1} {position[0]:.14f} "
        f"{position[1]:.14f} 0.0\n" for index, position in enumerate(tau)))
    (directory / "momfile").write_text(
        "1 1 1.0 0.0 0.0 1.0\n2 1 1.0 0.0 0.0 1.0\n")
    nn = [(0.5, S3 / 6.0), (-0.5, S3 / 6.0), (0.0, -S3 / 3.0)]
    (directory / "jfile").write_text(
        "".join(f"1 2 {dx:.14f} {dy:.14f} 0.0 1.0\n" for dx, dy in nn) +
        "".join(f"2 1 {-dx:.14f} {-dy:.14f} 0.0 1.0\n" for dx, dy in nn))
    mz = np.sqrt(1.0 - np.abs(field) ** 2)
    with (directory / "restart.in").open("w") as handle:
        handle.write("#" * 80 + "\n# File type: R\n# Simulation type: S\n")
        handle.write(f"# Number of atoms: {field.size:9d}\n")
        handle.write("# Number of ensembles:         1\n" + "#" * 80 + "\n")
        handle.write("  # iter     ens   iatom           |Mom|             M_x             M_y             M_z\n")
        for index, value in enumerate(field):
            handle.write(f"{0:8d}{1:8d}{index + 1:8d}  {1.0:16.8E}"
                         f"{value.real:24.16E}{value.imag:24.16E}"
                         f"{mz[index]:24.16E}\n")
    (directory / "inpsd.dat").write_text(f"""simid {simid}
ncell {BAND_N} {BAND_N} 1
BC P P 0
cell 1.0 0.0 0.0
     0.5 {S3 / 2.0:.14f} 0.0
     0.0 0.0 1.0
Sym 0
posfile ./posfile
posfiletype C
momfile ./momfile
exchange ./jfile
maptype 1
do_prnstruct 1
initmag 4
restartfile ./restart.in
ip_mode N
mode S
temp 0.0
damping 0.0
Nstep 1
timestep 1e-20
do_avrg N
do_oam_traj Y
oam_step 1
oam_gradient {gradient}
oam_axis 0 0 1
""")


def run_band_case(binary: str, directory: Path, field: np.ndarray,
                  gradient: str, simid: str) -> tuple[float | None, str]:
    write_band_case(directory, field, simid, gradient)
    executable = (shutil.which(binary) if os.path.sep not in binary
                  else str(Path(binary).expanduser().resolve()))
    if executable is None:
        raise FileNotFoundError(f"UppASD executable not found: {binary}")
    env = os.environ.copy()
    env.setdefault("OMP_NUM_THREADS", "1")
    completed = subprocess.run([executable], cwd=directory, text=True,
                               capture_output=True, env=env)
    log = completed.stdout + completed.stderr
    (directory / "uppasd.log").write_text(log)
    if completed.returncode:
        tail = "\n".join(log.splitlines()[-30:])
        raise RuntimeError(f"UppASD failed in {directory}:\n{tail}")
    output = directory / f"oam_traj.{simid}.out"
    if not output.exists():
        return None, log
    rows = [line.split() for line in output.read_text().splitlines()
            if line.strip() and not line.lstrip().startswith("#")]
    return (float(rows[0][2]) if rows else None), log


def run_band_check(binary: str, workdir: Path) -> int:
    cell = np.column_stack((BAND_N * C1[:2], BAND_N * C2[:2]))
    inv_cell = np.linalg.inv(cell)
    xy = coords_honeycomb(BAND_N)
    center = xy[:2].mean(axis=1)
    records = []
    spectral_probe = None
    for k0 in (1.0, 2.0):
        for ell in (0, 1):
            field, grad_x, grad_y = band_packet_with_gradient(xy, center, k0, ell)
            exact = band_exact_lambda(xy, cell, inv_cell, field, grad_x, grad_y)
            target = band_lswt_prediction(k0, ell)
            spectral, log = run_band_case(binary, workdir / f"s{k0:g}{ell}",
                                          field, "spectral", f"b55s{k0:g}{ell}")
            if spectral is None and spectral_probe is None:
                spectral_probe = log
            fem, _ = run_band_case(binary, workdir / f"f{k0:g}{ell}",
                                   field, "fem", f"b55f{k0:g}{ell}")
            records.append((k0, ell, spectral, exact, target, fem))
    if all(item[2] is None for item in records):
        if spectral_probe and "requires a build with USE_FFTW" in spectral_probe:
            print("[SKIP] B5.5 band bridge: binary lacks FFTW")
            return 0
        raise RuntimeError("B5.5 band bridge produced no spectral output")
    if any(item[2] is None for item in records):
        raise RuntimeError("B5.5 band bridge produced incomplete spectral output")
    for k0, ell, spectral, exact, target, fem in records:
        print(f"B5.5 k0={k0:.1f} l={ell}: spectral {spectral:+.6f} | "
              f"exact {exact:+.6f} | LSWT l-2<F> {target:+.6f} | FEM {fem:+.6f}")
    by_key = {(k0, ell): (spectral, exact, target)
              for k0, ell, spectral, exact, target, _ in records}
    errors = [abs(spectral - exact) for spectral, exact, _ in by_key.values()]
    shifts = [by_key[(k0, 1)][0] - by_key[(k0, 0)][0] for k0 in (1.0, 2.0)]
    lswt_errors = [abs(exact - target) for _, exact, target in by_key.values()]
    shift_error = max(abs(shift - 1.0) for shift in shifts)
    passed = (max(errors) <= 1.0e-5 and shift_error <= 1.0e-5
              and max(lswt_errors) <= 3.0e-3)
    return 0 if result("B5.5 band-packet bridge", passed,
                       f"max spectral error={max(errors):.3e}, "
                       f"max shift error={shift_error:.3e}, "
                       f"max LSWT gap={max(lswt_errors):.3e}") else 1


def write_band_dynamics_case(directory: Path, *, simid: str, ell: int,
                             gradient: str) -> None:
    """Write the undamped N=90 B5.5 band-packet conservation run."""

    directory.mkdir(parents=True, exist_ok=True)
    mkhoney.write(str(directory), C2, J=1.0, D=BAND_D,
                  kgrid=(30, 30), nphi=64, nr=16)
    coords = coords_honeycomb(BAND_N)
    center = coords[:2].mean(axis=1)
    field, _, _ = band_packet_with_gradient(coords, center, 1.0, ell)
    moments = np.vstack((field.real, field.imag,
                         np.sqrt(np.clip(1.0 - np.abs(field) ** 2, 0.0, None))))
    write_restart(directory, moments)
    origin = coords[:2].mean(axis=1)
    (directory / "inpsd.dat").write_text(f"""simid {simid}
ncell {BAND_N} {BAND_N} 1
BC P P 0
cell {C1[0]:.10f} {C1[1]:.10f} 0.0
     {C2[0]:.10f} {C2[1]:.10f} 0.0
     0.0 0.0 1.0
Sym 0
posfile ./posfile
posfiletype C
momfile ./momfile
exchange ./jfile
dm ./dmfile
maptype 1
do_prnstruct 0
initmag 4
restartfile ./restart.in
ip_mode N
mode S
temp 0.0
damping 0.0
Nstep 3001
timestep 1.0e-16
do_avrg N
do_tottraj N
do_oam_traj Y
oam_step 300
oam_buff 16
oam_gradient {gradient}
oam_axis 0 0 1
oam_sigma_max 0.6
oam_origin {origin[0]:.16e} {origin[1]:.16e} 0.0
do_chern N
do_oam_lswt N
""")


def run_band_dynamics(binary: str, workdir: Path) -> int:
    """Run the report-only conservative l=0/1 band-packet comparison."""

    trajectories = {}
    spectral_notice = ""
    for gradient in ("spectral", "fem"):
        for ell in (0, 1):
            simid = f"{'s' if gradient == 'spectral' else 'f'}{ell}"
            directory = workdir / simid
            write_band_dynamics_case(directory, simid=simid, ell=ell,
                                     gradient=gradient)
            run_binary(binary, directory)
            output = directory / f"oam_traj.{simid}.out"
            log = (directory / "uppasd.log").read_text()
            if not output.exists():
                if gradient == "spectral":
                    spectral_notice += log
                    trajectories[(gradient, ell)] = None
                    continue
                raise RuntimeError(f"No FEM OAM output in {directory}")
            rows = []
            for line in output.read_text().splitlines():
                if line.strip() and not line.lstrip().startswith("#"):
                    rows.append([float(value) for value in line.split()])
            trajectories[(gradient, ell)] = np.asarray(rows)
            data = trajectories[(gradient, ell)]
            print(f"band-dynamics {gradient} l={ell}: samples={len(data)}, "
                  f"lambda0={data[0, 2]:.7f}, lambda_range="
                  f"[{np.nanmin(data[:, 2]):.7f}, {np.nanmax(data[:, 2]):.7f}], "
                  f"Nm_relative_drift={np.max(np.abs(data[:, 3] / data[0, 3] - 1.0)):.3e}")

    if all(trajectories[("spectral", ell)] is None for ell in (0, 1)):
        if "requires a build with USE_FFTW" in spectral_notice:
            print("[SKIP] undamped band dynamics: binary lacks FFTW")
            return 0
        raise RuntimeError("spectral undamped band dynamics produced no output")
    if any(trajectories[("spectral", ell)] is None for ell in (0, 1)):
        raise RuntimeError("incomplete spectral undamped band dynamics output")

    spectral = [trajectories[("spectral", ell)] for ell in (0, 1)]
    if any(len(data) != len(spectral[0]) for data in spectral):
        raise RuntimeError("l=0 and l=1 spectral runs have different sample counts")
    shift_error = float(np.max(np.abs((spectral[1][:, 2] - spectral[0][:, 2]) - 1.0)))
    nm_drift = max(float(np.max(np.abs(data[:, 3] / data[0, 3] - 1.0)))
                   for data in spectral)
    spectral_ok = (all(np.max(np.abs(data[:, 2] - data[0, 2])) <= 3.0e-3
                       for data in spectral)
                   and abs(spectral[1][0, 2] - 1.0006) <= 3.0e-3
                   and shift_error <= 1.0e-3 and nm_drift <= 1.0e-4)
    result("undamped B5.5 spectral band dynamics", spectral_ok,
           f"lambda l1(t0)={spectral[1][0, 2]:.7f}, "
           f"max shift error={shift_error:.3e}, max Nm relative drift={nm_drift:.3e}")

    fem = [trajectories[("fem", ell)] for ell in (0, 1)]
    fem_bias = abs(float(np.mean(fem[1][:, 2])) - 0.880) <= 0.01
    fem_constant = all(np.max(np.abs(data[:, 2] - data[0, 2])) <= 3.0e-3
                       for data in fem)
    fem_ok = result("undamped B5.5 FEM control", fem_bias and fem_constant,
                    f"lambda l1 mean={np.mean(fem[1][:, 2]):.7f}, "
                    f"l1 range=[{np.min(fem[1][:, 2]):.7f}, {np.max(fem[1][:, 2]):.7f}]")
    print(f"Outputs retained in {workdir}")
    return 0 if spectral_ok and fem_ok else 1


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", default=os.environ.get("UPPASD", "uppasd"),
                        help="UppASD executable")
    parser.add_argument("--workdir", type=Path, default=Path("bridge_checks"),
                        help="generated run directory")
    parser.add_argument("--nstep", type=int, default=1000,
                        help="trajectory steps per packet")
    parser.add_argument("--force", action="store_true",
                        help="remove an existing generated work directory")
    parser.add_argument("--gradient", choices=("auto", "fem", "spectral"), default="fem",
                        help="trajectory OAM gradient method")
    parser.add_argument("--damping", type=float, default=0.0,
                        help="Gilbert damping for the generated packet")
    parser.add_argument("--temperature", type=float, default=0.0,
                        help="temperature for the generated packet")
    parser.add_argument("--band", action="store_true",
                        help="run the U5 B5.5 angular band-packet bridge check")
    parser.add_argument("--band-dynamics", action="store_true",
                        help="run the report-only undamped B5.5 band-packet conservation check")
    parser.add_argument("--clean", action="store_true",
                        help="run the one-step, grid-snapped clean control packet")
    parser.add_argument("--allow-known-bridge-gap", action="store_true",
                        help=argparse.SUPPRESS)
    args = parser.parse_args()

    workdir = args.workdir.resolve()
    if workdir.exists():
        if not args.force:
            parser.error(f"workdir exists; choose another directory or use --force: {workdir}")
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True)

    if args.band:
        return run_band_check(args.binary, workdir)
    if args.band_dynamics:
        return run_band_dynamics(args.binary, workdir)

    lswt_dir = workdir / "lswt"
    mkhoney.write(str(lswt_dir), C2, J=J, D=D, kgrid=(30, 30), nphi=64, nr=16,
                  extra="diamag_eps 1e-12")
    run_binary(args.binary, lswt_dir)
    k0, energy, fishman, mode = load_lswt_reference(lswt_dir)
    print(f"LSWT reference: k0={k0:.8f}, E={energy:.8f} meV, F_n={fishman:.8f}")

    if args.clean:
        n1 = 2 * round(k0 * CLEAN_N / (4.0 * np.pi))
        clean_k0 = 2.0 * np.pi * n1 / CLEAN_N
        _, clean_vectors = oracle_honey.evecs(np.array([clean_k0, 0.0]), J=J, D=D, S=1.0)
        clean_mode = clean_vectors[:, 0]
        print(f"Clean control: N={CLEAN_N}, n1={n1}, k0={clean_k0:.8f}, gradient={args.gradient}")
        clean_cases = {}
        for ell in (0, 1):
            simid = f"clean{ell}"
            directory = workdir / simid
            write_trajectory_case(directory, simid=simid, k0=clean_k0,
                                  mode=clean_mode, ell=ell, nstep=args.nstep,
                                  gradient=args.gradient, n=CLEAN_N, clean=True,
                                  damping=args.damping, temperature=args.temperature)
            run_binary(args.binary, directory)
            oam = load_oam(directory, simid)
            clean_cases[ell] = oam
            print(
                f"clean l={ell}: step1 lambda_centroid={oam[0, 2]:.8f}, "
                f"median lambda_centroid={np.nanmedian(oam[:, 2]):.8f}, "
                "phase residual=n/a"
            )
        clean_ok = True
        if args.gradient == "spectral":
            clean_ok = (abs(clean_cases[0][0, 2]) <= 1.0e-6
                        and abs(clean_cases[1][0, 2] - 1.0) <= 1.0e-4)
            result("U3 clean spectral targets", clean_ok,
                   f"lambda0={clean_cases[0][0, 2]:.8f}, "
                   f"lambda1={clean_cases[1][0, 2]:.8f}")
        else:
            print("[INFO] U3 clean FEM control is reported for comparison; spectral targets are not applied.")
        print(f"Outputs retained in {workdir}")
        return 0 if clean_ok else 1

    cases = {}
    for ell in (0, 1):
        simid = f"hb{ell}"
        directory = workdir / simid
        psi0 = write_trajectory_case(directory, simid=simid, k0=k0,
                                      mode=mode, ell=ell, nstep=args.nstep,
                                      gradient=args.gradient,
                                      damping=args.damping,
                                      temperature=args.temperature)
        run_binary(args.binary, directory)
        oam = load_oam(directory, simid)
        frequency, phase_residual, steps, projection = frequency_from_projection(
            directory, simid, psi0
        )
        cases[ell] = dict(oam=oam, frequency=frequency, phase_residual=phase_residual,
                          psi0=psi0, steps=steps, projection=projection)
        print(
            f"l={ell}: frequency={frequency:.8f} rad/ps, "
            f"phase residual={phase_residual:.3e}, "
            f"lambda_centroid step1={oam[0, 2]:.8f}, "
            f"lambda_centroid median={np.nanmedian(oam[:, 2]):.8f}, "
            f"N_m median={np.nanmedian(oam[:, 3]):.8f}"
        )

    expected_frequency = energy / HBAR_MEV_PS
    measured_frequency = cases[0]["frequency"]
    frequency_ok = abs(measured_frequency - expected_frequency) <= 0.01 * expected_frequency
    l_shift = np.nanmedian(cases[1]["oam"][:, 2]) - np.nanmedian(cases[0]["oam"][:, 2])
    bridge_value = np.nanmedian(cases[0]["oam"][:, 2])
    stable_ok = all(case["phase_residual"] < 0.15 for case in cases.values())

    # C14/B5.4 uses a full angular band packet and an analytic spatial
    # derivative. This isolates oracle/LSWT validation from the production
    # per-sublattice FEM diagnostic and its short-wavelength bias.
    bridge_results = [oracle_bridge.evaluate("periodic", N, k0, ell,
                                             D=D, sigma_k=0.18)
                      for ell in (0, 1)]
    bridge_errors = [result.lambda_exact_gradient - result.bridge_target
                     for result in bridge_results]
    bridge_shift = (bridge_results[1].lambda_exact_gradient -
                    bridge_results[0].lambda_exact_gradient)
    oracle_f_ok = abs(bridge_results[0].fishman_f - fishman) <= 0.01
    bridge_ok = (max(abs(error) for error in bridge_errors) <= 0.03 and
                 abs(bridge_shift - 1.0) <= 0.03 and oracle_f_ok)

    checks = [
        result("B5.4 frequency", frequency_ok,
               f"measured={measured_frequency:.8f}, expected={expected_frequency:.8f} rad/ps"),
        result("B5.4 particle-field bridge", bridge_ok,
               "; ".join([
                   f"F_n={bridge_results[0].fishman_f:.8f} (LSWT={fishman:.8f})",
                   f"lambda0={bridge_results[0].lambda_exact_gradient:.8f} "
                   f"target={bridge_results[0].bridge_target:.8f}",
                   f"lambda1={bridge_results[1].lambda_exact_gradient:.8f} "
                   f"target={bridge_results[1].bridge_target:.8f}",
                   f"shift={bridge_shift:.8f} (expected 1)"
               ])),
    ]
    print(
        "[INFO] Production fixed-spinor control targets: "
        "lambda(l=0)~=0, lambda(l=1)~=1, shift~=1; "
        f"observed phase_ok={stable_ok}, observed_shift={l_shift:.8f}, "
        f"lambda0={bridge_value:.8f}, gradient={args.gradient}; "
        "this fixed-spinor packet has no Berry term and is not used as the C14 bridge acceptance value"
    )
    print(f"Outputs retained in {workdir}")
    # The production fixed-spinor packet is intentionally not an acceptance
    # bridge packet.  B5.4 acceptance is the dynamics frequency plus the
    # independent, full-band particle-field convention oracle above.
    return 0 if frequency_ok and bridge_ok else 1


if __name__ == "__main__":
    sys.exit(main())
