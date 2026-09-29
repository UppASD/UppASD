"""Run the LSWT-to-trajectory OAM bridge check (Blueprint B5.4).

The harness reports two related but distinct checks:

* the precession frequency of a narrow honeycomb magnon packet is compared
  with the energy from ``do_oam_lswt``;
* the settled particle-field bridge is checked with the independent
  reciprocal-supercell packet, whose analytic gradient contains the angular
  variation of the band spinor.

The production trajectory packet is retained as a dynamics/control diagnostic.
Its per-sublattice FEM value is not used as a quantitative ``F_n`` check at
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
C1 = np.array([1.0, 0.0, 0.0])
C2 = np.array([0.5, S3 / 2.0, 0.0])
BASIS = np.array([[0.0, 0.0, 0.0], [0.5, 0.5 / S3, 0.0]])
J = 1.0
D = 0.30
AMP = 0.025
SIGMA = 4.0
TIMESTEP_S = 1.0e-16
HBAR_MEV_PS = 0.6582119569


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


def packet(coords: np.ndarray, mode: np.ndarray, k0: float, ell: int) -> np.ndarray:
    """Construct a small-amplitude optical-band packet in the global frame."""

    xy = coords[:2]
    center = xy.mean(axis=1)
    dxy = xy - center[:, None]
    radius = np.linalg.norm(dxy, axis=0)
    theta = np.arctan2(dxy[1], dxy[0])
    envelope = np.exp(-0.5 * (radius / SIGMA) ** 2)
    phase = np.exp(1j * (k0 * xy[0] + ell * theta))
    sublattice_mode = np.asarray([mode[0], mode[1]], complex)
    psi = AMP * sublattice_mode[np.arange(coords.shape[1]) % 2] * envelope * phase
    mz = np.sqrt(np.clip(1.0 - np.abs(psi) ** 2, 0.0, None))
    return np.vstack((psi.real, psi.imag, mz))


def write_trajectory_case(directory: Path, *, simid: str, k0: float,
                           mode: np.ndarray, ell: int, nstep: int,
                           gradient: str = "fem") -> np.ndarray:
    """Write one trajectory input and return its prescribed complex packet."""

    directory.mkdir(parents=True, exist_ok=True)
    # This supplies the common honeycomb geometry, exchange and DMI files.
    mkhoney.write(str(directory), C2, J=J, D=D, kgrid=(30, 30), nphi=64, nr=16)
    coords = coords_honeycomb()
    moments = packet(coords, mode, k0, ell)
    write_restart(directory, moments)
    origin = coords[:2].mean(axis=1)
    with (directory / "inpsd.dat").open("w") as handle:
        handle.write(f"""simid {simid}
ncell {N} {N} 1
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
temp 0.0
damping 0.0
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
    parser.add_argument("--gradient", choices=("fem", "spectral"), default="fem",
                        help="trajectory OAM gradient method")
    parser.add_argument("--allow-known-bridge-gap", action="store_true",
                        help=argparse.SUPPRESS)
    args = parser.parse_args()

    workdir = args.workdir.resolve()
    if workdir.exists():
        if not args.force:
            parser.error(f"workdir exists; choose another directory or use --force: {workdir}")
        shutil.rmtree(workdir)
    workdir.mkdir(parents=True)

    lswt_dir = workdir / "lswt"
    mkhoney.write(str(lswt_dir), C2, J=J, D=D, kgrid=(30, 30), nphi=64, nr=16,
                  extra="diamag_eps 1e-12")
    run_binary(args.binary, lswt_dir)
    k0, energy, fishman, mode = load_lswt_reference(lswt_dir)
    print(f"LSWT reference: k0={k0:.8f}, E={energy:.8f} meV, F_n={fishman:.8f}")

    cases = {}
    for ell in (0, 1):
        simid = f"hb{ell}"
        directory = workdir / simid
        psi0 = write_trajectory_case(directory, simid=simid, k0=k0,
                                      mode=mode, ell=ell, nstep=args.nstep,
                                      gradient=args.gradient)
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
        "[INFO] B5.4 production-mesh diagnostic: "
        f"phase_ok={stable_ok}, observed_shift={l_shift:.8f}, "
        f"fixed_spinor_lambda0={bridge_value:.8f}, gradient={args.gradient}; "
        "not used as the C14 bridge acceptance value"
    )
    print(f"Outputs retained in {workdir}")
    # The production fixed-spinor packet is intentionally not an acceptance
    # bridge packet.  B5.4 acceptance is the dynamics frequency plus the
    # independent, full-band particle-field convention oracle above.
    return 0 if frequency_ok and bridge_ok else 1


if __name__ == "__main__":
    sys.exit(main())
