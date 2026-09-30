"""Audit reference for U5 B5.5: run the production trajectory kernel on the angularly varying
honeycomb band packet and compare it with the exact-gradient value of l - 2<F>.

This is the script used in the audit of 877fcbb. It is a REFERENCE: port the logic into
run_bridge_checks.py (flag --band) and do not import from docs/.

usage: python3 b55_band_packet_reference.py /path/to/uppasd   (run from tests/SpinWaves/oam_lswt)

Audit output (FFTW Release build, 877fcbb):
  spectral k0=1.0 l=0: UppASD +0.000606 | exact +0.000606 | diff +5.3e-08
  spectral k0=1.0 l=1: UppASD +1.000602 | exact +1.000606 | diff -4.5e-06
  spectral k0=2.0 l=0: UppASD +0.040822 | exact +0.040822 | diff +2.4e-07
  spectral k0=2.0 l=1: UppASD +1.040821 | exact +1.040822 | diff -1.1e-06
  fem      k0=1.0 l=1: UppASD +0.882160 ;  fem k0=2.0 l=1: UppASD +0.601644
"""
import os, shutil, subprocess, sys
import numpy as np
import oracle_honey as H
import oracle_traj as O

EXE = sys.argv[1]
S3 = np.sqrt(3.0)
C1 = (1.0, 0.0, 0.0); C2 = (0.5, S3 / 2, 0.0)
BASIS = [(0.0, 0.0, 0.0), (0.5, 0.5 / S3, 0.0)]
N = 90                    # packet width ~ 1/dk ~ 7a; wrap amplitude ~ 1e-10
D = 1.0; BAND = 0         # upper (optical) band, J = S = 1
DK = 0.15; NPHI = 192; NR = 13
xy = O.build_coords(N, N, C1, C2, BASIS)      # basis fastest, UppASD order
ctr = xy[:2].mean(1)
cell = np.array([N * np.array(C1[:2]), N * np.array(C2[:2])]).T
inv = np.linalg.inv(cell)


def ring_gauge(k):
    """Band spinors around a ring in a smooth (parallel-transported) gauge."""
    phis = 2 * np.pi * np.arange(NPHI) / NPHI
    us = [H.evecs(np.array([k * np.cos(p), k * np.sin(p)]), D=D)[1][:, BAND] for p in phis]
    for j in range(1, NPHI):
        ov = np.vdot(us[j - 1], us[j]); us[j] = us[j] * np.conj(ov) / abs(ov)
    th = np.angle(np.vdot(us[-1], us[0]))
    return phis, [u * np.exp(1j * j * th / NPHI) for j, u in enumerate(us)]


def packet_with_grad(k0, ell):
    """psi(r) = sum_k w(|k|) e^{i l phi_k} u_s(k) e^{i k.r} and its exact gradient."""
    krs = np.linspace(k0 - 2.5 * DK, k0 + 2.5 * DK, NR)
    r = xy[:2] - ctr[:, None]; sub = np.arange(xy.shape[1]) % 2
    f = np.zeros(xy.shape[1], complex); gx = f.copy(); gy = f.copy(); prev = None
    for kr in krs:
        phis, us = ring_gauge(kr)
        if prev is not None:                       # radial continuity of the phi=0 phase
            ov = np.vdot(prev, us[0]); us = [u * np.conj(ov) / abs(ov) for u in us]
        prev = us[0]; w = np.exp(-(kr - k0) ** 2 / (2 * DK ** 2))
        for p, u in zip(phis, us):
            kv = kr * np.array([np.cos(p), np.sin(p)])
            t = w * np.exp(1j * ell * p) * u[sub] * np.exp(1j * (kv @ r))
            f += t; gx += 1j * kv[0] * t; gy += 1j * kv[1] * t
    s = 0.05 / np.abs(f).max()
    return f * s, gx * s, gy * s


def exact_lambda_centroid(f, gx, gy):
    """C3 lambda_L_centroid (site weights, circular centroid, minimum image) with exact gradient."""
    w = np.abs(f) ** 2; th = 2 * np.pi * ((inv @ xy[:2]) % 1)
    R = cell @ np.mod(np.arctan2((w * np.sin(th)).sum(1), (w * np.cos(th)).sum(1)) / (2 * np.pi), 1)
    d = inv @ (xy[:2] - R[:, None]); d -= np.round(d); lev = cell @ d
    return np.sum(np.imag(np.conj(f) * (lev[0] * gy - lev[1] * gx))) / w.sum()


def lswt_prediction(k0, ell):
    """l - 2<F>, <F> weighted by |w(k)|^2 k over the radial profile."""
    krs = np.linspace(k0 - 2.5 * DK, k0 + 2.5 * DK, NR)
    kk = np.linspace(0.0, krs[-1], 400)
    F = H.F_of_k(kk, BAND, D=D, nphi=384)
    Fk = np.interp(krs, kk, F)
    wt = np.exp(-(krs - k0) ** 2 / DK ** 2) * krs
    return ell - 2 * np.sum(wt * Fk) / wt.sum()


def run_uppasd(f, grad, tag):
    dd = f"b55_{tag}"; shutil.rmtree(dd, ignore_errors=True); os.makedirs(dd)
    tau = [np.array(b[:2]) for b in BASIS]
    open(f"{dd}/posfile", "w").write("".join(f"{i+1} {i+1} {t[0]:.14f} {t[1]:.14f} 0.0\n" for i, t in enumerate(tau)))
    open(f"{dd}/momfile", "w").write("1 1 1.0 0.0 0.0 1.0\n2 1 1.0 0.0 0.0 1.0\n")
    nn = [(0.5, S3 / 6), (-0.5, S3 / 6), (0.0, -S3 / 3)]
    open(f"{dd}/jfile", "w").write("".join(f"1 2 {v[0]:.14f} {v[1]:.14f} 0.0 1.0\n" for v in nn) +
                                   "".join(f"2 1 {-v[0]:.14f} {-v[1]:.14f} 0.0 1.0\n" for v in nn))
    mz = np.sqrt(1 - np.abs(f) ** 2); n = f.size
    with open(f"{dd}/restart.in", "w") as fh:
        fh.write("#" * 80 + f"\n# File type: R\n# Simulation type: S\n# Number of atoms: {n:9d}\n# Number of ensembles:         1\n" + "#" * 80 + "\n")
        fh.write("  # iter     ens   iatom           |Mom|             M_x             M_y             M_z\n")
        for i in range(n):
            fh.write(f"{0:8d}{1:8d}{i+1:8d}  {1.0:16.8E}{f[i].real:24.16E}{f[i].imag:24.16E}{mz[i]:24.16E}\n")
    open(f"{dd}/inpsd.dat", "w").write(f"""simid b55
ncell {N} {N} 1
BC P P 0
cell 1.0 0.0 0.0
     0.5 {S3/2:.14f} 0.0
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
oam_gradient {grad}
oam_axis 0 0 1
""")
    subprocess.run([EXE], cwd=dd, capture_output=True)
    c = np.loadtxt(f"{dd}/coord.b55.out")[:, 1:3].T
    assert np.abs(c - xy[:2]).max() < 1e-5, "coordinate ordering differs"
    val = np.loadtxt(f"{dd}/oam_traj.b55.out", comments="#", ndmin=2)[0][2]
    shutil.rmtree(dd)
    return val


for k0 in (1.0, 2.0):
    for ell in (0, 1):
        f, gx, gy = packet_with_grad(k0, ell)
        ex = exact_lambda_centroid(f, gx, gy); pred = lswt_prediction(k0, ell)
        sp = run_uppasd(f, "spectral", f"s{k0}{ell}"); fe = run_uppasd(f, "fem", f"f{k0}{ell}")
        print(f"k0={k0} l={ell}: spectral {sp:+.6f} | exact {ex:+.6f} (diff {sp-ex:+.1e}) | "
              f"LSWT l-2<F> {pred:+.6f} (exact-LSWT {ex-pred:+.1e}) | fem {fe:+.6f}")
