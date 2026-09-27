"""Reproduce the LSWT / Fishman-OAM audit findings on UppASD v6.1.0rc (b60b762).

Usage:  python3 run_lswt_checks.py --binary /path/to/uppasd
Needs mkhoney.py and oracle_honey.py next to this file (numpy only).

Model: FM honeycomb, NN J=1 mRy, Haldane NNN DMI D=0.1 mRy along z, spins along z,
easy-axis anisotropy K1=-0.05 mRy (gaps the Goldstone mode so check 2 can run
even without the diamag_eps fix).

Checks
  1. Goldstone mode: isotropic FM + DMI must have omega(Gamma) = 0 exactly.
     Unpatched: gapped (0.09 meV at D=0.1) because DM vectors are mounted on
     exchange-list bonds.  With only the DM fix it aborts instead (see 3).
  2. Cell-description invariance: the same lattice described with
     C2=(0.5,sqrt3/2) and C2=(1.5,sqrt3/2) must give identical F_n(k).
     Unpatched: differ by up to ~20x (polar mesh passes reduced coords).
  3. Plain isotropic FM (D=0): must run.  Unpatched: ERROR STOP (paraunitarity)
     because diamag_eps defaults to 0 and Goldstone Colpa is singular.
  4. Oracle: F_n(k) must equal the ring Berry phase / (4 pi) from an independent
     Wilson-loop LSWT (|ratio| ~ 0.996 at oam_nphi=128; sign convention unpinned).
"""
import argparse, os, subprocess, numpy as np
import mkhoney as M, oracle_honey as O

EXE = os.environ.get("UPPASD", "uppasd")
s3 = np.sqrt(3.0)

def run(d):
    r = subprocess.run([EXE], cwd=d, capture_output=True, text=True)
    return r.returncode, r.stdout + r.stderr

def load_F(d, band=1):
    f = np.loadtxt(f"{d}/oam_lswt.honey.out", comments="#")
    f = f[f[:, 1] == band]
    return f[:, 0], f[:, 2], f[:, 3]

def main():
    global EXE
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", default=EXE, help="UppASD executable")
    EXE = parser.parse_args().binary

    ok = True

    # 1. Goldstone
    M.write("chk_goldstone", [0.5, s3/2, 0.0], D=0.1)
    rc, log = run("chk_goldstone")
    if rc != 0:
        print("[1] Goldstone: run aborted (see check 3):", [l for l in log.splitlines() if "STOP" in l][:1])
        ok = False
    else:
        k, E, _ = load_F("chk_goldstone", band=2)
        print(f"[1] Goldstone omega(Gamma) = {E[0]:.4e} meV  ->", "PASS" if E[0] < 1e-2 else "FAIL")
        ok &= E[0] < 1e-2

    # 2 + 4. invariance and oracle
    M.write_aniso("chk_A", [0.5, s3/2, 0.0], K1=-0.05, D=0.1)
    M.write_aniso("chk_C", [1.5, s3/2, 0.0], K1=-0.05, D=0.1)
    for d in ("chk_A", "chk_C"):
        rc, log = run(d); assert rc == 0, log[-500:]
    k, _, FA = load_F("chk_A"); _, _, FC = load_F("chk_C")
    dev = np.max(np.abs(FA - FC)) / np.max(np.abs(FA))
    print(f"[2] cell-description invariance: max|F_A-F_C|/max|F_A| = {dev:.2e}  ->", "PASS" if dev < 1e-8 else "FAIL")
    ok &= dev < 1e-8
    Ft = O.F_of_k(k, 0, D=0.1)
    nz = np.abs(Ft) > 1e-10
    ratio = np.abs(FA[nz]) / np.abs(Ft[nz])
    print(f"[4] |F_uppasd|/|F_oracle| in [{ratio.min():.4f}, {ratio.max():.4f}]  ->",
          "PASS" if (ratio.min() > 0.99 and ratio.max() < 1.01) else "FAIL")
    ok &= ratio.min() > 0.99 and ratio.max() < 1.01

    # 3. plain isotropic FM
    M.write("chk_D0", [0.5, s3/2, 0.0], D=0.0)
    rc, log = run("chk_D0")
    print("[3] isotropic FM (D=0) runs:", "PASS" if rc == 0 else "FAIL (" + next((l for l in log.splitlines() if "STOP" in l), "abort") + ")")
    ok &= rc == 0

    print("ALL PASS" if ok else "SOME CHECKS FAILED")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
