"""Acceptance harness for the trajectory OAM (Blueprint B, WP5.2).

    UPPASD=/path/to/uppasd python3 run_traj_checks.py            # real run
    python3 run_traj_checks.py --selftest                         # harness only

--selftest writes each oam_traj file from the oracle instead of running
UppASD, which checks the harness itself.  It must pass before the harness is
used as evidence.

Output contract parsed here (CONVENTIONS_OAM.md, C6): oam_traj.<simid>.out,
'#' header lines, then whitespace columns
  step lambda_L_origin lambda_L_centroid N_m Lz_tot_hbar dSz_hbar balance R_x R_y sigma_psi
Non-finite values are written as NaN.
"""
import os, sys, shutil, subprocess
import numpy as np
import oracle_traj as O
import mkfixture_traj as F

EXE = os.environ.get("UPPASD", "uppasd")
SELFTEST = "--selftest" in sys.argv
WORK = os.path.abspath("traj_checks")
COLS = ["step", "lambda_L_origin", "lambda_L_centroid", "N_m", "Lz_tot_hbar",
        "dSz_hbar", "balance", "R_x", "R_y", "sigma_psi"]
N = 41
results = []


def reread(d):
    rs = np.loadtxt(f"{d}/restart.in", comments="#")
    return rs[:, 4:7].T


def run_case(name, lattice="square", bc=("P", "P"), shift=(0.0, 0.0), field=None, keys=None):
    d = os.path.join(WORK, name)
    shutil.rmtree(d, ignore_errors=True)
    keys = dict(keys or {})
    xy, _ = F.write_case(d, lattice=lattice, N=N, bc=bc, shift=shift, field=field, extra_keys=keys)
    m = reread(d)
    L = F.LATTICES[lattice]
    per = (bc[0] == "P", bc[1] == "P")
    simp, P = O.triangulate(N, N, 1, xy, L["C1"], L["C2"], per)
    origin = None
    if "oam_origin" in keys:
        origin = np.array([float(v) for v in keys["oam_origin"].split()[:2]])
    ref = O.evaluate(m, xy, simp, P, N, N, L["C1"], L["C2"], per, origin=origin,
                     weight=keys.get("oam_weight", "site"))
    fn = f"{d}/oam_traj.oamtest.out"
    if SELFTEST:
        with open(fn, "w") as f:
            f.write("# selftest\n# " + " ".join(COLS) + "\n")
            f.write("1 " + " ".join(repr(ref.get(c, float("nan"))) for c in COLS[1:]) + "\n")
        rc, log = 0, ""
    else:
        r = subprocess.run([EXE], cwd=d, capture_output=True, text=True)
        rc, log = r.returncode, r.stdout + r.stderr
        if rc == 0 and os.path.exists(f"{d}/coord.oamtest.out"):
            c = np.loadtxt(f"{d}/coord.oamtest.out")[:, 1:4].T
            assert np.abs(c - xy).max() < 1e-5, f"{name}: coordinate ordering differs from fixture"
    meas = None
    if rc == 0 and os.path.exists(fn):
        rows = [l.split() for l in open(fn) if l.strip() and not l.lstrip().startswith("#")]
        meas = dict(zip(COLS, [float(v) for v in rows[0]]))
    return ref, meas, rc, log


def check(label, ok, detail=""):
    results.append(ok)
    print(f"[{'PASS' if ok else 'FAIL'}] {label} {detail}")


def agree(ref, meas, key, tol=1e-6):
    a, b = ref[key], meas[key]
    if np.isnan(a) or np.isnan(b):
        return np.isnan(a) and np.isnan(b)
    return abs(a - b) <= tol * max(1.0, abs(a))


vort = lambda ell, amp=0.05, r0=5.0, c=None: (lambda xy: O.vortex(xy, xy[:2].mean(1) if c is None else c, ell, amp, r0))
os.makedirs(WORK, exist_ok=True)

# C1: angular index sweep, square, periodic
for ell in (-2, -1, 1, 2, 3):
    ref, meas, rc, log = run_case(f"ell{ell:+d}", field=vort(ell))
    ok = meas is not None and all(agree(ref, meas, k) for k in ("lambda_L_origin", "lambda_L_centroid", "N_m"))
    ok = ok and abs(meas["lambda_L_origin"] - ell) < 0.03 * abs(ell)  # FEM discretisation, r0=5
    check(f"C1 l={ell:+d}: matches oracle and l", ok, "" if meas is None else f"lambda={meas['lambda_L_origin']:.6f}")

# C2: amplitude independence
vals = []
for amp in (0.01, 0.05, 0.2):
    ref, meas, rc, log = run_case(f"amp{amp}", field=vort(1, amp))
    vals.append(np.nan if meas is None else meas["lambda_L_origin"])
vals = np.array(vals)
check("C2 amplitude independence", bool(np.all(np.isfinite(vals))) and float(np.ptp(vals)) < 1e-6, str(vals.tolist()))

# C3: rigid shift with fixed absolute origin.  The packet must carry net
# momentum (plane-wave boost), otherwise (R x P)_z = 0 and neither column moves.
def boosted(xy, k=0.3):
    v = O.vortex(xy, xy[:2].mean(1), 1, 0.05, 5.0)
    psi = (v[0] + 1j * v[1]) * np.exp(1j * k * xy[0])
    return np.stack([psi.real, psi.imag, v[2]])
out = {}
for tag, sh in (("a", (0.0, 0.0)), ("b", (0.37, 0.21))):
    ref, meas, rc, log = run_case(f"shift_{tag}", shift=sh, field=boosted, keys={"oam_origin": "0.0 0.0 0.0"})
    out[tag] = meas if (meas and agree(ref, meas, "lambda_L_origin") and agree(ref, meas, "lambda_L_centroid")) else None
ok = all(out.values()) and abs(out["a"]["lambda_L_centroid"] - out["b"]["lambda_L_centroid"]) < 1e-6 \
    and abs(out["a"]["lambda_L_origin"] - out["b"]["lambda_L_origin"]) > 1e-3
check("C3 centroid invariant, origin not, under rigid shift", ok)

# C4: saturated ferromagnet
ref, meas, rc, log = run_case("saturated", field=None)
check("C4 saturated FM -> NaN, no crash", rc == 0 and meas is not None and np.isnan(meas["lambda_L_origin"]))

# C5: packet straddling the periodic boundary
def straddle(xy):
    a = O.vortex(xy, (0.0, 20.0), 1, 0.05, 4.0); b = O.vortex(xy, (float(N), 20.0), 1, 0.05, 4.0)
    return np.where(np.abs(a[0] + 1j * a[1]) > np.abs(b[0] + 1j * b[1]), a, b)
ref, meas, rc, log = run_case("straddle", field=straddle)
ok = meas is not None and agree(ref, meas, "lambda_L_centroid") and min(abs(meas["R_x"]), abs(meas["R_x"] - N)) < 0.5
check("C5 boundary-straddling packet: centroid and lambda", ok)

# C6: site vs area weighting
r_s = run_case("w_site", field=vort(1), keys={"oam_weight": "site"})
r_a = run_case("w_area", field=vort(1), keys={"oam_weight": "area"})
ok = r_s[1] and r_a[1] and agree(r_a[0], r_a[1], "lambda_L_origin") \
    and abs(r_s[1]["lambda_L_origin"] - r_a[1]["lambda_L_origin"]) < 0.01
check("C6 site vs area agree within 1%, each matches oracle", bool(ok))

# C7: open vs periodic boundaries for a localised mode
r_o = run_case("bc_open", bc=("0", "0"), field=vort(1))
r_p = run_case("bc_per", bc=("P", "P"), field=vort(1))
ok = r_o[1] and r_p[1] and agree(r_o[0], r_o[1], "lambda_L_origin") \
    and abs(r_o[1]["lambda_L_origin"] - r_p[1]["lambda_L_origin"]) < 0.005
check("C7 open vs periodic within 0.5%", bool(ok))

# C8: oblique (hexagonal) cell
ref, meas, rc, log = run_case("hex", lattice="hex", field=vort(1))
ok = meas is not None and agree(ref, meas, "lambda_L_origin") and abs(meas["lambda_L_origin"] - 1) < 0.05
check("C8 hexagonal cell", ok)

# C9: delocalised field -> centroid column NaN, origin column finite
def noise(xy):
    rng = np.random.default_rng(7); v = np.zeros((3, xy.shape[1]))
    v[:2] = 0.05 * rng.normal(size=(2, xy.shape[1])); v[2] = np.sqrt(1 - v[0] ** 2 - v[1] ** 2); return v
ref, meas, rc, log = run_case("delocalised", field=noise)
check("C9 delocalised -> lambda_L_centroid NaN", meas is not None and np.isnan(meas["lambda_L_centroid"])
      and np.isfinite(meas["lambda_L_origin"]))

# C0: mesh diagnostic line (contract C7): periodic mesh tiles the cell exactly,
# open mesh drops the wrap cells.  Parsed from stdout.
if not SELFTEST:
    import re
    for name, bc, ncell in (("mesh_per", ("P", "P"), N * N), ("mesh_open", ("0", "0"), (N - 1) ** 2)):
        for lat in ("square", "hex"):
            ref, meas, rc, log = run_case(f"{name}_{lat}", lattice=lat, bc=bc, field=vort(1))
            L = F.LATTICES[lat]; cellA = abs(L["C1"][0] * L["C2"][1] - L["C1"][1] * L["C2"][0])
            mm = re.search(r"Mesh2D:\s*ntri=\s*(\d+)\s+total_area=\s*(\S+)\s+cell_area=\s*(\S+)\s+degenerate=\s*(\d+)", log)
            ok = bool(mm) and abs(float(mm.group(2)) - ncell * cellA) < 1e-8 * ncell and int(mm.group(4)) == 0
            check(f"C0 mesh {name} {lat}: total area = {ncell}*|C1xC2|, no degenerate", ok,
                  "" if not mm else f"(got {mm.group(2)})")

# C10: legacy key alias
if not SELFTEST:
    d = os.path.join(WORK, "legacy")
    shutil.rmtree(d, ignore_errors=True)
    F.write_case(d, field=vort(1), extra_keys={})
    txt = open(f"{d}/inpsd.dat").read().replace("do_oam_traj Y", "do_oam Y")
    open(f"{d}/inpsd.dat", "w").write(txt)
    r = subprocess.run([EXE], cwd=d, capture_output=True, text=True)
    check("C10 legacy do_oam -> trajectory OAM + deprecation warning",
          r.returncode == 0 and os.path.exists(f"{d}/oam_traj.oamtest.out") and "deprecat" in (r.stdout + r.stderr).lower())

print("ALL PASS" if all(results) else f"{results.count(False)} FAILED")
sys.exit(0 if all(results) else 1)
