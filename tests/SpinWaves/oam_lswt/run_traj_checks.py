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
import argparse, os, sys, shutil, subprocess
import numpy as np
import oracle_traj as O
import mkfixture_traj as F

_parser = argparse.ArgumentParser()
_parser.add_argument("--binary", default=os.environ.get("UPPASD", "uppasd"),
                     help="UppASD executable")
_parser.add_argument("--selftest", action="store_true",
                     help="write oracle outputs instead of running UppASD")
_args = _parser.parse_args()
EXE = _args.binary
SELFTEST = _args.selftest
WORK = os.path.abspath("traj_checks")
COLS = ["step", "lambda_L_origin", "lambda_L_centroid", "N_m", "Lz_tot_hbar",
        "dSz_hbar", "balance", "R_x", "R_y", "sigma_psi"]
N = 41
results = []


def reread(d):
    rs = np.loadtxt(f"{d}/restart.in", comments="#")
    rs = np.atleast_2d(rs)
    ensembles = int(np.max(rs[:, 1]))
    return np.stack([rs[rs[:, 1] == k, 4:7].T for k in range(1, ensembles + 1)], axis=2)


def run_case(name, lattice="square", bc=("P", "P"), shift=(0.0, 0.0), field=None,
             keys=None, layers=1, grid=None, ensemble_fields=None):
    d = os.path.join(WORK, name)
    shutil.rmtree(d, ignore_errors=True)
    keys = dict(keys or {})
    keys.setdefault("oam_gradient", "fem")
    ngrid = N if grid is None else grid
    xy, _ = F.write_case(d, lattice=lattice, N=ngrid, bc=bc, shift=shift, field=field,
                          extra_keys=keys, layers=layers, ensemble_fields=ensemble_fields)
    m = reread(d)
    L = F.LATTICES[lattice]
    per = (bc[0] == "P", bc[1] == "P")
    layer_n = ngrid * ngrid
    simp1, P1 = O.triangulate(ngrid, ngrid, 1, xy[:, :layer_n], L["C1"], L["C2"], per)
    simp = np.concatenate([simp1 + z * layer_n for z in range(layers)], axis=0)
    P = np.concatenate([P1 for _ in range(layers)], axis=0)
    origin = None
    if "oam_origin" in keys:
        origin = np.array([float(v) for v in keys["oam_origin"].split()[:2]])
    axis = None
    if "oam_axis" in keys:
        axis = np.array([float(v) for v in keys["oam_axis"].split()[:3]])
    gradient = keys.get("oam_gradient", "fem")
    if gradient == "auto":
        gradient = "spectral" if all(per) else "fem"
    if m.shape[2] == 1:
        ref = O.evaluate(m[:, :, 0], xy, simp, P, ngrid, ngrid, L["C1"], L["C2"], per,
                         origin=origin, weight=keys.get("oam_weight", "site"),
                         gradient=gradient, NA=1, axis=axis)
    else:
        ref = O.evaluate_ensembles([m[:, :, k] for k in range(m.shape[2])], xy, simp, P,
                                   ngrid, ngrid, L["C1"], L["C2"], per, origin=origin,
                                   weight=keys.get("oam_weight", "site"), gradient=gradient,
                                   NA=1, axis=axis)
    fn = f"{d}/oam_traj.oamtest.out"
    if SELFTEST:
        with open(fn, "w") as f:
            resolved = "spectral" if keys.get("oam_gradient") == "auto" and all(per) else gradient
            f.write(f"# selftest\n# oam_gradient = {resolved} (auto)\n# " + " ".join(COLS) + "\n")
            f.write("1 " + " ".join(repr(ref.get(c, float("nan"))) for c in COLS[1:]) + "\n")
        excluded = (ref.get("_excluded_origin", 0), ref.get("_excluded_centroid", 0))
        log = "excluded ensembles" if any(excluded) else ""
        rc = 0
    else:
        r = subprocess.run([EXE], cwd=d, capture_output=True, text=True)
        rc, log = r.returncode, r.stdout + r.stderr
        if rc == 0 and os.path.exists(f"{d}/coord.oamtest.out"):
            c = np.loadtxt(f"{d}/coord.oamtest.out")[:, 1:4].T
            assert np.abs(c - xy).max() < 1e-5, f"{name}: coordinate ordering differs from fixture"
        if rc == 0 and keys.get("oam_gradient") == "auto":
            header = open(fn).read() if os.path.exists(fn) else ""
            gradient = "spectral" if "oam_gradient = spectral (auto)" in header else "fem"
            if m.shape[2] == 1:
                ref = O.evaluate(m[:, :, 0], xy, simp, P, ngrid, ngrid, L["C1"], L["C2"], per,
                                 origin=origin, weight=keys.get("oam_weight", "site"),
                                 gradient=gradient, NA=1, axis=axis)
            else:
                ref = O.evaluate_ensembles([m[:, :, k] for k in range(m.shape[2])], xy, simp, P,
                                           ngrid, ngrid, L["C1"], L["C2"], per, origin=origin,
                                           weight=keys.get("oam_weight", "site"), gradient=gradient,
                                           NA=1, axis=axis)
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
    ref, meas, rc, log = run_case(f"shift_{tag}", shift=sh, field=boosted,
                                  keys={"oam_origin": "0.0 0.0 0.0"})
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

# C11: independent layers preserve lambda and add their magnetic weight.
r_one = run_case("multilayer_1", field=vort(1), layers=1)
r_two = run_case("multilayer_2", field=vort(1), layers=2)
ok = r_one[1] is not None and r_two[1] is not None
if ok:
    ok = all(agree(r_two[0], r_two[1], k) for k in ("lambda_L_origin", "lambda_L_centroid", "N_m"))
    ok = ok and abs(r_two[1]["lambda_L_origin"] - r_one[1]["lambda_L_origin"]) <= 1e-6
    ok = ok and abs(r_two[1]["N_m"] - 2.0 * r_one[1]["N_m"]) <= 1e-6 * max(1.0, abs(r_one[1]["N_m"]))
check("C11 two layers: lambda unchanged, N_m doubled", ok)

# C12: FFTW spectral gradients.  The first probe distinguishes a binary
# without USE_FFTW from an implementation failure; the former is explicitly
# outside this check's build requirement.
def high_k_boost(xy, k=1.5):
    center = xy[:2].mean(1)
    v = O.vortex(xy, center, 1, 0.05, 5.0)
    psi = (v[0] + 1j * v[1]) * np.exp(1j * k * xy[0])
    return np.stack([psi.real, psi.imag, np.sqrt(np.clip(1 - np.abs(psi) ** 2, 0, None))])


spectral_probe = run_case("spectral_probe", field=vort(1),
                          keys={"oam_gradient": "spectral"})
spectral_header = os.path.join(WORK, "spectral_probe", "oam_traj.oamtest.out")
has_fftw = SELFTEST or (
    spectral_probe[1] is not None
    and os.path.exists(spectral_header)
    and any("oam_gradient = spectral" in line for line in open(spectral_header))
)
if not has_fftw and not SELFTEST and "requires a build with USE_FFTW" in spectral_probe[3]:
    print("[SKIP] C12 spectral gradient: binary was not built with FFTW")
else:
    r_ell = run_case("spectral_ell1", field=vort(1),
                     keys={"oam_gradient": "spectral"})
    r_hex = run_case("spectral_hex", lattice="hex", field=vort(1),
                     keys={"oam_gradient": "spectral"})
    ok_a = all(r[1] is not None and all(agree(r[0], r[1], key, tol=1e-6)
                                        for key in ("lambda_L_origin", "lambda_L_centroid", "N_m"))
               for r in (r_ell, r_hex))
    check("C12a spectral ell+1 and hex match NumPy oracle", ok_a)

    r_high_s = run_case("spectral_high_k", field=high_k_boost,
                        keys={"oam_gradient": "spectral"})
    r_high_f = run_case("fem_high_k", field=high_k_boost)
    analytic_lambda = 1.0
    ok_b = (r_high_s[1] is not None and r_high_f[1] is not None
            and abs(r_high_s[1]["lambda_L_centroid"] - analytic_lambda) <= 1e-4
            and abs(r_high_f[1]["lambda_L_centroid"] - analytic_lambda) >= 0.1)
    check("C12b k.a=1.5 spectral analytic lambda; FEM bias", ok_b,
          "" if r_high_s[1] is None or r_high_f[1] is None else
          f"spectral={r_high_s[1]['lambda_L_centroid']:.8f}, "
          f"FEM={r_high_f[1]['lambda_L_centroid']:.8f}, analytic={analytic_lambda:.1f}")

    r_high_hex = run_case("spectral_high_k_hex", lattice="hex", grid=48,
                          field=lambda xy: high_k_boost(xy, k=3.0),
                          keys={"oam_gradient": "spectral", "oam_axis": "0 0 1"})
    ok_c = (r_high_hex[1] is not None
            and all(agree(r_high_hex[0], r_high_hex[1], key, tol=1e-6)
                    for key in ("lambda_L_origin", "lambda_L_centroid", "N_m"))
            and abs(r_high_hex[1]["lambda_L_centroid"] - 1.0) <= 1e-4)
    check("C12c hex k.a=3 Brillouin-zone fold", ok_c,
          "" if r_high_hex[1] is None else
          f"spectral={r_high_hex[1]['lambda_L_centroid']:.8f}, "
          f"oracle={r_high_hex[0]['lambda_L_centroid']:.8f}, analytic=1.0")

# C13: explicit lab-frame axis versus the C8 default frame.  This uses the
# original boosted packet; T4 restores it without subtracting its mean.
c13_keys = {"oam_origin": "0.0 0.0 0.0"}
r_c13_default = run_case("axis_default", field=boosted, keys=c13_keys)
r_c13_lab = run_case("axis_lab", field=boosted,
                     keys={**c13_keys, "oam_axis": "0.0 0.0 1.0"})
ok_default = (r_c13_default[1] is not None
              and all(agree(r_c13_default[0], r_c13_default[1], key)
                      for key in ("lambda_L_origin", "lambda_L_centroid", "N_m")))
ok_lab = (r_c13_lab[1] is not None
          and all(agree(r_c13_lab[0], r_c13_lab[1], key)
                  for key in ("lambda_L_origin", "lambda_L_centroid", "N_m")))
check("C13 default axis matches default-frame oracle", ok_default)
check("C13 explicit z axis matches lab-frame oracle", ok_lab)

# C14: norm-weighted ensemble aggregation.  The second packet has a different
# amplitude so arithmetic averaging of already-normalized lambda values would
# fail this check.
ensemble_vortices = run_case(
    "ensemble_vortices", field=None,
    ensemble_fields=[vort(1, amp=0.05), vort(2, amp=0.08)],
    keys={"oam_gradient": "fem"})
ensemble_meas = ensemble_vortices[1]
ensemble_ref = ensemble_vortices[0]
ensemble_keys = ("lambda_L_origin", "lambda_L_centroid", "N_m",
                 "Lz_tot_hbar", "dSz_hbar", "balance", "R_x", "R_y", "sigma_psi")
ensemble_ok = ensemble_meas is not None and all(
    (np.isnan(ensemble_ref[key]) and np.isnan(ensemble_meas[key])) or
    abs(ensemble_ref[key] - ensemble_meas[key]) <= 1.0e-10 * max(1.0, abs(ensemble_ref[key]))
    for key in ensemble_keys)
check("C14 norm-weighted ensemble aggregation", ensemble_ok)

ensemble_delocalised = run_case(
    "ensemble_delocalised", field=None,
    ensemble_fields=[vort(1, amp=0.05), noise],
    keys={"oam_gradient": "fem", "oam_axis": "0 0 1"})
deloc_meas = ensemble_delocalised[1]
deloc_ref = ensemble_delocalised[0]
deloc_log = ensemble_delocalised[3]
deloc_ok = (deloc_meas is not None and np.isfinite(deloc_meas["lambda_L_origin"])
            and np.isfinite(deloc_meas["lambda_L_centroid"])
            and agree(deloc_ref, deloc_meas, "lambda_L_centroid", tol=1.0e-10)
            and "excluded ensembles" in deloc_log)
check("C14 delocalised ensemble: centroid excludes one with warning", deloc_ok)

# C15: automatic gradient resolution.  A periodic FFTW binary must select
# spectral; open boundaries always select FEM.  The periodic FEM result is
# also the non-FFTW branch, which must not refuse.
auto_periodic = run_case("auto_periodic", field=vort(1),
                         keys={"oam_gradient": "auto"})
auto_open = run_case("auto_open", bc=("0", "0"), field=vort(1),
                     keys={"oam_gradient": "auto"})
auto_periodic_header = os.path.join(WORK, "auto_periodic", "oam_traj.oamtest.out")
auto_open_header = os.path.join(WORK, "auto_open", "oam_traj.oamtest.out")
auto_periodic_text = open(auto_periodic_header).read() if os.path.exists(auto_periodic_header) else ""
auto_open_text = open(auto_open_header).read() if os.path.exists(auto_open_header) else ""
auto_periodic_method = "spectral" if "oam_gradient = spectral (auto)" in auto_periodic_text else "fem"
auto_open_ok = (auto_open[1] is not None and "oam_gradient = fem (auto)" in auto_open_text
                and "disabled" not in auto_open[3].lower())
if auto_periodic_method == "spectral":
    auto_periodic_ok = (auto_periodic[1] is not None
                        and all(agree(auto_periodic[0], auto_periodic[1], key, tol=1.0e-10)
                                for key in ("lambda_L_origin", "lambda_L_centroid", "N_m")))
else:
    auto_periodic_ok = (auto_periodic[1] is not None
                        and "oam_gradient = fem (auto)" in auto_periodic_text
                        and "requires a build with USE_FFTW" not in auto_periodic[3])
check("C15 auto periodic resolves without refusal", auto_periodic_ok,
      f"resolved={auto_periodic_method}")
check("C15 auto open resolves to FEM", auto_open_ok)
check("C15 auto non-FFTW branch has no refusal", auto_periodic_ok and
      (auto_periodic_method == "fem" or "requires a build with USE_FFTW" not in auto_periodic[3]))

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
