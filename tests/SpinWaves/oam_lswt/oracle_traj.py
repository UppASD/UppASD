"""Independent reference for the UppASD trajectory OAM (pyswatter formulation).

This file is the oracle the Fortran implementation is pinned against.  It is
deliberately written without reference to the Fortran and must never be
"fixed" to agree with it: if they disagree, escalate.

Physics (CONVENTIONS_OAM.md, C1-C4):
    psi_i      = m_x,i + i m_y,i                 (frame: global, from <m> at init)
    grad psi   = linear-FEM triangle gradients, area-weighted to sites
    ell_z(i)   = Im[ psi_i^* ((x_i-x0) d_y psi_i - (y_i-y0) d_x psi_i) ]
    lambda_L   = sum_i ell_z(i) w_i / sum_i |psi_i|^2 w_i       (w_i = 1 or A_i)
    R          = |psi|^2 w-weighted centroid (circular mean along periodic axes)
    N_m        = sum_i (mmom_i/g) (1 - m_z,i)

Oracle identity used by the tests: for psi = g(r) exp(i l phi) about the
origin, (x d_y - y d_x) psi = i l psi, so lambda_L = l exactly for any
envelope and amplitude (up to discretisation).
"""
import json
import numpy as np


# ------------------------------------------------------------------ lattice --
def build_coords(N1, N2, C1, C2, basis):
    """UppASD atom ordering: atom = NA*(cell-1)+ia, cell index x fastest."""
    C1 = np.asarray(C1, float); C2 = np.asarray(C2, float)
    basis = np.asarray(basis, float).reshape(-1, 3)
    out = []
    for iy in range(N2):
        for ix in range(N1):
            for b in basis:
                out.append(ix * C1 + iy * C2 + b)
    return np.array(out).T                                    # (3, Natom)


def triangulate(N1, N2, NA, coords, C1, C2, periodic=(True, True)):
    """Port of delaunay_tri_tri with the WP1 fixes: minimum-image vertices,
    open boundaries drop wrap cells, counter-clockwise orientation.
    Returns simp (ntri,3) and local vertex positions P (ntri,3,2) relative
    to vertex 0, unwrapped."""
    C1 = np.asarray(C1, float)[:2]; C2 = np.asarray(C2, float)[:2]
    L1, L2 = N1 * C1, N2 * C2
    M = np.array([L1, L2]).T                                  # cell matrix (2x2)
    Minv = np.linalg.inv(M)

    def mi(d):                                                # minimum image, 2D
        s = Minv @ d
        s[0] -= np.round(s[0]) if periodic[0] else 0.0
        s[1] -= np.round(s[1]) if periodic[1] else 0.0
        return M @ s

    simp, P = [], []
    for iy in range(N2):
        for ix in range(N1):
            ixp, iyp = ix + 1, iy + 1
            if ixp == N1:
                if not periodic[0]:
                    continue
                ixp = 0
            if iyp == N2:
                if not periodic[1]:
                    continue
                iyp = 0
            for it in range(NA):
                i00 = NA * (iy * N1 + ix) + it
                i10 = NA * (iy * N1 + ixp) + it
                i01 = NA * (iyp * N1 + ix) + it
                i11 = NA * (iyp * N1 + ixp) + it
                r = {k: coords[:2, k] for k in (i00, i10, i01, i11)}
                rel = {k: mi(r[k] - r[i00]) for k in r}
                d1 = np.sum((rel[i00] - rel[i11]) ** 2)
                d2 = np.sum((rel[i10] - rel[i01]) ** 2)
                tris = ([(i00, i10, i11), (i00, i11, i01)] if d1 <= d2
                        else [(i00, i10, i01), (i10, i11, i01)])
                for t in tris:
                    p = np.array([rel[k] - rel[t[0]] for k in t])
                    a2 = (p[1, 0] - p[0, 0]) * (p[2, 1] - p[0, 1]) - (p[1, 1] - p[0, 1]) * (p[2, 0] - p[0, 0])
                    if a2 < 0:                                # enforce CCW
                        t = (t[0], t[2], t[1]); p = p[[0, 2, 1]]
                    simp.append(t); P.append(p)
    return np.array(simp), np.array(P)


def fem_gradient(psi, simp, P):
    """Linear-FEM gradient per triangle, area-weighted to sites."""
    x, y = P[:, :, 0], P[:, :, 1]
    A2 = (x[:, 1] - x[:, 0]) * (y[:, 2] - y[:, 0]) - (y[:, 1] - y[:, 0]) * (x[:, 2] - x[:, 0])
    A = 0.5 * A2
    b = np.stack([y[:, 1] - y[:, 2], y[:, 2] - y[:, 0], y[:, 0] - y[:, 1]], 1) / A2[:, None]
    c = np.stack([x[:, 2] - x[:, 1], x[:, 0] - x[:, 2], x[:, 1] - x[:, 0]], 1) / A2[:, None]
    pv = psi[simp]
    gx_t = np.sum(b * pv, 1); gy_t = np.sum(c * pv, 1)
    n = psi.size
    gx = np.zeros(n, complex); gy = np.zeros(n, complex); ws = np.zeros(n)
    for v in range(3):
        np.add.at(gx, simp[:, v], A * gx_t)
        np.add.at(gy, simp[:, v], A * gy_t)
        np.add.at(ws, simp[:, v], A)
    ok = ws > 0
    gx[ok] /= ws[ok]; gy[ok] /= ws[ok]
    return gx, gy, ws / 3.0, ok, A


def circular_centroid(coords, wt, N1, N2, C1, C2, periodic):
    """Centroid in reduced coordinates, circular mean along periodic axes."""
    C1 = np.asarray(C1, float)[:2]; C2 = np.asarray(C2, float)[:2]
    M = np.array([N1 * C1, N2 * C2]).T
    s = np.linalg.solve(M, coords[:2])                        # reduced in [0,1)
    W = wt.sum()
    red = []
    for ax in range(2):
        if periodic[ax]:
            th = 2 * np.pi * s[ax]
            ang = np.arctan2(np.sum(wt * np.sin(th)), np.sum(wt * np.cos(th)))
            red.append((ang / (2 * np.pi)) % 1.0)
        else:
            red.append(np.sum(wt * s[ax]) / W)
    return M @ np.array(red)


def evaluate(m, coords, simp, P, N1, N2, C1, C2, periodic, origin=None,
             weight="site", mmom=None, g=2.0, sigma_max=0.6):
    """Return the oam_traj column set for one configuration m (3,Natom)."""
    psi = m[0] + 1j * m[1]
    gx, gy, Ai, ok, _ = fem_gradient(psi, simp, P)
    w = np.where(ok, 1.0 if weight == "site" else Ai, 0.0)
    norm = np.sum(np.abs(psi) ** 2 * w)
    if origin is None:
        origin = coords[:2].mean(axis=1)
    mmom = np.ones(psi.size) if mmom is None else mmom
    Nm = float(np.sum(mmom / g * (1.0 - m[2])))
    out = dict(N_m=Nm)
    if norm < 1e-14:
        out.update(lambda_L_origin=float("nan"), lambda_L_centroid=float("nan"),
                   R_x=float("nan"), R_y=float("nan"), sigma_psi=float("nan"))
        return out

    def lam(x0):
        # positions relative to x0; for the centroid on a periodic cell use
        # minimum image so the lever arm is measured to the nearest image
        d = coords[:2] - np.asarray(x0)[:, None]
        if x0 is not origin:
            Mc = np.array([N1 * np.asarray(C1, float)[:2], N2 * np.asarray(C2, float)[:2]]).T
            s = np.linalg.solve(Mc, d)
            for ax in range(2):
                if periodic[ax]:
                    s[ax] -= np.round(s[ax])
            d = Mc @ s
        ell = np.imag(np.conj(psi) * (d[0] * gy - d[1] * gx))
        return float(np.sum(ell * w) / norm), d

    lo, _ = lam(origin)
    R = circular_centroid(coords, np.abs(psi) ** 2 * w, N1, N2, C1, C2, periodic)
    lc, dR = lam(R)
    sig = float(np.sqrt(np.sum(np.sum(dR ** 2, 0) * np.abs(psi) ** 2 * w) / norm))
    half = 0.5 * min(np.linalg.norm(N1 * np.asarray(C1)[:2]), np.linalg.norm(N2 * np.asarray(C2)[:2]))
    if sig > sigma_max * half:
        lc = float("nan")
    out.update(lambda_L_origin=lo, lambda_L_centroid=lc, R_x=float(R[0]), R_y=float(R[1]),
               sigma_psi=sig, Lz_tot_hbar=Nm * lc if np.isfinite(lc) else float("nan"),
               dSz_hbar=Nm)
    out["balance"] = out["dSz_hbar"] + out["Lz_tot_hbar"]
    return out


# ------------------------------------------------------------ test fields --
def vortex(coords, center, ell, amp, r0):
    """psi = amp (r/r0)^|l| exp(-r^2/2r0^2) exp(i l phi), m_z = +sqrt(1-|psi|^2)."""
    d = coords[:2] - np.asarray(center)[:, None]
    r = np.hypot(d[0], d[1]); phi = np.arctan2(d[1], d[0])
    g = (r / r0) ** abs(ell) * np.exp(-r ** 2 / (2 * r0 ** 2))
    g = amp * g / g.max()
    psi = g * np.exp(1j * ell * phi)
    return np.stack([psi.real, psi.imag, np.sqrt(np.clip(1 - np.abs(psi) ** 2, 0, None))])


def self_test():
    N = 41; C1 = (1, 0, 0); C2 = (0, 1, 0)
    xy = build_coords(N, N, C1, C2, [(0, 0, 0)])
    simp, P = triangulate(N, N, 1, xy, C1, C2, (True, True))
    A = 0.5 * np.abs((P[:, 1, 0] - P[:, 0, 0]) * (P[:, 2, 1] - P[:, 0, 1]) - (P[:, 1, 1] - P[:, 0, 1]) * (P[:, 2, 0] - P[:, 0, 0]))
    assert abs(A.sum() - N * N) < 1e-9 and A.min() > 0.49, "mesh area"
    c = xy[:2].mean(1)
    for ell in (-2, -1, 1, 2, 3):
        r = evaluate(vortex(xy, c, ell, 0.05, 6.0), xy, simp, P, N, N, C1, C2, (True, True))
        assert abs(r["lambda_L_origin"] - ell) < 0.06 and abs(r["lambda_L_centroid"] - ell) < 0.06, (ell, r)
    # packet straddling the periodic boundary: centroid must land on it
    m = vortex(xy, (0.0, 20.0), 1, 0.05, 4.0)
    m2 = vortex(xy, (41.0, 20.0), 1, 0.05, 4.0)
    m = np.where(np.abs(m[0] + 1j * m[1]) > np.abs(m2[0] + 1j * m2[1]), m, m2)
    r = evaluate(m, xy, simp, P, N, N, C1, C2, (True, True))
    assert min(abs(r["R_x"]), abs(r["R_x"] - 41)) < 0.5, r
    assert abs(r["lambda_L_centroid"] - 1) < 0.06, r
    print("oracle_traj self-test passed")


if __name__ == "__main__":
    self_test()
