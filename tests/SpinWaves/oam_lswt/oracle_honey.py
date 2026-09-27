"""Independent LSWT oracle for FM honeycomb + NNN (Haldane) DMI, spins along z.

Particle-conserving 2x2 magnon Hamiltonian (site-position Fourier convention):
  h(k) = d0(k) I + d(k).sigma,
  d_x - i d_y = -J S sum_delta exp(i k.delta)      (NN, A->B)
  d_z         = 2 D S nu sum_i sin(k.v_i)           (NNN, Haldane)
Eigenvectors depend only on d-hat, so F(k) depends only on D/J (and sign).
F_n(k) = ring average of O_n = -1/2 Im[u^dag d_phi u] in a gauge regular at k=0
       = (1/2) * Berry phase of the ring / (2 pi)   (Stokes)
The Berry phase is computed from a gauge-invariant Wilson loop and unwrapped
continuously from k=0, so no gauge fixing is involved.
"""
import numpy as np
s3 = np.sqrt(3.0)
A = np.array([0.0, 0.0]); B = np.array([0.5, 0.5/s3])
NN = [B-A, B-A-np.array([1.0,0]), B-A-np.array([0.5,s3/2])]
V  = [np.array([1.0,0]), np.array([-0.5,s3/2]), np.array([-0.5,-s3/2])]

def dvec(k, J=1.0, D=0.1, S=1.0, nu=+1):
    g = sum(np.exp(1j*k@d) for d in NN)
    dxy = -J*S*g
    dz = 2*D*S*nu*sum(np.sin(k@v) for v in V)
    return np.array([dxy.real, -dxy.imag, dz])

def evecs(k, **kw):
    d = dvec(k, **kw)
    h = np.array([[d[2], d[0]-1j*d[1]],[d[0]+1j*d[1], -d[2]]])
    w, u = np.linalg.eigh(h)
    return w[::-1], u[:, ::-1]          # band 1 = upper (optical), band 2 = lower

def ring_berry_phase(pts, band, **kw):
    us = [evecs(p, **kw)[1][:, band] for p in pts]
    prod = 1.0+0j
    for j in range(len(us)):
        prod *= np.vdot(us[j], us[(j+1) % len(us)])
    return -np.angle(prod)              # Berry phase gamma = -Im log prod

def F_of_k(kvals, band, mapping=None, nphi=512, **kw):
    """F(k) with branch continuity from k=0.  mapping: optional k->k_eff."""
    out, prev = [], 0.0
    for k in kvals:
        phis = 2*np.pi*np.arange(nphi)/nphi
        pts = [np.array([k*np.cos(p), k*np.sin(p)]) for p in phis]
        if mapping is not None:
            pts = [mapping(p) for p in pts]
        g = ring_berry_phase(pts, band, **kw) if k > 0 else 0.0
        g += 2*np.pi*np.round((prev-g)/(2*np.pi))
        prev = g
        out.append(g/(4*np.pi))         # F = gamma/(4 pi)
    return np.array(out)

def uppasd_polar_map(C1, C2):
    """What calculate_fishman_f_average actually hands to the Hamiltonian:
    q = k/(2pi) in Cartesian components, consumed as k_eff = 2pi q."""
    return lambda k: np.asarray(k)[:2]
