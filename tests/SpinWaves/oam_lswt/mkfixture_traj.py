"""Write an UppASD run directory holding a prescribed spin texture.

The texture is loaded through Initmag 4 (restart file) and measured at the
first step.  Nstep=1 with a 1e-20 s timestep and zero damping/temperature
keeps the configuration effectively frozen (drift ~1e-8). The restart
iteration column is 0 so rstep=0 and the first measured step is mstep=1, so the measured OAM is the OAM
of the prescribed field and can be compared with oracle_traj.evaluate().

Lattices: 'square' (C1=(1,0), C2=(0,1)) and 'hex' (C1=(1,0), C2=(1/2,sqrt3/2)),
one atom per cell, nearest-neighbour FM exchange.  ``layers`` makes a stack of
independent copies along C3; exchange remains in-plane.  Exchange only matters
for the (negligible) evolution during the single step.
"""
import os
import numpy as np
import oracle_traj as O

S3 = np.sqrt(3.0)
LATTICES = {
    "square": dict(C1=(1.0, 0.0, 0.0), C2=(0.0, 1.0, 0.0),
                   nn=[(1, 0), (-1, 0), (0, 1), (0, -1)]),
    "hex": dict(C1=(1.0, 0.0, 0.0), C2=(0.5, S3 / 2, 0.0),
                nn=[(1, 0), (-1, 0), (0.5, S3 / 2), (-0.5, -S3 / 2), (-0.5, S3 / 2), (0.5, -S3 / 2)]),
}


def write_case(dirname, *, lattice="square", N=41, bc=("P", "P"), shift=(0.0, 0.0),
               field=None, extra_keys=None, simid="oamtest", layers=1,
               ensemble_fields=None):
    """field: callable coords(3,Natom) -> m(3,Natom).  Returns coords."""
    if layers < 1:
        raise ValueError("layers must be positive")
    os.makedirs(dirname, exist_ok=True)
    L = LATTICES[lattice]
    # UppASD folds basis positions into the cell (reduced coordinates in [0,1));
    # do the same so the texture is generated on the coordinates UppASD uses.
    Mc = np.array([L["C1"][:2], L["C2"][:2]]).T
    red = np.linalg.solve(Mc, np.asarray(shift, float)) % 1.0
    shift = tuple(Mc @ red)
    basis = [(shift[0], shift[1], 0.0)]
    xy = O.build_coords(N, N, L["C1"], L["C2"], basis)
    coords = np.concatenate([xy + np.array([0.0, 0.0, float(z)])[:, None]
                             for z in range(layers)], axis=1)
    nsite = N * N * layers
    if ensemble_fields is None:
        ensemble_fields = [field(coords) if field is not None else
                           np.vstack([np.zeros((2, nsite)), np.ones(nsite)])]
    m = np.stack([(value(coords) if callable(value) else value) /
                  np.linalg.norm(value(coords) if callable(value) else value, axis=0)
                  for value in ensemble_fields], axis=2)
    nens = m.shape[2]
    with open(f"{dirname}/posfile", "w") as f:
        f.write(f"1 1 {shift[0]:.10f} {shift[1]:.10f} 0.0\n")
    with open(f"{dirname}/momfile", "w") as f:
        f.write("1 1 1.0 0.0 0.0 1.0\n")
    with open(f"{dirname}/jfile", "w") as f:
        for d in L["nn"]:
            f.write(f"1 1 {d[0]:.10f} {d[1]:.10f} 0.0 1.0\n")
    with open(f"{dirname}/restart.in", "w") as f:
        f.write("#" * 80 + "\n# File type: R\n# Simulation type: S\n"
                f"# Number of atoms: {nsite:9d}\n# Number of ensembles: {nens:9d}\n" + "#" * 80 + "\n")
        f.write("  # iter     ens   iatom           |Mom|             M_x             M_y             M_z\n")
        for k in range(nens):
            for i in range(nsite):
                f.write(f"{0:8d}{k+1:8d}{i+1:8d}  {1.0:16.8E}"
                        f"{m[0,i,k]:24.16E}{m[1,i,k]:24.16E}{m[2,i,k]:24.16E}\n")
    keys = {"do_oam_traj": "Y", "oam_step": "1"}
    keys.update(extra_keys or {})
    c1, c2 = L["C1"], L["C2"]
    with open(f"{dirname}/inpsd.dat", "w") as f:
        f.write(f"""simid {simid}
ncell {N} {N} {layers}
Mensemble {nens}
BC {bc[0]} {bc[1]} 0
cell {c1[0]:.10f} {c1[1]:.10f} 0.0
     {c2[0]:.10f} {c2[1]:.10f} 0.0
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
""")
        for k, v in keys.items():
            f.write(f"{k} {v}\n")
    return coords, m
