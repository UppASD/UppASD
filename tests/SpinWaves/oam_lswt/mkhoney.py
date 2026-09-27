"""Generate UppASD inputs: FM honeycomb, NN J, NNN Haldane-type DMI along z."""
import numpy as np, os, sys
s3 = np.sqrt(3.0)
def write(dirname, C2, J=1.0, D=0.1, hz=0.0, kgrid=(60,60), nphi=128, nr=32, extra=""):
    os.makedirs(dirname, exist_ok=True)
    C1 = np.array([1.0, 0.0, 0.0]); C2 = np.array(C2); C3 = np.array([0,0,1.0])
    A = np.array([0.0, 0.0, 0.0]); B = np.array([0.5, 0.5/s3, 0.0])
    # NN vectors A->B
    nn = [B - A, B - A - np.array([1.0,0,0]), B - A - np.array([0.5, s3/2, 0])]
    # NNN: three vectors at 120 deg (counter-clockwise set)
    v = [np.array([1.0,0,0]), np.array([-0.5, s3/2, 0]), np.array([-0.5,-s3/2,0])]
    with open(f"{dirname}/posfile","w") as f:
        f.write(f"1 1 {A[0]:.10f} {A[1]:.10f} 0.0\n2 2 {B[0]:.10f} {B[1]:.10f} 0.0\n")
    with open(f"{dirname}/momfile","w") as f:
        f.write("1 1 1.0 0.0 0.0 1.0\n2 1 1.0 0.0 0.0 1.0\n")
    with open(f"{dirname}/jfile","w") as f:
        for d in nn:
            f.write(f"1 2 {d[0]:.10f} {d[1]:.10f} 0.0 {J}\n")
            f.write(f"2 1 {-d[0]:.10f} {-d[1]:.10f} 0.0 {J}\n")
    with open(f"{dirname}/dmfile","w") as f:
        for i, sgn in ((1, +1), (2, -1)):
            for d in v:
                f.write(f"{i} {i} {d[0]:.10f} {d[1]:.10f} 0.0 0.0 0.0 {sgn*D}\n")
                f.write(f"{i} {i} {-d[0]:.10f} {-d[1]:.10f} 0.0 0.0 0.0 {-sgn*D}\n")
    with open(f"{dirname}/inpsd.dat","w") as f:
        f.write(f"""simid honey
ncell 24 24 1
BC P P 0
cell {C1[0]:.10f} {C1[1]:.10f} 0.0
     {C2[0]:.10f} {C2[1]:.10f} 0.0
     0.0 0.0 1.0
Sym 0
posfile ./posfile
momfile ./momfile
exchange ./jfile
dm ./dmfile
maptype 1
posfiletype C
initmag 3
ip_mode N
mode S
temp 0.0
Nstep 1
damping 0.1
timestep 1e-16
hfield 0.0 0.0 {hz}
do_avrg N
do_chern Y
do_oam_lswt Y
kgrid {kgrid[0]} {kgrid[1]} 1
oam_nphi {nphi}
oam_nr {nr}
{extra}
""")
if __name__ == "__main__":
    write("runs/hexA", [0.5, s3/2, 0.0])
    write("runs/hexB", [-0.5, s3/2, 0.0])

def write_aniso(dirname, C2, K1=0.05, **kw):
    write(dirname, C2, extra="anisotropy ./kfile", **kw)
    with open(f"{dirname}/kfile","w") as f:
        f.write(f"1 1 {K1} 0.0 0.0 0.0 1.0 0.0\n2 1 {K1} 0.0 0.0 0.0 1.0 0.0\n")
