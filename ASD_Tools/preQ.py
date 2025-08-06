# preQ.py
# Automagically creates a standardized path in reciprocal space from your `inpsd.dat`
# Also constructs the reduced q-mesh with proper symmetries for all your magnon DOS needs..
# (C) Anders Bergman, Uppsala University 2019

# Preprocessing script for UppASD simulations
# Reads the inpsd.dat input file to extract the lattice and basis for the system.
# That structural information is used to obtain the space-group, BZ symmetry points and
# suitable q-vector paths/grids for S(q,w) simulations using UppASD. The output q-vectors
# are given in direct coordinates.
#
# Dependencies: python>3.0, numpy, spglib, seekpath, tabulate. All available through `pip`.


# All important functionalities are provided by `spglib` and `seekpath`.
# References:
# 𝚂𝚙𝚐𝚕𝚒𝚋: a software library for crystal symmetry search”, Atsushi Togo and Isao Tanaka, https://arxiv.org/abs/1808.01590 (written at version 1.10.4)
# SeeKpath: Y. Hinuma, G. Pizzi, Y. Kumagai, F. Oba, I. Tanaka, Band structure diagram paths based on crystallography, Comp. Mat. Sci. 128, 140 (2017)
# UppASD: a simulation tool for atomistic spin dynamics and Monte Carlo simulations of Heisenberg spin systems. https://github.com/UppASD

import numpy as np
from tabulate import tabulate
from asd_io import (
    get_reciprocal_lattice,
    get_spacegroup,
    get_symmetry_points,
    get_kpath,
    read_inpsd,
)

# Read input
ifile = "inpsd.dat"
inpsd_data = read_inpsd(ifile)
lattice = inpsd_data["lattice"]
positions = inpsd_data["positions"]
numbers = inpsd_data["numbers"]
simid = inpsd_data["simid"]
mesh = inpsd_data["mesh"]
posfiletype = inpsd_data["posfiletype"]

# Convert positions to direct coordinates if needed
if posfiletype == "C":
    invlat = np.linalg.inv(lattice)
    invpos = np.copy(positions)
    for row in range(invpos.shape[0]):
        invpos[row, :] = np.matmul(invlat.T, positions[row, :])
    positions = np.copy(invpos)

# Print structure info
print("\nStructural data for UppASD simulation ", simid)
cell = (lattice, positions, numbers)
spacegroup = get_spacegroup(cell)
print("\nLattice:")
print(tabulate(lattice, floatfmt=".4f"))
print("\nAtomic positions:")
print(tabulate(positions, floatfmt=" .4f"))
print("\nSpacegroup:\n", spacegroup)

# Get symmetry points and primitive lattice
BZ, sympoints = get_symmetry_points(cell)
print("\nPrimitive lattice:")
print(tabulate(BZ, floatfmt=".4f"))

dictlist = []
for key, value in sympoints.items():
    str1 = " ".join("{: 4.4f}".format(e) for e in value)
    temp = [key, str1]
    dictlist.append(temp)

print("\nSymmetry points:")
print(tabulate(dictlist, floatfmt=" .4f"))

# Get k-path and print path string
import seekpath as spth
kpath_obj = spth.get_path(cell)
klist = list(sum(kpath_obj['path'], ()))
kpath = [ ' -> '.join(x) for x in zip(klist[0::2], klist[1::2]) ]
kstr = ', '.join('{}'.format(e) for e in kpath)
print("\nK-path written to 'qfile.kpath':")
print(kstr, '\n')

# Interpolate k-points and associate labels
mypath = spth.getpaths.get_explicit_k_path(cell, reference_distance=10)
mypath = mypath['explicit_kpoints_rel']
lpath = []
xpath = mypath[0,:]
for row in range(mypath.shape[0]-1):
    nwaves = (mypath[row+1,:] - mypath[row,:]) * mesh
    nsteps = np.int32(np.max(np.abs(nwaves)))
    lpath.append([k for k,v in kpath_obj['point_coords'].items() if (v == mypath[row]).all()][0])
    for iint in range(nsteps):
        ipart = (iint+1)/(nsteps)
        xpath = np.append(xpath, mypath[row,:]*(1.0-ipart) + mypath[row+1,:]*ipart, axis=0)
        if(iint < nsteps-1):
            lpath.append(' ')
lpath.append([k for k,v in kpath_obj['point_coords'].items() if (v == mypath[mypath.shape[0]-1]).all()][0])
np.savetxt('qfile.klabel', lpath, fmt='%s')

xpath = xpath.reshape(int(xpath.shape[0]/3),3)
# Save q-point mesh to file
nq = xpath.shape[0]
with open('qfile.kpath','w') as qf:
    print("         ", nq, file=qf)
    for i, gp in enumerate(xpath):
        print("%s     %s" % (' '.join(format(f, '10.8f') for f in gp), lpath[i]), file=qf)

# Save 2D q-point mesh to file
xpath2d = xpath[xpath[:,2]==0]
nq = xpath2d.shape[0]
with open('qfile.kpath2d','w') as qf:
    print("         ", nq, file=qf)
    for i, gp in enumerate(xpath):
        if gp[2]==0:
            print("%s     %s" % (' '.join(format(f, '10.8f') for f in gp), lpath[i]), file=qf)

# Reduced k-mesh (unchanged, can use modular routine if available)
import spglib as spg
mapping, grid = spg.get_ir_reciprocal_mesh(mesh, cell, is_shift=[0, 0, 0])
print(
    "\nMaximum k-point mesh size from `inpsd.dat`: \n",
    mesh[0],
    "x",
    mesh[1],
    "x",
    mesh[2],
)
irk_idx, mult = np.unique(mapping, return_counts=True)
print(
    "\nNumber of reduced k-points: ", irk_idx.size, " Written to file `qfile.reduced`."
)
with open("qfile.reduced", "w") as qf:
    print("         ", irk_idx.size, file=qf)
    for i, (ir_gp_id, gp) in enumerate(zip(irk_idx, grid[irk_idx])):
        print(
            "%s     %10.4f"
            % (" ".join(format(f, "10.8f") for f in gp / mesh), mult[i]),
            file=qf,
        )

# Done!