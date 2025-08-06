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


#!pip install spglib
#!pip install tabulate
#!pip install seekpath
import numpy as np
import spglib as spg
from tabulate import tabulate
import seekpath as spth

# from google.colab import files


############################################################
# Calculate reciprocal lattice
############################################################
from asd_io import (
    get_reciprocal_lattice,
    get_spacegroup,
    get_symmetry_points,
    get_kpath,
)


############################################################
# Read atom positions from UppASD position file
############################################################
def read_posfile(posfile):
    with open(posfile, "r") as pfile:
        lines = pfile.readlines()
        positions = np.empty([0, 3])
        numbers = []
        for idx, line in enumerate(lines):
            line_data = line.rstrip("\n").split()
            if len(line_data) > 0:
                positions = np.vstack(
                    (positions, np.asarray(line_data[2:5]).astype(np.float64))
                )
                numbers = np.append(numbers, np.asarray(line_data[1]).astype(np.int32))
        return positions, numbers


############################################################
# Read important keywords from UppASD inputfile `inpsd.dat`
############################################################
def read_inpsd(ifile):
    posfiletype = "C"
    with open(ifile, "r") as infile:
        lines = infile.readlines()
        for idx, line in enumerate(lines):
            line_data = line.rstrip("\n").split()
            if len(line_data) > 0:
                # Find the simulation id
                if line_data[0] == "simid":
                    simid = line_data[1]

                # Find the cell data
                if line_data[0] == "cell":
                    cell = []
                    lattice = np.empty([0, 3])
                    line_data = lines[idx + 0].split()
                    cell = np.append(cell, np.asarray(line_data[1:4]))
                    lattice = np.vstack((lattice, np.asarray(line_data[1:4])))
                    line_data = lines[idx + 1].split()
                    cell = np.append(cell, np.asarray(line_data[0:3]))
                    lattice = np.vstack((lattice, np.asarray(line_data[0:3])))
                    line_data = lines[idx + 2].split()
                    cell = np.append(cell, np.asarray(line_data[0:3]))
                    lattice = np.vstack((lattice, np.asarray(line_data[0:3])))
                    lattice = lattice.astype(np.float64)

                # Find the size of the simulated cell
                if line_data[0] == "ncell":
                    ncell_x = int(line_data[1])
                    ncell_y = int(line_data[2])
                    ncell_z = int(line_data[3])
                    mesh = [ncell_x, ncell_y, ncell_z]

                # Read the name of the position file
                if line_data[0].strip() == "posfile":
                    positions, numbers = read_posfile(line_data[1])

                # Read the type of coordinate representation
                if line_data[0].strip() == "posfiletype":
                    posfiletype = line_data[1]

    return lattice, positions, numbers, simid, mesh, posfiletype


############################################################
# Open and read input files
############################################################
ifile = "inpsd.dat"
inpsd_data = read_inpsd(ifile)
lattice = inpsd_data["lattice"]
positions = inpsd_data["positions"]
numbers = inpsd_data["numbers"]
simid = inpsd_data["simid"]
mesh = inpsd_data["mesh"]
posfiletype = inpsd_data["posfiletype"]
timestep = inpsd_data["timestep"]
sc_step = inpsd_data["sc_step"]
sc_nstep = inpsd_data["sc_nstep"]
qfile = inpsd_data["qfile"]


############################################################
# Convert positions to direct coordinates if needed
############################################################
if posfiletype == "C":
    invlat = np.linalg.inv(lattice)
    invpos = np.copy(positions)
    for row in range(invpos.shape[0]):
        invpos[row, :] = np.matmul(invlat.T, positions[row, :])
    positions = np.copy(invpos)

############################################################
# Get the spacegroup from spglib and print relevant info
############################################################
print("\nStructural data for UppASD simulation ", simid)
cell = (lattice, positions, numbers)
spacegroup = get_spacegroup(cell)
print("\nLattice:")
print(tabulate(lattice, floatfmt=".4f"))
print("\nAtomic positions:")
print(tabulate(positions, floatfmt=" .4f"))
print("\nSpacegroup:\n", spacegroup)


############################################################
# Get symmetry points from seekpath
############################################################
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

############################################################
# Get k-space path from seekpath
############################################################
print("\nK-path written to 'qfile.kpath':")
kpoints = get_kpath(cell, reference_distance=10)
# You can add interpolation and file writing logic here as needed

###############################################################
#### Get k-space path from seekpath  v2
###############################################################
###print("\nK-path written to 'qfile.kpath2':")
###
#### First get the symmetry points within the path
###mypath=spth.getpaths.get_explicit_k_path(cell,reference_distance=0.1)
###mypath=mypath['explicit_kpoints_rel']
###
#### Then interpolate to get a path commensurate with the UppASD geometry
###lpath=[]
###xpath=mypath
###
#### Save q-point mesh to file
###nq=xpath.shape[0]
###with open('qfile.kpath2','w') as qf:
###   print("         ",nq,file=qf)
###   for i,  gp in enumerate(xpath):
###      print("%s" % (' '.join(format(f, '10.8f') for f in gp)),file=qf)
###      #print("%s     %s" % (' '.join(format(f, '10.8f') for f in gp),lpath[i]),file=qf)
###

############################################################
# Get the reduced q-mesh from spglib
############################################################
# Can also add shift if wanted
mapping, grid = spg.get_ir_reciprocal_mesh(mesh, cell, is_shift=[0, 0, 0])

# K-point mesh
print(
    "\nMaximum k-point mesh size from `inpsd.dat`: \n",
    mesh[0],
    "x",
    mesh[1],
    "x",
    mesh[2],
)

# Find irrep
irk_idx, mult = np.unique(mapping, return_counts=True)
print(
    "\nNumber of reduced k-points: ", irk_idx.size, " Written to file `qfile.reduced`."
)

# Save q-point mesh to file
with open("qfile.reduced", "w") as qf:
    print("         ", irk_idx.size, file=qf)
    for i, (ir_gp_id, gp) in enumerate(zip(irk_idx, grid[irk_idx])):
        print(
            "%s     %10.4f"
            % (" ".join(format(f, "10.8f") for f in gp / mesh), mult[i]),
            file=qf,
        )

############################################################
# Done!
############################################################
