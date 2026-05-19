"""
This module provides functions to read and process UppASD simulation data files.
It includes functions to calculate reciprocal lattices, read atomic positions,
moments, and important keywords from UppASD input files.
Additionally, it provides functions to read and process spin Hamiltonian data.

Functions:
- get_reciprocal_lattice(lattice): Calculate the reciprocal lattice from a given lattice.
- read_posfile(posfile): Read atomic positions from a UppASD position file.
- read_momfile(momfile): Read atomic moments from a UppASD moment file.
- read_inpsd(ifile): Read important keywords from a UppASD input file `inpsd.dat`.
- get_uppasd_cell(infile): Convert UppASD geometry data to a cell tuple and get symmetry data.
- read_spin_hamiltonian(hamfile, atoms): Read spin Hamiltonian from UppASD `struct.simid.out` file.
- read_spin_hamiltonian_arr(hamfile, atoms): Read spin Hamiltonian from UppASD `struct.simid.out`
    file and store output as arrays.
"""
import numpy as np
import spglib as spg

############################################################
# Calculate reciprocal lattice
############################################################


def get_reciprocal_lattice(lattice):
    """Calculate reciprocal lattice"""
    avec_1 = lattice[0, :]
    avec_2 = lattice[1, :]
    avec_3 = lattice[2, :]
    volume = np.dot(avec_1, np.cross(avec_2, avec_3))

    kvec_1 = np.cross(avec_2, avec_3) / volume
    kvec_2 = np.cross(avec_3, avec_1) / volume
    kvec_3 = np.cross(avec_1, avec_2) / volume
    inv_lat = 2.0 * np.pi * np.asarray([kvec_1, kvec_2, kvec_3]).reshape(3, 3)

    return inv_lat


############################################################
# Read atom positions from UppASD position file
############################################################
def read_posfile(posfile):
    """Read atom positions from UppASD position file"""
    with open(posfile, "r", encoding="utf-8") as pfile:
        lines = pfile.readlines()
        positions = np.empty([0, 3])
        numbers = []
        for line in lines:
            line_data = line.rstrip("\n").split()
            if len(line_data) > 0:
                positions = np.vstack(
                    (positions, np.asarray(line_data[2:5]).astype(np.float64))
                )
                numbers = np.append(numbers, np.asarray(line_data[1]).astype(np.int32))
    return positions, numbers


############################################################
# Read atomic moments from UppASD moment file
############################################################


def read_momfile(momfile):
    """Read atomic moments from UppASD moment file"""

    moments = np.genfromtxt(momfile)[:, 2]

    return moments


############################################################
# Read important keywords from UppASD inputfile `inpsd.dat`
############################################################
def read_inpsd(ifile):
    """Read important keywords from UppASD inputfile `inpsd.dat`"""
    posfiletype = "C"
    with open(ifile, "r", encoding="utf-8") as infile:
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

    ############################################################
    # Convert positions to direct coordinates if needed
    ############################################################
    if posfiletype == "C":
        invlat = np.linalg.inv(lattice)
        invpos = np.copy(positions)
        for row in range(invpos.shape[0]):
            invpos[row, :] = np.matmul(invlat.T, positions[row, :])
        positions = np.copy(invpos)

    return lattice, positions, numbers, simid, mesh


def get_uppasd_cell(infile):
    """Convert UppASD geometry data to cell tuple"""
    lattice, positions, numbers, _, _ = read_inpsd(infile)

    cell = (lattice, positions, numbers)

    sym_data = spg.get_symmetry_dataset(cell)

    return cell, sym_data


def read_spin_hamiltonian(hamfile, atoms):
    """Read spin Hamiltonian from UppASD struct.simid.out
        - stores output as lists
    """

    hamdata = np.genfromtxt(hamfile)

    nnvec = []
    nndist = []
    nntype = []
    jxc = []

    for atom in atoms:
        locham = hamdata[hamdata[:, 0] == atom]
        i_list = []
        r_list = []
        d_list = []
        xc_list = []

        for coupling in locham:
            i_list.append(np.int32(coupling[3]))
            r_list.append(coupling[4:7])
            d_list.append(np.linalg.norm(coupling[4:7]))
            xc_list.append(coupling[7])

        nnvec.append(r_list)
        nntype.append(i_list)
        nndist.append(d_list)
        jxc.append(xc_list)

    return nnvec, nntype, nndist, jxc


def read_spin_hamiltonian_arr(hamfile, atoms):
    """Read spin Hamiltonian from UppASD struct.simid.out
        - stores output as arrays
    """

    hamdata = np.genfromtxt(hamfile)

    nnvec = []
    nndist = []
    nntype = []
    jxc = []

    for atom in atoms:
        locham = hamdata[hamdata[:, 0] == atom]
        i_arr = np.int32(locham[:, 3])
        r_arr = locham[:, 4:7]
        d_arr = np.linalg.norm(r_arr, axis=1)
        xc_arr = locham[:, 7]

        nnvec.append(r_arr)
        nntype.append(i_arr)
        nndist.append(d_arr)
        jxc.append(xc_arr)

    return nnvec, nntype, nndist, jxc
