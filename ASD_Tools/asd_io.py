# Structure and symmetry analysis routines
import spglib as spg
import seekpath as spth


def get_reciprocal_lattice(lattice):
    """Calculate reciprocal lattice from direct lattice vectors."""
    a1 = lattice[0, :]
    a2 = lattice[1, :]
    a3 = lattice[2, :]
    vol = np.dot(a1, np.cross(a2, a3))
    k1 = np.cross(a2, a3) / vol
    k2 = np.cross(a3, a1) / vol
    k3 = np.cross(a1, a2) / vol
    K = 2.0 * np.pi * np.asarray([k1, k2, k3]).reshape(3, 3)
    return K


def get_spacegroup(cell, symprec=1e-5):
    """Get spacegroup symbol and number from spglib."""
    return spg.get_spacegroup(cell, symprec=symprec)


def get_symmetry_points(cell):
    """Get symmetry points and reciprocal primitive lattice from seekpath."""
    kpath_obj = spth.get_path(cell)
    BZ = kpath_obj["reciprocal_primitive_lattice"]
    sympoints = kpath_obj["point_coords"]
    return BZ, sympoints


def get_kpath(cell, mesh=None, reference_distance=None):
    """Get explicit k-path from seekpath, optionally interpolated to mesh or reference_distance."""
    if reference_distance is not None:
        mypath = spth.getpaths.get_explicit_k_path(
            cell, reference_distance=reference_distance
        )
    else:
        mypath = spth.getpaths.get_explicit_k_path(cell)
    kpoints = mypath["explicit_kpoints_rel"]
    # Optionally interpolate to mesh
    if mesh is not None:
        # Interpolation logic can be added here if needed
        pass
    return kpoints


# asd_io.py
import numpy as np


def read_posfile(posfile):
    """Read atomic positions and numbers from UppASD position file."""
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


def read_inpsd(ifile):
    """Read important keywords from UppASD inputfile `inpsd.dat`. Returns lattice, positions, numbers, simid, mesh, posfiletype."""
    with open(ifile, "r") as infile:
        lines = infile.readlines()
        data = {
            "lattice": None,
            "positions": None,
            "numbers": None,
            "simid": None,
            "mesh": None,
            "posfiletype": "C",
            "timestep": None,
            "sc_step": None,
            "sc_nstep": None,
            "qfile": None,
        }
        for idx, line in enumerate(lines):
            line_data = line.rstrip("\n").split()
            if len(line_data) > 0:
                key = line_data[0].strip()
                if key == "simid":
                    data["simid"] = line_data[1]
                if key == "cell":
                    lattice = np.empty([0, 3])
                    line_data = lines[idx + 0].split()
                    lattice = np.vstack(
                        (lattice, np.asarray(line_data[1:4]).astype(np.float64))
                    )
                    line_data = lines[idx + 1].split()
                    lattice = np.vstack(
                        (lattice, np.asarray(line_data[0:3]).astype(np.float64))
                    )
                    line_data = lines[idx + 2].split()
                    lattice = np.vstack(
                        (lattice, np.asarray(line_data[0:3]).astype(np.float64))
                    )
                    data["lattice"] = lattice
                if key == "ncell":
                    if len(line_data) == 4:
                        ncell_x = int(line_data[1])
                        ncell_y = int(line_data[2])
                        ncell_z = int(line_data[3])
                    else:
                        ncell_x = int(line_data[1])
                        ncell_y = int(line_data[1])
                        ncell_z = int(line_data[1])
                    data["mesh"] = [ncell_x, ncell_y, ncell_z]
                if key == "posfile":
                    positions, numbers = read_posfile(line_data[1])
                    data["positions"] = positions
                    data["numbers"] = numbers
                if key == "posfiletype":
                    data["posfiletype"] = line_data[1]
                if key == "timestep":
                    try:
                        data["timestep"] = float(line_data[1])
                    except Exception:
                        data["timestep"] = line_data[1]
                if key == "sc_step":
                    try:
                        data["sc_step"] = int(line_data[1])
                    except Exception:
                        data["sc_step"] = line_data[1]
                if key == "sc_nstep":
                    try:
                        data["sc_nstep"] = int(line_data[1])
                    except Exception:
                        data["sc_nstep"] = line_data[1]
                if key == "qfile":
                    data["qfile"] = line_data[1]
        return data
