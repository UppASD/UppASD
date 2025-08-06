# Structure and symmetry analysis routines
import os
from collections import defaultdict
from string import Template

import numpy as np
import seekpath as spth
import spglib as spg
from pymatgen.core import Structure
from pymatgen.io.cif import CifParser
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer


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


def read_posfile(posfile):
    """Read atomic positions and numbers from UppASD position file."""
    with open(posfile, "r", encoding="utf-8") as pfile:
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
    with open(ifile, "r", encoding="utf-8") as infile:
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
                    except ValueError:
                        data["sc_step"] = line_data[1]
                if key == "sc_nstep":
                    try:
                        data["sc_nstep"] = int(line_data[1])
                    except ValueError:
                        data["sc_nstep"] = line_data[1]
                if key == "qfile":
                    data["qfile"] = line_data[1]
        return data

def extract_chemical_name_common(cif_file):
    """
    Extracts the _chemical_name_common value from a CIF file without using CifParser.

    Args:
        cif_file (str): Path to the CIF file.

    Returns:
        str: The value of _chemical_name_common, or None if not found.
    """
    with open(cif_file, "r", encoding="utf-8") as f:
        for line in f:
            if line.strip().startswith("_chemical_name_common"):
                parts = line.strip().split(None, 1)
                if len(parts) == 2:
                    value = parts[1].strip().strip("'\"")
                    return value
    return None

def extract_structure_from_cif(cif_file, primitive=False):
    """
    Reads a CIF file, extracts lattice parameters, and prepares a spglib cell tuple.

    Args:
        cif_file (str): Path to the CIF file.

    Returns:
        tuple: (lattice matrix, fractional positions, atomic numbers) for spglib.
    """
    parser = CifParser(cif_file)
    structure = parser.parse_structures(primitive=primitive)[0]
    # Extract chemical name using public API or fallback
    chemical_name = extract_chemical_name_common(cif_file)
    if chemical_name is not None:
        chemical_name = str(chemical_name)
    else:
        chemical_name = "unknown_"
    # Extract a_lat from structure.lattice.a
    a_lat = getattr(structure.lattice, 'a', None)

    # Extract lattice parameters
    lattice = structure.lattice
    a, b, c = lattice.a, lattice.b, lattice.c
    alpha, beta, gamma = lattice.alpha, lattice.beta, lattice.gamma
    print("Lattice parameters:")
    print(f"a = {a}, b = {b}, c = {c}")
    print(f"alpha = {alpha}, beta = {beta}, gamma = {gamma}")

    # Prepare spglib cell: (lattice, positions, numbers)
    cell = (
        lattice.matrix,
        [site.frac_coords for site in structure.sites],
        [site.specie.number for site in structure.sites],
    )
    print("\nSpglib cell:")
    print("Lattice matrix:")
    for row in cell[0]:
        print(row)
    print("Fractional positions:")
    for pos in cell[1]:
        print(pos)
    print("Atomic numbers:")
    print(cell[2])
    return cell, chemical_name, a_lat


def write_inpsd_dat(
    filename,
    simid,
    ncell,
    bc,
    cell_matrix,
    alat,
    posfile,
    momfile,
    exchangefile,
    qfile,
    extra_params=None,
):

    template_path = os.path.join(os.path.dirname(__file__), "inpsd_template.txt")
    with open(template_path, encoding="utf-8") as tf:
        template_str = tf.read()
    # Format cell_matrix as required for the template
    cell_matrix_str = "".join(f"     {row[0]:.3f}    {row[1]:.3f}    {row[2]:.3f}\n" for row in cell_matrix)
    # Format alat
    alat_str = f"{alat:.6e}"
    # Format extra_params
    extra_str = ""
    if extra_params:
        for k, v in extra_params.items():
            extra_str += f"{k} {v}\n"
    # Prepare substitution dictionary
    subs = {
        "simid": simid,
        "ncell_x": ncell[0],
        "ncell_y": ncell[1],
        "ncell_z": ncell[2],
        "bc_x": bc[0],
        "bc_y": bc[1],
        "bc_z": bc[2],
        "cell_matrix": cell_matrix_str,
        "alat": alat_str,
        "posfile": posfile,
        "momfile": momfile,
        "exchangefile": exchangefile,
        "qfile": qfile,
        "extra": extra_str,
    }
    content = Template(template_str).substitute(subs)
    with open(filename, "w", encoding="utf-8") as f:
        f.write(content)


def write_posfile(filename, sites):
    """
    Writes the posfile.

    Args:
        filename (str): Output file path.
        sites (list): List of tuples (index, atomic_number, x, y, z).
    """
    with open(filename, "w", encoding="utf-8") as f:
        print(sites)
        for idx, _, x, y, z in sites:
            # Print: idx (site number), idx (again), type_number, x, y, z
            f.write(f"{idx}   {idx}   {x:.6f}  {y:.6f}  {z:.6f}\n")


def write_momfile(filename, moments):
    """
    Writes the momfile.

    Args:
        filename (str): Output file path.
        moments (list): List of tuples (index, atomic_number, mx, my, mz).
    """
    print("Writing moments to file:", moments)
    with open(filename, "w") as f:
        for idx, _, m_mag, mx, my, mz in moments:
            f.write(f"{idx}  1   {m_mag}  {mx:.7f}     {my:.6f}   {mz:.6f}\n")


def find_real_space_neighbors(cell, cutoff):
    """
    Finds real space neighbors within a cutoff distance for each atom in the cell using pymatgen.

    Args:
        cell (tuple): (lattice, positions, numbers) as used by spglib.
        cutoff (float): Cutoff distance in Angstrom.

    Returns:
        dict: {index: [(neighbor_index, distance), ...], ...}
    """
    lattice, positions, numbers = cell
    # Build pymatgen Structure
    structure = Structure(lattice, numbers, positions)
    neighbors = {}
    for i, site in enumerate(structure.sites):
        nbrs = structure.get_neighbors(site, cutoff)
        print(f"Atom {i}:", nbrs)

        neighbors[i] = [(n.index, n.nn_distance) for n in nbrs]
    return neighbors

def get_symmetry_reduced_neighbor_shells(cell, cutoff=5.0, decimal=3):
    """
    Returns symmetry-reduced neighbor shells for a structure using pymatgen.

    Args:
        cell (tuple): (lattice, positions, numbers) as used by spglib.
        cutoff (float): Cutoff distance in Angstrom.
        decimal (int): Decimal places for grouping shells by distance.

    Returns:
        dict: {distance: [(i, j)], ...} where (i, j) are symmetry-unique pairs.
    """
    lattice, positions, numbers = cell
    structure = Structure(lattice, numbers, positions)
    sga = SpacegroupAnalyzer(structure)
    symm_struct = sga.get_symmetrized_structure()
    shells = defaultdict(list)
    # Only consider one representative from each symmetry-equivalent site
    for eq_indices in symm_struct.equivalent_indices:
        i = eq_indices[0]
        site = structure.sites[i]
        nbrs = structure.get_neighbors(site, cutoff)
        for n in nbrs:
            # Only keep unique pairs (i, j) with i <= j to avoid double counting
            j = n.index
            if i <= j:
                dist = round(n.nn_distance, decimal)
                shells[dist].append((i, j))
    return dict(shells)
