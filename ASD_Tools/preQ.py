#!/usr/bin/env python
# coding: utf-8

"""
UppASD Pre-Processing Module

This module provides the PreProcessor class for preparing UppASD simulations
by analyzing crystal structure, generating k-paths, and creating q-point meshes.

Automagically creates a standardized path in reciprocal space from your `inpsd.dat`
Also constructs the reduced q-mesh with proper symmetries for all your magnon DOS needs.

(C) Anders Bergman, Uppsala University 2019

Dependencies: python>3.0, numpy, spglib, seekpath, tabulate. All available through `pip`.

References:
- Spglib: a software library for crystal symmetry search, Atsushi Togo and Isao Tanaka
- SeeKpath: Y. Hinuma, G. Pizzi, Y. Kumagai, F. Oba, I. Tanaka, Band structure diagram paths based on crystallography, Comp. Mat. Sci. 128, 140 (2017)
- UppASD: a simulation tool for atomistic spin dynamics and Monte Carlo simulations of Heisenberg spin systems
"""

import numpy as np
import spglib as spg
import seekpath as spth
from tabulate import tabulate
from asd_io import (
    get_spacegroup,
    get_symmetry_points,
    read_inpsd,
)


class PreProcessor:
    """Pre-process UppASD simulations: analyze structure and generate k-paths."""

    def __init__(self, input_file="inpsd.dat"):
        """Initialize the PreProcessor with input file."""
        self.input_file = input_file

        # Input data dictionary
        self.inputs = {}

        # Structure and symmetry data
        self.structure = {
            "lattice": None,
            "positions": None,
            "numbers": None,
            "cell": None,
            "spacegroup": None,
            "BZ": None,
            "sympoints": None,
        }

        # K-path data
        self.kpath_data = {
            "kpath_obj": None,
            "xpath": None,
            "lpath": None,
            "xpath2d": None,
            "lpath2d": None,
        }

        # Mesh data
        self.mesh_data = {
            "mapping": None,
            "grid": None,
            "irk_idx": None,
            "mult": None,
        }

    def read_inputs(self):
        """Read input parameters from inpsd.dat."""
        inpsd_data = read_inpsd(self.input_file)
        self.inputs = {
            "lattice": inpsd_data["lattice"],
            "positions": inpsd_data["positions"],
            "numbers": inpsd_data["numbers"],
            "simid": inpsd_data["simid"],
            "mesh": inpsd_data["mesh"],
            "posfiletype": inpsd_data["posfiletype"],
        }

        # Store structure data
        self.structure["lattice"] = self.inputs["lattice"]
        self.structure["positions"] = self.inputs["positions"]
        self.structure["numbers"] = self.inputs["numbers"]

    def process_positions(self):
        """Convert positions to direct coordinates if needed."""
        if self.inputs["posfiletype"] == "C":
            lattice = self.structure["lattice"]
            positions = self.structure["positions"]
            invlat = np.linalg.inv(lattice)
            invpos = np.copy(positions)
            for row in range(invpos.shape[0]):
                invpos[row, :] = np.matmul(invlat.T, positions[row, :])
            self.structure["positions"] = np.copy(invpos)

    def analyze_structure(self):
        """Analyze crystal structure and determine space group."""
        lattice = self.structure["lattice"]
        positions = self.structure["positions"]
        numbers = self.structure["numbers"]

        self.structure["cell"] = (lattice, positions, numbers)
        self.structure["spacegroup"] = get_spacegroup(self.structure["cell"])
        self.structure["BZ"], self.structure["sympoints"] = get_symmetry_points(
            self.structure["cell"]
        )

    def generate_kpath(self):
        """Generate k-path using seekpath."""
        cell = self.structure["cell"]
        mesh = self.inputs["mesh"]

        # Get basic k-path
        self.kpath_data["kpath_obj"] = spth.get_path(cell)

        # Use get_explicit_k_path_orig_cell with reference_distance=2*pi/(N+1), N=max(mesh)
        N = max(mesh)
        ref_dist = 2.0 * np.pi / (N + 1)
        mypath_obj = spth.getpaths.get_explicit_k_path_orig_cell(
            cell, reference_distance=ref_dist
        )
        mypath = mypath_obj["explicit_kpoints_rel"]
        mylabels = mypath_obj["explicit_kpoints_labels"]

        xpath = []
        lpath = []
        for vec, label in zip(mypath, mylabels):
            xpath.append(np.array(np.float32(vec)))
            lpath.append([float(vec[0]), float(vec[1]), float(vec[2]), str(label)])

        self.kpath_data["xpath"] = np.array(xpath)
        self.kpath_data["lpath"] = lpath

        # Generate 2D k-path (z=0 plane)
        xpath = self.kpath_data["xpath"]
        lpath = self.kpath_data["lpath"]
        self.kpath_data["xpath2d"] = xpath[xpath[:, 2] == 0]
        self.kpath_data["lpath2d"] = [
            lpath[i] for i in range(len(xpath)) if xpath[i, 2] == 0
        ]

    def generate_reduced_mesh(self):
        """Generate reduced k-point mesh using spglib."""
        cell = self.structure["cell"]
        mesh = self.inputs["mesh"]

        mapping, grid = spg.get_ir_reciprocal_mesh(mesh, cell, is_shift=[0, 0, 0])
        irk_idx, mult = np.unique(mapping, return_counts=True)

        self.mesh_data["mapping"] = mapping
        self.mesh_data["grid"] = grid
        self.mesh_data["irk_idx"] = irk_idx
        self.mesh_data["mult"] = mult

    def print_structure_info(self):
        """Print structural information."""
        simid = self.inputs["simid"]
        lattice = self.structure["lattice"]
        positions = self.structure["positions"]
        spacegroup = self.structure["spacegroup"]
        BZ = self.structure["BZ"]
        sympoints = self.structure["sympoints"]

        print(f"\nStructural data for UppASD simulation {simid}")
        print("\nLattice:")
        print(tabulate(lattice, floatfmt=".4f"))
        print("\nAtomic positions:")
        print(tabulate(positions, floatfmt=" .4f"))
        print(f"\nSpacegroup:\n{spacegroup}")
        print("\nPrimitive lattice:")
        print(tabulate(BZ, floatfmt=".4f"))

        dictlist = []
        for key, value in sympoints.items():
            str1 = " ".join("{: 4.4f}".format(e) for e in value)
            temp = [key, str1]
            dictlist.append(temp)

        print("\nSymmetry points:")
        print(tabulate(dictlist, floatfmt=" .4f"))

    def print_kpath_info(self):
        """Print k-path information."""
        kpath_obj = self.kpath_data["kpath_obj"]

        klist = list(sum(kpath_obj["path"], ()))
        kpath = [" -> ".join(x) for x in zip(klist[0::2], klist[1::2])]
        kstr = ", ".join("{}".format(e) for e in kpath)
        print("\nK-path written to 'qfile.kpath':")
        print(kstr, "\n")

    def print_mesh_info(self):
        """Print mesh information."""
        mesh = self.inputs["mesh"]
        irk_idx = self.mesh_data["irk_idx"]

        print(
            f"\nMaximum k-point mesh size from `inpsd.dat`: \n{mesh[0]} x {mesh[1]} x {mesh[2]}"
        )
        print(
            f"\nNumber of reduced k-points: {irk_idx.size} Written to file `qfile.reduced`."
        )

    def write_output_files(self):
        """Write all output files."""
        self._write_kpath_file()
        self._write_kpath2d_file()
        self._write_reduced_mesh_file()

    def _write_kpath_file(self):
        """Write k-path to qfile.kpath."""
        lpath = self.kpath_data["lpath"]
        nq = len(lpath)

        with open("qfile.kpath", "w", encoding="utf-8") as f:
            f.write(f"         {nq}\n")
            for row in lpath:
                f.write(
                    f" {row[0]:10.8f}  {row[1]:10.8f}  {row[2]:10.8f}    {row[3]}\n"
                )

    def _write_kpath2d_file(self):
        """Write 2D k-path to qfile.kpath2d."""
        lpath2d = self.kpath_data["lpath2d"]
        nq2d = len(lpath2d)

        with open("qfile.kpath2d", "w", encoding="utf-8") as qf:
            qf.write(f"         {nq2d}\n")
            for row in lpath2d:
                qf.write(
                    f" {row[0]:10.8f}  {row[1]:10.8f}  {row[2]:10.8f}    {row[3]}\n"
                )

    def _write_reduced_mesh_file(self):
        """Write reduced mesh to qfile.reduced."""
        mesh = self.inputs["mesh"]
        irk_idx = self.mesh_data["irk_idx"]
        grid = self.mesh_data["grid"]
        mult = self.mesh_data["mult"]

        with open("qfile.reduced", "w", encoding="utf-8") as qf:
            qf.write(f"         {irk_idx.size}\n")
            for i, (ir_gp_id, gp) in enumerate(zip(irk_idx, grid[irk_idx])):
                gp_str = " ".join(f"{f/mesh[j]: 10.8f}" for j, f in enumerate(gp))
                qf.write(f"{gp_str}     {mult[i]:10.4f}\n")

    def run(self):
        """Execute the complete pre-processing workflow."""
        print("Reading input parameters...")
        self.read_inputs()

        print("Processing atomic positions...")
        self.process_positions()

        print("Analyzing crystal structure...")
        self.analyze_structure()

        print("Generating k-path...")
        self.generate_kpath()

        print("Generating reduced k-point mesh...")
        self.generate_reduced_mesh()

        print("Printing structure information...")
        self.print_structure_info()

        print("Printing k-path information...")
        self.print_kpath_info()

        print("Printing mesh information...")
        self.print_mesh_info()

        print("Writing output files...")
        self.write_output_files()

        print("Pre-processing complete!")


if __name__ == "__main__":
    # Execute pre-processing
    processor = PreProcessor()
    processor.run()
