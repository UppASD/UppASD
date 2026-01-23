#!/usr/bin/env python
"""
ASDImportSys.py

This module provides functions to parse .cif and .nml files and convert them into UppASD format.

    Parse .cif files and return them in UppASD format.

    Args:
        filename (str): Path to the .cif file.

    Returns:
        tuple: A tuple containing filenames, lattice, and alat.
    ...

    Read important keywords from UppASD input file `inpsd.dat`.

    Args:
        ifile (str): Path to the `inpsd.dat` file.

    Returns:
        tuple: A tuple containing lattice and alat.
    ...

    Parse .nml files and return them in UppASD format.

    Args:
        filename (str): Path to the .nml file.

    Returns:
        tuple: A tuple containing filenames, lattice_type, and alat.
    ...

    Read nml-file into dictionary with atom number, magnitude, and vector.

    Args:
        element_list (list): List of elements.
        path (str): Path to the directory containing .nml files.
        lattice (str): Lattice type.

    Returns:
        dict: A dictionary with moment information.
    ...

    Create momfile from dictionary.

    Args:
        momfile_dict (dict): Dictionary containing moment information.

    Returns:
        None
    ...
"""
# pylint: disable=invalid-name, no-name-in-module, no-member
# coding: utf-8
import os
import copy
import shutil
import subprocess
import tempfile

import f90nml
import numpy as np


def parse_cif(filename):
    """Parse .cif files and return them in UppASD format"""
    current_dir = os.getcwd()
    temp_dir = tempfile.TemporaryDirectory()
    temp_path = temp_dir.name

    shutil.copy2(filename, temp_path)

    subprocess.run("cif2cell -p uppasd -f " + filename, cwd=temp_path, shell=True, check=True)
    # subprocess.run("cif2cell -p uppasd -f FeCo_56273.cif", cwd=temp_path, shell=True, check=True)

    lattice, alat = read_inpsd(os.path.join(temp_path, "inpsd.dat"))

    shutil.copy2(os.path.join(temp_path, "posfile"), current_dir)
    shutil.copy2(os.path.join(temp_path, "momfile"), current_dir)

    temp_dir.cleanup()

    filenames = ["", "", "posfile", "momfile"]

    return filenames, lattice, alat


def read_inpsd(ifile):
    """Read important keywords from UppASD inputfile `inpsd.dat`"""
    with open(ifile, "r", encoding="utf-8") as infile:
        lines = infile.readlines()
        for idx, line in enumerate(lines):
            line_data = line.rstrip("\n").split()
            if len(line_data) > 0:
                # Find the cell data
                if line_data[0] == "cell":
                    lattice = np.empty([0, 3])
                    line_data = lines[idx + 1].split()
                    lattice = np.vstack((lattice, np.asarray(line_data[0:3])))
                    line_data = lines[idx + 2].split()
                    lattice = np.vstack((lattice, np.asarray(line_data[0:3])))
                    line_data = lines[idx + 3].split()
                    lattice = np.vstack((lattice, np.asarray(line_data[0:3])))
                    lattice = lattice.astype(np.float64)
                if line_data[0] == "alat":
                    alat = lines[idx].split()[-1]

    return lattice, alat


def parse_rs_lmto(filename):
    """:q"""
    current_dir = os.getcwd()
    path = "/".join(filename.split("/")[:-1])
    parser = f90nml.Parser()

    input_nml = parser.read(filename)

    elements = input_nml["atoms"]["label"]

    if current_dir != path:
        shutil.copy2(os.path.join(path, "jij.out"), current_dir)
        shutil.copy2(os.path.join(path, "dij.out"), current_dir)

    lattice_type = input_nml["lattice"]["crystal_sym"]
    alat = input_nml["lattice"]["alat"]
    filenames = ["./jij.out", "./dij.out", "./posfile", "./momfile"]

    if lattice_type == "data":
        raise NotImplementedError("data-file import not yet implemented!")

    momfile_dict = get_rs_lmto_moments(elements, path, lattice_type)
    write_momfile(momfile_dict)

    return filenames, lattice_type, alat


def get_rs_lmto_moments(element_list, path, lattice):
    """Read nml-file into dictonary with atom number, magnitude, vector
    TODO: The implementation for handling geometries with more than one
    particle in the unit cell is quite ugly as it is rigth now. Would like to
    read the information from the poosition file instead.
    """

    parser = f90nml.Parser()
    momfile_dict = {}

    for index, element in enumerate(element_list):
        data = parser.read(os.path.join(path, element + "_out.nml"))

        momfile_dict[element] = [[index + 1, 1]]
        momfile_dict[element].append(
            [
                np.sum(np.array(data["par"]["ql"][0])[:, 0])
                - np.sum(np.array(data["par"]["ql"][1])[:, 0])
            ]
        )
        momfile_dict[element].append(data["par"]["mom"])

    if lattice in ["bcc", "hcp"]:
        key = list(momfile_dict.keys())[0]
        momfile_dict[key + "2"] = copy.copy(momfile_dict[key])
        momfile_dict[key + "2"][0] = [2, 1]

    return momfile_dict


def write_momfile(momfile_dict):
    """Create momfile from dict."""
    momfile = open("momfile", "w", encoding="utf-8")
    momfile_string = ""

    for key in momfile_dict:
        for idx, inputs in enumerate(momfile_dict[key]):
            if idx == 0:
                momfile_string += "\t".join(str(inp) for inp in inputs) + "\t"
            else:
                momfile_string += (
                    "\t".join(f"{inp:2.5f}" for inp in inputs) + "\t"
                )
        momfile_string += "\n"

    momfile.write(momfile_string)
    momfile.close()
