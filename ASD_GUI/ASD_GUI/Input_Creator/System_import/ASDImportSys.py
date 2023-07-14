#!/usr/bin/env python
# coding: utf-8
import tempfile
import subprocess
import numpy as np
import shutil
import os
from PyQt6 import QtWidgets
import ASD_GUI.Input_Creator.System_import.SPRKKR_parser as SPRKKR_parser

def import_system(window, ASDInputGen):
    """Select import file."""
    dlg = QtWidgets.QFileDialog()
    dlg.setFileMode(QtWidgets.QFileDialog.FileMode.AnyFile)
    dlg.setDirectory('.')

    if dlg.exec():
        if window.sender() == window.InpImportCIFButton:
            filename = dlg.selectedFiles()[0]
            output_files, lattice = parse_cif(filename)
            cif = True
        if window.sender() == window.InpImportSPRKKRButton:
            filename = dlg.selectedFiles()[0]
            output_files, lattice = SPRKKR_parser.parse_sprkkr(filename)
            cif = False

        update_gui(window, ASDInputGen, output_files, lattice, cif)

def update_gui(window, ASDInputGen, output_files, lattice, cif):
    """ Update GUI with files and input from file-parser. """
    ASDInputGen.jfile = output_files[0]
    ASDInputGen.dmfile = output_files[1]
    ASDInputGen.posfile = output_files[2]
    ASDInputGen.momfile = output_files[3]

    if cif:
        ASDInputGen.UppASDKeywords['geometry']['posfiletype'] = 'D'
    
    gui_lines = np.array([line for line in window.findChildren(QtWidgets.QLineEdit)
                    if 'InpLineEditC' in line.objectName()]).reshape(3,3)

    for row, vector in enumerate(gui_lines):
        for column, element in enumerate(vector):
            element.setText(str(lattice[row, column]))

    return

def parse_cif(filename):
    """ Parse .cif files and return them in UppASD format"""
    current_dir = os.getcwd()
    temp_dir = tempfile.TemporaryDirectory()
    temp_path = temp_dir.name

    shutil.copy2(filename, temp_path)

    subprocess.run('cif2cell -p uppasd -f FeCo_56273.cif', cwd = temp_path, shell = True)

    lattice = read_inpsd(os.path.join(temp_path,'inpsd.dat'))
    shutil.copy2(os.path.join(temp_path, 'posfile'), current_dir)
    shutil.copy2(os.path.join(temp_path, 'momfile'), current_dir)

    temp_dir.cleanup()

    filenames = ['', '', 'posfile', 'momfile']

    return filenames, lattice

def read_inpsd(ifile):
    """Read important keywords from UppASD inputfile `inpsd.dat`"""
    with open(ifile, 'r', encoding="utf-8") as infile:
        lines = infile.readlines()
        for idx, line in enumerate(lines):
            line_data = line.rstrip('\n').split()
            if len(line_data) > 0:
                # Find the cell data
                if line_data[0] == 'cell':
                    lattice = np.empty([0, 3])
                    line_data = lines[idx+1].split()
                    lattice = np.vstack((lattice, np.asarray(line_data[0:3])))
                    line_data = lines[idx+2].split()
                    lattice = np.vstack((lattice, np.asarray(line_data[0:3])))
                    line_data = lines[idx+3].split()
                    lattice = np.vstack((lattice, np.asarray(line_data[0:3])))
                    lattice = lattice.astype(np.float64)

    return lattice