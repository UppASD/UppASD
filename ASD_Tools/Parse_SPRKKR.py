#!  /usr/bin/env python
################################################################################
# @author Jonathan Chico
# Script to plot the Jij's and Dij's obtained from SPRKKR. The script also creates
# Input files that are compatible with UppASD package for atomistic spin dynamics
# simulations.
# It works by reading the _JXC_Jij.dat raw data file created by SPRKKR versions 6.3,7.1 and 7.7
# (and if DMI is calculated it also reads the JXC_Dij_(x,y,z).dat files).
# The script generates legends with appropriate names, allows one to choose the
# energy range of plotting and allows one to  decide whether the parameters
# of an element should be plotted.
# For extra information the script will read the .sys file and the _SCF.out files
# to generate more appropriate labels, and UppASD input files.
# Tested with python 2.7 and 3.6
################################################################################

import glob
import string
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm as cm

plt.rc("text", usetex=True)
plt.rc("font", family="serif")


################################################################################
# Find the number of lines of the file
################################################################################
def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


################################################################################
# Control variables
################################################################################
maptype = 2  # Type of maptype used for UppASD input
Plot_Jij = "Y"  # Plot the exchange interactions
exchange_radius = 4.0  # Radius up to which the exchange interactions are plotted
exclude_list = [
    "Vc",
    "Nb",
    "Ta",
    "Se",
    "S",
]  # Exclusion list for elements to not be considered
font_size = 28  # Fontsize for plots
################################################################################
# Define the formats for printing the data
################################################################################
if maptype == 1:
    jfile_format = (
        "{:5.0f}  {:5.0f}  {: 4.10f}  {: 4.10f}  {: 4.10f}  {: 4.8f}  {:4.10f}\n"
    )
    dmfile_format = "{:5.0f}  {:5.0f}  {: 4.10f}  {: 4.10f}  {: 4.10f}  {: 4.8f}  {: 4.8f}  {: 4.8f}  {:4.10f}\n"
elif maptype == 2:
    jfile_format = (
        "{:5.0f}  {:5.0f}  {: 5.0f}  {: 5.0f}  {: 5.0f}  {: 4.8f}  {:4.10f}\n"
    )
    dmfile_format = "{:5.0f}  {:5.0f}  {: 5.0f}  {: 5.0f}  {: 5.0f}  {: 4.8f}  {: 4.8f}  {: 4.8f}  {:4.10f}\n"
posfile_format = "{:5.0f} {:5.0f}  {: 4.10f}  {: 4.10f}  {: 4.10f}\n"
momfile_format = "{:5.0f} {:5.0f}  {: 4.10f}  {: 4.10f}  {: 4.10f}  {: 4.10f}\n"
################################################################################
# Actual reading of the files
################################################################################
JXC_File_name = str(input("Input JXC file name: "))
# This is the name of the file that we want to plot
JXC_File = open(JXC_File_name + "_JJij.dat")

# Set a flag to check the start of a new type of data to be False
Data_start_found = False
# Intiate the counter
count = 0
# Declare an array that contains the labels of the atomic species found in the calculation
atom_types = []

# Check for certain key data in the file
while not Data_start_found:
    count = count + 1
    line = JXC_File.readline()
    data = str.split(line)
    # Find the name of the system to be used in the output name
    if len(data) > 3 and data[3] == "NQ":
        Num_NQ = int(data[5])
    if len(data) > 3 and data[3] == "NT":
        Num_NT = int(data[5])
    if len(data) > 0 and data[0] == "IT":
        Data_start_found = True

JXC_File.seek(0)
################################################################################
# Read the actual data for the Jijs
################################################################################
data = pd.read_csv(JXC_File, skiprows=count, header=None, delim_whitespace=True).values
non_zero_ind = np.where(data[:, 10] > 0)
bond_type1 = np.zeros([len(non_zero_ind[0]), 3], dtype=np.float64)
bond_type1[:, 0] = np.asarray(data[non_zero_ind[0], 7], dtype=np.float64)
bond_type1[:, 1] = np.asarray(data[non_zero_ind[0], 8], dtype=np.float64)
bond_type1[:, 2] = np.asarray(data[non_zero_ind[0], 9], dtype=np.float64)
bond_type2 = np.zeros([len(non_zero_ind[0]), 3], dtype=np.int)
bond_type2[:, 0] = np.asarray(data[non_zero_ind[0], 4], dtype=np.int)
bond_type2[:, 1] = np.asarray(data[non_zero_ind[0], 5], dtype=np.int)
bond_type2[:, 2] = np.asarray(data[non_zero_ind[0], 6], dtype=np.int)
itype = np.asarray(data[non_zero_ind[0], 0], dtype=np.int)
isite = np.asarray(data[non_zero_ind[0], 1], dtype=np.int)
jtype = np.asarray(data[non_zero_ind[0], 2], dtype=np.int)
jsite = np.asarray(data[non_zero_ind[0], 3], dtype=np.int)
mod_ij = np.asarray(data[non_zero_ind[0], 10], dtype=np.float64)
Jij_mRy = np.asarray(data[non_zero_ind[0], 11], dtype=np.float64)
del data
JXC_File.close()
################################################################################
# Check if the .sys file exists and try to get information from it
################################################################################
sys_file = glob.glob("*.sys")
if len(sys_file) > 0:
    file_size = file_len(sys_file[0])
    data_label = np.genfromtxt(
        sys_file[0], skip_header=file_size - Num_NT, usecols=(2,), dtype=str
    )
    data_label = data_label.reshape(data_label.size)
    labels = []
    trim_labels = []
    for ii in range(0, len(data_label)):
        curr_string = "$" + str(data_label[ii]) + "$"
        ind = str(data_label[ii]).find("_")
        if ind > 0:
            trim_labels.append(str(data_label[ii])[0:ind])
        else:
            trim_labels.append(str(data_label[ii]))
        labels.append(curr_string)
    ############################################################################
    # Read the sys file again to find positions of the atoms and the system name
    ############################################################################
    count = 0
    SYSFile = open(sys_file[0])
    icl = np.zeros(Num_NQ, dtype=np.int)
    iq = np.zeros(Num_NQ, dtype=np.int)
    pos = np.zeros([Num_NQ, 3], dtype=np.float32)
    Data_start_found = False
    while count < Num_NQ:
        line = SYSFile.readline()
        data = str.split(line)
        # Find the name of the system to be used in the output name
        if len(data) > 0 and data[0] == "IQ":
            Data_start_found = True
            line = SYSFile.readline()
            data = str.split(line)
        if Data_start_found:
            iq[count] = data[0]
            icl[count] = data[1]
            pos[count, 0] = data[2]
            pos[count, 1] = data[3]
            pos[count, 2] = data[4]
            count = count + 1
    num_non_exclude = 0
    ############################################################################
    # Print the posfile
    ############################################################################
    posfile = open("posfile.dat", "w")
    for ii in range(0, Num_NQ):
        if trim_labels[icl[ii] - 1] not in exclude_list:
            num_non_exclude += 1
            posfile.write(
                posfile_format.format(
                    ii + 1, icl[ii], pos[ii, 0], pos[ii, 1], pos[ii, 2]
                )
            )
    posfile.close()
    SYSFile.close()
else:
    labels = []
    icl = np.zeros(Num_NT, dtype=np.int)
    for ii in range(0, Num_NT):
        labels.append("$Dum_{" + str(ii + 1) + "}$")
        icl[ii] = ii + 1
    trim_labels = labels
################################################################################
# Print the jfile
################################################################################
jfile_name = "jfile.dat"
jfile = open(jfile_name, "w")
for ii in range(len(itype)):
    if (trim_labels[itype[ii] - 1] not in exclude_list) and (
        trim_labels[jtype[ii] - 1] not in exclude_list
    ):
        if maptype == 1:
            jfile.write(
                jfile_format.format(
                    isite[ii],
                    jsite[ii],
                    bond_type1[ii, 0],
                    bond_type1[ii, 1],
                    bond_type1[ii, 2],
                    Jij_mRy[ii],
                    mod_ij[ii],
                )
            )
        else:
            jfile.write(
                jfile_format.format(
                    isite[ii],
                    jsite[ii],
                    bond_type2[ii, 0],
                    bond_type2[ii, 1],
                    bond_type2[ii, 2],
                    Jij_mRy[ii],
                    mod_ij[ii],
                )
            )
jfile.close()
################################################################################
# Define the file names
################################################################################
Dij_File_name_x = glob.glob(JXC_File_name + "_Dij_x.dat")
Dij_File_name_y = glob.glob(JXC_File_name + "_Dij_y.dat")
Dij_File_name_z = glob.glob(JXC_File_name + "_Dij_z.dat")
################################################################################
# If the raw Dij files are present read them
################################################################################
if len(Dij_File_name_x) > 0 and len(Dij_File_name_y) > 0 and len(Dij_File_name_z) > 0:

    # Set a flag to check the start of a new type of data to be False
    Data_start_found = False
    # Intiate the counter
    count = 0
    # Declare an array that contains the labels of the atomic species found in the calculation
    atom_types = []

    Dij_File = open(Dij_File_name_x[0])
    # Check for certain key data in the file
    while not Data_start_found:
        count = count + 1
        line = Dij_File.readline()
        data = str.split(line)
        if len(data) > 0 and data[0] == "IT":
            Data_start_found = True
    Dij_File.seek(0)

    data_x = pd.read_csv(
        Dij_File_name_x[0],
        skiprows=count,
        header=None,
        delim_whitespace=True,
        usecols=[11],
    ).values
    data_y = pd.read_csv(
        Dij_File_name_y[0],
        skiprows=count,
        header=None,
        delim_whitespace=True,
        usecols=[11],
    ).values
    data_z = pd.read_csv(
        Dij_File_name_z[0],
        skiprows=count,
        header=None,
        delim_whitespace=True,
        usecols=[11],
    ).values
    data = pd.read_csv(
        Dij_File_name_x[0],
        skiprows=count,
        header=None,
        delim_whitespace=True,
        usecols=[10],
    ).values
    Dij_x = np.asarray(data_x[non_zero_ind[0]], dtype=np.float64)
    Dij_y = np.asarray(data_y[non_zero_ind[0]], dtype=np.float64)
    Dij_z = np.asarray(data_z[non_zero_ind[0]], dtype=np.float64)
    del data_x, data_y, data_z, data
    ############################################################################
    # Print the dmfile
    ############################################################################
    dmfile_name = "dmfile.dat"
    dmfile = open(dmfile_name, "w")
    for ii in range(len(itype)):
        if (trim_labels[itype[ii] - 1] not in exclude_list) and (
            trim_labels[jtype[ii] - 1] not in exclude_list
        ):
            if maptype == 1:
                dmfile.write(
                    dmfile_format.format(
                        isite[ii],
                        jsite[ii],
                        bond_type1[ii, 0],
                        bond_type1[ii, 1],
                        bond_type1[ii, 2],
                        float(Dij_x[ii]),
                        float(Dij_y[ii]),
                        float(Dij_z[ii]),
                        mod_ij[ii],
                    )
                )
            else:
                dmfile.write(
                    dmfile_format.format(
                        isite[ii],
                        jsite[ii],
                        bond_type2[ii, 0],
                        bond_type2[ii, 1],
                        bond_type2[ii, 2],
                        float(Dij_x[ii]),
                        float(Dij_y[ii]),
                        float(Dij_z[ii]),
                        mod_ij[ii],
                    )
                )
    dmfile.close()
    Dij_File.close()
################################################################################
# Find information regarding the SCF file
################################################################################
scffile = glob.glob("*_SCF.out")
if len(scffile) > 0:
    SCF_File = open(scffile[0])
    # Reading the SCF output file starting from the last line
    spin_mom = []
    for line in reversed(SCF_File.readlines()):
        data = str.split(line)
        if len(data) > 0:
            if data[0] == "sum":
                spin_mom.append(float(data[4]))
    momfile = open("momfile.dat", "w")
    for ii in range(0, Num_NT):
        if trim_labels[ii] not in exclude_list:
            ind = np.where(icl == (ii + 1))
            for jj in range(0, len(ind[0])):
                momfile.write(
                    momfile_format.format(
                        iq[ind[0][jj]], 1, spin_mom[Num_NT - 1 - ii], 0, 0, 1
                    )
                )
    momfile.close()
    SCF_File.close()
################################################################################
# Plot the exchange interactions
################################################################################
itype = itype - 1
jtype = jtype - 1
tol = 0.01
ind_cut = np.where(mod_ij <= exchange_radius)
mod_ij = mod_ij[ind_cut[0]]
Jij_mRy = Jij_mRy[ind_cut[0]]
itype = itype[ind_cut[0]]
jtype = jtype[ind_cut[0]]
bond_type1 = bond_type1[ind_cut[0], :]
if Plot_Jij == "Y":
    colors = cm.Paired(np.linspace(0, 1, 2 * num_non_exclude + 2))
    for ii in range(0, Num_NT):
        if trim_labels[ii] not in exclude_list:
            i_ind = np.where(itype == ii)
            fig = plt.figure()
            counter = 0
            for jj in range(0, Num_NT):
                j_ind = np.where(jtype[i_ind[0]] == jj)
                if len(spin_mom) > 0:
                    sign_i = np.sign(spin_mom[Num_NT - 1 - ii])
                    sign_j = np.sign(spin_mom[Num_NT - 1 - jj])
                    sign_ij = sign_i * sign_j
                else:
                    sign_ij = 1.0

                if trim_labels[jj] not in exclude_list:
                    counter += 1
                    plt.plot(
                        mod_ij[i_ind[0][j_ind[0]]],
                        sign_ij * Jij_mRy[i_ind[0][j_ind[0]]],
                        alpha=0.75,
                        lw=4,
                        c=colors[counter],
                        label=labels[ii] + "-" + labels[jj],
                    )
                    plt.scatter(
                        mod_ij[i_ind[0][j_ind[0]]],
                        sign_ij * Jij_mRy[i_ind[0][j_ind[0]]],
                        color=colors[counter],
                        alpha=0.75,
                        s=300,
                        lw=1.00,
                        edgecolor="black",
                    )

            ####################################################################
            # Plotting options
            ####################################################################
            extremum = max(abs(max(Jij_mRy[i_ind[0]])), abs(min(Jij_mRy[i_ind[0]])))
            extremum = extremum * 1.10
            plt.legend(fontsize=font_size, loc="upper right")
            plt.ylabel(r"$J_{ij}$ [mRy]", fontsize=font_size)
            plt.xlabel(r"$r_{ij}/a_{lat}$", fontsize=font_size)
            ax = plt.gca()
            ax.set_facecolor((1, 1, 1))
            ax.tick_params(axis="x", colors="black", labelsize=font_size, width=2)
            ax.tick_params(axis="y", colors="black", labelsize=font_size, width=2)
            ax.set_ylim(-extremum, extremum)
            plt.axhline(0, color="black", linestyle="--")
            for axis in ["top", "bottom", "left", "right"]:
                ax.spines[axis].set_linewidth(3)
            # We change the fontsize of minor ticks label
            plt.grid(False)
            fig.set_size_inches(18.5, 10.5)
            plt.savefig(
                "Jij_" + data_label[ii] + ".pdf",
                transparent=False,
                dpi=2540,
                bbox_inches="tight",
            )
            fig.clf()
            plt.clf()
            plt.cla()
            i_ind = np.where(itype == ii)
            i_ind_perp = np.where(abs(bond_type1[i_ind[0], 2]) > tol)
            i_ind_para = np.where(abs(bond_type1[i_ind[0], 2]) < tol)
            fig = plt.figure()
            counter = 0
            for jj in range(0, Num_NT):
                j_ind_perp = np.where(jtype[i_ind[0][i_ind_perp[0]]] == jj)
                j_ind_para = np.where(jtype[i_ind[0][i_ind_para[0]]] == jj)
                if len(spin_mom) > 0:
                    sign_i = np.sign(spin_mom[Num_NT - 1 - ii])
                    sign_j = np.sign(spin_mom[Num_NT - 1 - jj])
                    sign_ij = sign_i * sign_j
                else:
                    sign_ij = 1.0

                if trim_labels[jj] not in exclude_list:
                    dir_label = [r"r$_{ij}|_{z\neq 0}$", r"r$_{ij}|_{z=0}$"]
                    dir_modij = [
                        mod_ij[i_ind[0][i_ind_perp[0][j_ind_perp[0]]]],
                        mod_ij[i_ind[0][i_ind_para[0][j_ind_para[0]]]],
                    ]
                    dir_Jij_mRy = [
                        Jij_mRy[i_ind[0][i_ind_perp[0][j_ind_perp[0]]]],
                        Jij_mRy[i_ind[0][i_ind_para[0][j_ind_para[0]]]],
                    ]
                    for kk in range(0, 2):
                        counter += 1
                        plt.plot(
                            dir_modij[kk],
                            sign_ij * dir_Jij_mRy[kk],
                            alpha=0.75,
                            lw=4,
                            c=colors[counter],
                            label=labels[ii] + "-" + labels[jj] + ", " + dir_label[kk],
                        )
                        plt.scatter(
                            dir_modij[kk],
                            sign_ij * dir_Jij_mRy[kk],
                            color=colors[counter],
                            alpha=0.75,
                            s=300,
                            lw=1.00,
                            edgecolor="black",
                        )
            ####################################################################
            # Plotting options
            ####################################################################
            plt.legend(fontsize=font_size, loc="upper right")
            plt.ylabel(r"$J_{ij}$ [mRy]", fontsize=font_size)
            plt.xlabel(r"$r_{ij}/a_{lat}$", fontsize=font_size)
            ax = plt.gca()
            ax.set_facecolor((1, 1, 1))
            ax.tick_params(axis="x", colors="black", labelsize=font_size, width=2)
            ax.set_ylim(-extremum, extremum)
            ax.tick_params(axis="y", colors="black", labelsize=font_size, width=2)
            plt.axhline(0, color="black", linestyle="--")
            for axis in ["top", "bottom", "left", "right"]:
                ax.spines[axis].set_linewidth(3)
            # We change the fontsize of minor ticks label
            plt.grid(False)
            fig.set_size_inches(18.5, 10.5)
            plt.savefig(
                "Jij_" + data_label[ii] + "_decomposed.pdf",
                transparent=False,
                dpi=2540,
                bbox_inches="tight",
            )
            fig.clf()
            plt.clf()
            plt.cla()
            i_ind = np.where(itype == ii)
            i_ind_perp = np.where(abs(bond_type1[i_ind[0], 2]) > tol)
            i_ind_para = np.where(abs(bond_type1[i_ind[0], 2]) < tol)
            fig = plt.figure()
            counter = 0
            for jj in range(0, Num_NT):
                j_ind_perp = np.where(jtype[i_ind[0][i_ind_perp[0]]] == jj)
                j_ind_para = np.where(jtype[i_ind[0][i_ind_para[0]]] == jj)
                if len(spin_mom) > 0:
                    sign_i = np.sign(spin_mom[Num_NT - 1 - ii])
                    sign_j = np.sign(spin_mom[Num_NT - 1 - jj])
                    sign_ij = sign_i * sign_j
                else:
                    sign_ij = 1.0

                if trim_labels[jj] not in exclude_list:
                    dir_label = [r"r$_{ij}|_{z\neq 0}$", r"r$_{ij}|_{z=0}$"]
                    dir_modij = [
                        mod_ij[i_ind[0][i_ind_perp[0][j_ind_perp[0]]]],
                        mod_ij[i_ind[0][i_ind_para[0][j_ind_para[0]]]],
                    ]
                    dir_Jij_mRy = [
                        Jij_mRy[i_ind[0][i_ind_perp[0][j_ind_perp[0]]]],
                        Jij_mRy[i_ind[0][i_ind_para[0][j_ind_para[0]]]],
                    ]
                    for kk in range(0, 2):
                        counter += 1
                        plt.plot(
                            dir_modij[kk],
                            sign_ij * dir_Jij_mRy[kk] * dir_modij[kk] ** 2,
                            alpha=0.75,
                            lw=4,
                            c=colors[counter],
                            label=labels[ii] + "-" + labels[jj] + ", " + dir_label[kk],
                        )
                        plt.scatter(
                            dir_modij[kk],
                            sign_ij * dir_Jij_mRy[kk] * dir_modij[kk] ** 2,
                            color=colors[counter],
                            alpha=0.75,
                            s=300,
                            lw=1.00,
                            edgecolor="black",
                        )
            ####################################################################
            # Plotting options
            ####################################################################
            extremum = max(
                abs(max(Jij_mRy[i_ind[0]] * mod_ij[i_ind[0]] ** 2)),
                abs(min(Jij_mRy[i_ind[0]] * mod_ij[i_ind[0]] ** 2)),
            )
            extremum = extremum * 1.10
            plt.legend(fontsize=font_size, loc="upper right")
            plt.ylabel(
                r"$J_{ij} \left(\frac{r_{ij}}{a_{lat}}\right)^2$ [mRy]",
                fontsize=font_size,
            )
            plt.xlabel(r"$r_{ij}/a_{lat}$", fontsize=font_size)
            ax = plt.gca()
            ax.set_facecolor((1, 1, 1))
            ax.tick_params(axis="x", colors="black", labelsize=font_size, width=2)
            ax.set_ylim(-extremum, extremum)
            ax.tick_params(axis="y", colors="black", labelsize=font_size, width=2)
            plt.axhline(0, color="black", linestyle="--")
            for axis in ["top", "bottom", "left", "right"]:
                ax.spines[axis].set_linewidth(3)
            # We change the fontsize of minor ticks label
            plt.grid(False)
            fig.set_size_inches(18.5, 10.5)
            plt.savefig(
                "Jij_" + data_label[ii] + "_decomposed_RKKY.pdf",
                transparent=False,
                dpi=2540,
                bbox_inches="tight",
            )
