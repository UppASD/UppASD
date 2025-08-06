#!/usr/bin/env python
# coding: utf-8

from asd_io import (
    read_posfile,
    read_inpsd,
    get_reciprocal_lattice,
    get_spacegroup,
    get_symmetry_points,
    get_kpath,
)
import numpy as np
import matplotlib.pyplot as plt
import spglib as spg
import seekpath as spth
import matplotlib.cm as cmap
import os.path
from scipy import ndimage


# In[1]:


def is_close(A, B):
    check = True
    for a, b in zip(A, B):
        check = check and np.isclose(np.float64(a), np.float64(b), atol=1e-6)
    return check


hbar = 6.582119514e-13
ry_ev = 13.605693009

sigma_q = 0.5  # 0.5
sigma_w = 1.0  # 4.0


############################################################
# Set pyplot defaults
############################################################

font = {"family": "sans", "weight": "normal", "size": "14"}

plt.rc("font", **font)
plt.rc("lines", lw=2)


# In[3]:


# In[4]:


# plt.xkcd()


# In[5]:


# In[7]:


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
# Read adiabatic magnon spectra
############################################################
got_ams = os.path.isfile("ams." + simid + ".out")
got_ncams = os.path.isfile("ncams." + simid + ".out")
got_ncams_mq = os.path.isfile("ncams-q." + simid + ".out")
got_ncams_pq = os.path.isfile("ncams+q." + simid + ".out")
if got_ams:
    ams = np.loadtxt("ams." + simid + ".out")
    ams_dist_col = ams.shape[1] - 1
# else:
# got_ncams=os.path.isfile('ncams.'+simid+'.out')
if got_ncams:
    ams = np.loadtxt("ncams." + simid + ".out")
    ams_dist_col = ams.shape[1] - 1
if got_ncams_pq:
    ams_pq = np.loadtxt("ncams+q." + simid + ".out")
if got_ncams_mq:
    ams_mq = np.loadtxt("ncams-q." + simid + ".out")

# Set the x-vector
q_vecs = ams[:, 0]
q_min = np.min(q_vecs)
q_max = np.max(q_vecs)


if got_ams or got_ncams:
    emax_lswt = 1.10 * np.amax(ams[:, 1:ams_dist_col])
    # print('emax_lswt=',emax_lswt)


############################################################
# Read simulated dynamical structure factor (S(q,w))
############################################################
got_sqw = os.path.isfile("sqw." + simid + ".out")
if got_sqw:
    # sqw=np.loadtxt('sqw.'+simid+'.out')
    sqw = np.genfromtxt("sqw." + simid + ".out", usecols=(0, 4, 5, 6, 7, 8))
    nq = int(sqw[-1, 0])
    nw = int(sqw.shape[0] / nq)
    sqw_x = np.reshape(sqw[:, 2], (nq, nw))[:, :]
    sqw_y = np.reshape(sqw[:, 3], (nq, nw))[:, :]
    sqw_z = np.reshape(sqw[:, 4], (nq, nw))[:, :]
    sqw_t = sqw_x**2 + sqw_y**2


############################################################
# Read simulated dynamical structure factor (S(q,w)) intensity
############################################################
got_sqw_int = os.path.isfile("sqwintensity." + simid + ".out")
if got_sqw_int:
    sqw_t_int = np.genfromtxt("sqwintensity." + simid + ".out", usecols=(0, 4, 5, 6))
    nq_int = int(sqw_t_int[-1, 0])
    nw_int = int(sqw_t_int.shape[0] / nq_int)
    sqw_int = np.reshape(sqw_t_int[:, 2], (nq_int, nw_int))


############################################################
# Read linear spin wave theory structure factor (S(q,w)) intensity
############################################################
got_sqw_lint = os.path.isfile("ncsqw_intensity." + simid + ".out")
if got_sqw_lint:
    sqw_lt_int = np.genfromtxt(
        "ncsqw_intensity." + simid + ".out", usecols=(0, 4, 5, 6)
    )
    nq_lint = int(sqw_lt_int[-1, 0])
    nw_lint = int(sqw_lt_int.shape[0] / nq_lint)
    sqw_lint = np.reshape(sqw_lt_int[:, 2], (nq_lint, nw_lint))


############################################################
# Read simulated dynamical structure factor (S(q,w))
############################################################
got_sqw_tens = os.path.isfile("sqwtensa." + simid + ".out")
if got_sqw_tens:
    sqwt = np.genfromtxt(
        "sqwtensa." + simid + ".out", usecols=(0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)
    )
    nqt = int(sqwt[-1, 0])
    nwt = int(sqwt.shape[0] / nqt)
    sqw_tens = np.zeros((nqt, nwt, 3, 3))
    sqw_tens[:, :, 0, 0] = np.reshape(sqwt[:, 2], (nqt, nwt))
    sqw_tens[:, :, 0, 1] = np.reshape(sqwt[:, 3], (nqt, nwt))
    sqw_tens[:, :, 0, 2] = np.reshape(sqwt[:, 4], (nqt, nwt))
    sqw_tens[:, :, 1, 0] = np.reshape(sqwt[:, 5], (nqt, nwt))
    sqw_tens[:, :, 1, 1] = np.reshape(sqwt[:, 6], (nqt, nwt))
    sqw_tens[:, :, 1, 2] = np.reshape(sqwt[:, 7], (nqt, nwt))
    sqw_tens[:, :, 2, 0] = np.reshape(sqwt[:, 8], (nqt, nwt))
    sqw_tens[:, :, 2, 1] = np.reshape(sqwt[:, 9], (nqt, nwt))
    sqw_tens[:, :, 2, 2] = np.reshape(sqwt[:, 10], (nqt, nwt))
    # sqw_t=np.reshape(sqw[:,3],(nq,nw))
    # sqw_t=sqw_x**2+sqw_y**2


############################################################
# Read simulated dynamical structure factor (S(q,w))
############################################################
got_lswt_sqw_tens = os.path.isfile("ncsqw." + simid + ".out")
if got_lswt_sqw_tens:
    lswt_sqwt = np.genfromtxt(
        "ncsqw." + simid + ".out", usecols=(0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)
    )
    lswt_nqt = int(lswt_sqwt[-1, 0])
    lswt_nwt = int(lswt_sqwt.shape[0] / lswt_nqt)
    lswt_sqw_tens = np.zeros((lswt_nqt, lswt_nwt, 3, 3))
    lswt_sqw_tens[:, :, 0, 0] = np.reshape(lswt_sqwt[:, 2], (lswt_nqt, lswt_nwt))
    lswt_sqw_tens[:, :, 0, 1] = np.reshape(lswt_sqwt[:, 3], (lswt_nqt, lswt_nwt))
    lswt_sqw_tens[:, :, 0, 2] = np.reshape(lswt_sqwt[:, 4], (lswt_nqt, lswt_nwt))
    lswt_sqw_tens[:, :, 1, 0] = np.reshape(lswt_sqwt[:, 5], (lswt_nqt, lswt_nwt))
    lswt_sqw_tens[:, :, 1, 1] = np.reshape(lswt_sqwt[:, 6], (lswt_nqt, lswt_nwt))
    lswt_sqw_tens[:, :, 1, 2] = np.reshape(lswt_sqwt[:, 7], (lswt_nqt, lswt_nwt))
    lswt_sqw_tens[:, :, 2, 0] = np.reshape(lswt_sqwt[:, 8], (lswt_nqt, lswt_nwt))
    lswt_sqw_tens[:, :, 2, 1] = np.reshape(lswt_sqwt[:, 9], (lswt_nqt, lswt_nwt))
    lswt_sqw_tens[:, :, 2, 2] = np.reshape(lswt_sqwt[:, 10], (lswt_nqt, lswt_nwt))
    # sqw_t=np.reshape(sqw[:,3],(nq,nw))
    # sqw_t=sqw_x**2+sqw_y**2


############################################################
# Get the spacegroup from spglib and print relevant info
############################################################
cell = (lattice, positions, numbers)
spacegroup = get_spacegroup(cell)


# In[8]:


############################################################
# Get symmetry points from seekpath
############################################################
BZ, sympoints = get_symmetry_points(cell)


# In[9]:
############################################################
# Read the qpoint-file used for the simulation
############################################################

###############################################################
#### Extract symmetry points and their location for plotting
###############################################################

# if(qfile=='qfile.kpath' or qfile=='qfile.kpath2d'):
if qfile and "qfile.kpath" in qfile:
    with open(qfile, "r") as f:
        f.readline()
        qpts = f.readlines()

    axlab = []
    axidx = []
    for idx, row in enumerate(qpts):
        rs = row.split()
        # Expect: q_x q_y q_z label
        if len(rs) == 4:
            label = rs[3]
            axlab.append("$\\Gamma$" if label[0] == "G" else label)
            axidx.append(q_vecs[idx])
    # axidx now contains the x-positions for symmetry labels, axlab the labels

# In[12]:


############################################################
# Plot the AMS
############################################################
if got_ams or got_ncams:
    plt.figure(figsize=[8, 5])

    plt.plot(q_vecs[:], ams[:, 1:ams_dist_col])
    ala = plt.xticks()

    plt.xticks(axidx, axlab)
    # plt.xlabel('q')
    plt.ylabel("Energy (meV)")

    plt.autoscale(tight=True)
    plt.ylim(0)
    plt.grid(visible=True, which="major", axis="x")
    # plt.show()
    plt.savefig("ams.png")

    plt.close()

if got_ncams:

    plt.figure(figsize=[8, 5])

    plt.plot(q_vecs[:], ams[:, 1:ams_dist_col], label="E(q)")
    if got_ncams_pq:
        plt.plot(q_vecs[:], ams_pq[:, 1:ams_dist_col], label="E(q+q$_0$)")
    if got_ncams_mq:
        plt.plot(q_vecs[:], ams_mq[:, 1:ams_dist_col], label="E(q-q$_0$)")
    ala = plt.xticks()

    plt.xticks(axidx, axlab)
    # plt.xlabel('q')
    plt.legend()
    plt.ylabel("Energy (meV)")

    plt.autoscale(tight=True)
    plt.ylim(0)
    plt.grid(visible=True, which="major", axis="x")
    # plt.show()
    plt.savefig("ams_q.png")

    plt.close()


# In[ ]:

############################################################
# Plot the S(q,w)
############################################################
if got_sqw:
    fig = plt.figure(figsize=[8, 5])
    ax = plt.subplot(111)

    hbar = 4.135667662e-15
    emax = 0.5 * np.float64(hbar) / (np.float64(timestep) * np.float64(sc_step)) * 1e3
    sqw_temp = (sqw_x**2 + sqw_y**2) ** 0.5
    sqw_temp[:, 0] = sqw_temp[:, 0] / 100.0
    sqw_temp = sqw_temp.T / sqw_temp.T.max(axis=0)
    sqw_temp = ndimage.gaussian_filter1d(
        sqw_temp, sigma=sigma_q, axis=1, mode="constant"
    )
    sqw_temp = ndimage.gaussian_filter1d(
        sqw_temp, sigma=sigma_w, axis=0, mode="reflect"
    )
    # plt.imshow(sqw_temp, cmap=cmap.gist_ncar_r, interpolation='nearest',origin='lower',extent=[axidx_abs[0],axidx_abs[-1],0,emax])
    plt.imshow(
        sqw_temp,
        cmap=cmap.gist_ncar_r,
        interpolation="nearest",
        origin="lower",
        extent=[q_min, q_max, 0, emax],
    )
    plt.plot(q_vecs[:], ams[:, 1:ams_dist_col], "black", lw=0.5)
    ala = plt.xticks()

    plt.xticks(axidx, axlab)
    plt.xlabel("q")
    plt.ylabel("Energy (meV)")

    plt.grid(visible=True, which="major", axis="x")
    ax.set_aspect("auto")
    plt.autoscale(tight=True)
    # plt.show()
    plt.savefig("ams_sqw.png")


############################################################
# Plot the S(q,w) with full nc-LSWT support
############################################################
if got_sqw and got_ncams and got_ncams_pq and got_ncams_mq:
    fig = plt.figure(figsize=[8, 5])
    ax = plt.subplot(111)

    hbar = 4.135667662e-15
    emax = 0.5 * np.float64(hbar) / (np.float64(timestep) * np.float64(sc_step)) * 1e3
    sqw_temp = (sqw_x**2 + sqw_y**2) ** 0.5
    sqw_temp[:, 0] = sqw_temp[:, 0] / 100.0
    sqw_temp = sqw_temp.T / sqw_temp.T.max(axis=0)
    sqw_temp = ndimage.gaussian_filter1d(
        sqw_temp, sigma=sigma_q, axis=1, mode="constant"
    )
    sqw_temp = ndimage.gaussian_filter1d(
        sqw_temp, sigma=sigma_w, axis=0, mode="reflect"
    )
    imx_min = np.min(ams[:, 0] / ams[-1, 0] * axidx[-1])
    imx_max = np.max(ams[:, 0] / ams[-1, 0] * axidx[-1])
    # plt.imshow(sqw_temp, cmap=cmap.gist_ncar_r, interpolation='nearest',origin='lower',extent=[imx_min,imx_max,0,emax])
    plt.imshow(
        sqw_temp,
        cmap=cmap.gist_ncar_r,
        interpolation="nearest",
        origin="lower",
        extent=[q_min, q_max, 0, emax],
    )
    # plt.imshow(sqw_temp, cmap=cmap.gist_ncar_r, interpolation='nearest',origin='lower',extent=[axidx_abs[0],axidx_abs[-1],0,emax])
    # plt.imshow(sqw_temp, cmap=cmap.gist_ncar_r, interpolation='nearest',origin='lower',extent=[axidx_abs[0],axidx_abs[-1],0,emax])
    plt.plot(q_vecs[:], ams[:, 1:ams_dist_col], "black", lw=0.5)
    plt.plot(q_vecs[:], ams_pq[:, 1:ams_dist_col], "black", lw=0.5)
    plt.plot(q_vecs[:], ams_mq[:, 1:ams_dist_col], "black", lw=0.5)
    ala = plt.xticks()

    plt.xticks(axidx, axlab)
    # print('::::>',[np.min(ams[:,0]/ams[-1,0]*axidx_abs[-1]), np.max(ams[:,0]/ams[-1,0]*axidx_abs[-1])])
    plt.xlabel("q")
    plt.ylabel("Energy (meV)")

    plt.grid(visible=True, which="major", axis="x")
    ax.set_aspect("auto")
    plt.autoscale(tight=True)
    # plt.xlim([np.min(ams[:,0]/ams[-1,0]*axidx_abs[-1]), np.max(ams[:,0]/ams[-1,0]*axidx_abs[-1])])
    plt.savefig("ams_sqw_q.png")


############################################################
# Plot the simulated S(q,w) intensity
############################################################
if got_sqw_int:
    fig = plt.figure(figsize=[8, 5])
    ax = plt.subplot(111)

    hbar = 4.135667662e-15
    emax = 0.5 * np.float64(hbar) / (np.float64(timestep) * np.float64(sc_step)) * 1e3
    sqw_temp = sqw_int
    sqw_temp[:, 0] = sqw_temp[:, 0] / 100.0
    sqw_temp = sqw_temp.T / sqw_temp.T.max(axis=0)
    sqw_temp = ndimage.gaussian_filter1d(
        sqw_temp, sigma=sigma_q, axis=1, mode="constant"
    )
    sqw_temp = ndimage.gaussian_filter1d(
        sqw_temp, sigma=sigma_w, axis=0, mode="reflect"
    )
    plt.imshow(
        sqw_temp,
        cmap=cmap.gist_ncar_r,
        interpolation="nearest",
        origin="lower",
        extent=[axidx[0], axidx[-1], 0, emax],
    )
    plt.plot(
        ams[:, 0] / ams[-1, 0] * axidx[-1], ams[:, 1:ams_dist_col], "black", lw=1
    )
    ala = plt.xticks()

    plt.xticks(axidx, axlab)
    plt.xlabel("q")
    plt.ylabel("Energy (meV)")

    plt.autoscale(tight=False)
    ax.set_aspect("auto")
    plt.grid(visible=True, which="major", axis="x")
    # plt.show()
    plt.savefig("ams_sqw_int.png")


############################################################
# Plot the LSWT S(q,w) intensity
############################################################
if got_sqw_lint:
    fig = plt.figure(figsize=[8, 5])
    ax = plt.subplot(111)

    hbar = 4.135667662e-15
    emax = 0.5 * np.float64(hbar) / (np.float64(timestep) * np.float64(sc_step)) * 1e3
    sqw_temp = sqw_lint
    sqw_temp[:, 0] = sqw_temp[:, 0] / 100.0
    sqw_temp = sqw_temp.T / sqw_temp.T.max(axis=0)
    sqw_temp = ndimage.gaussian_filter1d(
        sqw_temp, sigma=sigma_q, axis=1, mode="constant"
    )
    sqw_temp = ndimage.gaussian_filter1d(
        sqw_temp, sigma=sigma_w, axis=0, mode="reflect"
    )
    plt.imshow(
        sqw_temp,
        cmap=cmap.gist_ncar_r,
        interpolation="nearest",
        origin="lower",
        extent=[axidx[0], axidx[-1], 0, emax_lswt],
    )
    plt.plot(
        ams[:, 0] / ams[-1, 0] * axidx[-1], ams[:, 1:ams_dist_col], "black", lw=1
    )
    ala = plt.xticks()

    plt.xticks(axidx, axlab)
    plt.xlabel("q")
    plt.ylabel("Energy (meV)")

    plt.autoscale(tight=False)
    ax.set_aspect("auto")
    plt.grid(visible=True, which="major", axis="x")
    # plt.show()
    plt.savefig("ams_ncsqw_int.png")


xyz = ("x", "y", "z")

############################################################
# Plot the S(q,w) on tensorial form (diagonal terms)
############################################################
if got_sqw_tens:
    fig = plt.figure(figsize=[16, 4])
    hbar = 4.135667662e-15
    emax = 0.5 * np.float64(hbar) / (np.float64(timestep) * np.float64(sc_step)) * 1e3

    plt.xticks(axidx, axlab)
    plt.xlabel("q")
    plt.ylabel("Energy (meV)")

    # plt.autoscale(tight=False)
    ax.set_aspect("auto")
    plt.grid(visible=True, which="major", axis="x")

    # sqw_x[:,0]=sqw_x[:,0]/100.0
    # sqw_x=sqw_x.T/sqw_x.T.max(axis=0)
    plt_idx = 130
    for ix in range(3):
        iy = ix
        sqw_temp = sqw_tens[:, :, ix, iy]
        sqw_temp[:, 0] = sqw_temp[:, 0] / 1e5
        sqw_temp = sqw_temp.T / sqw_temp.T.max(axis=0)
        plt_idx = plt_idx + 1
        ax = plt.subplot(plt_idx)
        stitle = "$S^{" + xyz[ix] + xyz[iy] + "}(q,\\omega)$"  # Fixed LaTeX syntax
        ax.set_title(stitle)

        sqw_temp = ndimage.gaussian_filter1d(
            sqw_temp, sigma=sigma_q, axis=1, mode="constant"
        )
        sqw_temp = ndimage.gaussian_filter1d(
            sqw_temp, sigma=sigma_w, axis=0, mode="reflect"
        )
        plt.imshow(
            sqw_temp,
            cmap=cmap.gist_ncar_r,
            interpolation="nearest",
            origin="lower",
            extent=[axidx[0], axidx[-1], 0, emax],
            aspect="auto",
        )
        plt.plot(
            ams[:, 0] / ams[-1, 0] * axidx[-1],
            ams[:, 1:ams_dist_col],
            "black",
            lw=1,
        )
        plt.xticks(axidx, axlab)

    plt.tight_layout()

    plt.savefig("sqw_diagonal.png")


############################################################
# Plot the S(q,w) on tensorial form (full tensor)
############################################################
if got_sqw_tens:
    fig = plt.figure(figsize=[16, 10])
    hbar = 4.135667662e-15
    emax = 0.5 * np.float64(hbar) / (np.float64(timestep) * np.float64(sc_step)) * 1e3

    plt.xticks(axidx, axlab)
    plt.xlabel("q")
    plt.ylabel("Energy (meV)")

    # plt.autoscale(tight=False)
    ax.set_aspect("auto")
    plt.grid(visible=True, which="major", axis="x")

    # sqw_x[:,0]=sqw_x[:,0]/100.0
    # sqw_x=sqw_x.T/sqw_x.T.max(axis=0)
    plt_idx = 330
    for ix in range(3):
        for iy in range(3):
            sqw_temp = sqw_tens[:, :, ix, iy]
            sqw_temp[:, 0] = sqw_temp[:, 0] / 1e5
            sqw_temp = sqw_temp.T / sqw_temp.T.max(axis=0)
            plt_idx = plt_idx + 1
            ax = plt.subplot(plt_idx)
            stitle = "$S^{" + xyz[ix] + xyz[iy] + "}(q,\\omega)$"  # Fixed LaTeX syntax
            ax.set_title(stitle)

            sqw_temp = ndimage.gaussian_filter1d(
                sqw_temp, sigma=sigma_q, axis=1, mode="constant"
            )
            sqw_temp = ndimage.gaussian_filter1d(
                sqw_temp, sigma=sigma_w, axis=0, mode="reflect"
            )
            plt.imshow(
                sqw_temp,
                cmap=cmap.gist_ncar_r,
                interpolation="nearest",
                origin="lower",
                extent=[axidx[0], axidx[-1], 0, emax],
                aspect="auto",
            )
            plt.plot(
                ams[:, 0] / ams[-1, 0] * axidx[-1],
                ams[:, 1:ams_dist_col],
                "black",
                lw=1,
            )
            plt.xticks(axidx, axlab)

    plt.tight_layout()

    plt.savefig("sqw_tensor.png")


############################################################
# Plot the LSWT S(q,w) on tensorial form (only diagonal terms)
############################################################
if got_lswt_sqw_tens:
    fig = plt.figure(figsize=[16, 4])
    hbar = 4.135667662e-15
    emax = 0.5 * np.float64(hbar) / (np.float64(timestep) * np.float64(sc_step)) * 1e3

    plt.xticks(axidx, axlab)
    plt.xlabel("q")
    plt.ylabel("Energy (meV)")

    # plt.autoscale(tight=False)
    # ax.set_aspect('auto')
    plt.grid(visible=True, which="major", axis="x")

    # sqw_x[:,0]=sqw_x[:,0]/100.0
    # sqw_x=sqw_x.T/sqw_x.T.max(axis=0)
    plt_idx = 130
    for ix in range(3):
        iy = ix
        sqw_temp = lswt_sqw_tens[:, :, ix, iy]
        sqw_temp = sqw_temp.T / (sqw_temp.T.max(axis=0) + 1.0e-20)
        plt_idx = plt_idx + 1
        ax = plt.subplot(plt_idx)
        stitle = "$S^{" + xyz[ix] + xyz[iy] + "}_{LSWT}$"
        ax.set_title(stitle)

        sqw_temp = ndimage.gaussian_filter1d(
            sqw_temp, sigma=sigma_q, axis=1, mode="constant"
        )
        sqw_temp = ndimage.gaussian_filter1d(
            sqw_temp, sigma=sigma_w, axis=0, mode="reflect"
        )
        plt.imshow(
            sqw_temp,
            cmap=cmap.gist_ncar_r,
            interpolation="nearest",
            origin="lower",
            extent=[axidx[0], axidx[-1], 0, emax_lswt],
            aspect="auto",
        )
        # plt.plot(ams[:,0]/ams[-1,0]*axidx[-1],ams[:,1:ams_dist_col],'black',lw=  1)
        plt.xticks(axidx, axlab)

    plt.tight_layout()

    plt.savefig("ncsqw_diagonal.png")


############################################################
# Plot the LSWT S(q,w) on tensorial form (full tensor)
############################################################
if got_lswt_sqw_tens:
    fig = plt.figure(figsize=[16, 10])
    hbar = 4.135667662e-15
    emax = 0.5 * np.float64(hbar) / (np.float64(timestep) * np.float64(sc_step)) * 1e3

    plt.xticks(axidx, axlab)
    plt.xlabel("q")
    plt.ylabel("Energy (meV)")

    # plt.autoscale(tight=False)
    # ax.set_aspect('auto')
    plt.grid(visible=True, which="major", axis="x")

    # sqw_x[:,0]=sqw_x[:,0]/100.0
    # sqw_x=sqw_x.T/sqw_x.T.max(axis=0)
    plt_idx = 330
    for ix in range(3):
        for iy in range(3):
            sqw_temp = lswt_sqw_tens[:, :, ix, iy]
            sqw_temp = sqw_temp.T / (sqw_temp.T.max(axis=0) + 1.0e-20)
            plt_idx = plt_idx + 1
            ax = plt.subplot(plt_idx)
            stitle = "$S^{" + xyz[ix] + xyz[iy] + "}_{LSWT}$"
            ax.set_title(stitle)

            sqw_temp = ndimage.gaussian_filter1d(
                sqw_temp, sigma=sigma_q, axis=1, mode="constant"
            )
            sqw_temp = ndimage.gaussian_filter1d(
                sqw_temp, sigma=sigma_w, axis=0, mode="reflect"
            )
            plt.imshow(
                sqw_temp,
                cmap=cmap.gist_ncar_r,
                interpolation="nearest",
                origin="lower",
                extent=[axidx[0], axidx[-1], 0, emax_lswt],
                aspect="auto",
            )
            # plt.plot(ams[:,0]/ams[-1,0]*axidx[-1],ams[:,1:ams_dist_col],'black',lw=  1)
            plt.xticks(axidx, axlab)

    plt.tight_layout()

    plt.savefig("ncsqw_tensor.png")
