#!/usr/bin/env python3
# coding: utf-8

# In[1]:


import numpy as np
import os.path
from scipy import ndimage


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
                positions = np.vstack((positions, np.asarray(line_data[2:5])))
                numbers = np.append(numbers, np.asarray(line_data[1]))
        return positions, numbers


############################################################
# Read important keywords from UppASD inputfile `inpsd.dat`
############################################################
def read_inpsd(ifile):
    with open(ifile, "r") as infile:
        lines = infile.readlines()
        for idx, line in enumerate(lines):
            line_data = line.rstrip("\n").split()
            if len(line_data) > 0:
                # Find the simulation id
                if line_data[0] == "simid":
                    simid = line_data[1][0:8]
                    print("simid: ", simid)

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

                # Find the size of the simulated cell
                if line_data[0] == "ncell":
                    ncell_x = int(line_data[1])
                    ncell_y = int(line_data[1])
                    ncell_z = int(line_data[1])
                    mesh = [ncell_x, ncell_y, ncell_z]

                if line_data[0] == "timestep":
                    timestep = line_data[1]
                    print("timestep: ", timestep)

                if line_data[0] == "sc_nstep":
                    sc_nstep = line_data[1]
                    print("sc_nstep: ", sc_nstep)

                if line_data[0] == "sc_step":
                    sc_step = line_data[1]
                    print("sc_step: ", sc_step)

                # Read the name of the position file
                if line_data[0].strip() == "posfile":
                    positions, numbers = read_posfile(line_data[1])

    return lattice, positions, numbers, simid, mesh, timestep, sc_step, sc_nstep


# In[7]:


############################################################
# Open and read input files
############################################################
ifile = "inpsd.dat"
lattice, positions, numbers, simid, mesh, timestep, sc_step, sc_nstep = read_inpsd(
    ifile
)


hbar = 4.135667662e-15
sigma_w = 2.0
emax = 0.5 * float(hbar) / (float(timestep) * float(sc_step)) * 1e3
############################################################
# Read simulated dynamical structure factor (S(q,w))
############################################################
got_sqw = os.path.isfile("sqw." + simid + ".out")
if got_sqw:
    # sqw=np.loadtxt('sqw.'+simid+'.out')
    sqw = np.genfromtxt("sqw." + simid + ".out", usecols=(0, 4, 5, 6, 7, 8))
    nq = int(sqw[-1, 0])
    nw = int(sqw[-1, 1])
    q = sqw[:, 0]
    w = sqw[:, 1]
    sqw_x = np.reshape(sqw[:, 2], (nq, nw))
    sqw_y = np.reshape(sqw[:, 3], (nq, nw))
    sqw_z = np.reshape(sqw[:, 4], (nq, nw))
    sqw_tot = np.reshape(sqw[:, 3], (nq, nw))
    sqw_perp = np.sqrt(sqw_x**2 + sqw_y**2)

    # sqw_x[:,0]=sqw_x[:,0]/100.0
    sqw_x = sqw_x.T / sqw_x.T.max(axis=0)
    sqw_x = ndimage.gaussian_filter1d(sqw_x, sigma=sigma_w, axis=0, mode="reflect")
    sqw_x = np.reshape(sqw_x.T, nq * nw)
    # sqw_y[:,0]=sqw_y[:,0]/100.0
    sqw_y = sqw_y.T / sqw_y.T.max(axis=0)
    sqw_y = ndimage.gaussian_filter1d(sqw_y, sigma=sigma_w, axis=0, mode="reflect")
    sqw_y = np.reshape(sqw_y.T, nq * nw)
    # sqw_z[:,0]=sqw_z[:,0]/100.0
    sqw_z = sqw_z.T / sqw_z.T.max(axis=0)
    sqw_z = ndimage.gaussian_filter1d(sqw_z, sigma=sigma_w, axis=0, mode="reflect")
    sqw_z = np.reshape(sqw_z.T, nq * nw)
    # sqw_perp[:,0]=sqw_perp[:,0]/100.0
    sqw_perp = sqw_perp.T / sqw_perp.T.max(axis=0)
    sqw_perp = ndimage.gaussian_filter1d(
        sqw_perp, sigma=sigma_w, axis=0, mode="reflect"
    )
    sqw_perp = np.reshape(sqw_perp.T, nq * nw)
    # sqw_tot[:,0]=sqw_tot[:,0]/100.0
    sqw_tot = sqw_tot.T / sqw_tot.T.max(axis=0)
    sqw_tot = ndimage.gaussian_filter1d(sqw_tot, sigma=sigma_w, axis=0, mode="reflect")
    sqw_tot = np.reshape(sqw_tot.T, nq * nw)

    # np.savetxt('sqw.normalized.dat',(q,w,np.reshape(sqw_x,nq*nw)))
    np.savetxt(
        "sqw_norm." + simid + ".dat",
        np.column_stack(
            (q.astype(int), w.astype(int), sqw_x, sqw_y, sqw_z, sqw_tot, sqw_perp)
        ),
        fmt="%7i %7i   %10.4e %10.4e %10.4e   %10.4e  %10.4e",
    )

    # sqw_x=ndimage.gaussian_filter1d(sqw_x,sigma=1,axis=1,mode='constant')
    # sqw_x=ndimage.gaussian_filter1d(sqw_x,sigma=5,axis=0,mode='reflect')


# In[ ]:


# In[ ]:
