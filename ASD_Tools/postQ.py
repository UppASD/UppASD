#!/usr/bin/env python
# coding: utf-8

"""
UppASD Post-Processing Module

This module provides the PostProcessor class for analyzing and visualizing
UppASD simulation results including magnon spectra and structure factors.
"""

import os.path

import matplotlib.cm as cmap
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage

from asd_io import (
    get_spacegroup,
    get_symmetry_points,
    read_inpsd,
)


class PostProcessor:
    """Post-process UppASD simulation results: read data, analyze, and plot."""

    def __init__(self, input_file="inpsd.dat"):
        """Initialize the PostProcessor with input file."""
        self.input_file = input_file

        # Input data dictionary
        self.inputs = {}

        # Simulation results dictionary
        self.results = {
            "ams": None,
            "ams_pq": None,
            "ams_mq": None,
            "sqw_x": None,
            "sqw_y": None,
            "sqw_z": None,
            "sqw_t": None,
            "sqw_int": None,
            "sqw_lint": None,
            "sqw_tens": None,
            "lswt_sqw_tens": None,
        }

        # Symmetry and plotting data
        self.symmetry = {
            "cell": None,
            "spacegroup": None,
            "BZ": None,
            "sympoints": None,
            "axlab": [],
            "axidx": [],
        }

        # Plot configuration
        self.plot_config = {
            "sigma_q": 0.5,
            "sigma_w": 1.0,
            "font": {"family": "sans", "weight": "normal", "size": "14"},
            "xyz": ("x", "y", "z"),
        }

        # Constants
        self.constants = {
            "hbar": 6.582119514e-13,
            "ry_ev": 13.605693009,
            "hbar_eV_s": 4.135667662e-15,
        }

        self._setup_plotting()

    def _setup_plotting(self):
        """Configure matplotlib plotting defaults."""
        plt.rc("font", **self.plot_config["font"])
        plt.rc("lines", lw=2)

    def read_inputs(self):
        """Read input simulation parameters from inpsd.dat."""
        inpsd_data = read_inpsd(self.input_file)
        self.inputs = {
            "lattice": inpsd_data["lattice"],
            "positions": inpsd_data["positions"],
            "numbers": inpsd_data["numbers"],
            "simid": inpsd_data["simid"],
            "mesh": inpsd_data["mesh"],
            "posfiletype": inpsd_data["posfiletype"],
            "timestep": inpsd_data["timestep"],
            "sc_step": inpsd_data["sc_step"],
            "sc_nstep": inpsd_data["sc_nstep"],
            "qfile": inpsd_data["qfile"],
        }

    def read_magnon_spectra(self):
        """Read adiabatic magnon spectra files."""
        simid = self.inputs["simid"]

        # Check which files exist
        got_ams = os.path.isfile(f"ams.{simid}.out")
        got_ncams = os.path.isfile(f"ncams.{simid}.out")
        got_ncams_mq = os.path.isfile(f"ncams-q.{simid}.out")
        got_ncams_pq = os.path.isfile(f"ncams+q.{simid}.out")

        if got_ams:
            self.results["ams"] = np.loadtxt(f"ams.{simid}.out")
        if got_ncams:
            self.results["ams"] = np.loadtxt(f"ncams.{simid}.out")
        if got_ncams_pq:
            self.results["ams_pq"] = np.loadtxt(f"ncams+q.{simid}.out")
        if got_ncams_mq:
            self.results["ams_mq"] = np.loadtxt(f"ncams-q.{simid}.out")

        # Set derived quantities
        if self.results["ams"] is not None:
            ams = self.results["ams"]
            self.results["ams_dist_col"] = ams.shape[1] - 1
            self.results["q_vecs"] = ams[:, 0]
            self.results["q_min"] = np.min(self.results["q_vecs"])
            self.results["q_max"] = np.max(self.results["q_vecs"])
            self.results["emax_lswt"] = 1.10 * np.amax(
                ams[:, 1 : self.results["ams_dist_col"]]
            )

    def read_structure_factors(self):
        """Read dynamical structure factor files."""
        simid = self.inputs["simid"]

        # Read S(q,w) data
        if os.path.isfile(f"sqw.{simid}.out"):
            sqw = np.genfromtxt(f"sqw.{simid}.out", usecols=(0, 4, 5, 6, 7, 8))
            nq = int(sqw[-1, 0])
            nw = int(sqw.shape[0] / nq)
            self.results["sqw_x"] = np.reshape(sqw[:, 2], (nq, nw))[:, :]
            self.results["sqw_y"] = np.reshape(sqw[:, 3], (nq, nw))[:, :]
            self.results["sqw_z"] = np.reshape(sqw[:, 4], (nq, nw))[:, :]
            self.results["sqw_t"] = (
                self.results["sqw_x"] ** 2 + self.results["sqw_y"] ** 2
            )

        # Read S(q,w) intensity
        if os.path.isfile(f"sqwintensity.{simid}.out"):
            sqw_t_int = np.genfromtxt(f"sqwintensity.{simid}.out", usecols=(0, 4, 5, 6))
            nq_int = int(sqw_t_int[-1, 0])
            nw_int = int(sqw_t_int.shape[0] / nq_int)
            self.results["sqw_int"] = np.reshape(sqw_t_int[:, 2], (nq_int, nw_int))

        # Read LSWT S(q,w) intensity
        if os.path.isfile(f"ncsqw_intensity.{simid}.out"):
            sqw_lt_int = np.genfromtxt(
                f"ncsqw_intensity.{simid}.out", usecols=(0, 4, 5, 6)
            )
            nq_lint = int(sqw_lt_int[-1, 0])
            nw_lint = int(sqw_lt_int.shape[0] / nq_lint)
            self.results["sqw_lint"] = np.reshape(sqw_lt_int[:, 2], (nq_lint, nw_lint))

        # Read S(q,w) tensor
        if os.path.isfile(f"sqwtensa.{simid}.out"):
            sqwt = np.genfromtxt(
                f"sqwtensa.{simid}.out", usecols=(0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)
            )
            nqt = int(sqwt[-1, 0])
            nwt = int(sqwt.shape[0] / nqt)
            sqw_tens = np.zeros((nqt, nwt, 3, 3))
            for i in range(3):
                for j in range(3):
                    sqw_tens[:, :, i, j] = np.reshape(
                        sqwt[:, 2 + i * 3 + j], (nqt, nwt)
                    )
            self.results["sqw_tens"] = sqw_tens

        # Read LSWT S(q,w) tensor
        if os.path.isfile(f"ncsqw.{simid}.out"):
            lswt_sqwt = np.genfromtxt(
                f"ncsqw.{simid}.out", usecols=(0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13)
            )
            lswt_nqt = int(lswt_sqwt[-1, 0])
            lswt_nwt = int(lswt_sqwt.shape[0] / lswt_nqt)
            lswt_sqw_tens = np.zeros((lswt_nqt, lswt_nwt, 3, 3))
            for i in range(3):
                for j in range(3):
                    lswt_sqw_tens[:, :, i, j] = np.reshape(
                        lswt_sqwt[:, 2 + i * 3 + j], (lswt_nqt, lswt_nwt)
                    )
            self.results["lswt_sqw_tens"] = lswt_sqw_tens

    def analyze_symmetry(self):
        """Analyze crystal symmetry and get Brillouin zone information."""
        lattice = self.inputs["lattice"]
        positions = self.inputs["positions"]
        numbers = self.inputs["numbers"]

        self.symmetry["cell"] = (lattice, positions, numbers)
        self.symmetry["spacegroup"] = get_spacegroup(self.symmetry["cell"])
        self.symmetry["BZ"], self.symmetry["sympoints"] = get_symmetry_points(
            self.symmetry["cell"]
        )

    def extract_symmetry_labels(self):
        """Extract symmetry point labels from qfile for plotting."""
        qfile = self.inputs.get("qfile")
        if qfile and "qfile.kpath" in qfile:
            with open(qfile, "r", encoding="utf-8") as f:
                f.readline()
                qpts = f.readlines()

            self.symmetry["axlab"] = []
            self.symmetry["axidx"] = []
            q_vecs = self.results.get("q_vecs")

            if q_vecs is not None:
                for idx, row in enumerate(qpts):
                    rs = row.split()
                    if len(rs) == 4:
                        label = rs[3]
                        self.symmetry["axlab"].append(
                            "$\\Gamma$" if label[0] == "G" else label
                        )
                        self.symmetry["axidx"].append(q_vecs[idx])

    def _calculate_emax(self):
        """Calculate maximum energy for plotting."""
        timestep = self.inputs["timestep"]
        sc_step = self.inputs["sc_step"]
        hbar = self.constants["hbar_eV_s"]
        return (
            0.5 * np.float64(hbar) / (np.float64(timestep) * np.float64(sc_step)) * 1e3
        )

    def _apply_gaussian_filter(self, data):
        """Apply Gaussian filtering to data for smoothing."""
        sigma_q = self.plot_config["sigma_q"]
        sigma_w = self.plot_config["sigma_w"]

        filtered = ndimage.gaussian_filter1d(
            data, sigma=sigma_q, axis=1, mode="constant"
        )
        filtered = ndimage.gaussian_filter1d(
            filtered, sigma=sigma_w, axis=0, mode="reflect"
        )
        return filtered

    def plot_magnon_spectra(self):
        """Plot adiabatic magnon spectra."""
        ams = self.results.get("ams")
        ams_pq = self.results.get("ams_pq")
        ams_mq = self.results.get("ams_mq")
        q_vecs = self.results.get("q_vecs")
        ams_dist_col = self.results.get("ams_dist_col")
        axidx = self.symmetry.get("axidx", [])
        axlab = self.symmetry.get("axlab", [])

        if ams is None or q_vecs is None or ams_dist_col is None:
            return

        # Basic AMS plot
        plt.figure(figsize=[8, 5])
        plt.plot(q_vecs[:], ams[:, 1:ams_dist_col])
        plt.xticks(axidx, axlab)
        plt.ylabel("Energy (meV)")
        plt.autoscale(tight=True)
        plt.ylim(0)
        plt.grid(visible=True, which="major", axis="x")
        plt.savefig("ams.png")
        plt.close()

        # Enhanced plot with +q and -q components if available
        if ams_pq is not None or ams_mq is not None:
            plt.figure(figsize=[8, 5])
            plt.plot(q_vecs[:], ams[:, 1:ams_dist_col], label="E(q)")
            if ams_pq is not None:
                plt.plot(q_vecs[:], ams_pq[:, 1:ams_dist_col], label="E(q+q$_0$)")
            if ams_mq is not None:
                plt.plot(q_vecs[:], ams_mq[:, 1:ams_dist_col], label="E(q-q$_0$)")
            plt.xticks(axidx, axlab)
            plt.legend()
            plt.ylabel("Energy (meV)")
            plt.autoscale(tight=True)
            plt.ylim(0)
            plt.grid(visible=True, which="major", axis="x")
            plt.savefig("ams_q.png")
            plt.close()



# In[ ]:

############################################################
# Plot the S(q,w)
############################################################
if got_sqw:
    fig = plt.figure(figsize=[8,5])
    ax=plt.subplot(111)
    
    
    hbar=4.135667662e-15
    emax=0.5*np.float64(hbar)/(np.float64(timestep)*np.float64(sc_step))*1e3
    sqw_temp=(sqw_x**2+sqw_y**2)**0.5
    sqw_temp[:,0]=sqw_temp[:,0]/100.0
    sqw_temp=sqw_temp.T/sqw_temp.T.max(axis=0)
    sqw_temp=ndimage.gaussian_filter1d(sqw_temp,sigma=sigma_q,axis=1,mode='constant')
    sqw_temp=ndimage.gaussian_filter1d(sqw_temp,sigma=sigma_w,axis=0,mode='reflect')
    #plt.imshow(sqw_temp, cmap=cmap.gist_ncar_r, interpolation='nearest',origin='lower',extent=[axidx_abs[0],axidx_abs[-1],0,emax])
    plt.imshow(sqw_temp, cmap=cmap.gist_ncar_r, interpolation='nearest',origin='lower',extent=[q_min, q_max,0,emax])
    plt.plot(q_vecs[:],ams[:,1:ams_dist_col],'black',lw=0.5)
    ala=plt.xticks()
    
    
    plt.xticks(axidx_abs,axlab)
    plt.xlabel('q')
    plt.ylabel('Energy (meV)')
    
    plt.grid(visible=True,which='major',axis='x')
    ax.set_aspect('auto')
    plt.autoscale(tight=True)
    #plt.show()
    plt.savefig('ams_sqw.png')


############################################################
# Plot the S(q,w) with full nc-LSWT support
############################################################
if got_sqw and got_ncams and got_ncams_pq and got_ncams_mq:
    fig = plt.figure(figsize=[8,5])
    ax=plt.subplot(111)
    
    
    hbar=4.135667662e-15
    emax=0.5*np.float64(hbar)/(np.float64(timestep)*np.float64(sc_step))*1e3
    sqw_temp=(sqw_x**2+sqw_y**2)**0.5
    sqw_temp[:,0]=sqw_temp[:,0]/100.0
    sqw_temp=sqw_temp.T/sqw_temp.T.max(axis=0)
    sqw_temp=ndimage.gaussian_filter1d(sqw_temp,sigma=sigma_q,axis=1,mode='constant')
    sqw_temp=ndimage.gaussian_filter1d(sqw_temp,sigma=sigma_w,axis=0,mode='reflect')
    imx_min=np.min(ams[:,0]/ams[-1,0]*axidx_abs[-1])
    imx_max=np.max(ams[:,0]/ams[-1,0]*axidx_abs[-1])
    #plt.imshow(sqw_temp, cmap=cmap.gist_ncar_r, interpolation='nearest',origin='lower',extent=[imx_min,imx_max,0,emax])
    plt.imshow(sqw_temp, cmap=cmap.gist_ncar_r, interpolation='nearest',origin='lower',extent=[q_min, q_max,0,emax])
    #plt.imshow(sqw_temp, cmap=cmap.gist_ncar_r, interpolation='nearest',origin='lower',extent=[axidx_abs[0],axidx_abs[-1],0,emax])
    #plt.imshow(sqw_temp, cmap=cmap.gist_ncar_r, interpolation='nearest',origin='lower',extent=[axidx_abs[0],axidx_abs[-1],0,emax])
    plt.plot(q_vecs[:],   ams[:,1:ams_dist_col],'black',lw=0.5)
    plt.plot(q_vecs[:],ams_pq[:,1:ams_dist_col],'black',lw=0.5)
    plt.plot(q_vecs[:],ams_mq[:,1:ams_dist_col],'black',lw=0.5)
    ala=plt.xticks()
    
    
    plt.xticks(axidx_abs,axlab)
    #print('::::>',[np.min(ams[:,0]/ams[-1,0]*axidx_abs[-1]), np.max(ams[:,0]/ams[-1,0]*axidx_abs[-1])])
    plt.xlabel('q')
    plt.ylabel('Energy (meV)')
    
    plt.grid(visible=True,which='major',axis='x')
    ax.set_aspect('auto')
    plt.autoscale(tight=True)
    #plt.xlim([np.min(ams[:,0]/ams[-1,0]*axidx_abs[-1]), np.max(ams[:,0]/ams[-1,0]*axidx_abs[-1])])
    plt.savefig('ams_sqw_q.png')



############################################################
# Plot the simulated S(q,w) intensity
############################################################
if got_sqw_int:
    fig = plt.figure(figsize=[8,5])
    ax=plt.subplot(111)
    
    
    
    hbar=4.135667662e-15
    emax=0.5*np.float64(hbar)/(np.float64(timestep)*np.float64(sc_step))*1e3
    sqw_temp=sqw_int
    sqw_temp[:,0]=sqw_temp[:,0]/100.0
    sqw_temp=sqw_temp.T/sqw_temp.T.max(axis=0)
    sqw_temp=ndimage.gaussian_filter1d(sqw_temp,sigma=sigma_q,axis=1,mode='constant')
    sqw_temp=ndimage.gaussian_filter1d(sqw_temp,sigma=sigma_w,axis=0,mode='reflect')
    plt.imshow(sqw_temp, cmap=cmap.gist_ncar_r, interpolation='nearest',origin='lower',extent=[axidx_abs[0],axidx_abs[-1],0,emax])
    plt.plot(ams[:,0]/ams[-1,0]*axidx_abs[-1],ams[:,1:ams_dist_col],'black',lw=1)
    ala=plt.xticks()
    
    
    plt.xticks(axidx_abs,axlab)
    plt.xlabel('q')
    plt.ylabel('Energy (meV)')
    
    plt.autoscale(tight=False)
    ax.set_aspect('auto')
    plt.grid(visible=True,which='major',axis='x')
    #plt.show()
    plt.savefig('ams_sqw_int.png')



############################################################
# Plot the LSWT S(q,w) intensity
############################################################
if got_sqw_lint:
    fig = plt.figure(figsize=[8,5])
    ax=plt.subplot(111)
    
    
    
    hbar=4.135667662e-15
    emax=0.5*np.float64(hbar)/(np.float64(timestep)*np.float64(sc_step))*1e3
    sqw_temp=sqw_lint
    sqw_temp[:,0]=sqw_temp[:,0]/100.0
    sqw_temp=sqw_temp.T/sqw_temp.T.max(axis=0)
    sqw_temp=ndimage.gaussian_filter1d(sqw_temp,sigma=sigma_q,axis=1,mode='constant')
    sqw_temp=ndimage.gaussian_filter1d(sqw_temp,sigma=sigma_w,axis=0,mode='reflect')
    plt.imshow(sqw_temp, cmap=cmap.gist_ncar_r, interpolation='nearest',origin='lower',extent=[axidx_abs[0],axidx_abs[-1],0,emax_lswt])
    plt.plot(ams[:,0]/ams[-1,0]*axidx_abs[-1],ams[:,1:ams_dist_col],'black',lw=1)
    ala=plt.xticks()
    
    
    plt.xticks(axidx_abs,axlab)
    plt.xlabel('q')
    plt.ylabel('Energy (meV)')
    
    plt.autoscale(tight=False)
    ax.set_aspect('auto')
    plt.grid(visible=True,which='major',axis='x')
    #plt.show()
    plt.savefig('ams_ncsqw_int.png')




xyz=('x','y','z')

############################################################
# Plot the S(q,w) on tensorial form (diagonal terms)
############################################################
if got_sqw_tens:
    fig = plt.figure(figsize=[16,4])
    hbar=4.135667662e-15
    emax=0.5*np.float64(hbar)/(np.float64(timestep)*np.float64(sc_step))*1e3

    plt.xticks(axidx_abs,axlab)
    plt.xlabel('q')
    plt.ylabel('Energy (meV)')
    
    #plt.autoscale(tight=False)
    ax.set_aspect('auto')
    plt.grid(visible=True,which='major',axis='x')

    #sqw_x[:,0]=sqw_x[:,0]/100.0
    #sqw_x=sqw_x.T/sqw_x.T.max(axis=0)
    plt_idx=130
    for ix in range(3):
        iy=ix
        sqw_temp=sqw_tens[:,:,ix,iy]
        sqw_temp[:,0]=sqw_temp[:,0]/1e5
        sqw_temp=sqw_temp.T/sqw_temp.T.max(axis=0)
        plt_idx=plt_idx+1
        ax=plt.subplot(plt_idx)
        stitle = f'$S^{{{xyz[ix] + xyz[iy]}}}(q,\\omega)$'
        ax.set_title(stitle)

        sqw_temp=ndimage.gaussian_filter1d(sqw_temp,sigma=sigma_q,axis=1,mode='constant')
        sqw_temp=ndimage.gaussian_filter1d(sqw_temp,sigma=sigma_w,axis=0,mode='reflect')
        plt.imshow(sqw_temp, cmap=cmap.gist_ncar_r, interpolation='nearest',origin='lower',extent=[axidx_abs[0],axidx_abs[-1],0,emax],aspect='auto')
        plt.plot(ams[:,0]/ams[-1,0]*axidx_abs[-1],ams[:,1:ams_dist_col],'black',lw=  1)
        plt.xticks(axidx_abs,axlab)

    plt.tight_layout()

    plt.savefig('sqw_diagonal.png')



############################################################
# Plot the S(q,w) on tensorial form (full tensor)
############################################################
if got_sqw_tens:
    fig = plt.figure(figsize=[16,10])
    hbar=4.135667662e-15
    emax=0.5*np.float64(hbar)/(np.float64(timestep)*np.float64(sc_step))*1e3

    plt.xticks(axidx_abs,axlab)
    plt.xlabel('q')
    plt.ylabel('Energy (meV)')
    
    #plt.autoscale(tight=False)
    ax.set_aspect('auto')
    plt.grid(visible=True,which='major',axis='x')

    #sqw_x[:,0]=sqw_x[:,0]/100.0
    #sqw_x=sqw_x.T/sqw_x.T.max(axis=0)
    plt_idx=330
    for ix in range(3):
        for iy in range(3):
            sqw_temp=sqw_tens[:,:,ix,iy]
            sqw_temp[:,0]=sqw_temp[:,0]/1e5
            sqw_temp=sqw_temp.T/sqw_temp.T.max(axis=0)
            plt_idx=plt_idx+1
            ax=plt.subplot(plt_idx)
            stitle = f'$S^{{{xyz[ix] + xyz[iy]}}}(q,\\omega)$'
            ax.set_title(stitle)

            sqw_temp=ndimage.gaussian_filter1d(sqw_temp,sigma=sigma_q,axis=1,mode='constant')
            sqw_temp=ndimage.gaussian_filter1d(sqw_temp,sigma=sigma_w,axis=0,mode='reflect')
            plt.imshow(sqw_temp, cmap=cmap.gist_ncar_r, interpolation='nearest',origin='lower',extent=[axidx_abs[0],axidx_abs[-1],0,emax],aspect='auto')
            plt.plot(ams[:,0]/ams[-1,0]*axidx_abs[-1],ams[:,1:ams_dist_col],'black',lw=  1)
            plt.xticks(axidx_abs,axlab)

    plt.tight_layout()

    plt.savefig('sqw_tensor.png')



############################################################
# Plot the LSWT S(q,w) on tensorial form (only diagonal terms)
############################################################
if got_lswt_sqw_tens:
    fig = plt.figure(figsize=[16,4])
    hbar=4.135667662e-15
    emax=0.5*np.float64(hbar)/(np.float64(timestep)*np.float64(sc_step))*1e3

    plt.xticks(axidx_abs,axlab)
    plt.xlabel('q')
    plt.ylabel('Energy (meV)')
    
    #plt.autoscale(tight=False)
    #ax.set_aspect('auto')
    plt.grid(visible=True,which='major',axis='x')

    #sqw_x[:,0]=sqw_x[:,0]/100.0
    #sqw_x=sqw_x.T/sqw_x.T.max(axis=0)
    plt_idx=130
    for ix in range(3):
        iy=ix
        sqw_temp=lswt_sqw_tens[:,:,ix,iy]
        sqw_temp=sqw_temp.T/(sqw_temp.T.max(axis=0)+1.0e-20)
        plt_idx=plt_idx+1
        ax=plt.subplot(plt_idx)
        stitle = f'$S^{{{xyz[ix] + xyz[iy]}}}_{{LSWT}}$'
        ax.set_title(stitle)

        sqw_temp=ndimage.gaussian_filter1d(sqw_temp,sigma=sigma_q,axis=1,mode='constant')
        sqw_temp=ndimage.gaussian_filter1d(sqw_temp,sigma=sigma_w,axis=0,mode='reflect')
        plt.imshow(sqw_temp, cmap=cmap.gist_ncar_r, interpolation='nearest',origin='lower',extent=[axidx_abs[0],axidx_abs[-1],0,emax_lswt],aspect='auto')
        #plt.plot(ams[:,0]/ams[-1,0]*axidx_abs[-1],ams[:,1:ams_dist_col],'black',lw=  1)
        plt.xticks(axidx_abs,axlab)

    plt.tight_layout()

    plt.savefig('ncsqw_diagonal.png')


############################################################
# Plot the LSWT S(q,w) on tensorial form (full tensor)
############################################################
if got_lswt_sqw_tens:
    fig = plt.figure(figsize=[16,10])
    hbar=4.135667662e-15
    emax=0.5*np.float64(hbar)/(np.float64(timestep)*np.float64(sc_step))*1e3

    plt.xticks(axidx_abs,axlab)
    plt.xlabel('q')
    plt.ylabel('Energy (meV)')
    
    #plt.autoscale(tight=False)
    #ax.set_aspect('auto')
    plt.grid(visible=True,which='major',axis='x')

    #sqw_x[:,0]=sqw_x[:,0]/100.0
    #sqw_x=sqw_x.T/sqw_x.T.max(axis=0)
    plt_idx=330
    for ix in range(3):
        for iy in range(3):
            sqw_temp=lswt_sqw_tens[:,:,ix,iy]
            sqw_temp=sqw_temp.T/(sqw_temp.T.max(axis=0)+1.0e-20)
            plt_idx=plt_idx+1
            ax=plt.subplot(plt_idx)
            stitle = f'$S^{{{xyz[ix] + xyz[iy]}}}_{{LSWT}}$'
            ax.set_title(stitle)

                        sqw_temp = self._apply_gaussian_filter(sqw_temp)
                        plt.imshow(
                            sqw_temp,
                            cmap=cmap.gist_ncar_r,
                            interpolation="nearest",
                            origin="lower",
                            extent=[axidx[0], axidx[-1], 0, emax_lswt],
                            aspect="auto",
                        )
                        plt.xticks(axidx, axlab)
                plt.tight_layout()
                plt.savefig("ncsqw_tensor.png")
                plt.close()

    def run(self):
        """Execute the complete post-processing workflow."""
        print("Reading input parameters...")
        self.read_inputs()

        print("Reading magnon spectra...")
        self.read_magnon_spectra()

        print("Reading structure factors...")
        self.read_structure_factors()

        print("Analyzing symmetry...")
        self.analyze_symmetry()

        print("Extracting symmetry labels...")
        self.extract_symmetry_labels()

        print("Plotting magnon spectra...")
        self.plot_magnon_spectra()

        print("Plotting S(q,w)...")
        self.plot_sqw()

        print("Plotting S(q,w) intensity...")
        self.plot_sqw_intensity()

        print("Plotting S(q,w) tensor (diagonal)...")
        self.plot_sqw_tensor(diagonal_only=True)

        print("Plotting S(q,w) tensor (full)...")
        self.plot_sqw_tensor(diagonal_only=False)

        print("Post-processing complete!")


if __name__ == "__main__":
    # Execute post-processing
    processor = PostProcessor()
    processor.run()
