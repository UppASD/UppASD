#!/usr/bin/env python
# coding: utf-8

"""
UppASD Post-Processing Module (refactored)

Provides PostProcessor class that safely reads UppASD outputs and
creates plots. All top-level work is encapsulated to avoid import-time
execution and to ensure missing files / missing inputs are handled
gracefully.

Refactor goals:
- safe defaults for missing inpsd keys
- presence flags for optional files
- no import-time plotting
- encapsulated plotting with defensive checks
"""

import os
import os.path
import numpy as np
import spglib as spg
import seekpath as spth
import matplotlib.cm as cmap
import matplotlib.pyplot as plt
from scipy import ndimage


def is_close(A, B, atol=1e-6):
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)
    return bool(np.allclose(A, B, atol=atol))


class PostProcessor:
    """Safe post-processor for UppASD results."""

    def __init__(self, input_file="inpsd.dat"):
        self.input_file = input_file

        # Inputs and results
        self.inputs = {}
        self.results = {}

        # plotting / defaults
        self.plot_config = {
            "sigma_q": 0.5,
            "sigma_w": 1.0,
            "font": {"family": "sans", "weight": "normal", "size": "14"},
        }
        plt.rc("font", **self.plot_config["font"])
        plt.rc("lines", lw=2)

        # constants
        self.constants = {
            "hbar": 6.582119514e-13,
            "ry_ev": 13.605693009,
            "hbar_eV_s": 4.135667662e-15,
        }

    def read_posfile(self, posfile):
        positions = np.empty((0, 3), dtype=float)
        numbers = []
        try:
            with open(posfile, "r", encoding="utf-8") as pfile:
                for line in pfile:
                    parts = line.split()
                    if len(parts) >= 5:
                        numbers.append(int(parts[1]))
                        positions = np.vstack((positions, np.asarray(parts[2:5], dtype=float)))
        except Exception:
            return None, None
        return positions, np.asarray(numbers, dtype=int)

    def read_inpsd(self):
        # conservative defaults
        defaults = {
            "lattice": None,
            "positions": None,
            "numbers": None,
            "simid": "0",
            "mesh": [1, 1, 1],
            "posfiletype": "C",
            "timestep": 1.0,
            "sc_step": 1,
            "sc_nstep": 1,
            "qfile": None,
        }

        data = {}
        try:
            with open(self.input_file, "r", encoding="utf-8") as f:
                lines = f.readlines()
        except Exception:
            self.inputs = defaults
            return

        lattice = None
        positions = None
        numbers = None
        simid = defaults["simid"]
        mesh = defaults["mesh"]
        timestep = defaults["timestep"]
        sc_step = defaults["sc_step"]
        sc_nstep = defaults["sc_nstep"]
        qfile = defaults["qfile"]

        for idx, line in enumerate(lines):
            parts = line.split()
            if not parts:
                continue
            key = parts[0]
            if key == "simid":
                simid = parts[1]
            elif key == "cell":
                # read 3 lines including this
                try:
                    l0 = np.asarray(lines[idx].split()[1:4], dtype=float)
                    l1 = np.asarray(lines[idx + 1].split()[0:3], dtype=float)
                    l2 = np.asarray(lines[idx + 2].split()[0:3], dtype=float)
                    lattice = np.vstack((l0, l1, l2))
                except Exception:
                    lattice = None
            elif key == "ncell":
                try:
                    n = int(parts[1])
                    mesh = [n, n, n]
                except Exception:
                    pass
            elif key == "timestep":
                try:
                    timestep = float(parts[1].replace("d", "e"))
                except Exception:
                    pass
            elif key == "sc_step":
                try:
                    sc_step = int(parts[1])
                except Exception:
                    pass
            elif key == "sc_nstep":
                try:
                    sc_nstep = int(parts[1])
                except Exception:
                    pass
            elif key == "qfile":
                qfile = parts[1]
            elif key.strip() == "posfile":
                posfile = parts[1]
                p, n = self.read_posfile(posfile)
                if p is not None:
                    positions = p
                    numbers = n

        self.inputs = {
            "lattice": lattice,
            "positions": positions,
            "numbers": numbers,
            "simid": simid,
            "mesh": mesh,
            "posfiletype": defaults["posfiletype"],
            "timestep": timestep,
            "sc_step": sc_step,
            "sc_nstep": sc_nstep,
            "qfile": qfile,
        }

    def read_magnon_spectra(self):
        simid = str(self.inputs.get("simid", "0"))
        got_ams = os.path.isfile(f"ams.{simid}.out")
        got_ncams = os.path.isfile(f"ncams.{simid}.out")
        got_ncams_pq = os.path.isfile(f"ncams+q.{simid}.out")
        got_ncams_mq = os.path.isfile(f"ncams-q.{simid}.out")

        self.results["got_ams"] = got_ams or got_ncams
        if got_ams:
            try:
                ams = np.loadtxt(f"ams.{simid}.out")
            except Exception:
                ams = None
        elif got_ncams:
            try:
                ams = np.loadtxt(f"ncams.{simid}.out")
            except Exception:
                ams = None
        else:
            ams = None

        self.results["ams"] = ams
        if ams is not None:
            if ams.ndim == 1:
                ams = ams.reshape(1, -1)
            self.results["ams_dist_col"] = int(ams.shape[1] - 1)
            self.results["q_vecs"] = ams[:, 0]
            self.results["q_min"] = float(np.min(ams[:, 0]))
            self.results["q_max"] = float(np.max(ams[:, 0]))
            try:
                self.results["emax_lswt"] = 1.10 * np.amax(ams[:, 1 : self.results["ams_dist_col"]])
            except Exception:
                self.results["emax_lswt"] = None
        else:
            self.results.setdefault("ams_dist_col", None)
            self.results.setdefault("q_vecs", None)
            self.results.setdefault("q_min", 0.0)
            self.results.setdefault("q_max", 1.0)
            self.results.setdefault("emax_lswt", None)

        if got_ncams_pq:
            try:
                self.results["ams_pq"] = np.loadtxt(f"ncams+q.{simid}.out")
            except Exception:
                self.results["ams_pq"] = None
        else:
            self.results["ams_pq"] = None

        if got_ncams_mq:
            try:
                self.results["ams_mq"] = np.loadtxt(f"ncams-q.{simid}.out")
            except Exception:
                self.results["ams_mq"] = None
        else:
            self.results["ams_mq"] = None

    def read_structure_factors(self):
        simid = str(self.inputs.get("simid", "0"))
        # sqw
        got_sqw = os.path.isfile(f"sqw.{simid}.out")
        self.results["got_sqw"] = got_sqw
        if got_sqw:
            try:
                sqw = np.genfromtxt(f"sqw.{simid}.out", usecols=(0, 4, 5, 6, 7, 8))
                nq = int(sqw[-1, 0])
                nw = int(sqw.shape[0] / nq)
                self.results["sqw_x"] = np.reshape(sqw[:, 2], (nq, nw))
                self.results["sqw_y"] = np.reshape(sqw[:, 3], (nq, nw))
                self.results["sqw_z"] = np.reshape(sqw[:, 4], (nq, nw))
                self.results["sqw_t"] = self.results["sqw_x"] ** 2 + self.results["sqw_y"] ** 2
            except Exception:
                self.results["got_sqw"] = False

        # sqw intensity
        got_sqw_int = os.path.isfile(f"sqwintensity.{simid}.out")
        self.results["got_sqw_int"] = got_sqw_int
        if got_sqw_int:
            try:
                sqw_t_int = np.genfromtxt(f"sqwintensity.{simid}.out", usecols=(0, 4, 5, 6))
                nq_int = int(sqw_t_int[-1, 0])
                nw_int = int(sqw_t_int.shape[0] / nq_int)
                self.results["sqw_int"] = np.reshape(sqw_t_int[:, 2], (nq_int, nw_int))
            except Exception:
                self.results["got_sqw_int"] = False

        # sqw tensor
        got_sqw_tens = os.path.isfile(f"sqwtensa.{simid}.out")
        self.results["got_sqw_tens"] = got_sqw_tens
        if got_sqw_tens:
            try:
                sqwt = np.genfromtxt(f"sqwtensa.{simid}.out", usecols=(0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13))
                nqt = int(sqwt[-1, 0])
                nwt = int(sqwt.shape[0] / nqt)
                tens = np.zeros((nqt, nwt, 3, 3))
                for i in range(3):
                    for j in range(3):
                        tens[:, :, i, j] = np.reshape(sqwt[:, 2 + i * 3 + j], (nqt, nwt))
                self.results["sqw_tens"] = tens
            except Exception:
                self.results["got_sqw_tens"] = False

        # LSWT tensors
        got_lswt = os.path.isfile(f"ncsqw.{simid}.out")
        self.results["got_lswt_sqw"] = got_lswt
        if got_lswt:
            try:
                lswt_sqwt = np.genfromtxt(f"ncsqw.{simid}.out", usecols=(0, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13))
                nqt = int(lswt_sqwt[-1, 0])
                nwt = int(lswt_sqwt.shape[0] / nqt)
                tens = np.zeros((nqt, nwt, 3, 3))
                for i in range(3):
                    for j in range(3):
                        tens[:, :, i, j] = np.reshape(lswt_sqwt[:, 2 + i * 3 + j], (nqt, nwt))
                self.results["lswt_sqw_tens"] = tens
            except Exception:
                self.results["got_lswt_sqw"] = False

    def analyze_symmetry(self):
        lattice = self.inputs.get("lattice")
        positions = self.inputs.get("positions")
        numbers = self.inputs.get("numbers")
        if lattice is None or positions is None or numbers is None:
            self.symmetry = {"cell": None, "spacegroup": None, "BZ": None, "sympoints": {}}
            return

        cell = (lattice, positions, numbers)
        try:
            sg = spg.get_spacegroup(cell)
            kpath_obj = spth.get_path(cell)
            BZ = np.asarray(kpath_obj.get("reciprocal_primitive_lattice", []))
            sympoints = kpath_obj.get("point_coords", {})
        except Exception:
            sg = None
            BZ = None
            sympoints = {}

        self.symmetry = {"cell": cell, "spacegroup": sg, "BZ": BZ, "sympoints": sympoints}

    def extract_symmetry_labels(self):
        qfile = self.inputs.get("qfile")
        q_vecs = self.results.get("q_vecs")
        self.symmetry.setdefault("axlab", [])
        self.symmetry.setdefault("axidx", [])
        if not qfile or q_vecs is None:
            return

        def _format_label(label):
            if label is None:
                return label
            s = str(label).strip()
            if not s:
                return s
            low = s.lower()
            if low == "g" or low == "gamma" or s == "Γ":
                return "$\\Gamma$"
            # sometimes seekpath returns strings like 'GAMMA' or 'Gamma'
            if low.startswith("gamma"):
                return "$\\Gamma$"
            return s

        try:
            if "qfile.kpath" in qfile:
                with open(qfile, "r", encoding="utf-8") as f:
                    f.readline()
                    qpts = f.readlines()
                for idx, row in enumerate(qpts):
                    rs = row.split()
                    if len(rs) >= 4:
                        label = _format_label(rs[3])
                        self.symmetry["axlab"].append(label)
                        self.symmetry["axidx"].append(q_vecs[idx])
            else:
                qpts = np.genfromtxt(qfile, skip_header=1, usecols=(0, 1, 2))
                for idx, row in enumerate(qpts):
                    for k, v in self.symmetry.get("sympoints", {}).items():
                        if is_close(v, row):
                            label = _format_label(k)
                            self.symmetry["axlab"].append(label)
                            self.symmetry["axidx"].append(q_vecs[idx])
                            break
        except Exception:
            # leave axlab/axidx as-is
            pass

    def _calculate_emax(self):
        timestep = float(self.inputs.get("timestep", 1.0))
        sc_step = float(self.inputs.get("sc_step", 1))
        hbar = self.constants.get("hbar_eV_s", 4.135667662e-15)
        return 0.5 * hbar / (timestep * sc_step) * 1e3

    def plot_magnon_spectra(self):
        ams = self.results.get("ams")
        q_vecs = self.results.get("q_vecs")
        ams_dist_col = self.results.get("ams_dist_col")
        if ams is None or q_vecs is None or ams_dist_col is None:
            print("AMS data incomplete; skipping magnon spectra plots.")
            return

        axlab = self.symmetry.get("axlab", [])
        axidx = self.symmetry.get("axidx", [])

        plt.figure(figsize=[8, 5])
        plt.plot(q_vecs, ams[:, 1:ams_dist_col])
        try:
            plt.xticks(axidx, axlab)
        except Exception:
            pass
        plt.ylabel("Energy (meV)")
        plt.autoscale(tight=True)
        plt.ylim(0)
        plt.grid(visible=True, which="major", axis="x")
        plt.savefig("ams.png")
        plt.close()

        # enhanced plot
        ams_pq = self.results.get("ams_pq")
        ams_mq = self.results.get("ams_mq")
        if ams_pq is not None or ams_mq is not None:
            plt.figure(figsize=[8, 5])
            plt.plot(q_vecs, ams[:, 1:ams_dist_col], label="E(q)")
            if ams_pq is not None:
                try:
                    plt.plot(q_vecs, ams_pq[:, 1:ams_dist_col], label="E(q+q0)")
                except Exception:
                    pass
            if ams_mq is not None:
                try:
                    plt.plot(q_vecs, ams_mq[:, 1:ams_dist_col], label="E(q-q0)")
                except Exception:
                    pass
            try:
                plt.xticks(axidx, axlab)
            except Exception:
                pass
            plt.legend()
            plt.ylabel("Energy (meV)")
            plt.autoscale(tight=True)
            plt.ylim(0)
            plt.grid(visible=True, which="major", axis="x")
            plt.savefig("ams_q.png")
            plt.close()

    def plot_sqw(self):
        if not self.results.get("got_sqw"):
            print("No S(q,w) data available; skipping plot_sqw.")
            return
        sqw_x = self.results.get("sqw_x")
        sqw_y = self.results.get("sqw_y")
        if sqw_x is None or sqw_y is None:
            print("S(q,w) arrays missing; skipping image plot.")
            return

        q_min = self.results.get("q_min", 0.0)
        q_max = self.results.get("q_max", 1.0)
        emax = self._calculate_emax()

        sqw_temp = np.sqrt(sqw_x ** 2 + sqw_y ** 2)
        try:
            sqw_temp[:, 0] = sqw_temp[:, 0] / 100.0
        except Exception:
            pass
        sqw_temp = sqw_temp.T
        maxcol = sqw_temp.max(axis=0) if sqw_temp.size else None
        if maxcol is not None and np.any(maxcol > 0):
            with np.errstate(invalid="ignore"):
                sqw_temp = sqw_temp / (maxcol + 1e-20)

        sqw_temp = ndimage.gaussian_filter1d(sqw_temp, sigma=self.plot_config["sigma_q"], axis=1, mode="constant")
        sqw_temp = ndimage.gaussian_filter1d(sqw_temp, sigma=self.plot_config["sigma_w"], axis=0, mode="reflect")

        fig = plt.figure(figsize=[8, 5])
        ax = plt.subplot(111)
        plt.imshow(sqw_temp, cmap=cmap.gist_ncar_r, interpolation="nearest", origin="lower", extent=[q_min, q_max, 0, emax])
        # overlay ams
        ams = self.results.get("ams")
        q_vecs = self.results.get("q_vecs")
        ams_dist_col = self.results.get("ams_dist_col")
        if ams is not None and q_vecs is not None and ams_dist_col is not None:
            try:
                plt.plot(q_vecs, ams[:, 1:ams_dist_col], "black", lw=0.5)
            except Exception:
                pass
        axlab = self.symmetry.get("axlab", [])
        axidx = self.symmetry.get("axidx", [])
        try:
            if axidx:
                plt.xticks(np.array(axidx, dtype=float), axlab)
        except Exception:
            pass
        plt.xlabel("q")
        plt.ylabel("Energy (meV)")
        plt.grid(visible=True, which="major", axis="x")
        ax.set_aspect("auto")
        plt.autoscale(tight=True)
        plt.savefig("ams_sqw.png")
        plt.close()

    def plot_sqw_intensity(self):
        if not self.results.get("got_sqw_int"):
            print("No S(q,w) intensity data; skipping plot_sqw_intensity.")
            return
        sqw_int = self.results.get("sqw_int")
        if sqw_int is None:
            print("S(q,w) intensity array missing; skipping.")
            return

        emax = self._calculate_emax()
        sqw_temp = sqw_int.copy()
        try:
            sqw_temp[:, 0] = sqw_temp[:, 0] / 100.0
        except Exception:
            pass
        sqw_temp = sqw_temp.T
        maxcol = sqw_temp.max(axis=0) if sqw_temp.size else None
        if maxcol is not None and np.any(maxcol > 0):
            sqw_temp = sqw_temp / (maxcol + 1e-20)

        sqw_temp = ndimage.gaussian_filter1d(sqw_temp, sigma=self.plot_config["sigma_q"], axis=1, mode="constant")
        sqw_temp = ndimage.gaussian_filter1d(sqw_temp, sigma=self.plot_config["sigma_w"], axis=0, mode="reflect")

        fig = plt.figure(figsize=[8, 5])
        ax = plt.subplot(111)
        axidx = self.symmetry.get("axidx", [])
        axlab = self.symmetry.get("axlab", [])
        axidx_abs = np.array(axidx, dtype=float) if axidx else np.array([0.0, 1.0])
        plt.imshow(sqw_temp, cmap=cmap.gist_ncar_r, interpolation="nearest", origin="lower", extent=[axidx_abs[0], axidx_abs[-1], 0, emax])
        plt.xticks(axidx_abs, axlab)
        plt.xlabel("q")
        plt.ylabel("Energy (meV)")
        plt.grid(visible=True, which="major", axis="x")
        plt.savefig("ams_sqw_int.png")
        plt.close()

    def plot_sqw_tensor(self, diagonal_only=True):
        if not self.results.get("got_sqw_tens"):
            print("No tensor S(q,w) data; skipping plot_sqw_tensor.")
            return
        tens = self.results.get("sqw_tens")
        if tens is None:
            print("Tensor data missing; skipping plot_sqw_tensor.")
            return

        emax = self._calculate_emax()
        axidx = self.symmetry.get("axidx", [])
        axlab = self.symmetry.get("axlab", [])
        axidx_abs = np.array(axidx, dtype=float) if axidx else np.array([0.0, 1.0])
        xyz = ("x", "y", "z")
        plt_idx = 130 if diagonal_only else 330
        fig = plt.figure(figsize=[16, 4 if diagonal_only else 10])
        for ix in range(3):
            for iy in ([ix] if diagonal_only else range(3)):
                sqw_temp = tens[:, :, ix, iy]
                try:
                    sqw_temp = sqw_temp.T / (sqw_temp.T.max(axis=0) + 1.0e-20)
                except Exception:
                    pass
                plt_idx += 1
                ax = plt.subplot(plt_idx)
                stitle = f'$S^{{{xyz[ix] + xyz[iy]}}}(q,\\omega)$'
                ax.set_title(stitle)
                sqw_temp = ndimage.gaussian_filter1d(sqw_temp, sigma=self.plot_config["sigma_q"], axis=1, mode="constant")
                sqw_temp = ndimage.gaussian_filter1d(sqw_temp, sigma=self.plot_config["sigma_w"], axis=0, mode="reflect")
                plt.imshow(sqw_temp, cmap=cmap.gist_ncar_r, interpolation="nearest", origin="lower", extent=[axidx_abs[0], axidx_abs[-1], 0, emax], aspect="auto")
                plt.xticks(axidx_abs, axlab)
        plt.tight_layout()
        outname = "sqw_diagonal.png" if diagonal_only else "sqw_tensor.png"
        plt.savefig(outname)
        plt.close()

    def run(self):
        print("Reading input parameters...")
        self.read_inpsd()
        print("Reading magnon spectra...")
        self.read_magnon_spectra()
        print("Reading structure factors...")
        self.read_structure_factors()
        print("Analyzing symmetry...")
        try:
            self.analyze_symmetry()
        except Exception:
            print("Warning: symmetry analysis failed; continuing.")
        print("Extracting symmetry labels...")
        try:
            self.extract_symmetry_labels()
        except Exception:
            print("Warning: extracting symmetry labels failed; continuing.")
        print("Plotting magnon spectra...")
        try:
            self.plot_magnon_spectra()
        except Exception:
            print("Warning: plot_magnon_spectra failed; continuing.")
        print("Plotting S(q,w)...")
        try:
            self.plot_sqw()
        except Exception:
            print("Warning: plot_sqw failed; continuing.")
        print("Plotting S(q,w) intensity...")
        try:
            self.plot_sqw_intensity()
        except Exception:
            print("Warning: plot_sqw_intensity failed; continuing.")
        print("Plotting S(q,w) tensor (diagonal)...")
        try:
            self.plot_sqw_tensor(diagonal_only=True)
        except Exception:
            print("Warning: plot_sqw_tensor (diagonal) failed; continuing.")
        print("Plotting S(q,w) tensor (full)...")
        try:
            self.plot_sqw_tensor(diagonal_only=False)
        except Exception:
            print("Warning: plot_sqw_tensor (full) failed; continuing.")
        # show what symmetry labels were detected so the user can verify
        try:
            print("Symmetry labels:", self.symmetry.get("axlab", []))
        except Exception:
            pass
        print("Post-processing complete!")


if __name__ == "__main__":
    p = PostProcessor()
    p.run()
