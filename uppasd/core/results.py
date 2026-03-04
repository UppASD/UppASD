"""
results.py

Lightweight container for UppASD simulation results.

Design principles:
- One canonical internal storage: self._tables
- Zero intelligence: no physics, no assumptions
- Explicit named accessors for common result types
- Introspection-friendly (__repr__, available, summary, describe)
"""

from pathlib import Path
from typing import Optional, Dict, Any
import numpy as np

from .output import read_all_outputs


class ASDResults:
    """
    Container for parsed UppASD output files.

    Parameters
    ----------
    workdir : str or Path
        Directory containing UppASD output files.
    """

    # ------------------------------------------------------------------
    # Construction
    # ------------------------------------------------------------------

    def __init__(self, workdir):
        self.workdir = Path(workdir)

        # Canonical internal storage for all parsed tables
        self._tables: Dict[str, Any] = {}

        tables, simid = read_all_outputs(self.workdir)
        self._tables.update(tables)
        self.simid = simid

    # ------------------------------------------------------------------
    # Introspection basics
    # ------------------------------------------------------------------

    @property
    def available(self):
        """
        Names of available result tables.
        """
        return tuple(self._tables.keys())

    def __repr__(self):
        tables = ", ".join(self.available) if self.available else "none"
        return (
            f"<ASDResults(workdir='{self.workdir}', "
            f"simid='{self.simid}', tables=[{tables}])>"
        )

    # ------------------------------------------------------------------
    # Dict-style access
    # ------------------------------------------------------------------

    def get(self, key):
        """
        Retrieve a parsed observable table by name.

        Returns None if not present.
        """
        return self._tables.get(key)

    def __getitem__(self, key):
        """
        Dict-style access: results["totenergy"]
        """
        try:
            return self._tables[key]
        except KeyError:
            raise KeyError(
                f"Observable '{key}' not found. "
                f"Available: {tuple(self._tables.keys())}"
            )

    # ------------------------------------------------------------------
    # Averaged / global observables
    # ------------------------------------------------------------------

    @property
    def averages(self) -> Optional[dict]:
        """Time-averaged magnetization and statistics."""
        return self._tables.get("averages")

    @property
    def cumulants(self) -> Optional[dict]:
        """Cumulants time series parsed from cumulants.<simid>.out."""
        return self._tables.get("cumulants")

    @property
    def thermodynamics(self) -> Optional[dict]:
        """Final thermodynamic summary parsed from cumulants.<simid>.json."""
        return self._tables.get("thermodynamics")

    @property
    def totenergy(self) -> Optional[dict]:
        """Total and decomposed energies vs iteration."""
        return self._tables.get("totenergy")

    # ------------------------------------------------------------------
    # Local / snapshot observables
    # ------------------------------------------------------------------

    @property
    def moment(self) -> Optional[dict]:
        """Local magnetic moments (multiple snapshots)."""
        return self._tables.get("moment")

    @property
    def befftot(self) -> Optional[dict]:
        """Local effective magnetic fields."""
        return self._tables.get("befftot")

    @property
    def restart(self) -> Optional[dict]:
        """Final spin configuration."""
        return self._tables.get("restart")

    @property
    def coord(self) -> Optional[dict]:
        """Expanded coordinates of the replicated system."""
        return self._tables.get("coord")

    @property
    def sq(self) -> Optional[dict]:
        """Static structure factor table parsed from `sq.<simid>.out`.

        Keys follow `read_sq` normalization: `iq, qx, qy, qz, qw, resxx, resyy, reszz, abs`.
        """
        return self._tables.get("sq")

    # ------------------------------------------------------------------
    # Scalar convenience accessors (safe, optional)
    # ------------------------------------------------------------------

    @property
    def final_energy(self):
        """
        Final total energy (last sampled step).
        """
        if self.totenergy is None:
            return None
        # tolerate different column naming conventions
        for key in ("Tot", "tot"):
            if key in self.totenergy:
                return self.totenergy[key][-1]
        return None

    @property
    def final_magnitude(self):
        """
        Final average magnetization magnitude.
        """
        if self.averages is None:
            return None
        for key in ("<M>", "m", "M"):
            if key in self.averages:
                return self.averages[key][-1]
        return None

    def final_thermo(self) -> Optional[dict]:
        """
        Return final thermodynamic observables as a compact dictionary.

        Preference order:
        1. ``thermodynamics`` table (from cumulants.<simid>.json)
        2. Last row of ``cumulants`` table (from cumulants.<simid>.out)

        Returns
        -------
        dict or None
            Keys: ``m``, ``binder``, ``chi``, ``cv``, ``energy``.
            Returns None if neither source is available.
        """
        if self.thermodynamics is not None:
            thermo = self.thermodynamics
            return {
                "m": float(thermo["m"]),
                "binder": float(thermo["binder"]),
                "chi": float(thermo["chi"]),
                "cv": float(thermo["cv"]),
                "energy": float(thermo["energy"]),
            }

        if self.cumulants is not None:
            cumulants = self.cumulants
            return {
                "m": float(cumulants["m"][-1]),
                "binder": float(cumulants["binder"][-1]),
                "chi": float(cumulants["chi"][-1]),
                "cv": float(cumulants["cv"][-1]),
                "energy": float(cumulants["energy"][-1]),
            }

        return None

    # ------------------------------------------------------------------
    # Didactic / user-facing inspection helpers
    # ------------------------------------------------------------------

    def summary(self):
        """
        Print a concise, didactic overview of available results.
        """
        print("UppASD results summary")
        print("-" * 32)
        print(f"Workdir : {self.workdir}")
        print(f"SimID   : {self.simid}")
        print()

        for key in self.available:
            table = self._tables[key]
            shape = self._infer_shape(table)
            print(f"• {key:<15} {shape}")

    def describe(self, verbose: bool = False):
        """
        Describe all available observables in more detail.

        Parameters
        ----------
        verbose : bool
            If True, also prints shapes of individual fields.
        """
        print("UppASD results description")
        print("-" * 36)

        for key in self.available:
            self.help(key, verbose=verbose)
            print()

    def help(self, key: str, verbose: bool = False):
        """
        Describe a single observable table.

        Parameters
        ----------
        key : str
            Observable name
        verbose : bool
            If True, print detailed field shapes
        """
        table = self[key]

        print(f"{key}")
        print("-" * len(key))

        kind = "dict" if isinstance(table, dict) else "array"
        print(f"Type  : {kind}")
        print(f"Shape : {self._infer_shape(table)}")

        if isinstance(table, dict):
            fields = list(table.keys())
            print(f"Fields: {fields}")

            if verbose:
                for name, val in table.items():
                    if hasattr(val, "shape"):
                        print(f"  {name:<10} shape={val.shape}")

    # ------------------------------------------------------------------
    # Reloading (interactive / long runs)
    # ------------------------------------------------------------------

    def reload(self):
        """
        Re-read output files from disk.
        """
        self._tables.clear()
        tables, simid = read_all_outputs(self.workdir)
        self._tables.update(tables)
        self.simid = simid

    # ------------------------------------------------------------------
    # Convenience converters
    # ------------------------------------------------------------------

    def coord_to_np(self) -> np.ndarray:
        """Return positions from `coord` as an (Ncoord, 3) numpy array.

        Raises
        -----
        RuntimeError: if `coord` table is not available.
        """
        coord = self.coord
        if coord is None:
            raise RuntimeError("coord table not available in results")

        return np.vstack((coord["x"], coord["y"], coord["z"])).T

    def restart_to_np(self, fallback_by_iatom: bool = True):
        """Return (pos_aligned, vecs) arrays aligned to the `restart` ordering.

        - `pos_aligned` (Nrestart, 3): positions taken from `coord` matched to
          each `restart` row by (replica, iatom). If no exact match is found
          and `fallback_by_iatom` is True, a row with the same `iatom` is used.
        - `vecs` (Nrestart, 3): reconstructed moment vectors = (mx,my,mz) * magnitude.

        Raises
        -----
        RuntimeError or KeyError on missing data or unrecoverable mapping failures.
        """
        coord = self.coord
        restart = self.restart
        if coord is None or restart is None:
            raise RuntimeError("coord and restart tables are required")

        # Build exact mapping (replica, iatom) -> index in coord arrays
        coord_pairs = [(int(r), int(i)) for r, i in zip(coord["replica"], coord["iatom"]) ]
        coord_map = {pair: idx for idx, pair in enumerate(coord_pairs)}

        def find_coord_index(replica, iatom):
            key = (int(replica), int(iatom))
            if key in coord_map:
                return coord_map[key]
            if not fallback_by_iatom:
                raise KeyError(f"No coord entry for replica/iatom {key}")

            matches = np.where(coord["iatom"] == int(iatom))[0]
            if matches.size == 1:
                return matches[0]
            if matches.size > 1:
                # prefer same replica if present among matches
                for m in matches:
                    if int(coord["replica"][m]) == int(replica):
                        return m
                return matches[0]

            raise KeyError(
                f"No coord entry for replica/iatom {key}; available sample: {coord_pairs[:20]}"
            )

        n = len(restart["iatom"])
        pos = np.empty((n, 3), dtype=float)
        vecs = np.empty((n, 3), dtype=float)

        for i in range(n):
            repl = restart["replica"][i]
            iatom = restart["iatom"][i]
            idx = find_coord_index(repl, iatom)
            pos[i, :] = (coord["x"][idx], coord["y"][idx], coord["z"][idx])
            mag = restart["magnitude"][i]
            vecs[i, :] = np.array((restart["mx"][i], restart["my"][i], restart["mz"][i])) * mag

        return pos, vecs

    # ------------------------------------------------------------------
    # Internals
    # ------------------------------------------------------------------

    @staticmethod
    def _infer_shape(table):
        """
        Infer a human-readable shape description.
        """
        if table is None:
            return "(none)"

        if isinstance(table, dict):
            # common UppASD conventions
            if {"mx", "my", "mz"} <= table.keys():
                return "(Nt × Natom × 3)"
            if {"bx", "by", "bz"} <= table.keys():
                return "(Nt × Natom × 3)"
            if {"rx", "ry"} <= table.keys():
                return "(Nt × 2)"
            if {"iter"} <= table.keys():
                return "(Nt, …)"
            return "(table)"

        if hasattr(table, "shape"):
            return str(table.shape)

        return "(unknown)"
