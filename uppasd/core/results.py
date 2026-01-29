"""
results.py

Lightweight container for UppASD simulation results.

Design principles:
- One canonical internal storage: self._tables
- Zero intelligence: no physics, no assumptions
- Explicit named accessors for each result type
- Introspection-friendly (__repr__, available)
"""

from pathlib import Path
from typing import Optional, Dict, Any

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
    # Introspection
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
    # Averaged / global observables
    # ------------------------------------------------------------------

    @property
    def averages(self) -> Optional[dict]:
        """Time-averaged magnetization and statistics."""
        return self._tables.get("averages")

    @property
    def cumulants(self) -> Optional[dict]:
        """Binder cumulants, susceptibility, heat capacity, etc."""
        return self._tables.get("cumulants")

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

    # ------------------------------------------------------------------
    # Scalar convenience accessors (optional but useful)
    # ------------------------------------------------------------------

    @property
    def final_energy(self):
        """
        Final total energy (last sampled step).
        """
        if self.totenergy is None:
            return None
        return self.totenergy["Tot"][-1]

    @property
    def final_magnitude(self):
        """
        Final average magnetization magnitude.
        """
        if self.averages is None:
            return None
        return self.averages["<M>"][-1]

    # ------------------------------------------------------------------
    # Reloading (useful for interactive / long runs)
    # ------------------------------------------------------------------

    def reload(self):
        """
        Re-read output files from disk.
        """
        self._tables.clear()
        tables, simid = read_all_outputs(self.workdir)
        self._tables.update(tables)
        self.simid = simid

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