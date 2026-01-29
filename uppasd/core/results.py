"""
results.py

Thin container for UppASD simulation results.

Design goals:
- Zero intelligence
- Zero assumptions
- Explicit access to parsed outputs
- Convenience properties for common observables
"""

from pathlib import Path
from typing import Optional, Dict

from .output import read_outputs


class ASDResults:
    """
    Lightweight results container.

    Parameters
    ----------
    workdir : str or Path
        UppASD run directory
    """

    def __init__(self, workdir):
        self.workdir: Path = Path(workdir)
        self.outputs: Dict[str, Optional[dict]] = {}
        self.reload()

    # ------------------------------------------------------------------
    # Magnetization convenience
    # ------------------------------------------------------------------

    @property
    def magnetization(self) -> Optional[dict]:
        return self.outputs.get("magnetization")

    @property
    def step(self):
        m = self.magnetization
        return None if m is None else m["step"]

    @property
    def mx(self):
        m = self.magnetization
        return None if m is None else m["mx"]

    @property
    def my(self):
        m = self.magnetization
        return None if m is None else m["my"]

    @property
    def mz(self):
        m = self.magnetization
        return None if m is None else m["mz"]

    @property
    def m(self):
        m = self.magnetization
        return None if m is None else m["m"]

    # ------------------------------------------------------------------
    # Energy convenience
    # ------------------------------------------------------------------

    @property
    def energy(self) -> Optional[dict]:
        return self.outputs.get("energy")

    @property
    def total_energy(self):
        e = self.energy
        return None if e is None else e["energy"]

    # ------------------------------------------------------------------
    # Metadata / helpers
    # ------------------------------------------------------------------

    @property
    def available(self):
        """List available output blocks."""
        return [k for k, v in self.outputs.items() if v is not None]

    # ------------------------------------------------------------------
    # Refresh (useful for interactive runs)
    # ------------------------------------------------------------------

    def reload(self):
        """Re-read output files from disk."""
        self.outputs = read_outputs(self.workdir)

    # ------------------------------------------------------------------

    def __repr__(self):
        avail = ", ".join(self.available) or "none"
        return f"ASDResults(workdir='{self.workdir}', outputs=[{avail}])"
