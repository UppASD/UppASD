"""
system.py

Definition of the SpinSystem class.

Responsibilities:
- Store lattice vectors, atomic positions, species, and magnetic moments
- Write UppASD structural input files (posfile, momfile)
- Perform basic consistency checks

NO physics logic.
NO defaults.
NO simulation state.
"""

from pathlib import Path
from typing import Sequence, List
import numpy as np


class SpinSystem:
    """
    Container for an atomistic spin system.
    Represents the UppASD *unit cell*.


    Parameters
    ----------
    cell : (3, 3) array-like
        Lattice vectors (Cartesian).
    positions : (N, 3) array-like
        Atomic positions (Cartesian).
    species : sequence of str
        Chemical species labels, length N.
    moments : (N, 3) array-like
        Magnetic moments for each atom.
    """

    def __init__(
        self,
        cell: Sequence[Sequence[float]],
        positions: Sequence[Sequence[float]],
        species: Sequence[str],
        moments: Sequence[Sequence[float]],
        ncell: Sequence[int] = None,
        bc: str = None,
    ):
        self.cell = np.asarray(cell, dtype=float)
        self.positions = np.asarray(positions, dtype=float)
        self.moments = np.asarray(moments, dtype=float)
        self.species = list(species)
        self.ncell = ncell
        self.bc = bc

        self._validate()

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    def _validate(self) -> None:
        """Perform basic consistency checks."""
        if self.cell.shape != (3, 3):
            raise ValueError("cell must be a (3, 3) array")

        natom = self.positions.shape[0]

        if self.positions.shape != (natom, 3):
            raise ValueError("positions must be an (N, 3) array")

        if self.moments.shape != (natom, 3):
            raise ValueError("moments must be an (N, 3) array")

        if len(self.species) != natom:
            raise ValueError(
                "species length must match number of atoms"
            )

    # ------------------------------------------------------------------
    # Introspection
    # ------------------------------------------------------------------

    @property
    def natom(self) -> int:
        """Number of atoms."""
        return self.positions.shape[0]

    # ------------------------------------------------------------------
    # Input block for inpsd.dat assembly
    # ------------------------------------------------------------------

    def input_block(self) -> dict:
        """
        Generate system block for inpsd.dat.
        
        Returns dict with keys in order: ncell, BC, cell, posfile, momfile.
        
        Note: ncell and BC have defaults (1 1 1 and P P P respectively)
        because Fortran requires them and they are fundamental structural parameters.
        """
        # Use ordered dict-like insertion to control key order
        block = {}
        
        # ncell: default to 1 1 1
        if self.ncell is not None:
            block["ncell"] = self.ncell
        else:
            block["ncell"] = [1, 1, 1]
        
        # BC: default to P P P (periodic in all directions)
        if self.bc is not None:
            block["BC"] = self.bc
        else:
            block["BC"] = "P P P"
        
        # cell
        if self.cell is not None:
            block["cell"] = self.cell
        
        # Files
        block["posfile"] = "./posfile"
        block["momfile"] = "./momfile"
        
        return block

    # ------------------------------------------------------------------
    # Writing UppASD input files
    # ------------------------------------------------------------------

    def write_posfile(self, filename) -> None:
        """
        Write UppASD posfile.

        Format:
            i_atom j_atom r_x r_y r_z
        
        Positions are written in lattice parameter units (fractional coordinates).
        """
        path = Path(filename)

        # Stable mapping: species label -> integer type (1-based)
        unique_species: List[str] = []
        for sp in self.species:
            if sp not in unique_species:
                unique_species.append(sp)

        type_map = {sp: i + 1 for i, sp in enumerate(unique_species)}

        # Convert Cartesian positions to fractional coordinates
        cell_inv = np.linalg.inv(self.cell).T
        fractional_positions = self.positions @ cell_inv

        with path.open("w") as f:
            # Atomic positions in fractional coordinates
            for i_atom, (frac_pos, sp) in enumerate(zip(fractional_positions, self.species), start=1):
                j_atom = type_map[sp]
                f.write(
                    f"{i_atom} {j_atom} {frac_pos[0]:.8f} {frac_pos[1]:.8f} {frac_pos[2]:.8f}\n"
                )

    # ------------------------------------------------------------------

    def write_momfile(self, filename) -> None:
        """
        Write UppASD momfile.

        Format:
            i_atom 1 m_mag m_x m_y m_z
        """
        path = Path(filename)

        with path.open("w") as f:
            for i_atom, m in enumerate(self.moments, start=1):
                m_mag = np.linalg.norm(m)
                f.write(
                    f"{i_atom} 1 {m_mag:.8f} {m[0]:.8f} {m[1]:.8f} {m[2]:.8f}\n"
                )

    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        species_set = sorted(set(self.species))
        return (
            f"SpinSystem(natom={self.natom}, "
            f"species={species_set})"
        )
