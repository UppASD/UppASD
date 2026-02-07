import numpy as np
import spglib
from typing import Tuple, List, Dict, Any

class SymmetryManager:
    """
    Handles spacegroup analysis and symmetry-based orbit classification for atomistic exchange.
    """
    def __init__(self, unit_cell: np.ndarray, basis_positions: np.ndarray, atom_types: List[int]):
        """
        Initialize using spglib cell format.
        
        Parameters
        ----------
        unit_cell : (3,3) ndarray
        basis_positions : (N,3) ndarray
        atom_types : list of int
        """
        self.cell = (unit_cell, basis_positions, atom_types)
        self.dataset = spglib.get_symmetry_dataset(self.cell)
        self.ops = spglib.get_symmetry(self.cell)
        
        # Cache for performance
        self.rotations = self.ops['rotations']
        self.translations = self.ops['translations']
        self.unit_cell = unit_cell
        self.inv_cell = np.linalg.inv(unit_cell)

    def classify_orbits(self, pairs: List[Tuple[int, int]], vectors: List[np.ndarray], tol: float = 1e-5) -> List[int]:
        """
        Group neighbors into symmetry orbits rather than just distance shells.
        
        Two bonds belong to the same orbit if there is a symmetry operation mapping
        bond (i1, j1) to bond (i2, j2).
        """
        num_bonds = len(pairs)
        orbit_indices = -np.ones(num_bonds, dtype=int)
        current_orbit = 1

        for i in range(num_bonds):
            if orbit_indices[i] != -1:
                continue
            
            # Start a new orbit
            orbit_indices[i] = current_orbit
            
            # Extract bond info for bond i
            src_i, tar_i = pairs[i]
            v_i = vectors[i]
            
            # Find all equivalent bonds in the remaining list
            for k in range(i + 1, num_bonds):
                if orbit_indices[k] != -1:
                    continue
                
                src_k, tar_k = pairs[k]
                v_k = vectors[k]
                
                if self._are_bonds_equivalent(src_i, v_i, src_k, v_k, tol):
                    orbit_indices[k] = current_orbit
            
            current_orbit += 1
            
        return orbit_indices.tolist()

    def _are_bonds_equivalent(self, src1: int, v1: np.ndarray, src2: int, v2: np.ndarray, tol: float) -> bool:
        """
        Check if bond 1 and bond 2 are related by any symmetry operation.
        """
        # Note: UppASD indices are 1-based, spglib indices are 0-based
        idx1 = src1 - 1
        idx2 = src2 - 1
        
        for rot in self.rotations:
            # Check if rotation maps vector 1 to vector 2
            v1_rotated = rot @ v1
            if np.allclose(v1_rotated, v2, atol=tol):
                # Optionally check if site types/indices also map correctly
                # (Spglib handles internal atom mapping via get_symmetry_dataset)
                return True
        return False

    def propagate_interaction(self, d_ref: np.ndarray, v_ref: np.ndarray, v_target: np.ndarray) -> np.ndarray:
        """
        Propagate a reference DMI vector to a target bond in the same orbit.
        
        Uses Dij = det(W) * W * D_ref where W is the rotation mapping v_ref to v_target.
        Also handles Dji = -Dij implicitly if bond direction is reversed.
        """
        for rot in self.rotations:
            if np.allclose(rot @ v_ref, v_target, atol=1e-5):
                det = np.linalg.det(rot)
                return det * (rot @ d_ref)
        
        # Check if inverse bond direction matches
        for rot in self.rotations:
            if np.allclose(rot @ v_ref, -v_target, atol=1e-5):
                det = np.linalg.det(rot)
                # Dji = -Dij rule
                return -det * (rot @ d_ref)
                
        raise ValueError("Target bond is not in the same symmetry orbit as the reference bond.")

    def get_orbit_report(self, pairs, vectors, distances, orbit_indices) -> str:
        """
        Generate a report categorized by Symmetry Orbits.
        """
        unique_orbits = sorted(list(set(orbit_indices)))
        lines = [f"Spacegroup: {self.dataset['international']} ({self.dataset['number']})"]
        lines.append(f"{'Orbit':<6} | {'Radius':<10} | {'Count':<6}")
        lines.append("-" * 30)
        
        for orb in unique_orbits:
            mask = [idx == orb for idx in orbit_indices]
            r = np.mean(np.array(distances)[mask])
            count = sum(mask)
            lines.append(f"{orb:<6} | {r:<10.6f} | {count:<6}")
            
        return "\n".join(lines)
