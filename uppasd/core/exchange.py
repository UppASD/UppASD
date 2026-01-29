import numpy as np
from collections import defaultdict
from typing import Iterable, Tuple, Dict, List, Union


ArrayLike = Union[List[float], np.ndarray]


class ExchangeShellTable:
    """
    Exchange interactions grouped by (iatom, jatom, shell).

    Internal representation:
        shells[(iatom, jatom, shell)] = list of (R_vec, Jij)

    Notes
    -----
    - iatom, jatom are 1-based indices (UppASD convention)
    - No symmetry is enforced: (i, j) and (j, i) are distinct
    - shell indices are user-defined integers
    """

    def __init__(self):
        self.shells: Dict[
            Tuple[int, int, int],
            List[Tuple[np.ndarray, float]]
        ] = defaultdict(list)

    # ------------------------------------------------------------------
    # Basic API
    # ------------------------------------------------------------------

    def add_bond(
        self,
        iatom: int,
        jatom: int,
        shell: int,
        R_vec: ArrayLike,
        Jij: float,
    ):
        """Add a single exchange bond."""
        if iatom < 1 or jatom < 1:
            raise ValueError("iatom and jatom must be 1-based (>= 1)")
        if not isinstance(shell, int):
            raise TypeError("shell index must be an integer")

        R_vec = np.asarray(R_vec, dtype=float)
        if R_vec.shape != (3,):
            raise ValueError("R_vec must be length-3")

        self.shells[(iatom, jatom, shell)].append((R_vec, float(Jij)))

    def add_shell(
        self,
        iatom: int,
        jatom: int,
        shell: int,
        vectors: Iterable[ArrayLike],
        Jij: Union[float, Iterable[float]],
    ):
        """
        Add a full shell of equivalent bonds.

        Jij may be:
            - scalar (applied to all vectors)
            - list/array matching vectors
        """
        vectors = [np.asarray(v, dtype=float) for v in vectors]

        if np.isscalar(Jij):
            Jij_vals = [float(Jij)] * len(vectors)
        else:
            Jij_vals = list(Jij)
            if len(Jij_vals) != len(vectors):
                raise ValueError("Length of Jij list must match vectors")

        for R, J in zip(vectors, Jij_vals):
            self.add_bond(iatom, jatom, shell, R, J)

    # ------------------------------------------------------------------
    # Query & manipulation
    # ------------------------------------------------------------------

    def shells_for_pair(self, iatom: int, jatom: int):
        return sorted(
            s for (i, j, s) in self.shells.keys()
            if i == iatom and j == jatom
        )

    def scale_shell(
        self,
        shell: int,
        factor: float,
        iatom: int = None,
        jatom: int = None,
    ):
        for (i, j, s), bonds in self.shells.items():
            if s != shell:
                continue
            if iatom is not None and i != iatom:
                continue
            if jatom is not None and j != jatom:
                continue

            for idx, (R, J) in enumerate(bonds):
                bonds[idx] = (R, J * factor)

    def remove_shell(self, iatom: int, jatom: int, shell: int):
        self.shells.pop((iatom, jatom, shell), None)

    # ------------------------------------------------------------------
    # Construction helpers
    # ------------------------------------------------------------------

    @classmethod
    def from_neighbors(
        cls,
        pairs: Iterable[Tuple[int, int]],
        vectors: Iterable[ArrayLike],
        distances: Iterable[float],
        shell_indices: Iterable[int],
        Jij_model,
    ):
        table = cls()
        for (i, j), R, r, shell in zip(
            pairs, vectors, distances, shell_indices
        ):
            table.add_bond(i, j, shell, R, Jij_model(r))
        return table

    # ------------------------------------------------------------------
    # I/O
    # ------------------------------------------------------------------

    def write_jfile(self, filename: str):
        lines = []
        for key in sorted(self.shells):
            iatom, jatom, shell = key
            for R, Jij in self.shells[key]:
                Rx, Ry, Rz = R
                lines.append(
                    f"{iatom:5d} {jatom:5d} "
                    f"{Rx:12.6f} {Ry:12.6f} {Rz:12.6f} "
                    f"{Jij:16.8e}\n"
                )

        with open(filename, "w") as f:
            f.writelines(lines)

    # ------------------------------------------------------------------

    def to_arrays(self):
        i_list, j_list, R_list, J_list = [], [], [], []

        for key in sorted(self.shells):
            i, j, _ = key
            for R, J in self.shells[key]:
                i_list.append(i)
                j_list.append(j)
                R_list.append(R.copy())
                J_list.append(J)

        return (
            np.array(i_list, dtype=int),
            np.array(j_list, dtype=int),
            np.array(R_list, dtype=float),
            np.array(J_list, dtype=float),
        )

    def __repr__(self):
        return f"ExchangeShellTable(n_shells={len(self.shells)})"


# ======================================================================
# DMI
# ======================================================================


class DMIShellTable:
    """
    Dzyaloshinskii–Moriya interactions grouped by (iatom, jatom, shell).

    Same conventions as ExchangeShellTable.
    """

    def __init__(self):
        self.shells: Dict[
            Tuple[int, int, int],
            List[Tuple[np.ndarray, np.ndarray]]
        ] = defaultdict(list)

    def add_bond(
        self,
        iatom: int,
        jatom: int,
        shell: int,
        R_vec: ArrayLike,
        D_vec: ArrayLike,
    ):
        if iatom < 1 or jatom < 1:
            raise ValueError("iatom and jatom must be 1-based (>= 1)")
        if not isinstance(shell, int):
            raise TypeError("shell index must be an integer")

        R_vec = np.asarray(R_vec, dtype=float)
        D_vec = np.asarray(D_vec, dtype=float)

        if R_vec.shape != (3,) or D_vec.shape != (3,):
            raise ValueError("R_vec and D_vec must be length-3")

        self.shells[(iatom, jatom, shell)].append((R_vec, D_vec))

    def scale_shell(
        self,
        shell: int,
        factor: float,
        iatom: int = None,
        jatom: int = None,
    ):
        for (i, j, s), bonds in self.shells.items():
            if s != shell:
                continue
            if iatom is not None and i != iatom:
                continue
            if jatom is not None and j != jatom:
                continue

            for idx, (R, D) in enumerate(bonds):
                bonds[idx] = (R, D * factor)

    def write_dmfile(self, filename: str):
        lines = []
        for key in sorted(self.shells):
            iatom, jatom, shell = key
            for R, D in self.shells[key]:
                Rx, Ry, Rz = R
                Dx, Dy, Dz = D
                lines.append(
                    f"{iatom:5d} {jatom:5d} "
                    f"{Rx:12.6f} {Ry:12.6f} {Rz:12.6f} "
                    f"{Dx:12.6e} {Dy:12.6e} {Dz:12.6e}\n"
                )

        with open(filename, "w") as f:
            f.writelines(lines)

    def to_arrays(self):
        i_list, j_list, R_list, D_list = [], [], [], []

        for key in sorted(self.shells):
            i, j, _ = key
            for R, D in self.shells[key]:
                i_list.append(i)
                j_list.append(j)
                R_list.append(R.copy())
                D_list.append(D.copy())

        return (
            np.array(i_list, dtype=int),
            np.array(j_list, dtype=int),
            np.array(R_list, dtype=float),
            np.array(D_list, dtype=float),
        )

    def __repr__(self):
        return f"DMIShellTable(n_shells={len(self.shells)})"
