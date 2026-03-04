import numpy as np
from typing import Dict, Tuple, Sequence, List


class AnisotropyTable:
    """
    On-site anisotropy table for UppASD kfile format.

    Stores per-atom anisotropy entries and writes them in the kfile format
    expected by UppASD:

        atom_index ani_type K2 K4 e_x e_y e_z 0

    Conventions:
    - `ani_type` == 1 => uniaxial
    - negative `K2` indicates easy-axis anisotropy
    """

    def __init__(self):
        # map: iatom -> (ani_type, K2, K4, ex, ey, ez)
        self.sites: Dict[int, Tuple[int, float, float, float, float, float]] = {}

    def add_site(
        self,
        iatom: int,
        ani_type: int,
        K2: float,
        K4: float,
        direction: Sequence[float],
    ) -> None:
        """Add or replace anisotropy for a given site (1-based index).

        Args:
            iatom: 1-based atom index
            ani_type: anisotropy type (1 = uniaxial)
            K2: second-order anisotropy constant
            K4: fourth-order anisotropy constant
            direction: length-3 sequence giving anisotropy axis (ex,ey,ez)
        """
        if iatom < 1:
            raise ValueError("iatom must be >= 1 (1-based indexing)")
        if not isinstance(ani_type, int):
            raise TypeError("ani_type must be an integer")

        dir_arr = np.asarray(direction, dtype=float)
        if dir_arr.shape != (3,):
            raise ValueError("direction must be length-3")

        ex, ey, ez = float(dir_arr[0]), float(dir_arr[1]), float(dir_arr[2])
        self.sites[int(iatom)] = (int(ani_type), float(K2), float(K4), ex, ey, ez)

    def write_kfile(self, filename: str) -> None:
        """Write the anisotropy `kfile`.

        Each line follows: iatom ani_type K2 K4 ex ey ez 0
        """
        lines: List[str] = []
        for iatom in sorted(self.sites.keys()):
            ani_type, K2, K4, ex, ey, ez = self.sites[iatom]
            lines.append(
                f"{iatom:5d} {ani_type:1d} {K2:12.6f} {K4:12.6f} "
                f"{ex:8.3f} {ey:8.3f} {ez:8.3f} 0\n"
            )

        with open(filename, "w", encoding="utf-8") as f:
            f.writelines(lines)

    def to_arrays(self):
        i_list, type_list, K2_list, K4_list, dir_list = [], [], [], [], []
        for iatom in sorted(self.sites.keys()):
            ani_type, K2, K4, ex, ey, ez = self.sites[iatom]
            i_list.append(iatom)
            type_list.append(ani_type)
            K2_list.append(K2)
            K4_list.append(K4)
            dir_list.append(np.array([ex, ey, ez], dtype=float))

        return (
            np.array(i_list, dtype=int),
            np.array(type_list, dtype=int),
            np.array(K2_list, dtype=float),
            np.array(K4_list, dtype=float),
            np.array(dir_list, dtype=float),
        )

    def __repr__(self):
        return f"AnisotropyTable(n_sites={len(self.sites)})"
