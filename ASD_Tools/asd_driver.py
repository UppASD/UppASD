"""
ASD Workflow Driver

This script orchestrates the full ASD workflow:
1. Structure input (from CIF or spglib tuple)
2. Symmetry analysis (spglib)
3. Neighbor shell extraction (pymatgen)
4. J assignment (interactive or config)
5. Export J to UppASD format
6. Setup input files (inpsd.dat, posfile, momfile)
7. Preprocess (preQ.py)
8. Run simulation (TBD)
9. Postprocess (postQ.py)
"""
import subprocess
import sys

from asd_io import (
    extract_structure_from_cif,
    get_symmetry_reduced_neighbor_shells,
    write_inpsd_dat,
    write_momfile,
    write_posfile,
)


# Step 1: Structure input (from CIF)
def step1_structure_input(cif_file, primitive=False):
    structure, name, a_lat = extract_structure_from_cif(cif_file, primitive=primitive)
    return structure, name, a_lat

# Step 2: Symmetry info (spglib)
def step2_symmetry_info(structure):
    import spglib
    spacegroup = spglib.get_spacegroup(structure)
    return spacegroup

# Step 3: Neighbor shells (pymatgen)
def step3_neighbor_shells(structure, cutoff=5.0):
    shells = get_symmetry_reduced_neighbor_shells(structure, cutoff=cutoff)
    return shells

# Step 4: Assign J (interactive, for now just assign J1=1.0, J2=0.0, ...)

# Exchange mapping class for shell assignment and export
class ASDExchangeMap:
    def __init__(self, shells):
        self.shells = shells  # {distance: [(i, j), ...]}
        self.shell_list = sorted(shells.keys())
        # Default: J1=1.0, J2=0.5, ...
        self.J_dict = {dist: 1.0 / (idx + 1) for idx, dist in enumerate(self.shell_list)}

    def print_shells(self):
        print("\nNeighbor shells:")
        for idx, dist in enumerate(self.shell_list):
            pairs = self.shells[dist]
            # Pick first pair as representative vector
            i, j = pairs[0]
            vec = f"({i} -> {j})"
            print(f"  Shell {idx+1}: distance={dist} Å, representative pair {vec}, n_pairs={len(pairs)}")

    def set_exchange_coupling(self, shell_idx, value):
        """
        Set the exchange coupling J for a given shell index (1-based).
        """
        if 1 <= shell_idx <= len(self.shell_list):
            dist = self.shell_list[shell_idx - 1]
            self.J_dict[dist] = value
        else:
            raise IndexError(f"Shell index {shell_idx} out of range.")

    def prompt_for_J(self):
        """
        Interactively prompt the user to assign J for each shell (CLI use).
        """
        print("\nAssign exchange J for each shell (press Enter for default value)")
        for idx, dist in enumerate(self.shell_list):
            default_J = round(1.0 / (idx + 1), 3)
            prompt = f"  J for shell {idx+1} at {dist} Å (default {default_J}): "
            val = input(prompt)
            try:
                J = float(val) if val.strip() else default_J
            except Exception:
                print(f"  Invalid input, using default {default_J}")
                J = default_J
            self.J_dict[dist] = J
            
    def get_J_dict(self):
        return self.J_dict

    def print_J(self):
        print("Assigned J values:")
        for idx, dist in enumerate(self.shell_list):
            print(f"  Shell {idx+1} (distance {dist} Å): J = {self.J_dict[dist]}")

# Step 5: Export J to UppASD format
def step5_export_J(shells, J_dict, filename="jfile"):
    """
    Export all neighbor pairs to jfile: i, j, r_ij(x,y,z), J_ij, |r_ij| (not symmetry reduced).
    i, j are 1-based indices.
    """
    with open(filename, "w", encoding="utf-8") as f:
        for dist in sorted(shells.keys()):
            J = J_dict[dist]
            for pair in shells[dist]:
                if len(pair) == 3:
                    i, j, vec = pair
                else:
                    i, j = pair
                    # reconstruct r_ij vector from structure (not available here, fallback to zeros)
                    vec = (0.0, 0.0, 0.0)
                x, y, z = vec
                f.write(f"{i+1:7d} {j+1:7d} {x:11.6f} {y:11.6f} {z:11.6f} {J:13.6f} {dist:11.6f}\n")

# Step 6: Setup input files (reuse uppasd_init)
def step6_setup_inputs(name, structure, a_lat):
    cell_matrix, positions, numbers = structure
    write_inpsd_dat(
        "inpsd.dat",
        name,
        (10, 10, 10),
        ("P", "P", "P"),
        cell_matrix,
        a_lat,
        "./posfile",
        "./momfile",
        "./jfile",
        "./qfile.kpath",
    )
    write_posfile("posfile", [(i+1, numbers[i], *positions[i]) for i in range(len(numbers))])
    write_momfile("momfile", [(i+1, numbers[i], 1.0, 0.0, 0.0, 1.0) for i in range(len(numbers))])

# Step 7: Preprocess (call preQ.py)
def step7_preprocess():
    subprocess.run([sys.executable, "preQ.py"])

# Step 8: Run simulation (TBD)
def step8_run_simulation():
    try:
        import uppasd.pyasd as asd
        print("Running UppASD simulation via asd.runuppasd()...")
        asd.runuppasd()
        print("UppASD simulation finished.")
    except ImportError:
        print("uppasd.pyasd module not found. Please ensure UppASD Python interface is installed.")
    except Exception as e:
        print(f"Error running UppASD simulation: {e}")

# Step 9: Postprocess (call postQ.py)
def step9_postprocess():
    subprocess.run([sys.executable, "postQ.py"])


import argparse


def main():
    parser = argparse.ArgumentParser(description="ASD Workflow Driver")
    parser.add_argument("cif_file", help="Input CIF file")
    parser.add_argument("--primitive", dest="primitive", action="store_true", default=True,
                        help="Use primitive cell (default: True)")
    parser.add_argument("--conventional", dest="primitive", action="store_false",
                        help="Use conventional cell instead of primitive")
    args = parser.parse_args()
    cif_file = args.cif_file
    primitive = args.primitive

    # 1. Structure input
    structure, name, a_lat = step1_structure_input(cif_file, primitive)
    cell_matrix, positions, numbers = structure
    cell_type = "primitive" if primitive else "conventional"
    print(f"Cell type: {cell_type}, number of sites: {len(numbers)}")
    # 2. Symmetry info
    spacegroup = step2_symmetry_info(structure)
    print(f"Spacegroup: {spacegroup}")
    # 3. Neighbor shells
    shells = step3_neighbor_shells(structure)
    print("Neighbor shells:")
    for dist, pairs in sorted(shells.items()):
        print(f"  Shell at {dist} Å: {pairs}")
    # 4. Assign J using ASDExchangeMap
    exchange_map = ASDExchangeMap(shells)
    exchange_map.print_shells()
    # Example: set J for shell 1 to 2.0 (user can do this in notebook/script)
    # exchange_map.set_exchange_coupling(1, 2.0)
    exchange_map.print_J()
    J_dict = exchange_map.get_J_dict()

    # 5. Export J
    # Prepare shells_with_vec for jfile export
    shells_with_vec = {}
    for dist, pairs in shells.items():
        new_pairs = []
        for pair in pairs:
            if len(pair) == 3:
                new_pairs.append(pair)
            else:
                i, j = pair
                cell_matrix, positions, numbers = structure
                r_i = positions[i]
                r_j = positions[j]
                vec = tuple(rj - ri for ri, rj in zip(r_i, r_j))
                new_pairs.append((i, j, vec))
        shells_with_vec[dist] = new_pairs
    step5_export_J(shells_with_vec, J_dict)

    # 6. Setup input files
    step6_setup_inputs(name, structure, a_lat)

    # 7. Preprocess using PreProcessor class
    try:
        from preQ import PreProcessor
        pre = PreProcessor()
        pre.run()
    except ImportError:
        print("preQ.py with PreProcessor class not found, skipping preprocessing.")

    # 8. Run simulation (TBD)
    step8_run_simulation()

    # 9. Postprocess using PostProcessor class
    try:
        from postQ import PostProcessor
        post = PostProcessor()
        post.run()
    except ImportError:
        print("postQ.py with PostProcessor class not found, skipping postprocessing.")

if __name__ == "__main__":
    main()
