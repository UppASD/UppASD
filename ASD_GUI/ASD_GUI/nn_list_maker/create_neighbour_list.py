""" Create neighbour lists. Currently not stored"""
""" - Input: """
"""      inpsd.dat   : UppASD input file """
"""      posfile     : The posfile as defined in `inpsd.dat` """
"""      momfile     : The momfile as defined in `inpsd.dat` """
"""      cutoff      : Range for which to extract neighbours """
""" - Output: """
"""      nnvec_list  : list of neighbouring vectors """
"""      nntype_list : list of neighbour types """
"""      clust_list  : list of which vectors belong to which symmetry cluster (shell)"""

#from structure import get_full_nnlist, sym_identify_clusters
import numpy as np

def create_neighbour_list():
    from read_uppasd import get_uppasd_cell
    fname = './inpsd.dat'

    cell, sym_data = get_uppasd_cell(fname)

    return cell, sym_data

def get_symmetry_ops_from_dataset(sym_data):
    """ Exctract symmetry operations from standardized dataset"""
    """ These are the operations (assuming only r) that will be used in the following"""
    symop_list = [(r, t) for r, t in zip(
        sym_data['cart_rotations'], sym_data['cart_translations'])]

    return symop_list

def get_symmetry_data_from_structure(cell):
    import spglib as spg

    # Create spglib dataset from original cell (using primitive_matrix which is in units of a_lat)
    sym_data = spg.get_symmetry_dataset(cell)

    # Exctract the lattice vectors from the cell structure
    (lattice, _, _) = cell

    # %%
    # Augment spglib dataset with transformed rotations which work on primitive cells
    sym_data['cart_rotations'] = []
    sym_data['cart_translations'] = []
    inv_mat = np.linalg.inv(lattice.T)
    for rotation, translation in zip(sym_data['rotations'], sym_data['translations']):
        sym_data['cart_rotations'].append(
            np.matmul(np.matmul(lattice.T, rotation), inv_mat))
        sym_data['cart_translations'].append(
            np.matmul(lattice.T, translation))
            #np.matmul(np.matmul(lattice.T, translation), inv_mat))

    sym_data['cart_rotations'] = np.array(sym_data['cart_rotations'])
    sym_data['cart_translations'] = np.array(sym_data['cart_translations'])

    return sym_data

def reduce_vectors_from_symmetry(cell, nnvec_list):

    """
    Reduce list of neighbouring vectors due to symmetry of lattice structure.

    - Input: 
            cell            :   tuple containing basis vectors, positions 
                                and numbering as numpy arrays. 
            nnvec_list      :   list of neighbouring vectors
    - Output:
            reduced_vectors :   resulting numpy array with sym. reduced vectors
    """

    sym_data = get_symmetry_data_from_structure(cell)

    rotation_matrices = sym_data['cart_rotations']

    rotation_vectors = np.zeros((len(nnvec_list)*len(rotation_matrices), 3))

    index = 0
    for vector in nnvec_list:
        for matrix in rotation_matrices:
            rotation_vectors[index] = (matrix@(vector.reshape(3,1))).T
            index += 1

    norms = np.linalg.norm(rotation_vectors, axis = 1)
    same_norm_indices = np.unique(norms, return_index= True)[1]
    reduced_vectors = np.sort(rotation_vectors[same_norm_indices])*(-1)
    reduced_vectors = np.where(reduced_vectors == -0, 0, reduced_vectors)

    return reduced_vectors