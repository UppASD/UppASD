"""Module for neighbour list creation and analysis"""
import numpy as np


def vector_in_list(i_vector, vector_list):
    """Checks if input vector is in list of vectors"""
    for j_vector in vector_list:
        if np.allclose(i_vector, j_vector, atol=1e-4, rtol=1e-3):
            return True
    return False

# %%


def get_full_nnlist(cell, c_atom, cutoff, in_cell_only=False):
    """ Routine for creating full neighbour list with cutoff (to later be symmetry reduced) """
    [lattice, positions, types] = cell
    l_scale = np.amax(lattice)
    l_range = np.int64(2*cutoff/l_scale+1.5)

    r_0 = np.matmul(lattice.T, positions[c_atom])
    shifts = np.mgrid[range(-l_range, l_range+1),
                      range(-l_range, l_range+1), range(-l_range, l_range+1)].reshape(3, -1).T
    r_list = []
    i_list = []
    t_list = []
    norm_list = []
    # asum_list = []
    for i_atom, coord in enumerate(positions):
        r_b = np.matmul(lattice.T, coord) - r_0
        # print('  ', i_atom, r_b)
        if in_cell_only:
            if np.linalg.norm(r_b) <= cutoff:
                r_list.append(r_b)
                i_list.append(i_atom)
                t_list.append(types[i_atom])
                norm_list.append(np.linalg.norm(r_b))
        else:
            vectors = np.matmul(shifts, lattice) + r_b
            filters = np.linalg.norm(
                vectors, axis=1).round(decimals=3) <= cutoff
            for hits in vectors[filters]:
                r_list.append(hits)
                i_list.append(i_atom)
                t_list.append(types[i_atom])
                norm_list.append(np.linalg.norm(hits))
                # asum_list.append(np.sum(np.abs(hits)))

    return np.array(r_list), np.array(t_list), np.array(norm_list)


# %%
def get_phonopy_nnlist(cell, c_atom, cutoff, symop_list):
    """Routine for creating neighbour list in phonopy fashion"""

    [lattice, positions, types] = cell
    r_0 = np.matmul(lattice.T, positions[c_atom])

    l_range = 1
    shifts = np.mgrid[range(-l_range, l_range+1),
                      range(-l_range, l_range+1), range(-l_range, l_range+1)].reshape(3, -1).T
    rfull_list = []
    ifull_list = []
    tfull_list = []
    nfull_list = []

    rsym_list = []
    isym_list = []
    tsym_list = []
    nsym_list = []

    for i_atom, coord in enumerate(positions):
        r_b = np.matmul(lattice.T, coord) - r_0
        vectors = np.matmul(shifts, lattice) + r_b
        norms = np.linalg.norm(vectors, axis=1)
        closest = np.argmin(norms)
        if norms[closest] < cutoff:
            vector = vectors[closest]
            r_b = vector
            rfull_list.append(r_b)
            ifull_list.append(i_atom)
            tfull_list.append(types[i_atom])
            nfull_list.append(np.linalg.norm(r_b))
            for rot, _ in symop_list:
                r_b = np.matmul(rot.T, vector)
                if not vector_in_list(r_b, rsym_list):
                    rsym_list.append(r_b)
                    isym_list.append(i_atom)
                    tsym_list.append(types[i_atom])
                    nsym_list.append(np.linalg.norm(r_b))

    print('Number of vectors:', len(rsym_list), np.amax(np.array(nsym_list)))
    return np.array(rsym_list), np.array(tsym_list), np.array(nsym_list), \
        np.amax(np.array(nsym_list))


# %%
def sort_vectors(in_arr):
    """Sort vectors according to magnitude (ascending)"""
    s_list = []
    for vector in in_arr:
        s_list.append(vector.tolist() + [np.linalg.norm(vector)])
    s_arr = np.array(s_list)
    s_arr = s_arr[s_arr[:, 3].argsort()]

    return s_arr[:, 0:3]


# %%
def sort_vectors_and_more(in_vecs, in_types, in_dists):
    """Routine to sort vectors according to magnitude (ascending)"""

    s_list = []
    for vector in in_vecs:
        s_list.append(np.linalg.norm(vector))
    s_arr = np.array(s_list).argsort()
    s_vecs = in_vecs[s_arr]
    s_types = in_types[s_arr]
    s_dists = in_dists[s_arr]

    return s_vecs, s_types, s_dists


# %%
# Routine to find nearest distance of a vector in a cell. Returns vector as a list.
# Works brute-force by making images of the cell in every direction and looking for
# the closest corner (r_0=[0 0 0])
# Assumes cubic lattice in direct coordinates for now, which should be ok for the application here.
# Notice the ambiguity with this approach where
# e.g. r=[0.5 0.5 0.5] is as close to r_0=[0 0 0] as r_0=[1 1 1]
# Assume cubic lattice in direct coordinates for now (should not matter)
def find_nearest_vector(i_vect):
    """Routine to find nearest distance of a vector in a cell. Returns vector as a list."""
    i_vect = np.array(i_vect)
    t_norm = np.linalg.norm(i_vect)
    o_vect = np.array(i_vect)
    for _x in range(-2, 3):
        for _y in range(-2, 3):
            for _z in range(-2, 3):
                t_vect = i_vect+np.array([_x, _y, _z])
                if np.linalg.norm(t_vect) <= t_norm:
                    t_norm = np.linalg.norm(t_vect)
                    o_vect = np.array(t_vect)

    return o_vect.tolist()

# %%


def get_shell_clust_list(types, dists):
    """Get list of clusters i.e. vectors with the same norm"""
    _, uniq_inv = np.unique(
        (dists+(types+1)*np.max(2*dists)).round(decimals=3), return_inverse=True)

    return np.max(uniq_inv), uniq_inv


# %%


def sym_identify_clusters(vecs, types, symops, sort=False, do_prn=False):
    """Finds symmetrically equivalent clusters for a full list of vectors"""

    if sort:
        sort_idx = np.linalg.norm(np.array(vecs), axis=1).argsort()
        vecs = np.array(vecs)[sort_idx]
        types = types[sort_idx]

    # vector = full list
    # Idea: loop through full list.
    # If vector in sym_list cycle; #
    # else if symop*vector not in sym_list:
    # add symop*vector to sym_list
    sym_vec_list = []
    sym_idx_list = []
    sym_norm_list = []
    sym_type_list = []
    uniq_vec_list = []
    for ivec, vec in enumerate(vecs):
        if not vector_in_list(vec, sym_vec_list):
            uniq_vec_list.append(vec)
            for rot, _ in symops:
                rot_vec = np.matmul(rot.T, vec)
                if not vector_in_list(rot_vec, sym_vec_list):
                    if vector_in_list(rot_vec, vecs):
                        sym_vec_list.append(rot_vec)
                        sym_idx_list.append(len(uniq_vec_list))
                        sym_norm_list.append(np.linalg.norm(rot_vec))
                        sym_type_list.append(types[ivec])

    if do_prn:
        for iclust, vector in enumerate(uniq_vec_list):
            print('------------------------------')
            print('Cluster         :', iclust+1)
            print('------------------------------')
            print('Number of sites :', np.sum(
                np.array(sym_idx_list, dtype=np.int32) == iclust+1))
            print('Vector          :', vector)
            print('Norm            :', np.linalg.norm(vector))

    clust_list = np.array(sym_idx_list) - 1

    return np.array(sym_vec_list), np.array(sym_type_list), clust_list
