/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 * SPDX-FileCopyrightText: 2026 UppASD contributors
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Reference model of UppASD multiscaleInterpolateInterfaces (src/model.f90).
 *
 * Updates emom, emom2 and emomM in place on the host in double precision.
 * Index arrays are 1-based. Arrays are column-major with Fortran shapes
 * emom(3, atom_count, ensemble_count) and mmom(atom_count, ensemble_count).
 *
 * @param atom_count      Number of atoms.
 * @param ensemble_count  Number of ensembles.
 * @param row_count       Number of interpolated atoms.
 * @param weight_count    Number of weights and neighbours.
 * @param indices         Row of each atom, 0 if not interpolated. atom_count elements.
 * @param first_neighbour Row offsets into weights and neighbours. row_count + 1 elements.
 * @param neighbours      Neighbour atom of each weight. weight_count elements.
 * @param weights         Interpolation weights. weight_count elements.
 * @param mmom            Moment magnitudes. atom_count * ensemble_count elements.
 * @param emom            Unit moments, updated in place. 3 * atom_count * ensemble_count elements.
 * @param emom2           Unit moments, updated in place. 3 * atom_count * ensemble_count elements.
 * @param emomM           Moments, updated only at interpolated atoms. 3 * atom_count * ensemble_count elements.
 */
void model_multiscale_interpolate_interfaces(int atom_count, int ensemble_count, int row_count,
                                             int weight_count, const int* indices,
                                             const int* first_neighbour, const int* neighbours,
                                             const double* weights, const double* mmom,
                                             double* emom, double* emom2, double* emomM);

#ifdef __cplusplus
}
#endif
