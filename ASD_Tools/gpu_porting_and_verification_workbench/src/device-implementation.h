/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 * SPDX-FileCopyrightText: 2026 UppASD contributors
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * Interpolates the moments of multiscale interface atoms on the device.
 *
 * Out-of-place counterpart of UppASD multiscaleInterpolateInterfaces. Each
 * interpolated atom receives normalized(sum_i weights[i] * emom_in[neighbours[i]]),
 * and likewise for emom2. Non-interpolated atoms are copied from the inputs.
 * Neighbours are read only from the inputs, also when they are interpolated, as in
 * the reference model src/model.f90.
 *
 * Index arrays are 1-based and use the LocalInterpolationInfo layout. Moment
 * arrays are column-major with Fortran shapes emom(3, atom_count, ensemble_count)
 * and mmom(atom_count, ensemble_count). All pointers are device pointers.
 * Outputs must not alias inputs.
 *
 * @param atom_count      Number of atoms.
 * @param ensemble_count  Number of ensembles.
 * @param indices         Row of each atom, 0 if not interpolated. atom_count elements.
 * @param first_neighbour Row offsets into weights and neighbours. Number of rows + 1 elements.
 * @param neighbours      Neighbour atom of each weight.
 * @param weights         Interpolation weights.
 * @param mmom            Moment magnitudes. atom_count * ensemble_count elements.
 * @param emom_in         Unit moments. 3 * atom_count * ensemble_count elements.
 * @param emom2_in        Unit moments. 3 * atom_count * ensemble_count elements.
 * @param emom_out        Interpolated emom_in. 3 * atom_count * ensemble_count elements.
 * @param emom2_out       Interpolated emom2_in. 3 * atom_count * ensemble_count elements.
 * @param emomM_out       emom_out * mmom, written only at interpolated atoms.
 *                        3 * atom_count * ensemble_count elements.
 */
void device_multiscale_interpolate_interfaces(const size_t atom_count, const size_t ensemble_count,
                                              const int* indices, const int* first_neighbour,
                                              const int* neighbours, const double* weights,
                                              const double* mmom, const double* emom_in,
                                              const double* emom2_in, double* emom_out,
                                              double* emom2_out, double* emomM_out);
#ifdef __cplusplus
}
#endif
