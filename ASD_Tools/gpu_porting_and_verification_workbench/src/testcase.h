/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include <stddef.h>
#include <stdint.h>

#include "common/datatypes.h"

// The index arrays are passed to the toolkit as int32_t, to the implementations as int and to
// Fortran as integer(c_int)
static_assert(_Generic((int32_t)0, int: true, default: false), "int32_t must be int");
// The toolkit math and verification functions (normalize_real_t, verify_real_t) exist only for
// real_t
static_assert(_Generic((real_t)0, double: true, default: false), "real_t must be double");

/**
 * Inputs and outputs of multiscaleInterpolateInterfaces.
 *
 * Index arrays are 1-based and use the LocalInterpolationInfo layout. Moment arrays are
 * column-major with Fortran shapes emom(3, atom_count, ensemble_count) and
 * mmom(atom_count, ensemble_count). The function that creates a testcase_t determines
 * whether its arrays are in host or device memory.
 */
typedef struct {
    size_t   atom_count;      /**< Number of atoms. */
    size_t   ensemble_count;  /**< Number of ensembles. */
    size_t   row_count;       /**< Number of interpolated atoms. */
    size_t   weight_count;    /**< Number of weights and neighbours. */
    int32_t* indices;         /**< Row of each atom, 0 if not interpolated. atom_count elements. */
    int32_t* first_neighbour; /**< Row offsets into weights and neighbours. row_count + 1 elements. */
    int32_t* neighbours;      /**< Neighbour atom of each weight. weight_count elements. */
    double*  weights;         /**< Interpolation weights. weight_count elements. */
    double*  mmom;            /**< Moment magnitudes. atom_count * ensemble_count elements. */
    double*  emom;            /**< Unit moments. 3 * atom_count * ensemble_count elements. */
    double*  emom2;           /**< Unit moments. 3 * atom_count * ensemble_count elements. */
    double*  emomM;           /**< emom * mmom. 3 * atom_count * ensemble_count elements. */
} testcase_t;

/** Command-line options of the test programs. */
typedef struct {
    size_t atom_count;     /**< --atoms N. Default 100000. */
    size_t ensemble_count; /**< --ensembles N. Default 3. */
    bool   verify;         /**< Cleared by --no-verify. Default true. */
} testcase_options_t;

/**
 * Parses the command line. Prints the usage and exits on an invalid argument.
 *
 * @param argc Argument count, as passed to main.
 * @param argv Arguments, as passed to main.
 * @return     The options.
 */
testcase_options_t testcase_parse_options(const int argc, char* const argv[]);

/**
 * Allocates an uninitialized test case in host memory.
 *
 * @param atom_count     Number of atoms.
 * @param ensemble_count Number of ensembles.
 * @param row_count      Number of interpolated atoms.
 * @param weight_count   Number of weights and neighbours.
 * @return               The test case. Free with testcase_host_destroy.
 */
testcase_t testcase_host_allocate(const size_t atom_count, const size_t ensemble_count,
                                  const size_t row_count, const size_t weight_count);

/**
 * Creates a test case in host memory.
 *
 * Atoms follow the UppASD order: real atoms, padding atoms, finite-difference nodes, and a second
 * padding block. Every second node is interpolated from 1 to 8 real atoms. Padding atoms are
 * interpolated from 2, 4 or 8 nodes, including interpolated ones of higher index (first block)
 * and lower index (second block). Neighbours are sorted, distinct and non-consecutive within a
 * row. Weights are positive and irregular, and sum to 1 for padding atoms. Moments are distinct,
 * with positive components.
 *
 * @param atom_count     Number of atoms. At least 88.
 * @param ensemble_count Number of ensembles.
 * @return               The test case. Free with testcase_host_destroy.
 */
testcase_t testcase_host_create(const size_t atom_count, const size_t ensemble_count);

/**
 * Copies a host test case into newly allocated host memory.
 *
 * @param src Host test case.
 * @return    The copy. Free with testcase_host_destroy.
 */
testcase_t testcase_host_clone(const testcase_t src);

/**
 * Frees a host test case and zeroes it.
 *
 * @param testcase Host test case.
 */
void testcase_host_destroy(testcase_t* testcase);
