/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include "testcase.h"

/**
 * Allocates an uninitialized test case in device memory.
 *
 * @param atom_count     Number of atoms.
 * @param ensemble_count Number of ensembles.
 * @param row_count      Number of interpolated atoms.
 * @param weight_count   Number of weights and neighbours.
 * @return               The test case. Free with testcase_device_destroy.
 */
testcase_t testcase_device_allocate(const size_t atom_count, const size_t ensemble_count,
                                    const size_t row_count, const size_t weight_count);

/**
 * Copies a host test case into newly allocated device memory.
 *
 * @param host Host test case.
 * @return     The device copy. Free with testcase_device_destroy.
 */
testcase_t testcase_device_clone_h2d(const testcase_t host);

/**
 * Copies a device test case into newly allocated host memory.
 *
 * @param device Device test case.
 * @return       The host copy. Free with testcase_host_destroy.
 */
testcase_t testcase_host_clone_d2h(const testcase_t device);

/**
 * Frees a device test case and zeroes it.
 *
 * @param testcase Device test case.
 */
void testcase_device_destroy(testcase_t* testcase);
