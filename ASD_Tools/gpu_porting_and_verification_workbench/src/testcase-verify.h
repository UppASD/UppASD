/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#pragma once

#include "testcase.h"

/**
 * Applies the reference model (src/model.f90) in place to emom, emom2 and emomM.
 *
 * @param testcase Host test case.
 */
void testcase_run_model(testcase_t* testcase);

/**
 * Compares emom, emom2 and emomM element-wise with verify_real_t. An array passes if its maximum
 * relative error is at most 4 * DBL_EPSILON. Prints the first and last few interpolated moments
 * of the input, reference and candidate, then the outcome for each array with the maximum
 * absolute and relative errors and their positions. Indexing starts from 0.
 *
 * @param input    Host test case with the values before the call.
 * @param expected Host test case with reference values.
 * @param actual   Host test case with candidate values.
 * @return         true if all elements match.
 */
bool testcase_compare(const testcase_t input, const testcase_t expected, const testcase_t actual);
