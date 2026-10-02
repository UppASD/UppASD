/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <stdlib.h>

#include "host-implementation.h"
#include "testcase-verify.h"
#include "testcase.h"

int
main(int argc, char* argv[])
{
    const testcase_options_t options = testcase_parse_options(argc, argv);

    testcase_t input = testcase_host_create(options.atom_count, options.ensemble_count);

    // The implementation writes emomM only at interpolated atoms
    testcase_t actual = testcase_host_clone(input);
    host_multiscale_interpolate_interfaces(input.atom_count, input.ensemble_count, input.indices,
                                           input.first_neighbour, input.neighbours, input.weights,
                                           input.mmom, input.emom, input.emom2, actual.emom,
                                           actual.emom2, actual.emomM);

    bool pass = true;
    if (options.verify) {
        testcase_t expected = testcase_host_clone(input);
        testcase_run_model(&expected);
        pass = testcase_compare(input, expected, actual);
        testcase_host_destroy(&expected);
    }

    testcase_host_destroy(&actual);
    testcase_host_destroy(&input);
    return pass ? EXIT_SUCCESS : EXIT_FAILURE;
}
