/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include <stdlib.h>

#include "device-implementation.h"
#include "testcase-device.h"
#include "testcase-verify.h"
#include "testcase.h"

int
main(int argc, char* argv[])
{
    const testcase_options_t options = testcase_parse_options(argc, argv);

    testcase_t input = testcase_host_create(options.atom_count, options.ensemble_count);

    testcase_t d_input = testcase_device_clone_h2d(input);
    // The implementation writes emomM only at interpolated atoms
    testcase_t d_output = testcase_device_clone_h2d(input);
    device_multiscale_interpolate_interfaces(d_input.atom_count, d_input.ensemble_count,
                                             d_input.indices, d_input.first_neighbour,
                                             d_input.neighbours, d_input.weights, d_input.mmom,
                                             d_input.emom, d_input.emom2, d_output.emom,
                                             d_output.emom2, d_output.emomM);

    bool pass = true;
    if (options.verify) {
        testcase_t actual   = testcase_host_clone_d2h(d_output);
        testcase_t expected = testcase_host_clone(input);
        testcase_run_model(&expected);
        pass = testcase_compare(input, expected, actual);
        testcase_host_destroy(&expected);
        testcase_host_destroy(&actual);
    }

    testcase_device_destroy(&d_output);
    testcase_device_destroy(&d_input);
    testcase_host_destroy(&input);
    return pass ? EXIT_SUCCESS : EXIT_FAILURE;
}
