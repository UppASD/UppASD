/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "testcase-verify.h"

#include <float.h>
#include <limits.h>
#include <stdio.h>

#include "common/errhandler.h"
#include "common/verify.h"

#include "model.h"

static void
print_moment(const testcase_t input, const size_t moment, const double* in, const double* expected,
             const double* actual)
{
    const size_t  base     = 3 * moment;
    const char*   labels[] = {"input", "reference", "candidate"};
    const double* arrays[] = {in, expected, actual};

    printf("  atom %zu, ensemble %zu\n", moment % input.atom_count, moment / input.atom_count);
    for (size_t i = 0; i < 3; ++i)
        printf("    %-9s (% .16e, % .16e, % .16e)\n", labels[i], arrays[i][base],
               arrays[i][base + 1], arrays[i][base + 2]);
}

// Prints the first and last few interpolated moments. The leading atoms are real atoms, which are
// not interpolated.
static void
print_samples(const char* label, const testcase_t input, const double* in, const double* expected,
              const double* actual)
{
    constexpr size_t sample_count = 3;

    const size_t moment_count = input.atom_count * input.ensemble_count;
    printf("%s samples (interpolated moments):\n", label);

    size_t printed = 0;
    for (size_t m = 0; m < moment_count && printed < sample_count; ++m) {
        if (input.indices[m % input.atom_count] != 0) {
            print_moment(input, m, in, expected, actual);
            ++printed;
        }
    }

    size_t last = moment_count;
    for (size_t found = 0; last > 0 && found < sample_count;)
        found += input.indices[--last % input.atom_count] != 0;

    printf("  ...\n");
    for (size_t m = last; m < moment_count; ++m)
        if (input.indices[m % input.atom_count] != 0)
            print_moment(input, m, in, expected, actual);
}

// Arrays have the Fortran shape (3, atom_count, ensemble_count)
static void
print_position(const char* label, const size_t index, const size_t atom_count)
{
    printf("%s index %zu: component %zu, atom %zu, ensemble %zu\n", label, index,
           index % 3, (index / 3) % atom_count, index / (3 * atom_count));
}

static bool
compare(const char* label, const size_t count, const size_t atom_count, const double* expected,
        const double* actual)
{
    // A few ulp of slack, e.g. for FMA contraction on the device
    constexpr long double tolerance = 4 * (long double)DBL_EPSILON;

    const verify_result_t result = verify_real_t(count, expected, actual);
    // NaN must fail. Keep the comparison as <=, because !(x > tolerance) passes NaN.
    const bool pass = result.max_rel_error <= tolerance;
    printf("%s: %s\n", label, pass ? "OK" : "FAIL");
    verify_print(result);
    print_position("Max absolute error", result.max_abs_error_index, atom_count);
    print_position("Max relative error", result.max_rel_error_index, atom_count);
    return pass;
}

void
testcase_run_model(testcase_t* testcase)
{
    ERRCHK(testcase != nullptr);
    ERRCHK(testcase->atom_count <= INT_MAX);
    ERRCHK(testcase->ensemble_count <= INT_MAX);
    ERRCHK(testcase->row_count <= INT_MAX);
    ERRCHK(testcase->weight_count <= INT_MAX);

    model_multiscale_interpolate_interfaces((int)testcase->atom_count,
                                            (int)testcase->ensemble_count,
                                            (int)testcase->row_count,
                                            (int)testcase->weight_count,
                                            testcase->indices,
                                            testcase->first_neighbour,
                                            testcase->neighbours,
                                            testcase->weights,
                                            testcase->mmom,
                                            testcase->emom,
                                            testcase->emom2,
                                            testcase->emomM);
}

bool
testcase_compare(const testcase_t input, const testcase_t expected, const testcase_t actual)
{
    ERRCHK(expected.atom_count == input.atom_count);
    ERRCHK(expected.ensemble_count == input.ensemble_count);
    ERRCHK(actual.atom_count == input.atom_count);
    ERRCHK(actual.ensemble_count == input.ensemble_count);

    const size_t count = 3 * input.atom_count * input.ensemble_count;

    print_samples("emom", input, input.emom, expected.emom, actual.emom);
    print_samples("emom2", input, input.emom2, expected.emom2, actual.emom2);
    print_samples("emomM", input, input.emomM, expected.emomM, actual.emomM);

    bool pass = true;
    pass &= compare("emom", count, input.atom_count, expected.emom, actual.emom);
    pass &= compare("emom2", count, input.atom_count, expected.emom2, actual.emom2);
    pass &= compare("emomM", count, input.atom_count, expected.emomM, actual.emomM);
    return pass;
}
