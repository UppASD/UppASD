/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#include "host-toolkit/reduce.h"

#include "common/errhandler.h"

void
reduce_sum(const size_t count, const real_t* data, real_t* result)
{
    ERRCHK(count > 0);
    ERRCHK(data != nullptr);
    ERRCHK(result != nullptr);

    real_t sum = 0;
    for (size_t i = 0; i < count; ++i)
        sum += data[i];

    *result = sum;
}
