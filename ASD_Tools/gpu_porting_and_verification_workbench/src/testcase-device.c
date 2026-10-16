/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "testcase-device.h"

#include "common/errhandler.h"
#include "device-toolkit/memory.h"

testcase_t
testcase_device_allocate(const size_t atom_count, const size_t ensemble_count,
                         const size_t row_count, const size_t weight_count)
{
    const size_t mmom_count = atom_count * ensemble_count;

    testcase_t testcase = {
        .atom_count     = atom_count,
        .ensemble_count = ensemble_count,
        .row_count      = row_count,
        .weight_count   = weight_count,
    };
    testcase.indices         = alloc_device_int32_t(atom_count);
    testcase.first_neighbour = alloc_device_int32_t(row_count + 1);
    testcase.neighbours      = alloc_device_int32_t(weight_count);
    testcase.weights         = alloc_device_double(weight_count);
    testcase.mmom            = alloc_device_double(mmom_count);
    testcase.emom            = alloc_device_double(3 * mmom_count);
    testcase.emom2           = alloc_device_double(3 * mmom_count);
    testcase.emomM           = alloc_device_double(3 * mmom_count);
    return testcase;
}

testcase_t
testcase_device_clone_h2d(const testcase_t host)
{
    const size_t mmom_count = host.atom_count * host.ensemble_count;

    testcase_t device = testcase_device_allocate(host.atom_count, host.ensemble_count,
                                                 host.row_count, host.weight_count);
    copy_h2d_int32_t(host.atom_count, host.indices, device.indices);
    copy_h2d_int32_t(host.row_count + 1, host.first_neighbour, device.first_neighbour);
    copy_h2d_int32_t(host.weight_count, host.neighbours, device.neighbours);
    copy_h2d_double(host.weight_count, host.weights, device.weights);
    copy_h2d_double(mmom_count, host.mmom, device.mmom);
    copy_h2d_double(3 * mmom_count, host.emom, device.emom);
    copy_h2d_double(3 * mmom_count, host.emom2, device.emom2);
    copy_h2d_double(3 * mmom_count, host.emomM, device.emomM);
    return device;
}

testcase_t
testcase_host_clone_d2h(const testcase_t device)
{
    const size_t mmom_count = device.atom_count * device.ensemble_count;

    testcase_t host = testcase_host_allocate(device.atom_count, device.ensemble_count,
                                             device.row_count, device.weight_count);
    copy_d2h_int32_t(device.atom_count, device.indices, host.indices);
    copy_d2h_int32_t(device.row_count + 1, device.first_neighbour, host.first_neighbour);
    copy_d2h_int32_t(device.weight_count, device.neighbours, host.neighbours);
    copy_d2h_double(device.weight_count, device.weights, host.weights);
    copy_d2h_double(mmom_count, device.mmom, host.mmom);
    copy_d2h_double(3 * mmom_count, device.emom, host.emom);
    copy_d2h_double(3 * mmom_count, device.emom2, host.emom2);
    copy_d2h_double(3 * mmom_count, device.emomM, host.emomM);
    return host;
}

void
testcase_device_destroy(testcase_t* testcase)
{
    ERRCHK(testcase != nullptr);

    free_device_double(&testcase->emomM);
    free_device_double(&testcase->emom2);
    free_device_double(&testcase->emom);
    free_device_double(&testcase->mmom);
    free_device_double(&testcase->weights);
    free_device_int32_t(&testcase->neighbours);
    free_device_int32_t(&testcase->first_neighbour);
    free_device_int32_t(&testcase->indices);
    *testcase = (testcase_t){};
}
