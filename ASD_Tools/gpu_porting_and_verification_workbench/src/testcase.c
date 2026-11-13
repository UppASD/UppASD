/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 *
 * SPDX-License-Identifier: GPL-3.0-or-later
 */

#include "testcase.h"

#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "common/errhandler.h"
#include "common/math-toolkit.h"
#include "host-toolkit/memory.h"

// Atom order of UppASD (multiscale.f90, createAtoms): real atoms, padding atoms, finite-difference
// nodes. A second padding block after the nodes is not in UppASD. It gives interpolated atoms
// with interpolated neighbours of lower index.
typedef struct {
    size_t real_count;
    size_t padding_count; // First padding block
    size_t node_count;
} layout_t;

static layout_t
layout_create(const size_t atom_count)
{
    return (layout_t){
        .real_count    = atom_count / 2,
        .padding_count = atom_count / 8,
        .node_count    = atom_count / 4,
    };
}

static size_t
node_begin(const layout_t layout)
{
    return layout.real_count + layout.padding_count;
}

static size_t
node_end(const layout_t layout)
{
    return node_begin(layout) + layout.node_count;
}

static bool
is_node(const layout_t layout, const size_t atom)
{
    return atom >= node_begin(layout) && atom < node_end(layout);
}

// Index of a padding atom across both padding blocks
static size_t
padding_index(const layout_t layout, const size_t atom)
{
    return atom < node_begin(layout) ? atom - layout.real_count
                                     : atom - node_end(layout) + layout.padding_count;
}

// Returns the number of neighbours of atom, 0 if it is not interpolated
static size_t
neighbour_count(const layout_t layout, const size_t atom)
{
    if (atom < layout.real_count)
        return 0;

    // Every second node is an interpolation node, with 1 to 8 real atoms as neighbours
    if (is_node(layout, atom)) {
        const size_t node = atom - node_begin(layout);
        return node % 2 == 0 ? 1 + (node / 2) % 8 : 0;
    }

    // Padding atoms have 2, 4 or 8 node corners, as in 1, 2 and 3 dimensions
    return (size_t)2 << (padding_index(layout, atom) % 3);
}

// Gap of 1 to 3 between neighbours i - 1 and i
static size_t
gap(const size_t key, const size_t i)
{
    return 1 + (key + i) % 3;
}

// Writes count sorted, distinct and non-consecutive 1-based neighbours and their weights. Varying
// gaps and weights from an irrational sequence avoid regular patterns. Padding neighbours are
// nodes, of which every second one is interpolated. weight is the index of the first weight.
static void
fill_neighbours(const layout_t layout, const size_t atom, const size_t count, const size_t weight,
                int32_t* neighbours, double* weights)
{
    const bool   padding    = !is_node(layout, atom);
    const size_t key        = padding ? padding_index(layout, atom) : (atom - node_begin(layout)) / 2;
    const size_t pool_begin = padding ? node_begin(layout) : 0;
    const size_t pool_count = padding ? layout.node_count : layout.real_count;

    size_t span = 0;
    for (size_t i = 1; i < count; ++i)
        span += gap(key, i);
    size_t neighbour = pool_begin + (5 * key) % (pool_count - span);

    double sum = 0;
    for (size_t i = 0; i < count; ++i) {
        if (i > 0)
            neighbour += gap(key, i);
        neighbours[i] = (int32_t)(neighbour + 1);
        weights[i]    = 0.5 + fmod(0.7548776662466927 * (double)(weight + i + 1), 1.0);
        sum += weights[i];
    }

    // Padding weights sum to 1, as in createPaddingInterpolationWeights
    if (padding)
        for (size_t i = 0; i < count; ++i)
            weights[i] /= sum;
}

[[noreturn]] static void
usage(const char* program)
{
    fprintf(stderr, "Usage: %s [--atoms N] [--ensembles N] [--no-verify]\n", program);
    exit(EXIT_FAILURE);
}

static size_t
parse_count(const char* value)
{
    const unsigned long long count = strtoull(value, nullptr, 10);
    ERRCHK(count <= SIZE_MAX);
    return (size_t)count;
}

testcase_options_t
testcase_parse_options(const int argc, char* const argv[])
{
    testcase_options_t options = {
        .atom_count     = 100'000,
        .ensemble_count = 3,
        .verify         = true,
    };
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "--atoms") == 0 && i + 1 < argc)
            options.atom_count = parse_count(argv[++i]);
        else if (strcmp(argv[i], "--ensembles") == 0 && i + 1 < argc)
            options.ensemble_count = parse_count(argv[++i]);
        else if (strcmp(argv[i], "--no-verify") == 0)
            options.verify = false;
        else
            usage(argv[0]);
    }
    return options;
}

testcase_t
testcase_host_allocate(const size_t atom_count, const size_t ensemble_count,
                       const size_t row_count, const size_t weight_count)
{
    const size_t mmom_count = atom_count * ensemble_count;

    testcase_t testcase = {
        .atom_count     = atom_count,
        .ensemble_count = ensemble_count,
        .row_count      = row_count,
        .weight_count   = weight_count,
    };
    testcase.indices         = alloc_host_int32_t(atom_count);
    testcase.first_neighbour = alloc_host_int32_t(row_count + 1);
    testcase.neighbours      = alloc_host_int32_t(weight_count);
    testcase.weights         = alloc_host_double(weight_count);
    testcase.mmom            = alloc_host_double(mmom_count);
    testcase.emom            = alloc_host_double(3 * mmom_count);
    testcase.emom2           = alloc_host_double(3 * mmom_count);
    testcase.emomM           = alloc_host_double(3 * mmom_count);
    return testcase;
}

testcase_t
testcase_host_create(const size_t atom_count, const size_t ensemble_count)
{
    // At least 22 real atoms and 22 nodes, the largest span of 8 neighbours with gaps up to 3
    ERRCHK(atom_count >= 88);
    ERRCHK(atom_count <= INT32_MAX / 8);
    ERRCHK(ensemble_count > 0);

    const layout_t layout = layout_create(atom_count);

    size_t row_count    = 0;
    size_t weight_count = 0;
    for (size_t atom = 0; atom < atom_count; ++atom) {
        const size_t count = neighbour_count(layout, atom);
        row_count += count != 0;
        weight_count += count;
    }

    testcase_t testcase = testcase_host_allocate(atom_count, ensemble_count, row_count,
                                                 weight_count);

    // Rows in ascending atom order, as in setupInterpolation (multiscalesetupsystem.f90)
    size_t row    = 0;
    size_t weight = 0;
    testcase.first_neighbour[0] = 1;
    for (size_t atom = 0; atom < atom_count; ++atom) {
        const size_t count = neighbour_count(layout, atom);
        if (count == 0) {
            testcase.indices[atom] = 0;
            continue;
        }

        testcase.indices[atom] = (int32_t)(row + 1);
        fill_neighbours(layout, atom, count, weight, &testcase.neighbours[weight],
                        &testcase.weights[weight]);
        weight += count;
        ++row;
        testcase.first_neighbour[row] = (int32_t)(weight + 1);
    }

    // Unit moments, as UppASD expects. Positive components keep the weighted sums away from zero
    // norm. Irrational steps make every moment distinct, so a wrong atom or ensemble offset fails.
    for (size_t m = 0; m < atom_count * ensemble_count; ++m) {
        testcase.mmom[m] = 1.0 + 0.1 * (double)m;

        double v[3];
        double v2[3];
        for (size_t k = 0; k < 3; ++k) {
            v[k]  = 1.0 + fmod(0.6180339887498949 * (double)(3 * m + k), 1.0);
            v2[k] = 1.0 + fmod(1.4142135623730951 * (double)(3 * m + k), 1.0);
        }
        normalize_real_t(3, v, &testcase.emom[3 * m]);
        normalize_real_t(3, v2, &testcase.emom2[3 * m]);

        for (size_t k = 0; k < 3; ++k)
            testcase.emomM[3 * m + k] = testcase.emom[3 * m + k] * testcase.mmom[m];
    }

    return testcase;
}

testcase_t
testcase_host_clone(const testcase_t src)
{
    const size_t mmom_count = src.atom_count * src.ensemble_count;

    testcase_t dst = testcase_host_allocate(src.atom_count, src.ensemble_count, src.row_count,
                                            src.weight_count);
    copy_h2h_int32_t(src.atom_count, src.indices, dst.indices);
    copy_h2h_int32_t(src.row_count + 1, src.first_neighbour, dst.first_neighbour);
    copy_h2h_int32_t(src.weight_count, src.neighbours, dst.neighbours);
    copy_h2h_double(src.weight_count, src.weights, dst.weights);
    copy_h2h_double(mmom_count, src.mmom, dst.mmom);
    copy_h2h_double(3 * mmom_count, src.emom, dst.emom);
    copy_h2h_double(3 * mmom_count, src.emom2, dst.emom2);
    copy_h2h_double(3 * mmom_count, src.emomM, dst.emomM);
    return dst;
}

void
testcase_host_destroy(testcase_t* testcase)
{
    ERRCHK(testcase != nullptr);

    free_host_double(&testcase->emomM);
    free_host_double(&testcase->emom2);
    free_host_double(&testcase->emom);
    free_host_double(&testcase->mmom);
    free_host_double(&testcase->weights);
    free_host_int32_t(&testcase->neighbours);
    free_host_int32_t(&testcase->first_neighbour);
    free_host_int32_t(&testcase->indices);
    *testcase = (testcase_t){};
}
