/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#include "host-toolkit/memory.h"

#include <stdint.h>
#include <stdio.h>
#include <string.h>

#include "common/errhandler.h"
#include "common/macros.h"

#define HPC_T float
#include "memory_def_template.h"
#undef HPC_T

#define HPC_T double
#include "memory_def_template.h"
#undef HPC_T

#define HPC_T int32_t
#include "memory_def_template.h"
#undef HPC_T

#define HPC_T uint64_t
#include "memory_def_template.h"
#undef HPC_T

#define HPC_T size_t
#include "memory_def_template.h"
#undef HPC_T
