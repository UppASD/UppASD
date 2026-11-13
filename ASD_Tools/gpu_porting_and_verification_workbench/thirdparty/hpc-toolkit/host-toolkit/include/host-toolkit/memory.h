/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#pragma once

#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "common/datatypes.h"
#include "common/errhandler.h"
#include "common/macros.h"

#ifdef __cplusplus
extern "C" {
#endif

#define HPC_T float
#include "host-toolkit/memory_decl_template.h"
#undef HPC_T

#define HPC_T double
#include "host-toolkit/memory_decl_template.h"
#undef HPC_T

#define HPC_T int32_t
#include "host-toolkit/memory_decl_template.h"
#undef HPC_T

#define HPC_T uint64_t
#include "host-toolkit/memory_decl_template.h"
#undef HPC_T

#define HPC_T size_t
#include "host-toolkit/memory_decl_template.h"
#undef HPC_T

#ifdef __cplusplus
}
#endif
