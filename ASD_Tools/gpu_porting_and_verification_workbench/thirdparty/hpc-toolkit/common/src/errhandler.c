/*
 * SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
 *
 * SPDX-License-Identifier: Apache-2.0
 */

#include <stdio.h>
#include <stdlib.h>

#include "common/errhandler.h"

static void
errhandler_default_callback(void)
{
    exit(EXIT_FAILURE);
}

static void (*errhandler_terminate)(void) = errhandler_default_callback;

void
errhandler_atexit(void (*callback)(void))
{
    errhandler_terminate = callback;
}

static void
print(const char* type, const char* file, const int line, const char* expr, const int result,
      const char* name, const char* description)
{
    fprintf(stderr, "%s: %s:%d; %s:%d; %s:%s\n", type, file, line, expr, result, name, description);
}

void
errhandler_log(const char* file, const int line, const char* expr, const int result,
               const char* name, const char* description)
{
    print("Log", file, line, expr, result, name, description);
}

void
errhandler_warn(const char* file, const int line, const char* expr, const int result,
                const char* name, const char* description)
{
    print("WARNING", file, line, expr, result, name, description);
}

void
errhandler_err(const char* file, const int line, const char* expr, const int result,
               const char* name, const char* description)
{
    print("ERROR", file, line, expr, result, name, description);
    errhandler_terminate();
}
