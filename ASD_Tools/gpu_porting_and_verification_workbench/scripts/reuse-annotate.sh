#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
#
# SPDX-License-Identifier: Apache-2.0

set -euo pipefail
cd "$(git rev-parse --show-toplevel)"

CSC="CSC – IT Center for Science"
UPPASD="UppASD contributors"

reuse download --all #Apache-2.0 GPL-3.0-or-later

# Derived from UppASD
reuse annotate --merge-copyrights --copyright "$UPPASD" --copyright "$CSC" --license GPL-3.0-or-later \
    tmp/gpuInterpolation.cpp \
    tmp/interpolation_kernels.cpp \
    src/device-implementation.cu \
    src/device-implementation.h \
    src/host-implementation.c \
    src/host-implementation.h \
    src/model.f90 \
    src/model.h

# Original work using GPLv3-licensed work
reuse annotate --merge-copyrights --copyright "$CSC" --license GPL-3.0-or-later \
    src/testcase.c \
    src/testcase.h \
    src/testcase-device.c \
    src/testcase-device.h \
    src/testcase-verify.c \
    src/testcase-verify.h \
    src/test-host.c \
    src/test-device.c

# Original work
reuse annotate --merge-copyrights --copyright "$CSC" --license Apache-2.0 \
    CMakeLists.txt \
    src/CMakeLists.txt \
    TODO.md \
    AGENTS.md \
    README.md \
    .gitlab-ci.yml \
    .gitmodules \
    .gitignore \
    scripts/reuse-annotate.sh

reuse lint
