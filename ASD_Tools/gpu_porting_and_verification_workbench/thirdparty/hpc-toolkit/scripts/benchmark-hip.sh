#!/usr/bin/env bash

# SPDX-FileCopyrightText: 2026 CSC – IT Center for Science
#
# SPDX-License-Identifier: Apache-2.0

module load PrgEnv-gnu # or PrgEnv-cray/PrgEnv-amd
module load craype-accel-amd-gfx90a
module load LUMI/25.03 partition/G rocm/6.3.4-extras

module load cray-python

python3 -m venv ~/rocprof-pyenv
source ~/rocprof-pyenv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install -r /pfs/lustrep3/appl/lumi/SW/LUMI-25.03/G/EB/rocm/6.3.4-extras/libexec/rocprofiler-compute/requirements.txt

# Generate kernel traces
rocprofv3 --stats --kernel-trace --truncate-kernels -- <APPLICATION>

# Perfetto HIP trace
rocprofv3 --hip-trace --output-format pftrace -- <APPLICATION>

# Print available metrics
rocprofv3 --list-metrics

# Collect metrics
rocprofv3 --input rocprof-counters.txt -- <APPLICATION>

# rocprof-compute
rocprof-compute profile --name profile -- <APPLICATION>

# rocprof-sys
rocprof-sys-run -- <APPLICATION>

# More information
# https://lumi-supercomputer.github.io/LUMI-EasyBuild-docs/r/rocm/#license-information
# https://462000265.lumidata.eu/paow-20260511/files/LUMI-paow-20260511-2_01_introduction-to-rocprof.pdf
