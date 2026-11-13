<!--
SPDX-FileCopyrightText: 2026 CSC – IT Center for Science

SPDX-License-Identifier: Apache-2.0
-->

# HPC Toolkit

A toolkit for common operations in high-performance computing applications.

# Building

```bash
mkdir build && cd build
cmake -DHPC_BUILD_TESTS=ON .. && cmake --build . --parallel
$SRUN ctest --output-on-failure
```

# Contributors

- Johannes Pekkilä (CSC – IT Center for Science)
