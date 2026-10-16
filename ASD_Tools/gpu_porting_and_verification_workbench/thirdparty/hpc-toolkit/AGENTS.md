<!--
SPDX-FileCopyrightText: 2026 CSC – IT Center for Science

SPDX-License-Identifier: Apache-2.0
-->

# AGENTS.md

A multi-platform HPC toolbox of reusable components for scientific applications.

## General

- Keep this file up to date; consult the user to resolve contradictions.
- You lack the project's full tacit context (requirements, style, goals): present options rather than decisions, and let the user make the final call.
- The list of tasks is in `TODO.md`. Keep this file updated.

## Commenting, answering, and running

- Write like scientific prose: every word must carry information; cut what doesn't, e.g. write "It is significant" rather than "It is significant, not negligible".
- Avoid jargon with rare or field-dependent meanings unless it is unambiguous in context.
- Avoid em dashes, and avoid overusing semicolons or colons. Prefer short, standard sentences, e.g. "CUDA must be used throughout. AMD portability is achieved with HIP." rather than one sentence joining both.
- Avoid wedge sentences, where a clause is nested inside the main clause. Split them into separate sentences instead.
- Do not build or run the project. Diagnose issues from source and the problem description only.
- Do not fetch content onto this machine, e.g. via curl, git clone, or WebFetch. Plain web search is fine.

## Conventions

- Code must be self-documenting: avoid comments except for TODO items, disclosing, sharp edges or pitfalls, and similar.
- Interface functions must be commented and in a format processable with standard linting tools.
- Code should use modern C23 conventions and best practices, for example:
    - const correctness
    - generics within reason
    - nullptr
    - designated initializers
    - type-generic math via `_Generic` wrappers over `<math.h>` (e.g. `HPC_SQRT`), not `<tgmath.h>`, which fails with Clang-based compilers such as ROCm's
- The code is designed to be copied into, or interfaced with, other independent projects.
    - Typically C, C++, or Fortran, but Python interfacing is also possible.
- Use `bytes` for byte-length fields/parameters and `count` for element-length ones. Never the ambiguous `size`.
- Prefer raw `count`/pointer parameters over toolkit-internal struct types (e.g. `buffer_device_double`) in function signatures. This keeps coupling between modules low.

### Composability

- Build the library from small, composable units. Base-level wrappers do error checking and operate on fundamental data types; higher-level functionality is layered on top of these, not written from scratch alongside them.
- This lets the caller pick which level of abstraction to operate at. For example, `alloc_host_double` (allocates a raw `double*`) is the base wrapper. `buffer_create_host_double`, which builds the `buffer_host_double` struct (count + data), calls `alloc_host_double` internally rather than duplicating the allocation. A caller who only needs a raw array can call `alloc_host_double` directly; a caller who wants the struct's bookkeeping calls `buffer_create_host_double`.

### Type-generic templates

- Per-type code is written once in a template file without an include guard. The file is instantiated once per type with `#define HPC_T <type>`, `#include`, `#undef HPC_T`.
- Each toolkit has a declaration template in `include/<toolkit>/` (e.g. `memory_decl_template.h`) and a definition template in `src/` (e.g. `memory_def_template.h`). Each template starts, after the license header, with `// Template: no include guard, included once per type`.
- Names put the operation first and end with the actual type. Templates write them with `HPC_T_NAME(name)`, which appends `_<HPC_T>`: `HPC_T_NAME(alloc_host)` gives e.g. `alloc_host_double`. Other examples are `copy_h2d_float`, `buffer_create_device_int32_t` and the struct `buffer_host_size_t`. `HPC_T_NAME` is defined and `#undef`'d in each template, using `HPC_CONCAT3` from `common/macros.h`.
- The type set is `float`, `double`, `int32_t`, `uint64_t` and `size_t`. The header and the source must instantiate the same list.
- Memory has no `real_t` variants. `real_t` remains in the compute modules, so passing a buffer's `.data` to them compiles only when `HPC_REAL_TYPE` matches.
- Pitfalls:
    - Template macros carry the `HPC_` prefix. Unprefixed names such as `T` clash with user macros and with C++ template parameters in CUDA, HIP and Thrust headers.
    - Template files include nothing. All includes come before the first `#define HPC_T`.
    - Device templates compile as C++17 and must be valid C++.
- In `device-toolkit`, host-side buffers are pinned (`buffer_pinned_T`), and every buffer copy checks that counts match. `realloc` exists only in `host-toolkit`.

## Performance

- Memory allocations must happen outside performance-critical functions.
    - Functions must not have side effects.
    - Memory allocation and deallocation are performed explicitly by the caller.
    - Functions must be concurrency-ready: data shared by concurrent functions must be immutable.
    - The location of memory allocation must be explicit (e.g., host, device).
- Deallocate in LIFO order (reverse of allocation), unless another constraint requires otherwise: allocate A, allocate B, then deallocate B, deallocate A.
- Functions return by value or `void` by default, and fail fast (e.g. via `ERRCHK`) on error. Return `hpc_status_t` only when the failure is meaningful and recoverable by the caller, and mark such functions `[[nodiscard]]` to enforce that the caller checks the result.

## Toolkits

The currently available toolkits are:
- `common`: functionality used by virtually every project, e.g., timing, benchmarking, error checking, and debug printing.
- `mpi-toolkit`: MPI-specific toolkit.
- `device-toolkit`: GPU-specific toolkit, supporting both AMD and Nvidia GPUs. CUDA must be used throughout. AMD portability is achieved with a `hip.h` header and HIP.
- `host-toolkit`: host implementations of `device-toolkit` functionality.

## Structure

The project is divided into independent toolkits.
- `common` is the only toolkit every other toolkit may depend on.
- No other dependencies exist between toolkits.
