<!--
SPDX-FileCopyrightText: 2026 CSC – IT Center for Science

SPDX-License-Identifier: Apache-2.0
-->

# AGENTS.md

## General

- Keep this file up to date; consult the user to resolve contradictions.
- You lack the project's full tacit context (requirements, style, goals): present options rather than decisions, and let the user make the final call.
- The list of tasks is in `TODO.md`. Keep this file updated.
- Never modify code in submodules; otherwise this would create licensing issues.
- Ensure the headers in files have correct licensing information and attribution.
- If in doubt or the situation is unclear, inform the user. Prefer Apache-2.0 but opt for GPLv3-or-later in unclear situations to be on the safe side.

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
- Use `bytes` for byte-length fields/parameters and `count` for element-length ones. Never the ambiguous `size`.
- Prefer raw `count`/pointer parameters over toolkit-internal struct types (e.g. `device_buffer_t`) in function signatures. This keeps coupling between modules low.

### Composability

- Build the library from small, composable units. Base-level wrappers do error checking and operate on fundamental data types; higher-level functionality is layered on top of these, not written from scratch alongside them.
- This lets the caller pick which level of abstraction to operate at. For example, `host_alloc_real_t` (allocates a raw `real_t*`) is the base wrapper. `host_buffer_create`, which builds the `host_buffer_t` struct (count + data), calls `host_alloc_real_t` internally rather than duplicating the allocation. A caller who only needs a raw array can call `host_alloc_real_t` directly; a caller who wants the struct's bookkeeping calls `host_buffer_create`.

## Performance

- Memory allocations must happen outside performance-critical functions.
    - Functions must not have side effects.
    - Memory allocation and deallocation are performed explicitly by the caller.
    - Functions must be concurrency-ready: data shared by concurrent functions must be immutable.
    - The location of memory allocation must be explicit (e.g., host, device).
- Deallocate in LIFO order (reverse of allocation), unless another constraint requires otherwise: allocate A, allocate B, then deallocate B, deallocate A.
- Functions return by value or `void` by default, and fail fast (e.g. via `ERRCHK`) on error. Return `hpc_status_t` only when the failure is meaningful and recoverable by the caller, and mark such functions `[[nodiscard]]` to enforce that the caller checks the result.

## Submodules

- `thirdparty/hpc-toolkit`: provides helper functions
