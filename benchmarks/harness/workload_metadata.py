"""Workload-metadata hooks (WP-02 section E).

Every case manifest names a `workload_metadata_method`: a parser registered
here that turns a generated run into the workload fields
`benchmark_record.v1.schema.json` requires -- `natom`,
`directed_interactions`, `mean_neighbors`, `max_neighbors` for neighbour
cases, or `natom`, `fft_grid`, `fft_grid_padded`, `fft_grid_points` for
dipole/FFT cases. Neither parser estimates a quantity from assumption where
real UppASD production input or output already supplies it:

* `neighbor_list_from_struct_output` reads UppASD's own
  ``struct.<simid>.out`` diagnostic dump, written by
  ``prn_exchange`` (source/Hamiltonian/printhamiltonian.f90) whenever a run
  enables ``do_prnstruct`` (1 or 4). The header carries the true `Natom` and
  maximum neighbour count directly; `directed_interactions` is the real
  count of (iatom, jatom) body lines, not an estimate.
* `fft_grid_from_replication` computes the dipole-FFT grid from the case's
  own supercell replication (``ncell`` -- N1/N2/N3 in
  source/Hamiltonian/dipolemanager.f90 are exactly the ncell dimensions,
  nothing FFT-specific) and the exact padding formula UppASD's own FFT
  backends use internally, `N_pad = 2*N - 1`
  (source/Hamiltonian/fftdipole_fftw.f90, fftdipole_mkl.f90). No UppASD
  build prints the padded grid anywhere, so this is the only way to obtain
  it without adding a new print statement to production code, which is out
  of scope here.
"""

from __future__ import annotations

import pathlib
import re

from harness import cases as cases_mod

_HEADER_NATOM_RE = re.compile(r"#\s*Number of atoms:\s*(\d+)")
_HEADER_MAXNEIGH_RE = re.compile(r"#\s*Maximum num of neighbours:\s*(\d+)")


class WorkloadMetadataError(ValueError):
    """A workload-metadata parser could not determine its fields."""


def count_basis_atoms(run_dir):
    posfile_name = cases_mod.read_keyword(run_dir, "posfile")
    if posfile_name is None:
        raise WorkloadMetadataError(f"{run_dir}: inpsd.dat has no posfile keyword")
    posfile_path = pathlib.Path(run_dir) / posfile_name
    if not posfile_path.is_file():
        raise WorkloadMetadataError(f"{run_dir}: posfile {posfile_name!r} does not exist")
    count = 0
    for line in posfile_path.read_text().splitlines():
        stripped = line.strip()
        if stripped and not stripped.startswith("#"):
            count += 1
    return count


def neighbor_list_from_struct_output(case, size, run_dir):
    """Parse ``struct.<simid>.out`` for real neighbour-workload metadata.

    Requires the run to have been generated with ``do_prnstruct`` set to 1
    or 4 and to have actually executed, so the struct file exists in
    ``run_dir``. Raises :class:`WorkloadMetadataError` if it is missing --
    it never falls back to guessing a neighbour count.
    """
    run_dir = pathlib.Path(run_dir)
    simid = cases_mod.read_simid(run_dir)
    struct_path = run_dir / f"struct.{simid}.out"
    if not struct_path.is_file():
        raise WorkloadMetadataError(
            f"{struct_path} does not exist -- rerun with do_prnstruct in "
            "{1, 4} so UppASD writes it"
        )

    natom = None
    max_neighbors = None
    directed_interactions = 0
    for line in struct_path.read_text().splitlines():
        if natom is None:
            match = _HEADER_NATOM_RE.search(line)
            if match:
                natom = int(match.group(1))
                continue
        if max_neighbors is None:
            match = _HEADER_MAXNEIGH_RE.search(line)
            if match:
                max_neighbors = int(match.group(1))
                continue
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        # One body line per directed (iatom, jatom) pair.
        directed_interactions += 1

    if natom is None or max_neighbors is None:
        raise WorkloadMetadataError(f"{struct_path}: missing header (Number of atoms / Maximum num of neighbours)")

    return {
        "natom": natom,
        "directed_interactions": directed_interactions,
        "mean_neighbors": directed_interactions / natom,
        "median_neighbors": None,
        "max_neighbors": max_neighbors,
    }


def fft_grid_from_replication(case, size, run_dir):
    """Derive the dipole-FFT grid from the case's own supercell replication.

    `fft_grid` is the physical [N1, N2, N3] the run was generated with;
    `fft_grid_padded` is the transform grid UppASD actually pads to
    internally. `natom` is the exact basis-atom count read from the run's
    own posfile, times the replication product -- not an assumption.
    """
    n1, n2, n3 = cases_mod.replication_to_ncell(case, size)
    fft_grid = [n1, n2, n3]
    fft_grid_padded = [2 * n1 - 1, 2 * n2 - 1, 2 * n3 - 1]
    fft_grid_points = fft_grid_padded[0] * fft_grid_padded[1] * fft_grid_padded[2]
    natom = count_basis_atoms(run_dir) * n1 * n2 * n3
    return {
        "natom": natom,
        "fft_grid": fft_grid,
        "fft_grid_padded": fft_grid_padded,
        "fft_grid_points": fft_grid_points,
    }


PARSERS = {
    "neighbor_list_from_struct_output": neighbor_list_from_struct_output,
    "fft_grid_from_replication": fft_grid_from_replication,
}


def expected_natom(case, size, run_dir):
    """Independently compute the atom count a run at ``(case, size)`` must report.

    Used by the WP-03 runner as a validity cross-check: it is derived from
    the case's own replication and the generated run's own posfile -- the
    same two real inputs every workload-metadata parser above already uses --
    never from whatever the executed run itself reports.
    """
    n1, n2, n3 = cases_mod.replication_to_ncell(case, size)
    return count_basis_atoms(run_dir) * n1 * n2 * n3


def compute_workload_metadata(case, size, run_dir):
    """Dispatch to the parser named by the case's `workload_metadata_method`."""
    method_name = case.manifest["workload_metadata_method"]
    try:
        parser = PARSERS[method_name]
    except KeyError:
        raise WorkloadMetadataError(
            f"case {case.id!r} names unknown workload_metadata_method {method_name!r}; "
            f"registered methods are {sorted(PARSERS)}"
        ) from None
    return parser(case, size, run_dir)
