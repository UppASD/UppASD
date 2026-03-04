"""
sweep.py

Temperature-sweep workflows for UppASD.

Responsibilities:
- Apply mode-aware temperature tokens to input
- Run one independent workspace per temperature
- Collect compact thermodynamic summaries
"""

import copy
from typing import Optional, Callable, Dict, Any

from uppasd.input.inputdata import ASDInput
from uppasd.input.inputdata import set_temperature_token, validate_temperature_token
from uppasd.core.exchange import ExchangeShellTable, DMIShellTable
from uppasd.core.anisotropy import AnisotropyTable
from uppasd.core.system import SpinSystem
from uppasd.core.results import ASDResults

from .simulator import ASDWorkspace, UppASDSimulator


def run_temperature_sweep(
    workdir_root: str,
    temperatures,
    system: SpinSystem,
    inp: ASDInput,
    exchange: Optional[ExchangeShellTable] = None,
    dmi: Optional[DMIShellTable] = None,
    ani: Optional[AnisotropyTable] = None,
    clean: bool = True,
    verbose: bool = True,
    progress_callback: Optional[Callable[[Dict[str, Any]], None]] = None,
):
    """
    Run independent measurement workflows across a temperature grid.

    For each temperature point this helper:
    1. Copies ``inp``
    2. Applies a consistent temperature token to initial/simulation blocks
    3. Prepares and runs a workspace
    4. Collects final thermodynamics via ``ASDResults.final_thermo()``

    Parameters
    ----------
    workdir_root : str
        Prefix for per-temperature workdirs.
    temperatures : iterable
        Sequence of temperature values.
    system : SpinSystem
        System definition.
    inp : ASDInput
        Base input template.
    exchange : ExchangeShellTable, optional
        Exchange table.
    dmi : DMIShellTable, optional
        DMI table.
    clean : bool, default=True
        Recreate per-temperature workspace if it exists.
    verbose : bool, default=True
        Print per-temperature progress messages.
    progress_callback : callable, optional
        Callback receiving progress payload dictionaries.

    Returns
    -------
    list[dict]
        One dictionary per temperature point.
    """
    sweep_results = []
    temp_list = [float(temp) for temp in temperatures]
    total = len(temp_list)

    for idx, temp_value in enumerate(temp_list, start=1):
        if verbose:
            print(f"[{idx}/{total}] Running temperature sweep point: T={temp_value:.6g} K")

        if progress_callback is not None:
            progress_callback({
                "stage": "start",
                "index": idx,
                "total": total,
                "temperature": temp_value,
            })

        inp_point = copy.deepcopy(inp)

        set_temperature_token(inp_point, temp_value)
        validate_temperature_token(inp_point)

        temp_label = str(temp_value).replace(".", "p")
        workdir = f"{workdir_root}_T{temp_label}"

        ws = ASDWorkspace(workdir, clean=clean)
        ws.prepare(system=system, inp=inp_point, exchange=exchange, dmi=dmi, ani=ani)

        sim = UppASDSimulator(ws)
        sim.initialize()
        sim.run_init_phase()
        sim.measure()
        sim.finalize()

        results = ASDResults(ws.path)
        row = {
            "temperature": temp_value,
            "workdir": str(ws.path),
        }

        thermo = results.final_thermo()
        if thermo is not None:
            row.update(thermo)

        if progress_callback is not None:
            progress_callback({
                "stage": "done",
                "index": idx,
                "total": total,
                "temperature": temp_value,
                "result": row,
            })

        # if verbose:
        #     print(f"[{idx}/{total}] Completed T={temp_value:.6g} K")

        sweep_results.append(row)

    return sweep_results
