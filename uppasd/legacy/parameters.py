from dataclasses import dataclass, field
from typing import Optional

@dataclass
class UppASDParameters:
    # --- Simulation control ---
    simid: str = "test"
    nstep: int = 10000
    timestep: float = 1.0e-15
    temperature: float = 300.0

    # --- Dynamics ---
    damping: float = 0.01
    gyromagnetic_ratio: float = 1.76e11

    # --- Flags ---
    do_relaxation: bool = False
    do_measurements: bool = True

    # --- Optional blocks ---
    mc: Optional[dict] = None
    llg: Optional[dict] = None
