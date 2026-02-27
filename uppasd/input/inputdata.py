"""
inputdata.py

Sparse, patch-based handling of UppASD input parameters.

Key design principles:
- Python writes ONLY explicitly specified keywords
- Fortran owns all defaults
- No schema mirroring
- Composable blocks, but no enforced hierarchy
"""

from pathlib import Path
from typing import Dict, Any, Optional
import yaml


class InputBlock:
    """
    A lightweight container for UppASD input keywords.

    Stores only explicitly set parameters.
    """

    def __init__(self, **kwargs):
        self._data: Dict[str, Any] = {}
        self.set(**kwargs)

    def __repr__(self):
        return f"InputBlock({self._data})"

    # ------------------------------------------------------------------
    # Core API
    # ------------------------------------------------------------------

    def set(self, **kwargs):
        """Set or update keyword values."""
        for key, value in kwargs.items():
            if value is None:
                continue
            self._data[key] = value

    def unset(self, *keys):
        """Remove keys if present."""
        for key in keys:
            self._data.pop(key, None)

    def clear(self):
        """Remove all keys."""
        self._data.clear()

    # ------------------------------------------------------------------
    # Introspection
    # ------------------------------------------------------------------

    def keys(self):
        return self._data.keys()

    def items(self):
        return self._data.items()

    def as_dict(self):
        return dict(self._data)

    # ------------------------------------------------------------------
    # YAML helpers
    # ------------------------------------------------------------------

    def apply_yaml(self, filename: str):
        """
        Apply a YAML patch.
        YAML keys are merged into this block.
        """
        with open(filename, "r") as f:
            data = yaml.safe_load(f)

        if data is None:
            return

        if not isinstance(data, dict):
            raise ValueError("YAML input must be a mapping")

        self.set(**data)

    def to_yaml(self, filename: str):
        """Write this block to YAML (for presets or inspection)."""
        with open(filename, "w") as f:
            yaml.safe_dump(self._data, f, sort_keys=True)

    # ------------------------------------------------------------------
    # Writing
    # ------------------------------------------------------------------

    def write(self, fh):
        """
        Write this block to an open file handle.
        Format: key value (space-separated, no equals sign)

        Multi-row values (2D arrays like cell) are written with:
        - First row: "key row1"
        - Subsequent rows: whitespace-padded to align
        """
        for key, value in self._data.items():
            formatted = self._format_value(value)
            if isinstance(formatted, list):
                # Multi-line values (e.g., cell matrix rows)
                # Write key only on first line, subsequent lines indented
                key_width = len(key)
                for i, line in enumerate(formatted):
                    if i == 0:
                        fh.write(f"{key} {line}\n")
                    else:
                        # Indent subsequent rows to align with first row data
                        fh.write(f"{' ' * (key_width + 1)}{line}\n")
            else:
                fh.write(f"{key} {formatted}\n")

    @staticmethod
    def _format_value(value):
        """
        Format values in a conservative UppASD-compatible way.

        Returns:
            - For 2D arrays (cell): list of space-separated row strings
            - For 1D arrays/lists: space-separated string
            - For scalars: the value as-is or formatted string
        """
        if isinstance(value, bool):
            return "T" if value else "F"
        if isinstance(value, (int, float)):
            return str(value)

        # Handle NumPy arrays
        if hasattr(value, "shape"):
            # 2D array (e.g., cell matrix)
            if len(value.shape) == 2:
                lines = []
                for row in value:
                    lines.append(" ".join(str(v) for v in row))
                return lines
            # 1D array
            elif len(value.shape) == 1:
                return " ".join(str(v) for v in value)

        # Handle Python lists/tuples
        if isinstance(value, (list, tuple)):
            # Check if it's a list of lists (2D)
            if value and isinstance(value[0], (list, tuple)):
                lines = []
                for row in value:
                    lines.append(" ".join(str(v) for v in row))
                return lines
            # 1D list
            return " ".join(str(v) for v in value)

        return str(value)


# ======================================================================
# Composite input container
# ======================================================================


class ASDInput:
    """
    A collection of InputBlocks written into a single inpsd.dat file.

    Blocks are written in insertion order.
    """

    def __init__(self):
        self.blocks: Dict[str, InputBlock] = {}

    def __repr__(self):
        return f"ASDInput(blocks={list(self.blocks.keys())})"

    def add_block(self, name: str, block: InputBlock):
        self.blocks[name] = block

    def block(self, name: str) -> InputBlock:
        """
        Get or create a block.
        """
        if name not in self.blocks:
            self.blocks[name] = InputBlock()
        return self.blocks[name]

    # ------------------------------------------------------------------
    # Writing
    # ------------------------------------------------------------------

    def write(self, filename: str):
        """
        Write all blocks to an UppASD input file.
        """
        path = Path(filename)
        with path.open("w") as f:
            for name, block in self.blocks.items():
                f.write(f"# --- {name} ---\n")
                block.write(f)
                f.write("\n")

    # ------------------------------------------------------------------
    # Convenience
    # ------------------------------------------------------------------

    def apply_yaml(self, filename: str, block: Optional[str] = None):
        """
        Apply a YAML patch to a specific block or all blocks.

        If block is None:
            YAML must contain a mapping of block_name -> dict
        """
        with open(filename, "r") as f:
            data = yaml.safe_load(f)

        if data is None:
            return

        if block is not None:
            self.block(block).set(**data)
            return

        if not isinstance(data, dict):
            raise ValueError("YAML must contain a mapping")

        for blk, values in data.items():
            if not isinstance(values, dict):
                raise ValueError(f"Block '{blk}' must map to a dict")
            self.block(blk).set(**values)


def _get_initial_block(inp: ASDInput) -> Optional[InputBlock]:
    """Return initial-phase block if present (supports both common names)."""
    if "initial" in inp.blocks:
        return inp.blocks["initial"]
    if "initialization" in inp.blocks:
        return inp.blocks["initialization"]
    return None


def _updated_sequence_temperature(seq, temperature: float):
    """Return seq with index 3 replaced by temperature, preserving container type."""
    if not isinstance(seq, (list, tuple)):
        raise TypeError("Initial-phase schedule must be list/tuple")
    if len(seq) < 4:
        raise ValueError("Initial-phase schedule must have at least 4 elements")

    items = list(seq)
    items[3] = float(temperature)
    return tuple(items) if isinstance(seq, tuple) else items


def set_temperature_token(inp: ASDInput, temperature: float):
    """
    Set temperature consistently for sweep workflows.

    Updates:
    - ``simulation.temperature``
    - ``initial.ip_temp`` when using scalar initial-phase syntax
    - ``initial.ip_nphase`` for schedule-based ``ip_mode='S'``
    - ``initial.ip_mcanneal`` for schedule-based ``ip_mode in {'M', 'H'}``

    Notes
    -----
    This helper does not invent missing initial-phase schedules.
    It updates whichever syntax is explicitly present in the input.
    """
    initial = _get_initial_block(inp)
    if initial is None:
        raise ValueError("Missing initial block ('initial' or 'initialization')")

    simulation = inp.block("simulation")
    simulation.set(temperature=float(temperature))

    data = initial.as_dict()
    ip_mode = str(data.get("ip_mode", "S")).strip().upper()

    if "ip_temp" in data:
        initial.set(ip_temp=float(temperature))
        return

    if ip_mode == "S":
        if "ip_nphase" not in data:
            raise ValueError("ip_mode='S' requires explicit ip_nphase in initial block")
        initial.set(ip_nphase=_updated_sequence_temperature(data["ip_nphase"], temperature))
        return

    if ip_mode in ("M", "H"):
        if "ip_mcanneal" not in data:
            raise ValueError("ip_mode='M'/'H' requires explicit ip_mcanneal in initial block")
        initial.set(ip_mcanneal=_updated_sequence_temperature(data["ip_mcanneal"], temperature))
        return

    raise ValueError(f"Unsupported ip_mode '{ip_mode}'. Expected 'S', 'M', or 'H'.")


def validate_temperature_token(inp: ASDInput) -> Dict[str, Any]:
    """
    Validate and report the active sweep temperature mapping.

    Returns
    -------
    dict
        ``{'ip_mode': ..., 'temperature_source': ..., 'initial_temperature': ..., 'simulation_temperature': ...}``
    """
    initial = _get_initial_block(inp)
    if initial is None:
        raise ValueError("Missing initial block ('initial' or 'initialization')")

    simulation = inp.block("simulation").as_dict()
    if "temperature" not in simulation:
        raise ValueError("simulation.temperature is not set")

    init_data = initial.as_dict()
    ip_mode = str(init_data.get("ip_mode", "S")).strip().upper()

    if "ip_temp" in init_data:
        initial_temperature = float(init_data["ip_temp"])
        temperature_source = "ip_temp"

        if ip_mode == "S" and "ip_nstep" not in init_data:
            raise ValueError("ip_mode='S' with ip_temp syntax requires explicit ip_nstep")
        if ip_mode in ("M", "H") and "ip_mcnstep" not in init_data:
            raise ValueError("ip_mode='M'/'H' with ip_temp syntax requires explicit ip_mcnstep")

    elif ip_mode == "S":
        if "ip_nphase" not in init_data:
            raise ValueError("ip_mode='S' requires explicit ip_nphase or ip_temp in initial block")
        initial_temperature = float(init_data["ip_nphase"][3])
        temperature_source = "ip_nphase"

    elif ip_mode in ("M", "H"):
        if "ip_mcanneal" not in init_data:
            raise ValueError("ip_mode='M'/'H' requires explicit ip_mcanneal or ip_temp in initial block")
        initial_temperature = float(init_data["ip_mcanneal"][3])
        temperature_source = "ip_mcanneal"
    else:
        raise ValueError(f"Unsupported ip_mode '{ip_mode}'. Expected 'S', 'M', or 'H'.")

    simulation_temperature = float(simulation["temperature"])

    if abs(initial_temperature - simulation_temperature) > 1e-12:
        raise ValueError(
            "Temperature token mismatch between initial-phase schedule and simulation.temperature"
        )

    return {
        "ip_mode": ip_mode,
        "temperature_source": temperature_source,
        "initial_temperature": initial_temperature,
        "simulation_temperature": simulation_temperature,
    }
