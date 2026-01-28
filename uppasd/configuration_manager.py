"""
Configuration management for UppASD measurement and simulation parameters.

This module provides a ConfigurationManager class that handles:
- Configuration validation against JSON schema
- Dependency checking between parameters
- Preset loading and composition
- Export to Fortran-compatible format

Key principle: Fortran provides sensible defaults for all parameters.
Python only handles user-provided overrides (deltas from defaults).
"""
from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Dict, List, Optional, Any
import warnings

try:
    import jsonschema
    HAS_JSONSCHEMA = True
except ImportError:
    HAS_JSONSCHEMA = False

try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False


logger = logging.getLogger(__name__)


class ConfigurationManager:
    """
    Manages UppASD measurement and simulation configuration.
    
    Handles configuration validation, preset composition, and export.
    Works with Fortran's default parameter values - only manages user overrides.
    
    Examples:
        Basic usage with presets:
        >>> cfg = ConfigurationManager()
        >>> cfg.add_preset('measurement_basic')
        >>> cfg.validate()
        >>> kwargs = cfg.to_fortran_kwargs()
        
        With custom parameters:
        >>> cfg = ConfigurationManager()
        >>> cfg.add_measurement(do_avrg='Y', avrg_step=50)
        >>> cfg.validate()
        
        Composition:
        >>> cfg = ConfigurationManager()
        >>> cfg.add_preset('measurement_spectroscopy')
        >>> cfg.add_field('microwave', mwf_freq=15.0)
        >>> cfg.validate()
    """
    
    # Load schema and presets at class initialization
    _SCHEMA_PATH = Path(__file__).parent / 'config_schema.json'
    _PRESETS_PATH = Path(__file__).parent / 'presets.json'
    
    _SCHEMA = None
    _PRESETS = None
    
    @classmethod
    def _load_schema(cls):
        """Load JSON schema from file."""
        if cls._SCHEMA is None:
            try:
                with open(cls._SCHEMA_PATH, 'r') as f:
                    cls._SCHEMA = json.load(f)
                logger.info(f"Loaded schema from {cls._SCHEMA_PATH}")
            except FileNotFoundError:
                logger.warning(f"Schema file not found at {cls._SCHEMA_PATH}")
                cls._SCHEMA = {}
        return cls._SCHEMA
    
    @classmethod
    def _load_presets(cls):
        """Load presets from file."""
        if cls._PRESETS is None:
            try:
                with open(cls._PRESETS_PATH, 'r') as f:
                    cls._PRESETS = json.load(f)
                logger.info(f"Loaded presets from {cls._PRESETS_PATH}")
            except FileNotFoundError:
                logger.warning(f"Presets file not found at {cls._PRESETS_PATH}")
                cls._PRESETS = {}
        return cls._PRESETS
    
    def __init__(self, base_config: Optional[Dict] = None):
        """
        Initialize ConfigurationManager.
        
        Args:
            base_config: Optional base configuration dict to start with
        """
        self.config = {}
        if base_config:
            self.merge(base_config)
        self.schema = self._load_schema()
        self.presets = self._load_presets()
    
    def add_preset(self, preset_name: str) -> 'ConfigurationManager':
        """
        Add a preset configuration.
        
        Args:
            preset_name: Name of preset (e.g., 'measurement_basic')
            
        Returns:
            Self for method chaining
            
        Raises:
            ValueError: If preset not found
        """
        if preset_name not in self.presets:
            raise ValueError(
                f"Preset '{preset_name}' not found. "
                f"Available: {list(self.presets.keys())}"
            )
        
        preset = self.presets[preset_name]
        # Extract config (everything except 'description')
        preset_config = {k: v for k, v in preset.items() if k != 'description'}
        self.merge(preset_config)
        logger.info(f"Added preset: {preset_name}")
        return self
    
    def add_measurement(self, **kwargs) -> 'ConfigurationManager':
        """
        Add measurement output configuration.
        
        Args:
            **kwargs: Measurement parameters (e.g., do_avrg='Y', avrg_step=100)
            
        Returns:
            Self for method chaining
        """
        if 'measurement' not in self.config:
            self.config['measurement'] = {}
        self.config['measurement'].update(kwargs)
        logger.debug(f"Added measurement config: {kwargs}")
        return self
    
    def add_spectroscopy(self, **kwargs) -> 'ConfigurationManager':
        """
        Add spectroscopy configuration.
        
        Args:
            **kwargs: Spectroscopy parameters (e.g., do_ams='Y', ams_nstep=500)
            
        Returns:
            Self for method chaining
        """
        if 'spectroscopy' not in self.config:
            self.config['spectroscopy'] = {}
        self.config['spectroscopy'].update(kwargs)
        logger.debug(f"Added spectroscopy config: {kwargs}")
        return self
    
    def add_protocol(self, **kwargs) -> 'ConfigurationManager':
        """
        Add protocol modification configuration.
        
        Args:
            **kwargs: Protocol parameters (e.g., do_tempexp='Y', tempexp_start=100)
            
        Returns:
            Self for method chaining
        """
        if 'protocol' not in self.config:
            self.config['protocol'] = {}
        self.config['protocol'].update(kwargs)
        logger.debug(f"Added protocol config: {kwargs}")
        return self
    
    def add_field(self, **kwargs) -> 'ConfigurationManager':
        """
        Add field/interaction configuration.
        
        Args:
            **kwargs: Field parameters (e.g., do_mwf='Y', mwf_freq=10.0)
            
        Returns:
            Self for method chaining
        """
        if 'field' not in self.config:
            self.config['field'] = {}
        self.config['field'].update(kwargs)
        logger.debug(f"Added field config: {kwargs}")
        return self
    
    def add_analysis(self, **kwargs) -> 'ConfigurationManager':
        """
        Add analysis configuration.
        
        Args:
            **kwargs: Analysis parameters (e.g., do_wl='Y')
            
        Returns:
            Self for method chaining
        """
        if 'analysis' not in self.config:
            self.config['analysis'] = {}
        self.config['analysis'].update(kwargs)
        logger.debug(f"Added analysis config: {kwargs}")
        return self
    
    def add_coupling(self, **kwargs) -> 'ConfigurationManager':
        """
        Add coupling configuration.
        
        Args:
            **kwargs: Coupling parameters (e.g., do_sld='Y')
            
        Returns:
            Self for method chaining
        """
        if 'coupling' not in self.config:
            self.config['coupling'] = {}
        self.config['coupling'].update(kwargs)
        logger.debug(f"Added coupling config: {kwargs}")
        return self
    
    def merge(self, other: Dict) -> 'ConfigurationManager':
        """
        Merge another configuration dict into this one.
        
        Args:
            other: Configuration dict to merge
            
        Returns:
            Self for method chaining
        """
        for key, value in other.items():
            if key in self.config and isinstance(self.config[key], dict):
                self.config[key].update(value)
            else:
                self.config[key] = value
        logger.debug(f"Merged config: {other}")
        return self
    
    def validate(self) -> List[str]:
        """
        Validate configuration against schema and check dependencies.
        
        Returns:
            List of validation issues (empty if valid)
        """
        issues = []
        
        # JSON schema validation (if available)
        if HAS_JSONSCHEMA and self.schema:
            try:
                jsonschema.validate(self.config, self.schema)
                logger.debug("Schema validation passed")
            except jsonschema.ValidationError as e:
                issues.append(f"Schema validation error: {e.message}")
                logger.error(f"Schema validation failed: {e}")
        
        # Dependency checks
        issues.extend(self._check_dependencies())
        
        if issues:
            logger.warning(f"Configuration validation issues: {issues}")
        else:
            logger.info("Configuration validation successful")
        
        return issues
    
    def _check_dependencies(self) -> List[str]:
        """Check for parameter dependencies and conflicts."""
        issues = []
        
        spec = self.config.get('spectroscopy', {})
        proto = self.config.get('protocol', {})
        
        # Check 1: do_sc requires qvectors
        if spec.get('do_sc') == 'Q':
            if 'qpoints' not in spec and 'qfile' not in spec:
                issues.append(
                    "WARNING: do_sc='Q' requires qvectors. "
                    "Provide qpoints or qfile in spectroscopy config."
                )
                logger.warning("Structure factor requires Q-vectors")
        
        # Check 2: do_tempexp and do_3tm are mutually exclusive
        if proto.get('do_tempexp') == 'Y' and proto.get('do_3tm') == 'Y':
            issues.append(
                "ERROR: Cannot enable both do_tempexp and do_3tm. "
                "These are mutually exclusive temperature protocols."
            )
            logger.error("Conflicting protocol settings")
        
        return issues
    
    def to_fortran_kwargs(self) -> Dict[str, Any]:
        """
        Convert to Fortran-compatible kwargs dictionary.
        
        This flattens the hierarchical structure into a single dict
        that can be passed to sim.measure() or sim.relax().
        
        Returns:
            Dictionary of all configuration parameters
        """
        flat = {}
        for category, params in self.config.items():
            if isinstance(params, dict):
                flat.update(params)
        logger.debug(f"Converted to Fortran kwargs: {flat}")
        return flat
    
    def to_json(self, indent: int = 2) -> str:
        """
        Export configuration as JSON string.
        
        Args:
            indent: JSON indentation level
            
        Returns:
            JSON string representation
        """
        return json.dumps(self.config, indent=indent)
    
    def to_yaml(self) -> str:
        """
        Export configuration as YAML string (if PyYAML available).
        
        Returns:
            YAML string representation
            
        Raises:
            ImportError: If PyYAML not installed
        """
        if not HAS_YAML:
            raise ImportError(
                "PyYAML not available. Install with: pip install pyyaml"
            )
        return yaml.dump(self.config, default_flow_style=False)
    
    @staticmethod
    def from_json(json_str: str) -> 'ConfigurationManager':
        """
        Create ConfigurationManager from JSON string.
        
        Args:
            json_str: JSON configuration string
            
        Returns:
            New ConfigurationManager instance
        """
        config = json.loads(json_str)
        return ConfigurationManager(config)
    
    @staticmethod
    def from_yaml(yaml_str: str) -> 'ConfigurationManager':
        """
        Create ConfigurationManager from YAML string (if PyYAML available).
        
        Args:
            yaml_str: YAML configuration string
            
        Returns:
            New ConfigurationManager instance
            
        Raises:
            ImportError: If PyYAML not installed
        """
        if not HAS_YAML:
            raise ImportError(
                "PyYAML not available. Install with: pip install pyyaml"
            )
        config = yaml.safe_load(yaml_str)
        return ConfigurationManager(config)
    
    def __repr__(self) -> str:
        """String representation."""
        categories = list(self.config.keys())
        total_params = sum(len(v) for v in self.config.values() if isinstance(v, dict))
        return f"ConfigurationManager({categories}, {total_params} params)"
    
    def __str__(self) -> str:
        """Human-readable string."""
        return self.to_json()

# ============================================================================
# ConfigurationBuilder: Fluent API for Configuration Construction
# ============================================================================


class ConfigurationBuilder:
    """Comprehensive fluent API for building UppASD simulation configurations.
    
    The single entry point for all file-based simulation setup. Handles:
    - Structural parameters (simid, ncell, BC, cell, posfile, etc.)
    - Simulation parameters (mode, temperature, steps, damping, etc.)
    - Measurement flags (plotenergy, do_avrg, do_sc, etc.)
    - Advanced features (spectroscopy, protocols, fields, analysis, coupling)
    
    Principle: Separate concerns clearly, compose via method chaining.
    
    Example (file-based workflow):
        >>> cfg = (ConfigurationBuilder()
        ...     .structure(simid='HeisWire', ncell=[1,1,100], BC=['0','0','P'])
        ...     .files(posfile='./posfile', exchange='./jfile', momfile='./momfile')
        ...     .simulate(mode='S', temperature=50, steps=2000, damping=0.1)
        ...     .measure(plotenergy=1, do_avrg='Y', avrg_step=100)
        ...     .build())
        >>> write_inpsd_file('inpsd.dat', cfg.to_fortran_kwargs())
        >>> with Simulator() as sim:
        ...     sim.relax(**cfg.relax_kwargs())
        ...     sim.measure(**cfg.measure_kwargs())
    """
    
    # Known measurement flags (validated set for sim.measure)
    MEASUREMENT_FLAGS = {
        'plotenergy', 'do_avrg', 'avrg_step',
        'do_cumu', 'cumu_step', 'cumu_buff',
        'do_tottraj', 'tottraj_step', 'do_proj_avrg',
        'do_sc', 'sc_step', 'sc_nstep', 'sc_emax', 'sc_eres',
        'sc_window_fun', 'sc_average', 'do_sc_local_axis',
        'do_sc_tens', 'qpoints', 'qfile', 'nc_qvect',
        'do_ams', 'do_diamag', 'do_magdos', 'magdos_freq', 'magdos_sigma',
        'do_stiffness', 'eta_max', 'eta_min', 'alat',
        'do_tempexp', 'tempexp_start', 'tempexp_end', 'tempexp_tau', 'tempexp_step',
        'gradient', 'temperature', 'skyno', 'stt', 'jvec', 'adibeta',
        'do_spintemp', 'use_vsl', 'do_jtensor', 'calc_jtensor', 'do_reduced',
    }
    
    # Known simulation parameters for sim.relax()
    SIMULATION_PARAMS = {
        'mode', 'ip_mode', 'temperature', 'temp', 'steps', 'nstep',
        'timestep', 'delta_t', 'damping', 'mplambda1'
    }
    
    # Structural parameters (written to inpsd, never passed to Simulator)
    STRUCTURAL_PARAMS = {
        'simid', 'ncell', 'BC', 'cell', 'Sym', 'posfile', 'exchange',
        'momfile', 'do_prnstruct', 'Mensemble', 'Initmag', 'mseed'
    }
    
    def __init__(self) -> None:
        """Initialize builder with empty configuration."""
        self._config = ConfigurationManager()
        self._structural = {}  # Structural parameters (separate track)
        self._simulation = {}  # Simulation parameters (separate track)
        logger.debug("ConfigurationBuilder initialized")
    
    # =========================================================================
    # STRUCTURAL SETUP (System definition, written to inpsd.dat)
    # =========================================================================
    
    def structure(
        self,
        simid: str,
        ncell: List[int],
        BC: Optional[List[str] | str] = None,
        cell: Optional[Any] = None,
        Sym: int = 1,
        **kwargs
    ) -> 'ConfigurationBuilder':
        """
        Set structural parameters (system definition).
        
        Parameters
        ----------
        simid : str
            Simulation ID
        ncell : list of int
            Supercell size [n1, n2, n3]
        BC : list of str or str, optional
            Boundary conditions (e.g., ['0', '0', 'P'])
        cell : ndarray, optional
            3x3 lattice vectors
        Sym : int, optional
            Symmetry flag. Default: 1
        **kwargs
            Additional structural parameters (Mensemble, Initmag, mseed, etc.)
        
        Returns
        -------
        self
            For method chaining
        """
        self._structural.update({
            'simid': simid,
            'ncell': ncell,
            'Sym': Sym,
        })
        if BC is not None:
            self._structural['BC'] = BC
        if cell is not None:
            self._structural['cell'] = cell
        self._structural.update(kwargs)
        logger.debug(f"Set structure: simid={simid}, ncell={ncell}")
        return self
    
    def files(
        self,
        posfile: str,
        exchange: str,
        momfile: Optional[str] = None,
        **kwargs
    ) -> 'ConfigurationBuilder':
        """
        Set input file paths.
        
        Parameters
        ----------
        posfile : str
            Path to posfile.dat
        exchange : str
            Path to jfile.dat (exchange interactions)
        momfile : str, optional
            Path to momfile.dat
        **kwargs
            Other file paths (anisfile, etc.)
        
        Returns
        -------
        self
            For method chaining
        """
        self._structural.update({
            'posfile': posfile,
            'exchange': exchange,
        })
        if momfile is not None:
            self._structural['momfile'] = momfile
        self._structural.update(kwargs)
        logger.debug(f"Set input files: pos={posfile}, ex={exchange}")
        return self
    
    # =========================================================================
    # SIMULATION SETUP (Dynamics control, passed to sim.relax())
    # =========================================================================
    
    def simulate(
        self,
        mode: str = 'M',
        temperature: Optional[float] = None,
        steps: Optional[int] = None,
        timestep: Optional[float] = None,
        damping: Optional[float] = None,
        **kwargs
    ) -> 'ConfigurationBuilder':
        """
        Set simulation parameters (dynamics control).
        
        Parameters
        ----------
        mode : str, optional
            Relaxation method ('M', 'S', or 'H'). Default: 'M'
        temperature : float, optional
            Temperature in Kelvin
        steps : int, optional
            Number of simulation steps
        timestep : float, optional
            Integration timestep (s)
        damping : float, optional
            LLG damping parameter
        **kwargs
            Other simulation parameters (ip_mode, etc.)
        
        Returns
        -------
        self
            For method chaining
        """
        if mode is not None:
            self._simulation['mode'] = mode
        if temperature is not None:
            self._simulation['temperature'] = temperature
        if steps is not None:
            self._simulation['steps'] = steps
        if timestep is not None:
            self._simulation['timestep'] = timestep
        if damping is not None:
            self._simulation['damping'] = damping
        self._simulation.update(kwargs)
        logger.debug(f"Set simulation: mode={mode}, T={temperature}, steps={steps}")
        return self
    
    # =========================================================================
    # MEASUREMENT SETUP (Output control, passed to sim.measure())
    # =========================================================================
    
    def measure(self, **kwargs) -> 'ConfigurationBuilder':
        """
        Add measurement output flags.
        
        Parameters
        ----------
        **kwargs
            Measurement parameters (plotenergy, do_avrg, avrg_step, etc.)
            Unknown flags are passed through (logged as warnings by Simulator)
        
        Returns
        -------
        self
            For method chaining
        """
        self._config.add_measurement(**kwargs)
        logger.debug(f"Added measurement: {list(kwargs.keys())}")
        return self
    
    # =========================================================================
    # ADVANCED FEATURES (Measurement details, spectroscopy, protocols, etc.)
    # =========================================================================
    
    def add_measurement(
        self,
        preset: Optional[str] = None,
        **kwargs
    ) -> ConfigurationBuilder:
        """Add measurement configuration.
        
        Parameters
        ----------
        preset : str, optional
            Preset name ('basic', 'thermal', 'spectroscopy', 'correlations')
        **kwargs
            Measurement-specific parameter overrides
            
        Returns
        -------
        self
            For method chaining
        """
        if preset:
            self._config.add_preset(f"measurement_{preset}")
            logger.debug(f"Added measurement preset: {preset}")
        if kwargs:
            self._config.add_measurement(**kwargs)
            logger.debug(f"Added measurement overrides: {kwargs}")
        return self
    
    def add_spectroscopy(self, **kwargs) -> ConfigurationBuilder:
        """Add spectroscopy configuration.
        
        Parameters
        ----------
        **kwargs
            Spectroscopy parameters (do_ams, do_magdos, etc.)
            
        Returns
        -------
        self
            For method chaining
        """
        self._config.add_spectroscopy(**kwargs)
        logger.debug(f"Added spectroscopy config: {kwargs}")
        return self
    
    def add_protocol(
        self,
        preset: Optional[str] = None,
        **kwargs
    ) -> ConfigurationBuilder:
        """Add protocol configuration.
        
        Parameters
        ----------
        preset : str, optional
            Preset name ('temperature_sweep', '3temperature', 'spinice')
        **kwargs
            Protocol-specific parameter overrides
            
        Returns
        -------
        self
            For method chaining
        """
        if preset:
            self._config.add_preset(f"protocol_{preset}")
            logger.debug(f"Added protocol preset: {preset}")
        if kwargs:
            self._config.add_protocol(**kwargs)
            logger.debug(f"Added protocol overrides: {kwargs}")
        return self
    
    def add_field(
        self,
        preset: Optional[str] = None,
        **kwargs
    ) -> ConfigurationBuilder:
        """Add field configuration.
        
        Parameters
        ----------
        preset : str, optional
            Preset name ('microwave')
        **kwargs
            Field-specific parameter overrides
            
        Returns
        -------
        self
            For method chaining
        """
        if preset:
            self._config.add_preset(f"field_{preset}")
            logger.debug(f"Added field preset: {preset}")
        if kwargs:
            self._config.add_field(**kwargs)
            logger.debug(f"Added field overrides: {kwargs}")
        return self
    
    def add_analysis(self, **kwargs) -> ConfigurationBuilder:
        """Add analysis configuration.
        
        Parameters
        ----------
        **kwargs
            Analysis parameters (do_wl, wl_nstep, etc.)
            
        Returns
        -------
        self
            For method chaining
        """
        self._config.add_analysis(**kwargs)
        logger.debug(f"Added analysis config: {kwargs}")
        return self
    
    def add_coupling(self, **kwargs) -> ConfigurationBuilder:
        """Add coupling configuration.
        
        Parameters
        ----------
        **kwargs
            Coupling parameters (do_sld, do_stt, etc.)
            
        Returns
        -------
        self
            For method chaining
        """
        self._config.add_coupling(**kwargs)
        logger.debug(f"Added coupling config: {kwargs}")
        return self
    
    def build(self) -> 'ConfigurationBuilder':
        """
        Validate and finalize the configuration.
        
        Returns
        -------
        self
            For continued method chaining after build
            
        Raises
        -------
        ValueError
            If validation fails
        """
        issues = self._config.validate()
        if issues:
            error_msg = "\n".join(issues)
            logger.error(f"Configuration validation failed:\n{error_msg}")
            raise ValueError(f"Configuration validation failed:\n{error_msg}")
        logger.info("Configuration built and validated successfully")
        return self
    
    # =========================================================================
    # EXTRACTION METHODS (for Simulator API and file writing)
    # =========================================================================
    
    def relax_kwargs(self) -> Dict[str, Any]:
        """
        Extract kwargs suitable for Simulator.relax().
        
        Returns only simulation parameters; measurement flags and structural
        parameters are filtered out.
        
        Returns
        -------
        kwargs : dict
            Arguments ready for sim.relax(**kwargs)
        """
        kwargs = {}
        
        # Map aliases: 'temperature' or 'temp', 'steps' or 'nstep', etc.
        param_map = {
            'mode': self._simulation.get('mode', 'M'),
            'temperature': self._simulation.get('temperature') or self._simulation.get('temp'),
            'steps': self._simulation.get('steps') or self._simulation.get('nstep'),
            'timestep': self._simulation.get('timestep') or self._simulation.get('delta_t'),
            'damping': self._simulation.get('damping') or self._simulation.get('mplambda1'),
            'ip_mode': self._simulation.get('ip_mode'),
        }
        
        # Add non-None values
        for key, value in param_map.items():
            if value is not None:
                if key == 'temperature' and value is None:
                    kwargs['temperature'] = 0.0
                elif value is not None:
                    kwargs[key] = value
        
        # Ensure temperature has a default
        if 'temperature' not in kwargs:
            kwargs['temperature'] = 0.0
        
        return kwargs
    
    def measure_kwargs(self) -> Dict[str, Any]:
        """
        Extract kwargs suitable for Simulator.measure().
        
        Returns only known measurement flags; simulation and structural
        parameters are filtered out.
        
        Returns
        -------
        kwargs : dict
            Arguments ready for sim.measure(**kwargs)
        """
        flat = self._config.to_fortran_kwargs()
        kwargs = {}
        
        for key, value in flat.items():
            if key in self.MEASUREMENT_FLAGS:
                kwargs[key] = value
        
        return kwargs
    
    def to_inpsd_dict(self) -> Dict[str, Any]:
        """
        Extract full configuration dict for write_inpsd_file().
        
        Combines structural, simulation, and measurement parameters into
        a single flat dict suitable for inpsd.dat writing.
        
        Returns
        -------
        config : dict
            Full configuration ready for write_inpsd_file()
        """
        result = {}
        
        # Add structural parameters
        result.update(self._structural)
        
        # Add simulation parameters
        result.update(self._simulation)
        
        # Add measurement flags (from ConfigurationManager)
        flat_cfg = self._config.to_fortran_kwargs()
        result.update(flat_cfg)
        
        return result
    
    def __repr__(self) -> str:
        """String representation."""
        n_struct = len(self._structural)
        n_sim = len(self._simulation)
        n_meas = len(self._config.config.get('measurement', {}))
        return f"ConfigurationBuilder(struct={n_struct}, sim={n_sim}, meas={n_meas})"