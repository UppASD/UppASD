"""
File I/O utilities for reading UppASD output files.

Provides convenient wrappers for reading standard UppASD output files
(averages, energy, trajectories, etc.) into Python data structures.

Inspired by the plotting routines in ASD_GUI/PLOT/ASDPlotsReading.py,
adapted for direct Python use in notebooks and analysis scripts.

Classes
-------
UppASDReader
    Main reader class for parsing UppASD output files.

Functions
---------
read_averages(filename)
    Read magnetization and energy averages.
read_energy(filename)
    Read total energy data.
read_trajectories(filenames)
    Read single-atom moment trajectories.
read_yaml(filename)
    Read UppASD YAML metadata file.

Examples
--------
>>> from uppasd.fileio import UppASDReader
>>> reader = UppASDReader(simid='mysim')
>>> avg_data = reader.read_averages()
>>> print(avg_data['time'])
>>> print(avg_data['magnetization'])
"""

import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple
import logging
import glob
import yaml

logger = logging.getLogger(__name__)


class UppASDReader:
    """
    Convenient reader for UppASD output files.
    
    Automatically discovers and reads standard output files based on
    simulation ID and directory.
    
    Parameters
    ----------
    simid : str, optional
        Simulation ID used in output file names (e.g., "mysim" for "averages.mysim.out")
        If not provided, will attempt to find YAML metadata file.
    workdir : str or Path, optional
        Working directory where output files are located.
        Default: current directory
    timestep : float, optional
        Integration timestep in seconds. Used to convert iterations to time.
        Default: 1.0 (use as-is if unknown)
    
    Attributes
    ----------
    simid : str
        Simulation identifier
    workdir : Path
        Working directory
    timestep : float
        Integration timestep
    
    Examples
    --------
    **Create reader with known simid:**
    
    >>> reader = UppASDReader(simid='relax_001')
    >>> avg = reader.read_averages()
    >>> ene = reader.read_energy()
    
    **Auto-detect from YAML:**
    
    >>> reader = UppASDReader()  # Looks for uppasd.*.yaml
    >>> files = reader.available_files()
    >>> print(files)
    """
    
    def __init__(
        self,
        simid: Optional[str] = None,
        workdir: Optional[str] = None,
        timestep: float = 1.0,
    ):
        """Initialize UppASD file reader."""
        self.workdir = Path(workdir) if workdir else Path.cwd()
        self.timestep = timestep
        self.simid = simid
        
        # Try to auto-detect simid and timestep from YAML if not provided
        if not self.simid:
            yaml_files = glob.glob(str(self.workdir / "uppasd.*.yaml"))
            if yaml_files:
                try:
                    self._load_from_yaml(yaml_files[0])
                except Exception as e:
                    logger.warning(f"Could not load YAML {yaml_files[0]}: {e}")
        
        if self.simid:
            logger.debug(f"UppASDReader initialized: simid={self.simid}, "
                        f"timestep={self.timestep}s")
        else:
            logger.debug(f"UppASDReader initialized: no simid (will search by pattern)")
    
    def _load_from_yaml(self, yaml_file: str) -> None:
        """Load simid and timestep from YAML metadata file."""
        with open(yaml_file, 'r', encoding='utf-8') as f:
            data = yaml.safe_load(f)
        
        if data and 'simid' in data:
            self.simid = data['simid']
        
        if data and 'siminfo' in data and 'timestep' in data['siminfo']:
            self.timestep = float(data['siminfo']['timestep'])
        
        logger.debug(f"Loaded from YAML: simid={self.simid}, "
                    f"timestep={self.timestep}s")
    
    def available_files(self) -> Dict[str, List[str]]:
        """
        List all detected UppASD output files in the working directory.
        
        Returns
        -------
        files : dict
            Dictionary with file categories as keys and file paths as values.
            Example: {
                'averages': ['averages.mysim.out'],
                'energy': ['totenergy.mysim.out'],
                'trajectory': ['trajectory.mysim.1.out', 'trajectory.mysim.2.out'],
                ...
            }
        """
        files = {}
        
        patterns = {
            'averages': f"averages.{self.simid}.out" if self.simid else "averages.*.out",
            'energy': f"totenergy.{self.simid}.out" if self.simid else "totenergy.*.out",
            'moments': f"moments.{self.simid}.out" if self.simid else "moments.*.out",
            'trajectory': f"trajectory.{self.simid}.*.out" if self.simid else "trajectory.*.out",
            'sqw': f"sqw.{self.simid}.out" if self.simid else "sqw.*.out",
            'ams': f"ams.{self.simid}.out" if self.simid else "ams.*.out",
        }
        
        for category, pattern in patterns.items():
            found = glob.glob(str(self.workdir / pattern))
            if found:
                files[category] = sorted(found)
        
        return files
    
    def read_averages(self, filename: Optional[str] = None) -> Dict[str, np.ndarray]:
        """
        Read magnetization and energy averages.
        
        Reads the ``averages.{simid}.out`` file containing time-averaged
        magnetization components, total magnetization, and other quantities.
        
        Parameters
        ----------
        filename : str, optional
            Explicit file path. If not provided, auto-detects based on simid.
        
        Returns
        -------
        data : dict
            Dictionary with columns as keys and numpy arrays as values:
            {
                'time': ndarray,              # Time in seconds
                'magnetization_x': ndarray,
                'magnetization_y': ndarray,
                'magnetization_z': ndarray,
                'magnetization_total': ndarray,
                ...
            }
        
        Examples
        --------
        >>> reader = UppASDReader(simid='mysim')
        >>> avg = reader.read_averages()
        >>> import matplotlib.pyplot as plt
        >>> plt.plot(avg['time'], avg['magnetization_total'])
        >>> plt.xlabel('Time (s)')
        >>> plt.ylabel('|M|')
        >>> plt.show()
        """
        if not filename:
            filename = self.workdir / f"averages.{self.simid}.out"
        else:
            filename = Path(filename)
        
        if not filename.exists():
            raise FileNotFoundError(f"Averages file not found: {filename}")
        
        try:
            data_df = pd.read_csv(filename, header=0, delim_whitespace=True, 
                                 escapechar="#")
            data = {}
            
            # Convert time column (first column, iterations) to time in seconds
            iterations = data_df.iloc[:, 0].values
            data['time'] = iterations * self.timestep
            data['iterations'] = iterations
            
            # Map standard column names
            col_mapping = {
                'M_x': 'magnetization_x',
                'M_y': 'magnetization_y',
                'M_z': 'magnetization_z',
                'M_abs': 'magnetization_total',
                'E_tot': 'energy_total',
            }
            
            for col in data_df.columns[1:]:
                col_clean = col.strip()
                key = col_mapping.get(col_clean, col_clean)
                data[key] = data_df[col_clean].values
            
            logger.debug(f"Read averages from {filename}: {len(data)} columns")
            return data
            
        except Exception as e:
            logger.error(f"Error reading averages file {filename}: {e}")
            raise
    
    def read_energy(self, filename: Optional[str] = None) -> Dict[str, np.ndarray]:
        """
        Read total energy data.
        
        Reads the ``totenergy.{simid}.out`` file containing energy values
        at each sampling step.
        
        Parameters
        ----------
        filename : str, optional
            Explicit file path. If not provided, auto-detects based on simid.
        
        Returns
        -------
        data : dict
            Dictionary with keys:
            {
                'time': ndarray,              # Time in seconds
                'iterations': ndarray,        # Iteration count
                'energy': ndarray,            # Total energy
                ...                           # Other columns
            }
        """
        if not filename:
            filename = self.workdir / f"totenergy.{self.simid}.out"
        else:
            filename = Path(filename)
        
        if not filename.exists():
            raise FileNotFoundError(f"Energy file not found: {filename}")
        
        try:
            data_df = pd.read_csv(filename, header=0, delim_whitespace=True,
                                 escapechar="#")
            data = {}
            
            iterations = data_df.iloc[:, 0].values
            data['time'] = iterations * self.timestep
            data['iterations'] = iterations
            
            # Standard energy column
            if 'E_tot' in data_df.columns:
                data['energy'] = data_df['E_tot'].values
            else:
                # Fallback: use second column
                data['energy'] = data_df.iloc[:, 1].values
            
            logger.debug(f"Read energy from {filename}: {len(data['energy'])} points")
            return data
            
        except Exception as e:
            logger.error(f"Error reading energy file {filename}: {e}")
            raise
    
    def read_trajectories(
        self,
        filenames: Optional[List[str]] = None
    ) -> Dict[int, Dict[str, np.ndarray]]:
        """
        Read single-atom moment trajectories.
        
        Reads ``trajectory.{simid}.*.out`` files containing 3D moment vectors
        for individual atoms over time.
        
        Parameters
        ----------
        filenames : list of str, optional
            Explicit list of trajectory files. If not provided, auto-detects
            from trajectory.{simid}.*.out pattern.
        
        Returns
        -------
        trajectories : dict
            Dictionary mapping atom index to trajectory data:
            {
                1: {'time': ndarray, 'm_x': ndarray, 'm_y': ndarray, 'm_z': ndarray},
                2: {...},
                ...
            }
        
        Examples
        --------
        >>> reader = UppASDReader(simid='mysim')
        >>> trajs = reader.read_trajectories()
        >>> atom_1 = trajs[1]
        >>> print(atom_1['m_x'].shape)  # (n_steps,)
        """
        if not filenames:
            pattern = self.workdir / f"trajectory.{self.simid}.*.out"
            filenames = sorted(glob.glob(str(pattern)))
        
        if not filenames:
            logger.warning(f"No trajectory files found matching pattern")
            return {}
        
        trajectories = {}
        
        for traj_file in filenames:
            try:
                # Extract atom index from filename
                # e.g., "trajectory.mysim.1.out" -> atom_id = 1
                fname = Path(traj_file).name
                parts = fname.split('.')
                if len(parts) >= 3:
                    try:
                        atom_id = int(parts[2])
                    except ValueError:
                        atom_id = len(trajectories) + 1
                else:
                    atom_id = len(trajectories) + 1
                
                # Read trajectory data (columns: time/iteration, m_x, m_y, m_z)
                data = pd.read_csv(traj_file, header=None, delim_whitespace=True,
                                  usecols=[0, 2, 3, 4])
                
                iterations = data.iloc[:, 0].values
                trajectories[atom_id] = {
                    'time': iterations * self.timestep,
                    'iterations': iterations,
                    'm_x': data.iloc[:, 1].values,
                    'm_y': data.iloc[:, 2].values,
                    'm_z': data.iloc[:, 3].values,
                }
                
                logger.debug(f"Read trajectory for atom {atom_id}: "
                            f"{len(iterations)} points from {traj_file}")
                
            except Exception as e:
                logger.warning(f"Error reading trajectory file {traj_file}: {e}")
        
        return trajectories
    
    def read_sqw(self, filename: Optional[str] = None) -> Dict[str, Any]:
        """
        Read dynamical structure factor S(q,ω).
        
        Reads the ``sqw.{simid}.out`` file.
        
        Parameters
        ----------
        filename : str, optional
            Explicit file path.
        
        Returns
        -------
        data : dict
            Dictionary with structure factor components:
            {
                'q_points': ndarray,          # Q-point values
                'frequencies': ndarray,       # Frequency values
                'Sxx': ndarray,               # S_xx(q,ω)
                'Syy': ndarray,               # S_yy(q,ω)
                'Szz': ndarray,               # S_zz(q,ω)
                'S_total': ndarray,           # Total S(q,ω)
            }
        """
        if not filename:
            filename = self.workdir / f"sqw.{self.simid}.out"
        else:
            filename = Path(filename)
        
        if not filename.exists():
            raise FileNotFoundError(f"SQW file not found: {filename}")
        
        try:
            sqw_data = np.genfromtxt(filename)
            
            # Extract q and omega dimensions from file metadata
            n_q = int(sqw_data[-1, 0])
            n_omega = int(sqw_data[-1, 4])
            
            data = {
                'q_points': sqw_data[:n_q * n_omega, 0],
                'frequencies': sqw_data[:n_q * n_omega, 4],
            }
            
            # Extract S components (columns 5 onwards)
            for i, label in enumerate(['Sxx', 'Syy', 'Szz', 'S_total']):
                if 4 + i + 1 < sqw_data.shape[1]:
                    data[label] = sqw_data[:, 4 + i + 1].reshape((n_q, n_omega))
            
            logger.debug(f"Read SQW from {filename}: "
                        f"{n_q} q-points, {n_omega} frequencies")
            return data
            
        except Exception as e:
            logger.error(f"Error reading SQW file {filename}: {e}")
            raise


# Convenience functions for direct use

def read_averages(
    simid: str,
    workdir: str = ".",
    timestep: float = 1.0,
) -> Dict[str, np.ndarray]:
    """
    Convenience function to read averages file.
    
    Equivalent to `UppASDReader(simid, workdir, timestep).read_averages()`.
    
    See `UppASDReader.read_averages` for details.
    """
    reader = UppASDReader(simid=simid, workdir=workdir, timestep=timestep)
    return reader.read_averages()


def read_energy(
    simid: str,
    workdir: str = ".",
    timestep: float = 1.0,
) -> Dict[str, np.ndarray]:
    """
    Convenience function to read energy file.
    
    Equivalent to `UppASDReader(simid, workdir, timestep).read_energy()`.
    
    See `UppASDReader.read_energy` for details.
    """
    reader = UppASDReader(simid=simid, workdir=workdir, timestep=timestep)
    return reader.read_energy()


def read_trajectories(
    simid: str,
    workdir: str = ".",
    timestep: float = 1.0,
) -> Dict[int, Dict[str, np.ndarray]]:
    """
    Convenience function to read trajectory files.
    
    Equivalent to `UppASDReader(simid, workdir, timestep).read_trajectories()`.
    
    See `UppASDReader.read_trajectories` for details.
    """
    reader = UppASDReader(simid=simid, workdir=workdir, timestep=timestep)
    return reader.read_trajectories()
