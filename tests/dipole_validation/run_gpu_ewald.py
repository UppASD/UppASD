#!/usr/bin/env python3
"""Run the CMake-defined complete-kernel GPU CTest target."""
import argparse
import subprocess
from pathlib import Path

p = argparse.ArgumentParser()
p.add_argument('--build', type=Path, default=Path('build_gpu'))
a = p.parse_args()
subprocess.run(['ctest', '--test-dir', str(a.build), '--output-on-failure', '-R', '^dipole-gpu-fft-convolution$'],
               check=True)
