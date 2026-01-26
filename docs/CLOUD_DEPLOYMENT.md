# Deploying UppASD Notebooks to the Cloud

This guide covers deploying UppASD notebooks to **Binder** and **Google Colab** for easy cloud-based access.

## Quick Launch Links

Once configured, users can launch notebooks with these badges:

**Binder:**
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/YOUR_USERNAME/UppASD/main?filepath=notebooks)

**Google Colab:**
[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/YOUR_USERNAME/UppASD/blob/main/notebooks/magnetic_structure_visualization.ipynb)

---

## 1. Dependency Management

### For Binder (Recommended for UppASD)

Binder builds a Docker container from your repository. It supports multiple configuration files:

#### Option A: requirements.txt (Simple Python dependencies)
```txt
numpy>=1.20
matplotlib>=3.5
scipy>=1.7
ase>=3.22  # Optional, for neighbor finding
```

#### Option B: environment.yml (Conda environment - RECOMMENDED)
```yaml
name: uppasd
channels:
  - conda-forge
  - defaults
dependencies:
  - python=3.11
  - numpy
  - scipy
  - matplotlib
  - gfortran  # For compiling Fortran code
  - cmake
  - pip
  - pip:
    - ase>=3.22
```

#### Option C: setup.py (Install UppASD package)
Binder will automatically run `pip install -e .` if it finds setup.py in the repo root.

**Recommended approach for UppASD:**
- Use `environment.yml` to set up conda environment with Fortran compiler
- Include `postBuild` script to compile and install UppASD

### For Google Colab

Colab comes with many packages pre-installed. Add installation cells at the top of each notebook:

```python
# First cell - Install UppASD from GitHub
!pip install git+https://github.com/YOUR_USERNAME/UppASD.git

# Or if using development version
!git clone https://github.com/YOUR_USERNAME/UppASD.git
%cd UppASD
!pip install -e .
%cd ..
```

---

## 2. File System Approach

### Current Approach: Intermediate Files ✅ **VIABLE**

**Your current file-based approach is perfectly fine for cloud/containers:**

```python
# This works in Binder and Colab
sim_dir = Path('./my_simulation')
sim_dir.mkdir(exist_ok=True)

nb.write_inpsd_file(sim_dir / 'inpsd.dat', config)
nb.write_posfile(sim_dir / 'posfile.dat', basis)
# ... run simulation
```

**Why it works:**
- Both Binder and Colab provide **ephemeral file systems**
- Files persist during the notebook session
- Users can download files before session ends
- Standard practice for simulation workflows

**Benefits:**
- ✅ Debuggable: Users can inspect input files
- ✅ Compatible: Works with existing UppASD file parsers
- ✅ Transparent: Clear what's being passed to simulator
- ✅ Flexible: Easy to modify files manually if needed

### Alternative: Pure API Approach (Optional)

For even cleaner cloud notebooks, you could add an API-only mode:

```python
# No intermediate files - everything in memory
from uppasd import Simulator

with Simulator.from_config(
    lattice=np.eye(3),
    basis=basis_array,
    exchange_matrix=exchange_dict,
    moments=moments_array,
    protocol='mc_sd'
) as sim:
    final_moments = sim.relax(temp=50, steps=2000)
```

**When to use each approach:**
- **File-based** (current): Better for teaching, debugging, compatibility
- **API-based**: Better for high-throughput, automation, minimal I/O

**Recommendation:** Keep file-based approach as primary, optionally add API shortcuts later.

---

## 3. Binder Configuration

### Step 1: Create Configuration Files

Create these files in your repository root:

#### `binder/environment.yml`
```yaml
name: uppasd-notebooks
channels:
  - conda-forge
  - defaults
dependencies:
  - python=3.11
  - numpy>=1.20
  - scipy>=1.7
  - matplotlib>=3.5
  - jupyter
  - notebook
  - ipywidgets  # For interactive widgets
  - gfortran_linux-64  # Fortran compiler (Linux)
  - cmake>=3.18
  - make
  - pip
  - pip:
    - ase>=3.22  # For neighbor finding
    - pyvista>=0.38  # For 3D visualization (optional)
```

#### `binder/postBuild` (executable script)
```bash
#!/bin/bash
# This script runs after environment is created

set -e  # Exit on error

echo "Building UppASD..."

# Compile UppASD
mkdir -p build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make -j2  # Use 2 cores
cd ..

# Install Python interface
pip install -e . --no-deps

echo "UppASD installation complete!"
```

Make it executable:
```bash
chmod +x binder/postBuild
```

#### `binder/apt.txt` (if needed for system packages)
```txt
gfortran
cmake
build-essential
```

### Step 2: Test Binder Build

1. Push configuration files to GitHub
2. Go to https://mybinder.org
3. Enter your repo URL: `https://github.com/YOUR_USERNAME/UppASD`
4. Click "launch"
5. Wait for build (first time: 5-10 minutes, subsequent: seconds)

### Step 3: Add Launch Badges to README

```markdown
## Try in the Cloud

### Interactive Notebooks on Binder
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/YOUR_USERNAME/UppASD/main?filepath=notebooks)

Launch a live Jupyter environment with all examples pre-configured.

### Individual Notebooks
- [Magnetic Structure Visualization](https://mybinder.org/v2/gh/YOUR_USERNAME/UppASD/main?filepath=notebooks/magnetic_structure_visualization.ipynb)
- [Multi-Species Setup](https://mybinder.org/v2/gh/YOUR_USERNAME/UppASD/main?filepath=notebooks/multi_species_interactive_setup.ipynb)
- [Interactive Setup](https://mybinder.org/v2/gh/YOUR_USERNAME/UppASD/main?filepath=notebooks/interactive_setup_notebook.ipynb)
```

---

## 4. Google Colab Configuration

### Step 1: Add Installation Cell to Each Notebook

Add this as the **first code cell** in each notebook:

```python
# === GOOGLE COLAB SETUP ===
# Detect if running in Colab
try:
    import google.colab
    IN_COLAB = True
    print("🌐 Running in Google Colab")
except ImportError:
    IN_COLAB = False
    print("💻 Running locally")

if IN_COLAB:
    # Install UppASD
    print("Installing UppASD from GitHub...")
    !git clone https://github.com/YOUR_USERNAME/UppASD.git /content/UppASD
    
    # Build and install
    %cd /content/UppASD
    !mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release && make -j2
    !pip install -e . --quiet
    
    # Install optional dependencies
    !pip install ase pyvista --quiet
    
    %cd /content
    print("✅ UppASD installed successfully!")
else:
    print("✅ Using local UppASD installation")

# Continue with normal imports
import numpy as np
from uppasd import notebook as nb
from pathlib import Path
```

### Step 2: Handle File Persistence in Colab

Colab files are ephemeral but can be saved to Google Drive:

```python
if IN_COLAB:
    # Optional: Mount Google Drive for persistent storage
    from google.colab import drive
    drive.mount('/content/drive')
    
    # Save results to Drive
    output_dir = Path('/content/drive/MyDrive/UppASD_Results')
    output_dir.mkdir(exist_ok=True)
else:
    output_dir = Path('./results')
    output_dir.mkdir(exist_ok=True)
```

### Step 3: Add Colab Launch Badges

Replace `YOUR_USERNAME` and `YOUR_REPO` in README:

```markdown
## Google Colab

Open notebooks directly in Google Colab (no installation required):

| Notebook | Launch |
|----------|--------|
| Magnetic Structure Visualization | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/YOUR_USERNAME/UppASD/blob/main/notebooks/magnetic_structure_visualization.ipynb) |
| Multi-Species Interactive Setup | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/YOUR_USERNAME/UppASD/blob/main/notebooks/multi_species_interactive_setup.ipynb) |
| Interactive Setup | [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/YOUR_USERNAME/UppASD/blob/main/notebooks/interactive_setup_notebook.ipynb) |
```

---

## 5. Best Practices

### A. Notebook Structure for Cloud Deployment

Structure notebooks with cloud-friendly patterns:

```python
# Cell 1: Cloud setup (Colab detection and installation)
# Cell 2: Standard imports
# Cell 3: Configuration (easy for users to modify)
# Cell 4+: Analysis workflow
# Last cell: Save/download results
```

### B. Resource Management

**Memory considerations:**
```python
# Keep supercell sizes reasonable for cloud
n1, n2, n3 = 4, 4, 3  # Good for Binder (192 atoms)
# n1, n2, n3 = 10, 10, 10  # Too large for free tier
```

**Computation time:**
```python
# Shorter simulations for cloud demos
config = nb.create_relaxation_protocol(
    'mc_sd',
    mc_steps=1000,  # Reduced from 2000
    sd_steps=1000   # Keep under 5 min for cloud
)
```

### C. Data Download Helpers

Add download functions for Colab users:

```python
def download_results(directory):
    """Package and download results."""
    import shutil
    from google.colab import files
    
    if IN_COLAB:
        # Create zip archive
        shutil.make_archive('simulation_results', 'zip', directory)
        files.download('simulation_results.zip')
        print("📦 Results downloaded!")
    else:
        print(f"💾 Results saved locally in {directory}")

# At end of notebook
download_results(sim_dir)
```

### D. Visualization Compatibility

**For PyVista (3D viz) in cloud:**
```python
if IN_COLAB:
    # Use static backend for Colab
    import pyvista as pv
    pv.set_jupyter_backend('static')
else:
    # Use interactive backend locally
    import pyvista as pv
    pv.set_jupyter_backend('client')  # or 'trame'
```

### E. Progress Indicators

Add progress feedback for long operations:

```python
from IPython.display import display, Markdown

display(Markdown("### ⏳ Building UppASD..."))
# ... installation
display(Markdown("### ✅ Ready to simulate!"))
```

---

## 6. Testing Your Cloud Setup

### Test Checklist

- [ ] Binder builds successfully (< 10 minutes)
- [ ] All dependencies install correctly
- [ ] UppASD compiles without errors
- [ ] Notebooks run start to finish
- [ ] File I/O works correctly
- [ ] Visualizations display properly
- [ ] Results can be downloaded
- [ ] Works in both Binder and Colab
- [ ] Links in README work

### Common Issues and Solutions

**Issue:** Binder build timeout
- **Solution:** Simplify dependencies, use mamba (faster than conda)
- Add to `environment.yml`: `- mamba` at top level

**Issue:** Fortran compiler not found
- **Solution:** Add `gfortran_linux-64` to conda dependencies

**Issue:** Colab can't find UppASD after install
- **Solution:** Restart runtime or check installation path
- Use `!pip show uppasd` to verify

**Issue:** Out of memory in cloud
- **Solution:** Reduce system size, use smaller supercells

---

## 7. Advanced: Docker Container (Optional)

For full reproducibility, create a Docker image:

**`Dockerfile`:**
```dockerfile
FROM jupyter/scipy-notebook:latest

USER root

# Install Fortran and build tools
RUN apt-get update && apt-get install -y \
    gfortran \
    cmake \
    build-essential \
    && rm -rf /var/lib/apt/lists/*

USER ${NB_UID}

# Copy UppASD source
COPY --chown=${NB_UID}:${NB_GID} . /home/jovyan/UppASD/

# Build and install UppASD
WORKDIR /home/jovyan/UppASD
RUN mkdir -p build && cd build && \
    cmake .. -DCMAKE_BUILD_TYPE=Release && \
    make -j4 && \
    cd .. && \
    pip install -e .

# Install optional dependencies
RUN pip install ase pyvista

WORKDIR /home/jovyan
EXPOSE 8888

CMD ["start-notebook.sh"]
```

**Build and push:**
```bash
docker build -t username/uppasd-notebooks .
docker push username/uppasd-notebooks
```

**Use in Binder:**
Create `Dockerfile` in root, and Binder will use it automatically.

---

## 8. Recommended Setup

For **UppASD**, I recommend:

### Primary: Binder
- ✅ Better for Fortran compilation
- ✅ Full conda environment control
- ✅ Persistent during session
- ✅ Better for computational work
- ✅ Free for open-source projects

### Secondary: Google Colab
- ✅ More familiar to many users
- ✅ Better GPU access (if needed later)
- ✅ Easy Google Drive integration
- ⚠️ Trickier for Fortran compilation
- ⚠️ More setup in each notebook

### File-Based Approach
- ✅ **Keep your current approach** - it works perfectly
- ✅ Standard practice for simulation tools
- ✅ Easy debugging and transparency
- ✅ Compatible with both platforms

---

## 9. Quick Start Checklist

1. **Create `binder/` directory with:**
   - `environment.yml` (conda dependencies)
   - `postBuild` (build script)
   - `apt.txt` (system packages if needed)

2. **Update notebooks with:**
   - Colab detection cell (first cell)
   - Installation instructions in markdown
   - Reasonable resource limits

3. **Add to README:**
   - Binder launch badge
   - Colab launch badges per notebook
   - Cloud usage instructions

4. **Test:**
   - Launch on Binder (wait for build)
   - Run all notebooks end-to-end
   - Verify in Colab

5. **Document:**
   - Expected runtime
   - Resource requirements
   - How to download results

---

## Summary

**Your Questions Answered:**

1. **Dependencies:** Use `environment.yml` for Binder, installation cell for Colab
2. **File approach:** ✅ Perfectly viable - standard practice for simulations
3. **Best practices:** 
   - Start with Binder (better for Fortran)
   - Add Colab support with detection cells
   - Keep resource usage reasonable
   - Provide download helpers
   - Test thoroughly before announcing

Your current file-based workflow is cloud-ready as-is! Just need configuration files and installation cells.
