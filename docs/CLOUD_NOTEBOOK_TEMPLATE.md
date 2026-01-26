# Cloud Notebook Template

This is a template for creating UppASD notebooks compatible with both local execution, Binder, and Google Colab.

## Usage

Add this cell as the **first code cell** in any UppASD notebook:

```python
# ========================================
# CLOUD SETUP (Binder & Google Colab)
# ========================================

import sys
from pathlib import Path

# Detect environment
try:
    import google.colab
    IN_COLAB = True
    ENV = "Google Colab"
except ImportError:
    IN_COLAB = False
    # Check if in Binder (has BINDER_SERVICE_HOST env var)
    import os
    if 'BINDER_SERVICE_HOST' in os.environ:
        ENV = "Binder"
    else:
        ENV = "Local"

print(f"🌍 Environment: {ENV}")

# Install UppASD if in Colab
if IN_COLAB:
    print("\n📦 Installing UppASD from GitHub...")
    
    # Clone repository
    !git clone --quiet https://github.com/UppASD/UppASD.git /content/UppASD 2>&1 | tail -5
    
    # Build and install
    import os
    os.chdir('/content/UppASD')
    
    print("🔨 Building UppASD (this takes 2-3 minutes)...")
    !mkdir -p build && cd build && cmake .. -DCMAKE_BUILD_TYPE=Release >/dev/null 2>&1 && make -j2 >/dev/null 2>&1
    
    print("📥 Installing Python interface...")
    !pip install -e . --quiet
    
    # Install optional dependencies
    !pip install ase --quiet
    
    os.chdir('/content')
    
    print("✅ UppASD installed successfully!\n")
else:
    print("✅ Using existing UppASD installation\n")

# Standard imports
import numpy as np
from uppasd import notebook as nb
from pathlib import Path
import matplotlib.pyplot as plt

# Verify installation
try:
    from uppasd import Simulator
    print("✅ UppASD Simulator available")
except ImportError:
    print("⚠️  Simulator not available (Fortran backend may not be compiled)")
    
print("✅ Helper functions available:", len([x for x in dir(nb) if not x.startswith('_')]), "functions")
print("\n" + "="*50)
print("Ready to begin!")
print("="*50 + "\n")
```

## For Notebooks with Visualization

If using PyVista for 3D visualization:

```python
# Additional setup for PyVista in Colab
if IN_COLAB:
    try:
        import pyvista as pv
        pv.set_jupyter_backend('static')  # Static images for Colab
        print("✅ PyVista configured for Colab (static backend)")
    except ImportError:
        print("⚠️  PyVista not installed. Install with: !pip install pyvista")
else:
    try:
        import pyvista as pv
        pv.set_jupyter_backend('client')  # Interactive for local/Binder
        print("✅ PyVista configured for interactive visualization")
    except ImportError:
        print("ℹ️  PyVista not available (optional)")
```

## For Saving Results

Add at the end of the notebook:

```python
# ========================================
# DOWNLOAD RESULTS (Google Colab)
# ========================================

def save_results(directory='./results'):
    """Package and download simulation results."""
    results_dir = Path(directory)
    
    if IN_COLAB:
        import shutil
        from google.colab import files
        
        if results_dir.exists():
            print("📦 Packaging results...")
            # Create zip archive
            archive_name = 'uppasd_results'
            shutil.make_archive(archive_name, 'zip', results_dir)
            
            # Download
            files.download(f'{archive_name}.zip')
            print(f"✅ Downloaded {archive_name}.zip")
        else:
            print(f"⚠️  Directory {results_dir} not found")
    else:
        print(f"💾 Results saved locally in: {results_dir.absolute()}")

# Uncomment to download results
# save_results('./relax_sim')
```

## Complete Example: magnetic_structure_visualization_colab.ipynb

See `notebooks/examples/magnetic_structure_visualization_colab.ipynb` for a complete working example with:

1. ✅ Cloud detection and setup
2. ✅ UppASD installation (Colab)
3. ✅ Standard notebook workflow
4. ✅ Result download helper
5. ✅ PyVista compatibility

## Testing Your Cloud Notebook

### Local Test
```bash
jupyter notebook your_notebook.ipynb
# Should work without modifications
```

### Binder Test
1. Push to GitHub
2. Go to https://mybinder.org
3. Enter repo URL
4. Wait for build
5. Run notebook

### Colab Test
1. Upload notebook to Colab
2. Run first cell (installation)
3. Wait 2-3 minutes
4. Continue with rest of notebook

## Best Practices

1. **Keep first cell focused** on environment setup only
2. **Show progress** during installation (users get anxious!)
3. **Test all three environments** (local, Binder, Colab)
4. **Document expected runtimes** for cloud users
5. **Provide download helpers** for Colab users
6. **Use reasonable system sizes** (cloud has memory limits)
7. **Add markdown cells** explaining cloud-specific behavior

## Troubleshooting

**Installation fails in Colab:**
- Check GitHub URL is correct
- Ensure repository is public
- Try clearing runtime and restarting

**"Module not found" after installation:**
- Restart runtime in Colab
- Check installation actually completed
- Verify with: `!pip show uppasd`

**Visualization doesn't work in Colab:**
- Use PyVista with 'static' backend
- Or use matplotlib instead
- Check display settings

**Out of memory:**
- Reduce supercell size (n1, n2, n3)
- Decrease number of simulation steps
- Use GPU runtime in Colab if available
