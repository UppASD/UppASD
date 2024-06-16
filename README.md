
<h1>UppASD</h1>

<b>Upp</b>sala <b>A</b>tomistic <b>S</b>pin <b>D</b>ynamics software
<!--![logo][logo]-->


[![Build binaries](https://github.com/UppASD/UppASD/actions/workflows/binaries.yml/badge.svg?branch=python_pip)](https://github.com/UppASD/UppASD/actions/workflows/binaries.yml)
[![Build wheels](https://github.com/UppASD/UppASD/actions/workflows/wheels.yml/badge.svg?branch=python_pip)](https://github.com/UppASD/UppASD/actions/workflows/wheels.yml)
<!---[![build status](https://gitlab.com/UppASD/UppASD/badges/master/pipeline.svg)](https://gitlab.com/UppASD/UppASD/pipelines)--->

The `UppASD` software package is a simulation suite to study magnetization dynamics by means of the atomistic version of the Landau-Lifshitz-Gilbert (LLG) equation.

***Applications:***
- [Critical temperature][11] determination
- [Magnon spectra][7] simulations
- [Micromagnetic exchange and DMI stiffness][6] calculations
- [Current][9], [field][8], and [temperature-driven][5] dynamics

***Features:***
- [LLG][12], Metropolis, and Heat-bath algorithms
- Generalized Hamiltonian including 
	-	bilinear and [biquadratic Heisenberg][10] exchange
	-	[Dzyaloshinskii-Moriya][15] exchange 
	-	magnetocrystalline anisotropies.
- Support for [arbitrary crystal structures][14] and [random alloys][13]
- [Spin-transfer][9] and spin-orbit torques
- [Temperature gradients][4]
- [Quantum statistics][3]

Detailed information about the method can be found in   
[**Atomistic Spin Dynamics: Foundations and Applications**  
O. Eriksson et. al,  Oxford University Press 2017][1]

---
---
---
# Installation:
## Binary
Download and unpack a binary from the [Release][Release] page. Installers for Linux and Windows are also available.

### Install as Python package (beta)
Python bindings for UppASD are available for installation using `pip` as follows
```python
pip install -i https://test.pypi.org/simple/ uppasd
```
The `pip` installation also provides the binary `uppasd`.

*Note: The pre-compiled binaries are not optimized, so building from source is recommended for production usage*

## Build from source
UppASD uses `cmake` for compiling the code. With `cmake` installed, UppASD can be compiled with
```
cmake -S . -B build   
cmake --build build
```
which results in a compiled binary `uppasd` locade in the `./bin/` directory.

---
# Examples and documentation

Examples are provided in the `./examples/` folder.

The code is documented in the [UppASD manual](https://uppasd.github.io/UppASD-manual/).

A tutorial with examples and exercises on atomistic spin-dynamics are contained in the [UppASD tutorial](https://uppasd.github.io/UppASD-tutorial/).

**Developers please look at the development guidelines in the `CONTRIBUTING.md` file, about how to make your contributions to UppASD.**

---
---
---

# Graphical User Interface
A `python` based `QT` GUI, named `asd_gui`, for the code is also available in the repository. 
The GUI allows for:
- Visualization of outputs via `VTK`.
- Plotting of several quantities via integrated `matplotlib` functionalities.
- Automatic generation of input files for `UppASD`.

## Installation Guide (pip)

The recommended way to use `asd_gui` is to install the offical version using `pip` as below:
```
pip install asd_gui
```
This will install the GUI as a Python module that is started by issuing the `asd_gui` command at the command prompt.


## Installation Guide (local)

For developing purposes, the GUI can also be installed from source using `pip`
```
cd ASD_GUI
pip install -e .
```
This install the GUI as an **editable** Python module that is started by issuing the `asd_gui` command at the command prompt.


---
(C) 2008-2024 [UppASD group][2]

[1]:https://global.oup.com/academic/product/atomistic-spin-dynamics-9780198788669
[2]:http://www.physics.uu.se/research/materials-theory/ongoing-research/uppasd/
[3]:https://journals.aps.org/prmaterials/abstract/10.1103/PhysRevMaterials.2.013802
[4]:https://journals.aps.org/prb/abstract/10.1103/PhysRevB.90.014434
[5]:https://www.nature.com/articles/ncomms12430
[6]:https://journals.aps.org/prb/abstract/10.1103/PhysRevB.92.214424
[7]:http://iopscience.iop.org/article/10.1088/0953-8984/27/24/243202/meta
[8]:https://journals.aps.org/prb/abstract/10.1103/PhysRevB.86.224401
[9]:https://www.nature.com/articles/srep25685
[10]:https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.111.127204
[11]:https://journals.aps.org/prb/abstract/10.1103/PhysRevB.93.214439
[12]:http://iopscience.iop.org/article/10.1088/0953-8984/20/31/315203
[13]:https://journals.aps.org/prb/abstract/10.1103/PhysRevB.94.214410
[14]:https://journals.aps.org/prb/abstract/10.1103/PhysRevB.92.094411
[15]:https://www.nature.com/articles/ncomms5815
[logo]:https://github.com/UppASD/UppASD/blob/master/docs/uppasd_rot.png
[Release]:https://github.com/UppASD/UppASD/releases