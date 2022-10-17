<h1>UppASD</h1>

<b>Upp</b>sala <b>A</b>tomistic <b>S</b>pin <b>D</b>ynamics software
<!--![logo][logo]-->


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
Quick-start:
```python
./setup_UppASD.sh   
make <profile>
```

Where `<profile>` is a suitable compiler profile, i.e. `ifort`, `gfortran`, `gfortran-osx` and so on.   

---
Quick-start alternative (using cmake):
```python
cmake -S . -B build   
cmake --build build
```

---
The binary is compiled to the `./bin/` directory. 

Examples are provided in `./examples/`

The code is documented in the [UppASD manual](https://uppasd.github.io/UppASD-manual/).

A tutorial with examples and exercises on atomistic spin-dynamics are contained in the [UppASD tutorial](https://uppasd.github.io/UppASD-tutorial/).

**Developers please look at the development guidelines in the `CONTRIBUTING.md` file, about how to make your contributions to UppASD.**

---

<h2>Graphical User Interface</h2>

---

A `python` based `QT` GUI for the code is also available at `./ASD_GUI/ASD_GUI.py`. 
This allows for:
- Visualization of outputs via `VTK`.
- Plotting of several quantities via integrated `matplotlib` functionalities.
- Automatic generation of input files for `UppASD`.

***Requirements***
- `Qt5`.
- `python3.6` or higher.
   - `pandas`
   - `numpy`
   - `PyYaml`
   - `matplotlib`
- `VTK7.0` or higher.

***Installation Guide (pip)***

The recommended way to install the prerequisites for the `ASD_GUI` is currently to use `pip` and `virtualenv` environments.

After installing `virtualenv` one can create virtual environment where to host the `ASD_GUI`. This can be done in the following way:

```
pip install virtualenv
virtualenv ASD_GUI_env 
source ASD_GUI_env/bin/activate
pip install numpy scipy matplotlib pandas pyqt5 vtk
```
This will generate a virtual environment named `ASD_GUI_env` which can be activated or deactivated to run the GUI.


***Installation Guide (anaconda)***

An alternative way to install the prerequisites for the `ASD_GUI` can be done via the `anaconda` framework and its environments.

After installing anaconda one can create virtual environment where to host the `ASD_GUI`. This can be done in the following way:

```
conda create --name ASD_GUI_env python vtk numpy scipy matplotlib yaml pyyaml pandas pyqt
source activate ASD_GUI_env
```
This will generate a conda environment named `ASD_GUI_env` which can be activated or deactivated to run the GUI.

---
(C) 2008-2022 [UppASD group][2]

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
