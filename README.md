<h1>UppASD</h1>

<b>Upp</b>sala <b>A</b>tomistic <b>S</b>pin <b>D</b>ynamics software
<!--![logo][logo]-->


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

The binary is compiled as `./source/sd`

Examples are provided in `./examples_revision_controlled/`

The manual is found at `./docs/UppASDmanual.pdf` 

---

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
   - `PYYaml`
   - `matplotlib`
   - `Enum`
- `VTK7.0` or higher.

---
***Downloading the source files for the GUI***

If you have downloaded the full UppASD repository, the GUI is accessible from the `./ASD_GUI/` directory within the repo. 

Alternatively it can also be downloaded stand-alone from https://github.com/UppASD/UppASD/releases/download/v5.1.0/ASD_GUI.tar.gz

---
***Installation Guide using pip (recommended)***

One way to install the prerequisites for the `ASD_GUI` can be done via `pip` and `virtualenv` or `venv` environments.

After installing `virtualenv` or `venv` one can create a virtual environment where to host the `ASD_GUI`. This can be done in the following way:

```
python3 -m pip install --user --upgrade pip
python3 -m pip install --user venv
python3 -m venv ~/ASD_GUI_env
source ~/ASD_GUI_env/bin/activate
python3 -m pip install numpy scipy matplotlib pandas pyqt5 vtk
```
This will generate a virtual environment named `ASD_GUI_env` which can be activated or deactivated to run the GUI.

---
***Installation Guide using anaconda***

One simplified way to install the prerequisites for the `ASD_GUI` can be done via the `anaconda` framework and its environments.

After installing anaconda one can create a virtual environment where to host the `ASD_GUI`. This can be done in the following way:

```
conda create --name ASD_GUI_env python
source activate ASD_GUI_env
conda install vtk numpy scipy matplotlib yaml pyyaml pandas pyqt

```
This will generate a virtual environment named `ASD_GUI_env` which can be activated or deactivated to run the GUI. 

Alternatively, the following variant of conda environment can be used.

```
conda create --name ASD_GUI_env python vtk numpy scipy matplotlib yaml pyyaml pandas jsoncpp=1.8.3 tbb=2020.2
conda activate ASD_GUI_env
```

---
***Running the GUI***

With all dependencies installed, according to the instructions above, the GUI can be started by 

```
python <your-path-to-the-gui>/ASD_GUI/ASD_GUI.py
```

where `<your-path-to-the-gui>` is the path to where you downloaded the repo or exctracted the stand-alone version.

It is most convenient to run the GUI from the directory where the output files to be analyzed resides.

---
(C) 2008-2018 [UppASD group][2]

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
[logo]:../docs/uppasd_rot.png
