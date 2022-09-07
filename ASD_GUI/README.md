---

<h2>The UppASD Graphical User Interface</h2>

---

A `python` based `QT` GUI for the code is provided with the `UppASD` software.
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
(C) 2008-2022 [UppASD group][1]


[1]:http://www.physics.uu.se/research/materials-theory/ongoing-research/uppasd/
[logo]:../docs/uppasd_rot.png
