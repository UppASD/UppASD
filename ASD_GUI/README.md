<h1>asd_gui</h1>

# Graphical User Interface
A `python` based `QT` GUI, named `asd_gui`, for the code is available in the repository. 
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
***Running the GUI***

With all dependencies installed, according to the instructions above, the GUI can be started by 

```
asd_gui
```
It is most convenient to run the GUI from the directory where the output files to be analyzed resides.

---
(C) 2008-2022 [UppASD group][1]


[1]:http://www.physics.uu.se/research/materials-theory/ongoing-research/uppasd/
[logo]:../docs/uppasd_rot.png
