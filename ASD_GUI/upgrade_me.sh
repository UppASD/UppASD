#!/bin/bash
# upgrade_me.sh
echo 'This script will restore the Qt requirements from Qt5 to Qt6 for the ASD_GUI'
echo 'Use this as a last resort, a better approach is to ask your system administrator'
echo 'to install an up-to-date Python environment.'

for i in setup.py bin/asd_gui `find . -name "*.py"` ; 
do 
   sed -i  "s/Qt5/Qt6/g" $i 
done
sed -i "0,/asd_gui_qt5/s/asd_gui_qt5/asd_gui/" setup.py
echo '-----------------------------------------------'

echo 'Downgrade complete.'
echo '-----------------------------------------------'
echo 'Install the GUI as follows:'

echo '> python -m pip install virtualenv'
echo '> python -m virtualenv ASD_GUI_env'
echo '> source ASD_GUI_env/bin/activate'
echo '> pip install numpy matplotlib pyyaml pandas pyqt6 vtk'
echo '> pip install . ' 
echo '-----------------------------------------------'
echo 'If successful, you can then invoke the GUI with' 
echo '> asd_gui'
echo 'If not, please consider reporting this at https://github.com/UppASD/UppASD/issues'
