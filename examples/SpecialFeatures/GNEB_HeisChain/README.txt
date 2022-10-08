This example demonstrates how Geodesic Nudged Elastic Band (GNEB) method can be 
used to reveal the mechanisms of magnetization switching in a Heisenberg chain 
with nearest-neighbor exchange (see 'jfile' file). Uniaxial anisotropy is 
included along the Y axis (see 'kfile' file), which is why there are two stable 
states in the system. One of the stable states (all magnetic moments point 
along the Y axis) is referred to as the initial state and the other one (all 
magnetic moments point opposite to the Y axis) is referred to as the final 
state. The minimum energy path (MEP) between the states is calculated. 

After the GNEB calculation has converged, various properties based on the MEP 
can be visualized.

******************VISUALIZE INITIAL AND FINAL CONFIGURATIONS*******************
The python script 'myMovie_if.py' is included in the folder. It visualizes the 
initial and the final configurations after the energy minimization. Change the 
command in the 28th line to 'momfiles=glob.glob("moment_if.*.in")' in order to 
visualize the initial and the final configurations before the energy 
minimization. The script requires 'ReadTimeData.py' module, which is also 
included in the folder. The script is based on the code created by Anders
Bergman.


*****************VISUALIZE CONFIGURATIONS ALONG TRANSITION PATH****************
The python script 'myMovie_path.py' is included in the folder. It sequentially 
shows configurations along the MEP between the initial and the final states. 
Change the command in the 28th line to 'momfiles=glob.glob("moment_path.*.in")' 
in order to visualize the path before relaxation. If the number of images has 
been changed (presently there are 16 images, see 'inpsd.dat' file), change the 
number in the 213th line of the script accordingly. The script requires 
'ReadTimeData.py' module, which is also included in the folder. The script is
based on the code created by Anders Bergman.


*****************************PLOT ENERGY ALONG MEP*****************************
The gnuplot script 'script_ene.plo' is included in the folder. It plots energy 
of the images and interpolated energy vs reaction coordinate. The energy 
maximum (Saddle Point (SP) energy) is also indicated. The result is saved in 
the 'mep.eps' file.
