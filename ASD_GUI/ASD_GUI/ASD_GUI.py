""" @package ASD_GUI
Main executable for the UppASD GUI. The GUI is designed to allow for pre-processing
and post-processing activities with ease.
The GUI has the following capabilities
    - VTK Visualization:
        - Restart files (Final snapshot of a simulation)
        - Moment files (Animation/video rendering of the data)
            - If one has a KMC file one can see the movement of the KMC particles.
            - If there is a cluster file one can see the location of the impurity cluster atoms.
        - Neighbour setup
            - Allows to see the Heisenberg exchange neighbours in real space, with colormapping
            being determined by the magnitude of the interactions.
            - Plot the DMI vectors between neighbours, with the colormap being determined
            by the magnitude of the vectors.
        - Site dependent energy, allowing for time-dependent rendering.
        - It also includes several general actors, such a colorbars, axes and time step visualization
        - It allows for on the fly modification of the glyphs for the magnetization,
        visualization of the magnetization density.
        - On the fly change of the size of the magnetization glyphs.
    - Matplotlib Plots:
        - Spin-spin correlation function.
            - Plots the colormap when the dynamical correlation function is used.
            - Plot the AMS obtained from linear spin wave theory.
        - Magnetization averages.
            - One can turn on/off component at will.
        - Single atom  trajectories in the unit sphere.
    - Input file generation
        - One can setup the basic features to create the inpsd.dat
        - Visual creation of the posfile and momfile
        - Creation of magnetic configurations for the restarfiles
            - Domain walls.
            - Skyrmions.
            - Helical spin spirals.
Author
----------
Jonathan Chico
"""
import sys
from PyQt6.QtWidgets import QApplication
from ASD_GUI.UI.ASDUIDriver import UppASDVizMainWindow

################################################################################
## @brief Main executable class to run the ASD_Visualizer
# @details It calls the ASDUIDriver which contains the UppASDVizMainWindow class
# containing the wrapper class defining the data needed to setup the GUI.
# @author Jonathan Chico
################################################################################
def main():
    # Open the Application Window
    app = QApplication(sys.argv)
    window = UppASDVizMainWindow()
    window.show()
    window.iren.Initialize() # Need this line to actually show the render inside Qt
    # Return
    sys.exit(app.exec())

if __name__ == '__main__':
    main()
