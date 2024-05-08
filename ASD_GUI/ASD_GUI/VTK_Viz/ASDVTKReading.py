""" @package ASDVTKReading
The ASDVTKReading clas encompasses all the reading routines necessary for the
several types of visualizations available in the ASD_visualizer.
These include:
    - Reading restart files
    - Reading moment files
    - Reading struct files
    - Reading coord files
The ASDVTKReading tries to make use of the linecache reader capabilities, which
should allow one to deal with large data sets in a much more efficient capacity
than the numpy and string options. However, they are still not at efficient as the pandas
read_csv routines, which are not included to avoid compatibility issues.
This routine also handles whether the structure is 2D or 3D, bear in mind, that 2D and 3D
in this context are not about a true 2D or 3D structure, if not what kind of density rendering
can be done, as for very thin samples, but with thickness different around 0, a 2D density rendering
approach is more appropiate than a full 3D render approach
When including any further visualization capabilities, requiring different files
it is recommended to include the reading routines inside this class.

Author
----------
Jonathan Chico
"""
################################################################################
# @brief Class dealing with the reading of data for the VTK mode.
# @ details Class for reading the data from ASD files, it makes use of the pandas framework
# to read the data needed for the visualization of magnetic moments, neighbours, etc.
# That is any data needed for the VTK visualization.
# @author Jonathan Chico
################################################################################
class ASDReading():
    ############################################################################
    # @brief Constructor for the ASDReading class.
    # @details Constructor for the ASDReading class. It contains a set of definitions
    # mainly dealing with the names of the data files.
    # @author Jonathan Chico
    ############################################################################
    def __init__(self):
        ASDReading.restart=[]
        ASDReading.posfiles=[]
        ASDReading.kmcfiles=[]
        ASDReading.enefiles=[]
        ASDReading.full_coord=[]
        ASDReading.full_mom=[]
        ASDReading.full_ene=[]
        ASDReading.full_KMC=[]
        ASDReading.structfiles=[]
        ASDReading.dmdatafiles=[]
        ASDReading.magnetization=[]
        ASDReading.not_read_pos=True
        ASDReading.not_read_mom=True
        ASDReading.not_read_ene=True
        ASDReading.not_read_neigh=True
        ASDReading.not_read_dmneigh=True
        return
    ############################################################################
    # @brief Function to get the file names for the different types of visualizations
    # @details Function that gets the file names needed for the different types
    # of visualization in the VTK mode, namely:
    #   - coord.*.out
    #   - restart.*.out
    #   - moment.*.out
    #   - restart.*.out
    #   - struct.*.out
    #   - dmstruct.*.out
    #   - localenergy.*.out
    #   - kmc_info.*.out
    #   - clus_info.*.out
    ############################################################################
    def getFileName(self,window):
        from PyQt6 import QtWidgets

        dlg = QtWidgets.QFileDialog()
        dlg.setFileMode(QtWidgets.QFileDialog.FileMode.ExistingFile)
        dlg.setDirectory('.')
        if dlg.exec():
            if window.sender()==window.actionCoordinate_File:
                ASDReading.posfiles=dlg.selectedFiles()[0]
                ASDReading.not_read_pos=True
            if window.sender()==window.actionMagnetization_File:
                ASDReading.magnetization=dlg.selectedFiles()[0]
                ASDReading.not_read_mom=True
            if window.sender()==window.actionStruct_File:
                ASDReading.structfiles=dlg.selectedFiles()[0]
                ASDReading.not_read_neigh=True
            if window.sender()==window.actionDM_File:
                ASDReading.dmdatafiles=dlg.selectedFiles()[0]
                ASDReading.not_read_dmneigh=True
            if window.sender()==window.actionKMC_File:
                ASDReading.kmcfiles=dlg.selectedFiles()[0]
            if window.sender()==window.actionEnergy_File:
                ASDReading.enefiles=dlg.selectedFiles()[0]
                ASDReading.not_read_ene=True
        return
    #--------------------------------------------------------------------------------
    # @brief Wrapper to read all the information needed for visualization
    # @details Wrapper handling the reading of all the files needed for the
    # visualization of the data in the VTK mode.
    # @author Jonathan Chico
    #--------------------------------------------------------------------------------
    def ReadingWrapper(self,mode,viz_type,file_names,window):
        import glob
        from ASD_GUI.UI import ASDInputWindows

        ASDReading.posfiles      = file_names[0]
        ASDReading.magnetization = file_names[1]
        ASDReading.kmcfiles      = file_names[2]
        ASDReading.structfiles   = file_names[3]
        ASDReading.enefiles      = file_names[4]
        ASDReading.dmdatafiles   = file_names[5]
        ASDReading.error_trap    = False
        #----------------------------------------------------------------------------
        # This selects which type of visualization one has moments, neigbours or energy
        #----------------------------------------------------------------------------
        #----------------------------------------------------------------------------
        # Magnetization (restart, moments) type of visualization
        #----------------------------------------------------------------------------
        if viz_type=='M':
            #------------------------------------------------------------------------
            # Find the restartfile
            #------------------------------------------------------------------------
            if len(ASDReading.magnetization)>0:
                ASDReading.MagFile = open(ASDReading.magnetization)
            else:
                print("No file name selected from menu. Trying to find a 'localenergy.*.out' file")
                ASDReading.MagFiles=glob.glob("restart.*.out")
                if len(ASDReading.MagFiles)>0:
                    ASDReading.MagFile= open(ASDReading.MagFiles[0])
                    window.Res_Info_Window=ASDInputWindows.Info_Window()
                    window.Res_Info_Window.FunMsg.setText("Information: no magnetic configuration given")
                    window.Res_Info_Window.InfoMsg.setText("File "+str(ASDReading.MagFiles[0])+" chosen as default.")
                    window.Res_Info_Window.show()
                else:
                    window.Res_Error_Window=ASDInputWindows.Error_Window()
                    window.Res_Error_Window.FunMsg.setText("I'm sorry, Dave. I'm afraid I can't do that.")
                    window.Res_Error_Window.ErrorMsg.setText("Error: No magnetic configuration file found!")
                    window.Res_Error_Window.show()
                    ASDReading.error_trap=True
                    print("Error: No magnetic configuration file selected!")
                    print("I'm sorry, Dave. I'm afraid I can't do that.")
        #----------------------------------------------------------------------------
        # Neighbour type of visualization
        #----------------------------------------------------------------------------
        elif viz_type=='N':
            if mode==1:
                #--------------------------------------------------------------------
                # Find the structfile
                #--------------------------------------------------------------------
                if len(ASDReading.structfiles)>0:
                    ASDReading.structFile=ASDReading.structfiles
                else:
                    print("No file name selected from menu. Trying to find a 'struct.*.out' file")
                    ASDReading.structfiles=glob.glob("struct.*.out")
                    if len(ASDReading.structfiles)>0:
                        ASDReading.structfiles=ASDReading.structfiles[0]
                        ASDReading.structFile= ASDReading.structfiles
                    else:
                        window.Neigh_Error_Window=ASDInputWindows.Error_Window()
                        window.Neigh_Error_Window.FunMsg.setText("I'm sorry, Dave. I'm afraid I can't do that.")
                        window.Neigh_Error_Window.ErrorMsg.setText("Error: No 'struct.*.out' file found!")
                        window.Neigh_Error_Window.show()
                        ASDReading.error_trap=True
                        print("No 'struct.*.out' file found!")
                        print("I'm sorry, Dave. I'm afraid I can't do that.")
            elif mode==2:
                #---------------------------------------------------------------
                # Find the dmdatafile
                #---------------------------------------------------------------
                if len(ASDReading.dmdatafiles)>0:
                    ASDReading.DMFile=ASDReading.dmdatafiles
                else:
                    print("No file name selected from menu. Trying to find a 'dmdata.*.out' file")
                    ASDReading.dmdatafiles=glob.glob("dmdata.*.out")
                    if len(ASDReading.dmdatafiles)>0:
                        ASDReading.dmdatafiles=ASDReading.dmdatafiles[0]
                        ASDReading.DMFile= ASDReading.dmdatafiles
                    else:
                        window.DMNeigh_Error_Window=ASDInputWindows.Error_Window()
                        window.DMNeigh_Error_Window.FunMsg.setText("I'm sorry, Dave. I'm afraid I can't do that.")
                        window.DMNeigh_Error_Window.ErrorMsg.setText("Error: No 'dmdata.*.out' file found!")
                        window.DMNeigh_Error_Window.show()
                        ASDReading.error_trap=True
                        print("No 'dmdata.*.out' file found!")
                        print("I'm sorry, Dave. I'm afraid I can't do that.")
        #-----------------------------------------------------------------------
        # Energy type of visualization
        #-----------------------------------------------------------------------
        elif viz_type=='E':
            #-------------------------------------------------------------------
            # Find the restartfile
            #-------------------------------------------------------------------
            if len(ASDReading.enefiles)>0:
                ASDReading.eneFile=ASDReading.enefiles
            else:
                print("No file name selected from menu. Trying to find a 'localenergy.*.out' file")
                ASDReading.enefiles=glob.glob("localenergy.*.out")
                if len(ASDReading.enefiles)>0:
                    ASDReading.enefiles=ASDReading.enefiles[0]
                    ASDReading.eneFile= ASDReading.enefiles
                else:
                    window.Ene_Error_Window=ASDInputWindows.Error_Window()
                    window.Ene_Error_Window.FunMsg.setText("I'm sorry, Dave. I'm afraid I can't do that.")
                    window.Ene_Error_Window.ErrorMsg.setText("Error: No 'localenergy.*.out' file found!")
                    window.Ene_Error_Window.show()
                    ASDReading.error_trap=True
                    print("No 'localenergy.*.out' file found!")
                    print("I'm sorry, Dave. I'm afraid I can't do that.")
        #-----------------------------------------------------------------------
        # Find the coordinate file
        #-----------------------------------------------------------------------
        if len(ASDReading.posfiles)>0:
            atomsFile = open(ASDReading.posfiles)
        else:
            print("No file name selected from menu. Trying to find a 'coord.*.out' file")
            ASDReading.posfiles=glob.glob("coord.*.out")
            if len(ASDReading.posfiles)>0:
                ASDReading.posfiles=ASDReading.posfiles[0]
                atomsFile = ASDReading.posfiles
            else:
                window.Coord_Error_Window=ASDInputWindows.Error_Window()
                window.Coord_Error_Window.FunMsg.setText("Sorry But Our Princess is in Another Castle.")
                window.Coord_Error_Window.ErrorMsg.setText("Error: No 'coord.*.out' file found!")
                window.Coord_Error_Window.show()
                ASDReading.error_trap=True
                print("No 'coord.*.out' file found!")
                print("Sorry But Our Princess is in Another Castle.")
        #-----------------------------------------------------------------------
        # Cluster coordinates file
        #-----------------------------------------------------------------------
        posfiles_c=glob.glob("clus_info.*.out")
        if len(posfiles_c)>0:
            ASDReading.cluster_flag=True
            #-------------------------------------------------------------------
            # Open the file if it exists
            #-------------------------------------------------------------------
            atomsFile_c = open(posfiles_c[0])
        else:
            ASDReading.cluster_flag=False
        #-----------------------------------------------------------------------
        # Check if the KMC file must be read
        #-----------------------------------------------------------------------
        if len(ASDReading.kmcfiles)>0:
            ASDReading.kmc_flag=True
            ASDReading.KMCFile=ASDReading.kmcfiles
        else:
            ASDReading.kmcfiles=glob.glob("kmc_info.*.out")
            if len(ASDReading.kmcfiles)>0:
                ASDReading.kmc_flag=True
                #---------------------------------------------------------------
                # Open the kmc file if it exists
                #---------------------------------------------------------------
                ASDReading.KMCFile = ASDReading.kmcfiles[0]
            else:
                ASDReading.kmc_flag=False
        #-----------------------------------------------------------------------
        # If not read already read the coordinates file
        #-----------------------------------------------------------------------
        if ASDReading.not_read_pos and not ASDReading.error_trap:
            #-------------------------------------------------------------------
            # Actually reading the files
            #-------------------------------------------------------------------
            (ASDReading.coord,ASDReading.nrAtoms,ASDReading.flag2D,                 \
            ASDReading.min_val)=self.readAtoms(atomsFile)
            ASDReading.not_read_pos=False
            #-------------------------------------------------------------------
            # Checking if the clusters are present
            #-------------------------------------------------------------------
            if ASDReading.cluster_flag:
                #---------------------------------------------------------------
                # Setting the coordinates of the impurity cluster
                #---------------------------------------------------------------
                (ASDReading.coord_c,ASDReading.nrAtoms_c,ASDReading.colors_clus,    \
                ASDReading.points_clus_imp,ASDReading.colors_imp,ASDReading.imp_nr)=\
                self.readAtoms_clus(atomsFile_c)
        #-----------------------------------------------------------------------
        # If the type of visualization is on the moments mode read the data about it
        #-----------------------------------------------------------------------
        if viz_type=='M' and ASDReading.not_read_mom and not ASDReading.error_trap:
            # Read the data for the vectors
            (ASDReading.moments,ASDReading.colors,ASDReading.number_time_steps,     \
            ASDReading.time_sep)=\
            self.readVectorsData(ASDReading.MagFile,0,ASDReading.nrAtoms,1)
            ASDReading.not_read_mom=False
            # Check if there are KMC files present
            if ASDReading.kmc_flag:
                (ASDReading.coord_KMC,self.nrKMCpar)=                               \
                self.readKMCData(ASDReading.KMCFile,0,1)
        #-----------------------------------------------------------------------
        # If the type of visualization is about the neigbours read that data instead
        #-----------------------------------------------------------------------
        elif viz_type=='N' and ASDReading.not_read_neigh and not ASDReading.error_trap and mode==1:
            (ASDReading.neighbours, ASDReading.Neigh_strength,ASDReading.curr_atom, \
            ASDReading.neigh_types,ASDReading.num_types_total) =                    \
            self.readNeighbours(ASDReading.structFile)
            #------------------------------------------------------------------------
            # Calculate neighbours to iAtom
            #------------------------------------------------------------------------
            (ASDReading.neighs,ASDReading.atomCenter,ASDReading.nTypes,             \
            ASDReading.neigh_colors,ASDReading.nNeighs,ASDReading.num_types,        \
            ASDReading.types_counters,ASDReading.types) =                           \
            self.setNeighbours(ASDReading.neighbours,0,ASDReading.coord,            \
            ASDReading.Neigh_strength,ASDReading.curr_atom,ASDReading.neigh_types)
            ASDReading.not_read_neigh=False
        #-----------------------------------------------------------------------
        # Read the DM neighbour data
        #-----------------------------------------------------------------------
        elif viz_type=='N' and ASDReading.not_read_dmneigh and not ASDReading.error_trap and mode==2:
            (ASDReading.neighbours,ASDReading.DM_vec,ASDReading.DM_strength,        \
            ASDReading.curr_atom,ASDReading.neigh_types,ASDReading.num_types_total)=\
            self.readDMNeighbours(ASDReading.DMFile)
            #------------------------------------------------------------------------
            # Calculate the neighbours to iAtom for the DM interactions
            #------------------------------------------------------------------------
            (ASDReading.neighs,ASDReading.atomCenter,ASDReading.nTypes,             \
            ASDReading.neigh_colors,ASDReading.DM_vectors,ASDReading.nNeighs,       \
            ASDReading.num_types,ASDReading.types_counters,ASDReading.types)=\
            self.setDMNeighbours(ASDReading.neighbours,0,ASDReading.coord,          \
            ASDReading.DM_vec,ASDReading.DM_strength,ASDReading.curr_atom,          \
            ASDReading.neigh_types)
            ASDReading.not_read_dmneigh=False
        #-----------------------------------------------------------------------
        # If it is an energy type of simulation and the file has not been read, read it
        #-----------------------------------------------------------------------
        elif viz_type=='E' and ASDReading.not_read_ene and not ASDReading.error_trap:
            # Read the data for the energy
            (ASDReading.energies, ASDReading.number_time_steps,ASDReading.time_sep)=\
            self.readEnergyData(ASDReading.eneFile,0,ASDReading.nrAtoms,1)
            ASDReading.not_read_ene=False
        #-----------------------------------------------------------------------
        # Store some global parameters
        #-----------------------------------------------------------------------
        ASDReading.mode=mode
        ASDReading.viz_type=viz_type
        return
    ############################################################################
    # @brief Read Location of Atoms
    # @details Function to read the position of the atoms in the system, as well
    # as calculating if the system should be treated like a 2D or 3D system.
    # @author Jonathan Chico
    ############################################################################
    def readAtoms(self,file_coord):
        from vtk import vtkPoints
        from vtk.util import numpy_support
        import numpy as np
        import pandas as pd
        #-----------------------------------------------------------------------
        # Define the parameters used to check if the system is 2D
        #-----------------------------------------------------------------------
        tol=2.00
        ASDReading.flag2D=False
        #-----------------------------------------------------------------------
        # Define vtkPoints as this type of arrays will be used to define the grid
        #-----------------------------------------------------------------------
        points = vtkPoints()
        #-----------------------------------------------------------------------
        # Check if the coordinate file is animated (gdisp)
        #-----------------------------------------------------------------------
        tmp_c=np.genfromtxt(fname=file_coord)
        if tmp_c.shape[1]==11:
            # Animated coordinate file (gdisp)
            coord = tmp_c[:,[2,3,4]]
            ASDReading.full_coord = coord
            ASDReading.nrAtoms=np.int32(tmp_c[-1,1])
        else:
            # Assume regular coord file
            coord = tmp_c[:,[1,2,3]]
            ASDReading.nrAtoms=len(coord)
        del tmp_c
        ASDReading.full_coord = coord
        ####-----------------------------------------------------------------------
        #### Read the data with pandas
        ####-----------------------------------------------------------------------
        ###coord=np.genfromtxt(fname=file_coord,usecols=[1,2,3])
        ####-----------------------------------------------------------------------
        #### Define the number of atoms in the system
        ####-----------------------------------------------------------------------
        ###ASDReading.nrAtoms=len(coord)
        ####-----------------------------------------------------------------------
        # Pass the numpy type arrays to vtk objects
        #-----------------------------------------------------------------------
        points.SetData(numpy_support.numpy_to_vtk(coord[0:ASDReading.nrAtoms]))
        #-----------------------------------------------------------------------
        # Data to check if one should consider the data to be rendered in 2D or 3D
        #-----------------------------------------------------------------------
        max_x=max(coord[0:ASDReading.nrAtoms,0])
        min_x=min(coord[0:ASDReading.nrAtoms,0])
        max_y=max(coord[0:ASDReading.nrAtoms,1])
        min_y=min(coord[0:ASDReading.nrAtoms,1])
        max_z=max(coord[0:ASDReading.nrAtoms,2])
        min_z=min(coord[0:ASDReading.nrAtoms,2])
        dist_x=np.sqrt((max_x-min_x)**2)
        dist_y=np.sqrt((max_y-min_y)**2)
        dist_z=np.sqrt((max_z-min_z)**2)
        min_dist=min(dist_x,dist_y,dist_z)
        #-----------------------------------------------------------------------
        # If the minimum distace is small set the flag to be 2D, this does not mean the
        # system is trully 2D, if not that the distance is small enough, such that the
        # 2D delaunay algorithm is used
        #-----------------------------------------------------------------------
        if min_dist<tol:
            ASDReading.flag2D=True
        ASDReading.min_val=min_z
        #-----------------------------------------------------------------------
        # Cleanup temporary arrays
        #-----------------------------------------------------------------------
        del coord
        return points, ASDReading.nrAtoms,ASDReading.flag2D,ASDReading.min_val
    ############################################################################
    # @biref Read Location of the cluster Atoms notice that this uses the clus_info file
    # @details Function to read the positions of the atoms in the cluster,
    # this also contains information about whether one should plot the inpurity
    # center atom.
    # @author Jonathan Chico
    ############################################################################
    def readAtoms_clus(self,file_clus):
        from vtk import vtkPoints, vtkUnsignedCharArray
        import numpy as np
        import pandas as pd
        tol=0.0001
        #-----------------------------------------------------------------------
        # Define vtkPoints arrays for the creation of the grid
        #-----------------------------------------------------------------------
        points_clus = vtkPoints()
        points_clus_imp = vtkPoints()
        #-----------------------------------------------------------------------
        # Define UnsignedCharArrays for the definitions of the colors to be used
        #-----------------------------------------------------------------------
        colors_clus = vtkUnsignedCharArray()
        colors_imp = vtkUnsignedCharArray()
        #-----------------------------------------------------------------------
        # Set the size of the arrays
        #-----------------------------------------------------------------------
        colors_clus.SetNumberOfComponents(3)
        colors_imp.SetNumberOfComponents(3)
        #-----------------------------------------------------------------------
        # Read the data with pandas
        #-----------------------------------------------------------------------
        data=pd.read_csv(file_clus,skiprows=1,header=None,delim_whitespace=True,    \
        usecols=[2,3,4,5]).values
        coord=data[:,0:3]
        itype=data[:,3]
        del data
        #-----------------------------------------------------------------------
        # Define the number of atoms in the cluster
        #-----------------------------------------------------------------------
        ASDReading.nrAtoms_clus=len(coord)
        #-----------------------------------------------------------------------
        # This will ensure that one can keep track of how many impurities one
        # actually has in the selected area
        #-----------------------------------------------------------------------
        imp_nr=0
        #-----------------------------------------------------------------------
        # Find the indices the correspond to the impurity atoms
        #-----------------------------------------------------------------------
        ind_type=np.where(itype==1)
        #-----------------------------------------------------------------------
        # Pass the numpy type data to vtk data
        #-----------------------------------------------------------------------
        for ii in range(0,ASDReading.nrAtoms_clus):
            points_clus.InsertPoint(ii,coord[ii,0],coord[ii,1],coord[ii,2])
            colors_clus.InsertTuple3(ii,0,0,0)
            for jj in range(0,len(ind_type[0])):
                #---------------------------------------------------------------
                # Data to display the center of the impurity cluster
                #---------------------------------------------------------------
                dist=np.sqrt((coord[ii,0]-coord[ind_type[0][jj],0])**2+             \
                             (coord[ii,1]-coord[ind_type[0][jj],1])**2+             \
                             (coord[ii,2]-coord[ind_type[0][jj],2])**2)
                #---------------------------------------------------------------
                # This is to ensure that the impurity cluster and the embedded cluster
                # are different structures that can be used differently
                #---------------------------------------------------------------
                if dist<tol:
                    colors_imp.InsertTuple3(imp_nr,51,160,44)
                    points_clus_imp.InsertPoint(imp_nr,coord[ii,0],coord[ii,1],coord[ii,2])
                    imp_nr=imp_nr+1
                elif dist<1:
                    colors_imp.InsertTuple3(imp_nr,106,61,154)
                    points_clus_imp.InsertPoint(imp_nr,coord[ii,0],coord[ii,1],coord[ii,2])
                    imp_nr=imp_nr+1
        #-----------------------------------------------------------------------
        # Cleanup the leftover data
        #-----------------------------------------------------------------------
        del coord
        return points_clus, ASDReading.nrAtoms_clus,colors_clus,points_clus_imp,colors_imp,imp_nr
    #--------------------------------------------------------------------------------
    # @brief Read the magnetic moments vectors, for both restart and moment visualization
    # @details Read the magnetic moments vectors, for both restart and moment visualization.
    # The number of total lines will be found and used to determine the number of frames,
    # these must be then passed back to the function
    # @author Jonathan Chico
    #--------------------------------------------------------------------------------
    def readVectorsData(self,file_mom,time,nrAtoms,temp_count):
        from vtk import vtkFloatArray
        from vtk.util import numpy_support
        import numpy as np
        import pandas as pd
        #----------------------------------------------------------------------------
        # Create a Double array which represents the vectors
        #----------------------------------------------------------------------------
        vectors = vtkFloatArray()
        colors_x = vtkFloatArray()
        colors_y = vtkFloatArray()
        colors_z = vtkFloatArray()
        #----------------------------------------------------------------------------
        # Define number of elements
        #----------------------------------------------------------------------------
        vectors.SetNumberOfComponents(3)
        colors_x.SetNumberOfComponents(1)
        colors_y.SetNumberOfComponents(1)
        colors_z.SetNumberOfComponents(1)
        #----------------------------------------------------------------------------
        # Check which kind of visualization is being done restartfile or momentfile
        #----------------------------------------------------------------------------
        if time==0:
            #------------------------------------------------------------------------
            # Find the format and type of file
            #------------------------------------------------------------------------
            (type_fmt,file_type)=self.check_format(file_mom)
            #------------------------------------------------------------------------
            # Check if the file is in the new format
            #------------------------------------------------------------------------
            if type_fmt=='new':
                #ASDReading.full_mom=pd.read_csv(file_mom,header=None,               \
                #delim_whitespace=True,skiprows=7,usecols=[4,5,6]).values
                ASDReading.full_mom=np.genfromtxt(fname=file_mom,usecols=[4,5,6])
                file_mom.seek(0)
                # Find how many different "times" there are
                ASDReading.number_time_steps=len(ASDReading.full_mom)/nrAtoms
                #--------------------------------------------------------------------
                # If there are more than one time step
                #--------------------------------------------------------------------
                if ASDReading.number_time_steps>1:
                    # Read the times
                    ASDReading.time_sep=np.genfromtxt(file_mom,usecols=[0])
                    #ASDReading.time_sep=pd.read_csv(file_mom,header=None,skiprows=7,\
                    #delim_whitespace=True,usecols=[0]).values
                    # Find the separations between different times
                    ASDReading.time_sep=np.unique(ASDReading.time_sep)
                    # If there is only one time check if there are several ensembles
                    if len(ASDReading.time_sep)==1:
                        # Read the ensembles
                        file_mom.seek(0)
                        ASDReading.time_sep=np.genfromtxt(file_mom,usecols=[1])
                        #ASDReading.time_sep=pd.read_csv(file_mom,header=None,       \
                        #skiprows=7,delim_whitespace=True,usecols=[1]).values
                        # Find how many different ensembles there are
                        ASDReading.time_sep=np.unique(ASDReading.time_sep)
                elif ASDReading.number_time_steps==1:
                    ASDReading.time_sep=0
            #------------------------------------------------------------------------
            # Read the file in the old format
            #------------------------------------------------------------------------
            elif type_fmt=='old':
                #--------------------------------------------------------------------
                # Read the restartfile
                #--------------------------------------------------------------------
                if file_type=='restart':
                    # Read the restartfile
                    ASDReading.full_mom=np.genfromtxt(file_mom,skip_header=1, \
                                                      usecols=[3,4,5])
                    #ASDReading.full_mom=pd.read_csv(file_mom,header=None,skiprows=1,\
                    #delim_whitespace=True,usecols=[3,4,5]).values
                    file_mom.seek(0)
                    # Find how many different "times" there are
                    ASDReading.number_time_steps=len(ASDReading.full_mom)/nrAtoms
                    if ASDReading.number_time_steps>1:
                        # Read the ensembles
                        ASDReading.time_sep=np.genfromtxt(file_mom,usecols=[0])
                        #ASDReading.time_sep=pd.read_csv(file_mom,header=None,       \
                        #delim_whitespace=True,usecols=[0]).values
                        # Find how many different ensembles there are
                        ASDReading.time_sep=np.unique(ASDReading.time_sep)
                    elif ASDReading.number_time_steps==1:
                        ASDReading.time_sep=0
                #--------------------------------------------------------------------
                # Read the momentfile
                #--------------------------------------------------------------------
                if file_type=='moment':
                    # Read the momentfile
                    ASDReading.full_mom=pd.read_csv(file_mom,header=None,           \
                    delim_whitespace=True,usecols=[2,3,4]).values
                    file_mom.seek(0)
                    # Find how many different "times" there are
                    ASDReading.number_time_steps=len(ASDReading.full_mom)/nrAtoms
                    if ASDReading.number_time_steps>1:
                        # Read the times
                        ASDReading.time_sep=pd.read_csv(file_mom,header=None,       \
                        delim_whitespace=True,usecols=[0]).values
                        # Find the separations between different times
                        ASDReading.time_sep=np.unique(ASDReading.time_sep)
                    elif ASDReading.number_time_steps==1:
                        ASDReading.time_sep=0
        else:
            ASDReading.number_time_steps=temp_count
        #----------------------------------------------------------------------------
        # Find the boundaries
        #----------------------------------------------------------------------------
        t_off = time*nrAtoms
        ### min_x=min(ASDReading.full_mom[t_off:t_off+nrAtoms,0])
        ### min_y=min(ASDReading.full_mom[t_off:t_off+nrAtoms,1])
        ### min_z=min(ASDReading.full_mom[t_off:t_off+nrAtoms,2])
        ### max_x=max(ASDReading.full_mom[t_off:t_off+nrAtoms,0])
        ### max_y=max(ASDReading.full_mom[t_off:t_off+nrAtoms,1])
        ### max_z=max(ASDReading.full_mom[t_off:t_off+nrAtoms,2])
        #----------------------------------------------------------------------------
        # Loop over all the atoms
        #----------------------------------------------------------------------------
        #print('Setting up vtk moment arrays')
        #vectors = numpy_support.numpy_to_vtk(ASDReading.full_mom)
        #colors_x = (ASDReading.full_mom[:,0] - min_x) / ( max_x - min_x)
        #colors_y = (ASDReading.full_mom[:,1] - min_y) / ( max_y - min_y)
        #colors_z = (ASDReading.full_mom[:,2] - min_z) / ( max_z - min_z)
        vectors = numpy_support.numpy_to_vtk(ASDReading.full_mom[t_off:t_off+nrAtoms,:])
        ### colors_x = (ASDReading.full_mom[t_off:t_off+nrAtoms,0] - min_x) / ( max_x - min_x + 1.0e-12)
        ### colors_y = (ASDReading.full_mom[t_off:t_off+nrAtoms,1] - min_y) / ( max_y - min_y + 1.0e-12)
        ### colors_z = (ASDReading.full_mom[t_off:t_off+nrAtoms,2] - min_z) / ( max_z - min_z + 1.0e-12)
        ### colors_x = numpy_support.numpy_to_vtk(colors_x)
        ### colors_y = numpy_support.numpy_to_vtk(colors_y)
        ### colors_z = numpy_support.numpy_to_vtk(colors_z)
        colors_x = numpy_support.numpy_to_vtk(ASDReading.full_mom[t_off:t_off+nrAtoms,0] )
        colors_y = numpy_support.numpy_to_vtk(ASDReading.full_mom[t_off:t_off+nrAtoms,1] )
        colors_z = numpy_support.numpy_to_vtk(ASDReading.full_mom[t_off:t_off+nrAtoms,2] )
        #for ii in range(0,nrAtoms):
        #    #------------------------------------------------------------------------
        #    # Pass the data from the numpy arrays to vtk data structures
        #    #------------------------------------------------------------------------
        #    #vectors.InsertTuple3(ii,ASDReading.full_mom[time*(nrAtoms)+ii,0],\
        #    #    ASDReading.full_mom[time*(nrAtoms)+ii,1],\
        #    #    ASDReading.full_mom[time*(nrAtoms)+ii,2])
        #    if ASDReading.flag2D:
        #        colors_x.InsertValue(ii,ASDReading.full_mom[time*(nrAtoms)+ii,0])
        #        colors_y.InsertValue(ii,ASDReading.full_mom[time*(nrAtoms)+ii,1])
        #        #colors_z.InsertValue(ii,ASDReading.full_mom[time*(nrAtoms)+ii,2])
        #    else:
        #        colors_x.InsertValue(ii,(ASDReading.full_mom[time*(nrAtoms)+ii,0]-min_x)/(max_x-min_x))
        #        colors_y.InsertValue(ii,(ASDReading.full_mom[time*(nrAtoms)+ii,1]-min_y)/(max_y-min_y))
        #        #colors_z.InsertValue(ii,(ASDReading.full_mom[time*(nrAtoms)+ii,2]-min_z)/(max_z-min_z))
        #print('     done.')
        #-----------------------------------------------------------------------
        # Pass the colors to an array
        #-----------------------------------------------------------------------
        colors=[]
        colors.append(colors_x)
        colors.append(colors_y)
        colors.append(colors_z)
        return vectors, colors,ASDReading.number_time_steps,ASDReading.time_sep

    #--------------------------------------------------------------------------------
    # Function to determine the type of file and the format of the file
    #--------------------------------------------------------------------------------
    def check_format(self,filename):

        try:
            file_data = open(filename)
        except:
            file_data=filename
        line= file_data.readline()
        data = str.split(line)
        if data[0][0]=='#':
            type_fmt='new'
            file_type=''
        else:
            type_fmt='old'
            if len(data)==1:
                file_type='restart'
            else:
                file_type='moment'
        file_data.seek(0)
        return type_fmt,file_type
    #--------------------------------------------------------------------------------
    # @brief Read the site dependent and time dependent energy
    # @detail Read the site dependent and time dependent energy. It works in the
    # same way that the visualization of the magnetic moments.
    # in the same way that for the moments one must calculate and then pass back the
    # number of time steps in the system.
    # @author Jonathan Chico
    #--------------------------------------------------------------------------------
    def readEnergyData(self,file_ene, time, nrAtoms,temp_count):
        from vtk import vtkFloatArray
        import numpy as np
        import pandas as pd
        #-----------------------------------------------------------------------
        # Create a Double array which represents the vectors
        #-----------------------------------------------------------------------
        ene_xc   = vtkFloatArray()
        ene_dm   = vtkFloatArray()
        ene_bq   = vtkFloatArray()
        ene_pd   = vtkFloatArray()
        ene_tot  = vtkFloatArray()
        ene_dip  = vtkFloatArray()
        ene_ani  = vtkFloatArray()
        ene_bqdm = vtkFloatArray()
        ene_bext = vtkFloatArray()
        ene_chir = vtkFloatArray()
        #-----------------------------------------------------------------------
        # Define number of elements
        #-----------------------------------------------------------------------
        ene_xc.SetNumberOfComponents(1)
        ene_dm.SetNumberOfComponents(1)
        ene_bq.SetNumberOfComponents(1)
        ene_pd.SetNumberOfComponents(1)
        ene_tot.SetNumberOfComponents(1)
        ene_dip.SetNumberOfComponents(1)
        ene_ani.SetNumberOfComponents(1)
        ene_bqdm.SetNumberOfComponents(1)
        ene_bext.SetNumberOfComponents(1)
        ene_chir.SetNumberOfComponents(1)
        #-----------------------------------------------------------------------
        # Find the number of time steps to perform the simulation
        #-----------------------------------------------------------------------
        if time==0:
            #-------------------------------------------------------------------
            # Read the data in the moment file to find the number of time steps
            #-------------------------------------------------------------------
            ASDReading.full_ene=pd.read_csv(file_ene,header=None,skiprows=1,        \
            delim_whitespace=True,usecols=[0,3,4,5,6,7,8,9,10,11,13]).values
            #-------------------------------------------------------------------
            # Defining the number of time steps
            #-------------------------------------------------------------------
            ASDReading.number_time_steps=len(ASDReading.full_ene)/nrAtoms
            temp_count=ASDReading.number_time_steps
            ASDReading.time_sep=np.unique(ASDReading.full_ene[:,0])
        else:
            ASDReading.number_time_steps=temp_count
        #-----------------------------------------------------------------------
        # Loop over all the atoms
        #-----------------------------------------------------------------------
        for ii in range(0,nrAtoms):
            #-------------------------------------------------------------------
            # Pass the data from the numpy arrays to vtk data structures
            #-------------------------------------------------------------------
            ene_tot.InsertValue(ii,ASDReading.full_ene[time*(nrAtoms)+ii,1])
            ene_xc.InsertValue(ii,ASDReading.full_ene[time*(nrAtoms)+ii,2])
            ene_ani.InsertValue(ii,ASDReading.full_ene[time*(nrAtoms)+ii,3])
            ene_dm.InsertValue(ii,ASDReading.full_ene[time*(nrAtoms)+ii,4])
            ene_pd.InsertValue(ii,ASDReading.full_ene[time*(nrAtoms)+ii,5])
            ene_bqdm.InsertValue(ii,ASDReading.full_ene[time*(nrAtoms)+ii,6])
            ene_bq.InsertValue(ii,ASDReading.full_ene[time*(nrAtoms)+ii,7])
            ene_dip.InsertValue(ii,ASDReading.full_ene[time*(nrAtoms)+ii,8])
            ene_bext.InsertValue(ii,ASDReading.full_ene[time*(nrAtoms)+ii,9])
            ene_chir.InsertValue(ii,ASDReading.full_ene[time*(nrAtoms)+ii,10])
        #-----------------------------------------------------------------------
        # Pass the energies to an array
        #-----------------------------------------------------------------------
        energies=[]
        energies.append(ene_tot)
        energies.append(ene_xc)
        energies.append(ene_dm)
        energies.append(ene_ani)
        energies.append(ene_bq)
        energies.append(ene_bqdm)
        energies.append(ene_pd)
        energies.append(ene_bext)
        energies.append(ene_dip)
        energies.append(ene_chir)
        return energies, ASDReading.number_time_steps,ASDReading.time_sep
    ############################################################################
    # @brief Read the data needed to visualize the time evolution of the KMC particles
    # @details Read the data needed to visualize the time evolution of the KMC particles,
    # it works in the same way that the visualization of the magnetic moments.
    # in the same way that for the moments one must calculate and then pass back the
    # number of time steps in the system.
    # @author Jonathan Chico
    ############################################################################
    def readKMCData(self,file_KMC,time,temp_nrKMCpar):
        from vtk import vtkPoints
        import numpy as np
        import pandas as pd
        
        coord_KMC = vtkPoints()
        #-----------------------------------------------------------------------
        # Read the file first to see how many KMC particles there are
        #-----------------------------------------------------------------------
        if time==0:
            #-------------------------------------------------------------------
            # Read the data in the KMC file to find the number of KMC particles
            #-------------------------------------------------------------------
            ASDReading.full_KMC=pd.read_csv(file_KMC,header=None,                   \
            delim_whitespace=True,usecols=[0,3,4,5]).values
            ASDReading.nrKMCpar=len(np.unique(ASDReading.full_KMC[:,0]))
            temp_nrKMCpar=ASDReading.nrKMCpar
        else:
            ASDReading.nrKMCpar=temp_nrKMCpar
        #-----------------------------------------------------------------------
        # Passing the helper arrays to vtk data structures
        #-----------------------------------------------------------------------
        for ii in range(0,ASDReading.nrKMCpar):
            coord_KMC.InsertPoint(ii,                                               \
            ASDReading.full_KMC[time*ASDReading.nrKMCpar+ii,1],                     \
            ASDReading.full_KMC[time*ASDReading.nrKMCpar+ii,2],                     \
            ASDReading.full_KMC[time*ASDReading.nrKMCpar+ii,3])
        return coord_KMC,ASDReading.nrKMCpar
    ############################################################################
    # Read and arrange the neighbours from the struct file
    # If one reads this file correctly one can generate a visual representation
    # of the neighbouring atoms as defined in the neighbour map
    # @author Jonathan Chico
    ############################################################################
    def readNeighbours(self,file):
        import numpy as np
        import pandas as pd
        # Read the data using pandas
        neigh_data=pd.read_csv(file,skiprows=5,header=None,delim_whitespace=True,   \
        usecols=[0,1,3,7]).values
        # Store the data in convenient arrays
        curr_atom=neigh_data[:,0]
        neighbours=neigh_data[:,1]
        neigh_types=neigh_data[:,2]
        Neigh_strength=neigh_data[:,3]
        num_types_total=len(np.unique(neigh_types))
        return neighbours,Neigh_strength,curr_atom,neigh_types,num_types_total
    ############################################################################
    # Read and arrange the neighbours from the dmdata file
    # If one reads this file correctly one can generate a visual representation
    # of the neighbouring atoms as defined in the dm neighbour map
    # @author Jonathan Chico
    ############################################################################
    def readDMNeighbours(self,file):
        import numpy as np
        import pandas as pd
        # Read the data using pandas
        neigh_data=pd.read_csv(file,skiprows=5,header=None,delim_whitespace=True,   \
        usecols=[0,1,3,7,8,9]).values
        # Store the data in convenient arrays
        curr_atom=neigh_data[:,0]
        neighbours=neigh_data[:,1]
        neigh_types=neigh_data[:,2]
        DM_vec=neigh_data[:,3:6]
        DM_strength=np.zeros(len(neigh_data),dtype=np.float64)
        for ii in range(0,len(neigh_data)):
            DM_strength[ii]=np.sqrt(DM_vec[ii].dot(DM_vec[ii]))
        num_types_total=len(np.unique(neigh_types))
        return neighbours,DM_vec,DM_strength,curr_atom,neigh_types,num_types_total
    ############################################################################
    # The previously defined set of data from the neighbours can then be
    # stored in vtk friendly arrays
    # Notice the variable iAtom, which defines which atom is being currently
    # visualized, that is the variable to be set by the slider
    # @author Anders Bergman
    ############################################################################
    def setNeighbours(self,neighbours,iAtom,coords,Neigh_strength,curr_atom,neigh_types):
        from vtk import vtkPoints,vtkFloatArray
        import numpy as np
        #-----------------------------------------------------------------------
        # Set the arrays for the neighbours and the center atom
        #-----------------------------------------------------------------------
        neighpoints = vtkPoints()
        atompoint   = vtkPoints()
        # Find the indices which correspond to the current atom
        ind=np.where(curr_atom==(iAtom+1))
        # Find the number of neighbours
        nNeighs=len(ind[0])
        # Find the number of types
        (types,types_counters)=np.unique(neigh_types[ind[0]],return_counts=True)
        num_types=len(types)
        # Find the positions of the current atom
        (x,y,z)=coords.GetPoint(iAtom)
        atompoint.InsertPoint(0,x,y,z)
        ntypes = vtkFloatArray()
        ntypes.SetNumberOfComponents(1)
        ntypes.InsertValue(0,1.25)
        colors = vtkFloatArray()
        colors.SetNumberOfComponents(1)
        colors.InsertValue(0,0)
        #-----------------------------------------------------------------------
        # Pass the data to vtk arrays
        #-----------------------------------------------------------------------
        for iNeigh in range(0,nNeighs):
            cNeigh=int(neighbours[ind[0][iNeigh]]-1)
            (xx,yy,zz)=coords.GetPoint(cNeigh)
            neighpoints.InsertPoint(iNeigh, xx,yy,zz)
            ntypes.InsertValue(iNeigh,neigh_types[ind[0][iNeigh]])
            colors.InsertValue(iNeigh,Neigh_strength[ind[0][iNeigh]])
        return neighpoints,atompoint, ntypes, colors,nNeighs,num_types,types_counters,types
    ############################################################################
    # @brief The previously defined set of data from the DM vectors neighbours can then be
    # stored in vtk friendly arrays
    # @detailsThe previously defined set of data from the DM vectors neighbours can then be
    # stored in vtk friendly arrays.
    #  Notice the variable iAtom, which defines which atom is being currently
    # visualized, that is the variable to be set by the slider
    # @author Jonathan Chico
    ############################################################################
    def setDMNeighbours(self,neighbours,iAtom,coords,DM_vec,DM_strength,curr_atom,neigh_types):
        from vtk import vtkPoints,vtkFloatArray
        import numpy as np
        #-----------------------------------------------------------------------
        # Set the arrays for the neighbours and the center atom
        #-----------------------------------------------------------------------
        neighpoints = vtkPoints()
        atompoint   = vtkPoints()
        DM_vectors = vtkFloatArray()
        #-----------------------------------------------------------------------
        # Define number of elements
        #-----------------------------------------------------------------------
        DM_vectors.SetNumberOfComponents(3)
        # Find the indices which correspond to the current atom
        ind=np.where(curr_atom==(iAtom+1))
        # Find the number of neighbours
        tol=1e-10
        ind_non_zero=np.where(DM_strength[ind[0]]>tol)
        nNeighs=len(ind[0][ind_non_zero[0]])
        # Find the number of types
        (types,types_counters)=np.unique(neigh_types[ind[0]],return_counts=True)
        num_types=len(types)
        # Find the positions of the current atom
        (x,y,z)=coords.GetPoint(iAtom)
        atompoint.InsertPoint(0,x,y,z)
        ntypes = vtkFloatArray()
        ntypes.SetNumberOfComponents(1)
        ntypes.InsertValue(0,1.25)
        colors = vtkFloatArray()
        colors.SetNumberOfComponents(1)
        colors.InsertValue(0,0)
        #-----------------------------------------------------------------------
        # Pass the data to vtk arrays
        #-----------------------------------------------------------------------
        for iNeigh in range(0,nNeighs):
            cNeigh=int(neighbours[ind[0][ind_non_zero[0]][iNeigh]]-1)
            (xx,yy,zz)=coords.GetPoint(cNeigh)
            neighpoints.InsertPoint(iNeigh, xx,yy,zz)
            ntypes.InsertValue(iNeigh,neigh_types[ind[0][ind_non_zero[0]][iNeigh]])
            colors.InsertValue(iNeigh,DM_strength[ind[0][ind_non_zero[0]][iNeigh]])
            DM_vectors.InsertTuple3(iNeigh,\
            DM_vec[ind[0][ind_non_zero[0]][iNeigh],0]/DM_strength[ind[0][ind_non_zero[0]][iNeigh]],\
            DM_vec[ind[0][ind_non_zero[0]][iNeigh],1]/DM_strength[ind[0][ind_non_zero[0]][iNeigh]],\
            DM_vec[ind[0][ind_non_zero[0]][iNeigh],2]/DM_strength[ind[0][ind_non_zero[0]][iNeigh]])
        return neighpoints,atompoint,ntypes,colors,DM_vectors,nNeighs,num_types,types_counters,types
