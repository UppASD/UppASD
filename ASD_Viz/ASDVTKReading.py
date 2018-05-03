#!/usr/bin/env vtkpython
################################################################################
# CLASS: ASDVTKReading
# @author Jonathan Chico (08/09/2017)
# @description
# Wrapper for the reading of UppASD data for the ASD_visualizer, the present
# class handles the reading of the data from the UppASD generated files.
# This class also deals with data pre-prossesing to tackle both 2D and 3D
# structures for density map visualizations.
################################################################################
import vtk
import sys
import time
import glob
import string
import os.path
import linecache
import numpy as np
import pandas as pd
from PyQt5 import QtWidgets

################################################################################
# Class for reading the data from ASD files this makes use of the linecache
# routines as they should allow for higher speed than the equivalent numpy
#  and/or string options.
################################################################################
class ASDReading():
    """The ASDVTKReading clas encompasses all the reading routines necessary for the
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
    """
    restart=[]
    moments=[]
    posfiles=[]
    kmcfiles=[]
    full_mom=[]
    full_KMC=[]
    structfiles=[]
    not_read_pos=True
    not_read_mom=True
    not_read_neigh=True
    def ReadingWrapper(self,mode,viz_type):
        # This selects which type of visualization one has
        # moments or neigbours
        if viz_type=='M':

            if mode==1:
                # restartfile
                if len(ASDReading.restart)>0:
                    ASDReading.MomentsFile = open(ASDReading.restart)
                else:
                    print("No file name selected from menu. Trying to find a 'restart.*.out' file")
                    ASDReading.restart=glob.glob("restart.*.out")
                    if len(ASDReading.restart)>0:
                        ASDReading.restart=ASDReading.restart[0]
                        ASDReading.MomentsFile = open(ASDReading.restart)
                    else:
                        print("No 'restart.*.out' file found!")
                        print("I'm sorry, Dave. I'm afraid I can't do that.")
                        sys.exit()

            elif mode==2:
                if len(ASDReading.moments)>0:
                    # MomentsFile
                    ASDReading.MomentsFile=ASDReading.moments
                else:
                    print("No file name selected from menu. Trying to find a 'moment.*.out' file")
                    ASDReading.moments=glob.glob("moment.*.out")
                    if len(ASDReading.moments)>0:
                        ASDReading.moments=ASDReading.moments[0]
                        ASDReading.MomentsFile=ASDReading.moments
                    else:
                        print("No 'moment.*.out' file found! ")
                        print("I'm sorry, Dave. I'm afraid I can't do that.")
                        sys.exit()

        # Neighbour mode
        elif viz_type=='N':
            # Structure file (for neighbour visualization)
            if len(ASDReading.structfiles)>0:
                ASDReading.structFile=open(ASDReading.structfiles)
                self.neigh_flag=True
            else:
                print("No file name selected from menu. Trying to find a 'struct.*.out' file")
                ASDReading.structfiles=glob.glob("struct.*.out")
                if len(ASDReading.structfiles)>0:
                    ASDReading.structfiles=ASDReading.structfiles[0]
                    self.neigh_flag=True
                    ASDReading.structFile = open(ASDReading.structfiles)
                else:
                    self.neigh_flag=False
                    print("No 'struct.*.out' file found!")
                    print("I'm sorry, Dave. I'm afraid I can't do that.")
                    sys.exit()

        if viz_type=='N' or viz_type=='M' and ASDReading.not_read_pos:
            # Coordinates file
            if len(ASDReading.posfiles)>0:
                atomsFile = open(ASDReading.posfiles)
            else:
                print("No file name selected from menu. Trying to find a 'coord.*.out' file")
                ASDReading.posfiles=glob.glob("coord.*.out")
                if len(ASDReading.posfiles)>0:
                    ASDReading.posfiles=ASDReading.posfiles[0]
                    atomsFile = open(ASDReading.posfiles)
                else:
                    print("No 'coord.*.out' file found!")
                    print("I'm sorry, Dave. I'm afraid I can't do that.")
                    sys.exit()

            # Cluster coordinates file
            posfiles_c=glob.glob("posfile_clus.dat")
            # Define a flag to check whether or not the file exists
            if len(posfiles_c)>0:
                ASDReading.cluster_flag=True
                # Open the file if it exists
                atomsFile_c = open(posfiles_c[0])
            else:
                ASDReading.cluster_flag=False

            # If the file has been defined open it
            if len(ASDReading.kmcfiles)>0:
                # KMC  file information
                ASDReading.kmc_flag=True
                ASDReading.KMCFile=ASDReading.kmcfiles
            else:
                ASDReading.kmcfiles=glob.glob("kmc_info.*.out")
                if len(ASDReading.kmcfiles)>0:
                    ASDReading.kmc_flag=True
                    # Open the kmc file if it exists
                    ASDReading.KMCFile = ASDReading.kmcfiles[0]
                else:
                    ASDReading.kmc_flag=False
            ASDReading.not_read_pos=False
            # Actually reading the files
            # Setting the coordinates of the full system
            (ASDReading.coord,ASDReading.nrAtoms,ASDReading.flag_2D,ASDReading.min_val)=\
            self.readAtoms(atomsFile)
            # Checking if the clusters are present
            if ASDReading.cluster_flag:
                # Setting the coordinates of the impurity cluster
                (ASDReading.coord_c,ASDReading.nrAtoms_c,ASDReading.colors_clus,\
                ASDReading.points_clus_imp,ASDReading.colors_imp,ASDReading.imp_nr)=\
                self.readAtoms_clus(atomsFile_c)

        # If the type of visualization is on the moments mode read the data about it
        if viz_type=='M' and ASDReading.not_read_mom:
            # Read the data for the vectors
            (ASDReading.rest_mom,ASDReading.colors_x,ASDReading.colors_y,ASDReading.colors_z,\
            ASDReading.number_time_steps)=\
            self.readVectorsData(ASDReading.MomentsFile,0,ASDReading.nrAtoms,mode,1)
            ASDReading.not_read_mom=False
            # Check if there are KMC files present
            if ASDReading.kmc_flag:
                (ASDReading.coord_KMC,self.nrKMCpar)=self.readKMCData(ASDReading.KMCFile,0,1)

            ASDReading.selected_points=ASDReading.coord
            ASDReading.selected_vectors=ASDReading.rest_mom
            ASDReading.selected_colors_x=ASDReading.colors_x
            ASDReading.selected_colors_y=ASDReading.colors_y
            ASDReading.selected_colors_z=ASDReading.colors_z
            ASDReading.reduced_nrAtoms=ASDReading.nrAtoms

        # If the type of visualization is about the neigbours read that data instead
        elif viz_type=='N' and ASDReading.not_read_neigh:
            ASDReading.neighbours = self.readNeighbours(ASDReading.structFile)
            # Calculate neighbours to iAtom
            (ASDReading.neighs, ASDReading.nTypes) = self.setNeighbours(ASDReading.neighbours,0,ASDReading.coord)
            ASDReading.not_read_neigh=False
        ASDReading.mode=mode
        ASDReading.viz_type=viz_type
        return

    ############################################################################
    # Read Location of Atoms
    ############################################################################
    def readAtoms(self,file_coord):
        # Define the parameters used to check if the system is 2D
        tol=0.75
        ASDReading.flag_2D=False
        # Define vtkPoints as this type of arrays will be used to define the grid
        points = vtk.vtkPoints()
        # Read the data with pandas
        coord=pd.read_csv(file_coord,header=None,delim_whitespace=True,usecols=[1,2,3]).as_matrix()
        #
        # Define the number of atoms in the system
        ASDReading.nrAtoms=len(coord)
        #
        # Pass the numpy type arrays to vtk objects
        for ii in range(0,ASDReading.nrAtoms):
            points.InsertPoint(ii,coord[ii,0],coord[ii,1],coord[ii,2])
        #
        # Data to check if one should consider the data to be rendered in 2D or 3D
        max_x=max(coord[:,0])
        min_x=min(coord[:,0])
        max_y=max(coord[:,1])
        min_y=min(coord[:,1])
        max_z=max(coord[:,2])
        min_z=min(coord[:,2])
        dist_x=np.sqrt((max_x-min_x)**2)
        dist_y=np.sqrt((max_y-min_y)**2)
        dist_z=np.sqrt((max_z-min_z)**2)
        min_dist=min(dist_x,dist_y,dist_z)
        #
        # If the minimum distace is small set the flag to be 2D, this does not mean the
        # system is trully 2D, if not that the distance is small enough, such that the
        # 2D delaunay algorithm is used
        if min_dist<tol:
            ASDReading.flag_2D=True
        ASDReading.min_val=min_z
        #
        # Cleanup temporary arrays
        del coord
        return points, ASDReading.nrAtoms,ASDReading.flag_2D,ASDReading.min_val

    ############################################################################
    # Read Location of the cluster Atoms notice that this uses a posfile
    ############################################################################
    def readAtoms_clus(self,file_clus):

        # Define vtkPoints arrays for the creation of the grid
        points_clus = vtk.vtkPoints()
        points_clus_imp = vtk.vtkPoints()

        # Define UnsignedCharArrays for the definitions of the colors to be used
        colors_clus = vtk.vtkUnsignedCharArray()
        colors_imp = vtk.vtkUnsignedCharArray()
        # Set the size of the arrays
        colors_clus.SetNumberOfComponents(3)
        colors_imp.SetNumberOfComponents(3)

        # Read the data with pandas
        coord=pd.read_csv(file_clus,header=None,delim_whitespace=True,usecols=[2,3,4]).as_matrix()
        # Define the number of atoms in the cluster
        ASDReading.nrAtoms_clus=len(coord)

        # This will ensure that one can keep track of how many impurities one
        # actually has in the selected area
        imp_nr=0
        # Pass the numpy type data to vtk data
        for ii in range(0,ASDReading.nrAtoms_clus):
            points_clus.InsertPoint(ii,coord[ii,0],coord[ii,1],coord[ii,2])
            colors_clus.InsertTuple3(ii,0,0,0)
            # Data to display the center of the impurity cluster
            dist=np.sqrt((coord[ii,0]-coord[1,0])**2+\
                         (coord[ii,1]-coord[1,1])**2+\
                         (coord[ii,2]-coord[1,2])**2)
            # This is to ensure that the impurity cluster and the embeded cluster
            # are different structures that can be used differently
            if dist<1 and ii==1:
                colors_imp.InsertTuple3(imp_nr,51,160,44)
                points_clus_imp.InsertPoint(imp_nr,coord[ii,0],coord[ii,1],coord[ii,2])
                imp_nr=imp_nr+1
            elif dist<1:
                colors_imp.InsertTuple3(imp_nr,106,61,154)
                points_clus_imp.InsertPoint(imp_nr,coord[ii,0],coord[ii,1],coord[ii,2])
                imp_nr=imp_nr+1

        # Cleanup the leftover data
        del coord
        return points_clus, ASDReading.nrAtoms_clus,colors_clus,points_clus_imp,colors_imp,imp_nr

    ############################################################################
    # Read vectors
    # We must know the time step and the number of atoms per time
    ############################################################################
    def readVectorsData(self,file_mom, time, nrAtoms,mode,temp_count):
        # Create a Double array which represents the vectors
        vectors = vtk.vtkFloatArray()
        colors_x = vtk.vtkFloatArray()
        colors_y = vtk.vtkFloatArray()
        colors_z = vtk.vtkFloatArray()

        # Define number of elemnts
        vectors.SetNumberOfComponents(3)
        colors_x.SetNumberOfComponents(1)
        colors_y.SetNumberOfComponents(1)
        colors_z.SetNumberOfComponents(1)

        # Check which kind of visualization is being done restartfile or momentfile
        if mode==1:
            # Read the data in the restartfile
            ASDReading.full_mom=pd.read_csv(file_mom,header=None,skiprows=1,delim_whitespace=True,usecols=[3,4,5]).as_matrix()
            # Pass the data to the numpy helper arrays
            ASDReading.number_time_steps=1
        elif mode==2:
            if time==0:
                # Read the data in the moment file to find the number of time steps
                ASDReading.full_mom=pd.read_csv(file_mom,header=None,delim_whitespace=True,usecols=[2,3,4]).as_matrix()
                # Defining the number of time steps
                ASDReading.number_time_steps=len(ASDReading.full_mom)/nrAtoms
                temp_count=ASDReading.number_time_steps
            else:
                ASDReading.number_time_steps=temp_count

        for ii in range(0,nrAtoms):
            # Pass the data from the numpy arrays to vtk data structures
            vectors.InsertTuple3(ii,ASDReading.full_mom[time*(nrAtoms)+ii,0],\
                ASDReading.full_mom[time*(nrAtoms)+ii,1],\
                ASDReading.full_mom[time*(nrAtoms)+ii,2])
            colors_x.InsertValue(ii,ASDReading.full_mom[time*(nrAtoms)+ii,0])
            colors_y.InsertValue(ii,ASDReading.full_mom[time*(nrAtoms)+ii,1])
            colors_z.InsertValue(ii,ASDReading.full_mom[time*(nrAtoms)+ii,2])

        return vectors, colors_x,colors_y,colors_z,ASDReading.number_time_steps

    ############################################################################
    # Read the data needed to visualize the time evolution of the KMC particles
    ############################################################################
    def readKMCData(self,file_KMC,time,temp_nrKMCpar):

        coord_KMC = vtk.vtkPoints()

        # Read the file first to see how many KMC particles there are
        if time==0:
            # Read the data in the KMC file to find the number of KMC particles
            ASDReading.full_KMC=pd.read_csv(file_KMC,header=None,delim_whitespace=True,usecols=[0,3,4,5]).as_matrix()
            ASDReading.nrKMCpar=len(np.unique(ASDReading.full_KMC[:,0]))
            temp_nrKMCpar=ASDReading.nrKMCpar
        else:
            ASDReading.nrKMCpar=temp_nrKMCpar

        # Passing the helper arrays to vtk data structures
        for ii in range(0,ASDReading.nrKMCpar):
            coord_KMC.InsertPoint(ii,ASDReading.full_KMC[time*ASDReading.nrKMCpar+ii,1],\
            ASDReading.full_KMC[time*ASDReading.nrKMCpar+ii,2],\
            ASDReading.full_KMC[time*ASDReading.nrKMCpar+ii,3])

        return coord_KMC,ASDReading.nrKMCpar

    ############################################################################
    # Read and arrange the neighbours from the struct file
    # If one reads this file correctly one can generate a visual representation
    # of the neighbouring atoms as defined in the neighbour map
    # @author Anders Bergman
    ############################################################################
    def readNeighbours(self,file):
        # Allocate neighbour array
        neighbours=[]

        # Read struct file to find max no. neighbours
        # Must read the first three lines to remove the unnecesary lines
        line = file.readline()
        line = file.readline()
        line = file.readline()

        # Read all data
        locNeigh=[]
        while line:
            data = string.split(line)
            # Find atom number
            if  data[0]=="Atom=" :
                iAtom=int(data[1])
                locNeigh=[]
                tNeigh=0
                aNeigh=0
            # Read neighbours (five per line)
            if  data[0]=="Shell=" :
                nNeigh=int(data[5])
                tNeigh=tNeigh+nNeigh
                iNeigh=1
                if(int(data[1])==1):
                    locNeigh.append(nNeigh)
                for iline in range(0,nNeigh/5+1) :
                    line = file.readline()
                    data = string.split(line)
                    lCount=0
                    while (iNeigh <= nNeigh) and (lCount<5) :
                        if int(data[lCount])>0 :
                            locNeigh.append(int(data[lCount]))
                            aNeigh=aNeigh+1
                        iNeigh=iNeigh+1
                        lCount=lCount+1
            if(line[5]=="-"):
                locNeigh[0]=aNeigh
                neighbours.append(locNeigh)
            line = file.readline()

        locNeigh[0]=tNeigh
        neighbours.append(locNeigh)
        return neighbours

    ############################################################################
    # The previously defined set of data from the neighbours can then be
    # stored in vtk friendly arrays
    # Notice the variable iAtom, which defines which atom is being currently
    # visualized, that is the variable to be set by the slider
    # @author Anders Bergman
    ############################################################################
    def setNeighbours(self,neighbours,iAtom,coords):

        points = vtk.vtkPoints()
        nNeighs = neighbours[iAtom][0]
        (x,y,z)=coords.GetPoint(iAtom)
        points.InsertPoint(0, x,y,z)

        ntypes = vtk.vtkFloatArray()
        ntypes.SetNumberOfComponents(1)
        ntypes.InsertValue(0,1.25)

        for iNeigh in range(1,nNeighs+1) :
            cNeigh=neighbours[iAtom][iNeigh]-1
            (xx,yy,zz)=coords.GetPoint(cNeigh)
            points.InsertPoint(iNeigh, xx,yy,zz)
            ntypes.InsertValue(iNeigh,1)

        return points , ntypes

    ############################################################################
    # Temporary 2D parser, need to find a better way to see if data is 2D or 3D
    ############################################################################
    def create2Dsystem(self,vectors,colors_x,colors_y,colors_z,coord,nrAtoms,min_val):
        tol=1e-5
        # Create vtkPoints arrays to create the grid
        selected_points = vtk.vtkPoints()
        # Create the vtkFloatArray for the selected vectors
        selected_vectors = vtk.vtkFloatArray()
        # Create vtkFloatArray for the colors to be used with a Lookup Table
        selected_colors_x = vtk.vtkFloatArray()
        selected_colors_y = vtk.vtkFloatArray()
        selected_colors_z = vtk.vtkFloatArray()

        # Define number of elemnts
        selected_vectors.SetNumberOfComponents(3)
        selected_colors_x.SetNumberOfComponents(1)
        selected_colors_y.SetNumberOfComponents(1)
        selected_colors_z.SetNumberOfComponents(1)

        # Need to find how many atoms are in the 2d slice
        ASDReading.reduced_nrAtoms=0
        for ii in range(0,nrAtoms):
            # Get the coordinates
            (x,y,z)=coord.GetPoint(ii)
            # If one has a 2D system
            if abs(z-min_val)<tol:
                (xx,yy,zz)=coord.GetPoint(ii)
                selected_points.InsertPoint(ASDReading.reduced_nrAtoms,xx,yy,0)
                mx,my,mz=vectors.GetTuple3(ii)
                selected_vectors.InsertTuple3(ASDReading.reduced_nrAtoms,mx,my,mz)
                mm_x=colors_x.GetValue(ii)
                mm_y=colors_y.GetValue(ii)
                mm_z=colors_z.GetValue(ii)
                selected_colors_x.InsertValue(ASDReading.reduced_nrAtoms,mm_x)
                selected_colors_y.InsertValue(ASDReading.reduced_nrAtoms,mm_y)
                selected_colors_z.InsertValue(ASDReading.reduced_nrAtoms,mm_z)
                ASDReading.reduced_nrAtoms=ASDReading.reduced_nrAtoms+1

        return selected_points, selected_vectors, selected_colors_x,selected_colors_y,\
        selected_colors_z,ASDReading.reduced_nrAtoms


    ############################################################################
    # Finding the file name of the coordinate file
    ############################################################################
    def getCoordFile(self):
        dlg = QtWidgets.QFileDialog()
        dlg.setFileMode(QtWidgets.QFileDialog.AnyFile)
        dlg.setDirectory('.')
        if dlg.exec_():
            ASDReading.posfiles=dlg.selectedFiles()[0]

    ############################################################################
    # Finding the file name of the restart file
    ############################################################################
    def getRestartFile(self):
        dlg = QtWidgets.QFileDialog()
        dlg.setFileMode(QtWidgets.QFileDialog.AnyFile)
        dlg.setDirectory('.')
        if dlg.exec_():
            ASDReading.restart=dlg.selectedFiles()[0]
    ############################################################################
    # Finding the file name of the moment file
    ############################################################################
    def getMomentFile(self):
        dlg = QtWidgets.QFileDialog()
        dlg.setFileMode(QtWidgets.QFileDialog.AnyFile)
        dlg.setDirectory('.')
        if dlg.exec_():
            ASDReading.moments=dlg.selectedFiles()[0]

    ############################################################################
    # finding the file name of the structuct file
    ############################################################################
    def getStructFile(self):
        dlg = QtWidgets.QFileDialog()
        dlg.setFileMode(QtWidgets.QFileDialog.AnyFile)
        dlg.setDirectory('.')
        if dlg.exec_():
            ASDReading.structfiles=dlg.selectedFiles()[0]

    ############################################################################
    # finding the file name of the KMC file
    ############################################################################
    def getKMCFile(self):
        dlg = QtWidgets.QFileDialog()
        dlg.setFileMode(QtWidgets.QFileDialog.AnyFile)
        dlg.setDirectory('.')
        if dlg.exec_():
            ASDReading.kmcfiles=dlg.selectedFiles()[0]
