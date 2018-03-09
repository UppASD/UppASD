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
from PyQt4 import QtCore, QtGui

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
    can be done, as for very thin samples, but with thickness different from 0, a 2D density rendering
    approach is more appropiate than a full 3D render approach
    When including any further visualization capabilities, requiring different files
    it is recommended to include the reading routines inside this class.
    """
    restart=[]
    moments=[]
    posfiles=[]
    kmcfiles=[]
    structfiles=[]
    def ReadingWrapper(self,mode,viz_type):
        # This selects which type of visualization one has
        # moments or neigbours
        if viz_type=='M':

            if mode==1:
                # restartfile
                if len(ASDReading.restart)>0:
                    ASDReading.MomentsFile = open(ASDReading.restart)
                else:
                    print "No file name selected from menu. Trying to find a 'restart.*.out' file"
                    ASDReading.restart=glob.glob("restart.*.out")
                    if len(ASDReading.restart)>0:
                        ASDReading.MomentsFile = open(ASDReading.restart[0])
                    else:
                        print "No 'restart.*.out' file found!"
                        print "I'm sorry, Dave. I'm afraid I can't do that."
                        sys.exit()

            elif mode==2:
                if len(ASDReading.moments)>0:
                    # MomentsFile
                    ASDReading.MomentsFile=ASDReading.moments
                else:
                    print "No file name selected from menu. Trying to find a 'moment.*.out' file"
                    ASDReading.moments=glob.glob("moment.*.out")
                    if len(ASDReading.moments)>0:
                        ASDReading.MomentsFile=ASDReading.moments[0]
                    else:
                        print "No 'moment.*.out' file found! "
                        print "I'm sorry, Dave. I'm afraid I can't do that."
                        sys.exit()

        # Neighbour mode
        elif viz_type=='N':
            # Structure file (for neighbour visualization)
            if len(ASDReading.structfiles)>0:
                ASDReading.structFile=open(ASDReading.structfiles)
                self.neigh_flag=True
            else:
                print "No file name selected from menu. Trying to find a 'struct.*.out' file"
                ASDReading.structfiles=glob.glob("struct.*.out")
                if len(ASDReading.structfiles)>0:
                    self.neigh_flag=True
                    ASDReading.structFile = open(ASDReading.structfiles[0])
                else:
                    self.neigh_flag=False
                    print "No 'struct.*.out' file found!"
                    print "I'm sorry, Dave. I'm afraid I can't do that."
                    sys.exit()

        if viz_type=='N' or viz_type=='M':
        	# Coordinates file
        	if len(ASDReading.posfiles)>0:
        	    atomsFile = open(self.posfiles)
        	else:
        	    print "No file name selected from menu. Trying to find a 'coord.*.out' file"
        	    ASDReading.posfiles=glob.glob("coord.*.out")
        	    if len(ASDReading.posfiles)>0:
        	        atomsFile = open(ASDReading.posfiles[0])
        	    else:
        	        print "No 'coord.*.out' file found!"
        	        print "I'm sorry, Dave. I'm afraid I can't do that."
        	        sys.exit()

        	# Cluster coordinates file
        	posfiles_c=glob.glob("posfile_clus.dat")
        	# Define a flag to check whether or not the file exists
        	if len(posfiles_c)>0:
        	    self.cluster_flag=True
        	    # Open the file if it exists
        	    atomsFile_c = open(posfiles_c[0])
        	else:
        	    self.cluster_flag=False

        	# If the file has been defined open it
        	if len(ASDReading.kmcfiles)>0:
        	    # KMC  file information
        	    ASDReading.kmc_flag=True
        	    ASDReading.KMCFile=ASDReading.kmcfiles
        	else:
        	    ASDReading.kmcfiles=glob.glob("kmc_info.*.out")
        	    if len(ASDReading.kmcfiles)>0:
        	        ASDReading.kmc_flag=True
        	        # Opet the kmc file if it exists
        	        ASDReading.KMCFile = ASDReading.kmcfiles[0]
        	    else:
        	        ASDReading.kmc_flag=False

        	# Actually reading the files
        	# Setting the coordinates of the full system
        	(ASDReading.coord,ASDReading.nrAtoms,ASDReading.flag_2D)=self.readAtoms(atomsFile)
        	# Checking if the clusters are present
        	if self.cluster_flag:
        	    # Setting the coordinates of the impurity cluster
        	    (ASDReading.coord_c,ASDReading.nrAtoms_c,ASDReading.colors_clus,ASDReading.points_clus_imp,\
        	    ASDReading.colors_imp,ASDReading.imp_nr)=\
        	    self.readAtoms_clus(atomsFile_c)

        # If the type of visualization is on the moments mode read the data about it
        if viz_type=='M':
            # Read the data for the vectors
            (ASDReading.rest_mom,ASDReading.colors_x,ASDReading.colors_y,ASDReading.colors_z,ASDReading.number_time_steps)=\
            self.readVectorsData(ASDReading.MomentsFile,0,ASDReading.nrAtoms,mode,1)

            # Check if there are KMC files present
            if ASDReading.kmc_flag:
                (ASDReading.coord_KMC,self.nrKMCpar)=self.readKMCData(ASDReading.KMCFile,0,1)
            # Check if the system should be represented as a 2D structure
            if ASDReading.flag_2D:
                # Creating a 2D data set
                (ASDReading.selected_points,ASDReading.selected_vectors,ASDReading.selected_colors_x,\
                ASDReading.selected_colors_y,ASDReading.selected_colors_z,ASDReading.reduced_nrAtoms)=\
                self.create2Dsystem(ASDReading.rest_mom,ASDReading.colors_x,ASDReading.colors_y,ASDReading.colors_z,\
                ASDReading.coord,ASDReading.nrAtoms)
            else:
                ASDReading.selected_points=ASDReading.coord
                ASDReading.selected_vectors=ASDReading.rest_mom
                ASDReading.selected_colors_x=ASDReading.colors_x
                ASDReading.selected_colors_y=ASDReading.colors_y
                ASDReading.selected_colors_z=ASDReading.colors_z
                ASDReading.reduced_nrAtoms=ASDReading.nrAtoms

        # If the type of visualization is about the neigbours read that data instead
        elif viz_type=='N':
            ASDReading.neighbours = self.readNeighbours(ASDReading.structFile)
            # Calculate neighbours to iAtom
            (ASDReading.neighs, ASDReading.nTypes) = self.setNeighbours(ASDReading.neighbours,0,ASDReading.coord)


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

        coord_x=[]
        coord_y=[]
        coord_z=[]
        # Read the data with numpy
        data=np.loadtxt(file_coord)
        coord_x=data[:,1]
        coord_y=data[:,2]
        coord_z=data[:,3]
        # Define the number of atoms in the system
        ASDReading.nrAtoms=len(coord_x)

        # Pass the numpy type arrays to vtk objects
        for ii in range(0,ASDReading.nrAtoms):
            points.InsertPoint(ii,coord_x[ii],coord_y[ii],coord_z[ii])

        # Data to check if one should consider the data to be rendered in 2D or 3D
        max_x=max(coord_x)
        min_x=min(coord_x)
        max_y=max(coord_y)
        min_y=min(coord_y)
        max_z=max(coord_z)
        min_z=min(coord_z)
        dist_x=np.sqrt((max_x-min_x)**2)
        dist_y=np.sqrt((max_y-min_y)**2)
        dist_z=np.sqrt((max_z-min_z)**2)
        min_dist=min(dist_x,dist_y,dist_z)

        # If the minimum distace is small set the flag to be 2D, this does not mean the
        # system is trully 2D, if not that the distance is small enough, such that the
        # 2D delaunay algorithm is used
        if min_dist<tol:
            ASDReading.flag_2D=True

        # Cleanup temporary arrays
        del coord_x
        del coord_y
        del coord_z
        return points, ASDReading.nrAtoms,ASDReading.flag_2D

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

        coord_x=[]
        coord_y=[]
        coord_z=[]
        # Read the data with numpy
        data=np.loadtxt(file_clus)
        coord_x=data[:,2]
        coord_y=data[:,3]
        coord_z=data[:,4]
        # Define the number of atoms in the cluster
        ASDReading.nrAtoms_clus=len(coord_x)

        # This will ensure that one can keep track of how many impurities one
        # actually has in the selected area
        imp_nr=0
        # Pass the numpy type data to vtk data
        for ii in range(0,ASDReading.nrAtoms_clus):
            points_clus.InsertPoint(ii,coord_x[ii],coord_y[ii],coord_z[ii])
            colors_clus.InsertTuple3(ii,0,0,0)
            # Data to display the center of the impurity cluster
            dist=np.sqrt((coord_x[ii]-coord_x[1])**2+\
                         (coord_y[ii]-coord_y[1])**2+\
                         (coord_z[ii]-coord_z[1])**2)
            # This is to ensure that the impurity cluster and the embeded cluster
            # are different structures that can be used differently
            if dist<1 and ii==1:
                colors_imp.InsertTuple3(imp_nr,51,160,44)
                points_clus_imp.InsertPoint(imp_nr,coord_x[ii],coord_y[ii],coord_z[ii])
                imp_nr=imp_nr+1
            elif dist<1:
                colors_imp.InsertTuple3(imp_nr,106,61,154)
                points_clus_imp.InsertPoint(imp_nr,coord_x[ii],coord_y[ii],coord_z[ii])
                imp_nr=imp_nr+1

        # Cleanup the leftover data
        del coord_x
        del coord_y
        del coord_z
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

        # Create numpy helper arrays
        mom_x=np.zeros(nrAtoms,dtype=np.float32)
        mom_y=np.zeros(nrAtoms,dtype=np.float32)
        mom_z=np.zeros(nrAtoms,dtype=np.float32)

        # Check which kind of visualization is being done restartfile or momentfile
        if mode==1:
            # Read the data in the restartfile
            mom_x=[]
            mom_y=[]
            mom_z=[]

            data=np.loadtxt(file_mom,skiprows=1)
            # Pass the data to the numpy helper arrays
            mom_x=data[:,3]
            mom_y=data[:,4]
            mom_z=data[:,5]
            ASDReading.number_time_steps=1

        elif mode==2:
            if time==0:
                open_file=open(file_mom)
                # Read the data in the moment file to find the number of time steps
                data=np.loadtxt(open_file,usecols=(0))
                # Defining the number of time steps
                ASDReading.number_time_steps=len(data)/nrAtoms
                temp_count=ASDReading.number_time_steps
                open_file.close()
            else:
                ASDReading.number_time_steps=temp_count

            ii=0
            # Read the first nrAtoms lines
            while ii<nrAtoms:
                linenum=int(ii+time*nrAtoms+1)
                line = linecache.getline(str(file_mom),linenum)
                data = string.split(line)
                mom_x[ii]=float(data[2])
                mom_y[ii]=float(data[3])
                mom_z[ii]=float(data[4])
                ii=ii+1

        for ii in range(0,nrAtoms):
            # Pass the data from the numpy arrays to vtk data structures
            vectors.InsertTuple3(ii,mom_x[ii],mom_y[ii],mom_z[ii])
            m_x=(mom_x[ii])
            m_y=(mom_y[ii])
            m_z=(mom_z[ii])
            colors_x.InsertValue(ii,m_x)
            colors_y.InsertValue(ii,m_y)
            colors_z.InsertValue(ii,m_z)

        # Cleanup the leftover data
        del mom_z
        del mom_y
        del mom_x
        return vectors, colors_x,colors_y,colors_z,ASDReading.number_time_steps

    ############################################################################
    # Read the data needed to visualize the time evolution of the KMC particles
    ############################################################################
    def readKMCData(self,file_KMC,time,temp_nrKMCpar):

        coord_KMC = vtk.vtkPoints()

        # Read the file first to see how many KMC particles there are
        if time==0:
            open_file=open(file_KMC)
            # Read the data in the KMC file to find the number of KMC particles
            data=np.loadtxt(open_file,usecols=(0))
            open_file.close()
            nrKMCpar=len(np.unique(data))
            temp_nrKMCpar=nrKMCpar
        else:
            nrKMCpar=temp_nrKMCpar

        # Create helper numpy arrays
        coord_x=np.zeros(nrKMCpar,dtype=np.float32)
        coord_y=np.zeros(nrKMCpar,dtype=np.float32)
        coord_z=np.zeros(nrKMCpar,dtype=np.float32)

        ii=0
        # Read the first nrKMCpar lines
        while ii<nrKMCpar:
            linenum=int(ii+time*nrKMCpar+1)
            line = linecache.getline(str(file_KMC),linenum)
            data = string.split(line)
            coord_x[ii]=float(data[3])
            coord_y[ii]=float(data[4])
            coord_z[ii]=float(data[5])
            ii=ii+1

        # Passing the helper arrays to vtk data structures
        for ii in range(0,nrKMCpar):
            coord_KMC.InsertPoint(ii,coord_x[ii],coord_y[ii],coord_z[ii])

        # Cleanup leftover data
        del coord_x
        del coord_y
        del coord_z
        return coord_KMC,nrKMCpar

    ############################################################################
    # Read and arrange the neighbours from the struct file
    # If one reads this file correctly one can generate a visual representation
    # of the neighbouring atoms as defined in the neighbour map
    # @author Anders Bergman
    ############################################################################
    def readNeighbours(self,file):
        # Local variables to check the dimensions of the arrays that will
        # contain the data
        maxNeigh = 0
        nAtom = 0

        # Read struct file to find max no. neighbours
        # Must read the first three lines to remove the unnecesary lines
        line = file.readline()
        line = file.readline()
        line = file.readline()
        # Read all data
        while line:
            data = string.split(line)
            # Find how many atoms there are in the system
            if  data[0]=="Atom=" :
                nAtom=nAtom+1
                iAtom=int(data[1])
            # Find how many neighbours a given atom has
            if  data[0]=="Shell=" :
                nNeigh=int(data[5])
                maxNeigh=max(nNeigh,maxNeigh)
            line = file.readline()

        # Allocate neighbour array
        neighbours=[]

        # Rewind file
        file.seek(0)

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
    def create2Dsystem(self,vectors,colors_x,colors_y,colors_z,coord,nrAtoms):
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
            if z==0:
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
        dlg = QtGui.QFileDialog()
        dlg.setFileMode(QtGui.QFileDialog.AnyFile)
        dlg.setDirectory('.')
        if dlg.exec_():
            ASDReading.posfiles=dlg.selectedFiles()[0]

    ############################################################################
    # Finding the file name of the restart file
    ############################################################################
    def getRestartFile(self):
        dlg = QtGui.QFileDialog()
        dlg.setFileMode(QtGui.QFileDialog.AnyFile)
        dlg.setDirectory('.')
        if dlg.exec_():
            ASDReading.restart=dlg.selectedFiles()[0]
    ############################################################################
    # Finding the file name of the moment file
    ############################################################################
    def getMomentFile(self):
        dlg = QtGui.QFileDialog()
        dlg.setFileMode(QtGui.QFileDialog.AnyFile)
        dlg.setDirectory('.')
        if dlg.exec_():
            ASDReading.moments=dlg.selectedFiles()[0]

    ############################################################################
    # finding the file name of the structuct file
    ############################################################################
    def getStructFile(self):
        dlg = QtGui.QFileDialog()
        dlg.setFileMode(QtGui.QFileDialog.AnyFile)
        dlg.setDirectory('.')
        if dlg.exec_():
            ASDReading.structfiles=dlg.selectedFiles()[0]

    ############################################################################
    # finding the file name of the KMC file
    ############################################################################
    def getKMCFile(self):
        dlg = qtgui.QFileDialog()
        dlg.setFileMode(QtGui.QFileDialog.AnyFile)
        dlg.setDirectory('.')
        if dlg.exec_():
            ASDReading.kmcfiles=dlg.selectedFiles()[0]
