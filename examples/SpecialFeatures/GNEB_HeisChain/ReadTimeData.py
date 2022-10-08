#Read atom- and direction data

import string
import vtk

# Read Location of Atoms
def readAtoms(file,Nmax):

    points = vtk.vtkPoints()
    nrAtoms=0
    # Read ahead
    line = file.readline()
    data = string.split(line)
    
    # Read all data
    while line:
	if nrAtoms <= Nmax:
	    a = int(data[0])
	    x, y, z = float(data[1]), float(data[2]), float(data[3])
	    #print "a ", a, " x ", x, " y ", y, " z ", z  
	    #points.InsertPoint(a, x, y, z)
	    points.InsertPoint(nrAtoms, x, y, z)
	nrAtoms=nrAtoms+1
	line = file.readline()
	data = string.split(line)
    return points, nrAtoms

# Read vectors
# We must know the time step and the number of atoms per time
def readVectors(file, time, nrAtoms,Nmax):
    # Create a Double array which represents the vectors
    vectors = vtk.vtkDoubleArray()
    # Define number of elemnts
    vectors.SetNumberOfComponents(3)

    i=0
    # Read all data for a certain time
    while i < nrAtoms:
        line = file.readline()
	if i < Nmax:
	    data = string.split(line)
	    t, a = int(data[0]), int(data[1])
	    x, y, z = float(data[2]), float(data[3]), float(data[4])
	    m = float(data[5])
	    #print "t ",t, " a ", a, " x ", x, " y ", y, " z ", z, " m ", m  
	    #vectors.InsertTuple3(a ,x, y, z)
	    vectors.InsertTuple3(i ,x, y, z)
        i=i+1
    return vectors
# Read vectors
# We must know the time step and the number of atoms per time
def readVectorsData(file, time, nrAtoms,Nmax):
    # Create a Double array which represents the vectors
    vectors = vtk.vtkDoubleArray()
#    mdata = vtk.vtkDoubleArray()
    # Define number of elemnts
    vectors.SetNumberOfComponents(3)
#    mdata.SetNumberOfComponents(1)
#    mdata.SetNumberOfValues(nrAtoms)

    i=0
    # Read all data for a certain time
    while i < nrAtoms:
        line = file.readline()
	if i < Nmax:
	    #print i
	    data = string.split(line)
	    t, a = int(data[0]), int(data[1])
	    x, y, z = float(data[2]), float(data[3]), float(data[4])
	    m = float(data[5])
#           m=z*0.98
#            m=(z+1.01)/2.05
#           m=(x+1.00)/3
#            if m <= 0.0 : 
#               m=float(0.001)
#            if m >= 1.0 : 
#               m=float(0.999)
#            m=(x+y+z)/1.8
            x=x*m*0.5#*0.95
            y=y*m*0.5#*0.95
            z=z*m*0.5#*0.95
#           x=x*m
#           y=y*m
#           z=z*m
#            if t == 2700:
#              print "t ",t, " a ", a, " x ", x, " y ", y, " z ", z, " m ", m  , (x*x+y*y+z*z)
#           print "t ",t, " a ", a, " x ", x, " y ", y, " z ", z, " m ", m  
	    #vectors.InsertTuple3(a ,x, y, z)
	    vectors.InsertTuple3(i ,x, y, z)
#            vectors.SetValue(i,m)
        i=i+1
    return vectors


# Read vectors as points on the sphere
# We must know the time step and the number of atoms per time
def readPoints(file, time, nrAtoms):

    points = vtk.vtkPoints()

    i=0
    # Read all data for a certain time
    while i < nrAtoms:
        line = file.readline()
        data = string.split(line)
        t, a = int(data[0]), int(data[1])
        x, y, z = float(data[2]), float(data[3]), float(data[4])
        m = float(data[5])
        points.InsertPoint(a, x, y, z)


        i=i+1
    return points
