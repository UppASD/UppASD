#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import matplotlib.pyplot as plt
import spglib as spg
import seekpath as spth
import matplotlib.cm as cmap
import os.path
from scipy import ndimage


# In[2]:

hbar=6.582119514e-13


############################################################
# Set pyplot defaults
############################################################

font = {'family' : 'sans',
        'weight' : 'normal',
        'size'   : '14'}

plt.rc('font', **font) 
plt.rc('lines', lw=2)


# In[3]:


# In[4]:


#plt.xkcd()


# In[5]:


############################################################
# Read atom positions from UppASD position file
############################################################
def read_posfile(posfile):
   with open(posfile,'r') as pfile:
      lines=pfile.readlines()
      positions=np.empty([0,3])
      numbers=[]
      for idx,line in enumerate(lines):
         line_data=line.rstrip('\n').split()
         if len(line_data)>0:
            positions=np.vstack((positions,np.asarray(line_data[2:5])))
            numbers=np.append(numbers,np.asarray(line_data[1]))
      return positions,numbers


# In[6]:


############################################################
# Read important keywords from UppASD inputfile `inpsd.dat`
############################################################
def read_inpsd(ifile):
   with open(ifile,'r') as infile:
      lines=infile.readlines()
      for idx,line in enumerate(lines):
         line_data=line.rstrip('\n').split()
         if len(line_data)>0:
             # Find the simulation id
             if(line_data[0]=='simid'):
                 simid=line_data[1]
                 print('simid: ',simid)

             # Find the cell data
             if(line_data[0]=='cell'):
                cell=[]
                lattice=np.empty([0,3])
                line_data=lines[idx+0].split()
                cell=np.append(cell,np.asarray(line_data[1:4]))
                lattice=np.vstack((lattice,np.asarray(line_data[1:4])))
                line_data=lines[idx+1].split()
                cell=np.append(cell,np.asarray(line_data[0:3]))
                lattice=np.vstack((lattice,np.asarray(line_data[0:3])))
                line_data=lines[idx+2].split()
                cell=np.append(cell,np.asarray(line_data[0:3]))
                lattice=np.vstack((lattice,np.asarray(line_data[0:3])))
                print('cell: ',cell)
                print('lattice: ',lattice)

             # Find the size of the simulated cell
             if(line_data[0]=='ncell'):
                ncell_x=int(line_data[1])
                ncell_y=int(line_data[1])
                ncell_z=int(line_data[1])
                mesh=[ncell_x,ncell_y,ncell_z]

             if(line_data[0]=='timestep'):
                 timestep=line_data[1]
                 print('timestep: ',timestep)

             if(line_data[0]=='sc_nstep'):
                 sc_nstep=line_data[1]
                 print('sc_nstep: ',sc_nstep)

             if(line_data[0]=='sc_step'):
                 sc_step=line_data[1]
                 print('sc_step: ',sc_step)

             # Read the name of the position file
             if(line_data[0].strip()=='posfile'):
                positions,numbers=read_posfile(line_data[1])

   return lattice,positions,numbers,simid,mesh,timestep,sc_step,sc_nstep


# In[7]:


############################################################
# Open and read input files
############################################################
ifile='inpsd.dat'
lattice,positions,numbers,simid,mesh,timestep,sc_step,sc_nstep=read_inpsd(ifile)

############################################################
# Read adiabatic magnon spectra
############################################################
got_ams=os.path.isfile('ams.'+simid+'.out')
if got_ams:
    ams=np.loadtxt('ams.'+simid+'.out')
    ams_dist_col=ams.shape[1]-1


############################################################
# Read simulated dynamical structure factor (S(q,w))
############################################################
got_sqw=os.path.isfile('sqw.'+simid+'.out')
if got_sqw:
    #sqw=np.loadtxt('sqw.'+simid+'.out')
    sqw=np.genfromtxt('sqw.'+simid+'.out',usecols=(0,4,5,6,7,8))
    nq=int(sqw[-1,0])
    nw=int(sqw[-1,1])
    sqw_x=np.reshape(sqw[:,2],(nq,nw))
    sqw_y=np.reshape(sqw[:,3],(nq,nw))
    sqw_z=np.reshape(sqw[:,4],(nq,nw))
    #sqw_t=np.reshape(sqw[:,3],(nq,nw))
    sqw_t=sqw_x**2+sqw_y**2



############################################################
# Get the spacegroup from spglib and print relevant info
############################################################
cell=(lattice,positions,numbers)
spacegroup = spg.get_spacegroup(cell, symprec=1e-5)



# In[8]:


############################################################
# Get symmetry points from seekpath
############################################################
kpath_obj=spth.get_path(cell)
BZ=np.asarray(kpath_obj['reciprocal_primitive_lattice'])
sympoints=kpath_obj['point_coords']


# In[9]:


############################################################
# Read the qpoint-file used for the simulation
############################################################
qpts=np.genfromtxt('qfile.kpath',skip_header=1,usecols=(0,1,2))


# In[10]:


############################################################
# Extract symmetry points and their location for plotting
############################################################

axlab=[]
axidx=[]
axidx_abs=[]

for idx,row in enumerate(qpts):
    for k,v in kpath_obj['point_coords'].items():
        if (v==row).all():
            axlab.append(k[0])
            axidx.append(ams[idx,ams_dist_col])
            axidx_abs.append(ams[idx,0])
    
#axlab=['$\Gamma$' if x=='G' else '${}$'.format(x) for x in axlab]
axlab=['$\Gamma$' if x=='G' else '{}'.format(x) for x in axlab]


# In[12]:


############################################################
# Plot the AMS
############################################################
plt.figure(figsize=[8,5])


plt.plot(ams[:,ams_dist_col],ams[:,1:ams_dist_col])
ala=plt.xticks()


plt.xticks(axidx,axlab)
#plt.xlabel('q')
plt.ylabel('Energy (meV)')

plt.autoscale(tight=True)
plt.grid(b=True,which='major',axis='x')
plt.show()


# In[ ]:

############################################################
# Plot the S(q,w)
############################################################
fig = plt.figure(figsize=[8,5])
ax=plt.subplot(111)



hbar=4.135667662e-15
emax=0.5*float(hbar)/(float(timestep)*float(sc_step))*1e3
sqw_x[:,0]=sqw_x[:,0]/100.0
sqw_x=sqw_x.T/sqw_x.T.max(axis=0)
sqw_x=ndimage.gaussian_filter1d(sqw_x,sigma=1,axis=1,mode='constant')
sqw_x=ndimage.gaussian_filter1d(sqw_x,sigma=5,axis=0,mode='reflect')
plt.imshow(sqw_x, cmap=cmap.gist_ncar_r, interpolation='nearest',origin='lower',extent=[axidx_abs[0],axidx_abs[-1],0,emax])
plt.plot(ams[:,0]/ams[-1,0]*axidx_abs[-1],ams[:,1:ams_dist_col],'r')
ala=plt.xticks()


plt.xticks(axidx_abs,axlab)
plt.xlabel('q')
plt.ylabel('Energy (meV)')

plt.autoscale(tight=False)
ax.set_aspect('auto')
plt.grid(b=True,which='major',axis='x')
plt.show()






# In[ ]:





# In[ ]:




