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

def is_close(A,B):
    check=True
    for a,b in zip(A,B):
        check=check and np.isclose(np.float64(a),np.float64(b),atol=1e-6)
    return check

# In[2]:

hbar=6.582119514e-13
ry_ev=13.605693009

sigma_q = 1.0
sigma_w = 1.0


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
   sc_emax = -1.0
   sc_eres = -1.0
   with open(ifile,'r') as infile:
      lines=infile.readlines()
      for idx,line in enumerate(lines):
         line_data=line.rstrip('\n').split()
         if len(line_data)>0:
             # Find the simulation id
             if(line_data[0]=='simid'):
                 simid=line_data[1]
                 #print('simid: ',simid)

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
                #print('cell: ',cell)
                #print('lattice: ',lattice)

             # Find the size of the simulated cell
             if(line_data[0]=='ncell'):
                ncell_x=int(line_data[1])
                ncell_y=int(line_data[1])
                ncell_z=int(line_data[1])
                mesh=[ncell_x,ncell_y,ncell_z]

             if(line_data[0]=='timestep'):
                 timestep=np.float64(line_data[1].replace('d','e'))
                 #print('timestep: ',timestep)

             if(line_data[0]=='sc_nstep'):
                 sc_nstep=line_data[1]
                 #print('sc_nstep: ',sc_nstep)

             if(line_data[0]=='sc_emax'):
                 sc_emax=np.float64(line_data[1])

             if(line_data[0]=='sc_eres'):
                 sc_eres=np.float64(line_data[1])
                 #print('sc_nstep: ',sc_nstep)

             if(line_data[0]=='sc_step'):
                 sc_step=line_data[1]

             if(line_data[0]=='qfile'):
                 qfile=line_data[1]
                 #print('sc_step: ',sc_step)

             # Read the name of the position file
             if(line_data[0].strip()=='posfile'):
                positions,numbers=read_posfile(line_data[1])

             if (sc_eres>0.0 and sc_emax>0.0):
                 sc_nstep = np.int64(sc_emax/sc_eres)
                 sc_step = np.int64(np.pi * hbar / (timestep * sc_emax * ry_ev))

   return lattice,positions,numbers,simid,mesh,timestep,sc_step,sc_nstep, qfile


# In[7]:


############################################################
# Open and read input files
############################################################
ifile='inpsd.dat'
lattice,positions,numbers,simid,mesh,timestep,sc_step,sc_nstep,qfile=read_inpsd(ifile)

############################################################
# Read adiabatic magnon spectra
############################################################
got_ams=os.path.isfile('ams.'+simid+'.out')
if got_ams:
    ams=np.loadtxt('ams.'+simid+'.out')
    ams_dist_col=ams.shape[1]-1
else: 
    got_ncams=os.path.isfile('ncams.'+simid+'.out')
    if got_ncams:
        ams=np.loadtxt('ncams.'+simid+'.out')
        ams_dist_col=ams.shape[1]-1
        got_ncams_pq=os.path.isfile('ncams+q.'+simid+'.out')
        if got_ncams_pq:
            ams_pq=np.loadtxt('ncams+q.'+simid+'.out')
        got_ncams_mq=os.path.isfile('ncams-q.'+simid+'.out')
        if got_ncams_mq:
            ams_mq=np.loadtxt('ncams-q.'+simid+'.out')
        

if got_ams or got_ncams:
    emax_lswt=1.10*np.amax(ams[:,1:ams_dist_col])
    #print('emax_lswt=',emax_lswt)


############################################################
# Read simulated dynamical structure factor (S(q,w))
############################################################
got_sqw=os.path.isfile('sqw.'+simid+'.out')
if got_sqw:
    #sqw=np.loadtxt('sqw.'+simid+'.out')
    sqw=np.genfromtxt('sqw.'+simid+'.out',usecols=(0,4,5,6,7,8))
    nq=int(sqw[-1,0])
    nw=int(sqw.shape[0]/nq)
    sqw_x=np.reshape(sqw[:,2],(nq,nw))
    sqw_y=np.reshape(sqw[:,3],(nq,nw))
    sqw_z=np.reshape(sqw[:,4],(nq,nw))
    sqw_t=sqw_x**2+sqw_y**2



############################################################
# Read simulated dynamical structure factor (S(q,w)) intensity
############################################################
got_sqw_int=os.path.isfile('sqwintensity.'+simid+'.out')
if got_sqw_int:
    sqw_t_int=np.genfromtxt('sqwintensity.'+simid+'.out',usecols=(0,4,5,6))
    nq_int=int(sqw_t_int[-1,0])
    nw_int=int(sqw_t_int.shape[0]/nq_int)
    sqw_int=np.reshape(sqw_t_int[:,2],(nq_int,nw_int))



############################################################
# Read linear spin wave theory structure factor (S(q,w)) intensity
############################################################
got_sqw_lint=os.path.isfile('ncsqw_intensity.'+simid+'.out')
if got_sqw_lint:
    sqw_lt_int=np.genfromtxt('ncsqw_intensity.'+simid+'.out',usecols=(0,4,5,6))
    nq_lint=int(sqw_lt_int[-1,0])
    nw_lint=int(sqw_lt_int.shape[0]/nq_lint)
    sqw_lint=np.reshape(sqw_lt_int[:,2],(nq_lint,nw_lint))



############################################################
# Read simulated dynamical structure factor (S(q,w))
############################################################
got_sqw_tens=os.path.isfile('sqwtensa.'+simid+'.out')
if got_sqw_tens:
    sqwt=np.genfromtxt('sqwtensa.'+simid+'.out',usecols=(0,4,5,6,7,8,9,10,11,12,13))
    nqt=int(sqwt[-1,0])
    nwt=int(sqwt.shape[0]/nqt)
    sqw_tens=np.zeros((nqt,nwt,3,3))
    sqw_tens[:,:,0,0]=np.reshape(sqwt[:,2] ,(nqt,nwt))
    sqw_tens[:,:,0,1]=np.reshape(sqwt[:,3] ,(nqt,nwt))
    sqw_tens[:,:,0,2]=np.reshape(sqwt[:,4] ,(nqt,nwt))
    sqw_tens[:,:,1,0]=np.reshape(sqwt[:,5] ,(nqt,nwt))
    sqw_tens[:,:,1,1]=np.reshape(sqwt[:,6] ,(nqt,nwt))
    sqw_tens[:,:,1,2]=np.reshape(sqwt[:,7] ,(nqt,nwt))
    sqw_tens[:,:,2,0]=np.reshape(sqwt[:,8] ,(nqt,nwt))
    sqw_tens[:,:,2,1]=np.reshape(sqwt[:,9] ,(nqt,nwt))
    sqw_tens[:,:,2,2]=np.reshape(sqwt[:,10],(nqt,nwt))
    #sqw_t=np.reshape(sqw[:,3],(nq,nw))
    #sqw_t=sqw_x**2+sqw_y**2


############################################################
# Read simulated dynamical structure factor (S(q,w))
############################################################
got_lswt_sqw_tens=os.path.isfile('ncsqw.'+simid+'.out')
if got_lswt_sqw_tens:
    lswt_sqwt=np.genfromtxt('ncsqw.'+simid+'.out',usecols=(0,4,5,6,7,8,9,10,11,12,13))
    lswt_nqt=int(lswt_sqwt[-1,0])
    lswt_nwt=int(lswt_sqwt.shape[0]/lswt_nqt)
    lswt_sqw_tens=np.zeros((lswt_nqt,lswt_nwt,3,3))
    lswt_sqw_tens[:,:,0,0]=np.reshape(lswt_sqwt[:,2] ,(lswt_nqt,lswt_nwt))
    lswt_sqw_tens[:,:,0,1]=np.reshape(lswt_sqwt[:,3] ,(lswt_nqt,lswt_nwt))
    lswt_sqw_tens[:,:,0,2]=np.reshape(lswt_sqwt[:,4] ,(lswt_nqt,lswt_nwt))
    lswt_sqw_tens[:,:,1,0]=np.reshape(lswt_sqwt[:,5] ,(lswt_nqt,lswt_nwt))
    lswt_sqw_tens[:,:,1,1]=np.reshape(lswt_sqwt[:,6] ,(lswt_nqt,lswt_nwt))
    lswt_sqw_tens[:,:,1,2]=np.reshape(lswt_sqwt[:,7] ,(lswt_nqt,lswt_nwt))
    lswt_sqw_tens[:,:,2,0]=np.reshape(lswt_sqwt[:,8] ,(lswt_nqt,lswt_nwt))
    lswt_sqw_tens[:,:,2,1]=np.reshape(lswt_sqwt[:,9] ,(lswt_nqt,lswt_nwt))
    lswt_sqw_tens[:,:,2,2]=np.reshape(lswt_sqwt[:,10],(lswt_nqt,lswt_nwt))
    #sqw_t=np.reshape(sqw[:,3],(nq,nw))
    #sqw_t=sqw_x**2+sqw_y**2



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

###############################################################
#### Extract symmetry points and their location for plotting
###############################################################

if(qfile=='qfile.kpath'):
    with open(qfile,'r') as f:
        f.readline()
        qpts=f.readlines()
    
    axlab=[]
    axidx=[]
    axidx_abs=[]
    for idx,row in enumerate(qpts):
        rs=row.split()
        if len(rs)==4:
            axlab.append(rs[3])
            axidx.append(ams[idx,ams_dist_col])
            axidx_abs.append(ams[idx,0])
    
    axlab=['$\Gamma$' if x[0]=='G' else '{}'.format(x) for x in axlab]
else:

   qpts=np.genfromtxt(qfile,skip_header=1,usecols=(0,1,2))
   
   axlab=[]
   axidx=[]
   axidx_abs=[]
   for idx,row in enumerate(qpts):
       for k,v in kpath_obj['point_coords'].items():
           if is_close(v,row):
               axlab.append(k[0])
               axidx.append(ams[idx,ams_dist_col])
               axidx_abs.append(ams[idx,0])
           elif idx==0:
               axlab.append(' ')
               axidx.append(0.0)
               axidx_abs.append(ams[idx,0])

       
   axlab=['$\Gamma$' if x[0]=='G' else '{}'.format(x) for x in axlab]

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
plt.ylim(0)
plt.grid(visible=True,which='major',axis='x')
#plt.show()
plt.savefig('ams.png')


# In[ ]:

############################################################
# Plot the S(q,w)
############################################################
if got_sqw:
    fig = plt.figure(figsize=[8,5])
    ax=plt.subplot(111)
    
    
    hbar=4.135667662e-15
    emax=0.5*np.float64(hbar)/(np.float64(timestep)*np.float64(sc_step))*1e3
    sqw_temp=(sqw_x**2+sqw_y**2)**0.5
    sqw_temp[:,0]=sqw_temp[:,0]/100.0
    sqw_temp=sqw_temp.T/sqw_temp.T.max(axis=0)
    sqw_temp=ndimage.gaussian_filter1d(sqw_temp,sigma=sigma_q,axis=1,mode='constant')
    sqw_temp=ndimage.gaussian_filter1d(sqw_temp,sigma=sigma_w,axis=0,mode='reflect')
    plt.imshow(sqw_temp, cmap=cmap.gist_ncar_r, interpolation='nearest',origin='lower',extent=[axidx_abs[0],axidx_abs[-1],0,emax])
    plt.plot(ams[:,0]/ams[-1,0]*axidx_abs[-1],ams[:,1:ams_dist_col],'black',lw=1)
    ala=plt.xticks()
    
    
    plt.xticks(axidx_abs,axlab)
    plt.xlabel('q')
    plt.ylabel('Energy (meV)')
    
    plt.autoscale(tight=False)
    ax.set_aspect('auto')
    plt.grid(visible=True,which='major',axis='x')
    #plt.show()
    plt.savefig('ams_sqw.png')


############################################################
# Plot the S(q,w) with full nc-LSWT support
############################################################
if got_sqw and got_ncams and got_ncams_pq and got_ncams_mq:
    fig = plt.figure(figsize=[8,5])
    ax=plt.subplot(111)
    
    
    hbar=4.135667662e-15
    emax=0.5*np.float64(hbar)/(np.float64(timestep)*np.float64(sc_step))*1e3
    sqw_temp=(sqw_x**2+sqw_y**2)**0.5
    sqw_temp[:,0]=sqw_temp[:,0]/100.0
    sqw_temp=sqw_temp.T/sqw_temp.T.max(axis=0)
    sqw_temp=ndimage.gaussian_filter1d(sqw_temp,sigma=sigma_q,axis=1,mode='constant')
    sqw_temp=ndimage.gaussian_filter1d(sqw_temp,sigma=sigma_w,axis=0,mode='reflect')
    plt.imshow(sqw_temp, cmap=cmap.gist_ncar_r, interpolation='nearest',origin='lower',extent=[axidx_abs[0],axidx_abs[-1],0,emax])
    plt.plot(ams[:,0]/ams[-1,0]*axidx_abs[-1],ams[:,1:ams_dist_col],'black',lw=1)
    plt.plot(ams[:,0]/ams[-1,0]*axidx_abs[-1],ams_pq[:,1:ams_dist_col],'black',lw=1)
    plt.plot(ams[:,0]/ams[-1,0]*axidx_abs[-1],ams_mq[:,1:ams_dist_col],'black',lw=1)
    ala=plt.xticks()
    
    
    plt.xticks(axidx_abs,axlab)
    plt.xlabel('q')
    plt.ylabel('Energy (meV)')
    
    plt.autoscale(tight=False)
    ax.set_aspect('auto')
    plt.grid(visible=True,which='major',axis='x')
    plt.savefig('ams_sqw_q.png')



############################################################
# Plot the simulated S(q,w) intensity
############################################################
if got_sqw_int:
    fig = plt.figure(figsize=[8,5])
    ax=plt.subplot(111)
    
    
    
    hbar=4.135667662e-15
    emax=0.5*np.float64(hbar)/(np.float64(timestep)*np.float64(sc_step))*1e3
    sqw_temp=sqw_int
    sqw_temp[:,0]=sqw_temp[:,0]/100.0
    sqw_temp=sqw_temp.T/sqw_temp.T.max(axis=0)
    sqw_temp=ndimage.gaussian_filter1d(sqw_temp,sigma=sigma_q,axis=1,mode='constant')
    sqw_temp=ndimage.gaussian_filter1d(sqw_temp,sigma=sigma_w,axis=0,mode='reflect')
    plt.imshow(sqw_temp, cmap=cmap.gist_ncar_r, interpolation='nearest',origin='lower',extent=[axidx_abs[0],axidx_abs[-1],0,emax])
    plt.plot(ams[:,0]/ams[-1,0]*axidx_abs[-1],ams[:,1:ams_dist_col],'black',lw=1)
    ala=plt.xticks()
    
    
    plt.xticks(axidx_abs,axlab)
    plt.xlabel('q')
    plt.ylabel('Energy (meV)')
    
    plt.autoscale(tight=False)
    ax.set_aspect('auto')
    plt.grid(visible=True,which='major',axis='x')
    #plt.show()
    plt.savefig('ams_sqw_int.png')



############################################################
# Plot the LSWT S(q,w) intensity
############################################################
if got_sqw_lint:
    fig = plt.figure(figsize=[8,5])
    ax=plt.subplot(111)
    
    
    
    hbar=4.135667662e-15
    emax=0.5*np.float64(hbar)/(np.float64(timestep)*np.float64(sc_step))*1e3
    sqw_temp=sqw_lint
    sqw_temp[:,0]=sqw_temp[:,0]/100.0
    sqw_temp=sqw_temp.T/sqw_temp.T.max(axis=0)
    sqw_temp=ndimage.gaussian_filter1d(sqw_temp,sigma=sigma_q,axis=1,mode='constant')
    sqw_temp=ndimage.gaussian_filter1d(sqw_temp,sigma=sigma_w,axis=0,mode='reflect')
    plt.imshow(sqw_temp, cmap=cmap.gist_ncar_r, interpolation='nearest',origin='lower',extent=[axidx_abs[0],axidx_abs[-1],0,emax_lswt])
    plt.plot(ams[:,0]/ams[-1,0]*axidx_abs[-1],ams[:,1:ams_dist_col],'black',lw=1)
    ala=plt.xticks()
    
    
    plt.xticks(axidx_abs,axlab)
    plt.xlabel('q')
    plt.ylabel('Energy (meV)')
    
    plt.autoscale(tight=False)
    ax.set_aspect('auto')
    plt.grid(visible=True,which='major',axis='x')
    #plt.show()
    plt.savefig('ams_ncsqw_int.png')




xyz=('x','y','z')

############################################################
# Plot the S(q,w) on tensorial form (diagonal terms)
############################################################
if got_sqw_tens:
    fig = plt.figure(figsize=[16,4])
    hbar=4.135667662e-15
    emax=0.5*np.float64(hbar)/(np.float64(timestep)*np.float64(sc_step))*1e3

    plt.xticks(axidx_abs,axlab)
    plt.xlabel('q')
    plt.ylabel('Energy (meV)')
    
    #plt.autoscale(tight=False)
    ax.set_aspect('auto')
    plt.grid(visible=True,which='major',axis='x')

    #sqw_x[:,0]=sqw_x[:,0]/100.0
    #sqw_x=sqw_x.T/sqw_x.T.max(axis=0)
    plt_idx=130
    for ix in range(3):
        iy=ix
        sqw_temp=sqw_tens[:,:,ix,iy]
        sqw_temp[:,0]=sqw_temp[:,0]/1e5
        sqw_temp=sqw_temp.T/sqw_temp.T.max(axis=0)
        plt_idx=plt_idx+1
        ax=plt.subplot(plt_idx)
        stitle='$S^{'+xyz[ix]+xyz[iy]+'}(q,\omega)}$'
        ax.set_title(stitle)

        sqw_temp=ndimage.gaussian_filter1d(sqw_temp,sigma=sigma_q,axis=1,mode='constant')
        sqw_temp=ndimage.gaussian_filter1d(sqw_temp,sigma=sigma_w,axis=0,mode='reflect')
        plt.imshow(sqw_temp, cmap=cmap.gist_ncar_r, interpolation='nearest',origin='lower',extent=[axidx_abs[0],axidx_abs[-1],0,emax],aspect='auto')
        plt.plot(ams[:,0]/ams[-1,0]*axidx_abs[-1],ams[:,1:ams_dist_col],'black',lw=  1)
        plt.xticks(axidx_abs,axlab)

    plt.tight_layout()

    plt.savefig('sqw_diagonal.png')



############################################################
# Plot the S(q,w) on tensorial form (full tensor)
############################################################
if got_sqw_tens:
    fig = plt.figure(figsize=[16,10])
    hbar=4.135667662e-15
    emax=0.5*np.float64(hbar)/(np.float64(timestep)*np.float64(sc_step))*1e3

    plt.xticks(axidx_abs,axlab)
    plt.xlabel('q')
    plt.ylabel('Energy (meV)')
    
    #plt.autoscale(tight=False)
    ax.set_aspect('auto')
    plt.grid(visible=True,which='major',axis='x')

    #sqw_x[:,0]=sqw_x[:,0]/100.0
    #sqw_x=sqw_x.T/sqw_x.T.max(axis=0)
    plt_idx=330
    for ix in range(3):
        for iy in range(3):
            sqw_temp=sqw_tens[:,:,ix,iy]
            sqw_temp[:,0]=sqw_temp[:,0]/1e5
            sqw_temp=sqw_temp.T/sqw_temp.T.max(axis=0)
            plt_idx=plt_idx+1
            ax=plt.subplot(plt_idx)
            stitle='$S^{'+xyz[ix]+xyz[iy]+'}(q,\omega)}$'
            ax.set_title(stitle)

            sqw_temp=ndimage.gaussian_filter1d(sqw_temp,sigma=sigma_q,axis=1,mode='constant')
            sqw_temp=ndimage.gaussian_filter1d(sqw_temp,sigma=sigma_w,axis=0,mode='reflect')
            plt.imshow(sqw_temp, cmap=cmap.gist_ncar_r, interpolation='nearest',origin='lower',extent=[axidx_abs[0],axidx_abs[-1],0,emax],aspect='auto')
            plt.plot(ams[:,0]/ams[-1,0]*axidx_abs[-1],ams[:,1:ams_dist_col],'black',lw=  1)
            plt.xticks(axidx_abs,axlab)

    plt.tight_layout()

    plt.savefig('sqw_tensor.png')



############################################################
# Plot the LSWT S(q,w) on tensorial form (only diagonal terms)
############################################################
if got_lswt_sqw_tens:
    fig = plt.figure(figsize=[16,4])
    hbar=4.135667662e-15
    emax=0.5*np.float64(hbar)/(np.float64(timestep)*np.float64(sc_step))*1e3

    plt.xticks(axidx_abs,axlab)
    plt.xlabel('q')
    plt.ylabel('Energy (meV)')
    
    #plt.autoscale(tight=False)
    #ax.set_aspect('auto')
    plt.grid(visible=True,which='major',axis='x')

    #sqw_x[:,0]=sqw_x[:,0]/100.0
    #sqw_x=sqw_x.T/sqw_x.T.max(axis=0)
    plt_idx=130
    for ix in range(3):
        iy=ix
        sqw_temp=lswt_sqw_tens[:,:,ix,iy]
        sqw_temp=sqw_temp.T/(sqw_temp.T.max(axis=0)+1.0e-20)
        plt_idx=plt_idx+1
        ax=plt.subplot(plt_idx)
        stitle='$S^{'+xyz[ix]+xyz[iy]+'}_{LSWT}$'
        ax.set_title(stitle)

        sqw_temp=ndimage.gaussian_filter1d(sqw_temp,sigma=sigma_q,axis=1,mode='constant')
        sqw_temp=ndimage.gaussian_filter1d(sqw_temp,sigma=sigma_w,axis=0,mode='reflect')
        plt.imshow(sqw_temp, cmap=cmap.gist_ncar_r, interpolation='nearest',origin='lower',extent=[axidx_abs[0],axidx_abs[-1],0,emax_lswt],aspect='auto')
        #plt.plot(ams[:,0]/ams[-1,0]*axidx_abs[-1],ams[:,1:ams_dist_col],'black',lw=  1)
        plt.xticks(axidx_abs,axlab)

    plt.tight_layout()

    plt.savefig('ncsqw_diagonal.png')


############################################################
# Plot the LSWT S(q,w) on tensorial form (full tensor)
############################################################
if got_lswt_sqw_tens:
    fig = plt.figure(figsize=[16,10])
    hbar=4.135667662e-15
    emax=0.5*np.float64(hbar)/(np.float64(timestep)*np.float64(sc_step))*1e3

    plt.xticks(axidx_abs,axlab)
    plt.xlabel('q')
    plt.ylabel('Energy (meV)')
    
    #plt.autoscale(tight=False)
    #ax.set_aspect('auto')
    plt.grid(visible=True,which='major',axis='x')

    #sqw_x[:,0]=sqw_x[:,0]/100.0
    #sqw_x=sqw_x.T/sqw_x.T.max(axis=0)
    plt_idx=330
    for ix in range(3):
        for iy in range(3):
            sqw_temp=lswt_sqw_tens[:,:,ix,iy]
            sqw_temp=sqw_temp.T/(sqw_temp.T.max(axis=0)+1.0e-20)
            plt_idx=plt_idx+1
            ax=plt.subplot(plt_idx)
            stitle='$S^{'+xyz[ix]+xyz[iy]+'}_{LSWT}$'
            ax.set_title(stitle)

            sqw_temp=ndimage.gaussian_filter1d(sqw_temp,sigma=sigma_q,axis=1,mode='constant')
            sqw_temp=ndimage.gaussian_filter1d(sqw_temp,sigma=sigma_w,axis=0,mode='reflect')
            plt.imshow(sqw_temp, cmap=cmap.gist_ncar_r, interpolation='nearest',origin='lower',extent=[axidx_abs[0],axidx_abs[-1],0,emax_lswt],aspect='auto')
            #plt.plot(ams[:,0]/ams[-1,0]*axidx_abs[-1],ams[:,1:ams_dist_col],'black',lw=  1)
            plt.xticks(axidx_abs,axlab)

    plt.tight_layout()

    plt.savefig('ncsqw_tensor.png')

