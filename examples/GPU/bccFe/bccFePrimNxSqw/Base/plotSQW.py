import numpy as np
import matplotlib.pyplot as pp
import matplotlib.colors as colors
import matplotlib as mpl
from matplotlib.collections import LineCollection
pp.style.use('default')
pp.rc("figure", facecolor="white")

# Load the adiabatic magnon spectra (AMS) from `ams.simid.out`
ams=np.genfromtxt("ams.bccFe000.out")
na=ams.shape
print(na)

# Load the S(q,w) intensity from `sqw.simid.out`
sqwa = np.genfromtxt("sqw.bccFe000.out")
qd=int(sqwa[sqwa.shape[0]-1,0])
ed=int(sqwa[sqwa.shape[0]-1,4])

# Define useful constants and parameters (like h for s(q,w))
h_mev=4.135667662e-12
hbar_mev=6.582119514e-13
timestep=1.0e-16
sc_step=10
sc_nstep=ed
delta_f=1.0/(sc_nstep*sc_step*timestep)/2
delta_w=2.0*np.pi/(sc_nstep*sc_step*timestep)/2

hf_scale=h_mev/(timestep*sc_step)
y_scale=hbar_mev*delta_w*sc_nstep/2

sqw_x=np.transpose((np.reshape(sqwa[:,5],(qd,ed))[:,0:int(ed/2)]))
sqw_y=np.transpose((np.reshape(sqwa[:,6],(qd,ed))[:,0:int(ed/2)]))
sqw_z=np.transpose((np.reshape(sqwa[:,7],(qd,ed))[:,0:int(ed/2)]))
sqw_a=np.transpose((np.reshape(sqwa[:,8],(qd,ed))[:,0:int(ed/2)]))
#sqwint=sqw_x+sqw_y
sqwint=np.log10(np.array(sqw_x+sqw_y))

# High-symmetry points: read the '#' labels in the qfile (q-point index = line number)
symlabels={'gamma':'$\\Gamma$','g':'$\\Gamma$'}
xticks=[]
xticklabels=[]
qfile=[l.split()[1] for l in open("inpsd.dat") if l.split()[:1]==['qfile']][0]
with open(qfile) as qf:
      qf.readline()
      for iq,line in enumerate(qf,start=1):
            if '#' in line:
                  lab=line.split('#')[1].strip()
                  xticks.append(iq)
                  xticklabels.append(symlabels.get(lab.lower(),'$'+lab+'$'))

xmax=qd
ymax=0
for x in range(1,na[1]-1):
      ymax=max(ymax,max(ams[:,x]))
ymax=1.1*ymax

#Plot adiabatic magnon spectrum 
fig_ams = pp.figure()
ax1 = fig_ams.add_subplot(111)
#ax1.set_title("Adiabatic magnon spectrum")    
ax1.set_xlabel('q-vector index')
ax1.set_ylabel('Energy (meV)')
#ax1.set_xticks([])
for x in range(1,na[1]-1):
      ax1.plot(ams[:,0],ams[:,x], c='k', linestyle='--')
ax1.axis('tight')
pp.xticks(xticks, xticklabels, fontsize='x-large')
#pp.show()
fig_ams.savefig('spinwave.eps', format='eps', dpi=100)
      
# Plot the dynamic correlation function $S(q,\omega)$ 
fig_sqw=pp.figure()
#pp.title('$S(Q,\omega)$')
pp.xlabel(r'$\mathbf{Q} \quad \mathrm{vector}$', fontsize='x-large')
pp.ylabel(r'$\mathrm{Energy \quad (meV})$', fontsize='x-large')
pp.ylim(0,ymax)
cmap1 = colors.LinearSegmentedColormap.from_list("", ["blue","white","red"])   
surf1=pp.imshow(sqwint,origin='lower',extent=[0.5,qd+0.5,0,y_scale],aspect='auto',cmap=cmap1)
fig_sqw.colorbar(surf1, shrink=0.5, aspect=5)
pp.xticks(xticks, xticklabels, fontsize='x-large')
#pp.show()
fig_sqw.savefig('dynstructfac.eps', format='eps', dpi=100, bbox_inches='tight')

# Plot the adiabatic spectra and the magnon dynamic correlation function $S(q,\omega)$ 
fig_amssqw, ax2 = pp.subplots(figsize=(6,4))
cmap2 = colors.LinearSegmentedColormap.from_list("", ["blue","white","red"])   
surf2=ax2.imshow(sqwint,origin='lower',extent=[0.5,qd+0.5,0,y_scale],aspect='auto',cmap=cmap2)
for x in range(1,na[1]-1):
      ax2.plot(ams[0:xmax,0],ams[0:xmax,x], c='k', linestyle='--')
ax2.set_xlim(0.5, xmax+0.5)
ax2.set_ylim(0, ymax)
ax2.set_xlabel(r'$\mathbf{Q} \quad \mathrm{vector}$', fontsize='x-large')
ax2.set_ylabel(r'$\mathrm{Energy \quad (meV})$', fontsize='x-large')
ax2.set_xticks(xticks, xticklabels, fontsize='x-large')
fig_amssqw.colorbar(surf2, shrink=0.5, aspect=5)
#cbar=fig_amssqw.colorbar(surf2, shrink=1.0, aspect=10)
#cbar.set_label(r'$\mathrm{Log}_{10}[S(\mathbf{Q},\omega)/S_0(\mathbf{Q},\omega)]$', fontsize='x-large')
#pp.show()
fig_amssqw.savefig('amssqw.eps', format='eps', dpi=100, bbox_inches='tight')
