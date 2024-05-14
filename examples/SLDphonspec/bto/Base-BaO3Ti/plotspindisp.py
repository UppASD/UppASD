import numpy as np
import matplotlib.pyplot as pp
import matplotlib.colors as colors
import matplotlib as mpl
from matplotlib.collections import LineCollection
pp.style.use('default')
pp.rc("figure", facecolor="white")

# Load the adiabatic magnon spectra (AMS) from `ams.simid.out`
ams=np.genfromtxt("ams.bao3ti00.out")
na=ams.shape
print(na)

xmax=na[0]
ymax=0
for x in range(1,na[1]-1):
      ymax=max(ymax,max(ams[:,x]))
ymax=1.1*ymax

#Plot adiabatic magnon spectrum 
fig_ams = pp.figure()
ax1 = fig_ams.add_subplot(111)
ax1.set_xlabel('q-vector index')
ax1.set_ylabel('Energy (meV)')
for x in range(1,na[1]):
      ax1.plot(ams[:,na[1]-1],ams[:,x], c='k', linestyle='--')
ax1.axis('tight')
fig_ams.savefig('spindispersion.png', format='png', dpi=100)
