import numpy as np
import matplotlib.pyplot as pp
import matplotlib.colors as colors
import matplotlib as mpl
from matplotlib.collections import LineCollection
pp.style.use('default')
pp.rc("figure", facecolor="white")

# Load the phonon spectra from `phondisp.simid.out`
ams=np.genfromtxt("ams.diamond0.out")
na=ams.shape
print(na)
phondisp=np.genfromtxt("phondisp.diamond0.out")
nb=phondisp.shape
print(nb)

xmax=na[0]
ymax=0
for x in range(1,na[1]-1):
      ymax=max(ymax,max(phondisp[:,x]))
ymax=1.1*ymax

#Plot phonon spectrum
fig_phondisp = pp.figure()
ax1 = fig_phondisp.add_subplot(111)
ax1.set_xlabel('q-vector index')
ax1.set_ylabel('Energy (meV)')
for x in range(1,nb[1]):
      ax1.plot(ams[:,na[1]-1],phondisp[:,x], c='k', linestyle='--')
ax1.axis('tight')
fig_phondisp.savefig('phonondispersion.png', format='png', dpi=100)
