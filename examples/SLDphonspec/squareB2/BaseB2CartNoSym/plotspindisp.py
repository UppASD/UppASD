import numpy as np
import matplotlib.pyplot as pp
import matplotlib.colors as colors
import matplotlib as mpl
from matplotlib.collections import LineCollection
pp.style.use('default')
pp.rc("figure", facecolor="white")

# Load the adiabatic magnon spectra (AMS) from `ams.simid.out`
ams=np.genfromtxt("ams.square2D.out")
na=ams.shape
print(na)

#Plot adiabatic magnon spectrum 
fig_ams3 = pp.figure()
ax3 = fig_ams3.add_subplot(111)
ax3.set_xlabel('q-vector index')
ax3.set_ylabel('Energy (meV)')
ax3.plot(ams[:,na[1]-1],ams[:,1], c='k', linestyle='-', marker='s', ms='4')
ax3.plot(ams[:,na[1]-1],ams[:,2], c='r', linestyle='--', marker='s', ms='4')
ax3.axis('tight')
pp.xticks((0.00, 3.16, 6.32, 11.35), ('$\Gamma$', '$X$', '$M$', '$\Gamma$'), fontsize='x-large')
#pp.show()
fig_ams3.savefig('spindispersion.png', format='png', dpi=100)
