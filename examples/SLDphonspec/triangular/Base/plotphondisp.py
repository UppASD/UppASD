import numpy as np
import matplotlib.pyplot as pp
import matplotlib.colors as colors
import matplotlib as mpl
from matplotlib.collections import LineCollection
pp.style.use('default')
pp.rc("figure", facecolor="white")

# Load the phonon spectra from `phondisp.simid.out`
ams=np.genfromtxt("ams.triang2D.out")
na=ams.shape
print(na)
phondisp=np.genfromtxt("phondisp.triang2D.out")
nb=phondisp.shape
print(nb)

#Plot phonon spectrum
fig_phondisp3 = pp.figure()
ax3 = fig_phondisp3.add_subplot(111)
#ax3.set_title("Phonon spectrum")
ax3.set_xlabel('q-vector index')
ax3.set_ylabel('Energy (meV)')
#ax3.set_xticks([])
ax3.plot(ams[:,na[1]-1],phondisp[:,1], c='k', linestyle='-', marker='s', ms='4')
ax3.plot(ams[:,na[1]-1],phondisp[:,2], c='r', linestyle='--', marker='s', ms='4')
ax3.plot(ams[:,na[1]-1],phondisp[:,3], c='g', linestyle='-.', marker='s', ms='4')
ax3.axis('tight')
pp.xticks((0.00, 10.00, 15.00, 23.660), ('$\Gamma$', '$K$', '$M$', '$\Gamma$'), fontsize='x-large')
#pp.show()
fig_phondisp3.savefig('phonondispersion.png', format='png', dpi=100)
