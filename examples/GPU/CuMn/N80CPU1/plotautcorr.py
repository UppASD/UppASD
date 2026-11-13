import numpy as np
import matplotlib.pyplot as pp
from matplotlib.pyplot import cm
pp.style.use('default')
pp.rc("figure", facecolor="white")

# Load the autocorrelation from `autocorr.simid.out`
autocorr=np.genfromtxt("autocorr.CuMn0000.out")
na=autocorr.shape
print(na)

# Load the waiting times t_w from acfile
tw=np.genfromtxt("acfile")
nb=tw.shape
print(nb)

# Rainbow color palette
color = iter(cm.rainbow(np.linspace(0, 1, na[1])))

# Plot the autocorrelation
fig_autocorr = pp.figure()
ax1 = fig_autocorr.add_subplot(111)
#ax1.set_title("Autocorrelation")
ax1.set_xlabel('time (simulation steps)')
ax1.set_ylabel('autocorrelation')
#ax1.set_xticks([])
for x in range(1,na[1]):
      # Plot the autocorrelation function for each waiting time
      c = next(color)
      ax1.plot(autocorr[:,0]-tw[x], autocorr[:,x], c=c, linestyle='--')
pp.xscale("log")
ax1.axis('tight')
pp.show()
fig_autocorr.savefig('autocorrelation.png', format='png', dpi=100)
