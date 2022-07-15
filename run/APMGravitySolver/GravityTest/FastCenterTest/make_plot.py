import numpy as na
import matplotlib.pyplot as plt
import glob


legend_added = False
for f in glob.glob('TestGravityCheckResults*'):
    Data = na.loadtxt(f)

    # Read in data from text file produced by code
    radius = Data[:,4]
    ForceTangentialComputed = Data[:,5]
    ForceRadialComputed = Data[:,6]
    ForceRadialTrue = Data[:,7]

    plt.loglog(radius, ForceRadialComputed, 'b', label='Frad', ls='None', marker='+')
    plt.loglog(radius, ForceRadialTrue, 'k', label='Frad_true', ls='None', marker='.')
    plt.loglog(radius, ForceTangentialComputed, 'm', label='Ftang', ls='None', marker='.')
    if not legend_added:
        plt.legend()
        legend_added = True

plt.xlabel('r (cells)')
plt.ylabel('Force')
plt.axis([5e-3,0.5,1e-4,1e3])
plt.savefig('GravityTest.png')
