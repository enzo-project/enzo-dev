from matplotlib import use; use('Agg')
import matplotlib.pyplot as plt
from yt.mods import *

pf = load('DD0001/data0001')
ray = pf.h.ortho_ray(0,(0,0))

plt.subplot(311)
plt.semilogy(ray['x'], ray['Density'], c='k', ls='-', marker='.', ms=5)
plt.ylabel(r'Density (cm$^{-3}$)')

plt.subplot(312)
plt.semilogy(ray['x'], ray['Temperature'], c='k', ls='-', marker='.', ms=5)
plt.ylabel('Temperature (K)')

plt.subplot(313)
plt.plot(ray['x'], 1e-5*ray['x-velocity'], c='k', ls='-', marker='.', ms=5)
plt.ylabel('Velocity (km/s)')
plt.xlabel('x')

plt.savefig('AMRZeldovichPancake.png')
