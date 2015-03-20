#
#  This script compares the position of the test particle to grid positions
#
#  There should be one grid per level and all grids above the root grid should be
#  centered on the test particle.
#
#  Additional useful scripts for analyzing this problem can be found in run/GravitySolver/TestOrbit
#
#  This script is written for yt3.0
#
import yt
from yt.utilities.linear_interpolators import *
import numpy as na
import matplotlib.pyplot as plt

delta_grid = []
grid_vol = []
time = []

ts = yt.DatasetSeries("*/*.hierarchy")
for ds in ts:
    
    if ds.index.num_grids != ds.index.max_level+1:
        print 'Something is wrong, there is not one grid per level',ds.index.num_grids,ds.index.max_level+1
        break

    ad = ds.all_data()

    # particle position
    xp = ad['particle_position_x'][1]
    yp = ad['particle_position_y'][1]
    zp = ad['particle_position_z'][1]
    
    # there should be 3 grids, one on each level, 0 through 2
    # compute the difference between the particle position & grid centers

    delta_grid_temp = []
    grid_vol_temp = []
    for gn in np.arange(ds.index.num_grids-1):
        le = ds.index.grid_left_edge[gn+1]
        re = ds.index.grid_right_edge[gn+1]
        dim = ds.index.grid_dimensions[gn+1]
        dx = (re[0]-le[0])/dim[0]
        center = (le+re)*0.5

        r = np.sqrt((center[0]-xp)**2.0 + (center[1]-yp)**2.0 + (center[2]-zp)**2.0)/dx
        ncells = dim[0]*dim[1]*dim[2]
        
        grid_vol_temp.append(ncells)
        delta_grid_temp.append(r)

    # ignore first value
    if ds.current_time > 0:
        delta_grid.append(delta_grid_temp)
        grid_vol.append(grid_vol_temp)
        time.append(ds.current_time)
    

# turn into NumPy arrays for convenient manipulation
delta_grid = na.array(delta_grid)
grid_vol = na.array(grid_vol)
time = na.array(time)

for gn in range(ds.parameters['MaximumRefinementLevel']):
   
    plt.plot(time,delta_grid[:,gn],'-o',label='Level '+str(gn+1)+' Grid')
    print 'Grid on Level '+str(gn+1)
    print '   mean grid offset from particle:',na.mean(delta_grid[:,gn]),'cells'
    print '   mean grid volume:',na.mean(grid_vol[:,gn]),'cells^3'

plt.xlabel('Orbit Number')
plt.ylabel('|Center of Grid - Particle Position| / Cell Width')
plt.xlim(1,21)
plt.legend(loc=0)
plt.savefig('GridPositions.png')
plt.close()
