from __future__ import print_function
import yt
import pylab
import numpy as numpy
### define simulation output directory and filename base
output_dir_base = 'DD'
datafile_base = 'data'

### load data
ts = yt.DatasetSeries.from_filenames("*/*.hierarchy")
for pf in ts:
    pylab.clf()
    print(pf.current_time)

    ### extract an ortho_ray (1D solution vector)
    ray = pf.ortho_ray(0, [0.5, 0.5])
    ray_sort = numpy.argsort(ray["x"])
    pylab.figure(1, figsize=(10,8))

    # Density Plot
    pylab.subplot(2,2,1)
    pylab.semilogy(ray['x'][ray_sort],ray['density'][ray_sort])
    pylab.xlabel('Position')
    pylab.ylabel('Density')

    # Temperature Plot
    pylab.subplot(2,2,2)
    pylab.semilogy(ray['x'][ray_sort],ray['Temperature'][ray_sort], 'b')
    pylab.xlabel('Position')
    pylab.ylabel('Temperature')

    # Mach Plot
    pylab.subplot(2,2,3)
    pylab.plot(ray['x'][ray_sort],ray['Mach'][ray_sort], 'k')
    pylab.xlabel('x')
    pylab.ylabel('Mach')
    
    # Mach Plot
    pylab.subplot(2,2,4)
    pylab.plot(ray['x'][ray_sort],ray[('gas','velocity_magnitude')][ray_sort], 'k')
    pylab.xlabel('x')
    pylab.ylabel('|v|')

    ### Save plot
    pylab.tight_layout()
    pylab.savefig('%s_thermal.png' % pf)

pylab.clf()

