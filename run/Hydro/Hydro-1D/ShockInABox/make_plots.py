from yt.mods import *
import pylab

### define simulation output directory and filename base
output_dir_base = 'DD'
datafile_base = 'data'

### load data
ts = TimeSeriesData.from_filenames("*/*.hierarchy")
for pf in ts:
    pylab.clf()
    print(pf.current_time)

    ### extract an ortho_ray (1D solution vector)
    ray = pf.h.ortho_ray(0, [0.5, 0.5])

    pylab.figure(1, figsize=(10,8))

    # Density Plot
    pylab.subplot(2,2,1)
    pylab.semilogy(ray['x'],ray['Density'], 'k')
    pylab.xlabel('Position')
    pylab.ylabel('Density')

    # Temperature Plot
    pylab.subplot(2,2,2)
    pylab.semilogy(ray['x'],ray['Temperature'], 'b')
    pylab.xlabel('Position')
    pylab.ylabel('Temperature')

    # Mach Plot
    pylab.subplot(2,2,3)
    pylab.plot(ray['x'],ray['Mach'], 'k')
    pylab.xlabel('x')
    pylab.ylabel('Mach')
    
    # Mach Plot
    pylab.subplot(2,2,4)
    pylab.plot(ray['x'],ray['VelocityMagnitude'], 'k')
    pylab.xlabel('x')
    pylab.ylabel('|v|')

    ### Save plot
    pylab.savefig('%s_thermal.png' % pf)

pylab.clf()

