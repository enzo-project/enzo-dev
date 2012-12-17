from yt.mods import *
import os
import sys
import pylab

def make_plot(pfname):
    pf = load(pfname)
    ### extract an ortho_ray (1D solution vector)
    ray = pf.h.ortho_ray(0, [0.5, 0.5])

    ### define fields vector
    fields = ('Density', 'x-velocity', 'TotalEnergy', 'Pressure' )

    ### make plot

    pylab.figure(1, figsize=(8,7))

    # Density Plot
    a = pylab.axes([0.09, 0.57, 0.38, 0.38])
    pylab.axhline(0,color='k',linestyle='dotted')
    pylab.plot(ray['x'],ray['Density'], 'ro', ms=4)

    pylab.xlabel('Position')
    pylab.ylabel('Density')

    # Velocity Plot
    a = pylab.axes([0.59, 0.57, 0.38, 0.38])
    pylab.axhline(0,color='k',linestyle='dotted')
    pylab.plot(ray['x'],ray['x-velocity'], 'ro', ms=4)

    pylab.xlabel('Position')
    pylab.ylabel('Velocity')

    # TotalEnergy Plot
    a = pylab.axes([0.59, 0.07, 0.38, 0.38])
    pylab.axhline(0,color='k',linestyle='dotted')
    pylab.plot(ray['x'],ray['TotalEnergy'], 'ro', ms=4)

    pylab.xlabel('Position')
    pylab.ylabel('Total Energy')

    ### Save plot
    pylab.savefig('%s.png' % pf)
    pylab.clf()
    
if __name__ == '__main__':
    for i in range(11):
        try: 
            make_plot('DD%04i/data%04i'% (i,i))
        except:
            break

    # To make a movie using avconv, uncomment the following 2 lines
    # os.system('avconv -r 10 -i data%04d.png -threads 8 -pass 1 -an -f webm -b 2000k movie.webm')
    # os.system('avconv -r 10 -i data%04d.png -threads 8 -pass 2 -an -f webm -b 2000k movie.webm')

