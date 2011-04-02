from yt.mods import *
import pylab

### define problem name
problem_name = 'Toro-4-ShockTube'


### define simulation output directory and filename base
output_dir_base = 'DD'
datafile_base = 'data'


### define output to be plotted
dumpid = '0001'


### construct input filename
filename = './' + output_dir_base + dumpid + '/' + datafile_base + dumpid
print "Plotting output file %s\n" % filename


### some more filenames
exact_solution_filename = './Toro-4-ShockTube_t=0.035_exact.txt'
png_filename = './' + problem_name + '.png'



### load data
pf = load(filename)


### define InternalEnergy field
def _InternalEnergy(field, data):
    return data['TotalEnergy'] - 0.5*data['x-velocity']*data['x-velocity']

add_field('InternalEnergy', function=_InternalEnergy, units=r'\rm{erg}/\rm{g}')

### define Pressure field
def _Pressure(field, data):
    return (data.pf['Gamma'] - 1.0) * data['Density'] * data['InternalEnergy']

add_field('Pressure', function=_Pressure, units=r'\rm{dyne}/\rm{cm}^{2}')



### extract an ortho_ray (1D solution vector)
ray = pf.h.ortho_ray(0, [0.5, 0.5])


### define fields vector
fields = ('Density', 'x-velocity', 'InternalEnergy', 'Pressure' )

### read exact solution
exact = pylab.csv2rec( exact_solution_filename, delimiter=' ', names=('x', 'Density', 'x-velocity', 'Pressure', 'InternalEnergy') )


### calculate difference norm

# first interpolate the exact solution onto the ray
ray_exact = {'x': ray['x'], 
             'Density': pylab.stineman_interp(ray['x'],exact['x'],exact['Density']),
             'x-velocity': pylab.stineman_interp(ray['x'],exact['x'],exact['x-velocity']),
             'Pressure': pylab.stineman_interp(ray['x'],exact['x'],exact['Pressure']),
             'InternalEnergy': pylab.stineman_interp(ray['x'],exact['x'],exact['InternalEnergy'])}


# now calculate the norm (first order, since we're dealing with
# conservation laws)

norm = {}
maxnorm = {}
for f in fields:
    delta = ray['dx'] * abs( ray[f] - ray_exact[f] )
    norm[f] = delta.sum()
    maxnorm[f] = delta.max()

### make plot

pylab.figure(1, figsize=(8,7))

# define error norm text label function
def error_label(norm, maxnorm, xpos, ypos):
    vexp =  pylab.floor( pylab.log10( norm ) )
    if norm == 0: vexp = 0
    vmant = norm / 10**vexp
    thistext = r'$\Vert$E$\Vert_1 = \, %4.2f' % vmant
    if vexp != 0:
        thistext += r' \times \, 10^{%2d}' % vexp
    thistext += '$\n'

    vexp =  pylab.floor( pylab.log10( maxnorm ) )
    if maxnorm == 0: vexp = 0
    vmant = maxnorm / 10**vexp
    thistext += r'$\Vert$E$\Vert_\infty \! = \, %4.2f' % vmant
    if vexp != 0:
        thistext += r' \times \, 10^{%2d}' % vexp
    thistext += '$'

    pylab.text(xpos, ypos, thistext, va='top')  


# Density Plot
a = pylab.axes([0.09, 0.57, 0.38, 0.38])
pylab.axhline(0,color='k',linestyle='dotted')
pylab.plot(exact['x'],exact['Density'])
pylab.plot(ray['x'],ray['Density'], 'ro', ms=4)
#pylab.plot(ray_exact['x'],ray_exact['Density'], 'g+', ms=4)

pylab.axis([0,1,0.0,35.0])
pylab.xlabel('Position')
pylab.ylabel('Density')

error_label(norm['Density'], maxnorm['Density'], 0.1, 30.0)


# Velocity Plot
a = pylab.axes([0.59, 0.57, 0.38, 0.38])
pylab.axhline(0,color='k',linestyle='dotted')
pylab.plot(exact['x'],exact['x-velocity'])
pylab.plot(ray['x'],ray['x-velocity'], 'ro', ms=4)
#pylab.plot(ray_exact['x'],ray_exact['x-velocity'], 'g+', ms=4)

pylab.axis([0,1,-10.0,25.0])
pylab.xlabel('Position')
pylab.ylabel('Velocity')

error_label(norm['x-velocity'], maxnorm['x-velocity'], 0.5, 22.0)


# Pressure Plot
a = pylab.axes([0.09, 0.07, 0.38, 0.38])
pylab.axhline(0,color='k',linestyle='dotted')
pylab.plot(exact['x'],exact['Pressure'])
pylab.plot(ray['x'],ray['Pressure'], 'ro', ms=4)
#pylab.plot(ray_exact['x'],ray_exact['Pressure'], 'g+', ms=4)

pylab.axis([0,1,0.0,1800.0])
pylab.xlabel('Position')
pylab.ylabel('Pressure')

error_label(norm['Pressure'], maxnorm['Pressure'], 0.05, 1600.0)


# InternalEnergy Plot
a = pylab.axes([0.59, 0.07, 0.38, 0.38])
pylab.axhline(0,color='k',linestyle='dotted')
pylab.plot(exact['x'],exact['InternalEnergy'])
pylab.plot(ray['x'],ray['InternalEnergy'], 'ro', ms=4)
#pylab.plot(ray_exact['x'],ray_exact['InternalEnergy'], 'g+', ms=4)

pylab.axis([0,1,0.0,350.0])
pylab.xlabel('Position')
pylab.ylabel('Internal Energy')

error_label(norm['InternalEnergy'], maxnorm['InternalEnergy'], 0.1, 100.0)


### Save plot
pylab.savefig(png_filename)
