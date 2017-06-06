from yt.mods import *
import pylab
import sys

### define problem name
problem_name = 'CRShockTube_'
#problem_name = 'SodShockTubeAMR'


### define simulation output directory and filename base
output_dir_base = 'DD'
datafile_base = 'data'


### define output to be plotted
dumpid = '0001'

problem_name += dumpid

### construct input filename
filename = './' + output_dir_base + dumpid + '/' + datafile_base + dumpid
print("Plotting output file %s\n" % filename)


### some more filenames
exact_solution_filename = 'analytic.txt'
png_filename = './' + problem_name + '.png'



### load data
pf = load(filename)


### define InternalEnergy field
def _InternalEnergy(field, data):
    return data['TotalEnergy']# - 0.5*data['x-velocity']*data['x-velocity']

add_field('InternalEnergy', function=_InternalEnergy, units=r'\rm{erg}/\rm{g}')

### define Pressure field
def _Pressure(field, data):
    return (data.pf['Gamma'] - 1.0) * data['Density'] * data['InternalEnergy']

add_field('Pressure', function=_Pressure, units=r'\rm{dyne}/\rm{cm}^{2}')

### define CR field
def _CRDensity(field, data):
    return data.pf['CREnergyDensity']

add_field('CREnergyDensity', function=_CRDensity, units=r'\rm{idk}')

### define CR Pressure
def _CRPressure(field,data):
    return (1.0/3.0)*data.pf['CREnergyDensity']

add_field('CRPressure',function=_CRPressure,units=r'\rm{idk}')

### extract an ortho_ray (1D solution vector)
ray = pf.h.ortho_ray(0, [0.5, 0.5])


### define fields vector
fields = ('Density', 'x-velocity', 'InternalEnergy', 'Pressure' )

### read exact solution
exact = pylab.csv2rec( exact_solution_filename, delimiter=' ', names=('x', 'Density', 'x-velocity', 'Pressure', 'InternalEnergy') )


### calculate difference norm

# first interpolate the exact solution onto the ray
#ray_exact = {'x': ray['x'], 
#             'Density': pylab.stineman_interp(ray['x'],exact['x'],exact['Density']),
#             'x-velocity': pylab.stineman_interp(ray['x'],exact['x'],exact['x-velocity']),
#             'Pressure': pylab.stineman_interp(ray['x'],exact['x'],exact['Pressure']),
#             'InternalEnergy': pylab.stineman_interp(ray['x'],exact['x'],exact['InternalEnergy'])}


# now calculate the norm (first order, since we're dealing with
# conservation laws)

#norm = {}
#maxnorm = {}
#for f in fields:
#    delta = ray['dx'] * abs( ray[f] - ray_exact[f] )
#    norm[f] = delta.sum()
#    maxnorm[f] = delta.max()


### make plot

pylab.figure(1, figsize=(8,7))

# define error norm text label function
#def error_label(norm, maxnorm, xpos, ypos):
#    vexp =  pylab.floor( pylab.log10( norm ) )
#    if norm == 0: vexp = 0
#    vmant = norm / 10**vexp
#    thistext = r'$\Vert$E$\Vert_1 = \, %4.2f' % vmant
#    if vexp != 0:
#        thistext += r' \times \, 10^{%2d}' % vexp
#    thistext += '$\n'
#
#    vexp =  pylab.floor( pylab.log10( maxnorm ) )
#    if maxnorm == 0: vexp = 0
#    vmant = maxnorm / 10**vexp
#    thistext += r'$\Vert$E$\Vert_\infty \! = \, %4.2f' % vmant
#    if vexp != 0:
#        thistext += r' \times \, 10^{%2d}' % vexp
#    thistext += '$'
#
#    pylab.text(xpos, ypos, thistext, va='top')  


# Density Plot
a = pylab.axes([0.09, 0.57, 0.38, 0.38])
pylab.plot((0,250,250,500),(1.0,1.0,0.2,0.2),color='k',linestyle='dotted')
pylab.plot(exact['x'],exact['Density'],color='r')
pylab.plot(ray['x'],ray['Density'], color='k')
#pylab.plot(ray_exact['x'],ray_exact['Density'], 'g+', ms=2)

pylab.axis([0,500,0,1.1])
pylab.xlabel('Position')
pylab.ylabel('Density')

#error_label(norm['Density'], maxnorm['Density'], 0.45, 1.0)


## Velocity Plot
a = pylab.axes([0.59, 0.57, 0.38, 0.38])
#pylab.axhline(0,color='k',linestyle='dotted')
pylab.plot(exact['x'],exact['x-velocity'],color='r')
pylab.plot(ray['x'],ray['x-velocity'], color='k')
##pylab.plot(ray_exact['x'],ray_exact['x-velocity'], 'g+', ms=2)
#
pylab.axis([0,500,0,500])
pylab.xlabel('Position')
pylab.ylabel('Velocity')
#
#error_label(norm['x-velocity'], maxnorm['x-velocity'], 0.4, 0.3)


# CR Energy Density Plot
#a = pylab.axes([0.59, 0.57, 0.38, 0.38])
#pylab.plot(ray['x'],ray['CREnergyDensity'], 'ro', ms=2)
#
#pylab.axis([0,500,0.0,150])
#pylab.xlabel('Position')
#pylab.ylabel('CR Energy Density')
#
# Pressure Plot
a = pylab.axes([0.09, 0.07, 0.38, 0.38])
pylab.plot((0,250,250,500),(2.0,2.0,0,0),color='k',linestyle='dotted')
pylab.plot(exact['x'],exact['Pressure'],color='r')
pylab.plot(ray['x'],ray['Pressure']/100000, color='k',linestyle='dashed')
pylab.plot(ray['x'],ray['CREnergyDensity']/300000,color='k',linestyle='dashed')
pylab.plot(ray['x'],ray['CREnergyDensity']/300000 + ray['Pressure']/100000,color='k')
#pylab.plot(ray_exact['x'],ray_exact['Pressure'], 'g+', ms=2)

pylab.axis([0,500,0,2.2])
pylab.xlabel('Position')
pylab.ylabel('Pressure/100,000')

#error_label(norm['Pressure'], maxnorm['Pressure'], 0.45, 1.0)


# InternalEnergy Plot
#a = pylab.axes([0.59, 0.07, 0.38, 0.38])
#pylab.axhline(0,color='k',linestyle='dotted')
#pylab.plot(exact['x'],exact['InternalEnergy'])
#pylab.plot(ray['x'],ray['InternalEnergy'], 'ro', ms=2)
#pylab.plot(ray_exact['x'],ray_exact['InternalEnergy'], 'g+', ms=2)

#pylab.axis([0,1,1.5,3.1])
#pylab.xlabel('Position')
#pylab.ylabel('Internal Energy')

#error_label(norm['InternalEnergy'], maxnorm['InternalEnergy'], 0.1, 2.95)


### Save plot
pylab.savefig(png_filename)
