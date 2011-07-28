from matplotlib import use; use('Agg')
from yt.mods import *
import pylab as plt

### define problem name
problem_name = 'BrioWu-MHD-1D'

### define simulation output directory and filename base
output_dir_base = 'DD'
datafile_base = 'data'

### define output to be plotted
dumpid = '0007'

### construct input filename
filename = './' + output_dir_base + dumpid + '/' + datafile_base + dumpid
print "Plotting output file %s\n" % filename

### some more filenames
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
fields = ('Density', 'x-velocity', 'By', 'InternalEnergy')
logme = (False, False, False, True)

count = 0
for f in fields:
    plt.subplot(2,2,count+1)
    if logme[count]:
        plt.semilogy(ray['x'], ray[f])
    else:
        plt.plot(ray['x'], ray[f])
    plt.ylabel(f)
    if count >= 2:
        plt.xlabel('x')
    count += 1
plt.savefig(png_filename)
