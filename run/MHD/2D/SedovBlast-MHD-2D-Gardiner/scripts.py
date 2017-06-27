from matplotlib import use; use('Agg')
from yt.mods import *
import pylab as plt

### define problem name
problem_name = 'SedovBlast-MHD-2D-Gardiner'

### define simulation output directory and filename base
output_dir_base = 'DD'
datafile_base = 'data'

### define output to be plotted
dumpid = '0009'

### construct input filename
filename = './' + output_dir_base + dumpid + '/' + datafile_base + dumpid
print("Plotting output file %s\n" % filename)

### some more filenames
png_filename = './' + problem_name

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



### define fields vector
fields = ('Density', 'Pressure', 'Bx', 'By')
pc = PlotCollection(pf, center=[0.0,0.0,0.0])

for f in fields:
    pc.add_slice(f, 2)

pc.save(png_filename)
