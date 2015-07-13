# Makes some projections through the gas.
from yt.mods import *
import sys
import os

# Plotting parameters
plot_folder = 'plots'

field = 'Density'
axis = 'z'

# Widths of plots in pc
widths = [100.0]

# Set up a folder for the plots
if not os.path.isdir(plot_folder):
   os.mkdir(plot_folder)

# Load each dataset and make a plot at each zoom level
for fn in sys.argv[1:]:
   pf = load(fn)
   
   for width in widths:
      prj = ProjectionPlot(pf, axis, field, weight_field='Density', width=width/pf['pc'])
      prj.annotate_velocity()
      prj.save('%s/%s_%f_pc.png' % (plot_folder, pf, width))
