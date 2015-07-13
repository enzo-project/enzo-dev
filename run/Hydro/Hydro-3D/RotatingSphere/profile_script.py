# Makes some radial plots of gas density and temperature.

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from yt.mods import *

import os
import sys

x_field = 'Radiuspc'
y_fields = ['Density', 'Temperature']

weight_field = 'CellMass'

x_min = 1.0e-1
x_max = 2.0e2

fns = sys.argv[1:]

plot_folder = 'profiles'
n_bins = 128

# Make a folder in which to save the profiles.
if not os.path.isdir(plot_folder):
   os.mkdir(plot_folder)

# Make plots for each dataset.
for fn in fns:
   pf = load(fn)

   profile = BinnedProfile1D(pf.h.all_data(), n_bins, x_field, x_min, x_max)
   profile.add_fields(y_fields, weight = weight_field)

   # Make the plot
   for y_field in y_fields:
      fig = plt.figure()
      ax = fig.add_subplot(111)

      ax.loglog(profile[x_field], profile[y_field])

      ax.set_xlabel(x_field)
      ax.set_ylabel(y_field)

      plt.savefig('%s/%s_%s.png' % (plot_folder, pf, y_field));

