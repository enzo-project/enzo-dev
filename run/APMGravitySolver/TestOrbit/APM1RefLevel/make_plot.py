import numpy as na
import sys
from yt.mods import *
import matplotlib as mpl
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
from math import *
from yt.mods import *
from yt.utilities.linear_interpolators import *

mpl.rcParams['figure.figsize'] = 6,6
mpl.rc('xtick', labelsize=19)
mpl.rc('ytick', labelsize=19)
mpl.rcParams.update({'lines.linewidth':2.0})
mpl.rcParams.update({'legend.numpoints': 1})
mpl.rcParams.update({'legend.fontsize': 19})
mpl.rcParams.update({'legend.frameon': False})

filename = "01.out"
zoom_box = 1
diff_levels = 1
subcycles = 0
solver = "APM"
dataname = "./DD0001/data0001"
size_dots = 2
number_to_skip = 1

# legend
if (diff_levels):
   line1 = "Particles in different levels"
else:
   line1 = "Particles in same level"

if (subcycles):
   line2 = "with subcycling"
else:
   line2 = "no subcycling"

# nb of orbits
pf = load(dataname)
n_orbits = int(pf.current_time)
if (n_orbits > 1):
   orbit_str = " orbits"
else:
   orbit_str = " orbit"

#textinfo1 = line1+"\n"+line2+"\n"+str(n_orbits)+orbit_str
#textinfo2 = solver+"\n"+str(nprocs)+proc_str
textinfo1 = line1+"\n"+line2+"\n"
textinfo2 = "\n"+solver

with open(filename, "r") as FileIn:
   xpos = []
   ypos = []
   xpos2 = []
   ypos2 = []
   for line in FileIn:
      if 'id=1' in line:
         lst = line.split()
         xpos.append(float(lst[1]))
         ypos.append(float(lst[2]))
      if 'id=0' in line:
         lst = line.split()
         xpos2.append(float(lst[1]))
         ypos2.append(float(lst[2]))

x1 = na.array(xpos)
y1 = na.array(ypos)
x2 = na.array(xpos2)
y2 = na.array(ypos2)

# plots orbit in x-y plane (should be a circle with zoom-in inset)
fig = plt.figure()
ax = fig.add_subplot(111, aspect='equal')
ax.scatter(x1[::number_to_skip], y1[::number_to_skip], c='b', s=size_dots, lw=0, marker='.').figure
ax.scatter(x2[::number_to_skip], y2[::number_to_skip], c='r', s=size_dots, lw=0, marker='.').figure
ax.set_xlabel('x',fontsize=21)
ax.set_ylabel('y',fontsize=21)
ax.set_xlim([0.155,0.79])
ax.set_ylim([0.182,0.818])
ax.xaxis.set_major_locator(mpl.ticker.FixedLocator([0.2,0.3,0.4,0.5,0.6,0.7]))
ax.yaxis.set_major_locator(mpl.ticker.FixedLocator([0.2,0.3,0.4,0.5,0.6,0.7]))

# legend
plt.text(0.16,0.735,textinfo1)
plt.text(0.16,0.20,textinfo2)

if (zoom_box):
   xbd = [0.709, 0.719]
   ybd = [0.619, 0.629]

   wh = na.where((xbd[0] < x1) & (xbd[1] > x1) & (ybd[0] < y1) & (ybd[1] > y1))
   iax = zoomed_inset_axes(ax, 10, loc=1)
   iax.scatter(x1[wh], y1[wh], marker='o', s=1, lw=0)
   iax.xaxis.set_ticks([])
   iax.xaxis.set_label('x')
   iax.yaxis.set_ticks([])
   iax.yaxis.set_label('y')
   iax.set_xlim(xbd)
   iax.set_ylim(ybd)
   mark_inset(ax, iax, loc1=2, loc2=4, fc='none', ec='0.5')

# Save
plt.savefig('Testorbit.png',bbox_inches='tight')
plt.close()
