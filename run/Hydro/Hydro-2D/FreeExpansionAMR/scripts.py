import sys
from matplotlib import use
from yt.mods import *

if len(sys.argv) > 1:
    num = int(sys.argv[1])
else:
    num = 20

amrfile = "DD%4.4d/output_%4.4d" % (num, num)
pf = load(amrfile)
pc = PlotCollection(pf, center=[0.5]*3)

# Make plots of density and gas energy
pc.add_slice('Density', 2)
pc.add_slice('Gas_Energy', 2)
pc.add_slice('VelocityMagnitude', 2)

# Make cuts of density and gas energy
pr1 = pc.add_ray([0.0,0.0,0.5], [1.0,0.0,0.5], "Density")
pr2 = pc.add_ray([0.0,0.0,0.5], [1.0,0.0,0.5], "Gas_Energy")
pr3 = pc.add_ray([0.0,0.0,0.5], [1.0,0.0,0.5], "VelocityMagnitude")

pc.save()
