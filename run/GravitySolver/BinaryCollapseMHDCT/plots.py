from yt.mods import *

for i in range(8):
    ds = load("DD%04d/data%04d" % (i,i))
    for ax in 'xyz':
        p = ProjectionPlot(ds, ax, 'Density', width=0.3).save()
