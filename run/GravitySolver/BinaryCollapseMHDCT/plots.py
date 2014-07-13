from yt.mods import *

for i in xrange(8):
    pf = load("DD%04d/data%04d" % (i,i))
    pc = PlotCollection(pf, center=[0.5]*3)
    for dim in xrange(3):
        pc.add_projection("Density", dim)
    pc.set_width(.3, '1')
    pc.save("%04d_" % i)


