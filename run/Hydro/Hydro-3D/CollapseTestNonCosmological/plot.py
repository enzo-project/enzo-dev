import yt.extensions.EnzoSimulation as ES
from matplotlib import pyplot
from yt.mods import *
import os

es = ES.EnzoSimulation("CollapseTestNonCosmological.enzo", initial_time=0)

output_file = 'rho_max.dat'
f = open(output_file,'w')
f.write('#t\trho_max\n')
f.close()

t = []
rho = []

output_dir = "frames"
if not os.path.exists(output_dir): os.mkdir(output_dir)

# Loop over all dataset in the requested time interval.
for output in es.allOutputs:

    pf = load(output['filename'])
    pc = PlotCollection(pf, center=[0.5]*3)
    pc.add_projection('Density', 0)
    pc.save("%s/%s" % (output_dir, pf.basename))

    rho_max = pf.h.find_max('Density')[0]
    t.append(pf.parameters['InitialTime'])
    rho.append(rho_max)

    f = open(output_file, 'a')
    f.write("%f\t%e\n" % (pf.parameters['InitialTime'], rho_max))
    f.close()

axes = pyplot.axes()
axes.semilogy(t, rho)
axes.set_xlabel('t [Code Units]')
axes.set_ylabel('$\\rho_{max}$ [g cm$^{-3}$]')
pyplot.savefig('rho_max.png')
