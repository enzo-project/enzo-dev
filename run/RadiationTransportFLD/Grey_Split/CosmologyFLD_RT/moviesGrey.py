# Creates 2D projection of the Grey_Radiation_Energy field
#
##########################################################

from yt.mods import *

min_output_number = 273
max_output_number = 273
skip = 1

field = "Grey_Radiation_Energy"
res = 64

frame_template = "frames/movie"+field+"_%04i.png"
basename_template = "DD%04i/DD%04i"

import pylab

for i in range(min_output_number, max_output_number+1, skip):
    pylab.clf()

    basename = basename_template % (i,i)

    pf = EnzoStaticOutput(basename)

    convD = pf.conversion_factors["Density"]
    convV = pf.conversion_factors["x-velocity"]
    print(convD)
    print(convV)

    proj = pf.h.proj(0, field, weight_field=field)

    frb = raven.FixedResolutionBuffer(proj, (0.0, 1.0, 0.0, 1.0), (res, res))

    cgsGrey = frb[field]*convD*convV*convV # convert values from code units to CGS

    z = pf["CosmologyCurrentRedshift"]
    yr = pf["InitialTime"] * pf["years"]

    pylab.imshow( na.log10(cgsGrey), interpolation='nearest', origin='lower')
    pylab.clim(-26,-14.4)
    pylab.title('z = %2.2f, t = %2.2e yr' % (z,yr))
    pylab.colorbar().set_label("Log GreyRadEnergy (erg/cm^3)")
    pylab.savefig(frame_template % (i), override=True)

