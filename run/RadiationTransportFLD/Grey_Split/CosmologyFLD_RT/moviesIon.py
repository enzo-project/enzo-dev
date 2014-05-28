# Can be modified to create 2D projection of the following fields:
# 1) IonizationFraction
# 2) LogIonizationFraction
# 3) NeutralFraction
# 4) LogNeutralFraction
#
#################################################################

from yt.mods import *

def NeutralFraction(field, data):
   return data["HI_Density"] / data["Density"]

add_field("NeutralFraction", function=NeutralFraction, units=r"\rho_{HI}/\rho")

def LogNeutralFraction(field, data):
   return na.log10(data["NeutralFraction"])

add_field("LogNeutralFraction", function=LogNeutralFraction, units=r"Log \rho_{HI}/\rho")

def IonizationFraction(field, data):
   return data["HII_Density"] / data["Density"]

add_field("IonizationFraction", function=IonizationFraction, units=r"\rho_{HII}/\rho")

def LogIonizationFraction(field, data):
   return na.log10(data["IonizationFraction"])

add_field("LogIonizationFraction", function=LogIonizationFraction, units=r"Log \rho_{HII}/\rho")

min_output_number = 273
max_output_number = 273
skip = 1

res = 64

field = "IonizationFraction"
#field = "LogIonizationFraction"
#field = "NeutralFraction"
#field = "LogNeutralFraction"

frame_template = "frames/movie"+field+"_%04i.png"
basename_template = "DD%04i/DD%04i"

import pylab

for i in range(min_output_number, max_output_number+1, skip):
    pylab.clf()

    basename = basename_template % (i,i)

    pf = load(basename)
    proj = pf.h.proj(0, field, weight_field=field)
    frb = raven.FixedResolutionBuffer(proj, (0.0, 1.0, 0.0, 1.0), (res, res))
    cgsfield = frb[field]

    z = pf["CosmologyCurrentRedshift"]
    yr = pf["InitialTime"] * pf["years"]

#    pylab.imshow( cgsfield, interpolation='nearest', origin='lower') # change to this if using Log of Ion/Neut
    pylab.imshow( na.log10(cgsfield), interpolation='nearest', origin='lower') # change to this if not using Log of Ion/Neut
    pylab.clim(-4,0)
    pylab.title('z = %2.2f, t = %2.2e yr' % (z,yr))
    pylab.colorbar().set_label(r"Log $\rho_{HII}/\rho$") # change label if change field
    pylab.savefig(frame_template % (i), override=True)

