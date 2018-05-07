import yt
from matplotlib import use; use('Agg')
import matplotlib.pyplot as plt
import numpy as na
import glob
# Calculate
# 1. Time evolution of the average temperature and ionization fraction
#    of the clump
# 2. Line cut of temperature through the clump center
# 3. Neutral fraction and temperature slices at the last output

first = 1
last = 15

Myr = 3.1557e13
kpc = 3.086e21

########################################################################
def _MyNeutralFrac(field, data):
    return data['H_p0_fraction'] / 0.75908798
yt.add_field("Neutral_Fraction", function=_MyNeutralFrac, take_log=True,
          units='')
########################################################################

center = [0.75757575, 0.5, 0.5]
radius = 0.121212
time = []
xe = []
temp = []

for i in range(first, last+1):
    amrfile = "DD%4.4d/data%4.4d" % (i,i)
    pf = yt.load(amrfile)

    sphere = pf.h.sphere(center, radius)

    avg_xe = 1.0 - sphere.quantities.weighted_average_quantity("Neutral_Fraction", "cell_volume")
    avg_temp = sphere.quantities.weighted_average_quantity("Temperature", "cell_volume")
    print("Time = ", pf.current_time.to('Myr'))
    time.append(pf.current_time.to('Myr'))
    xe.append(avg_xe)
    temp.append(avg_temp)

    del pf

time = na.array(time)
xe = na.array(xe)
temp = na.array(temp)

p = plt.subplot(211)
p.plot(time, xe, 'k-')
p.set_ylabel(r'$\bar{x}_e$ (clump)')
p = plt.subplot(212)
p.plot(time, temp, 'k-')
p.set_ylabel(r'$\bar{T}$ (clump)')
p.set_xlabel("Time (Myr)")
plt.savefig("ClumpEvo.png")

########################################################################
# Line cuts through the clump center
########################################################################
all_profiles = {}
fields = ["Neutral_Fraction", "Temperature"]
outputs = [1,3,5,10,15]

for outp in outputs:
    amrfile = "DD%4.4d/data%4.4d" % (outp, outp)
    pf = yt.load(amrfile)
    ray = pf.h.ortho_ray(0, (0.5,0.5))
    all_profiles[outp] = {}
    all_profiles[outp]['x'] = ray['x']
    for f in fields:
        all_profiles[outp][f] = ray[f]
    del pf

for f in fields:
    plt.clf()
    for outp in outputs:
        plt.semilogy(all_profiles[outp]['x'],
                     all_profiles[outp][f], label="%d Myr" % outp)
    plt.xlabel("x")
    plt.ylabel(f)
    plt.legend(loc='best')
    plt.savefig("%sLineEvo.png" % f)

########################################################################
# Some basic analysis on the final output
########################################################################

BASE = "./"
f = BASE + ("DD00%d/data00%d.hierarchy" % (last, last))

myfields = ["HI_kph", "Neutral_Fraction", "Temperature"]
ds = yt.load(f)
slc = yt.SlicePlot(ds, 2, center=[0.5,0.5,0.5], fields=myfields)
slc.save()

