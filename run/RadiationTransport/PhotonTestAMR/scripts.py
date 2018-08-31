import yt
from matplotlib import use; use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from yt.units.yt_array import YTQuantity, YTArray
# Calculate
# 1. Ionization front propagation and hydrodynamic response
# 2. Deviation from 1/r^2 in the photo-ionization field in the last output

first = 1
last = 25
NBINS=16
Myr = 3.1557e13
kpc = 3.086e21

########################################################################

def _MyRadius(field, data):
    center = data.get_field_parameter("center")
    dx = data["x"] - center[0]
    dy = data["y"] - center[1]
    dz = data["z"] - center[2]
    return np.sqrt(dx*dx + dy*dy + dz*dz)
yt.add_field("Radius", function=_MyRadius, take_log=False, units='')

def _MyNeutralFrac(field, data):
    return data['H_p0_fraction'] / 0.75908798
yt.add_field("Neutral_Fraction", function=_MyNeutralFrac, take_log=False,
          units='')
########################################################################

center = [1e-3]*3
time = []
radius = []

for i in range(first, last+1):
    amrfile = "DD%4.4d/data%4.4d" % (i,i)
    pf = yt.load(amrfile)

    x_bins_1d = 32
    r_min = pf.index.get_smallest_dx()
    r_max = pf.quan(1.0 - 1.0/64, 'code_length')
    
    sphere = pf.h.sphere(center, r_max)

    prof1d = yt.create_profile(sphere, 'radius', fields=["Neutral_Fraction", "HI_kph"],
                               n_bins=x_bins_1d,
                               units = {'radius':'code_length'})

    # Find the radius of the I-front (f_HI=0.5)
    res = np.abs(prof1d["Neutral_Fraction"] - 0.5)
    ir = np.where(res == res.min())[0]
    r = np.interp(0.5, prof1d["Neutral_Fraction"], prof1d.x.value)
    r = pf.quan(r, 'code_length')
    time.append(pf.current_time.to('s'))
    radius.append(r.to('cm'))

    del pf
    del prof1d

time = np.array(time)
radius = np.array(radius)

p = plt.subplot(111)
p.plot(time/Myr, radius/kpc, 'k-')
p.set_xlabel("Time (Myr)")
p.set_ylabel(r'$r_{\rm IF}$ (kpc)')
plt.savefig("IFrontRadius.png")

########################################################################
# Radial profiles for some outputs
########################################################################

all_profiles = {}
fields = ["Density", "Temperature", "H_p0_fraction", "H_p1_fraction"]
outputs = [3,5,10,15,25]
x_bins_1d = 20
r_min = 1.0/16
r_max = 1.0 - 1.0/16

for outp in outputs:
    amrfile = "DD%4.4d/data%4.4d" % (outp, outp)
    pf = yt.load(amrfile)
    sphere = pf.h.sphere(center, r_max)
    prof1d = yt.create_profile(sphere, 'radius', fields, extrema=dict(radius=(r_min,r_max)),
                               units = {'radius':'code_length'}, n_bins=NBINS)                                   
    all_profiles[outp] = prof1d
    del pf

for f in fields:
    plt.clf()
    for outp in outputs:
        plt.semilogy(all_profiles[outp].x.value,
                     all_profiles[outp][f], label="%d Myr" % outp)
    plt.xlabel("Radius")
    plt.ylabel(f)
    plt.legend()
    plt.savefig("%sEvo.png" % f)
print("HII = ", all_profiles[25]['H_p1_fraction'])
########################################################################
# Some basic analysis on the final output
########################################################################

pf = yt.load("DD%4.4d/data%4.4d" % (last,last))
MyFields = ['HI_kph', 'H_p0_fraction', 'El_fraction', 'temperature', 'density']
pc = yt.SlicePlot(pf,2,center=[0.5,0.5,1.0/64], fields=MyFields)

pc.save()
del pc

########################################################################
# Calculate deviation from 1/r^2 in the last output (inside 0.5*radius)
########################################################################

x_bins_1d = 20
r_min = 4*pf.index.get_smallest_dx()
r_max = 0.5*radius[-1]
r_max = YTQuantity(r_max, 'cm')
print("rmin, rmax = ", r_min, r_max)
sphere = pf.h.sphere(center, r_max)
fields = ['HI_kph']
prof1d = yt.create_profile(sphere, 'radius', fields,n_bins=NBINS,
                           extrema=dict(radius=(r_min,r_max.to('code_length'))),
                           units = {'radius':'code_length'})
print("np.log(prof1d.x.value[:-1]) = ", np.log(prof1d.x.value[:-1]))
print("np.log(prof1d['HI_kph'][:-1]) = ", np.log(prof1d['HI_kph'][:-1]))
coeff, residual, tr1, tr2, tr3 = \
       np.polyfit(np.log(prof1d.x.value[2:-1]),
                  np.log(prof1d['HI_kph'][2:-1]), 1, full=True)

print("="*72)
print("Inside 2*r_anyl: Radiation field slope = %f +/- %g" % \
      (coeff[0], residual))
print("="*72)
