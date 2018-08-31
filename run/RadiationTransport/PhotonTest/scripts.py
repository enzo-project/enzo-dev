import yt
from matplotlib import use; use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from yt.units.yt_array import YTQuantity, YTArray

# Calculate
# 1. Ionization front propagation and compare to Stroemgren solution
# 2. Slice of photo-ionization field
# 3. Deviation from 1/r^2 in the photo-ionization field in the last output

first = 1
last = 25
Myr = 3.1557e13
kpc = 3.086e21

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

center = [1e-3]*3
time = []
radius = []

for i in range(first, last+1):
    amrfile = "DD%4.4d/data%4.4d" % (i,i)
    pf = yt.load(amrfile)

    x_bins_1d = 32
    r_min = pf.index.get_smallest_dx()
    r_max = pf.arr(1.0 - 1.0/64, 'code_length')
    sphere = pf.h.sphere(center, r_max)
    xvals = np.linspace(r_min, r_max)
    #print("xvals = ", xvals)
    prof1d = yt.create_profile(sphere, 'radius', "Neutral_Fraction", n_bins=32,
                               units = {'radius':'code_length'})

    # Find the radius of the I-front (f_HI=0.5)
    res = np.abs(prof1d["Neutral_Fraction"] - 0.5)
    ir = np.where(res == res.min())[0]
    r = np.interp(0.5, prof1d["Neutral_Fraction"], prof1d.x.value)
    r = pf.quan(r, 'code_length')
    time.append(pf.current_time.to('s'))
    radius.append(r.to('cm'))
    #print("radius = ", r)
    del pf
    del prof1d
#print("time1 = ", time)
time = np.array(time)
radius = np.array(radius)
radius = YTArray(radius, 'cm')
#print("radius = ", radius.to('kpc'))
#print("time = ", time)

# Calculate analytic Stroemgren radius
alpha = 2.59e-13  # Recombination rate at T=1e4 K
nH = 1e-3    # Hydrogen density
lum = 5e48   # Ionizing photon luminosity
trec = 1.0 / (alpha * nH)  # Recombination time
rs = ((3 * lum) / (4*np.pi*alpha*nH**2))**(1.0/3)
anyl_radius = rs * (1.0 - np.exp(-time / trec))**(1.0/3)
#print("anyl_radius = ", anyl_radius)
anyl_radius = YTArray(anyl_radius, 'cm')
#print("anyl_radius = ", anyl_radius.to('kpc'))
ratio = radius / anyl_radius
error = np.abs(1-ratio)
imax = np.where(error == error.max())[0]

p = plt.subplot(211)
p.plot(time/Myr, radius.in_units('kpc'), 'ro')
p.plot(time/Myr, anyl_radius.to('kpc'), 'k-')
p.set_ylabel(r'$r_{\rm IF}$ (kpc)')
p = plt.subplot(212)
p.plot(time/Myr, ratio, 'k-')
p.set_ylabel(r'$r_{\rm IF} / r_{\rm anyl}$')
p.set_xlabel("Time (Myr)")
plt.savefig("StroemgrenRadius.png")

########################################################################
# Some basic analysis on the final output
########################################################################

pf = yt.load("DD%4.4d/data%4.4d" % (last,last))
myfields = ["HI_kph", "Neutral_Fraction", "El_fraction"]
pc = yt.SlicePlot(pf,2,center=[0.5,0.5,1.0/64], fields=myfields)
pc.save()
del pc

########################################################################
# Calculate deviation from 1/r^2 in the last output (inside 0.5*anyl_radius)
########################################################################

x_bins_1d = 32
r_min = 2*pf.index.get_smallest_dx()
r_max = YTQuantity(0.5*anyl_radius[-1], 'cm')
#print("anyl_radius = ", anyl_radius)
anyl_radius = pf.arr(anyl_radius, 'cm')
#print("anyl_radius = ", anyl_radius.to('code_length'))
#print("r_max = ", r_max)
#print("r_min = ", r_min)
r_max = pf.quan(r_max, 'cm')
r_max = r_max.in_units('code_length')
#print("r_max = ", r_max)
r_max = 100*r_min
sphere = pf.h.sphere(center, r_max)

prof1d = yt.create_profile(sphere, "radius", fields=["HI_kph"])
coeff, residual, tr1, tr2, tr3 = \
       np.polyfit(np.log(prof1d.x.value[:-1]),
                  np.log(prof1d['HI_kph'][:-1]), 1, full=True)

print("="*72)
print("Maximum error in ionization front radius = %g (at %f Myr)" % \
      (error.max(), time[imax]/Myr))
print("Average error in ionization front radius = %g" % (error.mean()))
print("Inside 2*r_anyl: Radiation field slope = %f +/- %g" % \
      (coeff[0], residual))
print("="*72)
