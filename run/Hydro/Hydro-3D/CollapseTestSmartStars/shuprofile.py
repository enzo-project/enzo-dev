import matplotlib
matplotlib.use('Agg')
import yt
import numpy as np
import sys
import matplotlib.pyplot as plt
from yt.units import G, kboltz, mh
from yt.units.yt_array import YTQuantity, YTArray
import glob
from scipy.signal import savgol_filter
yt.enable_parallelism()

mycolors = iter(['red', 'blue', 'black', 'green', 'magenta', 'grey'])

def _Accretion_Rate(field, data):  #M' = 4*pi*R*R*rho(R)*V(R)
    return -4.0*np.pi*data["radius"]*data["radius"]*data["density"]*data["radial_velocity"]

#Add Derived Rates
yt.add_field(("gas", "shell_accretion_rate"), units="Msun/yr",
             function=_Accretion_Rate, sampling_type="cell")

BASE = "./"
Mu = 3.0
CENTRE = [0.50, 0.50, 0.50]


Files = glob.glob(BASE + "DD00??/DD00??")
Files.sort()
Files = ["DD0000/DD0000",
         "DD0008/DD0008",
         "DD0012/DD0012",
         "DD0015/DD0015",
         "DD0017/DD0017",
         "DD0021/DD0021"]
print("Files = ", Files)


Fields = ["density", "temperature", "radial_velocity", "cell_mass",
          "shell_accretion_rate"]

YFields = []
FieldsFigure = {}
FieldsAxes = {}
for elem in Fields:
    #Setup Figures
    tmpfigure = plt.figure()
    FieldsFigure[elem] = tmpfigure
    FieldsAxes[elem] = tmpfigure.add_subplot(111)
    YFields.append(elem)

accrate = []
ages = []
cellwidth = 0.0
for f in Files:
    ff = f
    ds = yt.load(ff)
    dd = ds.all_data()
    NumAPs = 0
    flist = dir(ds.fields)
    if yt.is_root():
        print("CellWidth = %e cm" % (ds.index.get_smallest_dx().in_units("cm")))
        print("4*CellWidth = %e cm" % (4*ds.index.get_smallest_dx().in_units("cm")))
        print("CellWidth = %e pc" % (ds.index.get_smallest_dx().in_units("pc")))
        
        print("Current Time = %e s" % (ds.current_time.in_units("s")))
        print("Current Time = %e code" % (ds.current_time.in_units("code_time")))
        print("Current Time = %e yr" % (ds.current_time.in_units("yr")))
    try:
        NumAPs = len(dd["SmartStar", "level"])
        centre = dd["SmartStar", "particle_position"][0]
    except:
        centre = CENTRE
    cellwidth = ds.index.get_smallest_dx().in_units("cm")
    time = ds.current_time.in_units("yr")
    if(time > 0 and NumAPs > 0):
        timeindex = dd["SmartStar", "TimeIndex"][0].d
        accrate = dd["SmartStar", "AccretionRate"][0][1:int(timeindex+1)].in_units("Msun/yr")
        ages = dd["SmartStar", "AccretionRateTime"][0][1:int(timeindex+1)].in_units("yr")
        deltaT = np.ediff1d(ages, to_begin=1)
        masses = accrate*deltaT
        print("Pos = %f %f %f" % (dd["SmartStar", "particle_position"][0][0],
                                  dd["SmartStar", "particle_position"][0][1],
                                  dd["SmartStar", "particle_position"][0][2]))
        #print("ages = ", ages)
        #print("masses = ", masses)
    if(time > 0 and NumAPs == 0):
        continue
    print("NumAPs = ", NumAPs)
    print("Accreting Particle Position = ", centre)
  
    # Create a sphere of radius 1 kpc in the center of the box.
    my_sphere = ds.sphere(centre, (1.0, "pc"))
    color = next(mycolors)
    # Create a profile of the average density vs. radius.
    prof = yt.create_profile(my_sphere, "radius", YFields,
                             units = {'radius': 'cm'}, 
                             weight_field="cell_mass")

    Radius = prof.x[prof["density"] > 0.0]
    Density = prof["density"][prof["density"] > 0.0]
    Temperature = prof["temperature"][prof["density"] > 0.0]
    
    #Plot Density
    FieldsAxes["density"].loglog(Radius, Density, label="T = %2.2e yrs" %
                                 (ds.current_time.in_units("yr")), linewidth=2.0, color=color)
    
    #Add in analytic expression from Shu et al. (1977)
    T = ds.current_time.in_units("s")
    cs = np.sqrt(kboltz*Temperature/(Mu*mh))
    #cs = 0.166 * 1e5
    m0 = np.sqrt(133.0)
    csT = (cs*T).d
    #print("cs = ", cs/1e5)
    
    shu_density = np.power(cs, 1.5)*m0/(4*np.pi*G*np.sqrt(2.0*T)*np.power(Radius, 1.5))
    shu_density = shu_density[Radius.d < 2e16]
    shu_radius = Radius[Radius.d < 2e16]
    #print("shu_density = ", shu_density)
    FieldsAxes["density"].loglog(shu_radius, shu_density, ls='dotted', linewidth=2.0, color=color)
    #Radial Velocity Plots
    RadialVelocity = prof["radial_velocity"][prof["density"] > 0.0].in_units("km/s")
    EnclosedMass = prof["cell_mass"][prof["density"] > 0.0]
    
    FieldsAxes["radial_velocity"].semilogx(Radius, RadialVelocity,
                                           label="T = %2.2e yrs" %
                                           (ds.current_time.in_units("yr")),
                                           linewidth=2.0, color=color)
    #Add in analytic expression from Shu et al. (1977)
    # Need reasonably high resolution to get here.
    shu_radial = -cs*m0*np.sqrt(2*cs*T/Radius)/1e5
    shu_radial = shu_radial[Radius.d < 2e16]
    FieldsAxes["radial_velocity"].semilogx(shu_radius, shu_radial, ls='dotted', linewidth=2.0,
                                           color=color)
    #Accretion Plots
    AccretionRate = prof["shell_accretion_rate"][prof["density"] > 0.0]
    #print "Accretion Rate = ", AccretionRate
    if(np.sum(AccretionRate) > 0.0):
        FieldsAxes["shell_accretion_rate"].loglog(Radius, AccretionRate,
                                                  label="T = %2.2e yrs" %
                                                  (ds.current_time.in_units("yr"))
                                                  ,linewidth=2.0, color=color)
        #Add in analytic expression
        shu_accretion_rate = (m0*m0*np.power(cs, 3.0)/G).d
        shu_accretion_rate = YTArray(shu_accretion_rate,  'g/s', registry=ds.unit_registry)
        #print "Shu Accretion Rate = ",  shu_accretion_rate
        shu_accretion_rate = shu_accretion_rate[Radius.d < 2e16].in_units("Msun/yr")
        #print "Shu Accretion Rate = ",  shu_accretion_rate
        FieldsAxes["shell_accretion_rate"].loglog(shu_radius, shu_accretion_rate,
                                                  ls='dotted', linewidth=2.0, color=color)


    #ProjectionPlots
    prj = yt.ProjectionPlot(ds, "x", "number_density", weight_field="density", width=(0.025, 'pc'), center=centre,
                            axes_unit="cm")
    prj.set_zlim("number_density", 1e4, 1e8)
    prj.annotate_sphere(centre, radius=(0.0001, 'pc'),
                  circle_args={'color':'black', 'fill':'true'})
    prj.save("NumberDensity_Proj_%s.png" % (ds))


XS = cellwidth.d/2
XE = 1e17
#Add in accretion radius
#Add in minimum cell width
YCELLWIDTH = [1e-30, 1e-10]
XCELLWIDTH = np.zeros_like(YCELLWIDTH)
XCELLWIDTH[0] = ds.index.get_smallest_dx().in_units("cm")
XCELLWIDTH[1] = ds.index.get_smallest_dx().in_units("cm")
#print("XCELLWIDTH[0] = ", XCELLWIDTH[0])
FieldsAxes["density"].loglog(XCELLWIDTH, YCELLWIDTH, ls='dashed', label="Minimum Cellwidth", linewidth=2.0)
FieldsAxes["density"].loglog(XCELLWIDTH*4.0, YCELLWIDTH, ls='dashed', label="Accretion Radius", linewidth=2.0)
FieldsAxes["density"].legend()
FieldsAxes["density"].set_xlim(XS, XE)
FieldsAxes["density"].set_ylim(1e-19, 1e-11)
FieldsAxes["density"].set_xlabel("R [cm]", fontsize=18)
FieldsAxes["density"].set_ylabel("$\\rho$ [g cm$^{-3}$]", fontsize=18)
FieldsAxes["density"].xaxis.labelpad = -2
FieldsFigure["density"].savefig("DensityProfile.png")
FieldsFigure["density"].savefig("DensityProfile.pdf")

YCELLWIDTH = [-30, 0]
FieldsAxes["radial_velocity"].semilogx(XCELLWIDTH, YCELLWIDTH, ls='dashed', label="Minimum Cellwidth", linewidth=2.0)
FieldsAxes["radial_velocity"].semilogx(XCELLWIDTH*4.0, YCELLWIDTH, ls='dashed', label="Accretion Radius", linewidth=2.0)
FieldsAxes["radial_velocity"].legend(loc='best')
FieldsAxes["radial_velocity"].set_xlabel("R [cm]", fontsize=18)
FieldsAxes["radial_velocity"].set_ylabel("V$_r$ [km s$^{-1}$]", fontsize=18)
FieldsAxes["radial_velocity"].set_xlim(XS, XE)
FieldsAxes["radial_velocity"].set_ylim(-15, 0)
FieldsAxes["radial_velocity"].xaxis.labelpad = -2
FieldsFigure["radial_velocity"].savefig("RadialVelocityProfile.png")
FieldsFigure["radial_velocity"].savefig("RadialVelocityProfile.pdf")

YCELLWIDTH = [1e-10, 1e2]
FieldsAxes["shell_accretion_rate"].loglog(XCELLWIDTH, YCELLWIDTH, ls='dashed', label="Minimum Cellwidth", linewidth=2.0)
FieldsAxes["shell_accretion_rate"].loglog(XCELLWIDTH*4.0, YCELLWIDTH, ls='dashed', label="Accretion Radius", linewidth=2.0)
FieldsAxes["shell_accretion_rate"].legend()
FieldsAxes["shell_accretion_rate"].set_xlim(XS, XE)
FieldsAxes["shell_accretion_rate"].set_ylim(1e-6, 1e-2)
FieldsAxes["shell_accretion_rate"].set_xlabel("R [cm]", fontsize=18)
FieldsAxes["shell_accretion_rate"].set_ylabel("Accretion Rate [M$_{\odot}$ yr$^{-1}$]", fontsize=18)
FieldsAxes["shell_accretion_rate"].xaxis.labelpad = -2
FieldsFigure["shell_accretion_rate"].savefig("AccretionRateProfile.png")
FieldsFigure["shell_accretion_rate"].savefig("AccretionRateProfile.pdf")




plt.figure()

accrate =  savgol_filter(accrate, 33, 1,mode='nearest')
#plot the accretion rate over time
plt.plot(ages, accrate*1e4, lw=2.0, color='k')
plt.xlabel("Time [yr]", fontsize=18)
plt.ylabel("Accretion Rate [10$^4$ M$_{\odot}$ yr$^{-1}$]", fontsize=18)
#plt.ylim(1e-5, 1e-3)
plt.savefig("AccretionRateEvolution.png")
#Average Accretion Rate between 5000 and 12000 years
val_low = ages[ages > 5000][0]
index_low = np.asscalar(np.where(ages.d == val_low)[0])
val_high = ages[ages > 12000][0]
index_high = np.asscalar(np.where(ages.d == val_high)[0])

avg_array = accrate[index_low:index_high]
print("Average Accretion Rate = %e Msolar/yr" % (np.mean(avg_array)))

plt.figure()
masses = np.cumsum(masses)
masses =  savgol_filter(masses, 5, 1,mode='nearest')
#plot the mass against time
plt.plot(ages, masses, lw=2.0, color='k')
plt.xlabel("Time [yr]", fontsize=18)
plt.ylabel("Particle Mass [M$_{\odot}$]", fontsize=18)
plt.savefig("MassEvolution.png")

