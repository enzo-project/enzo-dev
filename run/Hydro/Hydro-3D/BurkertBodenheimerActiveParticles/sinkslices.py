import matplotlib
matplotlib.use('Agg')
import yt
import numpy as np
import sys
import glob
import matplotlib.pyplot as plt
from yt.units import G, kboltz, mh
from yt.units.yt_array import YTQuantity, YTArray
yt.enable_parallelism()
Mu = 3.0
CENTRE = [0.5, 0.5, 0.5]
#CENTRE = [0.50231934,  0.49035645,  0.49865723]
AccretionRadius = 4
BASE = "./GravPotential/"
fns = glob.glob(BASE + "DD014?/DD????.hierarchy")
fns.sort()
#print fns
WIDTH = 8000

ActiveParticle = "SmartStar"


for f in fns:
    ff = f
    ds = yt.load(ff)
    print("Stats = ", ds.index.get_smallest_dx().in_units("au"))
    dx = ds.index.get_smallest_dx().in_units("au")
    print("dx = ", dx.d)
    dd = ds.all_data()
    centre = CENTRE
    print("Time = ", ds.current_time.in_units("kyr"))
    #Now create a slice of the density field through the sink particle
    #slc = yt.SlicePlot(ds, "x", "density", center='c', width=(0.1, 'pc'))
    #slc.set_zlim('density', 1e-17, 1e-20)
    #slc.save("Slice_x_%s.png" % (ds))
    NumAPs = 0
    flist = dir(ds.fields)
    ap = flist[0]
    try:
        NumAPs = len(dd[ap, "level"])
        print("Number of %s active particles = %d" % (flist[0], NumAPs))
    except:
        print("No active particles found")
    
    #v,c = ds.find_max("velocity_magnitude")
    #print "Max Vel = ", v.in_units("km/s")
    #print "Loc = ", c
    #centre = CENTRE
    #v,c = ds.find_min("velocity_magnitude")
    #print "Min Vel = ", v.in_units("km/s")
    #print "Loc = ", c
    #prj = yt.ProjectionPlot(ds, "x", 'density', center='c',width=(WIDTH, 'au'),
    #                    weight_field=None)
    #prj.set_zlim('density', 1e-1, 1e2)
    #prj.annotate_scale(corner='upper_right')
    #prj.annotate_timestamp(corner='upper_left', redshift=False, draw_inset_box=True)
    
    #prj.save("Proj_x_%s.png" % (ds))

    #slc = yt.SlicePlot(ds, "y", "density", center='c', width=(0.25, 'pc'))
    #slc.set_zlim('density', 1e-17, 1e-20)
    #slc.save("Slice_y_%s.png" % (ds))
    #prj = yt.ProjectionPlot(ds, "y", 'density', center='c',width=(WIDTH, 'au'),
    #                    weight_field=None)
    #prj.set_zlim('density', 1e-1, 1e2)
    #prj.annotate_timestamp(corner='upper_left', redshift=False, draw_inset_box=True)
    #prj.annotate_scale(corner='upper_right')
    #prj.save("Proj_y_%s.png" % (ds))


    slc = yt.SlicePlot(ds, "z", "density", center=centre, width=(WIDTH, 'au'))
    #slc.set_zlim('density', 1e-13, 1e-19)
    if(NumAPs > 0):
        for i in range(NumAPs):
            pos = dd[ActiveParticle, "particle_position"][i]
            slc.annotate_sphere(pos, radius=(dx.d*AccretionRadius, 'au'),
                                circle_args={'color':'black', 'fill':True})
    slc.save("Slice_z_%s.png" % (ds))
    
    prj = yt.ProjectionPlot(ds, "z", 'number_density', center=centre, width=(WIDTH, 'au'),
                            weight_field='density')
    #prj.set_zlim('density', 1e-1, 1e2)
    prj.annotate_timestamp(corner='upper_left', redshift=False, draw_inset_box=True)
    #prj.annotate_scale(corner='upper_right')
    if(NumAPs > 0):
        for i in range(NumAPs):
            pos = dd[ActiveParticle, "particle_position"][i]
            prj.annotate_sphere(pos, radius=(dx.d*AccretionRadius,'au'),
                                circle_args={'color':'black', 'fill':True})
    prj.save("Proj_z_%s.png" % (ds))

    slc = yt.SlicePlot(ds, "z", "velocity_magnitude", center=centre, width=(WIDTH, 'au'))
    slc.set_unit("velocity_magnitude", "km/s")
    #slc.set_zlim('velocity_magnitude', 0.01, 100)
    if(NumAPs > 0):
        for i in range(NumAPs):
            pos = dd[ActiveParticle, "particle_position"][i]
            slc.annotate_sphere(pos, radius=(dx.d*AccretionRadius, 'au'),
                                circle_args={'color':'black', 'fill':True})
    slc.save("Slice_z_Vel_%s.png" % (ds))
    
    prj = yt.ProjectionPlot(ds, "z", 'velocity_magnitude', center=centre,width=(WIDTH, 'au'),
                            weight_field="velocity_magnitude")
    prj.set_unit("velocity_magnitude", "km/s")
    #prj.set_zlim('velocity_magnitude', 0.9, 2)
    prj.annotate_timestamp(corner='upper_left', redshift=False, draw_inset_box=True)
    #prj.annotate_scale(corner='upper_right')
    if(NumAPs > 0):
        for i in range(NumAPs):
            pos = dd[ActiveParticle, "particle_position"][i]
            prj.annotate_sphere(pos, radius = (dx.d*AccretionRadius, 'au'),
                                circle_args={'color':'black', 'fill':True})
    prj.save("Proj_z_Vel_%s.png" % (ds))
    #print dir(plt)
    #
    # Save the image.
    # Optionally, give a string as an argument
    # to name files with a keyword.
    #slc = yt.SlicePlot(ds, "z", "PotentialField", center='c', width=(WIDTH, 'au'))
    #print "Min/Max = ", np.nanmin(dd["PotentialField"]), np.nanmax(dd["PotentialField"])
    #slc.set_log('PotentialField', False)
    #if(NumAPs > 0):
    #    for i in range(NumAPs):
    #        pos = dd["SmartStar", "particle_position"][i]
    #        slc.annotate_sphere(pos, radius=(5, 'au'),
    #                            circle_args={'color':'black', 'fill':True})
    #slc.save("PotentialSlice_z_%s.png" % (ds))

   
    #prj = yt.ProjectionPlot(ds, "z", 'PotentialField', center='c',width=(WIDTH, 'au'),
    #                    weight_field=None)
    #prj.set_zlim('density', 1e-1, 1e2)
    #prj.annotate_timestamp(corner='upper_left', redshift=False, draw_inset_box=True)
    #prj.annotate_scale(corner='upper_right')
    #if(NumAPs > 0):
    #    for i in range(NumAPs):
    #        pos = dd["SmartStar", "particle_position"][i]
    #        prj.annotate_sphere(pos, radius=(5, 'au'),
     #                           circle_args={'color':'black', 'fill':True})
    #prj.save("Proj_z_%s.png" % (ds))
