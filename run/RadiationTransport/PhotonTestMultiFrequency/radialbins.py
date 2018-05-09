from matplotlib import use; use('Agg')
import yt
import matplotlib.pyplot as plt
import numpy as np


first = 0
last = 4
center = [0.5]*3

########################################################################
# Radial profiles for some outputs
########################################################################
YFields = []
Fields = ["Density", "H_fraction", "H_p1_fraction", "H_m1_fraction", "He_p1_fraction",
          "He_p2_fraction",  "HM_kph", "El_fraction", "HeI_kph", 
          "HeII_kph", "HI_kph", "H2II_kdiss", "H2I_kdiss", "H2_fraction"]
FieldsFigure = {}
FieldsAxes = {}
for elem in Fields:
    #Setup Figures
    tmpfigure = plt.figure()
    FieldsFigure[elem] = tmpfigure
    FieldsAxes[elem] = tmpfigure.add_subplot(111)
    YFields.append(elem)

print("Yfields = ", YFields)
outputs = [20]

for outp in outputs:
    amrfile = "DD%4.4d/data%4.4d" % (outp, outp)
    print("Load file %s" % (amrfile))
    ds = yt.load(amrfile)
    sphere = ds.h.sphere(center, (250, 'kpc'))
    rp = yt.create_profile(sphere, 'radius',  YFields, n_bins=8,
                           units = {'radius': 'kpc', 'HM_kph' : '1/s',
                                    'HI_kph' : '1/s','H2I_kdiss' : '1/s',
                                    'H2II_kdiss' : '1/s',
                                    'HeI_kph' : '1/s','HeII_kph' : '1/s'},
                           weight_field='density')
    
    for elem in YFields:
            FieldsAxes[elem].semilogy(rp.x.value, rp[elem], 
                                      label="T = %.2f" % (ds.current_time))
            FieldsAxes[elem].set_xlabel("Radius (kpc)")
            FieldsAxes[elem].set_ylabel(elem)
            FieldsAxes[elem].legend(loc='best')

    amrfile = "OT/DD%4.4d/data%4.4d" % (outp, outp)
    print("Load file %s" % (amrfile))
    ds = yt.load(amrfile)
    sphere = ds.h.sphere(center, (250, 'kpc'))
    rp = yt.create_profile(sphere, 'radius',  YFields, n_bins=8,
                           units = {'radius': 'kpc', 'HM_kph' : '1/s',
                                    'HI_kph' : '1/s','H2I_kdiss' : '1/s',
                                    'H2II_kdiss' : '1/s',
                                    'HeI_kph' : '1/s','HeII_kph' : '1/s'},
                           weight_field='density')
    
    for elem in YFields:
            FieldsAxes[elem].semilogy(rp.x.value, rp[elem], 
                                      label="T = %.2f, OT" % (ds.current_time))
            FieldsAxes[elem].set_xlabel("Radius (kpc)")
            FieldsAxes[elem].set_ylabel(elem)
            FieldsAxes[elem].legend(loc='best')
            
for elem in YFields:
    FieldsFigure[elem].savefig("%s_RadialProfile.png" % (elem))
    
