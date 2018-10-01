from matplotlib import use; use('Agg')
import yt
import matplotlib.pyplot as plt
import numpy as np


center = [0.5]*3
Width = 7.0 #in kpc

########################################################################
# Fields to examine
########################################################################
Fields = ["Density", "Temperature", "H_fraction", "H_p1_fraction", "He_p1_fraction", 
          "He_p2_fraction", "HI_kph", "HeI_kph", "HeII_kph"]

outputs = [1,12,23]

for outp in outputs:
    amrfile = "DD%4.4d/data%4.4d" % (outp, outp)
    print("Load file %s" % (amrfile))
    ds = yt.load(amrfile)
    
    for field in Fields:
        slc = yt.SlicePlot(ds, 1, field, center=center, width=(Width, 'kpc'))
        slc.save("%s_slice_%d.png" % (field, outp))
    
