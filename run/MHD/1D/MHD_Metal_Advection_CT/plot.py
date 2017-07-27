import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import yt
import tube

ds_list = [yt.load("DD%04d/data%04d"%(frame,frame)) for frame in [0,1]]
tube.tube(ds_list,   legend = True, delta=False, fields=['density','Metal_Density'],filename = 'tube.pdf',labels=['DD%04d'%n for n in [0,1]])
