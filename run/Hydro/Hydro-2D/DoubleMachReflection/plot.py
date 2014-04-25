from yt.mods import *

ts = TimeSeriesData.from_filenames("DD*/*.hierarchy")

for pf in ts:
    num = pf.basename[-4:]

    c = [2.0,0.5,0]

    slice = SlicePlot(pf, 'z', "Density", center=c,
                      width=((4.15, '1'), (1.15, '1')),
                      origin='domain', axes_unit=('1','1'))

    slice.set_zlim("Density",0.5, 20.0)

    slice.set_buff_size((1000, 300))
    
    slice.set_window_size(12)
    
    slice.save(mpl_kwargs={'bbox_inches':'tight'},name="DoubleMachTest_%s.png" % num)
