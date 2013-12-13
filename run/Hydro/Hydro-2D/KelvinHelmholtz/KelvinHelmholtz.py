from yt.mods import *

ts = TimeSeriesData.from_filenames("DD*/*.hierarchy")

for pf in ts:
    num = pf.basename[-4:]
    p = SlicePlot(pf,'z',"Density")
    p.set_log("Density", False)
    p.set_zlim("Density", 1,2)
    p.save('plot_%s.png' % num)
