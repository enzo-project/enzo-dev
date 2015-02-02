import yt
for i in range(8):
    ds = yt.load("DD%04d/data%04d" % (i,i))
    for ax in ("x", "y", "z"):
        p = yt.ProjectionPlot(ds, ax, 'density', width=(0.3, 'code_length')).save()
