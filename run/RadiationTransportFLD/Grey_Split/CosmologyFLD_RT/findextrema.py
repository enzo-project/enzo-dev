# this script calculates the min, max, mean value of the specified field
#
########################################################################

from yt.mods import *
pf = EnzoStaticOutput("DD0273/DD0273")
dd = pf.h.all_data()

# uncomment the wanted field for analysis
#field = "Emissivity"
#field = "Density"
#field = "Gas_Energy"
field = "Temperature"

# calculates for the extrema
extreme = dd.quantities["Extrema"](field)[0]
zmin = extreme[0]
zmax = extreme[1]

# calculates the mean value
mean = dd.quantities["WeightedAverageQuantity"](field,"CellVolume")
print(field)
print("min %16.16e" % zmin)
print("max %16.16e" % zmax)
print("mean %16.16e" % mean)
