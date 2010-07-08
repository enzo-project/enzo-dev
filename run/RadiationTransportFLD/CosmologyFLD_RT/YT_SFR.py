# This calculates the Star Formation Rate in units of Msun/yt/Mpc^3 (proper) in column 5
# Divide by (1+z)^3 to get it to comoving distances
#
########################################################################################


from yt.mods import *
from yt.extensions.StarAnalysis import *

filename = "DD0273/DD0273"
pf = load(filename, storage_filename = os.path.basename(filename) + ".yt") # only the last DD is required

all = pf.h.all_data()
sfr = StarFormationRate(pf, data_source=all, bins=100) # can add argument to change number of bins i.e. bins=100
sfr.write_out(name="2.0YTStarFormationRate.out") # outputs the file to YTStarFormationRate.out
