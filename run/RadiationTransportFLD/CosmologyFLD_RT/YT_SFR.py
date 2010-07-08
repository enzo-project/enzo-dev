# 1) Time (yrs)
# 2) Look-back time (yrs)
# 3) Redshift
# 4) Star formation rate in this bin per year (Msol/yr)
# 5) Star formation rate in this bin per year per Mpc**3 (Msol/yr/Mpc**3) (in proper Mpc)
# 6) Stars formed in this time bin (Msol)
# 7) Cumulative stars formed up to this time bin (Msol)
#
# Divide column 5 by (1+z)^3 to get it to comoving distances
#
########################################################################################


from yt.mods import *
from yt.extensions.StarAnalysis import *

filename = "DD0273/DD0273"
pf = load(filename, storage_filename = os.path.basename(filename) + ".yt") # only the last DD is required

all = pf.h.all_data()
sfr = StarFormationRate(pf, data_source=all, bins=100) # can add argument to change number of bins i.e. bins=100
sfr.write_out(name="2.0YTStarFormationRate.out") # outputs the file to YTStarFormationRate.out
