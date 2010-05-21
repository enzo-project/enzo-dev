from yt.mods import *
pf = EnzoStaticOutput("DD0273/DD0273")
v,c = pf.h.find_max("Emissivity")
print v,c
