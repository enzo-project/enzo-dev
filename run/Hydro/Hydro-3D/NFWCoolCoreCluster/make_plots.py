from yt.mods import *
from yt.config import ytcfg; ytcfg["yt","serialize"] = "False"

def profile(fn,imin=None,imax=None,rmax=None):
   if imin is None:
      imin=0
   if imax is None:
      imax=15
   if rmax is None:
      rmax=10.0
   pf = load(fn)
   pc = PlotCollection(pf)
   pc.add_projection("Density", 0)
   pc.add_projection("Density", 2)
   pc.set_width(4, 'kpc')
#   pc.add_profile_sphere(rmax, "mpc", ["Radiuskpc", "Density"],weight=None)
   pc.add_profile_sphere(rmax, "mpc", ["Radiuskpc", "Temperature"],weight="CellMassMsun")
   pc.save()
   return

if __name__=="__main__":
   import sys
   profile(sys.argv[1])
