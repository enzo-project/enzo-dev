# matplotlib-based plotting script for Iliev et al. test #2
# Daniel R. Reynolds, reynolds@smu.edu

# imports
from pylab import *
from yt.mods import *

# reduce log level of yt 
from yt.config import ytcfg
ytcfg["yt","loglevel"] = "50"

# set the total number of snapshots
te = 50

# set the solution tolerance
tol = 0.01

# load the reference solution
r_sol = np.load('r_sol.npy')

# set some constants
Ngammadot = 5.0e48     # ionization source strength [photons/sec]
aHII = 2.52e-13        # recombination rate coefficient
mp = 1.67262171e-24    # proton mass [g]
Myr = 3.15576e13       # duration of a Megayear [sec]
nH = 1.0e-3            # input hydrogen number density [cm^(-3)]
trec = 1.0/(aHII*nH)   # recombination time [sec]
rs0 = (3.0*Ngammadot/4/pi/aHII/nH/nH)**(1.0/3.0)   # Stromgren radius

# Define some derived fields
#   Neutral Hydrogen fraction (log plot)
def _xHI(field, data):
    return (data["HI_Density"]/data["Density"])
add_field("xHI", take_log=True, function=_xHI, 
          display_name="Neutral\; Fraction")

#   Ionized Hydrogen fraction (log plot)
def _xHII(field, data):
    return (data["HII_Density"]/data["Density"])
add_field("xHII", take_log=True, function=_xHII, 
          display_name="Ionized\; Fraction")

#   Radiation energy density (log plot)
def _logE(field, data):
    return (data["Grey_Radiation_Energy"])
def _convertlogE(data):
    return ( data.convert("MassUnits")/data.convert("TimeUnits") / 
             data.convert("TimeUnits")/data.convert("cm"))
add_field("logE", take_log=True, function=_logE, 
          convert_function=_convertlogE, 
          display_name="Radiation\; Energy\; Density", 
          units=r"\rm{erg}/\rm{cm}^3")

#   Temperature (log plot)
def _logT(field, data):
    mp = 1.67262171e-24
    kb = 1.3806504e-16
    gamma = 5.0/3.0
    tmp = (data["TotalEnergy"] * data["Density"] / 
           (2.0*data["Density"] - data["HI_Density"]))
    return ((gamma - 1.0)*mp/kb * tmp)
add_field("logT", take_log=True, function=_logT, 
          display_name="Temperature", units=r"\rm{K}")

#   Radius from domain center
def _radius(field, data):
    return (np.sqrt(data["x"]*data["x"] + data["y"]*data["y"] +
                    data["z"]*data["z"]))
def _convertradius(data):
    return (data.convert("cm"))
add_field("radius", take_log=False, function=_radius, 
          convert_function=_convertradius, 
          display_name="radius", units=r"\rm{cm}")



# initialize time-history outputs
#    row 1: time (t)
#    row 2: computed i-front radius
rdata = zeros( (2, te+1), dtype=float);


# loop over snapshots, loading values and times
for tstep in range(0,te+1):
    
    # load relevant information
    sdump = repr(tstep).zfill(4)
    pfile = 'DD' + sdump + '/data' + sdump
    pf = load(pfile)
    t = pf.current_time * pf["TimeUnits"]

    # determine if simulation was run with source in center or corner
    spherical = (2.0**(pf.domain_left_edge[0]/pf.domain_right_edge[0]+1.0) 
               * 2.0**(pf.domain_left_edge[1]/pf.domain_right_edge[1]+1.0) 
               * 2.0**(pf.domain_left_edge[2]/pf.domain_right_edge[2]+1.0))

    # compute I-front radius (assuming spherical)
    sp = pf.h.sphere([0.0, 0.0, 0.0], 1.0)
    HIIvolume = (sp["xHII"]*sp["CellVolumeCode"]*pf["cm"]**3).sum()*spherical
    radius = (3.0/4.0*HIIvolume/pi)**(1.0/3.0)
    
    # store data
    rdata[0][tstep] = t/trec
    rdata[1][tstep] = radius
    

# compute I-front radius comparison, error norm
r_err = (rdata[1] - r_sol)/rs0
r_err_norm = (np.sum(np.multiply(r_err,r_err))/te)**(0.5)
if (r_err_norm < tol):
    print 'Error of ',r_err_norm,' is below tolerance ',tol
    print 'PASS'
else:
    print 'Error of ',r_err_norm,' is above tolerance ',tol
    print 'FAIL'

