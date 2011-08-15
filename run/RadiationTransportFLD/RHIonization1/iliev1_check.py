# matplotlib-based plotting script for Iliev et al. test #1
# Daniel R. Reynolds, reynolds@smu.edu

# imports
from pylab import *
import numpy as np


# set the total number of snapshots
te = 50

# set the solution tolerance
tol = 0.004

# set some constants
Ngammadot = 5.0e48     # ionization source strength [photons/sec]
aHII = 2.52e-13        # recombination rate coefficient
mp = 1.67262171e-24    # proton mass [g]
Myr = 3.15576e13       # duration of a Megayear [sec]
nH = 1.0e-3            # input hydrogen number density [cm^(-3)]
trec = 1.0/(aHII*nH)   # recombination time [sec]
rs0 = (3.0*Ngammadot/4/pi/aHII/nH/nH)**(1.0/3.0)   # Stromgren radius

# initialize time-history outputs
#    row 1: time (t)
#    row 2: computed i-front radius
#    row 3: predicted i-front radius (rI)
#    row 4: stromgren sphere radius (rs)
rdata = zeros( (4, te+1), dtype=float);


# define some helpful functions
def get_params(file):
    """Returns t, vol, dUnit, tUnit, lUnit from a given parameter file"""
    import shlex
    f = open(file)
    for line in f:
        text = shlex.split(line)
        if ("InitialTime" in text):
            tval = float(text[len(text)-1])
        elif ("DensityUnits" in text):
            dUnit = float(text[len(text)-1])
        elif ("TimeUnits" in text):
            tUnit = float(text[len(text)-1])
        elif ("LengthUnits" in text):
            lUnit = float(text[len(text)-1])
        elif ("Gamma" in text):
            gamma = float(text[len(text)-1])
        elif ("DomainLeftEdge" in text):
            xL = float(text[len(text)-3])
            yL = float(text[len(text)-2])
            zL = float(text[len(text)-1])
        elif ("DomainRightEdge" in text):
            xR = float(text[len(text)-3])
            yR = float(text[len(text)-2])
            zR = float(text[len(text)-1])
    vol = (xR-xL)*(yR-yL)*(zR-zL)*lUnit*lUnit*lUnit
    tval = tval*tUnit
    return [tval, vol, dUnit, tUnit, lUnit]

def load_vals(tdump):
    """Returns t, vol, Eg, HIfrac, HIIfrac from a given data dump"""
    import h5py
    import numpy as np
    sdump = repr(tdump).zfill(4)
    pfile = 'DD' + sdump + '/data' + sdump
    hfile = pfile + '.cpu0000'
    tval, vol, dUnit, tUnit, lUnit = get_params(pfile)
    f = h5py.File(hfile,'r')
    Eg = f.get('/Grid00000001/Grey_Radiation_Energy')
    HI = f.get('/Grid00000001/HI_Density')
    HII = f.get('/Grid00000001/HII_Density')
    rho = f.get('/Grid00000001/Density')
    HIfrac = np.divide(HI,rho)
    HIIfrac = np.divide(HII,rho)
    Eg = np.add(Eg, 1.0e-30)
    Eg = np.multiply(Eg,dUnit*lUnit*lUnit/tUnit/tUnit)
    return [tval, vol, Eg, HIfrac, HIIfrac]



# loop over snapshots, loading values and times
for tstep in range(te+1):
    
    # load relevant information
    t, vol, Eg, xHI, xHII = load_vals(tstep)
    
    # compute volume element
    nx, ny, nz = Eg.shape
    dV = vol/nx/ny/nz
    
    # compute I-front radius (assuming spherical)
    HIIvolume = sum(xHII)*dV*8.0
    radius = (3.0/4.0*HIIvolume/pi)**(1.0/3.0)
    
    # compute analytical solution
    ranal = rs0*(1.0 - exp(-t/trec))**(1.0/3.0)
    
    # store data
    rdata[0][tstep] = t/trec
    rdata[1][tstep] = radius
    rdata[2][tstep] = ranal
    rdata[3][tstep] = rs0
    


# I-front radius comparison (skip left-most point)
r_err = (rdata[1][1:] - rdata[2][1:])/rs0
r_err_norm = (np.sum(np.multiply(r_err,r_err))/te)**(0.5)
if (r_err_norm < tol):
    print 'Error of ',r_err_norm,' is below tolerance ',tol
    print 'PASS'
else:
    print 'Error of ',r_err_norm,' is above tolerance ',tol
    print 'FAIL'
