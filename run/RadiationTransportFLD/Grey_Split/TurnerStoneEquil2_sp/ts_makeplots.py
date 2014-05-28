# matplotlib-based plotting script for Turner & Stone equilibration tests
# Daniel R. Reynolds, reynolds@smu.edu

# imports
from pylab import *
from scipy.integrate import *

# set the total number of snapshots
nt = 100

# set the graphics output type
pictype = '.png'

# define some helpful functions
def get_params(file):
    """Returns t, dUnit, tUnit, lUnit from a given parameter file"""
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
    return [tval, dUnit, tUnit, lUnit]

def load_vals(tdump):
    """Returns Eg, etot, tval from a given data dump"""
    import h5py
    sdump = repr(tdump).zfill(4)
    pfile = 'DD' + sdump + '/data' + sdump
    hfile = pfile + '.cpu0000'
    tval, dUnit, tUnit, lUnit = get_params(pfile)
    f = h5py.File(hfile,'r')
    Eg3D = f.get('/Grid00000001/Grey_Radiation_Energy')
    et3D = f.get('/Grid00000001/TotalEnergy')
    nx, ny, nz = Eg3D.shape
    Egval = sum(sum(sum(Eg3D)))/ny/nz/nx*dUnit*lUnit*lUnit/tUnit/tUnit
    etval = sum(sum(sum(et3D)))/ny/nz/nx*lUnit*lUnit/tUnit/tUnit
    result = [Egval, etval, tval*tUnit]
    return result


# initialize the data arrays
etot, times = [], []

# loop over snapshots, loading values and times
for it in range(0,nt+1):
    Egval, etotval, timeval = load_vals(it)
    etot.append(etotval)
    times.append(timeval)


# compute the true solution
e0 = etot[0]
def func(e,t):
    """Computes the RHS to the gas energy ODE"""
    import h5py
    c = 2.99792458e10
    k = 4.0e-8
    rho = 1.0e-7
    E = 1.0e12
    sig = 5.6704e-5
    gam = 5.0/3.0
    mu = 0.6
    mp = 1.67262171e-24
    kb = 1.3806504e-16
    dedt = c*k*E/rho - 4*k*sig/rho*((gam-1.0)*mu*mp*e/kb)**4
    return dedt
etrue = odeint(func,e0,times,atol=1.0e-12,rtol=1.0e-8)

# plot results
h = semilogy(times,etot,'r-',times,etrue,'b--')
xlabel('t')
ylabel('e')
title('Turner & Stone Equilibration History')
legend( ('computed', 'true') )
grid()
ax = axis()
axis([ ax[0], ax[1], 10.0**(13.2), 10.0**(15.2) ])
savefig('equil_history' + pictype)

