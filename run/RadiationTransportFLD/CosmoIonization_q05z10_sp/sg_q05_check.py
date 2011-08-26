# numpy/scipy-based error-checking script for Shapiro & Giroux q=0.5 test
# Daniel R. Reynolds, reynolds@smu.edu

# imports
from pylab import *
import numpy as np

# set the total number of snapshots
te = 18

# set the solution tolerance
tol = 0.01

# set some constants
q0 = 0.05              # deceleration parameter
Nph = 5.0e48           # ionization source strength [photons/sec]
alpha2 = 2.52e-13      # recombination rate coefficient
mp = 1.67262171e-24    # proton mass [g]
Myr = 3.15576e13       # duration of a Megayear [sec]

# initialize time-history outputs
#    row 1: i-front radius
#    row 2: stromgren sphere radius (rs)
#    row 3: redshift (z)
#    row 4: time (t)
#    row 5: i-front radius (analytical)
#    row 6: i-front velocity (analytical)
rdata = zeros( (6, te+1), dtype=float);


##########
# define some helpful functions
def get_params(file):
    """Returns z0, z, xR, t0, H0, dUnit, tUnit, lUnit from a parameter files"""
    import shlex
    f = open(file)
    for line in f:
        text = shlex.split(line)
        if ("CosmologyInitialRedshift" in text):
            z0 = float(text[len(text)-1])
        elif ("CosmologyCurrentRedshift" in text):
            z = float(text[len(text)-1])
        elif ("InitialTime" in text):
            t0 = float(text[len(text)-1])
        elif ("CosmologyHubbleConstantNow" in text):
            H0 = float(text[len(text)-1])
    f.close()
    f = open(file + '.rtmodule' )
    for line in f:
        text = shlex.split(line)
        if ("DensityUnits" in text):
            dUnit = float(text[len(text)-1])
        elif ("TimeUnits" in text):
            tUnit = float(text[len(text)-1])
        elif ("LengthUnits" in text):
            lUnit = float(text[len(text)-1])
    f.close()
    xR = lUnit
    t0 *= tUnit
    H0 *= 100*1e5/3.0857e24  # H0 units given 100km/s/Mpc, convert to 1/s 
    return [z0, z, xR, t0, H0, dUnit, tUnit, lUnit]


##########
def analytical_solution(q0,Nph,aval):
    """Analytical solution driver, returns rI, vI"""
    import h5py
    import numpy as np
    import scipy.integrate as sp
    z0, z, xR, t0, H0, dUnit, tUnit, lUnit = get_params('DD0000/data0000')
    f = h5py.File('DD0000/data0000.cpu0000','r')
    rho_data = f.get('/Grid00000001/Density')
    rho = rho_data[0][0][0]*dUnit
    del(rho_data)
    mp = 1.67262171e-24
    # initial nH: no need to scale by a, since a(z0)=1, but we do
    # need to accomodate for Helium in analytical soln
#    nH0 = rho/mp*0.76
    nH0 = rho/mp

    # We first set the parameter lamda = chi_{eff} alpha2 cl n_{H,0} t0, where
    #      chi_{eff} = correction for presence of He atoms [1 -- no correction]
    #      alpha2 = Hydrogen recombination coefficient [2.6e-13 -- case B]
    #      cl = the gas clumping factor [1 -- homogeneous medium]
    #      n_{H,0} = initial Hydrogen number density
    #      t0 = initial time
    alpha2 = 2.52e-13
    lamda = alpha2*nH0*t0
    
    # Compute the initial Stromgren radius, rs0 (proper, CGS units)
    rs0 = (Nph*3.0/4.0/pi/alpha2/nH0/nH0)**(1.0/3.0)  # no rescaling since a(z0)=1
    
    # We have the general formula for y(t):
    #    y(t) = (lamda/xi)exp(-tau(t)) integral_{1}^{a(t)} [da'
    #            exp(t(a'))/sqrt(1-2q0 + 2q0(1+z0)/a')] ,  where
    #    xi = H0*t0*(1+z0),
    #    H0 = Hubble constant
    #    tau(a) = (lamda/xi)*[F(a)-F(1)]/[3(2q0)^2(1+z0)^2/2],
    #    F(a) = [2(1-2q0) - 2q0(1+z0)/a]*sqrt(1-2q0+2q0(1+z0)/a)
    #
    # Here, a' is the variable of integration, not the time-derivative of a.
    F1 = (2.0*(1.0-2.0*q0) - 2.0*q0*(1.0+z0))*sqrt(1.0-2.0*q0+2.0*q0*(1.0+z0))
    xi = H0*t0*(1.0+z0)
    
    # set integration nodes/values (lots)
    inodes = 1000001
    if (aval == 1.0):
        numint = 0.0
    else:
        a = linspace(1,aval,inodes)
        integrand = zeros(inodes, dtype=float)
        arat = divide(2.0*q0*(1.0+z0), a)
        sqa = sqrt(add(1.0-2.0*q0, arat))
        afac = subtract(2*(1-2*q0), arat)
        arg1 = subtract(afac*sqa, F1)
        arg2 = exp(multiply((lamda/xi)/(6*q0*q0*(1+z0)*(1+z0)), arg1))
        integrand = divide(arg2,sqa)
        
        # perform numerical integral via composite Simpson's rule
        numint = sp.simps(integrand, a)
    tauval = (lamda/xi)*((2*(1-2*q0) - 2*q0*(1+z0)/aval)*sqrt(1-2*q0+2*q0*(1+z0)/aval)-F1)/(6*q0*q0*(1+z0)*(1+z0))
    y = lamda/xi*exp(-tauval)*numint;
    
    # extract the current Stromgren radius and velocity
    ythird = sign(y)*abs(y)**(1.0/3.0);
    rI = ythird/aval    # compute ratio rI/rS
    vI = (lamda/3)*aval/ythird*ythird*(1.0-y/aval**3);
    return [rI, vI]



##########
def load_vals(tdump):
    """Returns t, z, xR, nH, Eg, xHI, xHII from a given data dump"""
    import h5py
    import numpy as np
    sdump = repr(tdump).zfill(4)
    pfile = 'DD' + sdump + '/data' + sdump
    hfile = pfile + '.cpu0000'
    z0, z, xR, tval, H0, dUnit, tUnit, lUnit = get_params(pfile)
    f = h5py.File(hfile,'r')
    Eg = f.get('/Grid00000001/Grey_Radiation_Energy')
    HI = f.get('/Grid00000001/HI_Density')
    HII = f.get('/Grid00000001/HII_Density')
    rho = f.get('/Grid00000001/Density')
    # add floor values for happy numerics
    HI  = np.add(HI,  1.0e-10)
    HII = np.add(HII, 1.0e-10)
    Eg  = np.add(Eg,  1.0e-30)
    HIfrac = np.divide(HI,rho)
    HIIfrac = np.divide(HII,rho)
    Eg = np.multiply(Eg,dUnit*lUnit*lUnit/tUnit/tUnit)
    nH = rho[0][0][0]*dUnit/1.67262171e-24
    return [tval, z, xR, nH, Eg, HIfrac, HIIfrac]

##########




# loop over snapshots, loading values and times
for tstep in range(te+1):
    
    # load relevant information
    t, z, xR, nH, Eg, xHI, xHII = load_vals(tstep)

    # compute current Stromgren radius
    rs = (Nph*3.0/4.0/pi/alpha2/nH/nH)**(1.0/3.0)
    
    # store initial hydrogen number density
    if (tstep == 0):
        ti = t
        zi = z
        nHi = nH
        rsi = rs

    # compute volume element
    nx, ny, nz = Eg.shape
    dV = xR*xR*xR/nx/ny/nz
    
    # compute I-front radius (assuming spherical)
    HIIvolume = sum(xHII)*dV*8.0
    rloc = (3.0/4.0*HIIvolume/pi)**(1.0/3.0)
    
    # get analytical solutions for i-front position and velocity
    a = (1.0+zi)/(1.0+z)        # paper's version of a
    ranal, vanal = analytical_solution(q0,Nph,a)
    
    # store data
    rdata[0][tstep] = rloc
    rdata[1][tstep] = rs
    rdata[2][tstep] = z
    rdata[3][tstep] = t
    rdata[4][tstep] = ranal
    rdata[5][tstep] = vanal
    


# I-front radius/velocity plots vs analytical solutions
#   scaled i-front position
r_ratio = rdata[0]/rdata[1]
ranal_ratio = rdata[4]

#   i-front position comparison
r_err = []
for it in range(0, te+1):
    r_err.append( (r_ratio[it]-ranal_ratio[it])/(ranal_ratio[it]+r_ratio[it]+0.1) )

# compute the error norm
err_norm = (np.sum(np.multiply(r_err,r_err))/te)**(0.5)
if (err_norm < tol):
    print 'Error of ',err_norm,' is below tolerance ',tol
    print 'PASS'
else:
    print 'Error of ',err_norm,' is above tolerance ',tol
    print 'FAIL'

