# matplotlib-based plotting script for Iliev et al. test #2
# Daniel R. Reynolds, reynolds@smu.edu

# imports
from pylab import *
from yt.mods import *
from os import *

# reduce log level of yt 
from yt.config import ytcfg
ytcfg["yt","loglevel"] = "50"

# set the total number of snapshots
te = 19

# set the graphics output type
#pictype = 'pdf'
pictype = 'png'

# set some constants
q0 = 0.5               # deceleration parameter
Nph = 5.0e48           # ionization source strength [photons/sec]
alpha2 = 2.52e-13      # recombination rate coefficient
mp = 1.67262171e-24    # proton mass [g]


##########
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
add_field("logE", take_log=True, function=_logE, 
          display_name="Radiation\; Energy\; Density")

#   Radius from domain center
def _radius(field, data):
    return (np.sqrt(data["x"]*data["x"] + data["y"]*data["y"] +
                    data["z"]*data["z"]))
def _convertradius(data):
    return (data.convert("cm"))
add_field("radius", take_log=False, function=_radius, 
          convert_function=_convertradius, 
          display_name="radius", units=r"\rm{cm}")

##########

# define some helpful functions
def CompositeSimpson(fvals,xvals):
    """Performs the composite simpson numerical integral given a set of evenly-spaced nodes and function values (requires odd number of each)"""
    import numpy as np
    # determine width of each interval
    h = xvals[2] - xvals[0]
    # get total number of nodes
    n = fvals.size
    # initialize at end points
    integral = fvals[0] + fvals[n-1]
    # accumulate even-indexed interior contributions
    integral += 2.0*np.sum(fvals[2:n-1:2])
    # accumulate odd-indexed contributions
    integral += 4.0*np.sum(fvals[1:n:2])
    integral *= h/6.0
    return integral

def CompositeTrapezoid(fvals,xvals):
    """Performs the composite trapezoidal numerical integral given a set of evenly-spaced nodes and function values"""
    import numpy as np
    # determine width of each interval
    h = xvals[1] - xvals[0]
    # get total number of nodes
    n = fvals.size
    # accumulate even-indexed contributions
    integral = 2.0*h*np.sum(fvals) - h*fvals[0] - h*fvals[n-1]
    integral /= 2.0
    return integral

def analytical_solution(q0,Nph,pf0,pf):
    """Analytical solution driver, returns rI, vI"""
    import numpy as np
    #import scipy.integrate as sp
    z0 = pf["CosmologyInitialRedshift"]
    z  = pf["CosmologyCurrentRedshift"]
    H0 = pf0["CosmologyHubbleConstantNow"]*100*1e5/3.0857e24  # convert to CGS
    sph = pf0.h.sphere([0.0, 0.0, 0.0], 0.1)
    rho = sph["Density"].max()
    aval = (1.0+z0)/(1.0+z)     # paper's version of a
    mp = 1.67262171e-24
    # initial nH: no need to scale by a, since a(z0)=1, but we do
    # need to accomodate for Helium in analytical soln
    nH0 = rho/mp

    # set Enzo's TimeUnits (since yt doesn't have access to this)
    tUnit = 2.519445e17 / ( np.sqrt(pf0["CosmologyOmegaMatterNow"]) 
                          * pf0["CosmologyHubbleConstantNow"]
                          * (1.0 + pf0["CosmologyInitialRedshift"])**1.5 )
    # set current time [CGS]
    t0 = pf0["InitialTime"]*tUnit
    
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
    a = linspace(1,aval,inodes)  # uniformly spaced nodes
    integrand = zeros(inodes, dtype=float)
    arat = divide(2.0*q0*(1.0+z0), a)
    sqa = sqrt(add(1.0-2.0*q0, arat))
    afac = subtract(2*(1-2*q0), arat)
    arg1 = subtract(afac*sqa, F1)
    arg2 = exp(multiply((lamda/xi)/(6*q0*q0*(1+z0)*(1+z0)), arg1))
    integrand = divide(arg2,sqa)
    
    # perform numerical integral via composite Simpson's rule
    numint = CompositeSimpson(integrand, a)
    tauval = (lamda/xi)*((2*(1-2*q0) - 2*q0*(1+z0)/aval)*sqrt(1-2*q0+2*q0*(1+z0)/aval)-F1)/(6*q0*q0*(1+z0)*(1+z0))
    y = lamda/xi*exp(-tauval)*numint;
    
    # extract the current Stromgren radius and velocity
    ythird = sign(y)*abs(y)**(1.0/3.0);
    rI = ythird/aval    # compute ratio rI/rS
    vI = (lamda/3)*aval/ythird*ythird*(1.0-y/aval**3);
    return [rI, vI]

##########




# initialize time-history outputs
#    row 1: i-front radius
#    row 2: stromgren sphere radius (rs)
#    row 3: redshift (z)
#    row 4: time (t)
#    row 5: i-front radius (analytical)
#    row 6: i-front velocity (analytical)
rdata = zeros( (6, te+1), dtype=float);


# loop over snapshots, loading values and times
for tstep in range(0,te+1):
    
    # load relevant information
    sdump = repr(tstep).zfill(4)
    pfile = 'DD' + sdump + '/data' + sdump
    pf = load(pfile)
    xR = pf["DomainRightEdge"][0]*pf["cm"]
    z  = pf["CosmologyCurrentRedshift"]
    t  = pf["InitialTime"]
    sph = pf.h.sphere([0.0, 0.0, 0.0], 0.1)
    nH = sph["Density"].max()/mp

    # compute current Stromgren radius
    rs = (Nph*3.0/4.0/pi/alpha2/nH/nH)**(1.0/3.0)

    # store initial values of some constants
    if (tstep == 0):
        pf0 = pf
        ti = t
        zi = z
        nHi = nH
        rsi = rs

    # determine domain "center" for plotting
    xC = 0.5*(pf["DomainLeftEdge"][0] + pf["DomainRightEdge"][0])
    yC = 0.5*(pf["DomainLeftEdge"][1] + pf["DomainRightEdge"][1])
    zC = 0.5*(pf["DomainLeftEdge"][2] + pf["DomainRightEdge"][2])

    # determine if simulation was run with source in center or corner
    spherical = (2.0**(pf.domain_left_edge[0]/pf.domain_right_edge[0]+1.0) 
               * 2.0**(pf.domain_left_edge[1]/pf.domain_right_edge[1]+1.0) 
               * 2.0**(pf.domain_left_edge[2]/pf.domain_right_edge[2]+1.0))

    # compute I-front radius (assuming spherical)
    sp = pf.h.sphere([0.0, 0.0, 0.0], 1.0)
    HIIvolume = (sp["xHII"]*sp["CellVolumeCode"]*pf["cm"]**3).sum()*spherical
    rloc = (3.0/4.0*HIIvolume/pi)**(1.0/3.0)

    # get analytical solutions for i-front position and velocity
    ranal, vanal = analytical_solution(q0,Nph,pf0,pf)
    
    # store data
    rdata[0][tstep] = rloc
    rdata[1][tstep] = rs
    rdata[2][tstep] = z
    rdata[3][tstep] = t
    rdata[4][tstep] = ranal
    rdata[5][tstep] = vanal
    
    # generate 2D plots at certain times
    if (tstep == 1) or (tstep == 5) or (tstep == 10):
        
        # set time label
        tout = repr(tstep).zfill(2)
        
        # begin plot collection (center at (xC,yC,0))
        pc = PlotCollection(pf, [xC,yC,0.0])
        
        # xHI slice through z=0
        p = pc.add_slice("xHI",'z')
        p.modify["title"]('HI fraction, tstep =' + tout)

        p = pc.add_slice("xHII",'z')
        p.modify["title"]('HII fraction, tstep =' + tout)

        p = pc.add_slice("logE",'z')
        p.modify["title"]('radiation density, tstep =' + tout)

        pc.save('tstep' + tout, format=pictype)

        # rename generated files
        f1 = 'tstep' + tout + '_Slice_z_logE.' + pictype
        f2 = 'Econtour_' + tout + '.' + pictype
        rename(f1,f2)
        f1 = 'tstep' + tout + '_Slice_z_xHI.' + pictype
        f2 = 'HIcontour_' + tout + '.' + pictype
        rename(f1,f2)
        f1 = 'tstep' + tout + '_Slice_z_xHII.' + pictype
        f2 = 'HIIcontour_' + tout + '.' + pictype
        rename(f1,f2)

        # generate profile plots by averaging results from multiple rays
#         rays = np.array( ((1.0,0.0,0.0), (0.0,1.0,0.0), (0.0,0.0,1.0), 
#                           (-1.0,0.0,0.0), (0.0,-1.0,0.0), (0.0,0.0,-1.0), 
#                           (1.0,1.0,0.0), (1.0,0.0,1.0), (0.0,1.0,1.0), 
#                           (-1.0,1.0,0.0), (-1.0,0.0,1.0), (0.0,-1.0,1.0), 
#                           (1.0,-1.0,0.0), (1.0,0.0,-1.0), (0.0,1.0,-1.0), 
#                           (-1.0,-1.0,0.0), (-1.0,0.0,-1.0), (0.0,-1.0,-1.0), 
#                           (1.0,1.0,1.0), (-1.0,1.0,1.0), (1.0,-1.0,1.0), 
#                           (1.0,1.0,-1.0), (-1.0,-1.0,1.0), (-1.0,1.0,-1.0), 
#                           (1.0,-1.0,-1.0), (-1.0,-1.0,-1.0)) )
#         nrays = 26
        rays = np.array( ((1.0,0.0,0.0), (0.0,1.0,0.0), (0.0,0.0,1.0), 
                          (1.0,1.0,0.0), (1.0,0.0,1.0), (0.0,1.0,1.0), 
                          (1.0,1.0,1.0)) )
        nrays = 7
        nradii = 200
        rvals = np.linspace(0.0,1.0*pf["cm"]/rsi,nradii)
        HIprofile  = np.zeros(rvals.shape, dtype=float)
        HIIprofile = np.zeros(rvals.shape, dtype=float)

        # generate 1D profiles from ray emanating out from box center to corner
        for iray in range(0,nrays):
            rvec = rays[iray,:]
            rvec = rvec / sqrt(rvec[0]**2 + rvec[1]**2 + rvec[2]**2)
            r = pf.h.ray([0.0, 0.0, 0.0], rvec)
            HIprof  = log10(r["xHI"])
            HIIprof = log10(r["xHII"])
            Hradii  = r["radius"]/xR
        
            # sort results by radius (since that isn't quite working correctly from yt)
            ptype = [('r', float), ('xHI', float), ('xHII', float)]
            pdata = np.zeros(Hradii.shape, dtype=ptype);
            nrad = (Hradii.shape)[0]
            for irad in range(0,nrad):
                pdata[irad] = (Hradii[irad], HIprof[irad], HIIprof[irad])
            pdata = np.sort(pdata, order='r')
            
            # interpolate results into output arrays
            tmp = np.interp(rvals, pdata['r'], pdata['xHI'])
            HIprofile += tmp
            tmp = np.interp(rvals, pdata['r'], pdata['xHII'])
            HIIprofile += tmp
        HIprofile  /= nrays
        HIIprofile /= nrays
        
        # chemistry profiles
        figure()
        plot(rvals,HIprofile,'b-',rvals,HIIprofile,'r--')
        grid()
        xlabel('$r/L_{box}$')
        ylabel('log(xHI), log(xHII)')
        title('HI, HII Profiles, z = ' + repr(z))
        legend( ('xHI','xHII'), 'lower right' )
        axis([ 0.0, 1.2, -7.0, 1.0 ])
        savefig('profiles_' + tout + '.' + pictype)
        

# I-front radius plot vs analytical solution

#   scaled i-front position
r_ratio = rdata[0]/rdata[1]
ranal_ratio = rdata[4]

#   scaled redshift (cell centsteprs)
z_ratio = (1.0 + rdata[2])/(1.0+zi)

#   i-front position vs redshift plot
figure()
xdata = -log10(z_ratio)
plot(xdata,r_ratio,'b-',xdata,ranal_ratio,'r--')
xlabel('$-log[(1+z)/(1+z_i)]$')
ylabel('$r_I/r_S$')
title('r_i(t)/r_s(t) vs redshift, q_0 =' + repr(q0))
legend( ('computed', 'analytical'), loc=4 )
grid()
axis([ 0.0, 3.0, 0.0, 1.0 ])
savefig('radius.' + pictype)
