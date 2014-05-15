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
te = 50

# set the graphics output type
#pictype = 'pdf'
pictype = 'png'

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
    tmp = (data["Total_Energy"] * data["Density"] / 
           (2.0*data["Density"] - data["HI_Density"]))
    return ((gamma - 1.0)*mp/kb * tmp)
add_field("logT", take_log=True, function=_logT, 
          display_name="Temperature", units=r"\rm{K}")

#   Radius from domain center
def _radius(field, data):
    return (np.sqrt(data["x"]*data["x"] + data["y"]*data["y"] +
                    data["z"]*data["z"]))
add_field("radius", take_log=False, function=_radius, 
          display_name="radius", units=r"{r/L_{box}}")



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
    
    # generate 2D plots at certain times
    if (tstep == 1) or (tstep == 5) or (tstep == 10) or (tstep == 20) or (tstep == 50):
        
        # set time label
        if (tstep == 1):
            Myr = '010'
        elif (tstep == 5):
            Myr = '050'
        elif (tstep == 10):
            Myr = '100'
        elif (tstep == 20):
            Myr = '200'
        else:
            Myr = '500'
        
        # determine domain "center" for plotting
        xC = 0.5*(pf["DomainLeftEdge"][0] + pf["DomainRightEdge"][0])
        yC = 0.5*(pf["DomainLeftEdge"][1] + pf["DomainRightEdge"][1])
        zC = 0.5*(pf["DomainLeftEdge"][2] + pf["DomainRightEdge"][2])

        # begin plot collection
        pc = PlotCollection(pf, [xC,yC,0.0])
        
        # xHI slice through z=0
        p = pc.add_slice("xHI",'z')
        p.modify["title"]('HI fraction, t =' + Myr + ' Myr')

        p = pc.add_slice("xHII",'z')
        p.modify["title"]('HII fraction, t =' + Myr + ' Myr')

        p = pc.add_slice("logE",'z')
        p.modify["title"]('radiation density, t =' + Myr + ' Myr')

        p = pc.add_slice("logT",'z')
        p.modify["title"]('temperature, t =' + Myr + ' Myr')

        pc.save(Myr + 'Myr', format=pictype)

        # rename generated files
        f1 = Myr + 'Myr_Slice_z_logE.' + pictype
        f2 = 'Econtour_' + Myr + 'Myr.' + pictype
        rename(f1,f2)
        f1 = Myr + 'Myr_Slice_z_logT.' + pictype
        f2 = 'TempContour_' + Myr + 'Myr.' + pictype
        rename(f1,f2)
        f1 = Myr + 'Myr_Slice_z_xHI.' + pictype
        f2 = 'HIcontour_' + Myr + 'Myr.' + pictype
        rename(f1,f2)
        f1 = Myr + 'Myr_Slice_z_xHII.' + pictype
        f2 = 'HIIcontour_' + Myr + 'Myr.' + pictype
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
        rvals = np.linspace(0,1*pf["cm"]/rs0,nradii)
        HIprofile  = np.zeros(rvals.shape, dtype=float)
        HIIprofile = np.zeros(rvals.shape, dtype=float)
        Tprofile   = np.zeros(rvals.shape, dtype=float)

        # generate 1D profiles from ray emanating out from box center to corner
        for iray in range(0,nrays):
            rvec = rays[iray,:]
            rvec = rvec / sqrt(rvec[0]**2 + rvec[1]**2 + rvec[2]**2)
            r = pf.h.ray([0.0, 0.0, 0.0], rvec)
            HIprof  = log10(r["xHI"])
            HIIprof = log10(r["xHII"])
            Tprof   = log10(r["logT"])
            Hradii  = r["radius"]
        
            # sort results by radius (since that isn't quite working correctly from yt)
            ptype = [('r', float), ('xHI', float), ('xHII', float), ('T', float)]
            pdata = np.zeros(Hradii.shape, dtype=ptype);
            nrad = (Hradii.shape)[0]
            for irad in range(0,nrad):
                pdata[irad] = (Hradii[irad], HIprof[irad], HIIprof[irad], Tprof[irad])
            pdata = np.sort(pdata, order='r')
            
            # interpolate results into output arrays
            tmp = np.interp(rvals, pdata['r'], pdata['xHI'])
            HIprofile += tmp
            tmp = np.interp(rvals, pdata['r'], pdata['xHII'])
            HIIprofile += tmp
            tmp = np.interp(rvals, pdata['r'], pdata['T'])
            Tprofile += tmp
        HIprofile  /= nrays
        HIIprofile /= nrays
        Tprofile   /= nrays
        
        # chemistry profiles
        figure()
        plot(rvals,HIprofile,'b-',rvals,HIIprofile,'r--')
        grid()
        xlabel('$r/L_{box}$')
        ylabel('log(xHI), log(xHII)')
        title('HI, HII Profiles, t =' + Myr + ' Myr')
        legend( ('xHI','xHII'), 'lower right' )
        axis([ 0.0, 1.2, -7.0, 1.0 ])
        savefig('profiles_' + Myr + 'Myr.' + pictype)
        
        # Temperature profile
        figure()
        plot(rvals,Tprofile)
        grid()
        xlabel('$r/r_S$')
        ylabel('log(T) [K]')
        title('Temperature Profile, t =' + Myr + ' Myr')
        axis([ 0.0, 1.2, 3.5, 4.6 ])
        savefig('TempProfile_' + Myr + 'Myr.' + pictype)


# I-front radius
figure()
plot(rdata[0],rdata[1]/rs0,'b-')
xlabel('$t/t_{rec}$')
ylabel('$r_I/r_S$')
title('Propagation of HII Region')
axis([ 0.0, 4.25, 0.0, 1.2 ])
grid()
savefig('rad_vs_time.' + pictype)

