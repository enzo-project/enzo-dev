# SCRIPT TO MAKE INITIAL CONDITION FILES FOR COSMOLOGICAL SIMULATIONS
#
# Uses mpgrafic and degraf (if necessary)
#
# Written by John H. Wise (June 2008)
#-----------------------------------------------------------------------

import os, os.path, sys
from math import *

#-------------------- PARAMETERS --------------------
#enzo         = True              # On for Enzo ICs
nprocs       = 16                 # Number of processors
boxsize      = 17.831669               # Boxsize in comoving Mpc (not Mpc/h)
resolution   = 512                # Topgrid resolution
n_levels     = 0                 # Number of nested grids
inner_width  = 0.4375            # If using nested grids, width of finest grid
buffer_cells = 4                 # Number of cells between refined grids
seed         = 200905130         # Random seed (MUST be 9-digits)
name         = "LAE512"

center = [0.5, 0.5, 0.5]         # Center of interest (no shift = 0.5,0.5,0.5)
#center = [0.7717, 0.0405, 0.3435]   # Center of interest (no shift = 0.5,0.5,0.5)
LargeScaleCorrection = False     # If using noise from a low-resolution run
LargeScaleFile = "200902187_64.dat" # That noise file
OneDimPerFile = True             # Write one dimension per file

#-------------------- COSMOLOGY PARAMETERS --------------------
omega_m     = 0.279              # Omega matter
omega_v     = 0.721              # Omega lambda
omega_b     = 0.0462             # Omega baryon
h0          = 70.1               # Hubble constant
sigma8      = 0.817              # sigma_8
n_plawslope = 0.960              # Slope of power spectrum

#-----------------------------------------------------------------------

def write_grafic1inc(res):
    contents = """
!c  grafic1.inc
!c  User must set parameters then make grafic1
!c
!c  Parameters for offsetting baryon and CDM velocity fields by
!c  0.5 grid spacing.
!c  The origin of the baryon density field is (0,0,0).  The origin
!c  of the baryon velocity field will be offvelb*(0.5,0.5,0.5)*dx
!c  and similarly for offvelc.  offvelb/c should be +/-1 or 0.
	integer offvelb,offvelc
	parameter (offvelb=0,offvelc=0)
!c
!c  Initial density fluctuation amplitude on the ultimate subgrid scale.
!c  For large grids (> 10^7) one may wish to decrease this; for small
!c  refinement grids (< 10^6) and multiple levels one may wish to increase it.
	real sigstart
	parameter (sigstart=0.2)
!c
!c  Top grid mesh size.
!c
	integer np1,np2,np3
	parameter (np1=%d,np2=%d,np3=%d)
!c
!c np[1-3]:      particle lattice size for computing initial conditions.
!c                 All must be divisible by 2.  They need not be equal.
!c  The memory required by grafic1 is np1*np2*np3 words plus a few.
    """ % (res, res, res)
    open("grafic1.inc", "w").write(contents)

def run_mpgrafic():
    NoiseFile = "%9.9d_%d.dat" % (seed, resolution << n_levels)
    contents = """
4
%f %f %f
%f
%f
-%f
%f %f
-%g
1
1




1
%d
%s
%d
%s
    """ % \
    (omega_m, omega_v, h0, omega_b, n_plawslope, sigma8, 0.01, 100.0, \
     boxsize*(h0/100.0), seed, NoiseFile, LargeScaleCorrection, LargeScaleFile)
    open("params.in", "w").write(contents)
    cmdline = "mpirun -np %d ./mpgrafic < params.in" % (nprocs)
    os.system(cmdline)

def run_degraf():
    max_degrade = 2**n_levels
    if OneDimPerFile:
        _files = "'GridDensity' 'GridVelocities_x' 'GridVelocities_y' 'GridVelocities_z' 'ParticleVelocities_x' 'ParticleVelocities_y' 'ParticleVelocities_z'"
        _nfiles = 7
    else:
        _files = "'GridDensity' 'GridVelocities' 'ParticleVelocities'"
        _nfiles = 3
    for dim in range(3):
        center[dim] = fmod(center[dim]-0.5, 1.0)
        
    contents = """
&parameters
nfiles = %d
inputfiles = %s
degramax = %d
recurse = .true.
shift = %f %f %f
width_in = %f
buffercells = %d
/
    """ % \
    (_nfiles, _files, max_degrade, center[0], center[1], center[2],
     inner_width, buffer_cells)
    open("params.nml", "w").write(contents)
    cmdline = "mpirun -np %d ./degraf params.nml" % (nprocs)
    os.system(cmdline)

#------------------------------------------------------------------------

if len(sys.argv) == 2:
    name = sys.argv[1]
    if name.isdigit and len(name) == 9:
        seed = int(name)
    else:
        print ""
        print "usage: python %s [seed]" % sys.argv[0]
        print "--> seed is optional and must be 9 digits"
        print "--> seed = %s" % name
        sys.exit()

finest_res = resolution * 2**n_levels
dir_name = "data/%s" % (name)
NoiseFile = "%9.9d_%d.dat" % (seed, resolution << n_levels)

if not os.path.exists("data"):
    os.mkdir("data")
if not os.path.exists(dir_name):
    os.mkdir(dir_name)

# Run mpgrafic
os.chdir("mpgrafic-0.2/src")
if os.path.exists("mpgrafic"):
    os.remove("mpgrafic")
write_grafic1inc(finest_res)
os.chdir("..")
os.system("make clean")
os.system("make")
if not os.path.exists("src/mpgrafic"):
    print "mpgrafic does not exist.  Failure in make?  Must run ./configure first."
    os.chdir("..")
    sys.exit(1)
os.chdir("src")
run_mpgrafic()

# Move files
if OneDimPerFile:
    files = ['GridDensity', 'GridVelocities_x', 'GridVelocities_y',
             'GridVelocities_z', 'ParticleVelocities_x', 'ParticleVelocities_y',
             'ParticleVelocities_z']
else:
    files = ['GridDensity', 'GridVelocities', 'ParticleVelocities']
if (n_levels == 0):
    for f in files:
        os.rename("%s" % f, "../../%s/%s" % (dir_name, f))
else:
    os.remove(NoiseFile)
    for f in files:
        os.rename("%s" % f, "../../degraf/src/%s" % f)
os.chdir("../..")

# If necessary, degrade with degraf
if (n_levels > 0):
    os.chdir("degraf")
    if os.path.exists("src/degraf"):
        os.remove("src/degraf")
#    os.system("make clean")
    os.system("make")
    if not os.path.exists("src/degraf"):
        print "degraf does not exist.  Failure in make?  ", \
              "Must run ./configure first."
        os.chdir("..")
        sys.exit(1)
    os.chdir("src")
    run_degraf()

    for f in files:
        os.remove("%s" % (f))
    for i in range(n_levels+1):
        for f in files:
            os.rename("%s.%d" % (f, i),
                      "../../%s/%s.%d" % (dir_name, f, i))
    os.rename("enzo.params", "../../%s/enzo.params" % (dir_name))

