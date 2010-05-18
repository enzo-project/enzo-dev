"""
inits_sort.py -- A script to move particles to the correct grid for nested cosmology
initial conditions.

Purpose: When making nested initial conditions with 'inits.exe', depending on
the cosmology, some of the particles created for nested subgrids are given
positions outside that nested grid. Before running 'ring.exe', these particles
need to be moved to the correct grid. This script is therefore run between
'inits.exe' and 'ring.exe', and re-creates a corrected set of initial condition
files.

Method: ParticleMass files *must* be created when running inits for this script
to work. After the initial conditions are created, copy this script to the same
directory as the files. Edit this script (below) to match your settings, and run:
'python sort.py'. The new initial conditions are placed in a directory
named 'new_ICs'.

Requires: Python with Numpy and h5py.

Notes: Obviously, this only works for smaller datsets that fit into memory.

Stephen Skory
sskory@physics.ucsd.edu
December, 2009
"""

import h5py, sys, os,time
import numpy as na

# Put the files here, in list format, top grid first.
PP = ['ParticlePositions.0', 'ParticlePositions.1']
PV = ['ParticleVelocities.0', 'ParticleVelocities.1']
PM = ['ParticleMasses.0', 'ParticleMasses.1']

# Put the boundaries here in array format. Enter *only* the outer edges of
# each grid.
boundaries = []
# Top grid, which almost certainly the same as this.
boundaries.append((na.array([0.]*3), na.array([1.]*3)))
# Nested grid 1.
boundaries.append((na.array([0.375]*3), na.array([0.625]*3)))
# Nested grid 2.
# boundaries.append((na.array([0.45]*3), na.array([0.55]*3)))

# Grid refinement ratio. This is typically 2.
refinement = 2

#### END USER INPUT ####

tstart = time.time()

# Open the files
print 'Reading data...'
PPfp, PVfp, PMfp = [], [], []
for i,file in enumerate(PP):
    PPfp.append(h5py.File(file))
    PVfp.append(h5py.File(PV[i]))
    PMfp.append(h5py.File(PM[i]))

# Read in the data arrays.
PParr = []
PVarr = []
PMarr = []
for i,fp in enumerate(PPfp):
    PParr.append(fp[PP[i]][:,:])
    PVarr.append(PVfp[i][PV[i]][:,:])
    PMarr.append(PMfp[i][PM[i]][:])

tread = time.time()

# Put the particles into one big array, and then decide where to put them.
size = 0
for arr in PParr:
    size += arr.shape[1]
PPall = na.empty((3,size), dtype='float64')
PVall = na.empty((3,size), dtype='float64')
PMall = na.empty((1,size), dtype='float64')
# Tracks where they came from.
source = na.zeros(size, dtype='int8')
count = 0
for i,arr in enumerate(PParr):
    PPall[:,count:arr.shape[1]+count] = arr[:,:]
    PVall[:,count:arr.shape[1]+count] = PVarr[i][:,:]
    PMall[0,count:arr.shape[1]+count] = PMarr[i][:]
    source[count:arr.shape[1]+count] += i
    count += arr.shape[1]

# We're done with the original data.
print 'Calculating new locations...'
del PParr, PVarr, PMarr

# Make a mapping of where the particles should go.
dest = na.zeros(size, dtype='int8')
done = na.zeros(size, dtype='bool')
for i, b in enumerate(boundaries):
    # We go backwards, smallest to largest grid
    j = len(boundaries) - i - 1
    LE, RE = boundaries[j]
    is_inside = ( (PPall.T >= LE).all(axis=1) * \
                (PPall.T < RE).all(axis=1) )
    # We only change ones that haven't been assigned before.
    is_inside = is_inside * na.invert(done)
    dest += j * is_inside
    # Keep track of which particles have already been assigned.
    done = na.bitwise_or(is_inside, done)

# Clean up.
del is_inside, done

# Find the sizes of the final files using histogram.
bins = na.arange(len(boundaries) + 1)
sizes, bins = na.histogram(dest, bins)

# For the mass field, we need to fix the masses for particles that change
# assignments. A negative value of diff means that the particle went to a
# less refined grid, so it should drop in mass. And visa versa for positive
# diffs.
refinement = float(refinement)**3
diff = dest - source
multi = na.power(refinement, diff)
PMall *= multi

tcalc = time.time()

# With the places they go, make new files similar to the old ones, and
# save the data there.
print 'Writing new data to disk in the directory "new_ICs"...'
try:
    os.mkdir('new_ICs')
except OSError:
    pass
os.chdir('new_ICs')
for i, fp in enumerate(PPfp):
    # Choose the particles to save for this file.
    choose = (dest == i)
    this_size = sizes[i]
    # Positions
    new_PPfp = h5py.File(PP[i], 'w')
    PPds = new_PPfp.create_dataset(PP[i], data = PPall[:,choose])
    # Copy over the attributes, but making appropriate changes.
    for attr in PPfp[i][PP[i]].attrs.iteritems():
        if attr[0] == 'Component_Size' or attr[0] == 'Dimensions':
            PPds.attrs[attr[0]] = na.array([this_size])
        else:
            PPds.attrs[attr[0]] = attr[1]
    # save and close the files that we're done with.
    new_PPfp.flush()
    new_PPfp.close()
    PPfp[i].close()

    # Velocities.
    new_PVfp = h5py.File(PV[i], 'w')
    PVds = new_PVfp.create_dataset(PV[i], data = PVall[:,choose])
    for attr in PVfp[i][PV[i]].attrs.iteritems():
        if attr[0] == 'Component_Size' or attr[0] == 'Dimensions':
            PVds.attrs[attr[0]] = na.array([this_size])
        else:
            PVds.attrs[attr[0]] = attr[1]
    # save and close
    new_PVfp.flush()
    new_PVfp.close()
    PVfp[i].close()

    # Masses.
    new_PMfp = h5py.File(PM[i], 'w')
    PMds = new_PMfp.create_dataset(PM[i], data = PMall[:,choose])
    for attr in PMfp[i][PM[i]].attrs.iteritems():
        if attr[0] == 'Component_Size' or attr[0] == 'Dimensions':
            PMds.attrs[attr[0]] = na.array([this_size])
        else:
            PMds.attrs[attr[0]] = attr[1]
    # save and close
    new_PMfp.flush()
    new_PMfp.close()
    PMfp[i].close()

tdone = time.time()

# Figure out how many particles moved.
tosum = na.zeros(size, dtype='int8')
diff = (diff != 0)
tosum = tosum[diff]

# Figure out the timings.
full = tdone - tstart
IO = tread - tstart + tdone - tcalc
calc = tcalc - tread

print 'Moved %d particles to new grids in %f seconds.' % (tosum.size, full)
print 'Spent %f sec in disk IO and %f sec in calculations.' % (IO, calc)
