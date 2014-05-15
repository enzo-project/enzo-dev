# This script generates a random field with three components, and is meant
# to approximate a turbulent velocity field.
#
# The field is not divergence free, and thus is only suitable for compressable
# gas. For a more complete discussion of generating divergence free turbulent
# initial conditions, see Rogallo 1981 (1981STIN...8131508R).
#
# For standard Kolmogorov turbulence, the turbulent spectrum will have an index
# of -5/3. For compressible gas, however, the index will be steeper.

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import itertools

n_velocity_components = 3
grid_dim = 128
n = -3.0

# Min and max wave number in terms of the box size
k_min  = 2.0
k_max = float(grid_dim)

# Figure out the grid dimensions
grid_shape = []

for i in range(n_velocity_components):
   grid_shape.append(grid_dim)


# Initialize the RNG
np.random.seed()

# Generate the components in fourier space
turbulence_field = []

for component in range(n_velocity_components):
   F = np.zeros(grid_shape, dtype=complex)

   for index in itertools.product(range(grid_dim/2), repeat=n_velocity_components):
      k = np.sqrt(sum([a*b for a,b in zip(index, index)]))
   
      if k >= k_min and k <= k_max:
         amplitude = np.sqrt(k**n/(4.0 * np.pi * k * k))
   
         for region in itertools.product([-1, 1], repeat = n_velocity_components):
            new_index = tuple([a*b for a,b in zip(region, index)])
            phase = -1.0j*(2.0 * np.pi * np.random.random())
            F[new_index] = amplitude * np.exp(phase)

   f = np.fft.ifftn(F)
   turbulence_field.append(np.real(f))

# Write the perturbations to a file.
outf = open('turbulence.in', 'w')

# Write out a header.
outf.write('# Turbulent velocity field\n')
outf.write('# P(k) ~ k^%i\n' % n)
outf.write('# k_min: %d\n' % k_min)
outf.write('# k_max: %d\n' % k_max)
outf.write('# Dimensions: %i\n' % n_velocity_components)
outf.write('# Grid size:')

for component in range(n_velocity_components):
   outf.write(' %i' % grid_dim)

outf.write('\n')

# First non commented line is the number of dimensions
outf.write('%i\n' % n_velocity_components)

# Second non commented line is the grid size
for component in range(n_velocity_components):
   outf.write('%i ' % grid_dim)
outf.write('\n')

# Write out each index
for index in itertools.product(range(grid_dim), repeat=n_velocity_components):
   for comp in index:
      outf.write('%i ' % comp)
   for i in range(n_velocity_components - 1):
      outf.write('%e ' % turbulence_field[i][index])

   outf.write('%e\n' % turbulence_field[n_velocity_components-1][index])

outf.close()

# Make a plot of the field
for component in range(n_velocity_components):
   fig = plt.figure()
   ax = fig.add_subplot(111)

   if n_velocity_components == 2:
      ax.imshow(turbulence_field[component], interpolation='nearest')

   if n_velocity_components == 3:
      ax.imshow(turbulence_field[component][0], interpolation='nearest')
   
   plt.savefig('turbulence_%i.png' % component)
