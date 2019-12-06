# This script generates a density field with a white noise spectrum for k=0-20
# and amplitude 1%.
#

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import itertools

n_dimension = 3
grid_dim = 128
min_index = -int(grid_dim/2)
max_index =  int(grid_dim/2)
# Min and max wave number in terms of the box size
k_min = 1
k_max = 20 

# Figure out the grid dimensions
grid_shape = []

for i in range(n_dimension):
   grid_shape.append(grid_dim)


# Initialize the RNG
np.random.seed()

# Generate the components in fourier space
F = np.zeros(grid_shape, dtype=complex)

index_range = range(min_index, max_index)
for x in index_range:
   for y in index_range:
      for z in index_range:
         index = (x, y, z)
         k = np.linalg.norm(index)

         if k >= k_min and k <= k_max:
            amplitude = np.random.rand()
            phase =  -1.0j*(2.0 * np.pi * np.random.random())
            
            F[index] = amplitude*np.exp(phase)

f = np.fft.ifftn(F)
density_field = np.real(f)


# Write the perturbations to a file.
outf = open('white_noise.in', 'w')

# Write out a header.
outf.write('# white noise density field\n')
outf.write('# k_min: %d\n' % k_min)
outf.write('# k_max: %d\n' % k_max)
outf.write('# Dimensions: %i\n' % n_dimension)
outf.write('# Grid size:')

for component in range(n_dimension):
   outf.write(' %i' % grid_dim)

outf.write('\n')

# First non commented line is the dimensions
outf.write('%i\n' % n_dimension)

# Second non commented line is the grid size
for component in range(n_dimension):
   outf.write('%i ' % grid_dim)
outf.write('\n')

# Write out each index
for index in itertools.product(range(grid_dim), repeat=n_dimension):
   for comp in index:
      outf.write('%i ' % comp)
   outf.write('%e\n' % density_field[index])

outf.close()

# Make a plot of the field
fig = plt.figure()
ax = fig.add_subplot(111)
ax.imshow(density_field[0], interpolation='nearest')
plt.savefig('white_noise.png')


