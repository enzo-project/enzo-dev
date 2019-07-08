import h5py as h5
import numpy as np
from scipy import *
# This script generates a 3D Gaussian density field and then compute the real and imaginary portions of the wave function, writing it out to files that can be read in by Enzo's FDM Collapse problem.

rank = 3
dim  = 64

attr = {'Component_Rank':1, 'Component_Size':dim**rank, 'Dimensions':[dim]*rank, 'Rank':rank, 'TopGridDims':[dim]*rank, 'TopGridEnd':[dim]*rank, 'TopGridStart':[0]*rank}

# renormalize so mean is 1
dens = 1e-6*np.ones([dim]*rank)
xdim,ydim,zdim = dens.shape

# This should be sped up with meshgrid()
for i in range(xdim):
	for j in range(ydim):
		for k in range(zdim):
			dis = sqrt((i/xdim-0.5)**2 + (j/ydim-0.5)**2 + (k/zdim-0.5)**2)
			#if (dis<0.4):
			dens[i,j,k] += exp(- (dis/0.1)**2)

# write out new density to new file
f1 = h5.File('./GridDensity.new', 'w')
new_dens = f1.create_dataset('GridDensity.new', data=dens)

for a in attr:
	new_dens.attrs.create(a, attr[a])
#for a in dataset.attrs.keys():
#   new_dens.attrs.create(a, dataset.attrs[a])
f1.close()

# unit and constants
BoxLength = 10.0
hbar = 1.05457266e-27
mass_unit = 1.0
mass = 1e-22*1.6021772e-12/(2.99792458e10**2)*mass_unit 
print (mass_unit)

InitialRedshift = 100.
HubbleConstantNow = 0.704
OmegaMatterNow = 0.268

a0 = 1./(1+InitialRedshift)
InitialTime = 0.81651316219217
LengthUnits = 3.085678e24*BoxLength/HubbleConstantNow/(1 + InitialRedshift)
TimeUnits = 2.519445e17/sqrt(OmegaMatterNow)/HubbleConstantNow/(1 + InitialRedshift)**1.5
acoef = 1.5**(1./3.)*a0
coef = hbar/mass*TimeUnits/(LengthUnits**2)
print (coef)

# solve poisson equation for theta
N = dens.shape[-1]
# use fft to solve theta
LHS = -2./3.*(dens-1)/InitialTime/coef

klhs = np.fft.fftn(LHS)

G1d = N * np.fft.fftfreq(N) * 2*pi
kx, ky, kz = meshgrid(G1d, G1d, G1d, indexing='ij')

G2 = kx**2 + ky**2 + kz**2
G2[0,0,0] = 1.

thetak = klhs/(-G2)
thetak[0,0,0] = 0
theta = real(np.fft.ifftn(thetak))

#calculate wave function
repsi = sqrt(dens)*cos(theta)
impsi = sqrt(dens)*sin(theta)
#repsi = sqrt(dens)
#impsi = sqrt(dens)*0
# write out to new file
f1 = h5.File('./GridRePsi', 'w')
new_dens = f1.create_dataset('GridRePsi', data=repsi)
for a in attr:
	new_dens.attrs.create(a, attr[a])
f1.close()

f1 = h5.File('./GridImPsi', 'w')
new_dens = f1.create_dataset('GridImPsi', data=impsi)
for a in attr:
	new_dens.attrs.create(a, attr[a])
f1.close()
