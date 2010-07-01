# matplotlib-based plotting script for 1D streaming radiation datasets
# Daniel R. Reynolds, reynolds@smu.edu

import h5py
from pylab import *

# load first dataset, and put 1D radiation field into 'Eg'
f = h5py.File('DD0001/pc_amr_0001.cpu0000','r')
Eg = f.get('/Grid00000001/Grey_Radiation_Energy')

# set domain
N = Eg.shape
x = linspace(0.0, 1.0, N[0])

# plot field into figure handle 'h' and close HDF5 file
h = plot(x,Eg,'b-')
xlabel('x')
ylabel('E')
title('Streaming radiation history')
f.close()


# repeat process for next dataset
f = h5py.File('DD0002/pc_amr_0002.cpu0000','r')
Eg = f.get('/Grid00000001/Grey_Radiation_Energy')
plot(x,Eg,'g-')
f.close()

# repeat process for next dataset
f = h5py.File('DD0003/pc_amr_0003.cpu0000','r')
Eg = f.get('/Grid00000001/Grey_Radiation_Energy')
plot(x,Eg,'r-')
f.close()

# repeat process for next dataset
f = h5py.File('DD0004/pc_amr_0004.cpu0000','r')
Eg = f.get('/Grid00000001/Grey_Radiation_Energy')
plot(x,Eg,'c-')
f.close()

# repeat process for next dataset
f = h5py.File('DD0005/pc_amr_0005.cpu0000','r')
Eg = f.get('/Grid00000001/Grey_Radiation_Energy')
plot(x,Eg,'m-')
f.close()

# repeat process for next dataset
f = h5py.File('DD0006/pc_amr_0006.cpu0000','r')
Eg = f.get('/Grid00000001/Grey_Radiation_Energy')
plot(x,Eg,'y-')
f.close()

# repeat process for next dataset
f = h5py.File('DD0007/pc_amr_0007.cpu0000','r')
Eg = f.get('/Grid00000001/Grey_Radiation_Energy')
plot(x,Eg,'k-')
f.close()

# repeat process for next dataset
f = h5py.File('DD0008/pc_amr_0008.cpu0000','r')
Eg = f.get('/Grid00000001/Grey_Radiation_Energy')
plot(x,Eg,'b--')
f.close()

# repeat process for next dataset
f = h5py.File('DD0009/pc_amr_0009.cpu0000','r')
Eg = f.get('/Grid00000001/Grey_Radiation_Energy')
plot(x,Eg,'g--')
f.close()

# repeat process for last dataset
f = h5py.File('DD0010/pc_amr_0010.cpu0000','r')
Eg = f.get('/Grid00000001/Grey_Radiation_Energy')
plot(x,Eg,'r--')
f.close()

# finish off plot and save to file
#legend( ('t1', 't2', 't3', 't4', 't5', 't6', 't7', 't8', 't9', 't10') )
savefig('rad_snapshots.pdf')
