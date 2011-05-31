# numpy-based error-checking script for 1D streaming radiation datasets
# Daniel R. Reynolds, reynolds@smu.edu

import h5py
from pylab import *
import numpy as np

# set the solution tolerance
tol = 0.1

# load first dataset, and put 1D radiation field into 'Eg'
f = h5py.File('DD0001/data0001.cpu0000','r')
Eg = f.get('/Grid00000001/Grey_Radiation_Energy')

# set domain
N = Eg.shape
x = linspace(0.0, 1.0, N[0])

# set output flag
ret = 0

# compute error from analytical solution
Eg_anal = linspace(0.0, 1.0, N[0])
for i in range(0,N[0]):
    if (x[i] < 0.1):
        Eg_anal[i] = 1.0
    else:
        Eg_anal[i] = 0.0
Eg_err = Eg - Eg_anal
err_norm = (np.sum(np.multiply(Eg_err,Eg_err))/N[0])**(0.5)
if (err_norm > tol):
    ret += 1
f.close()


# repeat process for next dataset
f = h5py.File('DD0002/data0002.cpu0000','r')
Eg = f.get('/Grid00000001/Grey_Radiation_Energy')
Eg_anal = linspace(0.0, 1.0, N[0])
for i in range(0,N[0]):
    if (x[i] < 0.2):
        Eg_anal[i] = 1.0
    else:
        Eg_anal[i] = 0.0
Eg_err = Eg - Eg_anal
err_norm = (np.sum(np.multiply(Eg_err,Eg_err))/N[0])**(0.5)
if (err_norm > tol):
    ret += 1
f.close()


# repeat process for next dataset
f = h5py.File('DD0003/data0003.cpu0000','r')
Eg = f.get('/Grid00000001/Grey_Radiation_Energy')
Eg_anal = linspace(0.0, 1.0, N[0])
for i in range(0,N[0]):
    if (x[i] < 0.3):
        Eg_anal[i] = 1.0
    else:
        Eg_anal[i] = 0.0
Eg_err = Eg - Eg_anal
err_norm = (np.sum(np.multiply(Eg_err,Eg_err))/N[0])**(0.5)
if (err_norm > tol):
    ret += 1
f.close()


# repeat process for next dataset
f = h5py.File('DD0004/data0004.cpu0000','r')
Eg = f.get('/Grid00000001/Grey_Radiation_Energy')
Eg_anal = linspace(0.0, 1.0, N[0])
for i in range(0,N[0]):
    if (x[i] < 0.4):
        Eg_anal[i] = 1.0
    else:
        Eg_anal[i] = 0.0
Eg_err = Eg - Eg_anal
err_norm = (np.sum(np.multiply(Eg_err,Eg_err))/N[0])**(0.5)
if (err_norm > tol):
    ret += 1
f.close()


# repeat process for next dataset
f = h5py.File('DD0005/data0005.cpu0000','r')
Eg = f.get('/Grid00000001/Grey_Radiation_Energy')
Eg_anal = linspace(0.0, 1.0, N[0])
for i in range(0,N[0]):
    if (x[i] < 0.5):
        Eg_anal[i] = 1.0
    else:
        Eg_anal[i] = 0.0
Eg_err = Eg - Eg_anal
err_norm = (np.sum(np.multiply(Eg_err,Eg_err))/N[0])**(0.5)
if (err_norm > tol):
    ret += 1
f.close()


# repeat process for next dataset
f = h5py.File('DD0006/data0006.cpu0000','r')
Eg = f.get('/Grid00000001/Grey_Radiation_Energy')
Eg_anal = linspace(0.0, 1.0, N[0])
for i in range(0,N[0]):
    if (x[i] < 0.6):
        Eg_anal[i] = 1.0
    else:
        Eg_anal[i] = 0.0
Eg_err = Eg - Eg_anal
err_norm = (np.sum(np.multiply(Eg_err,Eg_err))/N[0])**(0.5)
if (err_norm > tol):
    ret += 1
f.close()


# repeat process for next dataset
f = h5py.File('DD0007/data0007.cpu0000','r')
Eg = f.get('/Grid00000001/Grey_Radiation_Energy')
Eg_anal = linspace(0.0, 1.0, N[0])
for i in range(0,N[0]):
    if (x[i] < 0.7):
        Eg_anal[i] = 1.0
    else:
        Eg_anal[i] = 0.0
Eg_err = Eg - Eg_anal
err_norm = (np.sum(np.multiply(Eg_err,Eg_err))/N[0])**(0.5)
if (err_norm > tol):
    ret += 1
f.close()


# repeat process for next dataset
f = h5py.File('DD0008/data0008.cpu0000','r')
Eg = f.get('/Grid00000001/Grey_Radiation_Energy')
Eg_anal = linspace(0.0, 1.0, N[0])
for i in range(0,N[0]):
    if (x[i] < 0.8):
        Eg_anal[i] = 1.0
    else:
        Eg_anal[i] = 0.0
Eg_err = Eg - Eg_anal
err_norm = (np.sum(np.multiply(Eg_err,Eg_err))/N[0])**(0.5)
if (err_norm > tol):
    ret += 1
f.close()


# repeat process for next dataset
f = h5py.File('DD0009/data0009.cpu0000','r')
Eg = f.get('/Grid00000001/Grey_Radiation_Energy')
Eg_anal = linspace(0.0, 1.0, N[0])
for i in range(0,N[0]):
    if (x[i] < 0.9):
        Eg_anal[i] = 1.0
    else:
        Eg_anal[i] = 0.0
Eg_err = Eg - Eg_anal
err_norm = (np.sum(np.multiply(Eg_err,Eg_err))/N[0])**(0.5)
if (err_norm > tol):
    ret += 1
f.close()


# repeat process for next dataset
f = h5py.File('DD0010/data0010.cpu0000','r')
Eg = f.get('/Grid00000001/Grey_Radiation_Energy')
Eg_anal = linspace(0.0, 1.0, N[0])
for i in range(0,N[0]):
    if (x[i] < 1.0):
        Eg_anal[i] = 1.0
    else:
        Eg_anal[i] = 0.0
Eg_err = Eg - Eg_anal
err_norm = (np.sum(np.multiply(Eg_err,Eg_err))/N[0])**(0.5)
if (err_norm > tol):
    ret += 1

# issue final success/failure statement
if (ret > 0):
    print 'Error tolerance ',tol,' exceeded for ',ret,'/10 datasets'
    print 'FAIL'
else:
    print 'PASS'


# end of script
