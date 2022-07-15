import sys
import numpy as np
import matplotlib.pyplot as plt
from yt.mods import *
import random
import glob

import matplotlib
matplotlib.rc('xtick', labelsize=19)
matplotlib.rc('ytick', labelsize=19)
matplotlib.rcParams.update({'font.size': 21})
matplotlib.rcParams.update({'lines.linewidth':2.0})
matplotlib.rcParams.update({'legend.numpoints': 1})
matplotlib.rcParams.update({'legend.fontsize': 19})
matplotlib.rcParams.update({'legend.frameon': False})

#--------------------------------------------------------------
sphere_type = 2
solver = 'APM'
data = './DD0000/data0000'
#--------------------------------------------------------------

#--------------------------------------------------------------
# constant density
def phi_anal0(r,Mstar,Rstar):
    if (r > Rstar):
        pot = -Mstar/r
    else:
        pot = -Mstar/(2*Rstar) * (3-r**2/Rstar**2)
    return pot

def acc_anal0(r,Mstar,Rstar):
    if (r > Rstar):
        acc = Mstar/r**2
    else:
        acc = Mstar*r/Rstar**3
    return acc

# 1/r**2
def phi_anal1(r,Mstar,Rstar):
    if (r > Rstar):
        pot = -Mstar/r
    else:
        pot = -Mstar/Rstar * (1 - np.log(r/Rstar))
    return pot

def acc_anal1(r,Mstar,Rstar):
    if (r > Rstar):
        acc = Mstar/r**2
    else:
        acc = Mstar/Rstar * 1/r
    return acc

# Plummer
def phi_anal2(r,Mstar,Rstar):
    if (r > Rstar):
        pot = -Mstar/r
    else:
        pot = -Mstar/Rstar * ( 2*np.sqrt(2)/np.sqrt(1.+r**2/Rstar**2) - 1 )
    return pot
def acc_anal2(r,Mstar,Rstar):
    if (r > Rstar):
        acc = Mstar/r**2
    else:
        acc = 2*np.sqrt(2)*Mstar/Rstar**3 * r * (1.+r**2/Rstar**2)**-1.5
    return acc

# Get Accelerations and Errors
def get_accelerations(data,Rho,Rstar):

    print("loading %s" % data)
    pf = load(data)
    pf.h.print_stats()
    ngrids_per_level = np.zeros(np.max(pf.h.ds.index.grid_levels)+1)
    # Get grids per level
    for grid in pf.h.ds.index.grids:
        ngrids_per_level[grid.Level] += 1

    #----------------------------------------------------
    # Gas from AccelerationField.out
    #----------------------------------------------------
    Files = glob.glob('AccelerationField.out*')
    procs = []
    level_gas = []
    is_ghost_gas = []
    i_gas = []
    j_gas = []
    k_gas = []
    x_gas = []
    y_gas = []
    z_gas = []
    r_gas = []
    ax_gas = []
    ay_gas = []
    az_gas = []
    arad_gas = []
    atan_gas = []

    for name in Files:
        Data = np.loadtxt(name)
        procs = procs + list(int(name[-4::])*np.ones(np.size(Data[:,0]),dtype=int))
        level_gas = level_gas + list(Data[:,0])
        is_ghost_gas = is_ghost_gas + list(Data[:,1])
        i_gas = i_gas + list(Data[:,2])
        j_gas = j_gas + list(Data[:,3])
        k_gas = k_gas + list(Data[:,4])
        x_gas = x_gas + list(Data[:,5])
        y_gas = y_gas + list(Data[:,6])
        z_gas = z_gas + list(Data[:,7])
        r_gas = r_gas + list(Data[:,8])
        ax_gas = ax_gas + list(Data[:,9])
        ay_gas = ay_gas + list(Data[:,10])
        az_gas = az_gas + list(Data[:,11])
        atan_gas = atan_gas + list(Data[:,12])
        arad_gas = arad_gas + list(Data[:,13])

    procs = np.array(procs)
    level_gas = np.array(level_gas)
    is_ghost_gas = np.array(is_ghost_gas)
    level_gas = level_gas.astype(int)
    is_ghost_gas = is_ghost_gas.astype(int)
    i_gas = np.array(i_gas)
    j_gas = np.array(j_gas)
    k_gas = np.array(k_gas)
    x_gas = np.array(x_gas)
    y_gas = np.array(y_gas)
    z_gas = np.array(z_gas)
    r_gas = np.array(r_gas)
    ax_gas = np.array(ax_gas)
    ay_gas = np.array(ay_gas)
    az_gas = np.array(az_gas)
    arad_gas = np.array(arad_gas)
    atan_gas = np.array(atan_gas)

    # Remove all ghost zones
    #Ind = ((is_ghost_gas == 0) | (level_gas > 0)).nonzero()[0]
    Ind = []
    for i in range(0,np.size(x_gas)):
        if (is_ghost_gas[i] == 0):
            Ind.append(i)
    Ind = np.array(Ind)

    procs = procs[Ind]
    level_gas = level_gas[Ind]
    is_ghost_gas = is_ghost_gas[Ind]
    i_gas = i_gas[Ind]
    j_gas = j_gas[Ind]
    k_gas = k_gas[Ind]
    x_gas = x_gas[Ind]
    y_gas = y_gas[Ind]
    z_gas = z_gas[Ind]
    r_gas = r_gas[Ind]
    ax_gas = ax_gas[Ind]
    ay_gas = ay_gas[Ind]
    az_gas = az_gas[Ind]
    arad_gas = arad_gas[Ind]
    atan_gas = atan_gas[Ind]

    ngas = np.size(x_gas)

    # RMS
    f_true = np.zeros(ngas)
    Error1 = np.zeros(ngas)
    Error2 = np.zeros(ngas)
    for i in range(0,ngas,1):
        if (sphere_type == 0):
            Mstar = 4/3. * np.pi * Rho * Rstar**3
            f_true[i] =  acc_anal0(r_gas[i],Mstar,Rstar)
        elif (sphere_type == 1):
            Mstar = 4 * np.pi * Rho *  Rstar**3
            f_true[i] =  acc_anal1(r_gas[i],Mstar,Rstar)
        elif (sphere_type == 2):
            Mstar = 4/3. * np.pi * Rho * Rstar**3 / (2.*np.sqrt(2.))
            f_true[i] =  acc_anal2(r_gas[i],Mstar,Rstar)

    Error1 = np.abs((arad_gas-f_true)/f_true)
    Error2 = np.abs(atan_gas/f_true)
    Mean1 = np.mean(Error1)
    Mean2 = np.mean(Error2)
    Rms1 = np.std(Error1)
    Rms2 = np.std(Error2)

    return (ngrids_per_level,procs,level_gas,r_gas,arad_gas,
            Error1,Error2,Mean1,Mean2,Rms1,Rms2)

# Values for analytical solution
Rho = 1.0
Rstar = 0.3

# Analytical solution
r_anal = np.arange(1e-3,1.1,1e-3)
n_anal = np.size(r_anal)
f_anal = np.zeros(n_anal)
for i in range(0,n_anal,1):
    if (sphere_type == 0):
        Mstar = 4/3. * np.pi * Rho * Rstar**3
        f_anal[i] =  acc_anal0(r_anal[i],Mstar,Rstar)
    elif (sphere_type == 1):
        Mstar = 4 * np.pi * Rho *  Rstar**3
        f_anal[i] =  acc_anal1(r_anal[i],Mstar,Rstar)
    elif (sphere_type == 2):
        Mstar = 4/3. * np.pi * Rho * Rstar**3 / (2.*np.sqrt(2.))
        f_anal[i] =  acc_anal2(r_anal[i],Mstar,Rstar)
    else:
        print("Wrong sphere type")
        sys.exit()

# Get accelerations, errors, etc...
(ngrids1,procs1,l1,r1,a1,er1,et1,mr1,mt1,rmsr1,rmst1) = get_accelerations(data,Rho,Rstar)

# output relevant values
print("Radial force error = ", mr1, " +/- ", rmsr1)
print("Tangential force error = ", mt1, " +/- ", rmst1)

# Plot
symbols = ['bx','ro','g.', 'c*']
if (sphere_type == 0):
    ylim_min = 1.1e-8
    ylim_max = 0.99
    legend_where = 'upper left'
elif  (sphere_type == 1):
    ylim_min = 1.1e-5
    ylim_max = 0.99
    legend_where = 'lower left'
else:
    ylim_min = 1.1e-6
    ylim_max = 0.99
    legend_where = 'upper left'

# All data (png)
# Radial acceleration
plt.figure(1)
plt.loglog(r_anal[0],f_anal[0],'k',label='$F_{\mathrm{anal}}$')
for i in range(0,int(max(l1)+1),1):
    Ind2 = np.where(l1 == i)
    txt = '$F_{\mathrm{rad, %i}}$ (%i grids)' %(i,ngrids1[i])
    plt.loglog(r1[Ind2],a1[Ind2],symbols[i],label=txt)
plt.loglog(r_anal,f_anal,'k')

# Limits
plt.xlim(6e-3,1.0)
if (sphere_type == 0):
    plt.ylim(2e-2,2.0)
elif  (sphere_type == 2):
    plt.ylim(2e-2,0.7)
else:
    plt.ylim(0.3,1e2)
plt.xlabel('$r$')
plt.ylabel('$\mathrm{Force}$')
plt.legend(loc=legend_where)
plt.vlines(Rstar,1e-8,1e99,'k')
plt.savefig('Force-'+str(sphere_type)+'-'+solver+'.png')

# Errors
plt.figure(2)
plt.subplots(2,1,sharex=True)
plt.subplots_adjust(hspace=0)
plt.subplot(211)
for i in range(0,int(max(l1)+1),1):
    k = int(max(l1))-i
    Ind2 = np.where(l1 == k)
    plt.loglog(r1[Ind2],er1[Ind2],symbols[k])
plt.vlines(Rstar,1e-8,1.0,'k')
plt.xlim(6e-3,1.0)
plt.ylim(ylim_min,ylim_max)
plt.ylabel('$\mathrm{Radial\,error}$')
plt.subplot(212)
for i in range(0,int(max(l1)+1),1):
    k = int(max(l1))-i
    Ind2 = np.where(l1 == k)
    plt.loglog(r1[Ind2],et1[Ind2],symbols[k])
plt.vlines(Rstar,1e-8,1.0,'k')
plt.xlim(6e-3,1.0)
plt.ylim(ylim_min,ylim_max)
plt.xlabel('$r$')
plt.ylabel('$\mathrm{Tangential\,error}$')
plt.savefig('Errors-'+str(sphere_type)+'-'+solver+'.png')
#----------------------------------------------------
