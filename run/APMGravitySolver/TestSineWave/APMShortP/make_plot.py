import sys
import numpy as np
import matplotlib.pyplot as plt
from yt.mods import *
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
period = 0.2
solver = 'APM'
data = './DD0001/data0001'
#--------------------------------------------------------------

#--------------------------------------------------------------
def anal_solution_acc(period,x):
    val = 4*np.pi*period/(2.*np.pi) * np.cos(x*2*np.pi/period)
    return val


def anal_solution_pot(period,x):
    val = -4*np.pi*(period/(2.*np.pi))**2 * np.sin(x*2*np.pi/period)
    return val
#--------------------------------------------------------------

print("loading %s" % data)
pf = load(data)
resol = pf.domain_dimensions[0]
pf.h.print_stats()
ngrids_per_level = np.zeros(np.max(pf.h.ds.index.grid_levels)+1)
# Get grids per level
for grid in pf.h.ds.index.grids:
    ngrids_per_level[grid.Level] += 1

#----------------------------------------------------
# Gas from AccelerationField.out
#----------------------------------------------------
Files = glob.glob('AccelerationField.out*_p0000')
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

for name in Files:
    Data = np.loadtxt(name)
    level_gas = level_gas + list(Data[:,0])
    is_ghost_gas = is_ghost_gas + list(Data[:,1])
    i_gas = i_gas + list(Data[:,2])
    j_gas = j_gas + list(Data[:,3])
    k_gas = k_gas + list(Data[:,4])
    x_gas = x_gas + list(Data[:,5]+0.5) # so coordinates go from 0 to 1 instead of -0.5 to 0.5
    y_gas = y_gas + list(Data[:,6]+0.5)
    z_gas = z_gas + list(Data[:,7]+0.5)
    r_gas = r_gas + list(Data[:,8])
    ax_gas = ax_gas + list(Data[:,9])
    ay_gas = ay_gas + list(Data[:,10])
    az_gas = az_gas + list(Data[:,11])

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

# Remove all ghost zones
Ind = []
for i in range(0,np.size(x_gas)):
    if (is_ghost_gas[i] == 0):
        Ind.append(i)
Ind = np.array(Ind)

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

ngas = np.size(x_gas)

# Analytical solution
x_anal = np.arange(-0.1,1.1,1e-3)
n_anal = np.size(x_anal)
f_anal = np.zeros(n_anal)
for i in range(0,n_anal,1):
    f_anal[i] = anal_solution_acc(period,x_anal[i])

# RMS
f_true = np.zeros(ngas)
for i in range(0,ngas,1):
    f_true[i] = anal_solution_acc(period,x_gas[i])

atan_gas = (ay_gas**2. + az_gas**2.)**0.5
Error = np.abs((ax_gas-f_true)/f_true)
Mean = np.mean(Error)
Rms = np.std(Error)
Error2 = np.abs(atan_gas/f_true)
Mean2 = np.mean(Error2)
Rms2 = np.std(Error2)

# Plot
symbols = ['bx','ro','g.', 'c*']

# limits
if (period == 0.2):
    xmin = np.min(x_gas)-1e-2
    xmax = np.max(x_gas)+1e-2
    ymin_anal = -0.5
    ymax_anal = 0.5
    if(resol == 32):
        ymin = 0.01
        ymax = 0.5
    else:
        ymin = 5e-3
        ymax = 0.2
else: # period == 1.0
    xmin = 0.28
    xmax = 0.72
    ymin_anal = -2.1
    ymax_anal = -0.5
    if (resol == 32):
        ymin = 4e-4
        ymax = 1e-2
    else:
        ymin = 7e-5
        ymax = 2e-3

# Accelerations all, png
plt.figure(1)
plt.plot(x_anal,f_anal,'k',label='$F_{\mathrm{anal}}$')
for i in range(0,int(max(level_gas)+1),1):
    Ind = np.where(level_gas == i)
    txt = '$F_{\mathrm{rad, %i}}$ (%i grids)' %(i,ngrids_per_level[i])
    plt.plot(x_gas[Ind],ax_gas[Ind],symbols[i],label=txt)

plt.xlim(xmin,xmax)
if (period == 0.2):
    plt.ylim(-0.5,0.5)
else:
    plt.ylim(-2.1,-0.5)
plt.xlabel('$x$')
plt.ylabel('$\mathrm{Force}$')
plt.savefig('Forces-'+str(period)+'-'+solver+'.png',bbox_inches='tight')

# Errors all, png
fig, ax1 = plt.subplots()
for i in range(0,int(max(level_gas)+1),1):
    Ind = np.where(level_gas == i)
    txt = '$F_{\mathrm{rad, %i}}$ (%i grids)' %(i,ngrids_per_level[i])
    ax1.semilogy(x_gas[Ind],Error[Ind],symbols[i],label=txt)
ax1.set_xlabel('$x$')
ax1.set_ylabel('$\mathrm{Force\,Error}$')
ax1.set_xlim(xmin,xmax)

ax2 = ax1.twinx()
ax2.plot(x_anal, f_anal, 'k--')
ax2.set_xlim(xmin,xmax)
ax2.set_ylim(ymin_anal,ymax_anal)
ax2.set_ylabel('$\mathrm{Force}$', color='k')
plt.savefig('Errors-'+str(period)+'-'+solver+'.png',bbox_inches='tight')
#----------------------------------------------------
