import sys
import numpy as N
import matplotlib.pyplot as plt

import matplotlib
matplotlib.rc('xtick', labelsize=19)
matplotlib.rc('ytick', labelsize=19)
matplotlib.rcParams.update({'font.size': 21})
matplotlib.rcParams.update({'lines.linewidth':2.0})
matplotlib.rcParams.update({'lines.markersize':10.0})
matplotlib.rcParams.update({'legend.numpoints': 1})
matplotlib.rcParams.update({'legend.fontsize': 19})
matplotlib.rcParams.update({'legend.frameon': False})
matplotlib.rcParams.update({'figure.autolayout': True})

#--------------------------------------------------------------
def get_offset(filename):
    FileIn = open(filename, "r")
    level = []
    time = []
    voffset = []
    root_cycle = []

    rec = -1
    for line in FileIn:
        if 'TopGrid' in line:
            rec += 1
        if 'Level' in line:
            lst = line.split()
            level.append(int(lst[2][0]))
            level.append(int(lst[2][0]))
        if 'PARTICLE' in line:
            lst = line.split()
            root_cycle.append(rec)
            time.append(float(lst[1]))
            vx = float(lst[5])
            vy = float(lst[6])
            vz = float(lst[7])
            dv2 = (vx - 1.0)**2 +  (vy - 1.0)**2 +  (vz - 1.0)**2
            dv = N.sqrt(dv2/3.)
            voffset.append(dv)
    FileIn.close()

    level = N.array(level)
    time = N.array(time)
    voffset = N.array(voffset)
    root_cycle = N.array(root_cycle)

    # only take the last point in a root grid cycle
    Ind = []
    Ind.append(0)
    for i in range(0,N.size(level),1):
        if (i<N.size(level)-1):
            if (root_cycle[i]<root_cycle[i+1]):
                Ind.append(i)
        else:
            Ind.append(i)

    Ind = N.array(Ind)
    level = level[Ind]
    time = time[Ind]
    voffset = voffset[Ind]
    root_cycle = root_cycle[Ind]

    return (level,time,voffset,root_cycle)
#--------------------------------------------------------------

filename1 = "./FastSubcycle/01.out"
filename2 = "./APMSubcycle/01.out"
fig_name  = "./SelfForce.png"

(level1,time1,voffset1,root_cycle1) = get_offset(filename1)
(level2,time2,voffset2,root_cycle2) = get_offset(filename2)

colors1 = ['bx','gx','mx','kx']
colors2 = ['b+','g+','m+','k+']
max_level1 = N.max(level1)
max_level2 = N.max(level2)
already_plot1 = N.zeros(max_level1+1)
already_plot2 = N.zeros(max_level2+1)

plt.figure(1)
for rec in range(0,N.max(level1)+1,1):
    Ind = N.where(level1 == rec)
    plt.plot(time1[Ind],voffset1[Ind]/1e-4,colors1[rec])
    if ((already_plot1[rec]==0) and (N.size(Ind)>0)):
        legen_str = 'level %i' %rec
        plt.plot(time1[Ind][0],voffset1[Ind][0]/1e-4,colors1[rec],label=legen_str)
        already_plot1[rec] = 1

for rec in range(0,N.max(level2)+1,1):
    Ind = N.where(level2 == rec)
    plt.plot(time2[Ind],voffset2[Ind]/1e-4,colors2[rec])
    if ((already_plot2[rec]==0) and (N.size(Ind)>0)):
        legen_str = 'level %i' %rec
        plt.plot(time2[Ind][0],voffset2[Ind][0]/1e-4,colors2[rec],label=legen_str)
        already_plot2[rec] = 1
plt.xlabel('$t$')
plt.ylabel('$\Delta v / 10^{-4}$')
plt.ylim(-0.1,1.5)
plt.xticks(N.arange(0.0, 0.411, 0.1))
plt.savefig(fig_name)
