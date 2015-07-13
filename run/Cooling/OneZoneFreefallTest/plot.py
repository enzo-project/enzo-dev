from matplotlib import pylab
import sys

from yt.mods import *

do_fH2 = True
do_t_cool = True

dust = True
if dust:
    keyword = 'with_dust'
else:
    keyword = 'without_dust'

es = EnzoSimulation('OneZoneFreefallTest.enzo', find_outputs=True)
es.get_time_series()

T = []
n = []
Z = []
fH2 = []
Tdust = []
t_cool = []
t_dyn = []

for pf in es:
    T.append(pf.h.grids[0]['Temperature'])
    n.append(pf.h.grids[0]['NumberDensity'])
    Z.append(pf.h.grids[0]['Metallicity'])
    if do_fH2: fH2.append(pf.h.grids[0]['H2I_Fraction'])
    if do_t_cool:
        t_cool.append(pf.h.grids[0]['Cooling_Time'])
        t_dyn.append(pf.h.grids[0]['DynamicalTime'])
    if dust: Tdust.append(pf.h.grids[0]['Dust_Temperature'])
    del pf

T = na.array(T)
n = na.array(n)
Z = na.array(Z)
fH2 = na.array(fH2)
Tdust = na.array(Tdust)
t_cool = na.array(t_cool)
t_dyn = na.array(t_dyn)

colors = ['black', 'purple', 'blue', 'green', 'orange', 'red']

met = na.round(na.log10(Z[0,0,:,0]))
for i in range(T.shape[2]):
    pylab.loglog(n[:, 0, i, 0], T[:, 0, i, 0], 
                 label='log (Z/Z$_{\odot}$) = %d' % met[i],
                 color=colors[i], linestyle='-')
    if dust:
        pylab.loglog(n[:, 0, i, 0], Tdust[:, 0, i, 0], 
                     color=colors[i], linestyle='--')
pylab.xlim(xmin=1.0)
pylab.ylim(1e0, 1e4)
pylab.xlabel('n [cm$^{-3}$]')
pylab.ylabel('T [K]')
pylab.legend(labelspacing=0.0, loc='lower right')
pylab.savefig('n-T_%s.png' % keyword)
pylab.clf()

if do_fH2:
    for i in range(T.shape[2]):
        pylab.loglog(n[:, 0, i, 0], fH2[:, 0, i, 0], 
                     label='log (Z/Z$_{\odot}$) = %d' % met[i],
                     color=colors[i])
    pylab.xlim(xmin=1.0)
    pylab.xlabel('n [cm$^{-3}$]')
    pylab.ylabel('f$_{H_{2}}$')
    pylab.ylim(1e-5, 1)
    pylab.legend(labelspacing=0.0, loc='lower right')
    pylab.savefig('n-fH2_%s.png' % keyword)
    pylab.clf()

if do_t_cool:
    for i in range(T.shape[2]):
        pylab.loglog(n[:, 0, i, 0], (t_cool[:, 0, i, 0]/t_dyn[:, 0, i, 0]), 
                     label='log (Z/Z$_{\odot}$) = %d' % met[i],
                     color=colors[i])
    pylab.xlim(xmin=1.0)
    pylab.xlabel('n [cm$^{-3}$]')
    pylab.ylabel('t$_{cool}$ / t$_{dyn}$')
    pylab.legend(labelspacing=0.0, loc='lower right')
    pylab.savefig('n-t_cool_t_dyn_%s.png' % keyword)
    pylab.clf()
