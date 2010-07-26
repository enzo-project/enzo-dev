from yt.mods import *
import pylab as pl

for dump in range(7):
    pf = EnzoStaticOutput('DD%04i/noh2D_%04i' % (dump,dump))
    dd = pf.h.all_data()
    t = pf['InitialTime']
    x = dd['x']
    y = dd['y']
    r = na.sqrt((x**2 + y**2))

    postshock_r = na.linspace(t*1./3., na.sqrt(2))
    postshock_den = (1 + t/postshock_r)

    pl.plot(r,dd['Density'],'b.',ms=4,alpha=0.2)

    pl.plot(r[(x==y)],dd['Density'][(x==y)],'r.',ms=6,alpha=0.8)

    pl.plot(postshock_r, postshock_den, 'k-')
    pl.legend(['Simulation All','Simulation Diagonal','Analytical'])
    pl.plot([0.0,t*1./3.],[16.0,16.0],'k-')
    pl.xlim(0.0,na.sqrt(2.0))
    pl.xlabel('r')
    pl.ylabel('Density')
    pl.savefig('%s_density.png'%pf)
    pl.clf()

    pc = PlotCollection(pf,center=[0.5,0.5,0.0])
    pc.add_slice('Density',2)
    pc.set_zlim(1.0,18.0,nticks=18)
    pc.save()
