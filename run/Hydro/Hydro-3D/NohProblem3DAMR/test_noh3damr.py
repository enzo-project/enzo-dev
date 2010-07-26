from yt.mods import *
from yt.extensions.enzo_test import YTStaticOutputTest
import pylab as pl
class TestShockImage(YTStaticOutputTest):
    name = "noh3damr_image"

    def run(self):
        # self.pf already exists
        sl = self.pf.h.slice(2, 0.5)
        frb = FixedResolutionBuffer(sl, (0.0, 1.0, 0.0, 1.0), (400, 400),antialias=False)
        self.result = frb["Density"]

    def compare(self, old_result):
        current_buffer = self.result
        old_buffer = old_result

        # We want our arrays to agree to some delta
        self.compare_array_delta(current_buffer, old_buffer, 5e-3)

    def plot(self):
        return []

class TestMaxValue(YTStaticOutputTest):
    name = "noh3damr_max"

    def run(self):
        # self.pf already exists
        self.result = self.pf.h.find_max("Density")[0]

    def compare(self, old_result):
        self.compare_value_delta(self.result, old_result, 5e-3)

    def plot(self):
        # There's not much to plot, so we just return an empty list.
        return []

class TestRadialDensity(YTStaticOutputTest):
    name = "noh3damr_radial"

    def run(self):
        # self.pf already exists
        dd = self.pf.h.all_data()
        t = self.pf['InitialTime']
        x = dd['x']
        y = dd['y']
        z = dd['z']
        r = na.sqrt((x**2 + y**2 + z**2))

        postshock_r = na.linspace(t*1./3., na.sqrt(3))
        postshock_den = (1 + t/postshock_r)**2

        diag_r = r[(x==y)*(y==z)]
        diag_den = dd['Density'][(x==y)*(y==z)]
        self.result = na.array(diag_den)

    def compare(self, old_result):
        current_buffer = na.array(self.result)
        old_buffer = na.array(old_result)

        # We want our arrays to agree to some delta
        self.compare_array_delta(current_buffer, old_buffer, 5e-3)

    def plot(self):
        dd = self.pf.h.all_data()
        t = self.pf['InitialTime']
        x = dd['x']
        y = dd['y']
        z = dd['z']
        r = na.sqrt((x**2 + y**2 + z**2))

        postshock_r = na.linspace(t*1./3., na.sqrt(3))
        postshock_den = (1 + t/postshock_r)**2

        all = pl.plot(r,dd['Density'],'b.',ms=4,alpha=0.2)
        diag_r = r[(x==y)*(y==z)]
        diag_den = dd['Density'][(x==y)*(y==z)]
        diag = pl.plot(diag_r,diag_den,'r.',ms=6,alpha=0.8)
        analy = pl.plot(postshock_r, postshock_den, 'k-')

        pl.legend([all, diag, analy],['Simulation All','Simulation Diagonal','Analytical'])
        
        im = pl.plot([0.0,t*1./3.],[64.0,64.0],'k-')
        pl.xlim(0.0,na.sqrt(3.0))
        pl.xlabel('r')
        pl.ylabel('Density')
        pl.savefig('%s_density.png'%self.pf)
        pl.clf()
        
        return ['%s_density.png'%self.pf]
        # There's not much to plot, so we just return an empty list.
