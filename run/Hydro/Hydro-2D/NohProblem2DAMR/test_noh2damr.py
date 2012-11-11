import pylab as pl
from yt.mods import *
from yt.testing import *
from yt.utilities.answer_testing.framework import \
     AnswerTestingTest, \
     requires_outputlog, \
     sim_dir_load

_pf_name = os.path.basename(os.path.dirname(__file__)) + ".enzo"
_dir_name = os.path.dirname(__file__)

class TestShockImage(AnswerTestingTest):
    _type_name = "noh2damr_image"
    _attrs = ()

    def __init__(self, pf):
        self.pf = pf
    
    def run(self):
        # self.pf already exists
        sl = self.pf.h.slice(2, 0.5)
        frb = FixedResolutionBuffer(sl, (0.0, 1.0, 0.0, 1.0), 
                                    (400, 400), antialias=False)
        return frb["Density"]

    def compare(self, new_result, old_result):
        assert_rel_equal(new_result, old_result, 3)

class TestRadialDensity(AnswerTestingTest):
    _type_name = "noh2damr_radial"
    _attrs = ()

    def __init__(self, pf):
        self.pf = pf
    
    def run(self):
        # self.pf already exists
        dd = self.pf.h.all_data()
        t = self.pf['InitialTime']
        x = dd['x']
        y = dd['y']
        r = na.sqrt((x**2 + y**2))

        postshock_r = na.linspace(t*1./3., na.sqrt(2))
        postshock_den = (1 + t/postshock_r)

        diag_r = r[(x==y)]
        diag_den = dd['Density'][(x==y)]
        return na.array(diag_den)

    def compare(self, new_result, old_result):
        assert_rel_equal(new_result, old_result, 3)

    def plot(self):
        dd = self.pf.h.all_data()
        t = self.pf['InitialTime']
        x = dd['x']
        y = dd['y']
        r = na.sqrt((x**2 + y**2))

        postshock_r = na.linspace(t*1./3., na.sqrt(2))
        postshock_den = (1 + t/postshock_r)

        all = pl.plot(r,dd['Density'],'b.',ms=4,alpha=0.2)
        diag_r = r[(x==y)]
        diag_den = dd['Density'][(x==y)]
        diag = pl.plot(diag_r,diag_den,'r.',ms=6,alpha=0.8)
        analy = pl.plot(postshock_r, postshock_den, 'k-')

        pl.legend([all, diag, analy],['Simulation All','Simulation Diagonal','Analytical'])
        im = pl.plot([0.0,t*1./3.],[16.0,16.0],'k-')
        pl.xlim(0.0,na.sqrt(2.0))
        pl.xlabel('r')
        pl.ylabel('Density')
        pl.savefig('%s_density.png'%self.pf)
        pl.clf()
        
        return ['%s_density.png'%self.pf]
        # There's not much to plot, so we just return an empty list.

@requires_outputlog(_dir_name, _pf_name)
def test_cooling_time():
    sim = sim_dir_load(_pf_name, path=_dir_name,
                       find_outputs=True)
    sim.get_time_series()
    for pf in sim:
        yield TestShockImage(pf)
        yield TestRadialDensity(pf)
