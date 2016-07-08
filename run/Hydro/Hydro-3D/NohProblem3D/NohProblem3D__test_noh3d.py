import pylab as pl
from yt.mods import *
from yt.testing import *
from yt.utilities.answer_testing.framework import \
     AnswerTestingTest, \
     sim_dir_load
from yt.frontends.enzo.answer_testing_support import \
     requires_outputlog

_pf_name = os.path.basename(os.path.dirname(__file__)) + ".enzo"
_dir_name = os.path.dirname(__file__)

class TestShockImage(AnswerTestingTest):
    _type_name = "noh3d_image"
    _attrs = ()

    def __init__(self, pf):
        self.ds = pf

    def run(self):
        sl = self.ds.slice(2, 0.5)
        frb = FixedResolutionBuffer(sl, (0.0, 1.0, 0.0, 1.0),
                                    (100, 100), antialias=False)
        dens = frb["Density"]
        return np.array([dens.mean(), dens.std(), dens.min(), dens.max()])

    def compare(self, new_result, old_result):
        tolerance = ytcfg.getint("yt", "answer_testing_tolerance")
        assert_allclose(new_result, old_result, rtol=10**-tolerance, atol=0)
    
class TestRadialDensity(AnswerTestingTest):
    _type_name = "noh3d_radial"
    _attrs = ()

    def __init__(self, pf):
        self.ds = pf
    
    def run(self):
        dd = self.ds.all_data()
        t = self.ds.parameters['InitialTime']
        x = dd['x']
        y = dd['y']
        z = dd['z']
        r = na.sqrt((x**2 + y**2 + z**2))

        postshock_r = na.linspace(t*1./3., na.sqrt(3))
        postshock_den = (1 + t/postshock_r)**2

        diag_r = r[(x==y)*(y==z)]
        diag_den = dd['Density'][(x==y)*(y==z)]
        return na.array(diag_den)

    def compare(self, new_result, old_result):
        tolerance = ytcfg.getint("yt", "answer_testing_tolerance")
        assert_allclose(new_result, old_result, rtol=10**-tolerance, atol=0)

    def plot(self):
        dd = self.ds.all_data()
        t = self.ds.parameters['InitialTime']
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
        pl.savefig('%s_density.png' % self.ds)
        pl.clf()
        
        return ['%s_density.png' % self.ds]

@requires_outputlog(_dir_name, _pf_name)
def test_noh3d():
    sim = sim_dir_load(_pf_name, path=_dir_name,
                       find_outputs=True)
    sim.get_time_series()
    pf = sim[-1]
    yield TestShockImage(pf)
    yield TestRadialDensity(pf)
