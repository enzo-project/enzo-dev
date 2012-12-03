from yt.mods import *
from yt.testing import *
from yt.utilities.answer_testing.framework import \
     AnswerTestingTest, \
     sim_dir_load
from yt.frontends.enzo.answer_testing_support import \
     requires_outputlog

_pf_name = os.path.basename(os.path.dirname(__file__)) + ".enzo"
_dir_name = os.path.dirname(__file__)

class TestFreeExpansionDistance(AnswerTestingTest):
    _type_name = "ExpansionDistance"
    _attrs = ()

    def __init__(self, pf):
        self.pf = pf
    
    def run(self):
        ray = self.pf.h.ray([0.0,0.0,0.5], [1.0,1.0,0.5])
        ray_length = np.sqrt(((ray.end_point - ray.start_point)**2).sum())
        ipos = na.argwhere(ray['VelocityMagnitude'] == 0.0)
        if len(ipos) > 0:
            ipos = ipos.min()
        else:
            ipos = -1
        return ray_length * ray['t'][ipos]

    def compare(self, new_result, old_result):
        tolerance = ytcfg.getint("yt", "answer_testing_tolerance")
        assert_allclose(new_result, old_result, 10**-tolerance, 0.0)

@requires_outputlog(_dir_name, _pf_name)
def test_collapse_max_value():
    sim = sim_dir_load(_pf_name, path=_dir_name, 
                       find_outputs=True)
    sim.get_time_series()
    pf = sim[-1]
    yield TestFreeExpansionDistance(pf)
