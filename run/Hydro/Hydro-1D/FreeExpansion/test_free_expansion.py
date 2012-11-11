from yt.mods import *
from yt.testing import *
from yt.utilities.answer_testing.framework import \
     AnswerTestingTest, \
     requires_outputlog, \
     sim_dir_load

_pf_name = os.path.basename(os.path.dirname(__file__)) + ".enzo"
_dir_name = os.path.dirname(__file__)

class TestFreeExpansionDistance(YTStaticOutputTest):
    _type_name = "ExpansionDistance"
    _attrs = ()

    def __init__(self, pf):
        self.pf = pf
    
    def run(self):
        ray = self.pf.h.ray([0.0,0.5,0.5], [1.0,0.5,0.5])
        ray_length = np.sqrt(((q.end_point - q.start_point)**2).sum())
        ipos = na.where(ray['VelocityMagnitude'] == 0.0)[0].argmin()
        return ray_length * ray['t'][ipos]

    def compare(self, new_result, old_result):
        assert_rel_equal(new_result, old_result, 2,
                         verbose=True)

@requires_outputlog(_dir_name, _pf_name)
def test_collapse_max_value():
    sim = sim_dir_load(_pf_name, path=_dir_name, 
                       find_outputs=True)
    sim.get_time_series()
    for pf in sim:
        yield TestFreeExpansionDistance(pf)
