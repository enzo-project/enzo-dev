from yt.mods import *
from yt.testing import *
from yt.utilities.answer_testing.framework import \
     AnswerTestingTest, \
     sim_dir_load
from yt.frontends.enzo.answer_testing_support import \
    requires_outputlog

_pf_name = os.path.basename(os.path.dirname(__file__)) + ".enzo"
_dir_name = os.path.dirname(__file__)

class TestAMRZeldovichMax(AnswerTestingTest):
    _type_name = "amr_zeldovich_max"
    _attrs = ()

    def __init__(self, sim):
        self.pf = sim
    
    def run(self):
        result = []
        for my_pf in self.pf:
            result.append(my_pf.h.find_max("Density")[0])
        return result

    def compare(self, new_result, old_result):
        for i in range(len(new_result)):
            assert_rel_equal(new_result[i], old_result[i], 4)

@requires_outputlog(_dir_name, _pf_name)
def test_collapse_max_value():
    sim = sim_dir_load(_pf_name, path=_dir_name, 
                       find_outputs=True)
    sim.get_time_series()
    
    yield TestAMRZeldovichMax(sim)
