from yt.mods import *
from yt.testing import *
from yt.utilities.answer_testing.api import AnswerTestingTest
from yt.utilities.answer_testing.framework import \
    requires_outputlog, \
    sim_dir_load

class TestCollapseMaxValue(AnswerTestingTest):
    _type_name = "MaxValue"
    _attrs = ()

    def __init__(self, sim):
        self.pf = sim
    
    def run(self):
        result = []
        for my_pf in pf:
            result.append(self.my_pf.h.find_max("Density")[0])
        return result

    def compare(self, new_result, old_result):
        for i in range(len(new_result)):
            assert_rel_equal(new_result[i], old_result[i], 2)

@requires_outputlog("./Hydro/Hydro-3D/CollapseTestNonCosmological", 
                    "CollapseTestNonCosmological.enzo"
def test_collapse_max_value():
    sim = sim_dir_load("CollapseTestNonCosmological.enzo", 
                       path="./Hydro/Hydro-3D/CollapseTestNonCosmological", 
                       find_outputs=True)
    sim.get_time_series()
    
    yield TestCollapseMaxValue(sim)
