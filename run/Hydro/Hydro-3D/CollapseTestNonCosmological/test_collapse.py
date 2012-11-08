from yt.mods import *
from yt.testing import *
from yt.utilities.answer_testing.api import YTStaticOutputTest

class TestCollapseMaxValue(YTStaticOutputTest):
    name = "collapse_max"

    def __init__(self, sim):
        self.pf = sim
    
    def run(self):
        result = []
        for my_pf in pf:
            result.append(self.my_pf.h.find_max("Density")[0])

    def compare(self, new_result, old_result):
        for i in range(len(new_result)):
            assert_rel_equal(new_result[i], old_result[i], 2)

def test_collapse_max_value():
    # Create a time series object.
    sim = simulation('CollapseTestNonCosmological.enzo', 'Enzo',
                     find_outputs=True)
    sim.get_time_series()
    
    yield TestCollapseMaxValue(sim)
