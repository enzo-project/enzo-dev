from yt.mods import *
from yt.testing import *
from yt.utilities.answer_testing.framework import \
     AnswerTestingTest, \
     sim_dir_load
from yt.frontends.enzo.answer_testing_support import \
     requires_outputlog

_pf_name = os.path.basename(os.path.dirname(__file__)) + ".enzo"
_dir_name = os.path.dirname(__file__)

class TestLVariation(AnswerTestingTest):
    _type_name = "deltaL"
    _attrs = ()

    def __init__(self, pf):
        self.pf = pf

    def run(self):
        ad = self.pf.h.all_data()
        return na.array(ad.quantities['TotalQuantity'](["AngularMomentumX", 
                                                        "AngularMomentumY", 
                                                        "AngularMomentumZ"]))

    def compare(self, new_result, old_result):
        tolerance = ytcfg.getint("yt", "answer_testing_tolerance")
        assert_allclose(new_result, old_result, rtol=10**-tolerance, atol=0)

@requires_outputlog(_dir_name, _pf_name)
def test_rotating_cylinder():
    sim = sim_dir_load(_pf_name, path=_dir_name,
                       find_outputs=True)
    sim.get_time_series()
    pf = sim[-1]
    yield TestLVariation(pf)
