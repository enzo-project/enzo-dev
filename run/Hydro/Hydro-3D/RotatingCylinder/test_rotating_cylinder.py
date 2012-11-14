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
        # self.pf already exists
        data = self.pf.h.all_data()
        AngMom = []
        AngMom.append(data.quantities["TotalQuantity"]("AngularMomentum")[0])
        AngMomInitial = AngMom[0]
        AngMomPercentageChange = []
        for i, item in enumerate(AngMom):
            AngMomPercentageChange.append(100.0*(item - AngMomInitial)/AngMomInitial)
        return max(AngMomPercentageChange)

    def compare(self, new_result, old_result):
        assert_allclose(new_result, old_result, rtol=1e-3, atol=0)

@requires_outputlog(_dir_name, _pf_name)
def test_rotating_cylinder():
    sim = sim_dir_load(_pf_name, path=_dir_name,
                       find_outputs=True)
    sim.get_time_series()
    for pf in sim:
        yield TestLVariation(pf)
