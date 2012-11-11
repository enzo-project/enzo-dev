from yt.mods import *
from yt.testing import *
from yt.utilities.answer_testing.framework import \
     AnswerTestingTest, \
     requires_outputlog, \
     sim_dir_load

_pf_name = os.path.basename(os.path.dirname(__file__)) + ".enzo"
_dir_name = os.path.dirname(__file__)
_fields = ('Density', 'HI_Fraction', 'HII_Fraction', 'HI_kph')

class TestAMRPhotonTest(AnswerTestingTest):
    _type_name = "photon_shadowing_image"
    _attrs = ("field", )

    def __init__(self, pf, field):
        self.pf = pf
        self.field = field

    def run(self):
        # self.pf already exists
        sl = self.pf.h.slice(2,0.5)
        frb = FixedResolutionBuffer(sl, (0,1,0,1), (200,200))
        return frb[self.field]

    def compare(self, new_result, old_result):
        assert_rel_equal(new_result, old_result, 3)

@requires_outputlog(_dir_name, _pf_name)
def test_cooling_time():
    sim = sim_dir_load(_pf_name, path=_dir_name,
                       find_outputs=True)
    sim.get_time_series()
    for pf in sim:
        for field in _fields:
            yield TestAMRPhotonTest(pf)
