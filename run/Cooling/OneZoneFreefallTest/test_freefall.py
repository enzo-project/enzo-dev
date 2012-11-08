from yt.mods import *
from yt.testing import *
from yt.utilities.answer_testing.api import AnswerTestingTest
from yt.utilities.answer_testing.framework import \
    requires_outputlog, \
    sim_dir_load

_rtol = 1.0e-2
_atol = 1.0e-7
_fields = ('Temperature', 'Dust_Temperature')

class FieldValuesTest(AnswerTestingTest):
    _type_name = "FreeFallFieldValues"
    _attrs = ("field", )

    def __init__(self, sim, field):
        self.pf = sim
        self.field = field
    
    def run(self):
        result = []
        for my_pf in self.pf:
            result.append(my_pf.h.all_data()[field])
        return result

    def compare(self, new_result, old_result):
        for i in range(len(new_result)):
            assert_allclose(new_result[i], old_result[i], _rtol, _atol)

@requires_outputlog(os.path.dirname(__file__),
                    "OneZoneFreefallTest.enzo")
def test_onezone_freefall():
    sim = sim_dir_load("OneZoneFreefallTest.enzo",
                       path="./Cooling/OneZoneFreefallTest",
                       find_outputs=True)
    sim.get_time_series()
    for field in _fields:
        yield FieldValuesTest(sim, field)
