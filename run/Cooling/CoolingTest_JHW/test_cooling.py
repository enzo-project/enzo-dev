from yt.mods import *
from yt.testing import *
from yt.utilities.answer_testing.framework import \
     FieldValuesTest, \
     requires_outputlog, \
     sim_dir_load

_fields = ("Cooling_Time",)
_pf_name = "CoolingTest_JHW.enzo"
_dir_name = os.path.dirname(__file__)

@requires_outputlog(_dir_name, _pf_name)
def test_cooling_time():
    sim = sim_dir_load(_pf_name, path=_dir_name,
                       find_outputs=True)
    sim.get_time_series()
    for pf in sim:
        for field in _fields:
            yield FieldValuesTest(pf, field)

