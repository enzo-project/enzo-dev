from yt.mods import *
from yt.testing import *
from yt.utilities.answer_testing.framework import \
     AllFieldValuesTest, \
     sim_dir_load
from yt.frontends.enzo.answer_testing_support import \
     requires_outputlog

_fields = ("particle_position_x", "particle_position_y")
_pf_name = os.path.basename(os.path.dirname(__file__)) + ".enzo"
_dir_name = os.path.dirname(__file__)

@requires_outputlog(_dir_name, _pf_name)
def test_orbit():
    sim = sim_dir_load(_pf_name, path=_dir_name,
                       find_outputs=True)
    sim.get_time_series()
    tolerance = ytcfg.getint("yt", "answer_testing_tolerance")
    pf = sim[-1]
    for field in _fields:
        yield AllFieldValuesTest(pf, field, decimals=tolerance)
