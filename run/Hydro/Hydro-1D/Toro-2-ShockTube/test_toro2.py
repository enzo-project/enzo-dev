from yt.mods import *
from yt.funcs import *
from yt.testing import *
from yt.frontends.enzo.answer_testing_support import \
    requires_outputlog, \
    ShockTubeTest
import os

_data_file = 'DD0001/data0001'
_solution_file = 'Toro-2-ShockTube_t=0.15_exact.txt'
_fields = ['Density','Pressure']
_les = [0.2, 0.7]
_res = [0.3, 0.8]
_rtol = 1.0e-1
_atol = 1.0e-7

# Verifies that OutputLog exists
@requires_outputlog(os.path.dirname(__file__), "Toro-2-ShockTube.enzo")
def test_toro2():
    test = ShockTubeTest(_data_file, _solution_file, _fields, 
                         _les, _res, _rtol, _atol)
    return test()
