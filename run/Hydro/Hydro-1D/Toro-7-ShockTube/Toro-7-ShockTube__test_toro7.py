from yt.mods import *
from yt.funcs import *
from yt.testing import *
from yt.frontends.enzo.answer_testing_support import \
    requires_outputlog, \
    ShockTubeTest
import os

_data_file = 'DD0001/data0001'
_solution_file = 'Toro-7-ShockTube_t=2.0_exact.txt'
_fields = ['Density','x-velocity','Pressure','ThermalEnergy']
_les = [0.0, 0.8]
_res = [0.6, 1.0]
_rtol = 1.0e-3
_atol = 1.0e-7

# Verifies that OutputLog exists
@requires_outputlog(os.path.dirname(__file__), "Toro-7-ShockTube.enzo")
def test_toro7():
    test = ShockTubeTest(_data_file, _solution_file, _fields, 
                         _les, _res, _rtol, _atol)
    return test()
