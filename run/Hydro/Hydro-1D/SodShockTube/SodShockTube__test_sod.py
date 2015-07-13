from yt.mods import *
from yt.funcs import *
from yt.testing import *
from yt.frontends.enzo.answer_testing_support import \
    requires_outputlog, \
    ShockTubeTest
import os

_data_file = 'DD0001/data0001'
_solution_file = 'SodShockTube_t=0.25_exact.txt'
_fields = ['Density','ThermalEnergy']
_les = [0.25, 0.85]
_res = [0.4, 0.9]
_rtol = 1.0e-2
_atol = 1.0e-7

# Verifies that OutputLog exists
@requires_outputlog(os.path.dirname(__file__), "SodShockTube.enzo")
def test_sod():
    test = ShockTubeTest(_data_file, _solution_file, _fields, 
                         _les, _res, _rtol, _atol)
    return test()
