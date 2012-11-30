from yt.mods import *
from yt.funcs import *
from yt.testing import *
from yt.frontends.enzo.answer_testing_support import \
    requires_outputlog, \
    ShockTubeTest, \
    standard_small_simulation
from yt.utilities.answer_testing.framework import \
    sim_dir_load, \
    VerifySimulationSameTest
import os


_data_file = 'DD0001/data0001'
_solution_file = 'Toro-6-ShockTube_t=2.0_exact.txt'
_fields = ['Density','x-velocity','Pressure','ThermalEnergy']
_les = [0.0]
_res = [1.0]
_rtol = 1.0e-6
_atol = 1.0e-7

# Verifies that OutputLog exists
@requires_outputlog(os.path.dirname(__file__), "Toro-6-ShockTube.enzo")
def test_toro6():
    test = ShockTubeTest(_data_file, _solution_file, _fields, 
                         _les, _res, _rtol, _atol)
    return test()

_base_fields = ('Density', 'TotalEnergy')

def test_almost_standard():
    sim = sim_dir_load("Toro-6-ShockTube.enzo",
                       path="./Hydro/Hydro-1D/Toro-6-ShockTube",
                       find_outputs=True)
    sim.get_time_series()
    yield VerifySimulationSameTest(sim)
    base_pf = sim[0]
    fields = [f for f in _base_fields if f in base_pf.h.field_list]
    # Only test the last output.
    pf = sim[-1]
    for test in standard_small_simulation(pf, fields): yield test

# Tests that OutputLog exists and fails otherwise
def test_exist():
    filename = os.path.dirname(__file__) + "/OutputLog"
    if not os.path.exists(filename):
        raise EnzoTestOutputFileNonExistent(filename)
