###
### MHDCT Test
###

import os
from yt.mods import *
from yt.testing import *
from yt.utilities.answer_testing.framework import \
    VerifySimulationSameTest, \
    sim_dir_load
from yt.frontends.enzo.answer_testing_support import \
    requires_outputlog, \
    standard_small_simulation

_base_fields = ("Density",
                "pressure",
                "x-velocity",
                "y-velocity",
                "Bx",
                "By")

@requires_outputlog(os.path.dirname(__file__), "MHDCTOrszagTangAMR.enzo") # Verifies that OutputLog exists
def test_standard():
    sim = sim_dir_load("MHDCTOrszagTangAMR.enzo",
                       path="./MHD/2D/MHDCTOrszagTangAMR",
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
    assert os.path.exists(filename)

@requires_outputlog(os.path.dirname(__file__), "MHDCTOrszagTangAMR.enzo") # Verifies that OutputLog exists
def test_DivB_CT():
    """ Make sure that Divergence of B is zero everywhere. """
    sim = sim_dir_load("MHDCTOrszagTangAMR.enzo",
                       path="./MHD/2D/MHDCTOrszagTangAMR",
                       find_outputs=True)
    sim.get_time_series()
    yield VerifySimulationSameTest(sim)
    # Only test the last output.
    pf = sim[-1]
    max_value = pf.h.find_max('DivB')
    max_value = float(max_value[0])
    assert (max_value < 1e-10)
