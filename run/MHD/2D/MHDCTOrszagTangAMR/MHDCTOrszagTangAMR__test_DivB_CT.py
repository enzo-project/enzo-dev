###
### MHDCT Test
###

import os
import numpy as np
from yt.testing import *
from yt.utilities.answer_testing.framework import \
    VerifySimulationSameTest, \
    sim_dir_load
from yt.frontends.enzo.answer_testing_support import \
    requires_outputlog, \
    standard_small_simulation

_base_fields = ("density",
                "pressure",
                "velocity_x",
                "velocity_y",
                "magnetic_field_x",
                "magnetic_field_y")

@requires_outputlog(os.path.dirname(__file__), "MHDCTOrszagTangAMR.enzo") # Verifies that OutputLog exists
def test_standard():
    sim = sim_dir_load("MHDCTOrszagTangAMR.enzo",
                       path="./MHD/2D/MHDCTOrszagTangAMR",
                       find_outputs=True)
    sim.get_time_series()
    yield VerifySimulationSameTest(sim)
    base_pf = sim[0]
    fields = [f for f in _base_fields if ("enzo", f) in base_pf.field_list]
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
    max_value = pf.find_max(('enzo', 'DivB'))
    max_value = float(max_value[0])
    assert (max_value < 1e-10)
