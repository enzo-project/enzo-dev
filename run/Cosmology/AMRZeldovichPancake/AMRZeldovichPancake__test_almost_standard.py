###
### This is a testing template
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
                "y-velocity",
                "z-velocity",
                "Gas_Energy",
                "particle_position_x",
                "particle_position_y",
                "particle_position_z")

@requires_outputlog(os.path.dirname(__file__), 
                    "AMRZeldovichPancake.enzo") # Verifies that OutputLog exists
def test_standard():
    sim = sim_dir_load("AMRZeldovichPancake.enzo",
                       path="./Cosmology/AMRZeldovichPancake",
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
