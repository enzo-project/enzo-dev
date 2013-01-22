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
                "Pressure",
                "x-velocity",
                "y-velocity",
                "Bx",
                "By")

_pf_name = os.path.basename(os.path.dirname(__file__)) + ".enzo"
_dir_name = os.path.dirname(__file__)

@requires_outputlog(_dir_name, _pf_name)
def test_standard():
    sim = sim_dir_load("MHDCTOrszagTang.enzo",
                       path="./MHD/2D/MHDCTOrszagTang",
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

# Defines a new derived field for divergence of B
def _DivB_CT(field,data):
    """ Doesn't seem to work right now..."""
    bx = data['BxF']
    by = data['ByF']
    bz = data['BzF']
    #print d.shape
    nx,ny,nz=bx.shape
    #print "ManualDivB: bx shape", bx.shape
    nx-=1 #because it came from bx, which is one too large
    return ((bx[1:nx+1,:,:]-bx[0:nx,:,:])/data['dx']+
            (by[:,1:ny+1,:]-by[:,0:ny,:])/data['dy']+
            (bz[:,:,1:nz+1]-bz[:,:,0:nz])/data['dz'])

# Tests that Div B = 0
@requires_outputlog(_dir_name, _pf_name)
def test_DivB_CT():
    """ Make sure that Divergence of B is zero everywhere. """
    add_field('DivB_CT', function = _DivB_CT,take_log=False,
              validators=[ValidateGridType()])
    sim = sim_dir_load("MHDCTOrszagTang.enzo",
                       path="./MHD/2D/MHDCTOrszagTang",
                       find_outputs=True)
    sim.get_time_series()
    yield VerifySimulationSameTest(sim)
    # Only test the last output.
    pf = sim[-1]
    max_value = pf.h.find_max('DivB_CT')
    max_value = float(max_value[0])
    assert (max_value < 1e-10)
