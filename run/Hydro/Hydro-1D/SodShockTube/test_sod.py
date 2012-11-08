import os
from yt.mods import *
from yt.funcs import *
from yt.testing import *
from yt.utilities.answer_testing.framework import \
    requires_outputlog, \
    small_patch_amr, \
    big_patch_amr, \
    data_dir_load

_solution_file = 'SodShockTube_t=0.25_exact.txt'
_fields = ['Density','ThermalEnergy']
_les = [0.25, 0.85]
_res = [0.4, 0.9]
_rtol = 1.0e-2
_atol = 1.0e-7

@requires_outputlog(os.path.dirname(__file__), "SodShockTube.enzo") # Verifies that OutputLog exists
def test_sod():
    # Read in the pf
    pf = load('DD0001/data0001')  
    exact = get_analytical_solution(_solution_file) 
   
    ad = pf.h.all_data()
    position = ad['x']
    for k in _fields:
        field = ad[k]
        for xmin, xmax in zip(_les, _res):
            mask = (position >= xmin)*(position <= xmax)
            exact_field = np.interp(position[mask], exact['pos'], exact[k]) 
            # yield test vs analytical solution 
            yield assert_allclose, field[mask], exact_field, _rtol, _atol

def get_analytical_solution(fname):
    # Reads in from file 
    pos, dens, vel, pres, inte = \
            np.loadtxt(fname, unpack=True)
    exact = {}
    exact['pos'] = pos
    exact['Density'] = dens
    exact['x-velocity'] = vel
    exact['Pressure'] = pres
    exact['ThermalEnergy'] = inte
    return exact

