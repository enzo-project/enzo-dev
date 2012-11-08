from yt.mods import *
from yt.funcs import *
from yt.testing import *
from yt.utilities.answer_testing.framework import \
    requires_pf, \
    small_patch_amr, \
    big_patch_amr, \
    data_dir_load

_solution_file = 'Toro-4-ShockTube_t=0.035_exact.txt'
_fields = ['Density','x-velocity','Pressure','ThermalEnergy']
_les = [0.]
_res = [1.]
_rtol = 1.0e-1
_atol = 1.0e-7

def test_shocktube():
    if not os.path.isfile('DD0001/data0001'):
        return
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

