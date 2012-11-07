from yt.mods import *
from yt.funcs import *
from yt.testing import *
from yt.utilities.answer_testing.framework import \
    requires_pf, \
    small_patch_amr, \
    big_patch_amr, \
    data_dir_load

def get_analytical_solution():
    # Reads in from file 
    return np.loadtxt('Toro-3-ShockTube_t=0.012_exact.txt', unpack=True)

def test_toro2():
    if not os.path.isfile('DD0001/data0001'):
        return
    # Read in the pf
    pf = load('DD0001/data0001')  
    pos, dens, vel, pres, inte = get_analytical_solution() 
    exact = {}
    exact['pos'] = pos
    exact['Density'] = dens
    exact['x-velocity'] = vel
    exact['Pressure'] = pres
    exact['ThermalEnergy'] = inte
   
    ad = pf.h.all_data()
    calc_position = ad['x']
    #for k in ['Density','x-velocity','Pressure','ThermalEnergy']:
    for k in ['Density']:
        calc_field = ad[k]
        for xmin, xmax in zip([0.1, 0.6], [0.25, 0.7]):
            mask = (calc_position >= xmin)*(calc_position <= xmax)
            exact_field = np.interp(calc_position[mask], exact['pos'], exact[k]) 
            # yield Test vs analytical solution (assert_relative_equal)
            yield assert_rel_equal, calc_field[mask], exact_field, 2
