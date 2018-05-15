import numpy as np
import os, glob
import yt
from yt.utilities.answer_testing.framework import \
    sim_dir_load, \
    assert_rel_equal

sim_dir = os.path.basename(os.getcwd())
test_data_dir = os.path.join(
    os.environ.get("COSMO_TEST_DATA_DIR", None), sim_dir)
if not os.path.exists(test_data_dir):
    os.makedirs(test_data_dir)
generate_answers = int(os.environ.get("COSMO_TEST_GENERATE",1))
tolerance       = os.environ.get("COSMO_TEST_MASS_TOLERANCE",8)

def test_max_density_halo_quantities():
    ds = yt.load('RD0009/RD0009')

    # Find the point of maximum density, center a sphere of radius
    # 1 Mpc around it, and sum the masses inside
    val,pos = ds.find_max('Density')
    sp = ds.sphere(pos,(1000.,'kpc'))
    ct = sp['creation_time']
    dm = (ct < 0)
    dm_mass = np.sum(sp['particle_mass'][dm]).in_units('Msun')
    gas_mass = np.sum(sp['cell_mass'].in_units('Msun'))

    # Also look at the radial profiles of density and temperature
    # within these spheres. The bin size is chosen to make the profiles
    # smooth and for each bin to be larger than the cell size.
    ptest0 = yt.create_profile(sp,"radius","density",n_bins=[20])
    ptest1 = yt.create_profile(sp,"radius","temperature",n_bins=[20])

    # Save the quantities to be compared
    data = {"dm_mass": dm_mass,
            "gas_mass": gas_mass,
            "max_position": pos,
            "density_profile": ptest0['density'],
            "temperature_profile": ptest1['temperature']}

    # save your results file
    filename = "max_density_halo_quantities.h5"
    yt.save_as_dataset(ds, filename, data)

    compare_filename = os.path.join(test_data_dir, filename)
    if generate_answers:
        os.rename(filename, compare_filename)
        return

    ds_comp          = yt.load(compare_filename)
    assert_rel_equal(data["dm_mass"], ds_comp.data["dm_mass"], tolerance)
    assert_rel_equal(data["gas_mass"], ds_comp.data["gas_mass"], tolerance)
    assert_rel_equal(data["max_position"], ds_comp.data["max_position"], tolerance)
    assert_rel_equal(data["density_profile"], ds_comp.data["density_profile"], tolerance)
    assert_rel_equal(data["temperature_profile"], ds_comp.data["temperature_profile"], tolerance)
    
    return
