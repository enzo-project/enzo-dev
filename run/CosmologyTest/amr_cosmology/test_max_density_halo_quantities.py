import numpy as np
import os
import yt
from yt.testing import \
    assert_rel_equal

test_data_dir = os.environ.get("COSMO_TEST_DATA_DIR", None)
compare_answers = int(os.environ.get("COSMO_TEST_COMPARE", 0))

# name your function test_<something>
def test_max_density_halo_quantities():
    ds = yt.load('RD0009/RD0009')

    val,pos = ds.find_max('Density')
    sp = ds.sphere(pos,(1000.,'kpc'))
    ct = sp['creation_time']
    dm = (ct < 0)
    dm_mass = np.sum(sp['particle_mass'][dm]).in_units('Msun')
    gas_mass = np.sum(sp['cell_mass'].in_units('Msun'))

    ptest0 = yt.create_profile(sp,"radius","density",n_bins=[20])
    ptest1 = yt.create_profile(sp,"radius","temperature",n_bins=[20])

    # example data
    data = {"dm_mass": dm_mass,
            "gas_mass": gas_mass,
            "max_position": pos,
            "density_profile": ptest0['density'],
            "temperature_profile": ptest1['temperature']}

    # save your results file
    filename = "max_density_halo_quantities.h5"
    yt.save_as_dataset(ds, filename, data)

    if compare_answers:
        compare_filename = os.path.join(test_data_dir, filename)
        ds_comp = yt.load(compare_filename)

        # assert quality to 8 decimals
        assert_rel_equal(data["dm_mass"], ds_comp.data["dm_mass"], 8)
        assert_rel_equal(data["gas_mass"], ds_comp.data["gas_mass"], 8)
        assert_rel_equal(data["max_position"], ds_comp.data["max_position"], 8)
        assert_rel_equal(data["density_profile"], ds_comp.data["density_profile"],8)
        assert_rel_equal(data["temperature_profile"], ds_comp.data["temperature_profile"],8)
