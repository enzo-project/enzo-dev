import yt
import matplotlib.pyplot as plt 
import numpy as np
import os
from yt.testing import assert_rel_equal

sim_dir = os.path.basename(os.getcwd())
test_data_dir = os.path.join(
    os.environ.get("COSMO_TEST_DATA_DIR", "."), sim_dir)
if not os.path.exists(test_data_dir):
    os.makedirs(test_data_dir)
generate_answers = int(os.environ.get("COSMO_TEST_GENERATE",0 ))

def test_phase():
    es = yt.simulation('amr_cosmology.enzo', 'Enzo')
    es.get_time_series(redshifts=[0])
    ds = es[-1]
    ad = ds.all_data()
    profile = ad.profile([("gas", "density")],
                         [("gas", "temperature"),
                          ("gas", "cell_mass")])
    profile1 = ad.profile([("gas", "density")],
                          [("gas", "temperature"),
                           ("gas", "cooling_time")],
                          weight_field=('gas', 'cell_mass'))
    density = profile.x
    temperature = profile['gas', 'temperature']
    cooling_time = profile1['gas', 'cooling_time']
    cell_mass = profile['gas', 'cell_mass']

    filename = 'phase_data.h5'
    data = {'density': density, 'temperature': temperature, 'cooling_time': cooling_time, 'cell_mass': cell_mass}
    yt.save_as_dataset(ds,filename,data)

    pp = yt.PhasePlot(ad, 'density', 'temperature', 'cell_mass')
    pp.set_unit('cell_mass', 'Msun')
    pp.save()
    pp1 = yt.PhasePlot(ad, 'density', 'temperature', 'cooling_time', weight_field='cell_mass')
    pp1.save()   

    compare_filename = os.path.join(test_data_dir, filename)
    if generate_answers:
        os.rename(filename, compare_filename)
        return

        # do the comparison
    ds_comp = yt.load(compare_filename)

 

    # assert quality to 8 decimals
    assert_rel_equal(data['density'], ds_comp.data['density'], 8)
    assert_rel_equal(data['temperature'], ds_comp.data['temperature'], 8)
    assert_rel_equal(data['cooling_time'], ds_comp.data['cooling_time'], 8)
    assert_rel_equal(data['cell_mass'], ds_comp.data['cell_mass'], 8)

