import glob
import yt
import matplotlib.pyplot as plt 
import numpy as np
import os
from yt_astro_analysis.halo_analysis import HaloCatalog
from yt.testing import assert_rel_equal
from numpy.testing import assert_equal

from yt.utilities.answer_testing.framework import \
     sim_dir_load

_pf_name = os.path.basename(os.path.dirname(__file__)) + ".enzo"
_dir_name = os.path.dirname(__file__)

sim_dir = os.path.basename(_dir_name)
test_data_dir = os.path.join(
    os.environ.get("COSMO_TEST_DATA_DIR", "."), sim_dir)
if not os.path.exists(test_data_dir):
    os.makedirs(test_data_dir)
generate_answers = int(os.environ.get("COSMO_TEST_GENERATE", 1))
tolerance        = os.environ.get("COSMO_TEST_MASS_TOLERANCE", 8)

def test_hmf():
    es = sim_dir_load(_pf_name, path=_dir_name)
    es.get_time_series()
    ds = es[-1]
    hc = HaloCatalog(
        data_ds=ds, finder_method='fof',
        output_dir=os.path.join(_dir_name, "halo_catalogs/catalog"))
    hc.create()
    masses = hc.data_ds.r[('all', 'particle_mass')].in_units('Msun')  
    h = ds.hubble_constant
    mtot = np.log10(masses*1.2) - np.log10(h)
    masses_sim = np.sort(mtot)
    sim_volume = ds.domain_width.in_units('Mpccm').prod()
    n_cumulative_sim = np.arange(len(mtot),0,-1)
    masses_sim,unique_indices = np.unique(masses_sim,return_index=True)
    
    n_cumulative_sim = n_cumulative_sim[unique_indices]/sim_volume
    filename = 'hmf.h5'
    save_filename = os.path.join(_dir_name, filename)
    data = {'masses': masses_sim, 'n_sim': n_cumulative_sim}
    yt.save_as_dataset(ds, save_filename, data)

    # make a plot
    fig = plt.figure(figsize=(8,8))
    plt.semilogy(masses_sim,n_cumulative_sim,'-')
    plt.ylabel('Cumulative Halo Number Density $\mathrm{Mpc}^{-3}$',fontsize=16)
    plt.xlabel('log Mass/$\mathrm{M}_{\odot}$',fontsize=16)
    plt.tick_params(labelsize=16)
    plt.savefig(os.path.join(_dir_name, 'hmf.png'),format='png')

    compare_filename = os.path.join(test_data_dir, filename)
    if generate_answers:
        os.rename(save_filename, compare_filename)
        return

    # do the comparison
    ds_comp = yt.load(compare_filename)

    # assert quality to 8 decimals
    assert_rel_equal(data['masses'], ds_comp.data['masses'], 8)
    assert_rel_equal(data['n_sim'], ds_comp.data['n_sim'], 8)

def test_dark_matter_mass():
    # gather most recent data set
    sim = sim_dir_load(_pf_name, path=_dir_name,
                       find_outputs=True)
    sim.get_time_series()
    ds    = sim[-1]
    data  = ds.all_data()

    # sum masses
    MDM   = np.sum(data[('all', 'particle_mass')][ data[('all', 'particle_type')] == 1 ].to('Msun'))

    output_data    = {('data', 'mass') : MDM}

    # save
    filename = "DM_mass_results.h5"
    save_filename = os.path.join(_dir_name, filename)
    yt.save_as_dataset(ds, save_filename, output_data)

    compare_filename = os.path.join(test_data_dir, filename)
    if generate_answers:
        os.rename(save_filename, compare_filename)
        return

    ds_comp          = yt.load(compare_filename)
    assert_rel_equal(output_data[('data', 'mass')], ds_comp.data[('data', 'mass')], tolerance)

def test_output_number():
    ds = yt.load(os.path.join(_dir_name, 'DD0000/DD0000'))

    DDnum = len(glob.glob(os.path.join(_dir_name, 'DD????/DD????')))
    RDnum = len(glob.glob(os.path.join(_dir_name, 'RD????/RD????')))

    output_data = {'number_of_files' : np.array([DDnum,RDnum])}

    filename = "outputnum_results.h5"
    save_filename = os.path.join(_dir_name, filename)
    yt.save_as_dataset(ds, save_filename, output_data)

    compare_filename = os.path.join(test_data_dir, filename)
    if generate_answers:
        os.rename(save_filename, compare_filename)
        return

    ds_comp = yt.load(compare_filename)
    assert_equal(output_data['number_of_files'],
                 ds_comp.data['number_of_files'])
