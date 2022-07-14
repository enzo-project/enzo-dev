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

def test_max_density_halo_quantities():
    ds = yt.load(os.path.join(_dir_name, 'RD0009/RD0009'))

    # Find the point of maximum density, center a sphere of radius
    # 1 Mpc around it, and sum the masses inside
    val,pos = ds.find_max(('enzo', 'Density'))
    sp = ds.sphere(pos,(1000.,'kpc'))
    ct = sp[('nbody', 'creation_time')]
    dm = (ct < 0)
    dm_mass = np.sum(sp[('all', 'particle_mass')][dm]).in_units('Msun')
    gas_mass = np.sum(sp[('gas', 'cell_mass')].in_units('Msun'))

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
    save_filename = os.path.join(_dir_name, filename)
    yt.save_as_dataset(ds, save_filename, data)

    compare_filename = os.path.join(test_data_dir, filename)
    if generate_answers:
        os.rename(save_filename, compare_filename)
        return

    ds_comp          = yt.load(compare_filename)
    assert_rel_equal(data["dm_mass"], ds_comp.data["dm_mass"], tolerance)
    assert_rel_equal(data["gas_mass"], ds_comp.data["gas_mass"], tolerance)
    assert_rel_equal(data["max_position"], ds_comp.data["max_position"], tolerance)
    assert_rel_equal(data["density_profile"], ds_comp.data["density_profile"], tolerance)
    assert_rel_equal(data["temperature_profile"], ds_comp.data["temperature_profile"], tolerance)

def test_dark_matter_mass():
    # gather most recent data set
    sim = sim_dir_load(_pf_name, path=_dir_name,
                       find_outputs=True)
    sim.get_time_series()
    ds    = sim[-1]
    data  = ds.all_data()

    # sum masses
    MDM   = np.sum(data[('all', 'particle_mass')][ data[('all', 'particle_type')] == 1 ].to('Msun'))

    output_data    = {(('data', 'mass') : MDM}

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

def test_individual_baryon_mass():
    # gather most recent data set
    sim = sim_dir_load(_pf_name, path=_dir_name,
                       find_outputs=True)

    if (sim.parameters['CosmologySimulationOmegaBaryonNow'] == 0.0):
        return

    sim.get_time_series()
    ds    = sim[-1]
    data  = ds.all_data()

    # sum masses
    Mstar = np.sum(data[('all', 'particle_mass')][ data[('all', 'particle_type')] == 2 ].to('Msun'))
    Mgas  = np.sum(data[('gas', 'cell_mass')].to('Msun'))

    output_data    = {'masses' : np.array([Mstar, Mgas])}

    # save
    filename = "gas_stars_mass_results.h5"
    save_filename = os.path.join(_dir_name, filename)
    yt.save_as_dataset(ds, save_filename, output_data)

    compare_filename = os.path.join(test_data_dir, filename)
    if generate_answers:
        os.rename(save_filename, compare_filename)
        return

    ds_comp          = yt.load(compare_filename)
    assert_rel_equal(output_data['masses'], ds_comp.data['masses'], tolerance)

def test_total_baryon_mass():
    # gather most recent data set
    sim = sim_dir_load(_pf_name, path=_dir_name,
                       find_outputs=True)

    if (sim.parameters['CosmologySimulationOmegaBaryonNow'] == 0.0):
        return

    sim.get_time_series()
    ds    = sim[-1]
    data  = ds.all_data()

    # sum masses
    Mstar = np.sum(data[('all', 'particle_mass')][ data[('all', 'particle_type')] == 2 ].to('Msun'))
    Mgas  = np.sum(data[('gas', 'cell_mass')].to('Msun'))

    output_data    = {'masses' : Mstar + Mgas}

    # save
    filename = "baryon_mass_results.h5"
    save_filename = os.path.join(_dir_name, filename)
    yt.save_as_dataset(ds, save_filename, output_data)

    compare_filename = os.path.join(test_data_dir, filename)
    if generate_answers:
        os.rename(save_filename, compare_filename)
        return

    ds_comp          = yt.load(compare_filename)
    assert_rel_equal(output_data['masses'], ds_comp.data['masses'], tolerance)

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

def test_phase():
    es = sim_dir_load(_pf_name, path=_dir_name)
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
    temperature = profile[('gas', 'temperature')]
    cooling_time = profile1[('gas', 'cooling_time')]
    cell_mass = profile[('gas', 'cell_mass')]

    filename = 'phase_data.h5'
    save_filename = os.path.join(_dir_name, filename)
    data = {('data', 'density'): density, ('data', 'temperature'): temperature,
            ('data', 'cooling_time'): cooling_time, ('data', 'cell_mass'): cell_mass}
    yt.save_as_dataset(ds, save_filename, data)

    pp = yt.PhasePlot(ad, ('gas', 'density'), ('gas', 'temperature'), ('gas', 'cell_mass'))
    pp.set_unit(('gas', 'cell_mass'), 'Msun')
    pp.save(_dir_name)
    pp1 = yt.PhasePlot(ad, ('gas', 'density'), ('gas', 'temperature'), ('gas', 'cooling_time'),
                       weight_field=('gas', 'cell_mass'))
    pp1.save(_dir_name)

    compare_filename = os.path.join(test_data_dir, filename)
    if generate_answers:
        os.rename(save_filename, compare_filename)
        return

        # do the comparison
    ds_comp = yt.load(compare_filename)

    # assert quality to 8 decimals
    assert_rel_equal(data[('data', 'density')], ds_comp.data[('data', 'density')], 8)
    assert_rel_equal(data[('data', 'temperature')], ds_comp.data[('data', 'temperature')], 8)
    assert_rel_equal(data[('data', 'cooling_time')], ds_comp.data[('data', 'cooling_time')], 8)
    assert_rel_equal(data[('data', 'cell_mass')], ds_comp.data[('data', 'cell_mass')], 8)
