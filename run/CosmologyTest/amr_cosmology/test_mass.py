"""
   Three functions to test (separately) the total dark matter mass,
   the total baryon mass, and the total stellar and total gas mass
   is the test cosmology simulations. Latter two only run if baryons
   are present.
"""
import yt
import numpy as np
from yt.utilities.answer_testing.framework import \
    sim_dir_load, \
    assert_rel_equal
import os, glob

test_data_dir   = os.environ.get("COSMO_TEST_DATA_DIR", None)
compare_answers = int(os.environ.get("COSMO_TEST_COMPARE",0))
tolerance       = os.environ.get("COSMO_TEST_MASS_TOLERANCE",8)

def test_dark_matter_mass():
    pf_name = glob.glob('./*.enzo')[0]

    # gather most recent data set
    sim = sim_dir_load(pf_name, path= './',
                       find_outputs=True)
    sim.get_time_series()
    ds    = sim[-1]
    data  = ds.all_data()

    # sum masses
    MDM   = np.sum(data['particle_mass'][ data['particle_type'] == 1 ].to('Msun'))

    output_data    = {'mass' : MDM}

    # save
    filename = "DM_mass_results.h5"
    yt.save_as_dataset(ds, filename, output_data)

    if compare_answers:
        compare_filename = os.path.join(test_data_dir, filename)
        ds_comp          = yt.load(compare_filename)

        # assert
        assert_rel_equal(output_data['mass'], ds_comp.data['mass'], tolerance)

    return


def test_individual_baryon_mass():
    pf_name = glob.glob('./*.enzo')[0]

    # gather most recent data set
    sim = sim_dir_load(pf_name, path= './',
                       find_outputs=True)

    if (sim.parameters['CosmologySimulationOmegaBaryonNow'] == 0.0):
        return

    sim.get_time_series()
    ds    = sim[-1]
    data  = ds.all_data()

    # sum masses
    Mstar = np.sum(data['particle_mass'][ data['particle_type'] == 2 ].to('Msun'))
    Mgas  = np.sum(data['cell_mass'].to('Msun'))

    output_data    = {'masses' : np.array([Mstar, Mgas])}

    # save
    filename = "gas_stars_mass_results.h5"
    yt.save_as_dataset(ds, filename, output_data)

    if compare_answers:
        compare_filename = os.path.join(test_data_dir, filename)
        ds_comp          = yt.load(compare_filename)

        # assert
        assert_rel_equal(output_data['masses'], ds_comp.data['masses'], tolerance)

    return

def test_total_baryon_mass():
    pf_name = glob.glob('./*.enzo')[0]

    # gather most recent data set
    sim = sim_dir_load(pf_name, path= './',
                       find_outputs=True)

    if (sim.parameters['CosmologySimulationOmegaBaryonNow'] == 0.0):
        return

    sim.get_time_series()
    ds    = sim[-1]
    data  = ds.all_data()

    # sum masses
    Mstar = np.sum(data['particle_mass'][ data['particle_type'] == 2 ].to('Msun'))
    Mgas  = np.sum(data['cell_mass'].to('Msun'))

    output_data    = {'masses' : Mstar + Mgas}

    # save
    filename = "baryon_mass_results.h5"
    yt.save_as_dataset(ds, filename, output_data)

    if compare_answers:
        compare_filename = os.path.join(test_data_dir, filename)
        ds_comp          = yt.load(compare_filename)

        # assert
        assert_rel_equal(output_data['masses'], ds_comp.data['masses'], tolerance)

    return

