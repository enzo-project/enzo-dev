import yt
import matplotlib.pyplot as plt 
import numpy as np
import os
from yt.analysis_modules.halo_mass_function.api import *
from yt.analysis_modules.halo_analysis.api import HaloCatalog
from yt.testing import assert_rel_equal

test_data_dir = os.environ.get("COSMO_TEST_DATA_DIR", None)
compare_answers = int(os.environ.get("COSMO_TEST_COMPARE", 0))

def test_hmf():
    es = yt.simulation('amr_cosmology.enzo','Enzo')
    es.get_time_series()
    ds = es[-1]
    hc = HaloCatalog(data_ds=ds,finder_method='fof',output_dir="halo_catalogs/catalog")
    hc.create()
    masses = hc.data_source['particle_mass'].in_units('Msun')  
    h = ds.hubble_constant
    mtot = np.log10(masses*1.2) - np.log10(h)
    masses_sim = np.sort(mtot)
    sim_volume = ds.domain_width.in_units('Mpccm').prod()
    n_cumulative_sim = np.arange(len(mtot),0,-1)
    masses_sim,unique_indices = np.unique(masses_sim,return_index=True)
    
    n_cumulative_sim = n_cumulative_sim[unique_indices]/sim_volume
    filename = 'data.h5' 
    data = {'masses': masses_sim, 'n_sim': n_cumulative_sim}
    yt.save_as_dataset(ds,filename,data)
    if compare_answers:
        compare_filename = os.path.join(test_data_dir, filename)
        ds_comp = yt.load(compare_filename)

        # assert quality to 8 decimals
        asswert_rel_equal(data['masses'], ds_comp.data['masses'], 8)
        asswert_rel_equal(data['n_sim'], ds_comp.data['n_sim'], 8)
    ds = yt.load('data.h5')
    masses = ds.data['masses']
    num = ds.data['n_sim']
    fig = plt.figure(figsize=(8,8))
    plt.semilogy(masses_sim,n_cumulative_sim,'-')
    plt.ylabel('Cumulative Halo Number Density $\mathrm{Mpc}^{-3}$',fontsize=16)
    plt.xlabel('log Mass/$\mathrm{M}_{\odot}$',fontsize=16)
    plt.tick_params(labelsize=16)
    plt.savefig('hmf.png',format='png')

