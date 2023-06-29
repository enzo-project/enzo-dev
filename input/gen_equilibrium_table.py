###############################################################################
#
# Use Grackle's python wrapper to make a table of equilibrium densities for 
# 9 chemical species on a square grid of density & temperature (with constant
# metallicity). This grid can be used to initialize the isolated galaxy
# problem type so that chemical species are closer to their equilibrium values.
#
# As is, this script makes a table spanning 1e-29 to 1e-23 g/cm^3
# and 1e4 to 1e6 K using 50 cells on a side. The constant metallicity is 
# 0.27 Zsun and the HM2012 UV background is used. The table is calculated for
# redshift 0. The table size and metallicity are included in the resulting HDF5 filename.
#
###############################################################################

import numpy as np
import yt

from pygrackle import \
    chemistry_data, \
    setup_fluid_container
from pygrackle.utilities.physical_constants import \
	cm_per_mpc as Mpc, \
	amu_cgs as m_u, \
	mass_sun_cgs as Msun, \
	keV_per_K as kboltz, \
	sec_per_Myr, \
	cm_per_kpc

# 1e-29 to 1e-23 g/cm**3
# 6.5e4 to 1.5e6 K

size_1d = 60 # square table; length of each side
Z = 0.3 # solar units

dens = np.logspace(np.log10(1e-31), np.log10(1e-23), size_1d) # cgs
temp = np.logspace(np.log10(3e3), np.log10(3e6), size_1d) # K
d, t = np.meshgrid(dens, temp)

# Set solver parameters
chem = chemistry_data()
chem.use_grackle = 1
chem.with_radiative_cooling = 0
chem.primordial_chemistry = 2
chem.metal_cooling = 1
chem.UVbackground = 1
chem.cmb_temperature_floor = 1
chem.grackle_data_file = b"/PATH/TO/GRACKLE/input/CloudyData_UVB=HM2012.h5"

chem.use_specific_heating_rate = 0
chem.use_volumetric_heating_rate = 0

# Set units
chem.comoving_coordinates = 0 # proper units
chem.a_units = 1.0
chem.a_value = 1.0
chem.density_units = 1.67e-27
chem.length_units = 1638.4 * cm_per_kpc
chem.time_units = sec_per_Myr
chem.velocity_units = chem.a_units \
    * (chem.length_units / chem.a_value) \
    / chem.time_units

# Call convenience function for setting up a fluid container.
# This container holds the solver parameters, units, and fields.
fc = setup_fluid_container(chem,
                           density=d.ravel(),
                           temperature=t.ravel(),
                           metal_mass_fraction=0.01295*Z,
                           converge=True,
                           tolerance=1e-7, # default: 0.01
                           max_iterations=10000000000)

# Save as fraction
dataset = {'HI'   :fc['HI']   /fc['density'],
           'HII'  :fc['HII']  /fc['density'],
           'HeI'  :fc['HeI']  /fc['density'],
           'HeII' :fc['HeII'] /fc['density'],
           'HeIII':fc['HeIII']/fc['density'],
           'HM'   :fc['HM']   /fc['density'],
           'H2I'  :fc['H2I']  /fc['density'],
           'H2II' :fc['H2II'] /fc['density'],
           'de'   :fc['de']   /fc['density'],
           'metal':fc['metal']/fc['density'],
           'density': dens,
           'temperature': temp}

field_types = {'HI'   :'table',
               'HII'  :'table',
               'HeI'  :'table',
               'HeII' :'table',
               'HeIII':'table',
               'HM'   :'table',
               'H2I'  :'table',
               'H2II' :'table',
               'de'   :'table',
               'metal':'table',
               'density': 'indexer',
               'temperature': 'indexer'}

yt.save_as_dataset(None, f'equilibrium_table_{size_1d}_0{Z*100:.0f}-Zsun.h5', dataset,
                   field_types = field_types)
