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

size_1d = 50 # square table; length of each side
z = 0.27 # solar units

dens = np.logspace(-29, -23, size_1d) # cgs
temp = np.logspace(4, 6, size_1d) # K
d, t = np.meshgrid(dens, temp)

current_redshift = 0.

# Set solver parameters
chem = chemistry_data()
chem.use_grackle = 1
chem.with_radiative_cooling = 0
chem.primordial_chemistry = 2
chem.metal_cooling = 1
chem.UVbackground = 1
chem.cmb_temperature_floor = 1
chem.h2_on_dust = 1
chem.grackle_data_file = "/home/claire/grackle-git/input/CloudyData_UVB=HM2012.h5"

chem.use_specific_heating_rate = 0
chem.use_volumetric_heating_rate = 0

# Set units
chem.comoving_coordinates = 0 # proper units
chem.a_units = 1.0
chem.a_value = 1.0 / (1.0 + current_redshift) / chem.a_units
chem.density_units = 1.67e-27
chem.length_units = 1638.4 * cm_per_kpc
chem.time_units = sec_per_Myr
chem.velocity_units = chem.a_units \
    * (chem.length_units / chem.a_value) \
    / chem.time_units

# Call convenience function for setting up a fluid container.
# This container holds the solver parameters, units, and fields.
fc = setup_fluid_container(chem,
                           density=d.flatten(),
                           temperature=t.flatten(),
                           metal_mass_fraction=0.02014*z,
                           converge=True)

# Save in cgs
dataset = {'HI'   :fc['HI']   *fc.chemistry_data.density_units,
           'HII'  :fc['HII']  *fc.chemistry_data.density_units,
           'HeI'  :fc['HeI']  *fc.chemistry_data.density_units,
           'HeII' :fc['HeII'] *fc.chemistry_data.density_units,
           'HeIII':fc['HeIII']*fc.chemistry_data.density_units,
           'HM'   :fc['HM']   *fc.chemistry_data.density_units,
           'H2I'  :fc['H2I']  *fc.chemistry_data.density_units,
           'H2II' :fc['H2II'] *fc.chemistry_data.density_units,
           'de'   :fc['de']   *fc.chemistry_data.density_units,
           'metal':fc['metal']*fc.chemistry_data.density_units,
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

yt.save_as_dataset(None, f'equilibrium_table_{size_1d}_0{z*100:.0f}-Zsun.h5', dataset,
                   field_types = field_types)
