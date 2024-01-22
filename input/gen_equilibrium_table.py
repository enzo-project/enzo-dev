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

import argparse

parser = argparse.ArgumentParser(description="Generate a square table of ionization fractions at equilibrium "
                                 "for a range of densities and temperatures. The assumed metallicity of the table is fixed. "
                                 "Users must specify a UV background. Included species are HI-II, HeI-III, "
                                 "HM (H minus), H2I-II, electrons, and metals to align with Grackle.",
                                 formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument(dest="uvb_file",
                    help="Path to your desired Grackle UVB file")
parser.add_argument("-s","--size", type=int, dest='size_1d', default=60,
                    help="Length of the square table in both dimensions")
parser.add_argument("-Z","--metallicity", type=float, dest='met', default=0.3,
                    help="Constant metallicity assumed for the table")
parser.add_argument("-d","--dens-bounds", nargs=2, type=float, dest='dens',
                    default=[1e-31,1e-23], help="Density bounds for the table in g/cc")
parser.add_argument("-t","--temp-bounds", nargs=2, type=float, dest='temp',
                    default=[3e3,3e6], help="Temperature bounds for the table in K")
parser.add_argument("--max-iter", type=int, dest="max_iter", default=1e10,
                    help="Maximum number of iterations for the table to converge")
parser.add_argument("--tol", type=float, default=1e-7,
                    help="Acceptable relative error between iterations")

args = parser.parse_args()
size_1d = args.size_1d # square table; length of each side
Z = args.met # solar units

dens = np.logspace(np.log10(args.dens[0]), np.log10(args.dens[1]), size_1d) # cgs
temp = np.logspace(np.log10(args.temp[0]), np.log10(args.temp[1]), size_1d) # K
d, t = np.meshgrid(dens, temp)

# Set solver parameters
chem = chemistry_data()
chem.use_grackle = 1
chem.with_radiative_cooling = 0
chem.primordial_chemistry = 2
chem.metal_cooling = 1
chem.UVbackground = 1
chem.cmb_temperature_floor = 1
chem.grackle_data_file = args.uvb_file

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
                           tolerance=args.tol, # default: 0.01
                           max_iterations=args.max_iter)

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
