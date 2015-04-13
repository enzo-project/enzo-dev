# find_ellipsoidal_refinement_mask.py
# Christine Simpson; April 2015
#
# This is a yt script that provides a worked example for how to select
# a ellipsoidal masked refinement region given a halo within a unigrid
# cosmological simulation.

final_data_dump = "DD0127/DD0127" 
initial_data_dump = "DD0000/DD0000"
selected_halo_index = 10

import explore_halo_modules as eh
from yt.mods import *
from yt.analysis_modules.halo_finding.api import *
import pylab

pf_final = load(final_data_dump)
halo_list = HaloFinder(pf_final) # This line can be replaced with a read call to a stored halo catalog

pf_initial = load(initial_data_dump)

halo_list.dump("HaloCatalog")

halo = halo_list[selected_halo_index]

halo_id_final = halo["particle_index"]
halo_radius = halo.virial_radius()
halo_pos = halo.maximum_density_location()

dd = pf_initial.h.all_data()

halo_initial_index = eh.match_halo_ids(halo_id_final,dd["particle_index"])
halo_initial_pos = np.array(zip(dd['particle_position_x'],
                                dd['particle_position_y'],
                                dd['particle_position_z']))[halo_initial_index]

# Here we print out a region point file
f = open('RegionPointFile'+str(selected_halo_index)+'.dat','w')
for i in range(len(halo_initial_pos)):
    line =  str(halo_initial_pos[i][0])+' '+str(halo_initial_pos[i][1])+' '+str(halo_initial_pos[i][2])+'\n'
    f.write(line)
f.close()

# If you want to specify an ellipsoidal matrix in the MUSIC config file, the 
# simplest way to do this is to run MUSIC with the region point file first (set 
# region = ellipsoid).  Once MUSIC reads the point file, it will spit out the 
# ellipsoid center and matrix (called xc and A).  These can be used to set 
# region_ellipsoid_center and region_ellipsoid_matrix.  Running MUSIC in either
# of these modes is equivalent, but it may be more convenient in some cases
# to specify the ellipsoidal matrix, depending on the situation.

