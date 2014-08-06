# YT script that will output the following by columns:
#
# 1) years after big bang
# 2) total emissivity contributing to feedback (actually emissivity density) values in the box in units erg/s/cm^3
# 3) redshift
# 4) number of cells with feedback sources
# 5) total mass of all star particles
# 6) total number of stars
#
# REMEMBER TO RENAME THE OLD STARFORMATIONRATE.TXT BEFORE RUNNING
#
#################################################################################################################

from yt.mods import *

def _ParticleAge(field, data):
   current_time = data.pf["InitialTime"]
   return (current_time - data["creation_time"])
def _convertParticleAge(data):
   return data.convert("years")
add_field("ParticleAge",
          function=_ParticleAge, particle_type=True, convert_function=_convertParticleAge)

def _StarParticleMass(field, data):
    return data["ParticleMassMsun"] * (data["creation_time"]>0)
add_field("StarParticleMass", function=_StarParticleMass, particle_type=True)

def _ZeroEmissivity(field, data):
    val = (data["Emissivity"] ==0.0).astype("float64")
    return val
add_field("ZeroEmissivity", function=_ZeroEmissivity)

def _TotEmissivity(field, data):
    val = data["Emissivity"]
    return val
add_field("TotEmissivity", function=_TotEmissivity)

def _StarParticles(field, data):
    return (data["creation_time"]>0)
add_field("StarParticles", function=_StarParticles, particle_type=True)


min_output_number = 180
max_output_number = 273
skip = 1
resolution = 64

# This will append the values into the file, so be sure to delete the old values before running script
def flush(fn, towrite):
   f = open(fn, "a")
   f.write(towrite)
   f.close()

for i in range(min_output_number, max_output_number+1, skip):

    pf = load("DD%04i/DD%04i" % (i,i))

    print ("loaded DD%04i" % i)

    dd = pf.h.all_data()

    redshift = pf.get_parameter("CosmologyCurrentRedshift")

    num_star = dd.quantities["TotalQuantity"]("StarParticles")[0]
    TotEmis = 0
    num_sources = 0

    z = pf["CosmologyCurrentRedshift"]
    yr = pf["InitialTime"] * pf["years"]

    if num_star == 0:
        yt.funcs.only_on_root(flush, "StarFormationRate.txt", "%12.12e %12.12e %12.12e %12.12i %12.12e %12.12i\n" % (yr, TotEmis, redshift, num_sources, 0, num_star))
    else:
        num_cells_zero = dd.quantities["TotalQuantity"]("ZeroEmissivity")[0]
        TotEmis = dd.quantities["TotalQuantity"]("TotEmissivity")[0]
        num_sources = resolution**3-num_cells_zero

        totalmass = dd.quantities["TotalQuantity"]("StarParticleMass")[0]

        yt.funcs.only_on_root(flush, "StarFormationRate.txt", "%12.12e %12.12e %12.12e %12.12i %12.12e %12.12i\n" % (yr, TotEmis, redshift, num_sources, totalmass, num_star))
