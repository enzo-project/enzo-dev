from yt.mods import *

# Using the parameters in the default parameter file
mach = 2.0
gamma = 1.4
sound_speed = na.sqrt(gamma) * (1.0/1.0) # see initializer for info
shock_pool_shock_density = 1.0 * ( ( (gamma + 1.0) * mach * mach) 
                                 / ( (gamma - 1.0) * mach * mach + 2.0))
shock_pool_shock_vel = sound_speed * mach 

for i in xrange(42):
    pf = load("DD%04i/DD%04i" % (i,i))
    pc = PlotCollection(pf, center=[0.5,0.5,0.5])
    pc.add_slice("Density", 2)
    p = pc.add_ray([0.0, 0.0, 0.5], [1.0, 1.0, 0.5], "Density")
    # We calculate the velocity of the shock, assuming mach number of 2.0, but
    # we divide out sqrt(2) to scale to 0..1, because it's across diagonals (45
    # deg angle for the shock)
    current_position = shock_pool_shock_vel * pf["InitialTime"] / na.sqrt(2.0)
    p.modify["line"]( (current_position, current_position),
                      (1.0, shock_pool_shock_density))
    pc.save("frames/%s" % pf)
