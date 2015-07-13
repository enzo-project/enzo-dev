#!/usr/bin/python
import sys
import os

kboltz = 1.3806504e-16  # erg K^-1
sec_per_Gyr  = 31.5576e15
mpc_per_cm    = 3.24077929e-25
cm_per_mpc    = 1.0 / mpc_per_cm
mh = 1.674534e-24  # g
mu = 0.6

def setup_shockbox(length_units, time_units, d1, T1, m, gamma=None):
    if gamma is None: gamma = 5./3.

    dens_units = mh
    temp_units = dens_units*(length_units/time_units)**2/kboltz
    p1 = (T1/temp_units)*d1/mu
    d2 = d1*((gamma+1.)*m*m)/((gamma-1.)*m*m + 2.);
    p2 = p1*(2.0*gamma*m*m - (gamma-1.))/(gamma+1.);
    c1 = (gamma*p1/(d1))**0.5;
    v1 = 0.0
    v2 = m*c1*(1.-d1/d2);
    shockspeed = 1.0 * c1 * m;

    vel_units = length_units/time_units

    lines = []
    lines.append('\n')
    lines.append('DensityUnits = %0.8e\n' % dens_units)
    lines.append('LengthUnits = %0.8e\n' % length_units)
    lines.append('TimeUnits = %0.8e\n' % time_units)
    lines.append('Gamma = %f\n' % gamma)
    lines.append('ShockInABoxLeftDensity = %0.8e\n' % (d1)) 
    lines.append('ShockInABoxLeftVelocity = %0.8e\n' % ((shockspeed - v1))) 
    lines.append('ShockInABoxLeftPressure = %0.8e\n' % (p1) ) 

    lines.append('ShockInABoxRightDensity = %0.8e\n' % (d2)) 
    lines.append('ShockInABoxRightVelocity = %0.8e\n' % ((shockspeed - v2))) 
    lines.append('ShockInABoxRightPressure = %0.8e\n' % (p2) ) 
    lines.append('\n')

    return lines
 
def get_header(d1, T1, T2, M):
    lines = []
    lines.append('# \n')
    lines.append('# Custom ShockInABox Problem.\n')
    lines.append('# Shock Mach Number: %f\n' % M)
    lines.append('# PreShock Temperature: %e\n' % T1)
    lines.append('# PostShock Temperature: %e\n' % T2)
    lines.append('# PreShock Density: %e\n' % d1)
    lines.append('# \n\n')
    return lines

def write_lines(outfile, lines):
    f = file(outfile,'w')
    f.writelines(lines)
    f.close()

def get_lines(infile):
    f = file(infile,'r')
    lines = f.readlines()
    f.close()
    return lines

def add_lines(infile, outfile, lines):
    orig_lines = get_lines(infile)
    of = file(outfile,'w')
    of.writelines(orig_lines)
    of.writelines(lines)
    of.close()

def TempJump(mach, Gamma=None):
    if Gamma is None:
        Gamma = 5.0/3.0
    M2 = mach*mach;
    TJ = (2.0*Gamma*M2-(Gamma-1.0))*((Gamma-1.0)*M2+2.0)/(M2*(Gamma+1.0)**2);
    return TJ

simtime = 1.0*sec_per_Gyr # 1 Gyr
boxsize = 1.0*cm_per_mpc # 1 Mpc box
pre_shock_den = 1.0e0 # part/cc
postT = 1.0e7 # K
mach = 3.0 
gas_gamma = 5./3.

infile = 'input_shock.enzo'
myname = 'CustomShockBox.enzo'

# No modification is needed below here.

preT = postT/TempJump(mach, gas_gamma)

header = get_header(pre_shock_den, preT, postT, mach)
inlines = get_lines(infile)
customlines = setup_shockbox(boxsize, simtime, pre_shock_den, preT, mach,
                             gamma=gas_gamma)

write_lines(myname, header+inlines+customlines)

