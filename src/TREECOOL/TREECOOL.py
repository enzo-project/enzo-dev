"""
Stephen Skory
sskory@physics.ucsd.edu
Jan 2010

Output the TREECOOL file for Gadget equilibrium cooling. This is a pythonification
of an IDL script provided by Pascal Paschos.
"""

import math

###############################
###### Start Editing ##########
###############################

settings = {
# 1 - Haardt and Madau (1996) quasar spectrum (alpha_q = 1.5)
# 2 - Haardt and Madau (1996) quasar spectrum (alpha_q = 1.8)
# 12 - Haardt and Madau (2001) quasar spectra (alpha_q = 1.57)
"RadiationFieldType"  : 12,
# This applies to all species (HI, HeI, HeII) uniformly in time.
# Values above 1.0 raise the photoionization rates.
"SetUVBAmplitude"     : 1.0,
# This is similar to the above factor, but it works only on HeII.
# It works on top of the above factor, so making both factors
# greater than 1 increases the HeII heating higher than the other two.
"SetHeIIHeatingScale" : 1.4,

##############
# HI and HeI #
##############
"Ramp"                : 1.0,
# Input the redshift at which to start up the background radiation.
# This is where the energy starts being non-zero, but is not at full strength.
"RadRedShiftOn"       : 7.0,
"RadRedShiftOff"      : 0.0,
# Input the redshift at which the background radiation is at full strength.
"RadRedShiftFullOn"   : 6.0,
"RadRedShiftDropOff"  : 0.0,
"RadRedShiftf0to3"    : 0.1,
########
# HeII #
########
"RampX"               : 1.0,
# Input the redshift at which to start up the background radiation.
# This is where the energy starts being non-zero, but is not at full strength.
# This covers HeII.
"XRadRedShiftOn"      : 7.0,
"XRadRedShiftOff"     : 0.0,
# Input the redshift at which the background radiation is at full strength.
"XRadRedShiftFullOn"  : 6.0,
"XRadRedShiftDropOff" : 0.0,
"XRadRedShiftf0to3"   : 0.1,

# Nearly always 0 unless you know what you're doing!
"zc_init"             : 0.0,
# Typically make this somewhat larger than your true
# maximum desired redshift. This will add a padding of zeros to the table
# for interpolation. This should generally be larger than (X)RadRedShiftOn.
"zc_final"            : 7.5,

# The number of lines written out in TREECOOL
"bins"                : 200
}

####################################
########### Stop Editing ###########
####################################

# Print current settings.
print "Current settings:\n"
for param in settings:
    print str(param) + " = " + str(settings[param])

ainvinit    = math.log10(1+settings["zc_init"])
ainvfinal   = math.log10(1+settings["zc_final"])
ainvarr     = []
zcarr       = []
for i in range(settings["bins"] + 1):
    ainvarr.append(float(i) / settings["bins"] * (ainvfinal-ainvinit))
    zcarr.append(10**ainvarr[i]-1.0)

fp = open("TREECOOL.mod", "w")

for ii,zc in enumerate(zcarr):
    if zc > settings["XRadRedShiftOn"] or zc < settings["XRadRedShiftOff"]:
        settings["RampX"] = 0.0
    else:
        if zc > settings["XRadRedShiftFullOn"]:
            settings["RampX"] = 0.5 - 0.5*math.tanh(15.0*(zc-0.5*\
                (settings["XRadRedShiftOn"] \
                + settings["XRadRedShiftFullOn"])))
        
        if zc < settings["XRadRedShiftDropOff"]:
            settings["RampX"] = ( zc - settings["XRadRedShiftDropOff"] + \
                settings["XRadRedShiftf0to3"] * \
                (settings["XRadRedShiftDropOff"]-zc)) / \
                (settings["XRadRedShiftDropOff"] - \
                settings["XRadRedShiftOff"])
    
    if zc > settings["RadRedShiftOn"] or zc < settings["RadRedShiftOff"]:
        settings["Ramp"] = 0.0
    else:
    
        if zc > settings["RadRedShiftFullOn"]:
            settings["Ramp"] = 0.5 - 0.5*math.tanh(15.0*(zc-0.5*\
                (settings["RadRedShiftOn"]+settings["RadRedShiftFullOn"])))
        if zc < settings["RadRedShiftDropOff"]:
            settings["Ramp"] = (zc - settings["RadRedShiftDropOff"]+\
            settings["RadRedShiftf0to3"] * \
            (settings["RadRedShiftDropOff"] - zc)) / \
            (settings["RadRedShiftDropOff"] - \
            settings["RadRedShiftOff"])
    
    exp_arg = -1.0 * (zc-2.3)**2
    
    #--------------------------------------------------------------------------
    #1) For the Haardt and Madau (1996) quasar spectrum (alpha_q = 1.5) 
    
    if (settings["RadiationFieldType"] == 1):
        k24 = 6.7e-13 * (1.0+zc)**(0.43) * math.exp(exp_arg/1.95)
        k25 = 6.3e-15 * (1.0+zc)**(0.51) * math.exp(exp_arg/2.35)
        k26 = 3.2e-13 * (1.0+zc)**(0.50) * math.exp(exp_arg/2.00)
        piHI   = 4.7e-24 * (1.0+zc)**(0.43) * math.exp(exp_arg/1.95)
        piHeI  = 8.2e-24 * (1.0+zc)**(0.50) * math.exp(exp_arg/2.00)
        piHeII = 1.6e-25 * (1.0+zc)**(0.51) * math.exp(exp_arg/2.35)
    
    #------------------------------------------------------------------
    #2) For the Haardt and Madau (1996) quasar spectrum (alpha_q = 1.8)
    
    if (settings["RadiationFieldType"] == 2):
        k24 = 5.6e-13 * (1.0+zc)**(0.43) * math.exp(exp_arg/1.95)
        k25 = 3.2e-15 * (1.0+zc)**(0.30) * math.exp(exp_arg/2.60)
        k26 = 4.8e-13 * (1.0+zc)**(0.43) * math.exp(exp_arg/1.95)
        piHI   = 3.9e-24 * (1.0+zc)**(0.43) * math.exp(exp_arg/1.95)
        piHeI  = 6.4e-24 * (1.0+zc)**(0.43) * math.exp(exp_arg/2.10)
        piHeII = 8.7e-26 * (1.0+zc)**(0.30) * math.exp(exp_arg/2.70)
    
    #-------------------------------------------------------------------------
    #12) For Haardt and Madau (2001) quasar spectra (alpha_q = 1.57)
    
    
    if (settings["RadiationFieldType"] == 12):
        k24 = 1.04e-12 * (1.0+zc)**(0.231) \
            * math.exp( -0.6818 * (zc-1.855)**(2.0) / \
            (1.0+0.1646 * (zc+0.3097)**(2.0)) )
        k25 = 1.84e-14 * (1.0+zc)**(-1.038 ) \
            * math.exp( -1.1640 * (zc-1.973)**(2.0) / \
            (1.0+0.1940 * (zc-0.6561)**(2.0)) )
        k26 = 5.79e-13 * (1.0+zc)**(0.278) \
            * math.exp( -0.8260 * (zc-1.973)**(2.0) / \
            (1.0+0.1730 * (zc+0.2880)**(2.0)) )
        piHI = 8.86e-24 * (1.0+zc)**(-0.0290) \
            * math.exp( -0.7055 * (zc-2.003)**(2.0) / \
            (1.0+0.1884 * (zc+0.2888)**(2.0)) )    
        piHeI = 5.86e-24 * (1.0+zc)**(0.1764) \
            * math.exp( -0.8029 * (zc-2.088)**(2.0) / \
            (1.0+0.1732 * (zc+0.1362)**(2.0)) )
        piHeII = 2.17e-25 * (1.0+zc)**(-0.2196) \
            * math.exp( -1.070 *(zc-1.782)**(2.0) / \
            (1.0+0.2124 * (zc-0.9213)**(2.0)) )
    
    k24    = k24 * settings["SetUVBAmplitude"] * settings["Ramp"]
    k25    = k25 * settings["SetUVBAmplitude"] * settings["RampX"]
    k26    = k26 * settings["SetUVBAmplitude"] * settings["Ramp"]
    piHI   = piHI * settings["SetUVBAmplitude"] * settings["Ramp"]
    piHeI  = piHeI * settings["SetUVBAmplitude"] * settings["Ramp"]
    piHeII = piHeII * (settings["SetUVBAmplitude"] * \
        settings["SetHeIIHeatingScale"]) * settings["RampX"]
    
    line = "%1.4f %1.7e %1.7e %1.7e %1.7e %1.7e %1.7e\n" % \
        (ainvarr[ii],k24,k26,k25,piHI,piHeI,piHeII)
    fp.write(line)

fp.close()
print '\nAll done. Saved to "TREECOOL.mod"'