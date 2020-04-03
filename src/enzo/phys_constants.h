#ifndef __phys_constants_h_
#define __phys_constants_h_
/***********************************************************************
/  
/ DEFINITION OF PHYSICAL CONSTANTS
/  
/ written by: Elizabeth Tasker (renamed by Daniel Reynolds)
/ date:       May, 2005 
/
/ Note: CGS units
/
*********************************************************************/

/* Physics constants */

/************************************************/

/* Planck's constant  [erg s] */

#define h_planck                        6.6260693E-27

/* Boltzmann's constant [cm2gs-2K-1] or [ergK-1] */

#define kboltz                          1.3806504e-16

/* Mass of hydrogen [g] */

#define mh                              1.67262171e-24   

/* Mass of an electron [g] */

#define me                              9.10938215e-28

/* Pi */

#define pi                              3.14159265358979323846

/* ergs per eV */

#define erg_eV                          1.60217653e-12

/* eV per erg */

#define eV_erg                          6.24150948e11

/* HI ionization cross section [cm^-2]*/

#define sigmaHI                         6.345e-18

/* Thompson cross section [cm^-2] */

#define sigma_thompson                  6.65E-25


/************************************************/

/* Astronomical constant */

/************************************************/

/* Speed of light [cms-1] */ 

#define clight                          2.99792458e10

/* Gravitational constant [cm3g-1s-2]*/

#define GravConst                       6.67428e-8

/* Solar mass [g] */

#define SolarMass                       1.9891e33

/* Solar Radius [cm] */

#define SolarRadius                     6.96e10

/* Solar Luminosity [erg/s] */

#define SolarLuminosity                 3.9e33



/* Megaparsec [cm] */

#define Mpc_cm                             3.0857e24

/* kiloparsec [cm] */

#define kpc_cm                             3.0857e21

/* parsec [cm] */

#define pc_cm                              3.0857e18

/* kilometer [cm] */

#define km_cm                              1.0E5

/* Megayear [s] */

#define Myr_s                              3.1556952E13

/* year [s] */

#define yr_s                               3.1556952E7


/************************************************/

/* Other constants - physical properties for convenience */

/************************************************/

/* Average energy [eV] of LW-band photon (11.2-13.6 eV) */

#define LW_photon_energy                    12.8

/* LW-band photon energy threshold [eV] */

#define LW_threshold_energy                 11.2

/* Average energy [eV] of FUV-band photon (5.6 - 11.2 eV) */

#define FUV_photon_energy                    8.4

/* FUV-band photon energy threshold [eV] */

#define FUV_threshold_energy                 5.6

/* Average energy [eV] of IR-band photon (0.75 - 5.6 eV) */

#define IR_photon_energy                     2.0

/* IR-band photon energy threshold [eV] */

#define IR_threshold_energy                 0.76

/* HI ionization threshold energy [eV] */

#define HI_ionizing_energy                13.5984

/* HeI ionization threshold energy [eV] */

#define HeI_ionizing_energy               24.5874

/* HeII ionization threshold energy [eV] */

#define HeII_ionizing_energy              54.42

#endif
