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



#endif
