/***********************************************************************
/
/  GLOBAL DATA DECLARATIONS FOR THE IMPLOSION PROBLEM
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:  Alexei Kritsuk, December 2004
/
/  PURPOSE:
/    This is global data that pertains only to the Implosion test problem.
/
************************************************************************/

#ifdef DEFINE_STORAGE
# define SPEXTERN
#else /* DEFINE_STORAGE */
# define SPEXTERN extern
#endif /* DEFINE_STORAGE */

/* Density, pressure and velocity of the surrounding state. */

SPEXTERN float ImplosionDensity;
SPEXTERN float ImplosionTotalEnergy;

/* Density and pressure of the "diamond" state. */

SPEXTERN float ImplosionDiamondDensity;
SPEXTERN float ImplosionDiamondTotalEnergy;

 
