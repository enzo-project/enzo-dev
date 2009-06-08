/***********************************************************************
/
/  GLOBAL DATA DECLARATIONS FOR THE SHOCK POOL PROBLEM
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:
/
/  PURPOSE:
/    This is global data that pertains only to the Shock Pool test problem.
/
************************************************************************/

#ifdef DEFINE_STORAGE
# define SPEXTERN
#else /* DEFINE_STORAGE */
# define SPEXTERN extern
#endif /* DEFINE_STORAGE */

/* Angle of the vector of propogation with respect to the x-axis.  */

SPEXTERN float ShockPoolAngle;

/* Speed of shock itself. */

SPEXTERN float ShockPoolShockSpeed;

/* Density, pressure and velocity of the pre shock state. */

SPEXTERN float ShockPoolDensity;
SPEXTERN float ShockPoolTotalEnergy;
SPEXTERN float ShockPoolVelocity[MAX_DIMENSION];

/* Density, pressure and velocity of the post shock state. */

SPEXTERN float ShockPoolShockDensity;
SPEXTERN float ShockPoolShockTotalEnergy;
SPEXTERN float ShockPoolShockVelocity[MAX_DIMENSION];
 
