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
# define SHOCK_POOL_EXTERN
#else /* DEFINE_STORAGE */
# define SHOCK_POOL_EXTERN extern
#endif /* DEFINE_STORAGE */

/* Angle of the vector of propogation with respect to the x-axis.  */

SHOCK_POOL_EXTERN float ShockPoolAngle;

/* Speed of shock itself. */

SHOCK_POOL_EXTERN float ShockPoolShockSpeed;

/* Density, pressure and velocity of the pre shock state. */

SHOCK_POOL_EXTERN float ShockPoolDensity;
SHOCK_POOL_EXTERN float ShockPoolTotalEnergy;
SHOCK_POOL_EXTERN float ShockPoolVelocity[MAX_DIMENSION];

/* Density, pressure and velocity of the post shock state. */

SHOCK_POOL_EXTERN float ShockPoolShockDensity;
SHOCK_POOL_EXTERN float ShockPoolShockTotalEnergy;
SHOCK_POOL_EXTERN float ShockPoolShockVelocity[MAX_DIMENSION];

/* Delay of shock */

SHOCK_POOL_EXTERN float ShockPoolDelay; 
