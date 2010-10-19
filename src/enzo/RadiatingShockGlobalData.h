/***********************************************************************
/
/  GLOBAL DATA DECLARATIONS FOR THE RADIATING SEDOV BLAST PROBLEM
/
/  written by: Brian O'Shea
/  date:       August 2007
/  modified1:  
/
/  PURPOSE:
/    This is global data that pertains only to the RadiatingShock Blast test problem.
/
************************************************************************/

#ifdef DEFINE_STORAGE
# define RADIATING_SHOCK_EXTERN
#else /* DEFINE_STORAGE */
# define RADIATING_SHOCK_EXTERN extern
#endif /* DEFINE_STORAGE */

/* Density and total energy in the ambient medium. */

RADIATING_SHOCK_EXTERN float RadiatingShockDensity;
RADIATING_SHOCK_EXTERN float RadiatingShockTotalEnergy;

/* Density and total energy in the "inner" state. */

RADIATING_SHOCK_EXTERN float RadiatingShockInnerTotalEnergy;




 
