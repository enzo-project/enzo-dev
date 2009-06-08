/***********************************************************************
/
/  GLOBAL DATA DECLARATIONS FOR THE SEDOV BLAST PROBLEM
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:  Alexei Kritsuk, January 2004
/
/  PURPOSE:
/    This is global data that pertains only to the Sedov Blast test problem.
/
************************************************************************/

#ifdef DEFINE_STORAGE
# define SPEXTERN
#else /* DEFINE_STORAGE */
# define SPEXTERN extern
#endif /* DEFINE_STORAGE */

/* Density and total energy in the ambient medium. */

SPEXTERN int   SedovBlastFullBox;
SPEXTERN float SedovBlastDensity;
SPEXTERN float SedovBlastTotalEnergy;

/* Density and total energy in the "inner" state. */

SPEXTERN float SedovBlastInnerTotalEnergy;


 
