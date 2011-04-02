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
# define SEDOV_BLAST_EXTERN
#else /* DEFINE_STORAGE */
# define SEDOV_BLAST_EXTERN extern
#endif /* DEFINE_STORAGE */

/* Density and total energy in the ambient medium. */

SEDOV_BLAST_EXTERN int   SedovBlastFullBox;
SEDOV_BLAST_EXTERN float SedovBlastDensity;
SEDOV_BLAST_EXTERN float SedovBlastTotalEnergy;

/* Density and total energy in the "inner" state. */

SEDOV_BLAST_EXTERN float SedovBlastInnerTotalEnergy;


 
