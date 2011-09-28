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
# define SEDEXTERN
#else /* DEFINE_STORAGE */
# define SEDEXTERN extern
#endif /* DEFINE_STORAGE */

/* Type of Sedov blast: 
 0) 2D with energy injected into multiple cells
 1) 3D with initial energy distribution from analytical start
 2) 3D with energy injected into multiple cells
 3) 3D with energy injected into a single cell
*/

SEDEXTERN int SedovBlastType;

/* Density and total energy in the ambient medium. */

SEDEXTERN int   SedovBlastFullBox;
SEDEXTERN float SedovBlastDensity;
SEDEXTERN float SedovBlastTotalEnergy;
SEDEXTERN float SedovBlastInputEnergy;

/* Density and total energy in the "inner" state. */

SEDEXTERN float SedovBlastInnerTotalEnergy;

/* Radius in cell lengths (on finest level) that energy is 
   injected into (type 0 and 2) and 
   the total number of cells this equates to. */
 
SEDEXTERN float SedovBlastEnergyZones;
SEDEXTERN int SedovBlastEnergyZonesUsed;
