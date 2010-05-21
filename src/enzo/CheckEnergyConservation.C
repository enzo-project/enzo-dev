/***********************************************************************
/
/  CHECK ENERGY CONSERVATION
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
 
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
 
 
int CheckEnergyConservation(HierarchyEntry *Grids[], int grid,
			    int NumberOfGrids, int level, float dt)
{
 
  /* Preserves the energy values:
      0 - baryon kinetic energy
      1 - baryon thermal energy
      2 - particle kinetic energy
      3 - potential energy
      4 - total energy.
      5 - Old dt
      6 - Layzer-Irvine energy term */
 
  static float EnergySum[MAX_DEPTH_OF_HIERARCHY][7];
  static float EnergySumOld[MAX_DEPTH_OF_HIERARCHY][7];
  static int   PrintHeader = TRUE;
 
  int i, l;
 
  /* If this is the first grid, zero the Energy values for this level
     (except for the cumulative LI term). */
 
  if (grid == 0) {
    for (i = 0; i < 6; i++) {
      EnergySumOld[level][i] = EnergySum[level][i];
      EnergySum[level][i] = 0;
    }
    EnergySum[level][5] = dt;
  }
 
  /* Compute energy sum in this grid. */
 
  if (Grids[grid]->GridData->ComputeEnergy(&EnergySum[level][0]) ==
      FAIL) {
    ENZO_FAIL("Error in grid->ComputeEnergy.\n");
  }
 
  /* If this is the last grid, sum it up and output result. */
 
  if (grid == NumberOfGrids-1) {
 
    /* Layzer-Irvine term. */
 
    FLOAT a = 1, ayed, aold, anew;
    FLOAT Time = Grids[grid]->GridData->ReturnTime();
    if (ComovingCoordinates) {
 
      CosmologyComputeExpansionFactor(Time, &a, &ayed);
      CosmologyComputeExpansionFactor(Time-EnergySumOld[level][5],&aold,&ayed);
      CosmologyComputeExpansionFactor(Time+dt, &anew, &ayed);
      EnergySum[level][6] += (anew - a)*(
		    EnergySum[level][0] +
		    EnergySum[level][2] +
		  (3.0*Gamma-4.0)*
		    EnergySum[level][1]
#ifdef UNUSED
		    0.5*(EnergySum[level][0]+EnergySumOld[level][0]) +
		    0.5*(EnergySum[level][2]+EnergySumOld[level][2]) +
		  (3.0*Gamma-4.0)*
		    0.5*(EnergySum[level][1]+EnergySumOld[level][1])
#endif /* UNUSED */
					);
    }
 
    /* Sum individual energy terms. */
 
    for (i = 0; i < 4; i++)
      EnergySum[level][4] += a*EnergySum[level][i];
    EnergySum[level][4] += EnergySum[level][6];
 
    /* Sum over levels. */
 
    float Energy[7];
    for (i = 0; i < 7; i++) {
      Energy[i] = 0;
      for (l = 0; l < MAX_DEPTH_OF_HIERARCHY; l++)
	Energy[i] += EnergySum[l][i];
    }
 
    /* Output results. */
 
    FILE *fptr;
    if ((fptr = fopen("amr_energy.out", "a")) == NULL) {
      ENZO_FAIL("error opening amr_energy.out\n");
    }
 
    if (PrintHeader) {
      fprintf(fptr, "# Time       ");
      if (ComovingCoordinates) fprintf(fptr, "a/ai     ");
      fprintf(fptr, "E_tot    E_kin (gas)  E_gas     E_kin (dm)    E_grav");
      if (ComovingCoordinates) fprintf(fptr,  "   LI term");
      fprintf(fptr, "\n");
      PrintHeader = FALSE;
    }
    fprintf(fptr, "%8.4"FSYM"  ", Time);
    if (ComovingCoordinates) fprintf(fptr, "%8.3"FSYM"  ", a);
    fprintf(fptr, "%8.3e  %8.3e  %8.3e  %8.3e  %8.3e", Energy[4], a*Energy[0],
	    a*Energy[1], a*Energy[2], a*Energy[3]);
    if (ComovingCoordinates) fprintf(fptr, "  %8.3e", Energy[6]);

    fprintf(fptr, "\n");
 
    fclose(fptr);
 
  }
 
  return SUCCESS;
}
