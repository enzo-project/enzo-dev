/***********************************************************************
/
/  RECALIBRATE ACCRETING MASS WHEN REQUESTED
/
/  written by: Ji-hoon Kim
/  date:       June, 2010
/  modified1: 
/
/  PURPOSE: This routine is implemented to recalibrate (correct) the 
/           accreting rate calculated in Star_CalculateMassAccretion,
/           if requested by MBHAccretingMassRatio = BONDI_ACCRETION_
/           CORRECT_NUMERICAL. This process is essential if the Bondi 
/           radius is bigger than the cell size where MBH is located.
/           
/  CAUTION: This is currently only for MBH Accretion.  
/         
/  INPUT:   BondiRadius and density from Star_CalculateMassAccretion.
/           Both are in their code units.  density is not bounded by
/           the Eddington limit. (yet!)
/
************************************************************************/
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
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
#include "CommunicationUtilities.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);

int RecalibrateAccretingMass(FLOAT star_pos[], LevelHierarchyEntry *LevelArray[], 
			     int level, float &BondiRadius, float &density,
			     float &RecalibrateAccretingMassRatio)
{

  const double pc = 3.086e18, Msun = 1.989e33, pMass = 1.673e-24, 
    gravConst = 6.673e-8, yr = 3.1557e7, Myr = 3.1557e13;

  float AvgVelocity[MAX_DIMENSION], MassEnclosed[2], Metallicity = 0, ColdGasMass = 0;
  float OneOverRSquaredSum, average_density_at_Bondi_radius = 0.0;
  int i, l, dim, count, FirstLoop = TRUE;
  LevelHierarchyEntry *Temp, *Temp2;
  FLOAT Time;

  /* This routine is invoked only when BONDI_ACCRETION_CORRECT_NUMERICAL */

  if (MBHAccretion <= 0 || 
      BondiRadius <= 0.0 || density <= 0.0 ||
      MBHAccretingMassRatio != BONDI_ACCRETION_CORRECT_NUMERICAL) 
    ENZO_FAIL("RecalibrateAM: Something is wrong!");
//    return SUCCESS;

//    printf("RecalibrateAM: Proc = %d, star_pos[] = %g, %g, %g\n", 
// 	  MyProcessorNumber, star_pos[0], star_pos[1], star_pos[2]);  

  /* Set the units. */

  Temp = LevelArray[level];
  Time = Temp->GridData->ReturnTime();
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits;
  double MassUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, &MassUnits, Time);

  for (dim = 0; dim < MAX_DIMENSION; dim++)
    AvgVelocity[dim] = 0.0;

  for (count = 0; count < 2; count++)
    MassEnclosed[count] = 0.0;
				   
  /* Get the cell width */

  int Rank, Dims[MAX_DIMENSION];
  float CellWidth;
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
  LevelArray[level]->GridData->ReturnGridInfo(&Rank, Dims, LeftEdge, RightEdge);
  CellWidth = (RightEdge[0] - LeftEdge[0]) / (Dims[0] - 2*NumberOfGhostZones);

  RecalibrateAccretingMassRatio = 1.0;

//   printf("RecalibrateAM1: dens = %g, ave_dens = %g, ME[0] = %g, ME[1] = %g, BondiR = %g, CellW = %g\n",
// 	 density, average_density_at_Bondi_radius, MassEnclosed[0], MassEnclosed[1], BondiRadius, CellWidth);  

  /* No need to proceed if the Bondi radius is smaller than the cell width */
  
  if (CellWidth >= BondiRadius)
    return SUCCESS;


  /***********************************************************************
                  FIND THE AVERAGE DENSITY AT BONDI RADIUS
  ***********************************************************************/

  /* Find the enclosed mass in the sphere of R_B, 
     and in the sphere of R_B + dx */

  for (count = 0; count < 2; count++) { 

    MassEnclosed[count] = 0;
    Metallicity = 0;
    ColdGasMass = 0;
    OneOverRSquaredSum = 0;
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      AvgVelocity[dim] = 0.0;
    for (l = 0; l < MAX_DEPTH_OF_HIERARCHY; l++) {
      Temp = LevelArray[l];
      while (Temp != NULL) {

	/* Zero under subgrid field */

	if (FirstLoop == 1) {
	  Temp->GridData->
	    ZeroSolutionUnderSubgrid(NULL, ZERO_UNDER_SUBGRID_FIELD);
	  Temp2 = LevelArray[l+1];
	  while (Temp2 != NULL) {
	    Temp->GridData->ZeroSolutionUnderSubgrid(Temp2->GridData, 
						     ZERO_UNDER_SUBGRID_FIELD);
	    Temp2 = Temp2->NextGridThisLevel;
	  }
	}

	/* Sum enclosed mass in this grid */

	if (Temp->GridData->GetEnclosedMass(star_pos, BondiRadius, MassEnclosed[count], 
					    Metallicity, ColdGasMass, 
					    AvgVelocity, OneOverRSquaredSum) == FAIL) {
	  ENZO_FAIL("Error in GetEnclosedMass.");
	}

	Temp = Temp->NextGridThisLevel;

      } // END: Grids

      if (l == MAX_DEPTH_OF_HIERARCHY-1)
	FirstLoop = 0;

    } // END: level

    CommunicationAllSumValues(&Metallicity, 1);
    CommunicationAllSumValues(&MassEnclosed[count], 1);
    CommunicationAllSumValues(&ColdGasMass, 1);
    CommunicationAllSumValues(&OneOverRSquaredSum, 1);
    CommunicationAllSumValues(AvgVelocity, 3);

    if (MassEnclosed[count] == 0) {
      ENZO_FAIL("RecalibrateAM: MassEnclosed = 0; something is wrong!\n");
    }

    Metallicity /= MassEnclosed[count];
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      AvgVelocity[dim] /= MassEnclosed[count];

    /* count=0 is at BondiRadius, count=1 is at BondiRadius+CellWidth */

    BondiRadius += CellWidth;

  }  // ENDFOR count

  /* Find the average density at BondiRadius, using 
     rho = dM/dV = { M(R_B+dx) - M(R_B) } / { 4pi(R_B^2)*dx } */

  average_density_at_Bondi_radius = (MassEnclosed[1] - MassEnclosed[0]) * Msun / MassUnits /
    (4 * PI * BondiRadius * BondiRadius * CellWidth);

  /* Find the correction ratio using "average_density_at_Bondi_radius" 
     and "density" (peak density at the site of MBH) */

  RecalibrateAccretingMassRatio = average_density_at_Bondi_radius / density;

//   printf("RecalibrateAM2: dens = %g, ave_dens = %g, ME[0] = %g, ME[1] = %g, BondiR = %g, CellW = %g, Recab = %g\n",
// 	 density, average_density_at_Bondi_radius, MassEnclosed[0], MassEnclosed[1], 
// 	 BondiRadius, CellWidth, RecalibrateAccretingMassRatio);  
	 
  return SUCCESS;

}
