/***********************************************************************
/
/  EXTERNAL BOUNDARY CLASS (IDENTIFY COMMONLY USED VARIABLES FROM THE LIST)
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:
/
/  PURPOSE:
/
/  NOTE:
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
#include "MHD2DTestGlobalData.h"
 
/* function prototypes */
 
int FindField(int f, int farray[], int n);
 
 
 
int ExternalBoundary::IdentifyPhysicalQuantities(int &DensNum, int &GENum,
						 int &Vel1Num, int &Vel2Num,
						 int &Vel3Num, int &TENum)
{
 
  DensNum = GENum = Vel1Num = Vel2Num = Vel3Num = TENum = 0;
 
  /* Find Density, if possible. */
 
  if ((DensNum = FindField(Density, BoundaryFieldType, NumberOfBaryonFields))
      < 0) {
    ENZO_FAIL("EBIPQ: Cannot find density.\n");
  }
 
  /* Find Total energy, if possible. */
 
  if ((TENum = FindField(TotalEnergy, BoundaryFieldType, NumberOfBaryonFields))
      < 0) {
    ENZO_FAIL("Cannot find total energy.\n");
  }
 
  /* Find gas energy, if possible. */
 
  if (DualEnergyFormalism == TRUE)
    if ((GENum = FindField(InternalEnergy, BoundaryFieldType,
			   NumberOfBaryonFields)) < 0) {
      ENZO_FAIL("Cannot find gas energy.\n");
    }
 
  /* Find Velocity1, if possible. */
 
  if ((Vel1Num = FindField(Velocity1, BoundaryFieldType, NumberOfBaryonFields))
      < 0) {
    ENZO_FAIL("Cannot find Velocity1.\n");
  }
 
  /* Find Velocity2, if possible. */
 
  if (BoundaryRank > 1)
    if ((Vel2Num = FindField(Velocity2, BoundaryFieldType,
			     NumberOfBaryonFields)) < 0) {
      ENZO_FAIL("Cannot find Velocity2.\n");
    }
 
  /* Find Velocity3, if possible. */
 
  if (BoundaryRank > 2)
    if ((Vel3Num = FindField(Velocity3, BoundaryFieldType,
			     NumberOfBaryonFields)) == 0) {
      ENZO_FAIL("Cannot find Velocity3.\n");

    }
 
  return SUCCESS;
}


int ExternalBoundary::IdentifyPhysicalQuantities(int &DensNum, int &GENum, int &Vel1Num, 
						 int &Vel2Num, int &Vel3Num, int &TENum,
						 int &B1Num, int &B2Num, int &B3Num, 
						 int &PhiNum)
{

  DensNum = GENum = Vel1Num = Vel2Num = Vel3Num = TENum = B1Num = B2Num = B3Num = PhiNum = 0;
    
  /* Find Density, if possible. */

  if ((DensNum = FindField(Density, BoundaryFieldType, NumberOfBaryonFields)) < 0) {
        ENZO_FAIL("Cannot find density.");
  }

  /* Find Total energy, if possible. */

  if ((TENum = FindField(TotalEnergy, BoundaryFieldType, NumberOfBaryonFields)) < 0) {
        ENZO_FAIL("Cannot find total energy.");
  }

  /* Find gas energy, if possible. */

  if (DualEnergyFormalism == TRUE) {
    if ((GENum = FindField(InternalEnergy, BoundaryFieldType,
			   NumberOfBaryonFields)) < 0) {
            ENZO_FAIL("Cannot find gas energy.");
    }
  }

  /* Find Velocity1, if possible. */
   
  if ((Vel1Num = FindField(Velocity1, BoundaryFieldType, NumberOfBaryonFields)) < 0) {
        ENZO_FAIL("Cannot find Velocity1.");
  }

  /* Find Velocity2, if possible. */

  if (BoundaryRank > 1 || HydroMethod == MHD_RK)
      if ((Vel2Num = FindField(Velocity2, BoundaryFieldType, 
			   NumberOfBaryonFields)) < 0) {
        ENZO_FAIL("Cannot find Velocity2.");
  }

  /* Find Velocity3, if possible. */

  if (BoundaryRank > 2 || HydroMethod == MHD_RK)
      if ((Vel3Num = FindField(Velocity3, BoundaryFieldType, 
                      NumberOfBaryonFields)) == 0) {
          ENZO_FAIL("Cannot find Velocity3.");
      }

  if (!UseMHD) {
    return SUCCESS;
  }

  /* Find magnetic field components */

  if ((B1Num = FindField(Bfield1, BoundaryFieldType, NumberOfBaryonFields)) < 0) {
        ENZO_FAIL("Cannot find Bfield1.");
  }

  if ((B2Num = FindField(Bfield2, BoundaryFieldType, NumberOfBaryonFields)) < 0) {
        ENZO_FAIL("Cannot find Bfield2.");
  }
   
  if ((B3Num = FindField(Bfield3, BoundaryFieldType, NumberOfBaryonFields)) < 0) {
        ENZO_FAIL("Cannot find Bfield3.");
  }

  if( HydroMethod == MHD_RK ){
      if ((PhiNum = FindField(PhiField, BoundaryFieldType, NumberOfBaryonFields)) < 0) {
            ENZO_FAIL("Cannot find Phi field.");
      }
  }

  return SUCCESS;
}
