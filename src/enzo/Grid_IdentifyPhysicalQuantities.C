/***********************************************************************
/
/  GRID CLASS (IDENTIFY CERTAIN COMMONLY USED VARIABLES FROM THE LIST)
/
/  written by: Greg Bryan
/  date:       April, 1995
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
 
/* function prototypes */
 
int FindField(int f, int farray[], int n);
  
int grid::IdentifyPhysicalQuantities(int &DensNum, int &GENum, int &Vel1Num,
				     int &Vel2Num, int &Vel3Num, int &TENum)
{
 
  DensNum = GENum = Vel1Num = Vel2Num = Vel3Num = TENum = 0;
 
  /* Find Density, if possible. */
 
  if ((DensNum = FindField(Density, FieldType, NumberOfBaryonFields)) < 0) {
    ENZO_FAIL("GIPQ: Cannot find density.");
  }
 
  /* Find Total energy, if possible. */
 
  if ((TENum = FindField(TotalEnergy, FieldType, NumberOfBaryonFields)) < 0) {
    if( EquationOfState != 1 )
      ENZO_FAIL("Cannot find total energy.");
  }
 
  /* Find gas energy, if possible. */
 
  if (DualEnergyFormalism == TRUE)
    if ((GENum = FindField(InternalEnergy, FieldType,
			   NumberOfBaryonFields)) < 0) {
      ENZO_FAIL("Cannot find gas energy.");
    }
 
  /* Find Velocity1, if possible. */
 
  if ((Vel1Num = FindField(Velocity1, FieldType, NumberOfBaryonFields)) < 0) {
    ENZO_FAIL("Cannot find Velocity1.");
  }
 
  /* Find Velocity2, if possible. */
 
  if (MaxVelocityIndex>1)
    if ((Vel2Num = FindField(Velocity2, FieldType,
			     NumberOfBaryonFields)) < 0) {
      ENZO_FAIL("Cannot find Velocity2.");
    }
 
  /* Find Velocity3, if possible. */
 
  if (MaxVelocityIndex>2)
    if ((Vel3Num = FindField(Velocity3, FieldType,
			     NumberOfBaryonFields)) == 0) {
      ENZO_FAIL("Cannot find Velocity3.");
    }
 
  return SUCCESS;
}

int grid::IdentifyPhysicalQuantities(int &DensNum, int &GENum, int &Vel1Num, 
				     int &Vel2Num, int &Vel3Num, int &TENum,
				     int &B1Num, int &B2Num, int &B3Num)
{

  DensNum = GENum = Vel1Num = Vel2Num = Vel3Num = TENum = 0;
    
  /* Find Density, if possible. */

  if ((DensNum = FindField(Density, FieldType, NumberOfBaryonFields)) < 0) {
        ENZO_FAIL("Cannot find density.");
  }

  /* Find Total energy, if possible. */

  if ((TENum = FindField(TotalEnergy, FieldType, NumberOfBaryonFields)) < 0) {
    if(EquationOfState != 1)
      ENZO_FAIL("Cannot find total energy.");
  }

  /* Find gas energy, if possible. */

  if (DualEnergyFormalism == TRUE) {
    if ((GENum = FindField(InternalEnergy, FieldType,
			   NumberOfBaryonFields)) < 0) {
      ENZO_FAIL("Cannot find gas energy.");
    }
  }

  /* Find Velocity1, if possible. */
   
  if ((Vel1Num = FindField(Velocity1, FieldType, NumberOfBaryonFields)) < 0) {
        ENZO_FAIL("Cannot find Velocity1.");
  }

  /* Find Velocity2, if possible. */

  if (MaxVelocityIndex>1)
    if ((Vel2Num = FindField(Velocity2, FieldType, 
			     NumberOfBaryonFields)) < 0) {
      ENZO_FAIL("Cannot find Velocity2.");
    }
  

  /* Find Velocity3, if possible. */
  if (MaxVelocityIndex>2)
    if ((Vel3Num = FindField(Velocity3, FieldType, 
			     NumberOfBaryonFields)) == 0) {
      ENZO_FAIL("Cannot find Velocity3.");
    }
  
  if (!UseMHD) {
    return SUCCESS;
  }

  /* Find magnetic field components */

  if ((B1Num = FindField(Bfield1, FieldType, NumberOfBaryonFields)) < 0) {
        ENZO_FAIL("Cannot find Bfield1.");
  }

  if ((B2Num = FindField(Bfield2, FieldType, NumberOfBaryonFields)) < 0) {
        ENZO_FAIL("Cannot find Bfield2.");
  }
  
  if ((B3Num = FindField(Bfield3, FieldType, NumberOfBaryonFields)) < 0) {
        ENZO_FAIL("Cannot find Bfield3.");
  }

  return SUCCESS;
}

int grid::IdentifyPhysicalQuantities(int &DensNum, int &GENum, int &Vel1Num, 
				     int &Vel2Num, int &Vel3Num, int &TENum,
				     int &B1Num, int &B2Num, int &B3Num, int &PhiNum)
{

  DensNum = GENum = Vel1Num = Vel2Num = Vel3Num = TENum = 0;
    
  /* Find Density, if possible. */

  if ((DensNum = FindField(Density, FieldType, NumberOfBaryonFields)) < 0) {
        ENZO_FAIL("Cannot find density.");
  }

  /* Find Total energy, if possible. */

  if ((TENum = FindField(TotalEnergy, FieldType, NumberOfBaryonFields)) < 0) {
    if( EquationOfState != 1)
      ENZO_FAIL("Cannot find total energy.");
  }

  /* Find gas energy, if possible. */

  if (DualEnergyFormalism == TRUE) {
    if ((GENum = FindField(InternalEnergy, FieldType,
			   NumberOfBaryonFields)) < 0) {
            ENZO_FAIL("Cannot find gas energy.");
    }
  }

  /* Find Velocity1, if possible. */
   
  if ((Vel1Num = FindField(Velocity1, FieldType, NumberOfBaryonFields)) < 0) {
        ENZO_FAIL("Cannot find Velocity1.");
  }

  /* Find Velocity2, if possible. */

  if (MaxVelocityIndex>1)
      if ((Vel2Num = FindField(Velocity2, FieldType, 
			   NumberOfBaryonFields)) < 0) {
        ENZO_FAIL("Cannot find Velocity2.");
  }

  /* Find Velocity3, if possible. */

  if (MaxVelocityIndex>2)
      if ((Vel3Num = FindField(Velocity3, FieldType, 
                      NumberOfBaryonFields)) == 0) {
          ENZO_FAIL("Cannot find Velocity3.");
      }

  if ( ! UseMHD ) {
    return SUCCESS;
  }

  /* Find magnetic field components */

  if ((B1Num = FindField(Bfield1, FieldType, NumberOfBaryonFields)) < 0) {
        ENZO_FAIL("Cannot find Bfield1.");
  }

  if ((B2Num = FindField(Bfield2, FieldType, NumberOfBaryonFields)) < 0) {
        ENZO_FAIL("Cannot find Bfield2.");
  }
   
  if ((B3Num = FindField(Bfield3, FieldType, NumberOfBaryonFields)) < 0) {
        ENZO_FAIL("Cannot find Bfield3.");
  }

  if ( HydroMethod == MHD_RK ){
      if ((PhiNum = FindField(PhiField, FieldType, NumberOfBaryonFields)) < 0) {
            ENZO_FAIL("Cannot find Phi field.");
      }
  }

  return SUCCESS;
}


int grid::IdentifyDrivingFields(int &Drive1Num, int &Drive2Num, int &Drive3Num)
{

  Drive1Num = Drive2Num = Drive3Num = 0;

  if ((Drive1Num = FindField(DrivingField1, FieldType, NumberOfBaryonFields)) < 0) {
        ENZO_FAIL("Cannot find DriveField1.");
  }

  if ((Drive2Num = FindField(DrivingField2, FieldType, NumberOfBaryonFields)) < 0) {
        ENZO_FAIL("Cannot find DriveField2.");
  }

  if ((Drive3Num = FindField(DrivingField3, FieldType, NumberOfBaryonFields)) < 0) {
        ENZO_FAIL("Cannot find DriveField3.");
  }

  return SUCCESS;
}

int grid::IdentifyPotentialField(int &PotenNum, int &Acce1Num, int &Acce2Num, int &Acce3Num)
{

  PotenNum = 0;

  if ((PotenNum = FindField(GravPotential, FieldType, NumberOfBaryonFields)) < 0) {
        ENZO_FAIL("Cannot find PotentialField.");
  }

  if ((Acce1Num = FindField(AccelerationField1, FieldType, NumberOfBaryonFields)) < 0) {
        ENZO_FAIL("Cannot find AccelerationField.");
  }

  if ((Acce2Num = FindField(AccelerationField2, FieldType, NumberOfBaryonFields)) < 0) {
        ENZO_FAIL("Cannot find AccelerationField.");
  }

  if ((Acce3Num = FindField(AccelerationField3, FieldType, NumberOfBaryonFields)) < 0) {
        ENZO_FAIL("Cannot find AccelerationField.");
  }

  return SUCCESS;
}
