/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE PROTOSTELLAR CORE COLLAPSE TEST) 
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Alexei Kritsuk, June 2005.
/
/  PURPOSE: Sets the initial core region.
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"



// function prototypes
int GetUnits(float *DensityUnits, float *LengthUnits, 
	     float *TemperatureUnits, float *TimeUnits, 
	     float *VelocityUnits, double *MassUnits, FLOAT Time);



int grid::RadHydroConstTestInitializeGrid(int NumChemicals,
					  float DensityConstant, 
					  float VxConstant, 
					  float VyConstant, 
					  float VzConstant, 
					  float IEConstant, 
					  float EgConstant, 
					  float HMassFrac, 
					  float InitFracHII, 
					  float InitFracHeII, 
					  float InitFracHeIII, 
					  int   local)
{
#ifdef TRANSFER
//   if (debug)
//     fprintf(stdout,"Entering grid::RadHydroConstTestInitializeGrid routine\n");

  // determine whether data should be allocated/initialized or not
  int NewData = TRUE;
  if ((ParallelRootGridIO == TRUE) && (local == 0))
    NewData = FALSE;

  // if grids allocated and already set up (i.e. restart), return
  if ((NumberOfBaryonFields > 5) && (BaryonField[5] != NULL))
    return SUCCESS;


  // create necessary baryon fields
  int RhoNum, TENum, IENum, V0Num, V1Num, V2Num, EgNum, DeNum, 
    HINum, HIINum, HeINum, HeIINum, HeIIINum;
  NumberOfBaryonFields = 0;
  FieldType[RhoNum = NumberOfBaryonFields++] = Density;
  FieldType[TENum = NumberOfBaryonFields++]  = TotalEnergy;
  if (DualEnergyFormalism) 
    FieldType[IENum = NumberOfBaryonFields++]  = InternalEnergy;
  FieldType[V0Num = NumberOfBaryonFields++]    = Velocity1;
  FieldType[V1Num = NumberOfBaryonFields++]    = Velocity2;
  FieldType[V2Num = NumberOfBaryonFields++]    = Velocity3;
  FieldType[EgNum = NumberOfBaryonFields++]    = RadiationFreq0;
  if (NumChemicals > 0) {
    FieldType[DeNum = NumberOfBaryonFields++]    = ElectronDensity;
    FieldType[HINum = NumberOfBaryonFields++]    = HIDensity;
    FieldType[HIINum = NumberOfBaryonFields++]   = HIIDensity;
  }
  if (NumChemicals == 3) {
    FieldType[HeINum = NumberOfBaryonFields++]   = HeIDensity;
    FieldType[HeIINum = NumberOfBaryonFields++]  = HeIIDensity;    
    FieldType[HeIIINum = NumberOfBaryonFields++] = HeIIIDensity;
  }

  // set the subgrid static flag (necessary??)
  SubgridsAreStatic = FALSE;  // no subgrids

  // Return if this doesn't concern us.
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  // Get various units
  double MassUnits=1;
  float DensityUnits=1, LengthUnits=1, TemperatureUnits=1, TimeUnits=1,
    VelocityUnits=1;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
    fprintf(stderr,"Error in GetUnits.\n");
    return FAIL;
  }
  if (debug) {
    fprintf(stdout,"  Internal Unit Conversion Factors:\n");
    fprintf(stdout,"         length = %g\n",LengthUnits);
    fprintf(stdout,"           mass = %lg\n",MassUnits);
    fprintf(stdout,"           time = %g\n",TimeUnits);
  }

  // compute size of fields
  int size = 1;
  for (int dim=0; dim<GridRank; dim++)  size *= GridDimension[dim];
 
  // allocate fields
  if (NewData == TRUE) {
//     printf("\n  P%"ISYM": Allocating %"ISYM" baryon fields of size %"ISYM" (%"ISYM"x%"ISYM"x%"ISYM")\n",
// 	   MyProcessorNumber, NumberOfBaryonFields, size, 
// 	   GridDimension[0], GridDimension[1], GridDimension[2]);

    this->AllocateGrids();
    
    // set fluid density, total energy, [internal energy,] velocities, 
    // radiation energy, electron density, chemical species
    int i;
    float TEConstant = IEConstant + 0.5*(VxConstant*VxConstant 
		       + VyConstant*VyConstant + VzConstant*VzConstant);
    float HIIConstant = InitFracHII*HMassFrac*DensityConstant;
    float HIConstant = HMassFrac*DensityConstant - HIIConstant;
    float HeIIConstant = InitFracHeII*DensityConstant*(1.0-HMassFrac);
    float HeIIIConstant = InitFracHeIII*DensityConstant*(1.0-HMassFrac);
    float HeIConstant = (1.0-HMassFrac)*DensityConstant-HeIIConstant-HeIIIConstant;
    float DeConstant = HIIConstant + 0.25*HeIIConstant + 0.5*HeIIIConstant;
    float eUnits = VelocityUnits*VelocityUnits;
    float EUnits = DensityUnits*eUnits;
    for (i=0; i<size; i++) {
      BaryonField[RhoNum][i] = DensityConstant/DensityUnits;
      BaryonField[TENum][i]  = TEConstant/eUnits;
      BaryonField[V0Num][i]  = VxConstant/VelocityUnits;
      BaryonField[V1Num][i]  = VyConstant/VelocityUnits;
      BaryonField[V2Num][i]  = VzConstant/VelocityUnits;
      BaryonField[EgNum][i]  = EgConstant/EUnits;
      if (NumChemicals > 0) {
	BaryonField[DeNum][i]  = DeConstant/DensityUnits;
	BaryonField[HINum][i]  = HIConstant/DensityUnits;
	BaryonField[HIINum][i] = HIIConstant/DensityUnits;
      }
      if (NumChemicals == 3) {
	BaryonField[HeINum][i]   = HeIConstant/DensityUnits;
	BaryonField[HeIINum][i]  = HeIIConstant/DensityUnits;
	BaryonField[HeIIINum][i] = HeIIIConstant/DensityUnits;
      }
    }
    if (DualEnergyFormalism)
      for (i=0; i<size; i++)
	BaryonField[IENum][i] = IEConstant/eUnits;
    
    if (debug) {
      printf("\n  Initializing constant fields using CGS values:\n");
      printf("        density = %g\n",DensityConstant);
      printf("   total energy = %g\n",TEConstant);
      if (DualEnergyFormalism)
	printf("   internal energy = %g\n",IEConstant);
      printf("     x-velocity = %g\n",VxConstant);
      printf("     y-velocity = %g\n",VyConstant);
      printf("     z-velocity = %g\n",VzConstant);
      printf("      radiation = %g\n",EgConstant);
      if (NumChemicals > 0) {
	printf("      electrons = %g\n",DeConstant);
	printf("            nHI = %g\n",HIConstant);
	printf("           nHII = %g\n",HIIConstant);
      }
      if (NumChemicals == 3) {
	printf("           nHeI = %g\n",HeIConstant);
	printf("          nHeII = %g\n",HeIIConstant);
	printf("         nHeIII = %g\n",HeIIIConstant);
      }

      printf("Corresponding scaled values:\n");
      printf("        density = %g\n",DensityConstant/DensityUnits);
      printf("   total energy = %g\n",TEConstant/eUnits);
      if (DualEnergyFormalism)
	printf("   internal energy = %g\n",IEConstant/eUnits);
      printf("     x-velocity = %g\n",VxConstant/VelocityUnits);
      printf("     y-velocity = %g\n",VyConstant/VelocityUnits);
      printf("     z-velocity = %g\n",VzConstant/VelocityUnits);
      printf("      radiation = %g\n",EgConstant/EUnits);
      if (NumChemicals > 0) {
	printf("      electrons = %g\n",DeConstant/DensityUnits);
	printf("            nHI = %g\n",HIConstant/DensityUnits);
	printf("           nHII = %g\n",HIIConstant/DensityUnits);
      }
      if (NumChemicals == 3) {
	printf("           nHeI = %g\n",HeIConstant/DensityUnits);
	printf("          nHeII = %g\n",HeIIConstant/DensityUnits);
	printf("         nHeIII = %g\n",HeIIIConstant/DensityUnits);
      }
    }

  } // end if NewData == TRUE

  return SUCCESS;

#else

  fprintf(stderr,"Error: TRANSFER must be enabled for this test!\n");
  return FAIL;

#endif

}
