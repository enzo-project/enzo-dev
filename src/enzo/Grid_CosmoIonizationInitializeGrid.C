/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A COLLAPSE TEST)
/
/  written by: Daniel Reynolds
/  date:       May 2008
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/
 
#include <stdio.h>
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
#include "CosmologyParameters.h"
 
// function prototypes
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
 

int grid::CosmoIonizationInitializeGrid(int NumChemicals,
					float VxConstant, 
					float VyConstant, 
					float VzConstant, 
					float IEConstant, 
					float EgConstant, 
					float InitialFractionHII, 
					float OmegaBaryonNow,
					int local)
{
#ifdef TRANSFER
//   if (MyProcessorNumber == ROOT_PROCESSOR)
//     printf("Entering grid::CosmoIonizationInitializeGrid routine\n");

  // determine whether data should be allocated/initialized or not
  int NewData = TRUE;
  if ((ParallelRootGridIO == TRUE) && (local == 0))
    NewData = FALSE;

  // if grids allocated and already set up (i.e. restart), return
  if ((NumberOfBaryonFields > 5) && (BaryonField[5] != NULL))
    return SUCCESS;

  // create necessary baryon fields
  int RhoNum, TENum, GENum, V0Num, V1Num, V2Num, EgNum, DeNum, 
    HINum, HIINum;
  NumberOfBaryonFields = 0;
  FieldType[RhoNum = NumberOfBaryonFields++]   = Density;
  FieldType[TENum = NumberOfBaryonFields++]    = TotalEnergy;
  if (DualEnergyFormalism) 
    FieldType[GENum = NumberOfBaryonFields++]  = InternalEnergy;
  FieldType[V0Num = NumberOfBaryonFields++]    = Velocity1;
  FieldType[V1Num = NumberOfBaryonFields++]    = Velocity2;
  FieldType[V2Num = NumberOfBaryonFields++]    = Velocity3;
  FieldType[EgNum = NumberOfBaryonFields++]    = RadiationFreq0;
  FieldType[DeNum = NumberOfBaryonFields++]    = ElectronDensity;
  FieldType[HINum = NumberOfBaryonFields++]    = HIDensity;
  FieldType[HIINum = NumberOfBaryonFields++]   = HIIDensity;

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
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    printf("  Internal Unit Conversion Factors:\n");
    printf("         length = %g\n",LengthUnits);
    printf("           mass = %lg\n",MassUnits);
    printf("           time = %g\n",TimeUnits);
  }

  // compute size of fields
  int size = 1;
  for (int dim=0; dim<GridRank; dim++)  size *= GridDimension[dim];
 
  // allocate fields
  if (NewData == TRUE) {
//     printf("\n  P%"ISYM": Allocating %"ISYM" baryon fields of size %"ISYM" (%"ISYM"x%"ISYM"x%"ISYM")\n",
// 	   MyProcessorNumber, NumberOfBaryonFields, size, 
// 	   GridDimension[0], GridDimension[1], GridDimension[2]);

    for (int field=0; field<NumberOfBaryonFields; field++)
      if (BaryonField[field] == NULL)
	BaryonField[field] = new float[size];
 
    // set field data
    int i;
    float TEConstant = (IEConstant + 0.5*(VxConstant*VxConstant + 
					  VyConstant*VyConstant + 
					  VzConstant*VzConstant));
    float mp = 1.67262171e-24;      // proton mass [g]
    float rhoConstant = OmegaBaryonNow/OmegaMatterNow*DensityUnits;
    float HIIConstant = rhoConstant*InitialFractionHII;
    float HIConstant = rhoConstant - HIIConstant;
    float DeConstant = HIIConstant;
    float eUnits = VelocityUnits*VelocityUnits;
    float EUnits = DensityUnits*eUnits;
    // initialize constant field
    for (i=0; i<size; i++) {
      BaryonField[TENum][i]    = TEConstant/eUnits;
      BaryonField[V0Num][i]    = VxConstant/VelocityUnits;
      BaryonField[V1Num][i]    = VyConstant/VelocityUnits;
      BaryonField[V2Num][i]    = VzConstant/VelocityUnits;
      BaryonField[EgNum][i]    = EgConstant/EUnits;
      BaryonField[RhoNum][i]   = rhoConstant/DensityUnits;
      BaryonField[DeNum][i]    = DeConstant/DensityUnits;
      BaryonField[HINum][i]    = HIConstant/DensityUnits;
      BaryonField[HIINum][i]   = HIIConstant/DensityUnits;
    }
    if (DualEnergyFormalism)
      for (i=0; i<size; i++)
	BaryonField[GENum][i] = IEConstant/eUnits;

    if (debug) {
      printf("\n  Initializing constant fields using CGS values:\n");
      printf("        density = %g\n",rhoConstant);
      printf("   total energy = %g\n",TEConstant);
      if (DualEnergyFormalism)
	printf("   internal energy = %g\n",IEConstant);
      printf("     x-velocity = %g\n",VxConstant);
      printf("     y-velocity = %g\n",VyConstant);
      printf("     z-velocity = %g\n",VzConstant);
      printf("      radiation = %g\n",EgConstant);
      printf("      electrons = %g\n",DeConstant);
      printf("            nHI = %g\n",HIConstant);
      printf("           nHII = %g\n",HIIConstant);
      
      printf("Corresponding scaled values:\n");
      printf("        density = %g\n",rhoConstant/DensityUnits);
      printf("   total energy = %g\n",TEConstant/eUnits);
      if (DualEnergyFormalism)
	printf("   internal energy = %g\n",IEConstant/eUnits);
      printf("     x-velocity = %g\n",VxConstant/VelocityUnits);
      printf("     y-velocity = %g\n",VyConstant/VelocityUnits);
      printf("     z-velocity = %g\n",VzConstant/VelocityUnits);
      printf("      radiation = %g\n",EgConstant/EUnits);
      printf("      electrons = %g\n",DeConstant/DensityUnits);
      printf("            nHI = %g\n",HIConstant/DensityUnits);
      printf("           nHII = %g\n",HIIConstant/DensityUnits);
    }

  } // end if NewData == TRUE

  return SUCCESS;

#else

  fprintf(stderr,"Error: TRANSFER must be enabled for this test!\n");
  return FAIL;

#endif

}
