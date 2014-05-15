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
					float HydrogenMassFraction,
					float InitialFractionHII, 
					float InitialFractionHeII, 
					float InitialFractionHeIII, 
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
    HINum, HIINum, HeINum, HeIINum, HeIIINum, kphHINum, kphHeINum, 
    kphHeIINum, gammaNum, kdissH2INum, etaNum;
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
  if ((NumChemicals == 3) || (MultiSpecies > 0)) {
    FieldType[HeINum   = NumberOfBaryonFields++] = HeIDensity;
    FieldType[HeIINum  = NumberOfBaryonFields++] = HeIIDensity;    
    FieldType[HeIIINum = NumberOfBaryonFields++] = HeIIIDensity;
  }
  // if using external chemistry/cooling, set rate fields
  if (RadiativeCooling) {
    FieldType[kphHINum = NumberOfBaryonFields++] = kphHI;
    FieldType[gammaNum = NumberOfBaryonFields++] = PhotoGamma;
    if (RadiativeTransferHydrogenOnly == FALSE) {
      FieldType[kphHeINum  = NumberOfBaryonFields++] = kphHeI;
      FieldType[kphHeIINum = NumberOfBaryonFields++] = kphHeII;
    }
    if (MultiSpecies > 1)
      FieldType[kdissH2INum = NumberOfBaryonFields++] = kdissH2I;
  }
  // if using the AMRFLDSplit solver, set a field for the emissivity
  if (ImplicitProblem == 6) 
    FieldType[etaNum = NumberOfBaryonFields++] = Emissivity0;


  // set the subgrid static flag (necessary??)
  SubgridsAreStatic = FALSE;

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
  if (debug && NewData) {
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
    float HIIConstant = rhoConstant*InitialFractionHII*HydrogenMassFraction;
    float HIConstant = rhoConstant*HydrogenMassFraction - HIIConstant;
    float HeIIConstant = rhoConstant*InitialFractionHeII*(1.0-HydrogenMassFraction);
    float HeIIIConstant = rhoConstant*InitialFractionHeIII*(1.0-HydrogenMassFraction);
    float HeIConstant = (1.0-HydrogenMassFraction)*rhoConstant - HeIIConstant - HeIIIConstant;
    float DeConstant = HIIConstant + 0.25*HeIIConstant + 0.5*HeIIIConstant;
    float eUnits = VelocityUnits*VelocityUnits;
    float EUnits = DensityUnits*eUnits;
    // initialize constant fields
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
      if ((NumChemicals == 3) || (MultiSpecies > 0)) {
	BaryonField[HeINum][i]   = HeIConstant/DensityUnits;
	BaryonField[HeIINum][i]  = HeIIConstant/DensityUnits;
	BaryonField[HeIIINum][i] = HeIIIConstant/DensityUnits;
      }
    }
    if (DualEnergyFormalism)
      for (i=0; i<size; i++)
	BaryonField[GENum][i] = IEConstant/eUnits;
    // if using external chemistry/cooling, set rate fields
    if (RadiativeCooling) {
      for (i=0; i<size; i++)  BaryonField[kphHINum][i] = 0.0;
      for (i=0; i<size; i++)  BaryonField[gammaNum][i] = 0.0;
      if (RadiativeTransferHydrogenOnly == FALSE) {
	for (i=0; i<size; i++)  BaryonField[kphHeINum][i]  = 0.0;
	for (i=0; i<size; i++)  BaryonField[kphHeIINum][i] = 0.0;
      }
      if (MultiSpecies > 1)
	for (i=0; i<size; i++)  BaryonField[kdissH2INum][i] = 0.0;
    }

    if (debug && NewData) {
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
      if ((NumChemicals == 3) || (MultiSpecies > 0)) {
	printf("           nHeI = %g\n",HeIConstant);
	printf("          nHeII = %g\n",HeIIConstant);
	printf("         nHeIII = %g\n",HeIIIConstant);
      }
      
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
      if ((NumChemicals == 3) || (MultiSpecies > 0)) {
	printf("           nHeI = %g\n",BaryonField[HeINum][1]);
	printf("          nHeII = %g\n",BaryonField[HeIINum][1]);
	printf("         nHeIII = %g\n",BaryonField[HeIIINum][1]);
      }
    }

  } // end if NewData == TRUE

  return SUCCESS;

#else

  fprintf(stderr,"Error: TRANSFER must be enabled for this test!\n");
  return FAIL;

#endif

}
