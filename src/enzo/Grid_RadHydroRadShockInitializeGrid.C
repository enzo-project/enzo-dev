/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE PROTOSTELLAR CORE COLLAPSE TEST) 
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  John Hayes, July 2007; cloned for the radiation shock problem.
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


int grid::RadHydroRadShockInitializeGrid(float DensityConstant, 
					 float TEConstant, 
					 float REConstant,
					 float VelocityConstant,
                                         int   ShockDir,
					 int   local)
{

#ifdef TRANSFER
//   if (MyProcessorNumber == ROOT_PROCESSOR)
//     fprintf(stdout,"Entering grid::RadHydroRadShockInitializeGrid routine\n");
 
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
  FieldType[RhoNum = NumberOfBaryonFields++]   = Density;
  FieldType[TENum = NumberOfBaryonFields++]    = TotalEnergy;
  if (DualEnergyFormalism) 
    FieldType[IENum = NumberOfBaryonFields++]  = InternalEnergy;
  FieldType[V0Num = NumberOfBaryonFields++]    = Velocity1;
  FieldType[V1Num = NumberOfBaryonFields++]    = Velocity2;
  FieldType[V2Num = NumberOfBaryonFields++]    = Velocity3;
  FieldType[EgNum = NumberOfBaryonFields++]    = RadiationFreq0;
  
  // set the subgrid static flag (necessary??)
  SubgridsAreStatic = FALSE;  // no subgrids
 
  // Return if this doesn't concern us.
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  // Get various units
  float DensityUnits=1, LengthUnits=1, TemperatureUnits=1, TimeUnits=1,
    VelocityUnits=1;
  double MassUnits=1;
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
    this->AllocateGrids();
 
    // set fluid density, total energy, [internal energy,] velocities, 
    // radiation energy
    int i, j, k;
    float IEConstant = TEConstant - 0.5 * VelocityConstant * VelocityConstant;
    float eUnits = VelocityUnits*VelocityUnits;
    float EUnits = DensityUnits*eUnits;
    for (i=0; i<size; i++) {
      BaryonField[RhoNum][i] = DensityConstant/DensityUnits;
      BaryonField[TENum][i]  = TEConstant/eUnits;
      BaryonField[EgNum][i]  = REConstant/EUnits;
    }
    
    if ( ShockDir == 0 )
      for (i=0; i<size; i++) {
        BaryonField[V0Num][i]  = VelocityConstant/VelocityUnits;
        BaryonField[V1Num][i]  = 0.0;
        BaryonField[V2Num][i]  = 0.0;
      }

    if ( ShockDir == 1 ) 
      for (i=0; i<size; i++) {
        BaryonField[V0Num][i]  = 0.0;
        BaryonField[V1Num][i]  = VelocityConstant/VelocityUnits;
        BaryonField[V2Num][i]  = 0.0;
      }

    if ( ShockDir == 2 ) 
      for (i=0; i<size; i++) {
        BaryonField[V0Num][i]  = 0.0;
        BaryonField[V1Num][i]  = 0.0;
        BaryonField[V2Num][i]  = VelocityConstant/VelocityUnits;
      }

    if (DualEnergyFormalism)
      for (i=0; i<size; i++)
	BaryonField[IENum][i] = IEConstant/eUnits;

    if (debug) {
      fprintf(stdout,"RadHydroRadShockInitializeGrid:\n");
      printf("           ShockDir = %"ISYM"\n",ShockDir);
      
      printf("    DensityConstant = %g\n",DensityConstant);    
      printf("         TEConstant = %g\n",TEConstant);    
      if (DualEnergyFormalism)
	printf("         IEConstant = %g\n",IEConstant);    
      printf("         REConstant = %g\n",REConstant);    
      
      printf("Corresponding scaled values:\n");
      printf("    Density = %g\n",DensityConstant/DensityUnits);
      printf("         TE = %g\n",TEConstant/eUnits);
      if (DualEnergyFormalism)
	printf("         IE = %g\n",IEConstant/eUnits);
      printf("         RE = %g\n",REConstant/EUnits);
    }

  } // end if NewData == TRUE

  return SUCCESS;

#else

  fprintf(stderr,"Error: TRANSFER must be enabled for this test!\n");
  return FAIL;

#endif

}
