/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE FREE-STREAMING RADIATION TEST) 
/
/  written by: Daniel Reynolds
/  date:       June 2009
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



int grid::FSMultiSourceInitializeGrid(float DensityConstant, 
				      float VxConstant, 
				      float VyConstant, 
				      float VzConstant, 
				      float TEConstant, 
				      float RadConstant, 
				      int   local)
{
#ifdef TRANSFER
//   if (debug)
//     fprintf(stdout,"Entering grid::FSMultiSourceInitializeGrid routine\n");

  // determine whether data should be allocated/initialized or not
  int NewData = TRUE;
  if ((ParallelRootGridIO == TRUE) && (local == 0))
    NewData = FALSE;


  // create necessary baryon fields
  int RhoNum, TENum, IENum, V0Num, V1Num, V2Num, RadNum;
  NumberOfBaryonFields = 0;
  FieldType[RhoNum = NumberOfBaryonFields++] = Density;
  FieldType[TENum = NumberOfBaryonFields++]  = TotalEnergy;
  FieldType[V0Num = NumberOfBaryonFields++]  = Velocity1;
  FieldType[V1Num = NumberOfBaryonFields++]  = Velocity2;
  FieldType[V2Num = NumberOfBaryonFields++]  = Velocity3;
  FieldType[RadNum = NumberOfBaryonFields++] = RadiationFreq0;

  // set the subgrid static flag (necessary??)
  SubgridsAreStatic = FALSE;  // no subgrids

  // Return if this doesn't concern us.
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  // Get various units
  double MassUnits = 1.0;
  float DensityUnits=1.0, LengthUnits=1.0, TemperatureUnits=1.0, 
    TimeUnits=1.0, VelocityUnits=1.0;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
    fprintf(stderr,"Error in GetUnits.\n");
    return FAIL;
  }
  if (debug) {
    fprintf(stdout,"  Internal Unit Conversion Factors:\n");
    fprintf(stdout,"         length = %g\n",LengthUnits);
    fprintf(stdout,"           mass = %g\n",MassUnits);
    fprintf(stdout,"           time = %g\n",TimeUnits);
  }

  // allocate fields
  if (NewData == TRUE) {

    // compute size of fields
    int size = 1;
    for (int dim=0; dim<GridRank; dim++)  size *= GridDimension[dim];
    
    // allocate the fields
    this->AllocateGrids();
    
    // set fluid density, total energy, [internal energy,] velocities, 
    // radiation energy, electron density, chemical species
    int i;
    float eUnits = VelocityUnits*VelocityUnits;
    float rUnits = DensityUnits*eUnits;
    for (i=0; i<size; i++) {
      BaryonField[RhoNum][i] = DensityConstant/DensityUnits;
      BaryonField[TENum][i]  = TEConstant/eUnits;
      BaryonField[V0Num][i]  = VxConstant/VelocityUnits;
      BaryonField[V1Num][i]  = VyConstant/VelocityUnits;
      BaryonField[V2Num][i]  = VzConstant/VelocityUnits;
      BaryonField[RadNum][i] = RadConstant/rUnits;
    }
    
    if (debug) {
      printf("\n  Initializing constant fields using CGS values:\n");
      printf("        density = %g\n",DensityConstant);
      printf("   total energy = %g\n",TEConstant);
      printf("     x-velocity = %g\n",VxConstant);
      printf("     y-velocity = %g\n",VyConstant);
      printf("     z-velocity = %g\n",VzConstant);
      printf("      radiation = %g\n",RadConstant);
      
      printf("Corresponding scaled values:\n");
      printf("        density = %g\n",DensityConstant/DensityUnits);
      printf("   total energy = %g\n",TEConstant/eUnits);
      printf("     x-velocity = %g\n",VxConstant/VelocityUnits);
      printf("     y-velocity = %g\n",VyConstant/VelocityUnits);
      printf("     z-velocity = %g\n",VzConstant/VelocityUnits);
      printf("      radiation = %g\n",RadConstant/rUnits);
    }

  } // end if NewData == TRUE

  return SUCCESS;

#else

  fprintf(stderr,"Error: TRANSFER must be enabled for this test!\n");
  return FAIL;

#endif

}
