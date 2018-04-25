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


#define DEFAULT_MU 0.6   // mean molecular mass for temperature field


// function prototypes
int GetUnits(float *DensityUnits, float *LengthUnits, 
	     float *TemperatureUnits, float *TimeUnits, 
	     float *VelocityUnits, double *MassUnits, FLOAT Time);



int grid::RadHydroStreamTestInitializeGrid(float DensityConstant, 
					   float EgConstant, 
					   int RadStreamDim,
					   int RadStreamDir,
					   int local)
{

#ifdef TRANSFER
//   if (debug)
//     fprintf(stdout,"Entering grid::RadHydroStreamTestInitializeGrid routine\n");

  // determine whether data should be allocated/initialized or not
  int NewData = TRUE;
  if ((ParallelRootGridIO == TRUE) && (local == 0))
    NewData = FALSE;

  // if grids allocated and already set up (i.e. restart), return
  if ((NumberOfBaryonFields > 5) && (BaryonField[5] != NULL))
    return SUCCESS;

  // create necessary baryon fields
  int RhoNum, TENum, IENum, V0Num, V1Num, V2Num, EgNum;
  NumberOfBaryonFields = 0;
  FieldType[RhoNum = NumberOfBaryonFields++] = Density;
  FieldType[TENum = NumberOfBaryonFields++]  = TotalEnergy;
  if (DualEnergyFormalism) 
    FieldType[IENum = NumberOfBaryonFields++] = InternalEnergy;
  FieldType[V0Num = NumberOfBaryonFields++] = Velocity1;
  FieldType[V1Num = NumberOfBaryonFields++] = Velocity2;
  FieldType[V2Num = NumberOfBaryonFields++] = Velocity3;
  FieldType[EgNum = NumberOfBaryonFields++] = RadiationFreq0;
  
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
    int i, j, k;
    float pi = 4.0*atan(1.0);
    float StBz = 5.6704e-5;
    float c = 2.99792458e10;
    float IEConstant = 1.0/(Gamma-1.0)/DEFAULT_MU*sqrt(sqrt((0.25*c*EgConstant/StBz)));
    float TEConstant = IEConstant;
    float eUnits = VelocityUnits*VelocityUnits;
    float EUnits = DensityUnits*eUnits;
    for (i=0; i<size; i++) {
      BaryonField[RhoNum][i] = DensityConstant/DensityUnits;
      BaryonField[TENum][i]  = TEConstant/eUnits;
      BaryonField[V0Num][i]  = 0.0;
      BaryonField[V1Num][i]  = 0.0;
      BaryonField[V2Num][i]  = 0.0;
      BaryonField[EgNum][i]  = EgConstant/EUnits;
    }
    if (DualEnergyFormalism)
      for (i=0; i<size; i++)
	BaryonField[IENum][i] = IEConstant/eUnits;

    if (debug) {
      printf("RadHydroStreamTestInitializeGrid:\n");
      printf("       RadStreamDim = %"ISYM"\n",MyProcessorNumber,RadStreamDim);
      printf("       RadStreamDir = %"ISYM"\n",MyProcessorNumber,RadStreamDir);
      printf("     MaxRadiationDt = %g\n",MyProcessorNumber,MaxRadiationDt);

      printf("    DensityConstant = %g\n",DensityConstant);    
      printf("         TEConstant = %g\n",TEConstant);    
      if (DualEnergyFormalism)
	printf("   internal energy = %g\n",IEConstant);
      printf("         EgConstant = %g\n",EgConstant);    
      
      printf("Corresponding scaled values:\n");
      printf("    Density = %g\n",DensityConstant/DensityUnits);
      printf("         TE = %g\n",TEConstant/eUnits);
      if (DualEnergyFormalism)
	printf("   internal energy = %g\n",IEConstant/eUnits);
      printf("         Eg = %g\n",EgConstant/EUnits);
    }

    // adjust Radiation Energy Density BC at streaming input edge
    int idx;
    if (RadStreamDim == 0) {
      if (RadStreamDir == 0) {
	if (GridLeftEdge[0] == DomainLeftEdge[0]) {
	  for (k=GridStartIndex[2]; k<=GridEndIndex[2]; k++)
	    for (j=GridStartIndex[1]; j<=GridEndIndex[1]; j++)
	      for (i=0; i<GridStartIndex[0]; i++) {
		idx = (k*GridDimension[1] + j)*GridDimension[0] + i;
		BaryonField[EgNum][idx] = 1.0;
	      }
	}
      }
      else if (RadStreamDir == 1) {
	if (GridRightEdge[0] == DomainRightEdge[0]) {
	  for (k=GridStartIndex[2]; k<=GridEndIndex[2]; k++)
	    for (j=GridStartIndex[1]; j<=GridEndIndex[1]; j++)
	      for (i=GridEndIndex[0]+1; i<GridDimension[0]; i++) {
		idx = (k*GridDimension[1] + j)*GridDimension[0] + i;
		BaryonField[EgNum][idx] = 1.0;
	      }
	}
      }
      else {
	fprintf(stderr,"RHStreamTest Error: illegal direction = %"ISYM"\n",
		RadStreamDir);
	return FAIL;
      }
    }
    else if (RadStreamDim == 1) {
      if (RadStreamDir == 0) {
	if (GridLeftEdge[1] == DomainLeftEdge[1]) {
	  for (k=GridStartIndex[2]; k<=GridEndIndex[2]; k++)
	    for (j=0; j<GridStartIndex[1]; j++)
	      for (i=GridStartIndex[0]; i<=GridEndIndex[0]; i++) {
		idx = (k*GridDimension[1] + j)*GridDimension[0] + i;
		BaryonField[EgNum][idx] = 1.0;
	      }
	}
      }
      else if (RadStreamDir == 1) {
	if (GridRightEdge[1] == DomainRightEdge[1]) {
	  for (k=GridStartIndex[2]; k<=GridEndIndex[2]; k++)
	    for (j=GridEndIndex[1]+1; j<GridDimension[1]; j++)
	      for (i=GridStartIndex[0]; i<=GridEndIndex[0]; i++) {
		idx = (k*GridDimension[1] + j)*GridDimension[0] + i;
		BaryonField[EgNum][idx] = 1.0;
	      }
	}
      }
      else {
	fprintf(stderr,"RHStreamTest Error: illegal direction = %"ISYM"\n",
		RadStreamDir);
	return FAIL;
      }
    }
    else if (RadStreamDim == 2) {
      if (RadStreamDir == 0) {
	if (GridLeftEdge[2] == DomainLeftEdge[2]) {
	  for (k=0; k<GridStartIndex[2]; k++)
	    for (j=GridStartIndex[1]; j<=GridEndIndex[1]; j++)
	      for (i=GridStartIndex[0]; i<=GridEndIndex[0]; i++) {
		idx = (k*GridDimension[1] + j)*GridDimension[0] + i;
		BaryonField[EgNum][idx] = 1.0;
	      }
	}
      }
      else if (RadStreamDir == 1) {
	if (GridRightEdge[2] == DomainRightEdge[2]) {
	  for (k=GridEndIndex[2]+1; k<GridDimension[2]; k++)
	    for (j=GridStartIndex[1]; j<=GridEndIndex[1]; j++)
	      for (i=GridStartIndex[0]; i<=GridEndIndex[0]; i++) {
		idx = (k*GridDimension[1] + j)*GridDimension[0] + i;
		BaryonField[EgNum][idx] = 1.0;
	      }
	}
      }
      else {
	fprintf(stderr,"RHStreamTest Error: illegal direction = %"ISYM"\n",
		RadStreamDir);
	return FAIL;
      }
    }
    else {
      fprintf(stderr,"RHStreamTest Error: illegal dimension = %"ISYM"\n",
	      RadStreamDim);
      return FAIL;    
    }
  } // end if NewData == TRUE

  return SUCCESS;

#else

  fprintf(stderr,"Error: TRANSFER must be enabled for this test!\n");
  return FAIL;  
  
#endif

}
