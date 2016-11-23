/***********************************************************************
/
/  GRID CLASS (FLAG CELLS TO BE REFINED BY THE JEAN'S CRITERION)
/
/  written by: Greg Bryan
/  date:       February, 1998
/  modified1:  Alexei Kritsuk, Feb. 2004 mods for isothermal EOS.
/
/  PURPOSE:
/
/  RETURNS:
/    number of flagged cells, or -1 on failure
/
************************************************************************/
 
#include <stdio.h>
#include <math.h>
#include "phys_constants.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "hydro_rk/EOS.h"
 
/* function prototypes */
 
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int QuantumGetUnits (float *DensityUnits, float *LengthUnits,
        float *TemperatureUnits, float *TimeUnits,
        float *VelocityUnits, double *MassUnits, FLOAT Time);
 
int grid::FlagCellsToBeRefinedByVelDiv()
{
  /* declarations */
 
  int i, j, k, dim;
  float *u, *v, *w;
  float dxcalc;
 
  /* error check */
 
  if (FlaggingField == NULL) {
    fprintf(stderr, "Flagging Field is undefined.\n");
    return -1;
  }

   /* Get density units. */

  float DensityUnits=1, LengthUnits=1, VelocityUnits=1, TimeUnits=1,
    TemperatureUnits=1, MassUnits=1;
  /* Get quantum coef. */
  if (QuantumGetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.");
  }
 
   FLOAT hmcoef = 5.9157166856e27*TimeUnits/pow(LengthUnits,2)/FDMMass;

   /* Compute expansion factor*/
   FLOAT a = 1.0, dadt;
    if (ComovingCoordinates){
      if (CosmologyComputeExpansionFactor(Time, &a, &dadt) == FAIL) {
    ENZO_FAIL("Error in ComputeExpansionFactor.\n");
     }
    }
 
  /* compute size */
 
  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* Find fields: density, total energy, velocity1-3. */
 
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }

  u = BaryonField[Vel1Num];
  v = BaryonField[Vel2Num];
  w = BaryonField[Vel3Num];




    /* allocate divergence of velocity */
  float *div = new float[size];
  int in = GridDimension[0];
  int jn = GridDimension[1];
  int kn = GridDimension[2];

  float dx = CellWidth[0][0];

  if (GridRank < 2) {
    v = new float[size];
    for (i = 0; i < size; i++)
      v[i] = 0;
  }
  if (GridRank < 3) {
    w = new float[size];
    for (i = 0; i < size; i++)
      w[i] = 0;
  }



  /* loop to compute velocity divergence */

  #define IDX(a,b,c) ( ((c)*jn + (b))*in + (a) )

  for (k=0; k<kn; k++){
  	for (j=0; j<jn; j++){
  		for (i=0; i<in; i++){
  			div[IDX(i,j,k)] = fabs((u[IDX(min(i+1,in-1),j,k)]-u[(IDX(i,j,k))]))/dx;
  			if (GridRank>1){
  				div[IDX(i,j,k)] += fabs((v[IDX(i,min(j+1,jn-1),k)]-v[(IDX(i,j,k))]))/dx;
  			}
  			if (GridRank>2){
  				div[IDX(i,j,k)] += fabs((w[IDX(i,j,min(k+1,kn-1))]-w[(IDX(i,j,k))]))/dx;
  			}
  			div[IDX(i,j,k)] = max(div[IDX(i,j,k)]*RefineByVelDivSafetyFactor, tiny_number);


  		}
  	}
  }

 
  /* Loop over grid. */
 
#ifdef UNUSED
  int j, k, index;
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++){
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      index = (j + k*GridDimension[1])*GridDimension[0];
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++){
              dxcalc = sqrt(hmcoef/a/div[i]);
	if (dx > dxcalc){
	  FlaggingField[index+i]++;}
        }
      }
    }
#endif /* UNUSED */
 
 // FLOAT CellWidthSquared = CellWidth[0][0]*CellWidth[0][0];
  for (i = 0; i < size; i++)
    {
      dxcalc = sqrt(hmcoef/a/div[i]);
	if (dx > dxcalc){
	  FlaggingField[i]++; 

      }    

    }
 
  /* clean up */
 
  /* Count number of flagged Cells. */
 
  int NumberOfFlaggedCells = 0;
  for (i = 0; i < size; i++) {
    FlaggingField[i] = (FlaggingField[i] >= 1)? 1 : 0;
    NumberOfFlaggedCells += FlaggingField[i];
  }

  delete [] div;
  if (GridRank < 2) delete [] v;
  if (GridRank < 3) delete [] w;
 
  return NumberOfFlaggedCells;
 
}
