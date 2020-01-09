/***********************************************************************
/
/  Check GravitatingMassField
/
/  written by: Passy JC
/  date:       April, 2013
/  modified1:
/
/           
/
************************************************************************/

#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "phys_constants.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
       	      float *TemperatureUnits, float *TimeUnits,
       	      float *VelocityUnits, double *MassUnits, FLOAT Time);

int grid::OutputGravitatingMassField(FILE *fptr, FILE *fptr2, int level)
{

  /* Return if this grid is not on this processor. */

  if (MyProcessorNumber != ProcessorNumber)
   return SUCCESS;

  /* declarations */

  int i,j,k,index,dim;
  int size_gmf = 1;
  int size_gmf2 = 1;
  int is_ghost_zone;

  float sum_gmf, sum_gmfp;

  /* Compute field size */

  for (dim = 0; dim < GridRank; dim++)
    size_gmf *= GravitatingMassFieldDimension[dim];

  for (dim = 0; dim < GridRank; dim++)
    size_gmf2 *= GravitatingMassFieldParticlesDimension[dim];


  /* Sanity check */

  if (GravitatingMassField == NULL)
    ENZO_FAIL("Error in grid->OutputGravitatingMassField: GravitatingMassField is NULL!\n");

  if (GravitatingMassFieldParticles == NULL)
    ENZO_FAIL("Error in grid->OutputGravitatingMassField: GravitatingMassFieldParticles is NULL!\n");

  /* Headers */

  fprintf(fptr,"# level = %"ISYM", time = %"ESYM"\n", level, Time);
  fprintf(fptr,"# GridLeftEdge = %"FSYM", %"FSYM", %"FSYM" \n",
	  GridLeftEdge[0],GridLeftEdge[1],GridLeftEdge[2]);
  fprintf(fptr,"# GridRightEdge = %"FSYM", %"FSYM", %"FSYM" \n",
	  GridRightEdge[0],GridRightEdge[1],GridRightEdge[2]);
  fprintf(fptr,"# GridStartIndex = %"ISYM" %"ISYM" %"ISYM"\n",
	  GridStartIndex[0],GridStartIndex[1],GridStartIndex[2]);
  fprintf(fptr,"# GridEndIndex = %"ISYM" %"ISYM" %"ISYM"\n",
	  GridEndIndex[0],GridEndIndex[1],GridEndIndex[2]);
  fprintf(fptr,"# GridDimension = %"ISYM" %"ISYM" %"ISYM"\n",
	  GridDimension[0],GridDimension[1],GridDimension[2]);
  fprintf(fptr,"# GravitatingMassFieldLeftEdge = %"FSYM", %"FSYM", %"FSYM" \n",
	  GravitatingMassFieldLeftEdge[0],GravitatingMassFieldLeftEdge[1],GravitatingMassFieldLeftEdge[2]);
  fprintf(fptr, "# size of GravitatingMassField = %"ISYM" = %"ISYM"x%"ISYM"x%"ISYM"\n",
	  size_gmf,
	  GravitatingMassFieldDimension[0],
	  GravitatingMassFieldDimension[1],
	  GravitatingMassFieldDimension[2]);
  fprintf(fptr,"# level, i, j, k, index, GMF[index]\n");
  fflush(fptr);

  //

  fprintf(fptr2,"# level = %"ISYM", time = %"ESYM"\n", level, Time);
  fprintf(fptr2,"# GridLeftEdge = %"FSYM", %"FSYM", %"FSYM" \n",
	  GridLeftEdge[0],GridLeftEdge[1],GridLeftEdge[2]);
  fprintf(fptr2,"# GridRightEdge = %"FSYM", %"FSYM", %"FSYM" \n",
	  GridRightEdge[0],GridRightEdge[1],GridRightEdge[2]);
  fprintf(fptr2,"# GridStartIndex = %"ISYM" %"ISYM" %"ISYM"\n",
	  GridStartIndex[0],GridStartIndex[1],GridStartIndex[2]);
  fprintf(fptr2,"# GridEndIndex = %"ISYM" %"ISYM" %"ISYM"\n",
	  GridEndIndex[0],GridEndIndex[1],GridEndIndex[2]);
  fprintf(fptr2,"# GridDimension = %"ISYM" %"ISYM" %"ISYM"\n",
	  GridDimension[0],GridDimension[1],GridDimension[2]);
  fprintf(fptr2,"# GravitatingMassFieldLeftEdge = %"FSYM", %"FSYM", %"FSYM" \n",
	  GravitatingMassFieldLeftEdge[0],GravitatingMassFieldLeftEdge[1],GravitatingMassFieldLeftEdge[2]);
  fprintf(fptr2, "# size of GravitatingMassField = %"ISYM" = %"ISYM"x%"ISYM"x%"ISYM"\n",
	  size_gmf2,
	  GravitatingMassFieldParticlesDimension[0],
	  GravitatingMassFieldParticlesDimension[1],
	  GravitatingMassFieldParticlesDimension[2]);
  fprintf(fptr2,"# level, i, j, k, index, GMF[index]\n");
  fflush(fptr2);


  sum_gmf = 0.0;
  sum_gmfp = 0.0;

  // GMF
  for (k = 0; k < GravitatingMassFieldDimension[2]; k++) {    
    for (j = 0; j < GravitatingMassFieldDimension[1]; j++) {
      for (i = 0; i < GravitatingMassFieldDimension[0]; i++) {


	index = i + GravitatingMassFieldDimension[0]*(j + GravitatingMassFieldDimension[1]*k);
	sum_gmf  += GravitatingMassField[index];
	sum_gmfp  += GravitatingMassFieldParticles[index];

	is_ghost_zone = 0;
	// Ghost zones
	if ((k < GridStartIndex[2]) || (k > GridEndIndex[2]) ||
            (j < GridStartIndex[1]) || (j > GridEndIndex[1]) ||
            (i < GridStartIndex[0]) || (i > GridEndIndex[0]))

          is_ghost_zone = 1;
	
	// For Level 0, flag an extra zone
	if (level == 0)
          if ((k == GridStartIndex[2]) || (k == GridEndIndex[2]) ||
              (j == GridStartIndex[1]) || (j == GridEndIndex[1]) ||
              (i == GridStartIndex[0]) || (i == GridEndIndex[0]))

            is_ghost_zone = 1;
	
	// GMF
	if (GravitatingMassField[index] > tiny_number) {
	  fprintf(fptr,"%"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ESYM"\n", 
		  level, i, j, k, index, GravitatingMassField[index]);
	  fflush(fptr);
	}

	// GMFP
	if (GravitatingMassFieldParticles[index] > tiny_number) {
	  fprintf(fptr2,"%"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ESYM"\n", 
		  level, i, j, k, index, GravitatingMassFieldParticles[index]);
	  fflush(fptr2);
	}


      }
    }
  }

  // Sum
  fprintf(fptr,"# total gmf: %"ESYM"\n",sum_gmf);
  fflush(fptr);
  fprintf(fptr2,"# total gmfp: %"ESYM"\n",sum_gmfp);
  fflush(fptr2);

  return SUCCESS;

}
