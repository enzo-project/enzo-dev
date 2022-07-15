/***********************************************************************
/
/  GRID CLASS (CHECK THE ACCELERATION FIELD)
/
/  written by: JC Passy
/  date:       March, 2013
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
 
int grid::OutputAccelerationField(FILE *fptr, int level)
{
 
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
 
  /* declarations */
 
  int dim,size=1;
  float dist[MAX_DIMENSION], Middle[MAX_DIMENSION], Width[MAX_DIMENSION];
  int i,j,k,index;

  /* Set origins. */
 
  for (dim = 0; dim < GridRank; dim++) {
    Middle[dim] = 0.5*(DomainLeftEdge[dim] + DomainRightEdge[dim]);
    Width[dim]  =      DomainRightEdge[dim] - DomainLeftEdge[dim] ;
  }

  /* Compute field size (in floats). */

  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  /* Diagnostic grid */
  
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
  fprintf(fptr,"# level, is_ghost_zone, i, j, k, x, y, z, r, ax, ay, az, atan, arad\n");

  /* Sanity check */

  for (dim = 0; dim < GridRank; dim++)
    //if (AccelerationField[dim] == NULL) {
    if (!AccelerationField[dim]) {
      fprintf(fptr,"AccelerationField is NULL!\n");
      return SUCCESS;
    }
  //ENZO_FAIL("Error in grid->OutputAccelerationField: AccelerationField is NULL!\n");

  /* Loop over all grid cells */
  
  float xpos,ypos,zpos,rpos,ax,ay,az,arad,atang,a2;
  int is_ghost_zone;  

  for (k = 0; k < GridDimension[2]; k++) {
    for (j = 0; j < GridDimension[1]; j++) {
      for (i = 0; i < GridDimension[0]; i++) {

	zpos = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
	ypos = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	xpos = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
    
	is_ghost_zone = 0;

	// Ghost zones
	if ((k < GridStartIndex[2]) || (k > GridEndIndex[2]) || 
	    (j < GridStartIndex[1]) || (j > GridEndIndex[1]) || 
	    (i < GridStartIndex[0]) || (i > GridEndIndex[0]))
	  
	  is_ghost_zone = 1;

	// For Level 0, flag an extra zone at the boundary of the domain
	if (level == 0)
	  if (
	      (zpos < DomainLeftEdge[2]+CellWidth[2][0]) || (zpos > DomainRightEdge[2]-CellWidth[2][0]) || 
	      (ypos < DomainLeftEdge[1]+CellWidth[1][0]) || (ypos > DomainRightEdge[1]-CellWidth[1][0]) || 
	      (xpos < DomainLeftEdge[0]+CellWidth[0][0]) || (xpos > DomainRightEdge[0]-CellWidth[0][0])
	      )
	  
	    is_ghost_zone = 1;

	// Distance form the center
	zpos -= Middle[2];
	ypos -= Middle[1];
	xpos -= Middle[0];

	index = i + GridDimension[0]*(j + GridDimension[1]*k);
	
	/* Calculate force components */	
	// With Zeus, it won't be perfect because of the interpolation
	if (HydroMethod == Zeus_Hydro && GravitySolverType == GRAVITY_SOLVER_FAST) {
	  
	  if (i < GridEndIndex[0])
	    ax = 0.5*(AccelerationField[0][index] + 
		      AccelerationField[0][index+1]);
	  else
	    ax = AccelerationField[0][index];

	  if (j < GridEndIndex[1])
	    ay = 0.5*(AccelerationField[1][index] + 
		      AccelerationField[1][index+GridDimension[0]]);
	  else
	    ay = AccelerationField[1][index];

	  if (k < GridEndIndex[2])
	    az = 0.5*(AccelerationField[2][index] + 
		      AccelerationField[2][index+GridDimension[0]*GridDimension[1]]);
	  else
	    az = AccelerationField[2][index];
	  
	} else {
	  
	  ax = AccelerationField[0][index];
	  ay = AccelerationField[1][index];
	  az = AccelerationField[2][index];

	}
	
	// If required, add external field for the APM solver
	if (GravitySolverType == GRAVITY_SOLVER_APM) 
	    if (ProblemType == 41 || ProblemType == 46) {
	      ax += AccelerationFieldExternalAPM[0][index];
	      ay += AccelerationFieldExternalAPM[1][index];
	      az += AccelerationFieldExternalAPM[2][index];
	    }
	
	rpos = pow(xpos*xpos + ypos*ypos + zpos*zpos, 0.5);
	
	arad = (ax*xpos + ay*ypos + az*zpos) / rpos;
	a2 = pow(ax,2.0) + pow(ay,2.0) + pow(az,2.0);
	atang = sqrt(max(a2-arad*arad, 0.0));

	/* Output results. */
	
	fprintf(fptr, "%"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM"\n",
		level,is_ghost_zone,
		i,j,k,
		xpos,ypos,zpos,rpos,
		ax, ay, az,
		atang, -arad);

      }
    }
  } //end loop grid dims

  return SUCCESS;

}
  
