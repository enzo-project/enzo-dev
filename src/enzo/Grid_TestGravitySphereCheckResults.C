/***********************************************************************
/
/  GRID CLASS (CHECK THE NUMERICAL GRAVITY FORCE AGAINST ANALYTIC)
/
/  written by: Greg Bryan
/  date:       September, 1995
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
#include "TestGravitySphereGlobalData.h"
 
float TestGravitySphereComputeRadialForce(float r, int GridRank);
 
int grid::TestGravitySphereCheckResults(FILE *fptr)
{
 
  char *TGSTangName     = "TangentialForce";
  char *TGSRadialName   = "RadialForce";
  char *TGSAnalyticName = "AnalyticForce";
 
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
 
  /* declarations */
 
  int dim, i;
  float dist[MAX_DIMENSION], Middle[MAX_DIMENSION];
  float r, fanalytic, fradial, ftang;
 
  /* Set constants. */
 
  for (dim = 0; dim < GridRank; dim++)
    Middle[dim] = 0.5*(DomainLeftEdge[dim] + DomainRightEdge[dim]);
  for (dim = GridRank; dim < MAX_DIMENSION; dim++)
    Middle[dim] = 0;
 
  /* Loop over particles, computing radial distance from the center and
     comparing the analyic force to the compute acceleration (determined
     by assuming the velocity is due soley to the analytic acceleration). */
 
  for (i = 0; i < NumberOfParticles; i++) {
 
    /* Compute distance and radial force. */
 
    r       = 0.0;
    fradial = 0.0;
    for (dim = 0; dim < GridRank; dim++) {
      dist[dim] = *(ParticlePosition[dim]+i) - Middle[dim];
      r += dist[dim]*dist[dim];
      fradial += (*(ParticleVelocity[dim] + i))/Time * dist[dim];
    }
 
    /* Normalize radial force components. */
 
    r = sqrt(r);
    fradial /= r;
 
    /* Compute tangential component. */
 
    ftang = 0.0;
    for (dim = 0; dim < GridRank; dim++)
      ftang += POW((*(ParticleVelocity[dim] + i))/Time, FLOAT(2));
    ftang = sqrt(max(ftang - fradial*fradial, tiny_number));
 
    /* Compute analytic acceleration. */
 
    fanalytic = TestGravitySphereComputeRadialForce(r, GridRank);
 
    /* Output results. */
 
    fprintf(fptr, "%"FSYM"  %e   %e   %e\n", r, ftang, -fradial, fanalytic);
 
  } // end loop over particles.
 
 
  /* If requested, set up the baryon field. */
 
  if (NumberOfBaryonFields > 0) {
 
    /* Create new fields to hold radial, tangential and analytic forces. */
 
    int size = 1;
    for (dim = 0; dim < GridRank; dim++)
      size *= GridDimension[dim];
 
    BaryonField[NumberOfBaryonFields++] = new float[size];
    BaryonField[NumberOfBaryonFields++] = new float[size];
    BaryonField[NumberOfBaryonFields++] = new float[size];
 
    DataLabel[NumberOfBaryonFields-3] = TGSTangName;
    DataLabel[NumberOfBaryonFields-2] = TGSRadialName;
    DataLabel[NumberOfBaryonFields-1] = TGSAnalyticName;
 
    float r, GridPos[MAX_DIMENSION];
    int n = 0;
 
    /* Loop over grid. */
 
    for (int k = 0; k < GridDimension[2]; k++)
      for (int j = 0; j < GridDimension[1]; j++)
	for (int i = 0; i < GridDimension[0]; i++, n++) {
 
	  /* Compute position */
 
	  GridPos[0] = *(CellLeftEdge[0]+i) + 0.5*(*(CellWidth[0]+i));
	  if (GridRank > 1)
	    GridPos[1] = *(CellLeftEdge[1]+j) + 0.5*(*(CellWidth[1]+j));
	  if (GridRank > 2)
	    GridPos[2] = *(CellLeftEdge[2]+k) + 0.5*(*(CellWidth[2]+k));
 
	  /* Compute distance and radial force. */
 
	  r       = 0.0;
	  fradial = 0.0;
	  for (dim = 0; dim < GridRank; dim++) {
	    dist[dim] = GridPos[dim] - Middle[dim];
	    r       += dist[dim]*dist[dim];
	    fradial += (*(BaryonField[2+dim] + n))/Time * dist[dim];
	  }
 
	  /* Normalize radial force components. */
 
	  r = sqrt(r);
	  fradial /= r;
 
	  /* Compute tangential component. */
 
	  ftang = 0.0;
	  for (dim = 0; dim < GridRank; dim++)
	    ftang += POW((*(BaryonField[2+dim] + n))/Time, FLOAT(2));
	  ftang = sqrt(max(ftang - fradial*fradial, tiny_number));
 
	  /* Compute analytic acceleration. */
 
	  fanalytic = TestGravitySphereComputeRadialForce(r, GridRank);
	
	  /* Output results. */
 
	  BaryonField[NumberOfBaryonFields-3][n] = ftang;
	  BaryonField[NumberOfBaryonFields-2][n] = -fradial;
	  BaryonField[NumberOfBaryonFields-1][n] = fanalytic;
 
	} // end: loop over grid
 
  }
 
  return SUCCESS;
}
 
 
 
 
 
/* This function compute the correct radial force given the radius and
   sphere type (and other parameters). */
 
float TestGravitySphereComputeRadialForce(float r, int GridRank)
{
 
  float fanalytic = 0;
 
  if (GridRank == 3) {
 
    /* 0) uniform sphere */
 
    if (TestGravitySphereType == 0) {
      if (r < TestGravitySphereRadius)
	fanalytic = 1.0/3.0*GravitationalConstant*
	  TestGravitySphereInteriorDensity*r;
      else
	fanalytic = 1.0/3.0*GravitationalConstant/(r*r)*
	  (TestGravitySphereInteriorDensity*POW(TestGravitySphereRadius, 3) +
	   TestGravitySphereExteriorDensity*POW(r                      , 3) -
	   TestGravitySphereExteriorDensity*POW(TestGravitySphereRadius, 3));
    }
 
    /* 1) r^-2 sphere */
 
    if (TestGravitySphereType == 1) {
      if (r < TestGravitySphereRadius)
	fanalytic = GravitationalConstant*POW(TestGravitySphereRadius, 2)*
	  TestGravitySphereInteriorDensity/r;
      else
	fanalytic = 1.0/3.0*GravitationalConstant/(r*r)*
	  (3.0*TestGravitySphereInteriorDensity*
	                                    POW(TestGravitySphereRadius, 3) +
	   TestGravitySphereExteriorDensity*POW(r                      , 3) -
	   TestGravitySphereExteriorDensity*POW(TestGravitySphereRadius, 3));
    }
 
    /* 2) r^-9/4 sphere */
 
    if (TestGravitySphereType == 2) {
      if (r < TestGravitySphereRadius)
	fanalytic = GravitationalConstant*POW(TestGravitySphereRadius, 2)*
	  TestGravitySphereInteriorDensity/
	  POW(r, float(5.0/4.0))/(float(3.0/4.0));
      else
	fanalytic = 1.0/3.0*GravitationalConstant/(r*r)*
	  (3.0*TestGravitySphereInteriorDensity/(float(3.0/4.0))*
                                 POW(TestGravitySphereRadius, float(11.0/4.0)) +
	   TestGravitySphereExteriorDensity*POW(r                      , 3) -
	   TestGravitySphereExteriorDensity*POW(TestGravitySphereRadius, 3));
    }
 
  }
 
  return fanalytic;
 
}
