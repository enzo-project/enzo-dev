/***********************************************************************
/
/  EXTERNAL BOUNDARY CLASS (SETS INFLOW BOUNDARY CONDITIONS FOR WENGEN TEST)
/
/  written by: Tom Abel
/  date:       October 2010
/  modified1:
/
/  PURPOSE:
/
/  BEWARE: Unfortuantely has a number of hard coded values since these routines
/          do not have access to all problem generating ones. Your mileage will vary.
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "MHD2DTestGlobalData.h"
#include "hydro_rk/EOS.h"
 
/* function prototypes */
 
int FindField(int f, int farray[], int n);
 
 
int ExternalBoundary::SetWengenCollidingFlowBoundary(FLOAT time, FLOAT CellLeftEdge[],
					    FLOAT CellWidth[])
{
  /* declarations */
 
  int i, j, dim, index;
  int NumberOfZones[MAX_DIMENSION], Offset[MAX_DIMENSION];
 
  /* Compute size of entire mesh. */
 
  int size = 1;
  for (dim = 0; dim < BoundaryRank; dim++)
    size = size*BoundaryDimension[dim];
 
  /* Find fields: density, total energy, velocity1-3. */
 
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum, B1Num, B2Num, B3Num, PhiNum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				       Vel3Num, TENum, B1Num, B2Num, B3Num, PhiNum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }
 
  /* set the appropriate BoundaryValues on the left side */
  //  fprintf(stdout, "boundary ints: %i %i %i \n", BoundaryDimension[0],BoundaryDimension[1],BoundaryRank);
  for (dim = 0; dim < BoundaryRank; dim++)
    if (BoundaryDimension[dim] != 1) {
 
      /* If the BoundaryValue fields are missing, create them. */
 
      for (int field = 0; field < NumberOfBaryonFields; field++) {
	if (BoundaryValue[field][dim][0] == NULL)
	  BoundaryValue[field][dim][0] =
	    new float[size/BoundaryDimension[dim]];
	if (BoundaryValue[field][dim][1] == NULL)
	  BoundaryValue[field][dim][1] =
	    new float[size/BoundaryDimension[dim]];
      }
      /* Compute quantities needed for boundary face loop (below). */
 
      int dim1, dim2;
      FLOAT x,y;
      float vx,vy;
      dim1 = (dim == 0) ? 1 : 0;
      dim2 = dim1 + 1;
      dim2 = (dim2 == dim) ? dim2+1 : dim2;
      for (i = 0; i < 3; i++) {
	NumberOfZones[i] = max(BoundaryDimension[i] - 2*NumberOfGhostZones,1);
	Offset[i]        = min(NumberOfGhostZones, BoundaryDimension[i]) - 1;
      }
      //      fprintf(stdout, "ints: %i %i %i \n", dim, dim1, dim2);
      //      fprintf(stdout, "ints: %i %i %i \n", NumberOfZones[dim], NumberOfZones[dim1], NumberOfZones[dim2]);
      /* Loop over the boundary face. */
 
      for (i = 0; i < BoundaryDimension[dim1]; i++)
	for (j = 0; j < BoundaryDimension[dim2]; j++) {
 	  x = (float(i-Offset[dim1]) -0.5 )*
	    ((DomainRightEdge[dim1]-DomainLeftEdge[dim1]) /
	    float(NumberOfZones[dim1])) ;
	  y = (float(j-Offset[dim2]) -0.5)*
	    ((DomainRightEdge[dim]-DomainLeftEdge[dim]) /
	     float(NumberOfZones[dim])) ;

	  /* Set the field values. */
	  //	  fprintf(stdout, "lower: %g %g \n", x ,y);

	  index = j*BoundaryDimension[dim1] + i;
	  float pres, eintl, eintu, h, cs, dpdrho, dpde,ramp,rhot;
	  float rho, vx, vy, bx,by, f, vxu, vxl, vyu, vyl, Bxl, Byl, Bxu, Byu, etotl;
	  rho = LowerDensity;
	  vxl = LowerVelocityX;
	  vxu = UpperVelocityX;
	  vyl = LowerVelocityY;
	  vyu = UpperVelocityY;
	  Bxl = LowerBx;
	  Byl = LowerBy;
	  Bxu = UpperBx;
	  Byu = UpperBy;

	  // lower y boundary
	  pres = 0.112611*0.112611* rho; // isothermal sound speed = 0.112611
	  pres = (EOSType > 0) ? EOSSoundSpeed*EOSSoundSpeed*rho : // 
	    EOSSoundSpeed*EOSSoundSpeed*rho ; 
	  EOS(pres, rho, eintl, h, cs, dpdrho, dpde, 0, 1); // compute eintl
	  ramp =  1./(1.+exp(-2/RampWidth*(y-0.5)));
	  f = cos(2.*M_PI*x*10.)*exp(-fabs(y-0.5)*10.)*cos(2.*M_PI*x*3);
	  vx = f * (vxl + ramp*(vxu-vxl)) ;
	  vy = vyl + ramp*(vyu - vyl);
	  bx = (Bxl+ ramp*(Bxu-Bxl))  ;
	  by = (Byl+ ramp*(Byu-Byl))  ;

	  etotl = eintl + 0.5*(vx*vx + vy*vy) + 0.5*(bx*bx+by*by)/rho;

	  *(BoundaryValue[DensNum][dim][0] + index) = rho;
	  *(BoundaryValue[TENum][dim][0]   + index) = etotl;
	  *(BoundaryValue[Vel1Num][dim][0] + index) = vx ;
	  if (BoundaryRank > 1)
	    *(BoundaryValue[Vel2Num][dim][0] + index) = vy;
	  if (BoundaryRank > 2)
	    *(BoundaryValue[Vel3Num][dim][0] + index) = 0.;

	  if (HydroMethod == MHD_RK) {
	    *(BoundaryValue[B1Num][dim][0] + index) = bx ;
	    *(BoundaryValue[B2Num][dim][0] + index) = by;
	    *(BoundaryValue[B3Num][dim][0] + index) = 0;
	  }

	  // upper y boundary

	  y = DomainRightEdge[dim] + (float(j-Offset[dim2]) + 0.5)*
	    ((DomainRightEdge[dim]-DomainLeftEdge[dim]) /
	     float(NumberOfZones[dim])) ;

	  //	  fprintf(stdout, "upper: %g %g \n", x,y);

	  pres = 0.112611*0.112611* rho; // isothermal sound speed = 0.112611
	  pres = (EOSType > 0) ? EOSSoundSpeed*EOSSoundSpeed*rho : // 
	    EOSSoundSpeed*EOSSoundSpeed*rho ; 
	  EOS(pres, rho, eintl, h, cs, dpdrho, dpde, 0, 1); // compute eintl
	  ramp =  1./(1.+exp(-2/RampWidth*(y-0.5)));
	  f = cos(2.*M_PI*x*10.)*exp(-fabs(y-0.5)*10.)*cos(2.*M_PI*x*3);
	  vx = f * (vxl + ramp*(vxu-vxl)) ;
	  vy = vyl + ramp*(vyu - vyl);
	  bx = (Bxl+ ramp*(Bxu-Bxl))  ;
	  by = (Byl+ ramp*(Byu-Byl))  ;

	  etotl = eintl + 0.5*(vx*vx + vy*vy) + 0.5*(bx*bx+by*by)/rho;
 
	  *(BoundaryValue[DensNum][dim][1] + index) = rho;
	  *(BoundaryValue[TENum][dim][1]   + index) = etotl;
	  *(BoundaryValue[Vel1Num][dim][1] + index) = vx ;
	  if (BoundaryRank > 1)
	    *(BoundaryValue[Vel2Num][dim][1] + index) = vy;
	  if (BoundaryRank > 2)
	    *(BoundaryValue[Vel3Num][dim][1] + index) = 0.;

	  if (HydroMethod == MHD_RK) {
	    *(BoundaryValue[B1Num][dim][1] + index) = bx;
	    *(BoundaryValue[B2Num][dim][1] + index) = by;
	    *(BoundaryValue[B3Num][dim][1] + index) = 0;
	  }

 
    } // end loop over boundary directions

    } 
  return SUCCESS;
 
}
