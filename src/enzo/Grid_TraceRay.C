#define DEBUG 1
/***********************************************************************
/
/  GRID CLASS (WALK PHOTON PACKAGES ACROSS GRID)
/
/  written by: Tom Abel
/  date:       August, 2003
/  modified1:
/
/  PURPOSE: Trace a ray across the grid for NumberOfSegments ray segments
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
#include "ExternalBoundary.h"
#include "Fluxes.h"
#include "GridList.h"
#include "Grid.h"

int grid::TraceRay(int NumberOfSegments,
		   FLOAT r,
		   FLOAT x, FLOAT y, FLOAT z,
		   FLOAT ux, FLOAT uy, FLOAT uz,
		   FLOAT dr[],
		   long cindex[],
		   int ci[],
		   int cj[],
		   int ck[]) {
  int i,j,k;
  FLOAT oldr, drx, dry, drz;
  i=j=k=0;

  // source position
  FLOAT sx = (x - r*ux);
  FLOAT sy = (y - r*uy);
  FLOAT sz = (z - r*uz);

  // Current Cell ? 
  i = int((x-GridLeftEdge[0])/CellWidth[0][0]);
  j = int((y-GridLeftEdge[1])/CellWidth[1][0]);
  k = int((z-GridLeftEdge[2])/CellWidth[2][0]);

  // in floating point:
  FLOAT fx = GridLeftEdge[0] + (FLOAT) (i) * CellWidth[0][0];
  FLOAT fy = GridLeftEdge[1] + (FLOAT) (j) * CellWidth[1][0];
  FLOAT fz = GridLeftEdge[2] + (FLOAT) (k) * CellWidth[2][0];

  // in HEALPIX only that z & y component can have a zero directional vector
  if (uz == 0) uz = 1.e-20;
  if (uy == 0) uy = 1.e-20;

  // on cell boundaries the index will change in negative directions
  if (x == fx) i += (int) (sign(ux)-1)/2;
  if (y == fy) j += (int) (sign(uy)-1)/2;
  if (z == fz) k += (int) (sign(uz)-1)/2; 
  
  oldr = r;
  long count = 0;
  while (count < NumberOfSegments) {
    drx = dry = drz = 1;
    // now the radius from the current left edge and the next one 
    drx = (GridLeftEdge[0] + (FLOAT (i) + (sign(ux)+1)/2)*CellWidth[0][0] - sx)/ux;
    dry = (GridLeftEdge[1] + (FLOAT (j) + (sign(uy)+1)/2)*CellWidth[1][0] - sy)/uy;
    drz = (GridLeftEdge[2] + (FLOAT (k) + (sign(uz)+1)/2)*CellWidth[2][0] - sz)/uz;
    
 
    if (drx <= min(dry, drz)) {
      r = drx;
      i += sign(ux);
    }
    if (dry <= min(drx, drz)) {
      r = dry;
      j += sign(uy);
    }
    if (drz <= min(drx, dry)) {
      r = drz;
      k += sign(uz);
    }
    dr[count] = r-oldr;

    ci[count] = i;
    cj[count] = j;
    ck[count] = k;

    cindex[count] = (((k+GridStartIndex[2])*GridDimension[1]+
		      (j+GridStartIndex[1]))*GridDimension[0]+GridStartIndex[0]+i);
      //    cindex[count] = GRIDINDEX(i, j, k);

   if (dr[count] < 0 && DEBUG) {
      fprintf(stderr,"%"ISYM" tr: %"ISYM" %"ISYM" %"ISYM" | %"ISYM" - %"ISYM" %"ISYM" %"ISYM" | r %"GSYM" dr %"GSYM"  (%"GSYM" %"GSYM" %"GSYM") | uxyz: %"GSYM" %"GSYM" %"GSYM" \n", 
  	      count,    i,j,k,	cindex[count], ci[count], cj[count], ck[count], r, dr, drx,dry,drz,x,y,z); 

    }

    oldr = r;
    count++;

  } // end while (count < NumberOfSegments ...    

  
  return SUCCESS;
}
