/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR SEDOV BLAST WAVE TEST in 3D) 
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Alexei Kritsuk, March 2005.
/
/  PURPOSE: Sets the total energy in the initial explosion region
/           based on precomputed analytical 1D-spherical solution at 
/           time SedovBlastInitialTime. 
/
/  RETURNS: FAIL or SUCCESS
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
#include "SedovBlastGlobalData.h"

int grid::SedovBlastInitializeGrid3D(char * SedovBlastFileName)
{

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if (GridRank != 3)
    ENZO_FAIL("GridRank != 3!\n");

  /* declarations */

  int size = 1, dim;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];


  /* Open initial data file and read. */

  FILE *fptr;
  if ((fptr = fopen(SedovBlastFileName, "r")) == NULL) {
    ENZO_VFAIL("Cannot open SedovBlast Initial Data File %s\n", 
	   SedovBlastFileName)
  }

  int nl = 0, nlines;
  char line[MAX_LINE_LENGTH];
  if (fgets(line, MAX_LINE_LENGTH, fptr) == NULL)
    ERROR_MESSAGE;
  sscanf(line, "%"ISYM, &nlines);
  if (debug)
    printf("GSBIG: %"ISYM" lines\n", nlines);
  fgets(line, MAX_LINE_LENGTH, fptr); // skip variable names
  float * rad = new float[nlines];
  float * den = new float[nlines];
  float * pre = new float[nlines];
  float * vel = new float[nlines];
  float dummy;
  int idummy;
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
    if (nl > nlines)
      ERROR_MESSAGE;
    sscanf(line, "%"ISYM" %"ISYM" %"ISYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM,
	   &idummy, // nt   - skip
	   &idummy, // l    - skip
	   &idummy, // lev  - skip
	   rad+nl,  // x    = radial coordinate
	   vel+nl,  // xdot = radial velocity
	   den+nl,  // rho  = density
	   &dummy,  // tev  - skip temperature (eV)
	   pre+nl   // p    = pressure
	   );

    /* Normalization corrections for Explosion Energy != 1.
       P ~ E, V ~ sqrt(E)
    */

    pre[nl] *= 2.0259759e-16;
    vel[nl] *= 1.4233680e-8;
    if (debug)
      printf("%"GSYM" %"GSYM" %"GSYM" %"GSYM"\n", 
	     rad[nl], 
	     den[nl], 
	     pre[nl], 
	     vel[nl]);
    nl++;
  }
  fclose(fptr);

  if (nl != nlines)
    ERROR_MESSAGE;

  /* get prepared to "interpolate" */

  float r0 = rad[0];
  float dr = (rad[nlines-1] - r0)/(nlines-1);

  /* Set fields as a function of radius: x^2+y^2+z^2 = r^2. 
     To avoid interpolation here, provide as detailed initial data
     table as needed for the highest resolution subgrid.
   */

  int index, jndex, kndex, i;
  float zonex, zoney, zonez, radius;
  for (i = 0; i < size; i++) {
    index  = i % GridDimension[0];
    jndex  = (i-index) % (GridDimension[0]*GridDimension[1]);
    kndex  = (i-index - jndex)/(GridDimension[0]*GridDimension[1]);
    jndex /= GridDimension[0];
    zonex  = *(CellLeftEdge[0] + index) + 0.5*(*(CellWidth[0] + index));
    zoney  = *(CellLeftEdge[1] + jndex) + 0.5*(*(CellWidth[1] + jndex));
    zonez  = *(CellLeftEdge[2] + kndex) + 0.5*(*(CellWidth[2] + kndex));
    if (SedovBlastFullBox) {
      zonex -= 0.5*CellWidth[0][0]*(GridEndIndex[0] - GridStartIndex[0] + 1) + 
	GridLeftEdge[0];
      zoney -= 0.5*CellWidth[1][0]*(GridEndIndex[1] - GridStartIndex[1] + 1) + 
	GridLeftEdge[1];
      zonez -= 0.5*CellWidth[2][0]*(GridEndIndex[2] - GridStartIndex[2] + 1) + 
	GridLeftEdge[2];
    }
    radius = sqrt(zonex*zonex + zoney*zoney + zonez*zonez);

    int num = nint((radius - r0)/dr);
    if (num < nlines) {

      BaryonField[0][i] = den[num];
      BaryonField[1][i] = pre[num]/den[num]/(Gamma-1.0) + 
	vel[num]*vel[num]/2.0; 
      BaryonField[2][i] = vel[num]*zonex/radius;
      BaryonField[3][i] = vel[num]*zoney/radius; 
      BaryonField[4][i] = vel[num]*zonez/radius;
    }
  }

  /* clean up */

  delete [] rad;
  delete [] den;
  delete [] pre;
  delete [] vel;


  return SUCCESS;
}
