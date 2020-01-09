/***********************************************************************
/
/  DEPOSIT 1D TSC 'PARTICLES' (EITHER PERIODIC OR 'PILE-UP')
/
/  written by: Greg Bryan
/  date:       March, 1995
/  modified1:
/
/  PURPOSE:
/
/  NOTE: 
/
************************************************************************/

#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"

/* Periodic version. */

void DepositPositionsPeriodicTSC1D(FLOAT *Position[], float *Mass, 
				   int Number, float *Field, FLOAT LeftEdge[],
				   int Dimension[], float CellSize)
{
  int i0, i0p, i0m, n;
  float dx, wx0, wxm, wxp, xpos;

  for (n = 0; n < Number; n++) {

    /* compute index of central cell */

    xpos = ( (*(Position[0] + n)) - LeftEdge[0] ) / CellSize;
    i0   = int(xpos);

    /* compute the weights */

    dx   = xpos - float(i0) - 0.5;

    wxm  = 0.5*(0.5 - dx)*(0.5 - dx);
    wxp  = dx + wxm;
    wx0  = 1.0 - wxp - wxm;

    /* check for off-edge particles. */

    if (i0 < 0            ) i0 += Dimension[0];
    if (i0 >= Dimension[0]) i0 -= Dimension[0];

    /* determine offsets */

    i0m  = i0 - 1;
    i0p  = i0 + 1;

    /* wrap indexes */

    if (i0m < 0)              i0m += Dimension[0];
    if (i0p >= Dimension[0])  i0p -= Dimension[0];

    /* deposit mass */

    *(Field + i0m) += wxm*(*(Mass + n));
    *(Field + i0 ) += wx0*(*(Mass + n));
    *(Field + i0p) += wxp*(*(Mass + n));

  } // next particle

}

/* Particles which extend past the edge of the grid are deposit at the
   edge in this version, causing mass to pile-up. */

void DepositPositionsPileUpTSC1D(FLOAT *Position[], float *Mass, 
				 int Number, float *Field, FLOAT LeftEdge[],
				 int EffectiveStart[], int EffectiveDim[], 
				 float CellSize)
{
  int i0, i0p, i0m, n;
  float dx, wx0, wxm, wxp, xpos;

  for (n = 0; n < Number; n++) {

    /* compute index of central cell */

    xpos = ( (*(Position[0] + n)) - LeftEdge[0] ) / CellSize;
    i0   = int(xpos);

    /* compute the weights */

    dx   = xpos - float(i0) - 0.5;

    wxm  = 0.5*(0.5 - dx)*(0.5 - dx);
    wxp  = dx + wxm;
    wx0  = 1.0 - wxp - wxm;

    /* check for off-edge particles. */

    if (i0 < EffectiveStart[0] || i0 >= EffectiveDim[0]) 
      continue;

    /* determine offsets */

    i0m  = i0 - 1;
    i0p  = i0 + 1;

    /* wrap indexes */

    if (i0m <  EffectiveStart[0]) i0m = EffectiveStart[0];
    if (i0p >= EffectiveDim[0]  ) i0p = EffectiveDim[0]-1;

    /* deposit mass */

    *(Field + i0m) += wxm*(*(Mass + n));
    *(Field + i0 ) += wx0*(*(Mass + n));
    *(Field + i0p) += wxp*(*(Mass + n));

  } // next particle

}
