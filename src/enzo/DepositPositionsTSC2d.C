/***********************************************************************
/
/  DEPOSIT 2D TSC 'PARTICLES' (EITHER PERIODIC OR 'PILE-UP')
/
/  written by: Greg Bryan
/  date:       May, 1995
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

void DepositPositionsPeriodicTSC2D(float *Position[], float *Mass, 
				   int Number, float *Field, float LeftEdge[],
				   int Dimension[], float CellSize)
{
  int i0, i0p, i0m, j0, j0p, j0m, n;
  float dx, dy, wx0, wxm, wxp, wy0, wym, wyp, xpos, ypos;

  for (n = 0; n < Number; n++) {

    /* compute index of central cell */

    xpos = ( (*(Position[0] + n)) - LeftEdge[0] ) / CellSize;
    ypos = ( (*(Position[1] + n)) - LeftEdge[1] ) / CellSize;
    i0   = int(xpos);
    j0   = int(ypos);

    /* compute the weights */

    dx   = xpos - float(i0) - 0.5;
    dy   = ypos - float(j0) - 0.5;

    wxm  = 0.5*(0.5 - dx)*(0.5 - dx);
    wxp  = dx + wxm;
    wx0  = 1.0 - wxp - wxm;

    wym  = 0.5*(0.5 - dy)*(0.5 - dy);
    wyp  = dy + wym;
    wy0  = 1.0 - wyp - wym;

    /* check for off-edge particles. */

    if (i0 < 0            ) i0 += Dimension[0];
    if (j0 < 0            ) j0 += Dimension[1];
    if (i0 >= Dimension[0]) i0 -= Dimension[0];
    if (j0 >= Dimension[1]) j0 -= Dimension[1];

    /* determine offsets */

    i0m  = i0 - 1;
    i0p  = i0 + 1;
    j0m  = j0 - 1;
    j0p  = j0 + 1;

    /* wrap indexes */

    if (i0m < 0)              i0m += Dimension[0];
    if (j0m < 0)              j0m += Dimension[1];
    if (i0p >= Dimension[0])  i0p -= Dimension[0];
    if (j0p >= Dimension[1])  j0p -= Dimension[1];

    /* deposit mass */

    *(Field + j0m*Dimension[0] + i0m) += wym*wxm*(*(Mass + n));
    *(Field + j0m*Dimension[0] + i0 ) += wym*wx0*(*(Mass + n));
    *(Field + j0m*Dimension[0] + i0p) += wym*wxp*(*(Mass + n));
    *(Field + j0 *Dimension[0] + i0m) += wy0*wxm*(*(Mass + n));
    *(Field + j0 *Dimension[0] + i0 ) += wy0*wx0*(*(Mass + n));
    *(Field + j0 *Dimension[0] + i0p) += wy0*wxp*(*(Mass + n));
    *(Field + j0p*Dimension[0] + i0m) += wyp*wxm*(*(Mass + n));
    *(Field + j0p*Dimension[0] + i0 ) += wyp*wx0*(*(Mass + n));
    *(Field + j0p*Dimension[0] + i0p) += wyp*wxp*(*(Mass + n));

  } // next particle

}

/* Particles which extend past the edge of the grid are deposit at the
   edge in this version, causing mass to pile-up. */

void DepositPositionsPileUpTSC2D(FLOAT *Position[], float *Mass, 
				 int Number, float *Field, FLOAT LeftEdge[],
				 int EffectiveStart[], int EffectiveDim[], 
				 int Dimension[], float CellSize)
{
  int i0, i0p, i0m, j0, j0p, j0m, n;
  float dx, dy, wx0, wxm, wxp, wy0, wym, wyp, xpos, ypos;

  for (n = 0; n < Number; n++) {

    /* compute index of central cell */

    xpos = ( (*(Position[0] + n)) - LeftEdge[0] ) / CellSize;
    ypos = ( (*(Position[1] + n)) - LeftEdge[1] ) / CellSize;
    i0   = int(xpos);
    j0   = int(ypos);

    /* compute the weights */

    dx   = xpos - float(i0) - 0.5;
    dy   = ypos - float(j0) - 0.5;

    wxm  = 0.5*(0.5 - dx)*(0.5 - dx);
    wxp  = dx + wxm;
    wx0  = 1.0 - wxp - wxm;

    wym  = 0.5*(0.5 - dy)*(0.5 - dy);
    wyp  = dy + wym;
    wy0  = 1.0 - wyp - wym;

    /* check for off-edge particles. */

    if (i0 < EffectiveStart[0] || i0 >= EffectiveDim[0] ||
        j0 < EffectiveStart[1] || j0 >= EffectiveDim[1]) 
      continue;

    /* determine offsets */

    i0m  = i0 - 1;
    i0p  = i0 + 1;
    j0m  = j0 - 1;
    j0p  = j0 + 1;

    /* wrap indexes */

    if (i0m <  EffectiveStart[0]) i0m = EffectiveStart[0];
    if (j0m <  EffectiveStart[1]) j0m = EffectiveStart[1];
    if (i0p >= EffectiveDim[0]  ) i0p = EffectiveDim[0]-1;
    if (j0p >= EffectiveDim[1]  ) j0p = EffectiveDim[1]-1;

    /* deposit mass */

    *(Field + j0m*Dimension[0] + i0m) += wym*wxm*(*(Mass + n));
    *(Field + j0m*Dimension[0] + i0 ) += wym*wx0*(*(Mass + n));
    *(Field + j0m*Dimension[0] + i0p) += wym*wxp*(*(Mass + n));
    *(Field + j0 *Dimension[0] + i0m) += wy0*wxm*(*(Mass + n));
    *(Field + j0 *Dimension[0] + i0 ) += wy0*wx0*(*(Mass + n));
    *(Field + j0 *Dimension[0] + i0p) += wy0*wxp*(*(Mass + n));
    *(Field + j0p*Dimension[0] + i0m) += wyp*wxm*(*(Mass + n));
    *(Field + j0p*Dimension[0] + i0 ) += wyp*wx0*(*(Mass + n));
    *(Field + j0p*Dimension[0] + i0p) += wyp*wxp*(*(Mass + n));

  } // next particle

}
