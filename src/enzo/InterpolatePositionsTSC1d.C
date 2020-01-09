/***********************************************************************
/
/  INTERPOLATE FIELD TO 1D TSC 'PARTICLES' (EITHER PERIODIC OR 'PILE-UP')
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
#include "typedefs.h"

/* Periodic version. */

void InterpolatePositionsPeriodicTSC1D(FLOAT *Position[], int Number,
                                       float *SumField,
                                       float *Field, FLOAT LeftEdge[],
				                               int Dimension[], FLOAT CellSize)
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

    /* wrap central index */

    if (i0 < 0)              i0 += Dimension[0];
    if (i0 >= Dimension[0])  i0 -= Dimension[0];

    /* determine offsets */

    i0m  = i0 - 1;
    i0p  = i0 + 1;

    /* wrap indexes */

    if (i0m < 0)              i0m += Dimension[0];
    if (i0p >= Dimension[0])  i0p -= Dimension[0];

    /* interpolate from Field to SumField */

    *(SumField + n) += wxm*(*(Field + i0m)) +
                       wx0*(*(Field + i0 )) +
		       wxp*(*(Field + i0p));

  } // next particle

}


/* Particles which extend past the edge of the grid are deposit at the
   edge in this version, causing mass to pile-up. */

void InterpolatePositionsPileUpTSC1D(FLOAT *Position[], int Number,
				     float *SumField,
				     float *Field, FLOAT LeftEdge[],
				     int EffectiveDim[], FLOAT CellSize)
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

    /* fix off-edge central index */

    if (i0 < 0)                i0 = 0;
    if (i0 >= EffectiveDim[0]) i0 = EffectiveDim[0]-1;

    /* determine offsets */

    i0m  = i0 - 1;
    i0p  = i0 + 1;

    /* fix off-edge indexes */

    if (i0m < 0)                i0m = 0;
    if (i0p >= EffectiveDim[0]) i0p = EffectiveDim[0]-1;

    /* interpolate from Field to SumField */

    *(SumField + n) += wxm*(*(Field + i0m)) +
                       wx0*(*(Field + i0 )) +
		       wxp*(*(Field + i0p));

  } // next particle

}
