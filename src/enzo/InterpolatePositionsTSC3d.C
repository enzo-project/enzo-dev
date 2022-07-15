/***********************************************************************
/
/  INTERPOLATE FIELD TO 3D TSC 'PARTICLES' (EITHER PERIODIC OR 'PILE-UP')
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

/* Use FORTRAN or C++ version? */

#define NO_USE_FORTRAN

/* External defines. */

#ifdef USE_FORTRAN
extern "C" void FORTRAN_NAME(tsc3d_period)(
                                      float *posx, float *posy, float *posz,
                                      int *npositions, float *sumfield,
                                      float *field, float leftedge[],
                                      int *dim1, int *dim2, int *dim3,
                                      float *cellsize);
extern "C" void FORTRAN_NAME(tsc3d_pile)(float *posx, float *posy, float *posz,
                                         int *npositions, float *sumfield,
                                         float *field, float leftedge[],
                                         int *edim1, int *edim2, int *edim3,
                                         int *dim1, int *dim2, int *dim3,
                                         float *cellsize);
#endif /* USE_FORTRAN */

/* Periodic version. */

void InterpolatePositionsPeriodicTSC3D(FLOAT *Position[], int Number,
				       float *SumField,
				       float *Field, FLOAT LeftEdge[],
				       int Dimension[], FLOAT CellSize)
{

#ifdef USE_FORTRAN

/* Call FORTRAN routine to do all the hard work. */

FORTRAN_NAME(tsc3d_period)(Position[0], Position[1], Position[2],
                           &Number, SumField, Field, LeftEdge,
                           &Dimension[0], &Dimension[1], &Dimension[2],
                           &CellSize);

#else /* USE_FORTRAN */

  int dim12, i0, i0p, i0m, j0, j0p, j0m, k0, k0p, k0m, n;
  float dx, dy, dz, wx0, wxm, wxp, wy0, wym, wyp, wz0, wzm, wzp,
        xpos, ypos, zpos;

  for (n = 0; n < Number; n++) {

    /* compute index of central cell */

    xpos = ( (*(Position[0] + n)) - LeftEdge[0] ) / CellSize;
    ypos = ( (*(Position[1] + n)) - LeftEdge[1] ) / CellSize;
    zpos = ( (*(Position[2] + n)) - LeftEdge[2] ) / CellSize;

    i0   = int(xpos);
    j0   = int(ypos);
    k0   = int(zpos);

    /* compute the weights */

    dx   = xpos - float(i0) - 0.5;
    dy   = ypos - float(j0) - 0.5;
    dz   = zpos - float(k0) - 0.5;

    wxm  = 0.5*(0.5 - dx)*(0.5 - dx);
    wxp  = dx + wxm;
    wx0  = 1.0 - wxp - wxm;

    wym  = 0.5*(0.5 - dy)*(0.5 - dy);
    wyp  = dy + wym;
    wy0  = 1.0 - wyp - wym;

    wzm  = 0.5*(0.5 - dz)*(0.5 - dz);
    wzp  = dz + wzm;
    wz0  = 1.0 - wzp - wzm;

    /* wrap central index */

    if (i0 < 0)              i0 += Dimension[0];
    if (j0 < 0)              j0 += Dimension[1];
    if (k0 < 0)              k0 += Dimension[2];

    if (i0 >= Dimension[0])  i0 -= Dimension[0];
    if (j0 >= Dimension[1])  j0 -= Dimension[1];
    if (k0 >= Dimension[2])  k0 -= Dimension[2];

    /* determine offsets */

    i0m  = i0 - 1;
    i0p  = i0 + 1;
    j0m  = j0 - 1;
    j0p  = j0 + 1;
    k0m  = k0 - 1;
    k0p  = k0 + 1;

    /* wrap indexes */

    if (i0m < 0)              i0m += Dimension[0];
    if (j0m < 0)              j0m += Dimension[1];
    if (k0m < 0)              k0m += Dimension[2];
    if (i0p >= Dimension[0])  i0p -= Dimension[0];
    if (j0p >= Dimension[1])  j0p -= Dimension[1];
    if (k0p >= Dimension[2])  k0p -= Dimension[2];

    /* interpolate from Field to SumField */

    dim12 = Dimension[0]*Dimension[1];

    *(SumField + n) +=
      wzm*wym*wxm*(*(Field + k0m*dim12 + j0m*Dimension[0] + i0m)) +
      wzm*wym*wx0*(*(Field + k0m*dim12 + j0m*Dimension[0] + i0 )) +
      wzm*wym*wxp*(*(Field + k0m*dim12 + j0m*Dimension[0] + i0p)) +
      wzm*wy0*wxm*(*(Field + k0m*dim12 + j0 *Dimension[0] + i0m)) +
      wzm*wy0*wx0*(*(Field + k0m*dim12 + j0 *Dimension[0] + i0 )) +
      wzm*wy0*wxp*(*(Field + k0m*dim12 + j0 *Dimension[0] + i0p)) +
      wzm*wyp*wxm*(*(Field + k0m*dim12 + j0p*Dimension[0] + i0m)) +
      wzm*wyp*wx0*(*(Field + k0m*dim12 + j0p*Dimension[0] + i0 )) +
      wzm*wyp*wxp*(*(Field + k0m*dim12 + j0p*Dimension[0] + i0p));

    *(SumField + n) +=
      wz0*wym*wxm*(*(Field + k0 *dim12 + j0m*Dimension[0] + i0m)) +
      wz0*wym*wx0*(*(Field + k0 *dim12 + j0m*Dimension[0] + i0 )) +
      wz0*wym*wxp*(*(Field + k0 *dim12 + j0m*Dimension[0] + i0p)) +
      wz0*wy0*wxm*(*(Field + k0 *dim12 + j0 *Dimension[0] + i0m)) +
      wz0*wy0*wx0*(*(Field + k0 *dim12 + j0 *Dimension[0] + i0 )) +
      wz0*wy0*wxp*(*(Field + k0 *dim12 + j0 *Dimension[0] + i0p)) +
      wz0*wyp*wxm*(*(Field + k0 *dim12 + j0p*Dimension[0] + i0m)) +
      wz0*wyp*wx0*(*(Field + k0 *dim12 + j0p*Dimension[0] + i0 )) +
      wz0*wyp*wxp*(*(Field + k0 *dim12 + j0p*Dimension[0] + i0p));

    *(SumField + n) +=
      wzp*wym*wxm*(*(Field + k0p*dim12 + j0m*Dimension[0] + i0m)) +
      wzp*wym*wx0*(*(Field + k0p*dim12 + j0m*Dimension[0] + i0 )) +
      wzp*wym*wxp*(*(Field + k0p*dim12 + j0m*Dimension[0] + i0p)) +
      wzp*wy0*wxm*(*(Field + k0p*dim12 + j0 *Dimension[0] + i0m)) +
      wzp*wy0*wx0*(*(Field + k0p*dim12 + j0 *Dimension[0] + i0 )) +
      wzp*wy0*wxp*(*(Field + k0p*dim12 + j0 *Dimension[0] + i0p)) +
      wzp*wyp*wxm*(*(Field + k0p*dim12 + j0p*Dimension[0] + i0m)) +
      wzp*wyp*wx0*(*(Field + k0p*dim12 + j0p*Dimension[0] + i0 )) +
      wzp*wyp*wxp*(*(Field + k0p*dim12 + j0p*Dimension[0] + i0p));

  } // next particle

#endif /* USE_FORTRAN */

}


/* Particles which extend past the edge of the grid are deposit at the
   edge in this version, causing mass to pile-up. */

void InterpolatePositionsPileUpTSC3D(FLOAT *Position[], int Number,
				     float *SumField,
				     float *Field, FLOAT LeftEdge[],
				     int EffectiveDim[], int Dimension[],
				     FLOAT CellSize)
{

#ifdef USE_FORTRAN

/* Call FORTRAN routine to do all the hard work. */

FORTRAN_NAME(tsc3d_pile)(Position[0], Position[1], Position[2],
                         &Number, SumField, Field, LeftEdge,
                         &EffectiveDim[0], &EffectiveDim[1], &EffectiveDim[2],
                         &Dimension[0], &Dimension[1], &Dimension[2],
                         &CellSize);

#else /* USE_FORTRAN */

  int dim12, i0, i0p, i0m, j0, j0p, j0m, k0, k0p, k0m, n;
  float dx, dy, dz, wx0, wxm, wxp, wy0, wym, wyp, wz0, wzm, wzp,
        xpos, ypos, zpos;

  for (n = 0; n < Number; n++) {

    /* compute index of central cell */

    xpos = ( (*(Position[0] + n)) - LeftEdge[0] ) / CellSize;
    ypos = ( (*(Position[1] + n)) - LeftEdge[1] ) / CellSize;
    zpos = ( (*(Position[2] + n)) - LeftEdge[2] ) / CellSize;

    i0   = int(xpos);
    j0   = int(ypos);
    k0   = int(zpos);

    /* compute the weights */

    dx   = xpos - float(i0) - 0.5;
    dy   = ypos - float(j0) - 0.5;
    dz   = zpos - float(k0) - 0.5;

    wxm  = 0.5*(0.5 - dx)*(0.5 - dx);
    wxp  = dx + wxm;
    wx0  = 1.0 - wxp - wxm;

    wym  = 0.5*(0.5 - dy)*(0.5 - dy);
    wyp  = dy + wym;
    wy0  = 1.0 - wyp - wym;

    wzm  = 0.5*(0.5 - dz)*(0.5 - dz);
    wzp  = dz + wzm;
    wz0  = 1.0 - wzp - wzm;

    /* fix off-edge central index */

    if (i0 < 0)                i0 = 0;
    if (j0 < 0)                j0 = 0;
    if (k0 < 0)                k0 = 0;
    if (i0 >= EffectiveDim[0]) i0 = EffectiveDim[0]-1;
    if (j0 >= EffectiveDim[1]) j0 = EffectiveDim[1]-1;
    if (k0 >= EffectiveDim[2]) k0 = EffectiveDim[2]-1;

    /* determine offsets */

    i0m  = i0 - 1;
    i0p  = i0 + 1;
    j0m  = j0 - 1;
    j0p  = j0 + 1;
    k0m  = k0 - 1;
    k0p  = k0 + 1;

    /* fix off-edge indexes */

    if (i0m < 0)                i0m = 0;
    if (j0m < 0)                j0m = 0;
    if (k0m < 0)                k0m = 0;
    if (i0p >= EffectiveDim[0]) i0p = EffectiveDim[0]-1;
    if (j0p >= EffectiveDim[1]) j0p = EffectiveDim[1]-1;
    if (k0p >= EffectiveDim[2]) k0p = EffectiveDim[2]-1;

    dim12 = Dimension[0]*Dimension[1];

    /* Interpolate from field.  This is split to improve performance. */

    *(SumField + n) +=
      wzm*wym*wxm*(*(Field + k0m*dim12 + j0m*Dimension[0] + i0m)) +
      wzm*wym*wx0*(*(Field + k0m*dim12 + j0m*Dimension[0] + i0 )) +
      wzm*wym*wxp*(*(Field + k0m*dim12 + j0m*Dimension[0] + i0p)) +
      wzm*wy0*wxm*(*(Field + k0m*dim12 + j0 *Dimension[0] + i0m)) +
      wzm*wy0*wx0*(*(Field + k0m*dim12 + j0 *Dimension[0] + i0 )) +
      wzm*wy0*wxp*(*(Field + k0m*dim12 + j0 *Dimension[0] + i0p)) +
      wzm*wyp*wxm*(*(Field + k0m*dim12 + j0p*Dimension[0] + i0m)) +
      wzm*wyp*wx0*(*(Field + k0m*dim12 + j0p*Dimension[0] + i0 )) +
      wzm*wyp*wxp*(*(Field + k0m*dim12 + j0p*Dimension[0] + i0p));

    *(SumField + n) +=
      wz0*wym*wxm*(*(Field + k0 *dim12 + j0m*Dimension[0] + i0m)) +
      wz0*wym*wx0*(*(Field + k0 *dim12 + j0m*Dimension[0] + i0 )) +
      wz0*wym*wxp*(*(Field + k0 *dim12 + j0m*Dimension[0] + i0p)) +
      wz0*wy0*wxm*(*(Field + k0 *dim12 + j0 *Dimension[0] + i0m)) +
      wz0*wy0*wx0*(*(Field + k0 *dim12 + j0 *Dimension[0] + i0 )) +
      wz0*wy0*wxp*(*(Field + k0 *dim12 + j0 *Dimension[0] + i0p)) +
      wz0*wyp*wxm*(*(Field + k0 *dim12 + j0p*Dimension[0] + i0m)) +
      wz0*wyp*wx0*(*(Field + k0 *dim12 + j0p*Dimension[0] + i0 )) +
      wz0*wyp*wxp*(*(Field + k0 *dim12 + j0p*Dimension[0] + i0p));

    *(SumField + n) +=
      wzp*wym*wxm*(*(Field + k0p*dim12 + j0m*Dimension[0] + i0m)) +
      wzp*wym*wx0*(*(Field + k0p*dim12 + j0m*Dimension[0] + i0 )) +
      wzp*wym*wxp*(*(Field + k0p*dim12 + j0m*Dimension[0] + i0p)) +
      wzp*wy0*wxm*(*(Field + k0p*dim12 + j0 *Dimension[0] + i0m)) +
      wzp*wy0*wx0*(*(Field + k0p*dim12 + j0 *Dimension[0] + i0 )) +
      wzp*wy0*wxp*(*(Field + k0p*dim12 + j0 *Dimension[0] + i0p)) +
      wzp*wyp*wxm*(*(Field + k0p*dim12 + j0p*Dimension[0] + i0m)) +
      wzp*wyp*wx0*(*(Field + k0p*dim12 + j0p*Dimension[0] + i0 )) +
      wzp*wyp*wxp*(*(Field + k0p*dim12 + j0p*Dimension[0] + i0p));

  } // next particle

#endif /* USE_FORTRAN */

}
