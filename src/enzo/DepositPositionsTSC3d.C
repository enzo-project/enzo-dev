/***********************************************************************
/
/  DEPOSIT 3D TSC 'PARTICLES' (EITHER PERIODIC OR 'PILE-UP')
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
extern "C" void FORTRAN_NAME(dep_pile_tsc3d)
      (float *posx, float *posy, float *posz, float *mass, int *npositions,
       float *field, float leftedge[],
       int *estart1, int *estart2, int *estart3,
       int *edim1, int *edim2, int *edim3,
       int *dim1, int *dim2, int *dim3, float *cellsize);
#endif /* USE_FORTRAN */

/* Periodic version. */

void DepositPositionsPeriodicTSC3D(float *Position[], float *Mass, 
				   int Number, float *Field, float LeftEdge[],
				   int Dimension[], float CellSize)
{

  int dim12, i0, i0p, i0m, j0, j0p, j0m, k0, k0p, k0m, n;
  float dx, dy, dz, wx0, wxm, wxp, wy0, wym, wyp, wz0, wzm, wzp, 
        xpos, ypos, zpos;

  for (n = 0; n < Number; n++) {

    /* compute index of central cell */

    xpos = ( Position[0][n] - LeftEdge[0] ) / CellSize;
    ypos = ( Position[1][n] - LeftEdge[1] ) / CellSize;
    zpos = ( Position[2][n] - LeftEdge[2] ) / CellSize;

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

    /* check for off-edge particles. */

    if (i0 < 0            ) i0 += Dimension[0];
    if (j0 < 0            ) j0 += Dimension[1];
    if (k0 < 0            ) k0 += Dimension[2];
    if (i0 >= Dimension[0]) i0 -= Dimension[0];
    if (j0 >= Dimension[1]) j0 -= Dimension[1];
    if (k0 >= Dimension[2]) k0 -= Dimension[2];

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

    /* deposit mass */

    dim12 = Dimension[0]*Dimension[1];

    Field[k0m*dim12 + j0m*Dimension[0] + i0m] += wzm*wym*wxm*Mass[n];
    Field[k0m*dim12 + j0m*Dimension[0] + i0 ] += wzm*wym*wx0*Mass[n];
    Field[k0m*dim12 + j0m*Dimension[0] + i0p] += wzm*wym*wxp*Mass[n];
    Field[k0m*dim12 + j0 *Dimension[0] + i0m] += wzm*wy0*wxm*Mass[n];
    Field[k0m*dim12 + j0 *Dimension[0] + i0 ] += wzm*wy0*wx0*Mass[n];
    Field[k0m*dim12 + j0 *Dimension[0] + i0p] += wzm*wy0*wxp*Mass[n];
    Field[k0m*dim12 + j0p*Dimension[0] + i0m] += wzm*wyp*wxm*Mass[n];
    Field[k0m*dim12 + j0p*Dimension[0] + i0 ] += wzm*wyp*wx0*Mass[n];
    Field[k0m*dim12 + j0p*Dimension[0] + i0p] += wzm*wyp*wxp*Mass[n];

    Field[k0 *dim12 + j0m*Dimension[0] + i0m] += wz0*wym*wxm*Mass[n];
    Field[k0 *dim12 + j0m*Dimension[0] + i0 ] += wz0*wym*wx0*Mass[n];
    Field[k0 *dim12 + j0m*Dimension[0] + i0p] += wz0*wym*wxp*Mass[n];
    Field[k0 *dim12 + j0 *Dimension[0] + i0m] += wz0*wy0*wxm*Mass[n];
    Field[k0 *dim12 + j0 *Dimension[0] + i0 ] += wz0*wy0*wx0*Mass[n];
    Field[k0 *dim12 + j0 *Dimension[0] + i0p] += wz0*wy0*wxp*Mass[n];
    Field[k0 *dim12 + j0p*Dimension[0] + i0m] += wz0*wyp*wxm*Mass[n];
    Field[k0 *dim12 + j0p*Dimension[0] + i0 ] += wz0*wyp*wx0*Mass[n];
    Field[k0 *dim12 + j0p*Dimension[0] + i0p] += wz0*wyp*wxp*Mass[n];

    Field[k0p*dim12 + j0m*Dimension[0] + i0m] += wzp*wym*wxm*Mass[n];
    Field[k0p*dim12 + j0m*Dimension[0] + i0 ] += wzp*wym*wx0*Mass[n];
    Field[k0p*dim12 + j0m*Dimension[0] + i0p] += wzp*wym*wxp*Mass[n];
    Field[k0p*dim12 + j0 *Dimension[0] + i0m] += wzp*wy0*wxm*Mass[n];
    Field[k0p*dim12 + j0 *Dimension[0] + i0 ] += wzp*wy0*wx0*Mass[n];
    Field[k0p*dim12 + j0 *Dimension[0] + i0p] += wzp*wy0*wxp*Mass[n];
    Field[k0p*dim12 + j0p*Dimension[0] + i0m] += wzp*wyp*wxm*Mass[n];
    Field[k0p*dim12 + j0p*Dimension[0] + i0 ] += wzp*wyp*wx0*Mass[n];
    Field[k0p*dim12 + j0p*Dimension[0] + i0p] += wzp*wyp*wxp*Mass[n];

  } // next particle

}

/* Particles which extend past the edge of the grid are deposit at the
   edge in this version, causing mass to pile-up. */

void DepositPositionsPileUpTSC3D(FLOAT *Position[], float *Mass, 
				 int Number, float *Field, FLOAT LeftEdge[],
				 int EffectiveStart[], int EffectiveDim[], 
				 int Dimension[], float CellSize)
{

  //  fprintf(stderr, "DepositTSC: %d %d %d\n", Dimension[0], Dimension[1], Dimension[2]);

#ifdef USE_FORTRAN

/* Call FORTRAN routine to do all the hard work. */

FORTRAN_NAME(dep_pile_tsc3d)(
                 Position[0], Position[1], Position[2], Mass, &Number, 
                 Field, LeftEdge, 
                 &EffectiveStart[0], &EffectiveStart[1], &EffectiveStart[2],
                 &EffectiveDim[0], &EffectiveDim[1], &EffectiveDim[2],
                 &Dimension[0], &Dimension[1], &Dimension[2], &CellSize);

#else /* USE_FORTRAN */

  int dim12, i0, i0p, i0m, j0, j0p, j0m, k0, k0p, k0m, n;
  float dx, dy, dz, wx0, wxm, wxp, wy0, wym, wyp, wz0, wzm, wzp, 
        xpos, ypos, zpos;

  for (n = 0; n < Number; n++) {

    /* compute index of central cell */

    xpos = ( Position[0][n] - LeftEdge[0] ) / CellSize;
    ypos = ( Position[1][n] - LeftEdge[1] ) / CellSize;
    zpos = ( Position[2][n] - LeftEdge[2] ) / CellSize;

    i0   = int(xpos);
    j0   = int(ypos);
    k0   = int(zpos);
    
    /* check for off-edge particles. */

    if (i0 < EffectiveStart[0] || i0 >= EffectiveDim[0] ||
        j0 < EffectiveStart[1] || j0 >= EffectiveDim[1] ||
        k0 < EffectiveStart[2] || k0 >= EffectiveDim[2]) {
    
      fprintf(stdout,"Reject particle %"ISYM" in DepositPositionsTSC3d.C: off-edge.\n", n);
      continue;
    }

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

    /* determine offsets */

    i0m  = i0 - 1;
    i0p  = i0 + 1;
    j0m  = j0 - 1;
    j0p  = j0 + 1;
    k0m  = k0 - 1;
    k0p  = k0 + 1;

    /* wrap indexes */

    if (i0m <  EffectiveStart[0]) i0m = EffectiveStart[0];
    if (j0m <  EffectiveStart[1]) j0m = EffectiveStart[1];
    if (k0m <  EffectiveStart[2]) k0m = EffectiveStart[2];
    if (i0p >= EffectiveDim[0]  ) i0p = EffectiveDim[0]-1;
    if (j0p >= EffectiveDim[1]  ) j0p = EffectiveDim[1]-1;
    if (k0p >= EffectiveDim[2]  ) k0p = EffectiveDim[2]-1;

    /* deposit mass */

    dim12 = Dimension[0]*Dimension[1];

    Field[k0m*dim12 + j0m*Dimension[0] + i0m] += wzm*wym*wxm*Mass[n];
    Field[k0m*dim12 + j0m*Dimension[0] + i0 ] += wzm*wym*wx0*Mass[n];
    Field[k0m*dim12 + j0m*Dimension[0] + i0p] += wzm*wym*wxp*Mass[n];
    Field[k0m*dim12 + j0 *Dimension[0] + i0m] += wzm*wy0*wxm*Mass[n];
    Field[k0m*dim12 + j0 *Dimension[0] + i0 ] += wzm*wy0*wx0*Mass[n];
    Field[k0m*dim12 + j0 *Dimension[0] + i0p] += wzm*wy0*wxp*Mass[n];
    Field[k0m*dim12 + j0p*Dimension[0] + i0m] += wzm*wyp*wxm*Mass[n];
    Field[k0m*dim12 + j0p*Dimension[0] + i0 ] += wzm*wyp*wx0*Mass[n];
    Field[k0m*dim12 + j0p*Dimension[0] + i0p] += wzm*wyp*wxp*Mass[n];

    Field[k0 *dim12 + j0m*Dimension[0] + i0m] += wz0*wym*wxm*Mass[n];
    Field[k0 *dim12 + j0m*Dimension[0] + i0 ] += wz0*wym*wx0*Mass[n];
    Field[k0 *dim12 + j0m*Dimension[0] + i0p] += wz0*wym*wxp*Mass[n];
    Field[k0 *dim12 + j0 *Dimension[0] + i0m] += wz0*wy0*wxm*Mass[n];
    Field[k0 *dim12 + j0 *Dimension[0] + i0 ] += wz0*wy0*wx0*Mass[n];
    Field[k0 *dim12 + j0 *Dimension[0] + i0p] += wz0*wy0*wxp*Mass[n];
    Field[k0 *dim12 + j0p*Dimension[0] + i0m] += wz0*wyp*wxm*Mass[n];
    Field[k0 *dim12 + j0p*Dimension[0] + i0 ] += wz0*wyp*wx0*Mass[n];
    Field[k0 *dim12 + j0p*Dimension[0] + i0p] += wz0*wyp*wxp*Mass[n];

    Field[k0p*dim12 + j0m*Dimension[0] + i0m] += wzp*wym*wxm*Mass[n];
    Field[k0p*dim12 + j0m*Dimension[0] + i0 ] += wzp*wym*wx0*Mass[n];
    Field[k0p*dim12 + j0m*Dimension[0] + i0p] += wzp*wym*wxp*Mass[n];
    Field[k0p*dim12 + j0 *Dimension[0] + i0m] += wzp*wy0*wxm*Mass[n];
    Field[k0p*dim12 + j0 *Dimension[0] + i0 ] += wzp*wy0*wx0*Mass[n];
    Field[k0p*dim12 + j0 *Dimension[0] + i0p] += wzp*wy0*wxp*Mass[n];
    Field[k0p*dim12 + j0p*Dimension[0] + i0m] += wzp*wyp*wxm*Mass[n];
    Field[k0p*dim12 + j0p*Dimension[0] + i0 ] += wzp*wyp*wx0*Mass[n];
    Field[k0p*dim12 + j0p*Dimension[0] + i0p] += wzp*wyp*wxp*Mass[n];

  } // next particle

#endif /* USE_FORTRAN */

}
