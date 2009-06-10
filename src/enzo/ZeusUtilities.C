/*=======================================================================
  //////////////////////////  ZEUS_UTILITIES  \\\\\\\\\\\\\\\\\\\\\\\\\\\
  ======================================================================= */

#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "fortran.def"

#define IDX(a,b,c) (( (c)*jdim + (b) )*idim + (a))

/* ===================================================================== */

int vanlr1_zc(float *q, int idim, int jdim, int kdim, int is, int ie, 
	      int j, int k, float dx[], float dt, float u[], float qstar[])
{

  /* Definitions */

  int i;
  float dq1[MAX_ANY_SINGLE_DIRECTION], dq2[MAX_ANY_SINGLE_DIRECTION];

  for (i = is-2; i <= ie; i++)
    dq1[i] = (q[IDX(i+1,j,k)] - q[IDX(i,j,k)])*dt/dx[i];

  for (i = is-1; i <= ie; i++) {
    if (dq1[i]*dq1[i-1] > 0.0)
      dq2[i] = 2.0*dq1[i]*dq1[i-1]/(dq1[i-1] + dq1[i]);
    else
      dq2[i] = 0.0;
  }

  dq2[is-2] = dq1[is-2];
  dq2[ie+1] = dq1[ie];
  for (i = is-1; i <= ie+1; i++) {
    if (u[i] > 0.0)
      qstar[i] = q[IDX(i-1,j,k)] + (dx[i-1]/dt - u[i])*0.5*dq2[i-1];
    else
      qstar[i] = q[IDX(i  ,j,k)] - (dx[i  ]/dt + u[i])*0.5*dq2[i  ];
  }

//  qstar[is-1] = 0.5*(q[IDX(is-1,j,k)] + q[IDX(is-2,j,k)]);
//  qstar[ie+1] = 0.5*(q[IDX(ie+1,j,k)] + q[IDX(ie  ,j,k)]);

  return SUCCESS;
}



/* ===================================================================== */

int vanlr1_fc(float *q, int idim, int jdim, int kdim, int is, int ie, 
	      int j, int k, float dx[], float dt, float u[], float qstar[])
{

  /* Definitions */

  int i;
  float dq1[MAX_ANY_SINGLE_DIRECTION], dq2[MAX_ANY_SINGLE_DIRECTION];

  for (i = is-1; i <= ie+1; i++) 
    dq1[i] = (q[IDX(i+1,j,k)] - q[IDX(i,j,k)])*dt/dx[i];

  for (i = is-1; i <= ie; i++) {
    if (dq1[i]*dq1[i+1] > 0.0)
      dq2[i] = 2.0*dq1[i]*dq1[i+1]/(dq1[i+1] + dq1[i]);
    else
      dq2[i] = 0.0;
  }

  dq2[is-2] = dq1[is-1];
  dq2[ie+1] = dq1[ie+1];
  for (i = is-1; i <= ie+1; i++) {
    if (u[i] > 0.0)
      qstar[i] = q[IDX(i  ,j,k)] + (dx[i]/dt - u[i])*0.5*dq2[i-1];
    else
      qstar[i] = q[IDX(i+1,j,k)] - (dx[i+1]/dt + u[i])*0.5*dq2[i];
  }

//  qstar[is-1] = 0.5*(q[IDX(is-1,j,k)] + q[IDX(is  ,j,k)])
//  qstar[ie+1] = 0.5*(q[IDX(ie+1,j,k)] + q[IDX(ie+2,j,k)])

  return SUCCESS;
}


/* ===================================================================== */

int vanlr2_zc(float *q, int idim, int jdim, int kdim, int js, int je, 
	      int i, int k, float dy[], float dt, float v[], float qstar[])
{

  /* Definitions */

  int j;
  float dq1[MAX_ANY_SINGLE_DIRECTION], dq2[MAX_ANY_SINGLE_DIRECTION];

  for (j = js-2; j <= je; j++)
    dq1[j] = (q[IDX(i,j+1,k)] - q[IDX(i,j,k)])*dt/dy[j];

  for (j = js-1; j <= je; j++) {
    if (dq1[j]*dq1[j-1] > 0.0)
      dq2[j] = 2.0*dq1[j]*dq1[j-1]/(dq1[j-1] + dq1[j]);
    else
      dq2[j] = 0.0;
  }

  dq2[js-2] = dq1[js-2];
  dq2[je+1] = dq1[je];
  for (j = js-1; j <= je+1; j++) {
    if (v[j] > 0.0)
      qstar[j] = q[IDX(i,j-1,k)] + (dy[j-1]/dt - v[j])*0.5*dq2[j-1];
    else
      qstar[j] = q[IDX(i,j  ,k)] - (dy[j  ]/dt + v[j])*0.5*dq2[j  ];
  }

//  qstar[js-1] = 0.5*(q[IDX(i,js-1,k)] + q[IDX(i,js-2,k)]);
//  qstar[je+1] = 0.5*(q[IDX(i,je+1,k)] + q[IDX(i,je  ,k)]);

  return SUCCESS;
}



/* ===================================================================== */

int vanlr2_fc(float *q, int idim, int jdim, int kdim, int js, int je, 
	      int i, int k, float dy[], float dt, float v[], float qstar[])
{

  /* Definitions */

  int j;
  float dq1[MAX_ANY_SINGLE_DIRECTION], dq2[MAX_ANY_SINGLE_DIRECTION];

  for (j = js-1; j <= je+1; j++)
    dq1[j] = (q[IDX(i,j+1,k)] - q[IDX(i,j,k)])*dt/dy[j];

  for (j = js-1; j <= je; j++) {
    if (dq1[j]*dq1[j+1] > 0.0)
      dq2[j] = 2.0*dq1[j]*dq1[j+1]/(dq1[j+1] + dq1[j]);
    else
      dq2[j] = 0.0;
  }

  dq2[js-2] = dq1[js-1];
  dq2[je+1] = dq1[je+1];
  for (j = js-1; j <= je+1; j++) {
    if (v[j] > 0.0)
      qstar[j] = q[IDX(i,j  ,k)] + (dy[j  ]/dt - v[j])*0.5*dq2[j-1];
    else
      qstar[j] = q[IDX(i,j+1,k)] - (dy[j+1]/dt + v[j])*0.5*dq2[j  ];
  }
  
//  qstar[js-1] = 0.5*(q[IDX(i,js-1,k)] + q[IDX(i,js  ,k)]);
//  qstar[je+1] = 0.5*(q[IDX(i,js+1,k)] + q[IDX(i,je+2,k)]);

  return SUCCESS;
}


/* ===================================================================== */

int vanlr3_zc(float *q, int idim, int jdim, int kdim, int ks, int ke, 
	      int i, int j, float dz[], float dt, float w[], float qstar[])
{

  /* Definitions */

  int k;
  float dq1[MAX_ANY_SINGLE_DIRECTION], dq2[MAX_ANY_SINGLE_DIRECTION];

  for (k = ks-2; k <= ke; k++)
    dq1[k] = (q[IDX(i,j,k+1)] - q[IDX(i,j,k)])*dt/dz[k];

  for (k = ks-1; k <= ke; k++) {
    if (dq1[k]*dq1[k-1] > 0.0)
      dq2[k] = 2.0*dq1[k]*dq1[k-1]/(dq1[k-1] + dq1[k]);
    else
      dq2[k] = 0.0;
  }

  dq2[ks-2] = dq1[ks-2];
  dq2[ke+1] = dq1[ke];
  for (k = ks-1; k <= ke+1; k++) {
    if (w[k] > 0.0)
      qstar[k] = q[IDX(i,j,k-1)] + (dz[k-1]/dt - w[k])*0.5*dq2[k-1];
    else
      qstar[k] = q[IDX(i,j,k  )] - (dz[k  ]/dt + w[k])*0.5*dq2[k  ];
  }

//  qstar[ks-1] = 0.5*(q[IDX(i,j,ks-1)] + q[IDX(i,j,ks-2)]);
//  qstar[ke+1] = 0.5*(q[IDX(i,j,ke+1)] + q[IDX(i,j,ke  )]);

  return SUCCESS;
}



/* ===================================================================== */

int vanlr3_fc(float *q, int idim, int jdim, int kdim, int ks, int ke, 
	      int i, int j, float dz[], float dt, float w[], float qstar[])
{

  /* Definitions */

  int k;
  float dq1[MAX_ANY_SINGLE_DIRECTION], dq2[MAX_ANY_SINGLE_DIRECTION];

  for (k = ks-1; k <= ke+1; k++)
    dq1[k] = (q[IDX(i,j,k+1)] - q[IDX(i,j,k)])*dt/dz[k];

  for (k = ks-1; k <= ke; k++) {
    if (dq1[k]*dq1[k+1] > 0.0)
      dq2[k] = 2.0*dq1[k]*dq1[k+1]/(dq1[k+1] + dq1[k]);
    else
      dq2[k] = 0.0;
  }

  dq2[ks-2] = dq1[ks-1];
  dq2[ke+1] = dq1[ke+1];
  for (k = ks-1; k <= ke+1; k++) {
    if (w[k] > 0.0)
      qstar[k] = q[IDX(i,j,k  )] + (dz[k  ]/dt - w[k])*0.5*dq2[k-1];
    else
      qstar[k] = q[IDX(i,j,k+1)] - (dz[k+1]/dt + w[k])*0.5*dq2[k  ];
  }

//  qstar[ks-1] = 0.5*(q[IDX(i,j,ks-1)] + q[IDX(i,j,ks  )]);
//  qstar[ke+1] = 0.5*(q[IDX(i,j,ke+1)] + q[IDX(i,j,ke+2)]);

  return SUCCESS;
}

