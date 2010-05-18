/***********************************************************************
/
/  FOF ROUTINE: MERGE PARTICLE USE FRIEND-OF-FRIEND ALGORITHM
/
/  written by: Peng Wang
/  date:       Januaray, 2009
/  modified1:
/
/  PURPOSE:
/
/  DESCRIPTION:
/
************************************************************************/

#include <math.h>
#include "macros_and_parameters.h"

int fofi(FLOAT *x, FLOAT *y, FLOAT *z, const int &np, const FLOAT &l, const int &ip,
         int *group, int &ng)
{

  FLOAT l2 = l*l;

  for (int i = 0; i < np; i++) {
    if (i == ip || group[i] >= 0) continue;
    FLOAT r2 = pow(x[i] - x[ip],2) + pow(y[i] - y[ip], 2) + pow(z[i] - z[ip], 2);
    if (r2 <= l2) {
      group[i] = group[ip];
      fofi(x, y, z, np, l, i, group, ng);
    }
  }

}


int fof(FLOAT *x, FLOAT *y, FLOAT *z, const int &np, const FLOAT &l,
        int *group, int &ng)
{

  if (np <= 1)
    return 1;

  for (int i = 0; i < np; i++)
    group[i] = -1;

  for (int i = 0; i < np; i++) {
    if (group[i] >= 0) continue;
    group[i] = ng;
    fofi(x, y, z, np, l, i, group, ng);
    ng++;
  }

  return 1;
}
