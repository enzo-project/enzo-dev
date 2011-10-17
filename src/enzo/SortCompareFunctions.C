/***********************************************************************
/
/  FUNCTIONS FOR QUICKSORT CALLS
/
/  written by: John Wise
/  date:       May, 2009
/  modified:   
/
/  PURPOSE: 
/
************************************************************************/

#include <stdlib.h>
#include "macros_and_parameters.h"
#include "typedefs.h"

Eint32 compare_grid(const void *a, const void *b)
{
  struct particle_data *ia = (struct particle_data*) a;
  struct particle_data *ib = (struct particle_data*) b;
  if (ia->grid - ib->grid < 0)
    return -1;
  else if (ia->grid - ib->grid > 0)
    return 1;
  return 0;
}

/***********************************************************************/

Eint32 compare_proc(const void *a, const void *b)
{
  struct particle_data *ia = (struct particle_data*) a;
  struct particle_data *ib = (struct particle_data*) b;
  if (ia->proc - ib->proc < 0)
    return -1;
  else if (ia->proc - ib->proc > 0)
    return 1;
  return 0;
}

/************************************************************************/

Eint32 compare_star_grid(const void *a, const void *b)
{
  struct star_data *ia = (struct star_data*) a;
  struct star_data *ib = (struct star_data*) b;
  if (ia->grid - ib->grid < 0)
    return -1;
  else if (ia->grid - ib->grid > 0)
    return 1;
  return 0;
}
/***********************************************************************/

Eint32 compare_star_proc(const void *a, const void *b)
{
  struct star_data *ia = (struct star_data*) a;
  struct star_data *ib = (struct star_data*) b;
  if (ia->proc - ib->proc < 0)
    return -1;
  else if (ia->proc - ib->proc > 0)
    return 1;
  return 0;
}

/***********************************************************************/

Eint32 compare_flt(const void *a, const void *b)
{
  if (*(float*)a - *(float*)b < 0)
    return -1;
  else if (*(float*)a - *(float*)b > 0)
    return 1;
  return 0;
}

/***********************************************************************/

Eint32 compare_hkey(const void *a, const void *b)
{
  struct hilbert_data *ia = (struct hilbert_data*) a;
  struct hilbert_data *ib = (struct hilbert_data*) b;
  if (ia->hkey - ib->hkey < 0)
    return -1;
  else if (ia->hkey - ib->hkey > 0)
    return 1;
  return 0;
}
