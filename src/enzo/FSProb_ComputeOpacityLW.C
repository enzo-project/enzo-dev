/*****************************************************************************
 *                                                                           *
 * Copyright 2009 John H. Wise                                               *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Free-streaming Radiation Implicit Problem Class
/  Compute spatially-dependent opacity for Lyman-Werner radiation
/
/  written by: John Wise
/  date:       September, 2009
/  modified:   
/
/  PURPOSE: If desired (kappa_h2on == 1), calculate the opacity to 
/           Lyman-Werner radiation due to self-shielding.
/
************************************************************************/
#ifdef TRANSFER
#include "FSProb.h"

int FSProb::ComputeOpacityLW(float *H2Density)
{

  const float H2ISigma = 3.71e-18;  // cm^-2
  const float mh = 1.673e-24;
  double factor = H2ISigma * (DenUnits/mh);

  int i, j, k, dim, size;

  if (kappa_h2on == 0)
    return SUCCESS;

  for (dim = 0, size = 1; dim < rank; dim++)
    size *= ArrDims[dim];

  float *opacity = kappa->GetData(0);
  for (i = 0; i < size; i++)
    opacity[i] = H2Density[i] * factor;

  return SUCCESS;

}
#endif // TRANSFER
