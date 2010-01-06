/************************************************************************
/
/   CALCULATE PHOTO-IONIZATION CROSS-SECTION
/
/   written by: John Wise
/   date: July, 2007
/   modified1:
/
/   PURPOSE: For a given type of photon and energy, calculate the
/            appropriate cross-section.  Uses fits from Verner et al. 
/            (1996).
/
/   INPUTS:  type: photon type
/            energy: photon energy in eV
/
/   RETURNS: cross-section in cm^2
/
/************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "ExternalBoundary.h"
#include "Fluxes.h"
#include "GridList.h"
#include "Grid.h"
#include "CosmologyParameters.h"

FLOAT FindCrossSection(int type, float energy)
{

  float sigma;
  float e_th, e_max, e0, sigma0, ya, P, yw, y0, y1;

  switch (type) {

    // HI
  case 0:
    e_th = 13.6;
    e_max = 5e4;
    e0 = 4.298e-1;
    sigma0 = 5.475e4;
    ya = 32.88;
    P = 2.963;
    yw = y0 = y1 = 0.0;
    break;

    // HeI
  case 1:
    e_th = 24.59;
    e_max = 5e4;
    e0 = 13.61;
    sigma0 = 9.492e2;
    ya = 1.469;
    P = 3.188;
    yw = 2.039;
    y0 = 0.4434;
    y1 = 2.136;
    break;
  
    // HeII
  case 2:
    e_th = 54.42;
    e_max = 5e4;
    e0 = 1.720;
    sigma0 = 1.369e4;
    ya = 32.88;
    P = 2.963;
    yw = y0 = y1 = 0.0;
    break;

    // Lyman-Werner
  case 3:
    e_th = 11.18;
    e_max = 13.6;
    sigma0 = 3.71;  
    break;

  } // ENDSWITCH

  float x, y, fy;
  if (type != 3) {
    x = energy/e0 - y0;
    y = sqrt(x*x + y1*y1);
    fy = ((x-1.0)*(x-1.0) + yw*yw) * pow(y, 0.5*P-5.5) * 
      pow((1.0 + sqrt(y/ya)), -P);
  } else
    fy = 1.0;

  sigma = sigma0 * fy * 1e-18;

  return sigma;

}
