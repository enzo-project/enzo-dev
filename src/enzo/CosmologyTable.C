/***********************************************************************
/
/  WRITES COSMOLOGY PARAMETERS TO AN OUTPUT FILE
/
/  written by: Britton Smith
/  date:       November, 2017
/  modified1:
/
/  PURPOSE: get a(t) from a table
/
/  NOTE:
/
************************************************************************/

#include <string.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "CosmologyParameters.h"

FLOAT myageint(FLOAT a)
{
  FLOAT a2 = a * a;
  float OmegaCurvatureNow = 1 - OmegaMatterNow -
    OmegaLambdaNow - OmegaRadiationNow;
  return pow(((OmegaMatterNow / a) + (OmegaRadiationNow / a2) +
              (OmegaLambdaNow * a2) + OmegaCurvatureNow), -0.5);
}

int InitializeCosmologyTable()
{

  if ((CosmologyTableLoga != NULL) || (CosmologyTableLogt != NULL)) {
    ENZO_FAIL("Cosmology tables have already been initialized.\n");
  }

  CosmologyTableLoga = new FLOAT[CosmologyTableNumberOfBins];
  CosmologyTableLogt = new FLOAT[CosmologyTableNumberOfBins];

  FLOAT loga_i, dloga, coef, a_1, a_12, a_2, dadt_1, dadt_12, dadt_2;
  loga_i = CosmologyTableLogaInitial;
  dloga = (CosmologyTableLogaFinal - loga_i) /
    (CosmologyTableNumberOfBins - 1);

  // Use Simpson's Rule to integrate
  a_1 = pow(10, (loga_i - dloga));
  dadt_1 = myageint(a_1);
  for (int i = 0; i < CosmologyTableNumberOfBins; i++) {
    CosmologyTableLoga[i] = i * dloga + loga_i;
    a_2 = pow(10, (i * dloga + loga_i));
    a_12  = 0.5 * (a_1 + a_2);
    coef = (a_2 - a_1) / 6.;
    dadt_2  = myageint(a_2);
    dadt_12 = myageint(a_12);

    CosmologyTableLogt[i] = coef * (dadt_1 + 4 * dadt_12 + dadt_2);
    if (i > 0) {
      CosmologyTableLogt[i] += CosmologyTableLogt[i-1];
    }

    a_1 = a_2;
    dadt_1 = dadt_2;
  }

  for (int i = 1; i < CosmologyTableNumberOfBins; i++) {
    CosmologyTableLogt[i] = log10(CosmologyTableLogt[i]);
  }

  return SUCCESS;
}
