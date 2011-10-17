/***********************************************************************
/
/  GRID CLASS (Initialize gadget cooling information)
/
/  written by: Brian O'Shea
/  date:       June 2002
/  modified1:
/
/  PURPOSE:
/    calculate tables of cooling rates vs. temperature for
/    various species
/
/  RETURNS:
/    void - no returns
/
/  NOTE:  This routine is based on code in the GADGET code, written
/    by Volker Springel and Naoki Yoshida.  The subroutines are
/    essentially identical to the Gadget cooling routines, with an
/    Enzo-friendly wrapper around them.  As of the addition of this
/    comment to the code (Feb. 2004) the cooling code is private and
/    NOT FOR PUBLIC RELEASE until Volker makes the cooling code public
/    in Gadget.  See the Gadget web page for more details:
/
/      http://www.mpa-garching.mpg.de/~volker/gadget/
/
/      --B.W.O'Shea, Feb. 2004
/
************************************************************************/

#include <math.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"
#include "Gadget.h"


void GadgetMakeCoolingTable(void)
{

  int i;
  float Tfact;
  float T;

  /* Tmin, Tmax are somewhat arbitrary, but for logged temperature,
     these are reasonable bounds */
  Tmin=0.1; 
  Tmax=9.0;
  XH = HYDROGEN_MASSFRAC;

  yhelium = (1 - XH) / (4 * XH);

  mhboltz = PROTONMASS / BOLTZMANN;


  /* minimum internal energy for neutral gas */

  if(MINGASTEMP > 0.0)
    Tmin = log10(0.1 * MINGASTEMP);
  else
    Tmin = 1.0;

  deltaT = (Tmax - Tmin) / NCOOLTAB;

  ethmin = POW(10.0, Tmin) * (1. + yhelium) / ((1. + 4. * yhelium) *
           mhboltz * GAMMA_MINUS1);

  for(i = 0; i <= NCOOLTAB; i++)
    {
      BetaH0[i] =
	BetaHep[i] =
	Betaff[i] =
	AlphaHp[i] = AlphaHep[i] = AlphaHepp[i] = Alphad[i] = GammaeH0[i] = GammaeHe0[i] = GammaeHep[i] = 0;

      T = POW(10.0, Tmin + deltaT * i);

      Tfact = 1.0 / (1 + sqrt(T / 1.0e5));

      if(118348 / T < 70)
	BetaH0[i] = 7.5e-19 * exp(-118348 / T) * Tfact;

      if(473638 / T < 70)
	BetaHep[i] = 5.54e-17 * POW(T, -0.397) * exp(-473638 / T) * Tfact;

      Betaff[i] = 1.43e-27 * sqrt(T) * (1.1 + 0.34 * exp(-(5.5 - log10(T)) * (5.5 - log10(T)) / 3));

      AlphaHp[i] = 8.4e-11 * POW(T / 1000, -0.2) / (1. + POW(T / 1.0e6, 0.7)) / sqrt(T);

      AlphaHep[i] = 1.5e-10 * POW(T, -0.6353);

      AlphaHepp[i] = 4. * AlphaHp[i];

      if(470000 / T < 70)
	Alphad[i] = 1.9e-3 * POW(T, -1.5) * exp(-470000 / T) * (1. + 0.3 * exp(-94000 / T));

      if(157809.1 / T < 70)
	GammaeH0[i] = 5.85e-11 * sqrt(T) * exp(-157809.1 / T) * Tfact;

      if(285335.4 / T < 70)
	GammaeHe0[i] = 2.38e-11 * sqrt(T) * exp(-285335.4 / T) * Tfact;

      if(631515.0 / T < 70)
	GammaeHep[i] = 5.68e-12 * sqrt(T) * exp(-631515.0 / T) * Tfact;
    }
}
