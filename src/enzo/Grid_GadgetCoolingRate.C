/***********************************************************************
/
/  GRID CLASS (SOLVE THE GADGET COOLING/HEATING RATE EQUATIONS)
/
/  written by: Brian O'Shea
/  date:       June 2002
/  modified1:
/
/  PURPOSE:
/              return cooling rate, given temperature, density,
/              redshift, guess at electron density
/
/  RETURNS:
/    heating/cooling rate in cgs units
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

#include <stdio.h>
#include <math.h>

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

void Gadgetfind_abundances_and_rates(float logT, float rho, float *ne_guess);

float grid::GadgetCoolingRate(float logT, float rho, float *nelec, float redshift)
{

  float Lambda, Heat;
  float LambdaExc, LambdaIon, LambdaRec, LambdaFF, LambdaCmptn = 0.0;
  float LambdaExcH0, LambdaExcHep, LambdaIonH0, LambdaIonHe0, LambdaIonHep;
  float LambdaRecHp, LambdaRecHep, LambdaRecHepp, LambdaRecHepd;
  float T;

  if(logT <= Tmin)
    logT = Tmin + 0.5 * deltaT;	/* floor at Tmin */

  nHcgs = XH * rho / PROTONMASS;	/* hydrogen number dens in cgs units */

  if(logT < Tmax)
    {
      Gadgetfind_abundances_and_rates(logT, rho, nelec);

      /* Compute cooling and heating rate (cf KWH Table 1) in units of nH**2 */

      T = POW(10.0, logT);

      LambdaExcH0 = bH0 * ne * nH0;
      LambdaExcHep = bHep * ne * nHep;
      LambdaExc = LambdaExcH0 + LambdaExcHep;	/* excitation */

      LambdaIonH0 = 2.18e-11 * geH0 * ne * nH0;
      LambdaIonHe0 = 3.94e-11 * geHe0 * ne * nHe0;
      LambdaIonHep = 8.72e-11 * geHep * ne * nHep;
      LambdaIon = LambdaIonH0 + LambdaIonHe0 + LambdaIonHep;	/* ionization */

      LambdaRecHp = 1.036e-16 * T * ne * (aHp * nHp);
      LambdaRecHep = 1.036e-16 * T * ne * (aHep * nHep);
      LambdaRecHepp = 1.036e-16 * T * ne * (aHepp * nHepp);
      LambdaRecHepd = 6.526e-11 * ad * ne * nHep;
      LambdaRec = LambdaRecHp + LambdaRecHep + LambdaRecHepp + LambdaRecHepd;

      LambdaFF = bff * (nHp + nHep + 4 * nHepp) * ne;

      Lambda = LambdaExc + LambdaIon + LambdaRec + LambdaFF;

      if(ComovingCoordinates)
	{
	  LambdaCmptn = 5.65e-36 * ne * (T - 2.73 * (1. + redshift)) * POW(1. + redshift, 4.) / nHcgs;

	  Lambda += LambdaCmptn;
	}
      else
	LambdaCmptn = 0;

      Heat = 0;
      if(J_UV != 0)
	Heat += (nH0 * epsH0 + nHe0 * epsHe0 + nHep * epsHep) / nHcgs;
    }
  else				/* here we're outside of tabulated rates, T>Tmax K */
    {
      /* at high T (fully ionized); only free-free and Compton cooling are present.
         Assumes no heating. */

      Heat = 0;

      LambdaExcH0 = LambdaExcHep = LambdaIonH0 = LambdaIonHe0 = LambdaIonHep =
	LambdaRecHp = LambdaRecHep = LambdaRecHepp = LambdaRecHepd = 0;

      /* very hot: H and He both fully ionized */
      nHp = 1.0;
      nHep = 0;
      nHepp = yhelium;
      ne = nHp + 2.0 * nHepp;
      *nelec = ne;		/* note: in units of the hydrogen number density */

      T = POW(10.0, logT);
      LambdaFF =
	1.42e-27 * sqrt(T) * (1.1 + 0.34 * exp(-(5.5 - logT) * (5.5 - logT) / 3)) * (nHp + 4 * nHepp) * ne;

      if(ComovingCoordinates)
	{

	  /* add inverse Compton cooling off the microwave background */
	  LambdaCmptn = 5.65e-36 * ne * (T - 2.73 * (1. + redshift)) * POW(1. + redshift, 4.) / nHcgs;
	}
      else
	LambdaCmptn = 0;

      Lambda = LambdaFF + LambdaCmptn;
    }

  return (Heat - Lambda);

}
