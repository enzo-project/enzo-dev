/***********************************************************************
/
/  GRID CLASS (GADGET inits stuff)
/
/  written by: Brian O'Shea
/  date:       June 2002
/  modified1:
/
/  PURPOSE:
/     analytically calculate ionization parameters
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

#define UVALPHA         1.0

void grid::GadgetIonizeParamsFunction(float redshift)
{

  int i, nint;
  float gint, eint, t, tinv, fac, eps;
  float at, beta, s;
  float pi;
  float Jold=-1.0;


  J_UV = 0.;
  gJHe0 = gJHep = gJH0 = 0.;
  epsHe0 = epsHep = epsH0 = 0.;


  if(ComovingCoordinates==1)	/* analytically compute params from power law J_nu */
    {  /* ionize params */
      if(redshift >= 6)
	J_UV = 0.;
      else
	{
	  if(redshift >= 3)
	    J_UV = 4e-22 / (1 + redshift);
	  else
	    {
	      if(redshift >= 2)
		J_UV = 1e-22;
	      else
		J_UV = 1.e-22 * POW(3.0 / (1 + redshift), 3.0);
	    }
	}

      if(J_UV == Jold)
	return;


      Jold = J_UV;

      if(J_UV == 0)
	return;


      a0 = 6.30e-18;
      planck = 6.6262e-27;
      ev = 1.6022e-12;
      e0_H = 13.6058 * ev;
      e0_He = 24.59 * ev;
      e0_Hep = 54.4232 * ev;

      gint = 0.0;
      eint = 0.0;
      nint = 5000;
      at = 1. / ((float) nint);

      for(i = 1; i <= nint; i++)
	{
	  t = (float) i;
	  t = (t - 0.5) * at;
	  tinv = 1. / t;
	  eps = sqrt(tinv - 1.);
	  fac = exp(4. - 4. * atan(eps) / eps) / (1. - exp(-2. * M_PI / eps)) * POW(t, UVALPHA + 3.);
	  gint += fac * at;
	  eint += fac * (tinv - 1.) * at;
	}

      gJH0 = a0 * gint / planck;
      epsH0 = a0 * eint * (e0_H / planck);
      gJHep = gJH0 * POW(e0_H / e0_Hep, UVALPHA) / 4.0;
      epsHep = epsH0 * POW((e0_H / e0_Hep), UVALPHA - 1.) / 4.0;

      at = 7.83e-18;
      beta = 1.66;
      s = 2.05;

      gJHe0 = (at / planck) * POW((e0_H / e0_He), UVALPHA) *
	(beta / (UVALPHA + s) + (1. - beta) / (UVALPHA + s + 1));
      epsHe0 = (e0_He / planck) * at * POW(e0_H / e0_He, UVALPHA) *
	(beta / (UVALPHA + s - 1) + (1 - 2 * beta) / (UVALPHA + s) - (1 - beta) / (UVALPHA + s + 1));

      pi = M_PI;
      gJH0 *= 4. * pi * J_UV;
      gJHep *= 4. * pi * J_UV;
      gJHe0 *= 4. * pi * J_UV;
      epsH0 *= 4. * pi * J_UV;
      epsHep *= 4. * pi * J_UV;
      epsHe0 *= 4. * pi * J_UV;
    }
}
