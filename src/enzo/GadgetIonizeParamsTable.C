/***********************************************************************
/
/  GRID CLASS (GADGET inits stuff)
/
/  written by: Brian O'Shea
/  date:       June 2002
/  modified1:
/
/  PURPOSE:
/     get ionization parameters
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

void GadgetIonizeParamsTable(float redshift)
{

  int i, ilow;
  float logz, dzlow, dzhi;


  J_UV=gJH0=gJHep=gJHe0=epsH0=epsHep=epsHe0=0.;

  if(ComovingCoordinates==0)
    {
      gJHe0 = gJHep = gJH0 = 0;
      epsHe0 = epsHep = epsH0 = 0;
      J_UV = 0;
      return;
    }

  logz = log10(redshift + 1.0);
  ilow = 0;
  for(i = 0; i < nheattab; i++)
    {
      if(inlogz[i] < logz)
	ilow = i;
      else
	break;
    }

  dzlow = logz - inlogz[ilow];
  dzhi = inlogz[ilow + 1] - logz;

  if(logz > inlogz[nheattab - 1] || gH0[ilow] == 0 || gH0[ilow + 1] == 0 || nheattab == 0)
    {
      gJHe0 = gJHep = gJH0 = 0;
      epsHe0 = epsHep = epsH0 = 0;
      J_UV = 0;
      return;
    }
  else
    J_UV = 1.e-21;		/* irrelevant as long as it's not 0 */

  gJH0 = JAMPL * POW(10., (dzhi * log10(gH0[ilow]) + dzlow * log10(gH0[ilow + 1])) / (dzlow + dzhi));
  gJHe0 = JAMPL * POW(10., (dzhi * log10(gHe[ilow]) + dzlow * log10(gHe[ilow + 1])) / (dzlow + dzhi));
  gJHep = JAMPL * POW(10., (dzhi * log10(gHep[ilow]) + dzlow * log10(gHep[ilow + 1])) / (dzlow + dzhi));
  epsH0 = JAMPL * POW(10., (dzhi * log10(eH0[ilow]) + dzlow * log10(eH0[ilow + 1])) / (dzlow + dzhi));
  epsHe0 = JAMPL * POW(10., (dzhi * log10(eHe[ilow]) + dzlow * log10(eHe[ilow + 1])) / (dzlow + dzhi));
  epsHep = JAMPL * POW(10., (dzhi * log10(eHep[ilow]) + dzlow * log10(eHep[ilow + 1])) / (dzlow + dzhi));

  return;
}
