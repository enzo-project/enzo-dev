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


void GadgetIonizeParamsTable(float redshift);
void GadgetIonizeParamsFunction(float redshift);

void GadgetIonizeParams(float redshift)
{

#ifdef IONIZETABLE
  GadgetIonizeParamsTable(redshift);
#else
  GadgetIonizeParamsFunction(redshift);
#endif

}
