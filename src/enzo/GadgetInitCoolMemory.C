/***********************************************************************
/
/  GRID CLASS (SET UP FOR GADGET COOLING)
/
/  written by: Brian O'Shea
/  date:       June 2002
/  modified1:
/
/  PURPOSE:
/    allocates memory for the cooling table data
/
/  RETURNS:
/    nothing - void
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

void GadgetInitCoolMemory(void)
{

  if(debug) printf("in GadgetInitCoolMemory\n");

  /* the executive decision was made to use 'new' rather
     than 'malloc' because it's easier.  As of 27 June 2002
     we don't actually clean up after this - ie, no 'delete' - 
     because it should only be called once, and the memory
     will be freed when the simulation finishes anyhow.
  */

  BetaH0 = new float [(NCOOLTAB + 1)];
  BetaHep = new float [(NCOOLTAB + 1)];
  AlphaHp = new float [(NCOOLTAB + 1)];
  AlphaHep = new float [(NCOOLTAB + 1)];
  Alphad = new float [(NCOOLTAB + 1)];
  AlphaHepp = new float [(NCOOLTAB + 1)];
  GammaeH0 = new float [(NCOOLTAB + 1)];
  GammaeHe0 = new float [(NCOOLTAB + 1)];
  GammaeHep = new float [(NCOOLTAB + 1)];
  Betaff = new float [(NCOOLTAB + 1)];
  
}
