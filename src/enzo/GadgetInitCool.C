/***********************************************************************
/
/  GRID CLASS (SET UP GADGET COOLING STUFF)
/
/  written by: Brian O'Shea
/  date:       June 2002
/  modified1:
/
/  PURPOSE:
/      initializes the cooling tables and some global variables 
/      for gadget equilibrium cooling
/
/  RETURNS:
/    doesn't return anything - void
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

void GadgetInitCoolMemory(void);
void GadgetMakeCoolingTable(void);
void GadgetReadIonizeParams(char *fname);
void GadgetIonizeParams(float redshift);

void grid::GadgetInitCool(float redshift)
{

  if(debug) printf("in GadgetInitCool\n");
  
  GadgetInitCoolMemory();
  GadgetMakeCoolingTable();
  
#ifdef IONIZETABLE
  GadgetReadIonizeParams("TREECOOL");
#endif 
  GadgetIonizeParams(redshift);
  
  if(debug) printf("GadgetInitCool:  Cooling tables initialized.  Done\n");

}
