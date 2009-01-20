/***********************************************************************
/
/  Initialize GADGET Equilibrium cooling data
/
/  written by: Brian O'Shea
/  date:       September 2003
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
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

#include <string.h>
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "CosmologyParameters.h"
#include "Gadget.h"


int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
void GadgetInitCoolMemory(void);
void GadgetMakeCoolingTable(void);
void GadgetReadIonizeParams(char *fname);
void GadgetIonizeParamsTable(float redshift);

int InitializeGadgetEquilibriumCoolData(FLOAT Time)
{

  if(debug) fprintf(stderr,"in InitializeGadgetEquilibriumCoolData\n");

  FLOAT a = 1.0, dadt, redshift;
  
  if (CosmologyComputeExpansionFactor(Time, &a, &dadt) 
      == FAIL) {
    fprintf(stderr, "Error in CosmologyComputeExpansionFactors.\n");
    return FAIL;
  }

  redshift = 1.0 / a - 1.0;

  GadgetInitCoolMemory();
  GadgetMakeCoolingTable();  
  GadgetReadIonizeParams("TREECOOL");
  GadgetIonizeParamsTable(redshift);
  
  if(debug) printf("InitializeGadgetEquilibriumCoolData:  Cooling tables initialized.  Done\n");

  return SUCCESS;
}
