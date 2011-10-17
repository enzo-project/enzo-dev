/***********************************************************************
/
/  GRID CLASS (initialize gadget cooling routines)
/
/  written by: Brian O'Shea
/  date:       June 2002
/  modified1:
/
/  PURPOSE:
/    given the name of a file, read the ionization parameters
/    and store them in an array
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

void GadgetReadIonizeParams(char *fname)
{

  int i;
  FILE *fdcool;

  if(!(fdcool = fopen(fname, "r")))
    {
      fprintf(stderr," Cannot read ionization table in file `%s'\n", fname);
      /*  endrun(456); */  /* need to add appropriate exit/error checking stuff */
    }

  for(i = 0; i < TABLESIZE; i++)
    gH0[i] = 0;

  for(i = 0; i < TABLESIZE; i++)
    if(fscanf(fdcool, "%"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM,
	      &inlogz[i], &gH0[i], &gHe[i], &gHep[i], &eH0[i], &eHe[i], &eHep[i]) == EOF)
      break;

  fclose(fdcool);

  /*  nheattab is the number of entries in the table */

  for(i = 0, nheattab = 0; i < TABLESIZE; i++)
    if(gH0[i] != 0.0)
      nheattab++;
    else
      break;

  //if(ProcessorNumber == 0){

  if (MyProcessorNumber == ROOT_PROCESSOR){

    printf("GadgetReadIonizeParams: read ionization table with %"ISYM" entries in file `%s'.\n\n", nheattab, fname);
    for(i=0; i < TABLESIZE ; i++)
      //      printf("%"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM,
      printf("%e %e %e %e %e %e %e\n",
	     inlogz[i], gH0[i], gHe[i], gHep[i], eH0[i], eHe[i], eHep[i]);
    fflush(stdout);
  }
   /*  need to put this into the files as well...  */

}
