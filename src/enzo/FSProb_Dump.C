/*****************************************************************************
 *                                                                           *
 * Copyright 2009 Daniel R. Reynolds                                         *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Free-streaming Radiation Implicit Problem Class
/  Dump routine
/  
/  written by: Daniel Reynolds
/  date:       March, 2009
/  modified:   
/
/  PURPOSE: Outputs entire problem state to stderr and files.  Useful 
/           upon solve failure.
/
************************************************************************/
#ifdef TRANSFER
#include "FSProb.h"


int FSProb::Dump(EnzoVector *ucur)
{

  if (debug) 
    fprintf(stderr,"\n\nDumping FSProb:\n");

  if (debug) {
    fprintf(stderr,"  told = %g\n",told);
    fprintf(stderr,"  tnew = %g\n",tnew);
    fprintf(stderr,"  dt = %g\n",dt);
    fprintf(stderr,"  theta = %g\n",theta);
    fprintf(stderr,"  kappa = %g\n",kappa);
    fprintf(stderr,"  NGammaDot = %g\n", NGammaDot);
    fprintf(stderr,"  EtaRadius = %g\n", EtaRadius);
    fprintf(stderr,"  EtaCenter = %g %g %g\n", 
	    EtaCenter[0], EtaCenter[1], EtaCenter[2]);
    fprintf(stderr,"  LimType = %"ISYM"\n",LimType);
    fprintf(stderr,"  dt_suggest = %g\n",dt_suggest);
    fprintf(stderr,"  aUnits = %g\n",aUnits);
    fprintf(stderr,"  EUnits = %g\n",EUnits);
    fprintf(stderr,"  EScale = %g\n",EScale);
    fprintf(stderr,"  LenUnits = %g\n",LenUnits);
    fprintf(stderr,"  TimeUnits = %g\n",TimeUnits);
  }

  if (debug) {
    fprintf(stderr,"Dumping FS module parameters to file FSdump.params\n");
    FILE *fptr = fopen("FSdump.params", "w");
    this->WriteParameters(fptr);
    fclose(fptr);
  }

  char *ofile = new char[12];
  strcpy(ofile,"Ef.vec");
  if (debug) 
    fprintf(stderr,"  writing FS radiation to %s\n",ofile);
  sol->writeall(ofile,0);

  return SUCCESS;
}
#endif
