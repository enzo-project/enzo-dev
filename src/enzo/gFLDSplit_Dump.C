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
/  Gray Flux-Limited Diffusion Split Implicit Problem Class 
/  Dump routine
/  
/
/  written by: Daniel Reynolds
/  date:       July 2009
/  modified1:  
/
/  PURPOSE: Outputs entire problem state to stderr and files.  Useful 
/           upon solve failure.
/
************************************************************************/
#ifdef TRANSFER
#include "gFLDSplit.h"


int gFLDSplit::Dump(EnzoVector *ucur)
{

  if (debug) 
    fprintf(stderr,"\n\nDumping gFLDSplit:\n");

  if (debug) {
    fprintf(stderr,"  told = %g\n",told);
    fprintf(stderr,"  tnew = %g\n",tnew);
    fprintf(stderr,"  dt = %g\n",dt);
    fprintf(stderr,"  dtrad = %g\n",dtrad);
    fprintf(stderr,"  dtchem = %g\n",dtchem);
    fprintf(stderr,"  maxsubcycles = %g\n",maxsubcycles);
    fprintf(stderr,"  maxchemsub = %g\n",maxchemsub);
    fprintf(stderr,"  theta = %g\n",theta);
    fprintf(stderr,"  Nchem = %"ISYM"\n",Nchem);
    fprintf(stderr,"  Model = %"ISYM"\n",Model);
    fprintf(stderr,"  ESpectrum = %"ISYM"\n",ESpectrum);
    fprintf(stderr,"  aUnits = %g\n",aUnits);
    fprintf(stderr,"  ErUnits = %g\n",ErUnits);
    fprintf(stderr,"  ecUnits = %g\n",ecUnits);
    fprintf(stderr,"  NiUnits = %g\n",NiUnits);
    fprintf(stderr,"  ErScale = %g\n",ErScale);
    fprintf(stderr,"  ecScale = %g\n",ecScale);
    fprintf(stderr,"  NiScale = %g\n",NiScale);
    fprintf(stderr,"  DenUnits = %g\n",DenUnits);
    fprintf(stderr,"  LenUnits = %g\n",LenUnits);
    fprintf(stderr,"  TimeUnits = %g\n",TimeUnits);
    fprintf(stderr,"  VelUnits = %g\n",VelUnits);
  }

  if (debug) {
    fprintf(stderr,"Dumping FLD module parameters to file RTdump.params\n");
    FILE *fptr = fopen("RTdump.params", "w");
    this->WriteParameters(fptr);
    fclose(fptr);
  }


  float rmstmp, inftmp;
  for (int ns=0; ns<=Nchem+1; ns++) {
    rmstmp = U0->rmsnorm_component(ns);
    inftmp = U0->infnorm_component(ns);
    if (debug) 
      fprintf(stderr,"\n  U0(%"ISYM"): rms = %g, max = %g\n",ns,rmstmp,inftmp);

    rmstmp = ucur->rmsnorm_component(ns);
    inftmp = ucur->infnorm_component(ns);
    if (debug) 
      fprintf(stderr,"  u(%"ISYM"): rms = %g, max = %g\n",ns,rmstmp,inftmp);
  }
  

  char *ofile = new char[12];
  char *tmp_str = new char[3];
  for (int ns=0; ns<=Nchem+1; ns++) {

    sprintf(tmp_str,"%"ISYM,ns);

    strcpy(ofile,"u0_");
    strcat(strcat(ofile,tmp_str),".vec");
    if (debug) 
      fprintf(stderr,"  writing U0(%"ISYM") to %s\n",ns,ofile);
    U0->writeall(ofile,ns);

    strcpy(ofile,"u_");
    strcat(strcat(ofile,tmp_str),".vec");
    if (debug) 
      fprintf(stderr,"  writing u(%"ISYM") to %s\n",ns,ofile);
    ucur->writeall(ofile,ns);
  }

  // // output opacities (create temporary vector to do output)
  // int empty = 1;
  // EnzoVector Opacities = EnzoVector(LocDims[0], LocDims[1], LocDims[2], 
  // 				    GhDims[0][0], GhDims[0][1], GhDims[1][0], 
  // 				    GhDims[1][1], GhDims[2][0], GhDims[2][1], 
  // 				    1, NBors[0][0], NBors[0][1], NBors[1][0], 
  // 				    NBors[1][1], NBors[2][0], NBors[2][1], empty);
  // Opacities.SetData(0, OpacityE);
  // rmstmp = Opacities.rmsnorm_component(0);
  // inftmp = Opacities.infnorm_component(0);
  // if (debug) 
  //   fprintf(stderr,"\n  Energy mean opacity: rms = %g, max = %g\n",rmstmp,inftmp);

  // strcpy(ofile,"opacityE");
  // strcat(ofile,".vec");
  // if (debug) 
  //   fprintf(stderr,"  writing Energy mean opacity to %s\n",ofile);
  // Opacities.writeall(ofile,0);

  return SUCCESS;
}
#endif
