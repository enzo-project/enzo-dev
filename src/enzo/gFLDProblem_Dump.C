/*****************************************************************************
 *                                                                           *
 * Copyright 2006 Daniel R. Reynolds                                         *
 * Copyright 2006 Laboratory for Computational Astrophysics                  *
 * Copyright 2006 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Gray Flux-Limited Diffusion Implicit Problem Class dump routine
/  
/
/  written by: Daniel Reynolds
/  date:       April, 2007
/  modified1:  
/
/  PURPOSE: Outputs entire problem state to stderr and files.  Useful 
/           upon solve failure.
/
************************************************************************/
#ifdef TRANSFER
#include "gFLDProblem.h"


int gFLDProblem::Dump(EnzoVector *ucur)
{

  if (debug) 
    fprintf(stderr,"\n\nDumping gFLDProblem:\n");

  if (debug) {
    fprintf(stderr,"  told = %g\n",told);
    fprintf(stderr,"  tnew = %g\n",tnew);
    fprintf(stderr,"  dt = %g\n",dt);
    fprintf(stderr,"  theta = %g\n",theta);
    fprintf(stderr,"  LimType = %"ISYM"\n",LimType);
    fprintf(stderr,"  Nchem = %"ISYM"\n",Nchem);
    fprintf(stderr,"  Model = %"ISYM"\n",Model);
    fprintf(stderr,"  AnalyticChem = %"ISYM"\n",AnalyticChem);
    fprintf(stderr,"  AprxJacobian = %"ISYM"\n",approx_jac);
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
    fprintf(stderr,"  TempUnits = %g\n",TempUnits);
    fprintf(stderr,"  VelUnits = %g\n",VelUnits);
    fprintf(stderr,"  MassUnits = %g\n",MassUnits);
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

    rmstmp = rhs0->rmsnorm_component(ns);
    inftmp = rhs0->infnorm_component(ns);
    if (debug) 
      fprintf(stderr,"  rhs0(%"ISYM"): rms = %g, max = %g\n",ns,rmstmp,inftmp);

    rmstmp = rhs->rmsnorm_component(ns);
    inftmp = rhs->infnorm_component(ns);
    if (debug) 
      fprintf(stderr,"  rhs(%"ISYM"): rms = %g, max = %g\n",ns,rmstmp,inftmp);

    for (int ns2=0; ns2<=Nchem+1; ns2++) {
      rmstmp = (L[ns])->rmsnorm_component(ns2);
      inftmp = (L[ns])->infnorm_component(ns2);
      if (debug)  fprintf(stderr,"  Lblock(%"ISYM",%"ISYM"): rms = %g, max = %g\n",
			  ns,ns2,rmstmp,inftmp);
    }
  }
  

  char *ofile = new char[12];
  char *tmp_str = new char[3];
  char *tmp_str2 = new char[3];
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

    strcpy(ofile,"rhs0_");
    strcat(strcat(ofile,tmp_str),".vec");
    if (debug) 
      fprintf(stderr,"  writing rhs0(%"ISYM") to %s\n",ns,ofile);
    rhs0->writeall(ofile,ns);

    strcpy(ofile,"rhs_");
    strcat(strcat(ofile,tmp_str),".vec");
    if (debug) 
      fprintf(stderr,"  writing rhs(%"ISYM") to %s\n",ns,ofile);
    rhs->writeall(ofile,ns);
 
    for (int ns2=0; ns2<=Nchem+1; ns2++) {
      sprintf(tmp_str2,"%"ISYM,ns2);
      strcpy(ofile,"Lblock_");
      strcat(ofile,tmp_str);
      strcat(ofile,tmp_str2);
      strcat(ofile,".vec");
      if (debug) 
	fprintf(stderr,"  writing Lblock(%"ISYM",%"ISYM") to %s\n",ns,ns2,ofile);
      (L[ns])->writeall(ofile,ns2);
    }
  }

  // output opacities (create temporary vector to do output)
  int empty = 1;
  EnzoVector Opacities = EnzoVector(LocDims[0], LocDims[1], LocDims[2], 
				    GhDims[0][0], GhDims[0][1], GhDims[1][0], 
				    GhDims[1][1], GhDims[2][0], GhDims[2][1], 
				    3, NBors[0][0], NBors[0][1], NBors[1][0], 
				    NBors[1][1], NBors[2][0], NBors[2][1], empty);
  Opacities.SetData(0, OpacityP);
  Opacities.SetData(1, OpacityE);
  rmstmp = Opacities.rmsnorm_component(0);
  inftmp = Opacities.infnorm_component(0);
  if (debug) 
    fprintf(stderr,"\n  Planck mean opacity: rms = %g, max = %g\n",rmstmp,inftmp);
  rmstmp = Opacities.rmsnorm_component(1);
  inftmp = Opacities.infnorm_component(1);
  if (debug) 
    fprintf(stderr,"  Energy mean opacity: rms = %g, max = %g\n",rmstmp,inftmp);

  strcpy(ofile,"opacityP");
  strcat(ofile,".vec");
  if (debug) 
    fprintf(stderr,"  writing Planck mean opacity to %s\n",ofile);
  Opacities.writeall(ofile,0);
  strcpy(ofile,"opacityE");
  strcat(ofile,".vec");
  if (debug) 
    fprintf(stderr,"  writing Energy mean opacity to %s\n",ofile);
  Opacities.writeall(ofile,1);

  // output current residual (store in tmp3)
  this->nlresid(tmp3,ucur);
  for (int ns=0; ns<=Nchem+1; ns++) {
    rmstmp = tmp3->rmsnorm_component(ns);
    inftmp = tmp3->infnorm_component(ns);
    if (debug)
      fprintf(stderr,"  f(%"ISYM"): rms = %g, max = %g\n",ns,rmstmp, inftmp);
    sprintf(tmp_str,"%"ISYM,ns);
    strcpy(ofile,"fu_");
    strcat(strcat(ofile,tmp_str),".vec");
    if (debug) 
      fprintf(stderr,"  writing f(%"ISYM") to %s\n",ns,ofile);
    tmp3->writeall(ofile,ns);
  }
  delete[] ofile;
  delete[] tmp_str;

  return SUCCESS;
}
#endif
