/*****************************************************************************
 *                                                                           *
 * Copyright 2010 Daniel R. Reynolds                                         *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Gray Flux-Limited Diffusion Split Implicit Problem Class 
/  FillRates routine.
/
/  written by: Daniel Reynolds
/  date:       June 2010
/
/  PURPOSE: Fills the photo-ionization and photo-heating arrays used by
/           chemistry and cooling routines using time-averaged internal
/           values for the radiation.
/ 
/  NOTE: In order to save on memory, the photo-heating rates are 
/        combined into a single rate, and scaled by the current number 
/        density of HI, to be later unpacked by rescaling back with HI.  
/        This loses accuracy in the case that during chemistry 
/        subcycling the chemistry changes significantly, since we retain 
/        the initial rate scaling but use updated HI values in rescaling.
/
************************************************************************/
#ifdef TRANSFER
#include "gFLDSplit.h"


int gFLDSplit::FillRates(EnzoVector *u, EnzoVector *u0, float *phHI, 
			 float *phHeI, float *phHeII, float *photogamma, 
			 float *dissH2I)
{

//   if (debug)
//     printf("Entering gFLDSplit::FillRates routine\n");

  // get local mesh description
  int usz[4], ghXl, ghXr, ghYl, ghYr, ghZl, ghZr;
  u->size(&usz[0], &usz[1], &usz[2], &usz[3], 
	  &ghXl, &ghXr, &ghYl, &ghYr, &ghZl, &ghZr);
  if (usz[0] != LocDims[0]) 
    ENZO_FAIL("FillRates error: x0 vector dims do not match");
  if (usz[1] != LocDims[1]) 
    ENZO_FAIL("FillRates error: x1 vector dims do not match");
  if (usz[2] != LocDims[2]) 
    ENZO_FAIL("FillRates error: x2 vector dims do not match");
  if (usz[3] != (2+Nchem)) 
    ENZO_FAIL("FillRates error: nspecies dims do not match");
  if ((usz[0]+ghXl+ghXr) != ArrDims[0]) 
    ENZO_FAIL("FillRates error: x0 vector sizes do not match");
  if ((usz[1]+ghYl+ghYr) != ArrDims[1]) 
    ENZO_FAIL("FillRates error: x1 vector sizes do not match");
  if ((usz[2]+ghZl+ghZr) != ArrDims[2]) 
    ENZO_FAIL("FillRates error: x2 vector sizes do not match");

  // set some physical constants
  float c = 2.99792458e10;        // speed of light [cm/s]
  float hp = 6.6260693e-27;       // Planck's constant [ergs*s]
  float mp = 1.67262171e-24;      // mass of a proton [g]
  float ev2erg = 1.60217653e-12;  // conversion constant from eV to ergs
  float dom = DenUnits*a*a*a/mp;
  float tbase1 = TimeUnits;
  float xbase1 = LenUnits/a/aUnits;
  float dbase1 = DenUnits*a*a*a*aUnits*aUnits*aUnits;
  float coolunit = aUnits*aUnits*aUnits*aUnits*aUnits * xbase1*xbase1
    * mp*mp / tbase1/tbase1/tbase1 / dbase1;
  float rtunits = ev2erg/TimeUnits/coolunit/dom;
  // float ErUn = ErUnits;  // original
  float ErUn = (ErUnits+ErUnits0)*0.5;   // arithmetic mean
  // float ErUn = sqrt(ErUnits*ErUnits0);   // geometric mean
  // float ErUn = 2.0*ErUnits*ErUnits0/(ErUnits+ErUnits0);  // harmonic mean
  
  // access radiation energy density array
  float *Er = u->GetData(0);

  // compute the size of the fields
  int i, dim, size=1;
  for (dim=0; dim<rank; dim++)  size *= ArrDims[dim];

  // fill HI photo-ionization rate
  float pHIconst = c*TimeUnits*intSigESigHInu/hp/intSigE;
  for (i=0; i<size; i++)  phHI[i] = Er[i]*ErUn*pHIconst;

  // fill HeI and HeII photo-ionization rates
  float pHeIconst  = c*TimeUnits*intSigESigHeInu/hp/intSigE;
  float pHeIIconst = c*TimeUnits*intSigESigHeIInu/hp/intSigE;
  if (RadiativeTransferHydrogenOnly == FALSE) {
    for (i=0; i<size; i++)  phHeI[i]  = Er[i]*ErUn*pHeIconst;
    for (i=0; i<size; i++)  phHeII[i] = Er[i]*ErUn*pHeIIconst;
  }
   
  // fill photo-heating rate
  float phScale    = c*TimeUnits/intSigE/VelUnits/VelUnits/mp/rtunits;
  float GHIconst   = phScale*(intSigESigHI   - 13.6*ev2erg/hp*intSigESigHInu);
  float GHeIconst  = phScale*(intSigESigHeI  - 24.6*ev2erg/hp*intSigESigHeInu);
  float GHeIIconst = phScale*(intSigESigHeII - 54.4*ev2erg/hp*intSigESigHeIInu);
  if (Nchem == 1)
    for (i=0; i<size; i++)  photogamma[i] = Er[i]*ErUn*GHIconst;
  if (Nchem == 3) {
    float *HI = U0->GetData(2);
    float *HeI = U0->GetData(3);
    float *HeII = U0->GetData(4);
    for (i=0; i<size; i++)  photogamma[i] = Er[i]*ErUn *
	      (GHIconst*HI[i] + GHeIconst*HeI[i] + GHeIIconst*HeII[i])/HI[i];
  }

  // fill H2 dissociation rate (none for grey FLD problems)
  if (MultiSpecies > 1) 
    for (i=0; i<size; i++)  dissH2I[i] = 0.0;

  // return success
  return SUCCESS;
}
#endif
