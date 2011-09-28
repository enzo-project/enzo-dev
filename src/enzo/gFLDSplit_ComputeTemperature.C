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
/  Temperature calculation routine.
/
/  written by: Daniel Reynolds
/  date:       July 2009
/
/  modified1:  Elizabeth Tasker, May 2011: DEFAULT_MU changed to parameter Mu
/
/  PURPOSE: Takes in EnzoVector and returns temperature field at 
/           desired time.  This routine will be called repeatedly, so 
/           it should NOT result in any net allocation of memory.
/           
/           Note: in this module, we only use the Temperature for 
/           computing the emissivity in the LTE regime.  To disable 
/           LTE-based emissivity in the radiation field, we set the 
/           'Temperature' to 0 for non-LTE physics.
/
************************************************************************/
#ifdef TRANSFER
#include "gFLDSplit.h"



/* default constants */
#define MIN_TEMP 1.0     // minimum temperature [K]


int gFLDSplit::ComputeTemperature(float *TempArr, EnzoVector *u) 
{

//   if (debug)
//     printf("Entering gFLDSplit::ComputeTemperature routine\n");

  // get local mesh description
  int usz[4], ghXl, ghXr, ghYl, ghYr, ghZl, ghZr;
  u->size(&usz[0], &usz[1], &usz[2], &usz[3], 
	  &ghXl, &ghXr, &ghYl, &ghYr, &ghZl, &ghZr);
  if (usz[0] != LocDims[0]) 
    ENZO_FAIL("Temperature error: x0 vector dims do not match");
  if (usz[1] != LocDims[1]) 
    ENZO_FAIL("Temperature error: x1 vector dims do not match");
  if (usz[2] != LocDims[2]) 
    ENZO_FAIL("Temperature error: x2 vector dims do not match");
  if (usz[3] != (2+Nchem)) 
    ENZO_FAIL("Temperature error: nspecies dims do not match");
  if ((usz[0]+ghXl+ghXr) != ArrDims[0]) 
    ENZO_FAIL("Temperature error: x0 vector sizes do not match");
  if ((usz[1]+ghYl+ghYr) != ArrDims[1]) 
    ENZO_FAIL("Temperature error: x1 vector sizes do not match");
  if ((usz[2]+ghZl+ghZr) != ArrDims[2]) 
    ENZO_FAIL("Temperature error: x2 vector sizes do not match");


  // set some physical constants
  float mp=1.67262171e-24;    // proton mass [g]
  float kb=1.3806504e-16;     // boltzmann constant [erg/K]

  // extract fluid energy array
  float *ec = u->GetData(1);

  // Compute the size of the fields
  int size=1;
  int i, dim;
  for (dim=0; dim<rank; dim++)  size *= ArrDims[dim];


  // compute 'Temperature' based on LTE/non-LTE physics (see Note at top)
  if (Nchem == 0) {


    ////////////////////////////
    // Compute the internal energy first, storing in Temperature array, and
    // converting to physical units
    if (DualEnergyFormalism) {
      // just set TempArr to the internal energy contained in eh and ec
      for (i=0; i<size; i++)  
	TempArr[i] = ecUnits*ec[i] + VelUnits*VelUnits*eh[i];
    }
    else {
      // compute internal energy by subtracting off kinetic component
      for (i=0; i<size; i++) {
	TempArr[i] = ecUnits*ec[i] + VelUnits*VelUnits*(eh[i] 
	  	   - 0.5*(vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]));
	// correct for near-zero temperature 
	TempArr[i] = max(TempArr[i], tiny_number);
      }
    }
    
    // special case for Lowrie & Edwards radiating shock
    if ( ProblemType == 405 ) {
      for (i=0; i<size; i++)
	TempArr[i] = max(TempArr[i]/2.218056e12/kb*1.60219e-12, MIN_TEMP);
    } 
    // special case for the astrophysical radiating shock
    else if ( ProblemType == 404 ) {
      for (i=0; i<size; i++)
	TempArr[i] = max((Gamma-1.0)*0.5*mp*TempArr[i]/kb, MIN_TEMP);
    } 
    // general LTE case
    else {
      for (i=0; i<size; i++)
	TempArr[i] = max((Gamma-1.0)*Mu*mp*TempArr[i]/kb, MIN_TEMP);
    }
  }
  //  Chemistry case: non-LTE physics => 0 'Temperature'
  else {
    for (i=0; i<size; i++) TempArr[i] = 0.0;
  }

  // return success
  return SUCCESS;
}
#endif
