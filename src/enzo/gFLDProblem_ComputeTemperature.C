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
/  Gray Flux-Limited Diffusion Implicit Problem Class temperature 
/  calculation routine (nearly identical to 
/  Grid_ComputeTemperatureField, but allows input arguments based on u 
/  fields, as opposed to just extracting all information out of the grid.
/
/  written by: Daniel Reynolds
/  date:       October, 2006
/  modified1:  June 12, 2007 by John Hayes; added case for computing temperature
/              appropriate to Marshak wave test problems of type published by
/              Su & Olson, JCP, 62 (1999).
/  modified2:  September 11, 2007 by John Hayes, added some flexibility in how
/              the mean molecular weight is evaluated so that radiating shock
/              problems can be accommodated properly.
/  modified3:  Elizabeth Tasker, May 2011: changed DEFAULT_MU to parameter Mu
/              (default still 0.6).
/
/  PURPOSE: Takes in EnzoVector and returns temperature field at 
/           desired time.  This routine will be called repeatedly, so 
/           it should NOT result in any net allocation of memory.
/
************************************************************************/
#ifdef TRANSFER
#include "gFLDProblem.h"



/* default constants */
#define MIN_TEMP 1.0     // minimum temperature [K]


int gFLDProblem::ComputeTemperature(float *TempArr, float time, 
				    FLOAT a, EnzoVector *u) 
{

//   if (debug)
//     printf("Entering gFLDProblem::ComputeTemperature routine\n");

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
  float rc=7.56e-15;          // radiation constant [erg/cm^3/K^4]
  float Cv, everg, mmw;


  // extract fluid energy array
  float *ec = u->GetData(1);

  // Compute the size of the fields
  int size=1;
  int i, dim;
  for (dim=0; dim<rank; dim++)  size *= ArrDims[dim];


  float maxval, minval, avgval;


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

  ////////////////////////////
  // compute the temperature, depending on whether we are using chemistry or not
  //   no chemistry: compute using approximation of fully ionized gas
  if (Nchem == 0) {
    // 
    // Check for Su-Olson Marshak model for mat-rad coupling first. If not
    // selected, do ordinary ideal gas. Marshak problems are designated by
    // Model numbers in the range 20 <= Model <= 29.  Currently, we have:
    // 
    // Model 20: 
    //   Grey Su & Olson #1 (Su & Olson, 1996, JQSRT vol 56, #3, pp 337-351)
    //   Grey Su & Olson #2 (Su & Olson, 1997, Ann. Nucl. Energy v.24, pp 1035-1055)
    //   Grey Olson, Auer, & Hall (Olson et. al 2000, JQSRT vol 64, pp 619-634)
    //
    if ( Model >= 20 && Model <= 29) {
      float SuOlsonGreyEps = MarshakParms[0];
      for (i=0; i<size; i++)
        TempArr[i] = sqrt( sqrt( SuOlsonGreyEps*TempArr[i]/rc ) );
    } 
    //  STANDARD temperature computation for fully ionized gas
    else {
      if ( ProblemType == 405 ) {
	// special case for Lowrie & Edwards radiating shock
	Cv    = 2.218056e12;
	everg = 1.60219e-12;
	mmw   = everg / (Gamma-1.0) / Cv / mp;
	Cv    = 2.218056e12 * kb / everg ;
	for (i=0; i<size; i++)
	  TempArr[i] = max(TempArr[i]/Cv, MIN_TEMP);
      } else if ( ProblemType == 404 ) {
	// special case for the astrophysical radiating shock
	mmw = 0.5;
      } else {
	mmw = Mu;
      }
      if ( ProblemType != 405 ) {
	for (i=0; i<size; i++)
	  TempArr[i] = max((Gamma-1.0)*mmw*mp*TempArr[i]/kb, MIN_TEMP);
      }
    }
  } 
  //  Chemistry case: compute temperature with self-consistent MU value
  else {

    float mu, num_density, ne, nHI, nHII, nHeI, nHeII, nHeIII;
    // extract chemistry arrays
    float **ni = new float *[Nchem];
    for (int i=1; i<=Nchem; i++)  ni[i-1] = u->GetData(1+i);

    //   If performing Ionization test 0 or 1, decouple the gas energy from 
    //   the chemistry: use variable Mu (default 0.6) as in standard approach
    if ((Model == 4) || (Model == 5)) {
      for (i=0; i<size; i++)
        TempArr[i] = max((Gamma-1.0)*Mu*mp*TempArr[i]/kb, MIN_TEMP);
    }
    else {
      // Hydrogen only
      if ((Nchem == 1) && (HFrac == 1.0)) {
	for (i=0; i<size; i++) {
	  // compute mean molecular weight
	  nHI = ni[0][i]*NiScale;
	  nHII = rho[i]*HFrac - nHI;
	  ne = nHII;
	  num_density = (nHI + nHII + ne);
	  mu = rho[i]/num_density;
	  
	  // compute temperature
	  TempArr[i] = max((Gamma-1.0)*mu*mp*TempArr[i]/kb, MIN_TEMP);
	}
      }
      // Hydrogen and Helium
      else if (Nchem == 3) {
	for (i=0; i<size; i++) {
	  // compute mean molecular weight
	  nHI = ni[0][i]*NiScale;
	  nHII = rho[i]*HFrac - nHI;
	  nHeI = ni[1][i]*NiScale;
	  nHeII = ni[2][i]*NiScale;
	  nHeIII = rho[i]*(1.0-HFrac) - nHeI - nHeII;
	  ne = nHII + nHeII/4.0 + nHeIII/2.0;
	  num_density = (0.25*(nHeI+nHeII+nHeIII) + nHI + nHII + ne);
	  mu = rho[i]/num_density;
	  
	  // compute temperature
	  TempArr[i] = max((Gamma-1.0)*mu*mp*TempArr[i]/kb, MIN_TEMP);
	}
      }
      // otherwise return error
      else
	ENZO_FAIL("ComputeTemperature error: Nchem must be {0,1,3}");
    }

    // clean up a bit
    delete[] ni;

  } // end chemistry case

  // return success
  return SUCCESS;
}
#endif
