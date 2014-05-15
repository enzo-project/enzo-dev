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
/  Radiation spectrum evaluation routine 
/
/  written by: Daniel Reynolds
/  date:       July 2009
/  modified1:  
/
/  PURPOSE: Takes in a frequency and returns the assumed radiation 
/  energy spectrum at that frequency.
/
************************************************************************/
#ifdef TRANSFER
#include "gFLDSplit.h"

 
 

float gFLDSplit::RadiationSpectrum(float nu)
{

  // set necessary constants
  float h = 6.6260693e-27;          // Planck's constant [ergs*s]
  float kb = 1.3806504e-16;         // Boltzmann's constant [ergs/K]
  float pi = 4.0*atan(1.0);         // pi
  float c = 2.99792458e10;          // speed of light [cm/s]
  float ev2erg = 1.60217653e-12;    // conversion constant from eV to ergs
  float nu0 = hnu0_HI*ev2erg/h;     // ionization threshold of Hydrogen (hz)
  float nu1 = 2.5*nu0;              // ionization of Wolf Reyet stars + HeliumI
  float nu2 = 4.0*nu0;              // ionization threshold of HeliumII
  float nu3 = 100.0*nu0;            // parameter used to characterize PopII SED cutoff frequency
  float GBbrk = 2.5;                // parameters in Ricotti 2002 fig.4 fit
  float GXray = 2.0e-3;             // parameters in Ricotti 2002 fig.4 fit
  float sigma;

  // check that frequency is within the allowed range
  if (nu < nu0) {
    fprintf(stderr,"gFLDSplit::RadiationSpectrum Error: (nu = %g) < (nu0 = %g)\n",
	    nu,nu0);
    return -1.0;
  }

  // evaluate the radiation spectrum based on the internal ESpectrum parameter
  switch (ESpectrum) {

  case 1:
    // T = 1e5 K blackbody spectrum
    sigma = 8.0*pi*h*POW(nu/c,3)/(exp(h*nu/kb/1e5)-1.0);
    break;

  // SED with photons above 4 Ryd truncated
  case 2:
    if (nu < nu1)
      sigma = 1/nu0/POW(nu/nu0, 1.8);
    else if (nu < nu2)
      sigma = 1/nu0/GBbrk/POW(nu/nu0, 1.8);
    else if (nu >= nu2)
      sigma=0;
    break;

  default:
    // simple power law spectrum with power -1.5
    sigma = POW((nu/nu0),-1.5);
  }

  return sigma;
}

#endif
