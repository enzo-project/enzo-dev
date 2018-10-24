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
/  Gray Flux-Limited Diffusion Implicit Problem Class 
/  Radiation Spectrum Evaluation routine 
/
/  written by: Daniel Reynolds
/  date:       March, 2007
/  modified1:  
/
/  PURPOSE: Takes in a frequency and returns the assumed radiation 
/  energy spectrum at that frequency.
/
************************************************************************/
#ifdef TRANSFER
#include "gFLDProblem.h"
#include "phys_constants.h"


float gFLDProblem::RadiationSpectrum(float nu)
{

  // set necessary constants
  float nu0 = hnu0_HI*erg_eV/h_planck;     // ionization threshold of Hydrogen (hz)
  float sigma;

  // check that frequency is within the allowed range
  if (nu < nu0) {
    fprintf(stderr,"gFLDProblem::RadiationSpectrum Error: (nu = %g) < (nu0 = %g)\n",
	    nu,nu0);
    return -1.0;
  }

  // evaluate the radiation spectrum based on the internal ESpectrum parameter
  switch (ESpectrum) {

  case 1:
    // T = 1e5 K blackbody spectrum
    sigma = 8.0*pi*h_planck*POW(nu/clight,3)/(exp(h_planck*nu/kboltz/1e5)-1.0);
    break;

  // Add new spectrum choices here
  // case 2:
  //   ...
  //   break;

  default:
    // simple power law spectrum with power -1.5
    sigma = POW((nu/nu0),-1.5);
  }

  return sigma;
}

#endif
