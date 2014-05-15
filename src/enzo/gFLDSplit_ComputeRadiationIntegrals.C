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
/  Radiation Spectrum Evaluation routine 
/
/  written by: Daniel Reynolds
/  date:       July 2009
/  modified1:  
/
/  PURPOSE: Computes the Radiation Energy spectrum integrals
/                  int_{nu0_HI}^{inf} sigma_E(nu) dnu
/                  int_{nu0_HI}^{inf} sigma_E(nu)*sigma_HI(nu) dnu
/                  int_{nu0_HI}^{inf} sigma_E(nu)*sigma_HI(nu)/nu dnu
/                  int_{nu0_HeI}^{inf} sigma_E(nu)*sigma_HeI(nu) dnu
/                  int_{nu0_HeI}^{inf} sigma_E(nu)*sigma_HeI(nu)/nu dnu
/                  int_{nu0_HeII}^{inf} sigma_E(nu)*sigma_HeII(nu) dnu
/                  int_{nu0_HeII}^{inf} sigma_E(nu)*sigma_HeII(nu)/nu dnu
/           where nu0_* is the ionization threshold of the relevant 
/           species, sigma_E(nu) is the spectrum of the grey radiation 
/           energy density, sigma_HI(nu) is the ionization cross section 
/           of HI, sigma_HeI(nu) is the ionization cross section of HeI, 
/           and sigma_HeII(nu) is the ionization cross section of HeII.
/
/           These are computed using a simple 4th-order accurate 
/           numerical quadrature rule (composite Simpson's), on the 
/           re-mapped indefinite integral
/                  int_{nu0}^{inf} f(nu) dnu
/                = int_0^1 nu0/(x^2)*f(nu0/x) dx
/           The composite Simpson's rule computes this integral as 
/           the sum of small quadrature intervals [xl,xr]_i, where 
/                  [0,1] = Union_i [xl,xr]_i,
/           and xm = (xl+xr)/2, via the quadrature formula
/                  int_{xl}^{xr} g(x) dx 
/                = (xr-xl)/6*(g(xl) + 4*g(xm) + g(xr))
/           Note: since 1/0 = infty, we begin the mapped quadrature 
/           integration at a small positive value.
/
/           We further note that these integrals are typically 
/           computed once and stored in the gFLDSplit object, upon 
/           initialization of said object.  We compute these here to 
/           allow for control over the assumed radiation energy 
/           density spectrum, as stored in 
/           gFLDSplit_RadiationSpectrum.C, as opposed to fixing the 
/           spectrum, calculating these same integrals outside of the 
/           simulation, and hard-coding the physics routines to use 
/           the fixed integral values.
/
************************************************************************/
#ifdef TRANSFER
#include "gFLDSplit.h"

 
 

int gFLDSplit::ComputeRadiationIntegrals()
{

  if (debug)  printf("Entering gFLDSplit::ComputeRadiationIntegrals\n");

  // set necessary constants
  float h = 6.6260693e-27;               // Planck's constant [ergs*s]
  float ev2erg = 1.60217653e-12;         // conversion constant from eV to ergs
  float nu0_HI = hnu0_HI*ev2erg/h;       // ionization threshold of HI (hz)
  float nu0_HeI = hnu0_HeI*ev2erg/h;     // ionization threshold of HeI (hz)
  float nu0_HeII = hnu0_HeII*ev2erg/h;   // ionization threshold of HeII (hz)
  float epsilon = 1.0;                   // floating point roundoff
  float nuMoray = 5.0777776e+15;         // 21 eV photon frequency (Hz)
  while ((1.0 + epsilon) > 1.0)  epsilon*=0.5;

  // set integration parameters
  float FreqH = 1e-4;        // normalized width of quadrature bins
  float Llimit, Ulimit;           // limits of integration, shift away 
  Llimit = POW(epsilon,0.35);     // from [0,1] to avoid overflow/underflow
  Ulimit = 1.0-POW(epsilon,0.47); //  errors in calculations.

  // initialize integrals
  intSigE = intSigESigHI = intSigESigHeI = intSigESigHeII = 0.0;
  intSigESigHInu = intSigESigHeInu = intSigESigHeIInu = 0.0;

  // initialize quadrature points and values
  float xl, xm, xr, nu_l, nu_m, nu_r, fl, fm, fr;
  float fl_E, fm_E, fr_E, fl_ni, fm_ni, fr_ni, fl_nu, fm_nu, fr_nu;


  // for monochromatic problems, integration uses delta function
  //    monochromatic at nu0_HI
  if (ESpectrum == -1) {
    // evaluation point
    nu_m = nu0_HI*(1.0 + 2.0*epsilon);

    // integral( delta_{nu0}(nu) )
    intSigE = 1.0;

    // integral( delta_{nu0}(nu) * sig_{HI}(nu) )
    intSigESigHI = this->CrossSections(nu_m,0);

    // integral( delta_{nu0}(nu) * sig_{HI}(nu) / nu )
    intSigESigHInu = this->CrossSections(nu_m,0) / nu_m;

    if (debug) {
      printf("  Computed Radiation Integrals:\n");
      printf("    intSigE          = %22.16e\n",intSigE);
      printf("    intSigESigHI     = %22.16e\n",intSigESigHI);
      printf("    intSigESigHInu   = %22.16e\n",intSigESigHInu);
    }

    // (integrals using sigHeI and sigHeII are zero at nu0_HI)
    return SUCCESS;
  }

  //    monochromatic at nu0_HeI
  if (ESpectrum == -2) {
    // evaluation point
    nu_m = nu0_HeI*(1.0 + 2.0*epsilon);

    // integral( delta_{nu0}(nu) )
    intSigE = 1.0;

    // integral( delta_{nu0}(nu) * sig_{HI}(nu) )
    intSigESigHI = this->CrossSections(nu_m,0);

    // integral( delta_{nu0}(nu) * sig_{HI}(nu) / nu )
    intSigESigHInu = this->CrossSections(nu_m,0) / nu_m;

    // integral( delta_{nu0}(nu) * sig_{HeI}(nu) )
    intSigESigHeI = this->CrossSections(nu_m,1);

    // integral( delta_{nu0}(nu) * sig_{HeI}(nu) / nu )
    intSigESigHeInu = this->CrossSections(nu_m,1) / nu_m;

    if (debug) {
      printf("  Computed Radiation Integrals:\n");
      printf("    intSigE          = %22.16e\n",intSigE);
      printf("    intSigESigHI     = %22.16e\n",intSigESigHI);
      printf("    intSigESigHInu   = %22.16e\n",intSigESigHInu);
      printf("    intSigESigHeI    = %22.16e\n",intSigESigHeI);
      printf("    intSigESigHeInu  = %22.16e\n",intSigESigHeInu);
    }

    // (integrals using sigHeII are zero at nu0_HeI)
    return SUCCESS;
  }

  //    monochromatic at nu0_HeII
  if (ESpectrum == -3) {
    // evaluation point
    nu_m = nu0_HeII*(1.0 + 2.0*epsilon);

    // integral( delta_{nu0}(nu) )
    intSigE = 1.0;

    // integral( delta_{nu0}(nu) * sig_{HI}(nu) )
    intSigESigHI = this->CrossSections(nu_m,0);

    // integral( delta_{nu0}(nu) * sig_{HI}(nu) / nu )
    intSigESigHInu = this->CrossSections(nu_m,0) / nu_m;

    // integral( delta_{nu0}(nu) * sig_{HeI}(nu) )
    intSigESigHeI = this->CrossSections(nu_m,1);

    // integral( delta_{nu0}(nu) * sig_{HeI}(nu) / nu )
    intSigESigHeInu = this->CrossSections(nu_m,1) / nu_m;

    // integral( delta_{nu0}(nu) * sig_{HeII}(nu) )
    intSigESigHeII = this->CrossSections(nu_m,2);

    // integral( delta_{nu0}(nu) * sig_{HeII}(nu) / nu )
    intSigESigHeIInu = this->CrossSections(nu_m,2) / nu_m;

    if (debug) {
      printf("  Computed Radiation Integrals:\n");
      printf("    intSigE          = %22.16e\n",intSigE);
      printf("    intSigESigHI     = %22.16e\n",intSigESigHI);
      printf("    intSigESigHInu   = %22.16e\n",intSigESigHInu);
      printf("    intSigESigHeI    = %22.16e\n",intSigESigHeI);
      printf("    intSigESigHeInu  = %22.16e\n",intSigESigHeInu);
      printf("    intSigESigHeII   = %22.16e\n",intSigESigHeII);
      printf("    intSigESigHeIInu = %22.16e\n",intSigESigHeIInu);
    }
    return SUCCESS;
  }

  // Geoffrey So 2/25/2014 adding 21 eV for Moray/FLD comparison
  if (ESpectrum == -10) {
    // evaluation point
    // changing nu0_HI to 5.0777e15 Hz for 21 eV photons
    nu_m = nuMoray*(1.0 + 2.0*epsilon);

    // integral( delta_{nu0}(nu) )
    intSigE = 1.0;

    // integral( delta_{nu0}(nu) * sig_{HI}(nu) )
    intSigESigHI = this->CrossSections(nu_m,0);

    // integral( delta_{nu0}(nu) * sig_{HI}(nu) / nu )
    intSigESigHInu = this->CrossSections(nu_m,0) / nu_m;

    if (debug) {
      printf("  Computed Radiation Integrals:\n");
      printf("    intSigE          = %22.16e\n",intSigE);
      printf("    intSigESigHI     = %22.16e\n",intSigESigHI);
      printf("    intSigESigHInu   = %22.16e\n",intSigESigHInu);
    }

    // (integrals using sigHeI and sigHeII are zero at nu0_HI)
    return SUCCESS;
  }


  //////////////////////////////////////////////////////////////
  // Compute intSigE, intSigESigHI and intSigESigHInu integrals

  //   left end of integration
  //      can't start at 0, shift over a bit
  xr = Llimit;
  nu_r = nu0_HI/xr;
  //      get function values at this location
  fr_E = this->RadiationSpectrum(nu_r);
  fr_ni = this->CrossSections(nu_r,0);
  fr_nu = 1.0/nu_r;
  if ((fr_E == -1.0) || (fr_ni == -1.0)) 
    ENZO_FAIL("ComputeRadiationIntegrals Error in evaluating spectrum");
  
  //   iterate over intervals
  for (int i=1; i<1e9; i++) {

    //      set quadrature points in interval
    xl = xr;  // cannot start at 0, so shift over a bit
    xr = min(xl+FreqH,Ulimit);
    xm = 0.5*(xl+xr);

    //      copy left subinterval function value, location, etc
    nu_l = nu_r;
    fl_E = fr_E;
    fl_ni = fr_ni;
    fl_nu = fr_nu;

    //      evaluate sigma_E(), sigma_HI, 1/nu at remapped quad. pts.
    nu_m = nu0_HI/xm;
    fm_E = this->RadiationSpectrum(nu_m);
    fm_ni = this->CrossSections(nu_m,0);
    fm_nu = 1.0/nu_m;
    if ((fm_E == -1.0) || (fm_ni == -1.0)) 
      ENZO_FAIL("ComputeRadiationIntegrals Error in evaluating spectrum");
  
    nu_r = nu0_HI/xr;
    fr_E = this->RadiationSpectrum(nu_r);
    fr_ni = this->CrossSections(nu_r,0);
    fr_nu = 1.0/nu_r;
    if ((fr_E == -1.0) || (fr_ni == -1.0)) 
      ENZO_FAIL("ComputeRadiationIntegrals Error in evaluating spectrum");
  
    //      compute integrals using re-scaled function values
    fl = fl_E/xl/xl;
    fm = fm_E/xm/xm;
    fr = fr_E/xr/xr;
    intSigE += nu0_HI*(xr-xl)/6.0*(fl + 4.0*fm + fr);
    
    fl = fl_E*fl_ni/xl/xl;
    fm = fm_E*fm_ni/xm/xm;
    fr = fr_E*fr_ni/xr/xr;  
    intSigESigHI += nu0_HI*(xr-xl)/6.0*(fl + 4.0*fm + fr);
    
    fl = fl_E*fl_ni*fl_nu/xl/xl;
    fm = fm_E*fm_ni*fm_nu/xm/xm;
    fr = fr_E*fr_ni*fr_nu/xr/xr;
    intSigESigHInu += nu0_HI*(xr-xl)/6.0*(fl + 4.0*fm + fr);

    //      quit if we have finished interval
    if (xr >= (Ulimit-10.0*epsilon))  break;
  
  }
  


  //////////////////////////////////////////////////////////////
  // Compute intSigESigHeI and intSigESigHeInu integrals

  //   left end of integration
  //      can't start at 0, shift over a bit
  xr = Llimit;
  nu_r = nu0_HeI/xr;
  //      get function values at this location
  fr_E = this->RadiationSpectrum(nu_r);
  fr_ni = this->CrossSections(nu_r,1);
  fr_nu = 1.0/nu_r;
  if ((fr_E == -1.0) || (fr_ni == -1.0)) 
    ENZO_FAIL("ComputeRadiationIntegrals Error in evaluating spectrum");
  
  //   iterate over intervals
  for (int i=1; i<1e9; i++) {

    //      set quadrature points in interval
    xl = xr;  // cannot start at 0, so shift over a bit
    xr = min(xl+FreqH,Ulimit);
    xm = 0.5*(xl+xr);

    //      copy left subinterval function value, location, etc
    nu_l = nu_r;
    fl_E = fr_E;
    fl_ni = fr_ni;
    fl_nu = fr_nu;

    //      evaluate sigma_E(), sigma_HeI, 1/nu at remapped quad. pts.
    nu_m = nu0_HeI/xm;
    fm_E = this->RadiationSpectrum(nu_m);
    fm_ni = this->CrossSections(nu_m,1);
    fm_nu = 1.0/nu_m;
    if ((fm_E == -1.0) || (fm_ni == -1.0)) 
      ENZO_FAIL("ComputeRadiationIntegrals Error in evaluating spectrum");
  
    nu_r = nu0_HeI/xr;
    fr_E = this->RadiationSpectrum(nu_r);
    fr_ni = this->CrossSections(nu_r,1);
    fr_nu = 1.0/nu_r;
    if ((fr_E == -1.0) || (fr_ni == -1.0)) 
      ENZO_FAIL("ComputeRadiationIntegrals Error in evaluating spectrum");
  
    //      compute integrals using re-scaled function values
    fl = fl_E*fl_ni/xl/xl;
    fm = fm_E*fm_ni/xm/xm;
    fr = fr_E*fr_ni/xr/xr;
    intSigESigHeI += nu0_HeI*(xr-xl)/6.0*(fl + 4.0*fm + fr);
    
    fl = fl_E*fl_ni*fl_nu/xl/xl;
    fm = fm_E*fm_ni*fm_nu/xm/xm;
    fr = fr_E*fr_ni*fr_nu/xr/xr;
    intSigESigHeInu += nu0_HeI*(xr-xl)/6.0*(fl + 4.0*fm + fr);

    //      quit if we have finished interval
    if (xr >= (Ulimit-10.0*epsilon))  break;
  
  }
  

  //////////////////////////////////////////////////////////////
  // Compute intSigESigHeII and intSigESigHeIInu integrals

  //   left end of integration
  //      can't start at 0, shift over a bit
  xr = Llimit;
  nu_r = nu0_HeII/xr;
  //      get function values at this location
  fr_E = this->RadiationSpectrum(nu_r);
  fr_ni = this->CrossSections(nu_r,2);
  fr_nu = 1.0/nu_r;
  if ((fr_E == -1.0) || (fr_ni == -1.0)) 
    ENZO_FAIL("ComputeRadiationIntegrals Error in evaluating spectrum");
  
  //   iterate over intervals
  for (int i=1; i<1e9; i++) {

    //      set quadrature points in interval
    xl = xr;  // cannot start at 0, so shift over a bit
    xr = min(xl+FreqH,Ulimit);
    xm = 0.5*(xl+xr);

    //      copy left subinterval function value, location, etc
    nu_l = nu_r;
    fl_E = fr_E;
    fl_ni = fr_ni;
    fl_nu = fr_nu;

    //      evaluate sigma_E(), sigma_HeII, 1/nu at remapped quad. pts.
    nu_m = nu0_HeII/xm;
    fm_E = this->RadiationSpectrum(nu_m);
    fm_ni = this->CrossSections(nu_m,2);
    fm_nu = 1.0/nu_m;
    if ((fm_E == -1.0) || (fm_ni == -1.0)) 
      ENZO_FAIL("ComputeRadiationIntegrals Error in evaluating spectrum");
  
    nu_r = nu0_HeII/xr;
    fr_E = this->RadiationSpectrum(nu_r);
    fr_ni = this->CrossSections(nu_r,2);
    fr_nu = 1.0/nu_r;
    if ((fr_E == -1.0) || (fr_ni == -1.0)) 
      ENZO_FAIL("ComputeRadiationIntegrals Error in evaluating spectrum");
  
    //      compute integrals using re-scaled function values
    fl = fl_E*fl_ni/xl/xl;
    fm = fm_E*fm_ni/xm/xm;
    fr = fr_E*fr_ni/xr/xr;
    intSigESigHeII += nu0_HeII*(xr-xl)/6.0*(fl + 4.0*fm + fr);
    
    fl = fl_E*fl_ni*fl_nu/xl/xl;
    fm = fm_E*fm_ni*fm_nu/xm/xm;
    fr = fr_E*fr_ni*fr_nu/xr/xr;
    intSigESigHeIInu += nu0_HeII*(xr-xl)/6.0*(fl + 4.0*fm + fr);

    //      quit if we have finished interval
    if (xr >= (Ulimit-10.0*epsilon))  break;
  
  }
  

  if (debug) {
    printf("  Computed Radiation Integrals:\n");
    printf("    intSigE          = %22.16e\n",intSigE);
    printf("    intSigESigHI     = %22.16e\n",intSigESigHI);
    printf("    intSigESigHInu   = %22.16e\n",intSigESigHInu);
    printf("    intSigESigHeI    = %22.16e\n",intSigESigHeI);
    printf("    intSigESigHeInu  = %22.16e\n",intSigESigHeInu);
    printf("    intSigESigHeII   = %22.16e\n",intSigESigHeII);
    printf("    intSigESigHeIInu = %22.16e\n",intSigESigHeIInu);
  }

  return SUCCESS;
}

#endif
