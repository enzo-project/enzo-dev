/***********************************************************************
/
/  SETS THE MULTI-SPECIES RATES BASED ON THE EXTERNAL RADIATION FIELD
/
/  written by: Greg Bryan
/  date:       October, 1996
/  modified1:  Pascal Paschos, Robert Harkness
/  date:       1st July 2002
/  modified2:  Pascal Paschos
/  date:       August, 2002	
/  modified3:  Robert Harkness - Killed 32-bit IBM C++ bug
/  date:       15 September 2002
/  modified4:  Pascal Paschos & Robert Harkness - parameter controls
/  date:       March, 2006
/  modified5:  Robert Harkness - removed C++ I/O
/  date:       February 29th, 2008
/  modified6:  Ji-hoon Kim - changes to work with non-cosmological runs
/  date:       November, 2009
/
/  PURPOSE:    
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

#include <string.h>
#include <stdio.h>
#include <math.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "CosmologyParameters.h"

/* function prototypes */

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int RadiationFieldLymanWernerTable(float Redshift,float *J21);

int RadiationFieldCalculateRates(FLOAT Time)
{
  /* Return if there is no radiation (rates should be all zero). */

  if (RadiationFieldType == 0)
    return SUCCESS;

  /* Set units. */

/*
  FILE *fp;

  fp = fopen("uvb_data", "a");

  fprintf(fp, "AdjustUVBackground = %"ISYM"\n", AdjustUVBackground);
  fprintf(fp, "SetUVBAmplitude = %"ISYM"\n", SetUVBAmplitude);
  fprintf(fp, "SetHeIIHeatingScale = %"FSYM"\n", SetHeIIHeatingScale);
  fprintf(fp, "RadiationFieldType = %"ISYM"\n", RadiationFieldType);
*/

/*
  static ofstream output("uvb_data");
  output<<"AdjustUVBackground = "<<AdjustUVBackground<<'\n';
  output<<"SetUVBAmplitude = "<<SetUVBAmplitude<<'\n';
  output<<"SetHeIIHeatingScale = "<<SetHeIIHeatingScale<<'\n';
  output<<"RadiationFieldType = "<<RadiationFieldType<<'\n';
*/

  /*
  if (!ComovingCoordinates) {   
    ENZO_FAIL("RadiationField only defined for cosmology.\n");
  }  
  */

  FLOAT a = 1.0, dadt;
  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1, 
    VelocityUnits = 1, TimeUnits = 1, aUnits = 1;
  float Redshift;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }

  if (ComovingCoordinates) {
    CosmologyComputeExpansionFactor(Time, &a, &dadt);
    aUnits = 1.0/(1.0 + InitialRedshift);
    Redshift = 1.0/(a*aUnits) - 1;
  } else {  
    Redshift = RadiationFieldRedshift;   
    CoolData.RadiationRedshiftOn = RadiationFieldRedshift+0.2;
    CoolData.RadiationRedshiftOff = 0.0;
    CoolData.RadiationRedshiftFullOn = RadiationFieldRedshift+0.1;
    CoolData.RadiationRedshiftDropOff = 0.0;
  }

    
  double tbase1 = TimeUnits;
  double xbase1 = LengthUnits/(a*aUnits);
  double dbase1 = DensityUnits*POW(a*aUnits, 3);
  double mh     = 1.67e-24;
  double CoolingUnits = (POW(aUnits, 5) * xbase1*xbase1 * mh*mh) /
                        (POW(tbase1, 3) * dbase1);

  /* ------------------------------------------------------------------ */
  /* First, calculate the ramp value, a number between 0 and 1 which
     is used as an external control to the radiation. 
     (Only used if RadiationFieldType = 1 to 4 or equal to 12). */

  float Ramp = 0;

  if (Redshift < CoolData.RadiationRedshiftOn && 
      Redshift >= CoolData.RadiationRedshiftOff) {

    if (Redshift > CoolData.RadiationRedshiftFullOn)
      Ramp = 0.5 - 0.5*tanh(15.0*(Redshift - 0.5*
	    (CoolData.RadiationRedshiftOn+CoolData.RadiationRedshiftFullOn)));
    else if (Redshift < CoolData.RadiationRedshiftDropOff)
      Ramp = (Redshift - CoolData.RadiationRedshiftOff + CoolData.f0to3*
	                 (CoolData.RadiationRedshiftDropOff - Redshift)) /
             (CoolData.RadiationRedshiftDropOff - 
	      CoolData.RadiationRedshiftOff);
    else
      Ramp = 1.0;

  }

  float RampX = 0.0;  // this is unused below
  float XRadRedShiftOn = 7.0;
  float XRadRedShiftOff = 0.0 ;
  float XRadRedShiftFullOn = 6.0 ;
  float XRadRedShiftDropOff = 0.0 ;
  float XRadRedShiftf0to3 = 0.1 ;
  
  if (Redshift < XRadRedShiftOn &&
      Redshift >= XRadRedShiftOff) {
    
    if (Redshift > XRadRedShiftFullOn)
      RampX = 0.5 - 0.5*tanh(15.0*(Redshift - 0.5*
				   (XRadRedShiftOn + XRadRedShiftFullOn)));
    else if (Redshift < XRadRedShiftDropOff)
      RampX = (Redshift - XRadRedShiftOff + XRadRedShiftf0to3*
	       (XRadRedShiftDropOff - Redshift)) /
	(XRadRedShiftDropOff -
	 XRadRedShiftOff);
    else
      RampX = 1.0;
    
  }

  float exp_arg = -1.0 * POW(Redshift-2.3, 2);

  /* Above z = 2.3, increase the width of the Gaussian by a factor of
     3 for HI and HeI. */

  float beta2 = (AdjustUVBackgroundHighRedshift && Redshift > 2.3) ? 9.0 : 1.0;

  /* ------------------------------------------------------------------ */
  /* 1) For the Haardt and Madau (1996) quasar spectrum (alpha_q = 1.5) */

  if (RadiationFieldType == 1) {

    RateData.k24 = 6.7e-13 * POW(1.0+Redshift, 0.43) * exp(exp_arg/1.95/beta2)
                     * TimeUnits * Ramp;
    RateData.k25 = 6.3e-15 * POW(1.0+Redshift, 0.51) * exp(exp_arg/2.35) 
                     * TimeUnits * Ramp;
    RateData.k26 = 3.2e-13 * POW(1.0+Redshift, 0.50) * exp(exp_arg/2.00/beta2) 
                     * TimeUnits * Ramp;
    CoolData.piHI   = 4.7e-24 * POW(1.0+Redshift, 0.43) * exp(exp_arg/1.95/beta2) 
                     / CoolingUnits * Ramp;
    CoolData.piHeI  = 8.2e-24 * POW(1.0+Redshift, 0.50) * exp(exp_arg/2.00/beta2) 
                     / CoolingUnits * Ramp;
    CoolData.piHeII = 1.6e-25 * POW(1.0+Redshift, 0.51) * exp(exp_arg/2.35) 
                     / CoolingUnits * Ramp;
  }   
    
  /* ------------------------------------------------------------------ */
  /* 2) For the Haardt and Madau (1996) quasar spectrum (alpha_q = 1.8) */

  if (RadiationFieldType == 2) {
    RateData.k24 = 5.6e-13 * POW(1.0+Redshift, 0.43) * exp(exp_arg/1.95/beta2)
                 * TimeUnits * Ramp;
    RateData.k25 = 3.2e-15 * POW(1.0+Redshift, 0.30) * exp(exp_arg/2.60)
                 * TimeUnits * Ramp;
    RateData.k26 = 4.8e-13 * POW(1.0+Redshift, 0.43) * exp(exp_arg/1.95/beta2)
                 * TimeUnits * Ramp;
    CoolData.piHI   = 3.9e-24 * POW(1.0+Redshift, 0.43) * exp(exp_arg/1.95/beta2)
                 / CoolingUnits * Ramp;
    CoolData.piHeI  = 6.4e-24 * POW(1.0+Redshift, 0.43) * exp(exp_arg/2.10/beta2)
                 / CoolingUnits * Ramp;
    CoolData.piHeII = 8.7e-26 * POW(1.0+Redshift, 0.30) * exp(exp_arg/2.70)
                 / CoolingUnits * Ramp;
  }

  /* ------------------------------------------------------------------ */
  /* 3) This used to be a modified version of (1) but with the HeII heating rate
     multiplied by 1.8. Since the parameter is now controlled externally this option
     is not needed altogether. Instead, we are going to use it as the 
     HM-12 background reconstructed and modified to match observations 
     by Kirkman & Tytler (2005). */ 
    

  if (RadiationFieldType == 3) {
    RateData.k24 = 1.04e-12 * POW(1.0+Redshift, 0.231)
                 * exp( -0.6818 * POW(Redshift-1.855, 2.0) /
             (1.0+0.1646 * POW(Redshift+0.3097, 2.0)) )
      * TimeUnits * Ramp;
    RateData.k25 = 1.84e-14 * POW(1.0+Redshift, -1.038)
                 * exp( -1.1640 * POW(Redshift-1.973, 2.0) /
                 (1.0+0.1940 * POW(Redshift-0.6561, 2.0)) )
                 * TimeUnits * Ramp;
    RateData.k26 = 5.79e-13 * POW(1.0+Redshift, 0.278)
                 * exp( -0.8260 * POW(Redshift-1.973, 2.0) /
                 (1.0+0.1730 * POW(Redshift+0.2880, 2.0)) )
                 * TimeUnits * Ramp;
    CoolData.piHI   = 8.86e-24 * POW(1.0+Redshift, -0.0290)
                 * exp( -0.7055 * POW(Redshift-2.003, 2.0) /
                 (1.0+0.1884 * POW(Redshift+0.2888, 2.0)) )
                 / CoolingUnits * Ramp;
    CoolData.piHeI  = 5.86e-24 * POW(1.0+Redshift, 0.1764)
                 * exp( -0.8029 * POW(Redshift-2.088, 2.0) /
                 (1.0+0.1732 * POW(Redshift+0.1362, 2.0)) )
                 / CoolingUnits * Ramp;
    CoolData.piHeII = 2.17e-25 * POW(1.0+Redshift, -0.2196)
                 * exp( -1.070 * POW(Redshift-1.782, 2.0) /
                 (1.0+0.2124 * POW(Redshift-0.9213, 2.0)) )
                 / CoolingUnits * Ramp;

    float KP_mod = 1.0 ; 

    if (Redshift >= 3 ) {
      KP_mod = 1.3 ;
    }
    if (Redshift > 2 && Redshift < 3) {
      KP_mod = (1.0+(Redshift-2.0)*0.3) ;
    }

    RateData.k24 *= KP_mod ;
    RateData.k25 *= KP_mod ;
    RateData.k26 *= KP_mod ;
    CoolData.piHI *= KP_mod ; 
    CoolData.piHeI *= KP_mod ;
    CoolData.piHeII *= KP_mod ;

                                                                            
}

  /* ------------------------------------------------------------------ */
  /* 4) For the Haardt and Madau (2001) quasar spectrum (alpha_q = 1.57) 
        with X-ray Compton heating from Madau & Efstathiou (9902080). */

  if (RadiationFieldType == 4) {

    RateData.k24 = 1.04e-12 * POW(1.0+Redshift, 0.231)
                 * exp( -0.6818 * POW(Redshift-1.855, 2.0) /
             (1.0+0.1646 * POW(Redshift+0.3097, 2.0)) )
      * TimeUnits * Ramp;
    RateData.k25 = 1.84e-14 * POW(1.0+Redshift, -1.038)
                 * exp( -1.1640 * POW(Redshift-1.973, 2.0) /
                 (1.0+0.1940 * POW(Redshift-0.6561, 2.0)) )
                 * TimeUnits * Ramp;
    RateData.k26 = 5.79e-13 * POW(1.0+Redshift, 0.278)
                 * exp( -0.8260 * POW(Redshift-1.973, 2.0) /
                 (1.0+0.1730 * POW(Redshift+0.2880, 2.0)) )
                 * TimeUnits * Ramp;
    CoolData.piHI   = 8.86e-24 * POW(1.0+Redshift, -0.0290)
                 * exp( -0.7055 * POW(Redshift-2.003, 2.0) /
                 (1.0+0.1884 * POW(Redshift+0.2888, 2.0)) )
                 / CoolingUnits * Ramp;
    CoolData.piHeI  = 5.86e-24 * POW(1.0+Redshift, 0.1764)
                 * exp( -0.8029 * POW(Redshift-2.088, 2.0) /
                 (1.0+0.1732 * POW(Redshift+0.1362, 2.0)) )
                 / CoolingUnits * Ramp;
    CoolData.piHeII = 2.17e-25 * POW(1.0+Redshift, -0.2196)
                 * exp( -1.070 * POW(Redshift-1.782, 2.0) /
                 (1.0+0.2124 * POW(Redshift-0.9213, 2.0)) )
                 / CoolingUnits * Ramp;

    /* This is sigma_thompson * c * (effective <h \mu>/<m_e c^2>) *
       U_xray * 1eV.  U_xray is the energy density of XRB in , <h \mu> is the
       average photon energy in keV, corrected for relativistic effects. 
       Eq.(4) and Eq.(11) of Madau & Efstathiou (1999) */

    float RedshiftXrayCutoff = 5.0;
    
    CoolData.comp_xray = 6.65e-25 * 3.0e10 * 
                        (31.8*POW(1.0+Redshift, 0.3333)/511.0) * 
                        (6.3e-5 * 1.6e-12) * 
                        POW(1.0 + Redshift, 4) * 
                        exp(-POW(Redshift/RedshiftXrayCutoff, 2)) / 
                        CoolingUnits; 
    /*
    CoolData.comp_xray = 6.65e-25 * 3.0e10 * 
                        (1.0/511.0e3) * 
                        (4.0 * 1.38e-16/1.6e-12) *
                        (6.3e-5 * 1.6e-12) * 
                        POW(1.0 + Redshift, 4) * 
                        exp(-POW(Redshift/RedshiftXrayCutoff, 2)) / 
                        CoolingUnits; */

    /* The effective temperature (in K).  Eq.(10) of Madau & Efstathiou (1999) 
       with U_xray(z=0) = 6.3e-5 and U_cmb(z=0) = 0.256 eV/cm3  */

    CoolData.temp_xray = 31.8e3*POW(1.0+Redshift, 0.3333)*1.6e-12/
                         (4.0*1.38e-16) *
                         6.3e-5 * POW(1.0 + Redshift, 4) * 
                         exp(-POW(Redshift/RedshiftXrayCutoff, 2)) /
                         (0.256 * (1+Redshift));  

    /*
    CoolData.temp_xray = 31.8e3*POW(1.0+Redshift, 0.3333)*1.6e-12/
       (4.0*1.38e-16); */  // <-- this seems to be wrong, Ji-hoon Kim
  }


  /* ------------------------------------------------------------------ */
  /* 5) Featureless power-law spectrum (see calc_rates.src). */

  if (RadiationFieldType == 5)
    ;

  /* ------------------------------------------------------------------ */
  /* 8) An absorbed (hard) quasar-like spectrum plus molecular H constant
     photo-dissociation  (the rates are calculated in calc_rates.src). */

  if (RadiationFieldType == 8) {

    /* Insert redshift-dependent term here */

    /* molecular hydrogen constant photo-dissociation
       rate is 1.13e-8 * F_LW  (flux in Lyman-Werner bands)
       Note: this is hard-coded to F_LW = 1e-21; flux normalization controls
       the hard-radiation component (i.e. > 13.6 eV)*/

    RateData.k31 = 1.13e8 * 1.0e-21 * TimeUnits;
  }

  /* ------------------------------------------------------------------ */
  /* 9) molecular hydrogen constant photo-dissociation only! 
     rate is 1.13e-8 * F_LW  (flux in Lyman-Werner bands) */

  if (RadiationFieldType == 9) {

    /* Redshift dependent background. */
    if (TabulatedLWBackground) {

      /* Get J_21 of LW background from table. */
      float LW_J21;
      if (RadiationFieldLymanWernerTable(Redshift, &LW_J21) == FAIL) {
	ENZO_FAIL("Error in RadiationFieldLymanWernerTable.\n");
      }

      /* molecular hydrogen constant photo-dissociation
	 rate is 1.13e8 * F_LW  (flux in Lyman-Werner bands) 
	 F = 4 Pi J, so rate is 1.42e9 * J_LW. */
      RateData.k31 = 1.42e9 * LW_J21 * 1.0e-21 * TimeUnits;

    }

    /* Constant background. */
    else {
      RateData.k31 = 1.13e8 * CoolData.f3 * TimeUnits;
    }

  }

  /* ------------------------------------------------------------------ */
  /* 10 & 11) - internally-computed radiation field.  Most of the rates
     are calculated in RadiationFieldUpdate, but the comp_xray rate is
     calculated here, from the X-ray energy density and temperature. */

  if (RadiationFieldType >= 10 && RadiationFieldType <= 11) {
    CoolData.comp_xray = 8.0/3.0*0.665e6/(9.1e0*3.0e0) *
      RadiationData.ComptonXrayEnergyDensity*1.38e2*
      (double(1.0e-30)/CoolingUnits);
    CoolData.temp_xray = RadiationData.ComptonXrayTemperature;
  }

  /* ------------------------------------------------------------------ */
  /* 12) For the Haardt and Madau (2001) QSO+GAL (alpha_q = 1.57)        */

  if (RadiationFieldType == 12) {

   
    RateData.k24 = 1.04e-12 * POW(1.0+Redshift, 0.231)
                 * exp( -0.6818 * POW(Redshift-1.855, 2.0) /
             (1.0+0.1646 * POW(Redshift+0.3097, 2.0)) )
      * TimeUnits * Ramp;
    RateData.k25 = 1.84e-14 * POW(1.0+Redshift, -1.038)
                 * exp( -1.1640 * POW(Redshift-1.973, 2.0) /
                 (1.0+0.1940 * POW(Redshift-0.6561, 2.0)) )
                 * TimeUnits * Ramp;
    RateData.k26 = 5.79e-13 * POW(1.0+Redshift, 0.278)
                 * exp( -0.8260 * POW(Redshift-1.973, 2.0) /
                 (1.0+0.1730 * POW(Redshift+0.2880, 2.0)) )
                 * TimeUnits * Ramp;
    CoolData.piHI   = 8.86e-24 * POW(1.0+Redshift, -0.0290)
                 * exp( -0.7055 * POW(Redshift-2.003, 2.0) /
                 (1.0+0.1884 * POW(Redshift+0.2888, 2.0)) )
                 / CoolingUnits * Ramp;
    CoolData.piHeI  = 5.86e-24 * POW(1.0+Redshift, 0.1764)
                 * exp( -0.8029 * POW(Redshift-2.088, 2.0) /
                 (1.0+0.1732 * POW(Redshift+0.1362, 2.0)) )
                 / CoolingUnits * Ramp;
    CoolData.piHeII = 2.17e-25 * POW(1.0+Redshift, -0.2196)
                 * exp( -1.070 * POW(Redshift-1.782, 2.0) /
                 (1.0+0.2124 * POW(Redshift-0.9213, 2.0)) )
                 / CoolingUnits * Ramp;
    }

 /* 13) Fits to the (unpublished) Haardt & Madau (2005) model that
        ships with Cloudy v8.00 (Michael Kuhlen, 06/02/2010).

        Like RadiationFieldType=4, this model includes both a quasar
        and a galaxy contribution, but no X-ray Compton heating. */

  if (RadiationFieldType == 13) {

    RateData.k24 = 4.68e-12 * POW(1.0+Redshift, -0.592) * 
      exp( -0.7156 * POW(Redshift-2.292, 2.0) / 
	   (1.0 + 0.2009 * POW(Redshift+0.548, 2.0)) ) *
      TimeUnits * Ramp;

    RateData.k25 = 3.67e-14 * POW(1.0+Redshift, -0.258) * 
      exp( -1.1378 * POW(Redshift-1.875, 2.0) / 
	   (1.0 + 0.1727 * POW(Redshift-0.607, 2.0)) ) *
      TimeUnits * Ramp;

    RateData.k26 = 1.11e-11 * POW(1.0+Redshift, -1.451) * 
      exp( -0.7552 * POW(Redshift-2.954, 2.0) / 
	   (1.0 + 0.2625 * POW(Redshift+0.918, 2.0)) ) *
      TimeUnits * Ramp;

    if (MultiSpecies > 1) {

      RateData.k27 = 9.37e-10 * POW(1.0+Redshift, 2.526) *
	exp( 12.7560 * POW(Redshift+2.707, 2.0) /
	     (1.0 + -0.0301 * POW(Redshift+111.200, 2.0)) ) *
	TimeUnits * Ramp;
      
      RateData.k28 = 9.84e-12 * POW(1.0+Redshift, 2.115) * 
	exp( -1.9137 * POW(Redshift-1.748, 2.0) / 
	     (1.0 + 0.3896 * POW(Redshift+3.679, 2.0)) ) * 
	TimeUnits * Ramp;
      
      RateData.k29 = 6.87e-12 * POW(1.0+Redshift, -0.456) * 
	exp( -0.7601 * POW(Redshift-2.234, 2.0) / 
	     (1.0 + 0.1955 * POW(Redshift+0.655, 2.0)) ) * 
	TimeUnits * Ramp;
      
      RateData.k30 = 2.61e-12 * POW(1.0+Redshift, -2.198) * 
	exp( -0.5824 * POW(Redshift-3.811, 2.0) / 
	     (1.0 + 0.3241 * POW(Redshift+0.838, 2.0)) ) * 
	TimeUnits * Ramp;
      
      /* LymanSawtoothSuppressionFactor is supposed to account for the
	 suppression of LW flux due to Lyman-series absorption (giving
	 a sawtooth pattern), a la Haiman & Abel, & Rees (2000).  This
	 is really only applicable when there is plenty of neutral
	 hydrogen is around, so I'm scaling it with Ramp. */
      float LymanSawtoothSuppressionFactor = 0.1 + 0.9 * Ramp;

      RateData.k31 = 2.36e-12 * POW(1.0+Redshift, 1.185) * 
	exp( -0.2173 * POW(Redshift-3.211, 2.0) / 
	     (1.0 + 1.1852 * POW(Redshift+0.156, 2.0)) ) * 
	TimeUnits * Ramp * LymanSawtoothSuppressionFactor;

    }

    CoolData.piHI = 4.09e-23 * POW(1.0+Redshift, -0.746) * 
      exp( -0.7469 * POW(Redshift-2.419, 2.0) / 
	   (1.0 + 0.2120 * POW(Redshift+0.686, 2.0)) )
      / CoolingUnits * Ramp;

    CoolData.piHeII = 1.28e-24 * POW(1.0+Redshift, -0.504) *
      exp( -1.0742 * POW(Redshift-1.889, 2.0) /
	   (1.0 + 0.1919 * POW(Redshift-0.518, 2.0)) )
      / CoolingUnits * Ramp;

    CoolData.piHeI = 4.86e-22 * POW(1.0+Redshift, -2.302) * 
      exp( -0.5250 * POW(Redshift-3.900, 2.0) / 
	   (1.0 + 0.3452 * POW(Redshift+0.673, 2.0)) )
      / CoolingUnits * Ramp;

 }

  /* ------------------------------------------------------------------ */
  /* 14) molecular hydrogen photo-dissociation only with a fit from
     the Wise & Abel (2005) model, updated for WMAP7, reionization at
     z=6.8.  Only valid for 6<z<50.  At z>50, set to tiny.  At z<6,
     set it to CoolData.f3.

     rate is 1.13e-8 * F_LW  (flux in Lyman-Werner bands) */

  if (RadiationFieldType == 14) {
    float logJ;
    if (Redshift > 50.0)
      RateData.k31 = tiny_number;
    else if (Redshift > 6.0) {
      logJ = -23.56688 + 4.56213e-1 * (1.0+Redshift) -
	2.67982e-2 * POW(1.0+Redshift, 2.0) + 
	5.88234e-4 * POW(1.0+Redshift, 3.0) -
	5.05576e-6 * POW(1.0+Redshift, 4.0);
      RateData.k31 = 1.13e8 * POW(10.0,logJ) * TimeUnits;  //*4.0*M_PI
    } else
      RateData.k31 = 1.13e8 * CoolData.f3 * TimeUnits;
  }


  /* ------------------------------------------------------------------ */
  /* 15) tabulated photoionization and photoheating rates from 
     Haardt and Madau 2012.  rates are read into memory in 
     InitializeHM12Photorates.C and interpolated here. */

  if (RadiationFieldType == 15) {

    double zhi, delta, frac;
    int ilo, ihi;


    /* first find indices that bracket the input redshift */
    if ( Redshift <= RateData.HM12RedshiftLo ) {
      ilo = 0;
      ihi = ilo + 1;
      frac = 0.0;
    }
    else if ( Redshift >= RateData.HM12RedshiftHi ) {
      ihi = RateData.HM12NumberOfRedshiftBins - 1;
      ilo = ihi - 1;
      frac = 1.0;
    }
    else {
      ihi = 0;
      zhi = RateData.HM12Redshifts[ihi];
      while ( zhi < Redshift ) {
	ihi = ihi + 1;
	zhi = RateData.HM12Redshifts[ihi];	  
      }
      ilo = ihi-1;
      delta = RateData.HM12Redshifts[ihi] - RateData.HM12Redshifts[ilo];
      frac = ( Redshift - RateData.HM12Redshifts[ilo] ) / delta;
    }
    
    /* now interpolate the rates */
    delta = RateData.HM12GH1[ihi] - RateData.HM12GH1[ilo];
    RateData.k24 = POW( 10.0, RateData.HM12GH1[ilo] + frac * delta );
    RateData.k24 *= TimeUnits * Ramp;
    
    delta = RateData.HM12GHe2[ihi] - RateData.HM12GHe2[ilo];
    RateData.k25 = POW( 10.0, RateData.HM12GHe2[ilo] + frac * delta );
    RateData.k25 *= TimeUnits * Ramp;

    delta = RateData.HM12GHe1[ihi] - RateData.HM12GHe1[ilo];
    RateData.k26 = POW( 10.0, RateData.HM12GHe1[ilo] + frac * delta );
    RateData.k26 *= TimeUnits * Ramp;

    delta = RateData.HM12GhH1[ihi] - RateData.HM12GhH1[ilo];
    CoolData.piHI = POW( 10.0, RateData.HM12GhH1[ilo] + frac * delta );
    CoolData.piHI = CoolData.piHI / CoolingUnits * Ramp;

    delta = RateData.HM12GhHe1[ihi] - RateData.HM12GhHe1[ilo];
    CoolData.piHeI = POW( 10.0, RateData.HM12GhHe1[ilo] + frac * delta );
    CoolData.piHeI = CoolData.piHeI / CoolingUnits * Ramp;
                 
    delta = RateData.HM12GhHe2[ihi] - RateData.HM12GhHe2[ilo];
    CoolData.piHeII = POW( 10.0, RateData.HM12GhHe2[ilo] + frac * delta );
    CoolData.piHeII = CoolData.piHeII / CoolingUnits * Ramp;

    /*
    if (MyProcessorNumber == ROOT_PROCESSOR) {
      printf( "Altay: TimeUnits = %"ESYM"\n", TimeUnits);
      printf( "Altay: CoolingUnits = %"ESYM"\n", CoolingUnits);
      printf( "Altay: Redshift = %"ESYM"\n", Redshift);
      printf( "Altay: RadiationRedshiftOn = %"ESYM"\n", CoolData.RadiationRedshiftOn);
      printf( "Altay: RadiationRedshiftFullOn = %"ESYM"\n", CoolData.RadiationRedshiftFullOn);
      printf( "Altay: Ramp = %"ESYM"\n", Ramp);
      printf( "Altay: k24 = %"ESYM"\n",RateData.k24); 
      printf( "Altay: k25 = %"ESYM"\n",RateData.k25); 
      printf( "Altay: k26 = %"ESYM"\n",RateData.k26); 
      printf( "Altay: piHI = %"ESYM"\n",CoolData.piHI);
      printf( "Altay: piHeI = %"ESYM"\n",CoolData.piHeI);
      printf( "Altay: piHeII = %"ESYM"\n",CoolData.piHeII);
    }
    */
    
  }


/* ------------------------------------------------------------------ */
  if (RadiationFieldType < 0 || RadiationFieldType > 15) {
    ENZO_VFAIL("RadiationFieldType %"ISYM" not recognized.\n", 
	    RadiationFieldType)
   }

  if (AdjustUVBackground < 0 || AdjustUVBackground > 2 ) {
   ENZO_VFAIL("AdjustUVBackground Type %"ISYM" not recognized.\n",
            AdjustUVBackground)
  }

  if (AdjustUVBackground == 0 ) {
/* All rates are computed at their theoretical values */
    return SUCCESS;
  }

#ifdef NO_HEII_UVB
  if (MyProcessorNumber == ROOT_PROCESSOR)
    printf("WARNING: Not using a radiation background for HeIII\n");
  RateData.k25 = 0.0;
  CoolData.piHeII = 0.0;
#endif /* NO_HEII_UVB */

  if (RadiationFieldType == 1 || RadiationFieldType == 2 || 
      RadiationFieldType == 3 || RadiationFieldType == 4 || 
      RadiationFieldType == 12 || RadiationFieldType == 13) {
  
    if (AdjustUVBackground == 1 ) {
      /* Increases the HM-type HeII photoheating rates by the default value of 1.8 */
      SetUVBAmplitude = 1.0 ;
      SetHeIIHeatingScale = 1.8 ;
      CoolData.piHeII *= (SetUVBAmplitude * SetHeIIHeatingScale);
    }
    if (AdjustUVBackground == 2) {

      /* Alters the HM-type rates by the user supplied values.
	 If some are not supplied then it uses the default 
	 values defined in the previous if statement and set
	 in the SetDefaultGlobalValues.C  */

/*
  fprintf(fp, "In Loop with RadiationFieldType = %"ISYM"\n", RadiationFieldType);
  fprintf(fp, "In Loop with AdjustUVBackground = %"ISYM"\n", AdjustUVBackground);
  fprintf(fp, "In Loop with SetUVBAmplitude = %"ISYM"\n", SetUVBAmplitude);
  fprintf(fp, "In Loop with SetHeIIHeatingScale= %"FSYM"\n", SetHeIIHeatingScale);
*/

/*
  output<<"In Loop with RadiationFieldType ="<<RadiationFieldType<<'\n';
  output<<"In Loop with AdjustUVBackground ="<<AdjustUVBackground<<'\n';
  output<<"In Loop with SetUVBAmplitude ="<<SetUVBAmplitude<<'\n';
  output<<"In Loop with SetHeIIHeatingScale="<<SetHeIIHeatingScale<<'\n';
*/

      RateData.k24    *= SetUVBAmplitude ;
      RateData.k25    *= SetUVBAmplitude ;
      RateData.k26    *= SetUVBAmplitude ;
      CoolData.piHI   *= SetUVBAmplitude ;
      CoolData.piHeI  *= SetUVBAmplitude ;
      CoolData.piHeII *= (SetUVBAmplitude * SetHeIIHeatingScale);
    }

  }

/*
  fprintf(fp, "Redshift = %"FSYM"\n", Redshift);
  fprintf(fp, "RateData.k24 = %"FSYM"\n", RateData.k24/TimeUnits/Ramp);
  fprintf(fp, "RateData.k25 = %"FSYM"\n", RateData.k25/TimeUnits/Ramp);
  fprintf(fp, "RateData.k26 = %"FSYM"\n", RateData.k26/TimeUnits/Ramp);
  fprintf(fp, "CoolData.piHI = %"FSYM"\n", CoolData.piHI*CoolingUnits/Ramp);
  fprintf(fp, "CoolData.piHeI = %"FSYM"\n", CoolData.piHeI*CoolingUnits/Ramp);
  fprintf(fp, "CoolData.piHeII = %"FSYM"\n", CoolData.piHeII*CoolingUnits/Ramp);
*/

/*
  output<<"Redshift = "<<Redshift<<'\n';
  output<<"RateData.k24 = "<<RateData.k24/TimeUnits/Ramp<<'\n';
  output<<"RateData.k25 = "<<RateData.k25/TimeUnits/Ramp<<'\n';
  output<<"RateData.k26 = "<<RateData.k26/TimeUnits/Ramp<<'\n';
  output<<"CoolData.piHI = "<<CoolData.piHI*CoolingUnits/Ramp<<'\n';
  output<<"CoolData.piHeI = "<<CoolData.piHeI*CoolingUnits/Ramp<<'\n';
  output<<"CoolData.piHeII = "<<CoolData.piHeII*CoolingUnits/Ramp<<'\n';
  output<<"-----------------------------------------------------------"<<'\n';
*/

/*
  fclose(fp);
*/


  return SUCCESS;
}
