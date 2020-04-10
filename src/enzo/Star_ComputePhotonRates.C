/***********************************************************************
/
/  COMPUTE STELLAR PHOTON EMISSION RATES
/
/  written by: John Wise
/  date:       November, 2005
/  modified1:
/
/  ---------- SPECIES --------
/  0 : HI
/  1 : HeI
/  2 : HeII
/  3 : Lyman-Werner (H2)
/
************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "phys_constants.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"

#include "IndividualStarProperties.h"

float ReturnValuesFromSpectrumTable(float ColumnDensity, float dColumnDensity, int mode);

int Star::ComputePhotonRates(const float TimeUnits, int &nbins, float E[], double Q[])
{


  int i;
  double L_UV, cgs_convert, _mass;
  float x, x2, x3, EnergyFractionLW, MeanEnergy, XrayLuminosityFraction;
  float Mform, EnergyFractionHeI, EnergyFractionHeII;

  /* for individual star */
  float Z, M, tau;
  float DensityUnits, LengthUnits, TemperatureUnits, tunits, VelocityUnits;
  float R;

  if ( (ABS(this->type) == IndividualStar) || (ABS(this->type) == IndividualStarPopIII)){
    _mass = this->BirthMass;
  } else if (this->Mass < 0.1){  // Not "born" yet
    _mass = this->FinalMass;
  } else {
    _mass = this->Mass;
  }
  x = log10((float)(_mass));
  x2 = x*x;
  x3 = x*x*x;



  switch(ABS(this->type)) {

    /* Luminosities from Schaerer (2002) */

  case IndividualStarPopIII:
  {
    if (PopIIIRadiationModel == 0){
      nbins = (PopIIIHeliumIonization &&
             !RadiativeTransferHydrogenOnly) ? 3 : 1;
#ifdef TRANSFER
      if (!RadiativeTransferOpticallyThinH2) nbins++;
#endif


      E[0] = 28.0;
      E[1] = 30.0;
      E[2] = 58.0;
      E[3] = LW_photon_energy;
      E[4] = IR_photon_energy;
      _mass = max(min((float)(_mass), 500), 5);
      if (_mass > 9 && _mass <= 500) {
        Q[0] = pow(10.0, 43.61 + 4.9*x   - 0.83*x2);
        Q[1] = pow(10.0, 42.51 + 5.69*x  - 1.01*x2);
        Q[2] = pow(10.0, 26.71 + 18.14*x - 3.58*x2);
        Q[3] = pow(10.0, 44.03 + 4.59*x  - 0.77*x2);
        Q[4] = 0.0;
      } else if (_mass > 5 && _mass <= 9) {
        Q[0] = pow(10.0, 39.29 + 8.55*x);
        Q[1] = pow(10.0, 29.24 + 18.49*x);
        Q[2] = pow(10.0, 26.71 + 18.14*x - 3.58*x2);
        Q[3] = pow(10.0, 44.03 + 4.59*x  - 0.77*x2);
        Q[4] = 0.0;
      } // ENDELSE
      else {
        for (i = 0; i < nbins; i++) Q[i] = 0.0;
      }

    } else if (PopIIIRadiationModel == 1){

      if( (IndividualStarFUVHeating || IndividualStarLWRadiation) &&
          (M >= IndividualStarOTRadiationMass)){

        nbins = 8 ; // + LW and + FUV
      }

      /* Fits to Heger and Woosley 2010

          https://github.com/aemerick/galaxy_analysis/blob/master/physics_data/Heger_Woosley_2010/Heger_woosley_radiation_fits.ipynb

          3 degree polynomial fits over 2 mass range leads to
          typical errors in the 0.1% range and maximal error
          around 5% . This fit produces large errors in the IR photon
          rate because model is not monotonic here at high mass. But fit captures
          general trend and approximately similar photon count integrated
          over mass range */

      // Stars > the maximum mass will be treated at 100 and scaled
      // later. Stars below (10Msun) will just be extrapolated using fit
      float fit_mass = min((float)(_mass), 100.0);
      const float mass_cut = 35.0; // separation in polynomials

      /* These fits return the TOTAL photons. Need to divide by lifetime */

      if ( fit_mass < mass_cut ){
        Q[0] = pow(10.0,  -1.2618e+00*x3 + 4.8613e+00*x2 -4.4976e+00*x + 6.3507e+01);
        Q[1] = pow(10.0,  -1.0090e-01*x3 -2.2907e-01*x2 + 3.3247e+00*x + 5.8788e+01);
        Q[2] = pow(10.0,  3.7096e+00*x3 -1.6701e+01*x2 + 2.8367e+01*x + 4.3733e+01);
        Q[3] = pow(10.0,  -1.8868e+00*x3 + 7.6565e+00*x2 -8.9860e+00*x + 6.5457e+01);
        Q[4] = pow(10.0,  -2.7242e+00*x3 + 1.1538e+01*x2 -1.5151e+01*x + 6.8400e+01);
        Q[7] = pow(10.0,  -2.3146e+00*x3 + 9.6461e+00*x2 -1.2168e+01*x + 6.7769e+01);
      } else if (fit_mass <= 100.0) {
        Q[0] = pow(10.0,  1.6675e+00*x3 -9.1524e+00*x2 + 1.7875e+01*x + 5.1583e+01);
        Q[1] = pow(10.0,  4.1553e-01*x3 -2.8421e+00*x2 + 7.5444e+00*x + 5.6596e+01);
        Q[2] = pow(10.0,  3.1015e-01*x3 -3.1192e+00*x2 + 1.0218e+01*x + 5.1874e+01);
        Q[3] = pow(10.0,  4.7327e+00*x3 -2.5093e+01*x2 + 4.5244e+01*x + 3.5422e+01);
        Q[4] = pow(10.0,  -3.3190e+01*x3 + 1.7386e+02*x2 -2.9988e+02*x + 2.3327e+02);
        Q[7] = pow(10.0,  -1.3126e+00*x3 + 6.2667e+00*x2 -8.4472e+00*x + 6.6389e+01);
      } else {
        for (i=0;i<nbins;i++) Q[i]=0.0;
      }
      E[0] = pow(10.0,  7.1694e-04*x3 -5.2935e-02*x2 + 2.5684e-01*x + 1.0847e+00);
      E[1] = pow(10.0,  -6.5004e-03*x3 + 5.8389e-03*x2 + 9.5175e-02*x + 1.3787e+00);
      E[2] = pow(10.0,  -9.2050e-04*x3 -5.7004e-03*x2 + 5.5792e-02*x + 1.7217e+00);
      E[3] = LW_photon_energy;
      E[4] = IR_photon_energy;
      E[7] = 0.8*FUV_photon_energy;

      if (_mass > 100.0){
        // scale massive stars assuming fixed photon rate
        // per unit mass (not a bad assumption for massive stars)
        for (i = 0; i < nbins; i++) Q[i] *= _mass / 100.0;
      }

      // convert from total photons to rate
      float inv_lifetime = 1.0 / (this->LifeTime*TimeUnits);
      for(i=0;i<nbins;i++) Q[i] *= inv_lifetime; // now photons / second
    }

  }
  break;

  case PopIII:

      nbins = (PopIIIHeliumIonization &&
             !RadiativeTransferHydrogenOnly) ? 3 : 1;
#ifdef TRANSFER
      if (!RadiativeTransferOpticallyThinH2) nbins++;
#endif


      E[0] = 28.0;
      E[1] = 30.0;
      E[2] = 58.0;
      E[3] = LW_photon_energy;
      E[4] = IR_photon_energy;
      _mass = max(min((float)(_mass), 500), 5);
      if (_mass > 9 && _mass <= 500) {
        Q[0] = pow(10.0, 43.61 + 4.9*x   - 0.83*x2);
        Q[1] = pow(10.0, 42.51 + 5.69*x  - 1.01*x2);
        Q[2] = pow(10.0, 26.71 + 18.14*x - 3.58*x2);
        Q[3] = pow(10.0, 44.03 + 4.59*x  - 0.77*x2);
        Q[4] = 0.0;
      } else if (_mass > 5 && _mass <= 9) {
        Q[0] = pow(10.0, 39.29 + 8.55*x);
        Q[1] = pow(10.0, 29.24 + 18.49*x);
        Q[2] = pow(10.0, 26.71 + 18.14*x - 3.58*x2);
        Q[3] = pow(10.0, 44.03 + 4.59*x  - 0.77*x2);
        Q[4] = 0.0;
      } // ENDELSE
      else {
        for (i = 0; i < nbins; i++) Q[i] = 0.0;
      }

    break;

  case IndividualStar:
    /* Compute HI and HeI ionizing photon rates using individual star model */

    nbins = 3;

    M   = this->BirthMass;   // interpolate grids on initial ZAMS mass
    Z   = this->Metallicity;
    tau = this->LifeTime;
    R   = this->ReturnRadius();
    tau = tau * (TimeUnits); // convert to cgs


    /* set ionizing radiation photon rates */
    if (M >= IndividualStarIonizingRadiationMinimumMass){

      this->ComputeIonizingRates(Q[0], Q[1], Q[2]);

      // compute average energy by integrating over the black body spectrum
      ComputeAverageEnergy(E[0],   (HI_ionizing_energy / eV_erg), this->Teff);
      ComputeAverageEnergy(E[1],  (HeI_ionizing_energy / eV_erg), this->Teff);
      ComputeAverageEnergy(E[2], (HeII_ionizing_energy / eV_erg), this->Teff);

      // convert to eV from cgs
      E[0] = E[0] * eV_erg; // HI
      E[1] = E[1] * eV_erg; // HeI
      E[2] = E[2] * eV_erg; // HeII

      // Functions above return the ionizing flux at stellar surface.
      // Convert to ionizing photon rate
      Q[0] = Q[0] * 4.0 * pi * R*R;
      Q[1] = Q[1] * 4.0 * pi * R*R;
      Q[2] = Q[2] * 4.0 * pi * R*R;

    }

    /* compute optically thin rates */
    if( (IndividualStarFUVHeating || IndividualStarLWRadiation) &&
        (M >= IndividualStarOTRadiationMass)){

      nbins = 8 ; // + LW and + FUV

      E[3] = LW_photon_energy; // LW-band radiation (11.2 - 13.6 eV)

      // For now, no IR, Xray, or spectrum types (no need to actually zero)
      //  E[5] = 0.0; E[6] = 0.0;
      //  Q[5] = 0.0; Q[6] = 0.0;

      E[4] = IR_photon_energy; // IR-band radiation (2.0 eV)

      // FUV in the 5.6 - 11.2 eV band --- 11.2 to 13.6 is in LW
      //      - set to 8.4 eV for convenience
      E[7] = FUV_photon_energy ;  // FUV radiation - average of 5.6 to 11.2 eV range

      if(IndividualStarLWRadiation){
        float l_lw;
        this->ComputeLWLuminosity(l_lw);
        Q[3] = l_lw / (E[3] / eV_erg);
      }

      if(IndividualStarIRRadiation){
        float l_ir;
        this->ComputeIRLuminosity(l_ir);
        Q[4] = l_ir / (E[4] / eV_erg);
      }

      if(IndividualStarFUVHeating){
        float l_fuv;
        this->ComputeFUVLuminosity(l_fuv);
        Q[7] = l_fuv / (E[7] / eV_erg); // photon rate
      }

    } // end optically thin radiation

    break;

    /* Average energy from Schaerer (2003) */

  case PopII:
    nbins = (StarClusterHeliumIonization &&
	     !RadiativeTransferHydrogenOnly) ? 3 : 1;
#ifdef TRANSFER
    if (!RadiativeTransferOpticallyThinH2 &&
	MultiSpecies > 1) nbins++;
#endif
    EnergyFractionLW   = 1.288;
    EnergyFractionHeI  = 0.2951;
    EnergyFractionHeII = 2.818e-4;
    E[0] = 21.62; // eV (good for a standard, low-Z IMF)
    E[1] = 30.0;
    E[2] = 60.0;
    E[3] = LW_photon_energy;
    Q[0] = StarClusterIonizingLuminosity * _mass;
    if (StarClusterHeliumIonization) {
      Q[1] = EnergyFractionHeI * Q[0];
      Q[2] = EnergyFractionHeII * Q[0];
      Q[0] *= 1.0 - EnergyFractionHeI - EnergyFractionHeII;
    } else {
      Q[1] = 0.0;
      Q[2] = 0.0;
    }
    Q[3] = EnergyFractionLW * Q[0];
    break;

    /* Approximation to the multi-color disk and power law of an
       accreting BH (Kuhlen & Madau 2004; Alvarez et al. 2009) */

  case BlackHole:
    nbins = 1;
    XrayLuminosityFraction = 0.43;
    EnergyFractionLW = 1.51e-3;
    MeanEnergy = 93.0;  // eV
    E[0] = 460.0;
    E[1] = 0.0;
    E[2] = 0.0;
    E[3] = LW_photon_energy;
    Q[0] = 1.12e66 * PopIIIBHLuminosityEfficiency * XrayLuminosityFraction *
      this->last_accretion_rate / E[0];
//    Below is wrong!
//    Q[0] = 3.54e58 * PopIIIBHLuminosityEfficiency * XrayLuminosityFraction *
//      this->DeltaMass / E[0];
    Q[1] = 0.0;
    Q[2] = 0.0;
    Q[3] = EnergyFractionLW * (E[0]/MeanEnergy) * Q[0];
    break;

    /* Average quasar SED by Sazonov et al.(2004), where associated
       spectral temperature is 2 keV, for accreting massive BH */

  case MBH:
    nbins = 1;
    XrayLuminosityFraction = 1.0;
    E[0] = 2000.0; //2keV
    E[1] = 0.0;
    E[2] = 0.0;
    E[3] = 0.0;
    // 1.99e33g/Ms * (3e10cm/s)^2 * 6.24e11eV/ergs = 1.12e66 eV/Ms
    Q[0] = 1.12e66 * MBHFeedbackRadiativeEfficiency * XrayLuminosityFraction *
      this->last_accretion_rate / E[0];
    Q[1] = 0.0;
    Q[2] = 0.0;
    Q[3] = 0.0;

#define NOT_HII_REGION_TEST
#ifdef HII_REGION_TEST
    Q[0] = 1.0e45 * MBHFeedbackRadiativeEfficiency * XrayLuminosityFraction / E[0];
#endif

//    fprintf(stdout, "star::ComputePhotonRates: this->last_accretion_rate = %g, Q[0]=%g\n",
//    	    this->last_accretion_rate, Q[0]);

#ifdef TRANSFER

    if (RadiativeTransferTraceSpectrum == TRUE) {
      nbins = 1;
      E[0] = ReturnValuesFromSpectrumTable(0.0, 0.0, 3); //##### mean energy if column density=0
      E[1] = 0.0;
      E[2] = 0.0;
      E[3] = 0.0;

      Q[0] = 1.12e66 * MBHFeedbackRadiativeEfficiency *
	this->last_accretion_rate / E[0];
      Q[1] = 0.0;
      Q[2] = 0.0;
      Q[3] = 0.0;

      //better check the initial mean energy when tracing spectrum
      if (MyProcessorNumber == ROOT_PROCESSOR)
	fprintf(stdout, "star::CPP: check initial mean E of photon SED: E[0] = %g\n", E[0]);
    }

#endif

    break;

  case SimpleSource:
    nbins = 1;
    // radiating particle that ramps with time, independant of mass
    E[0] = 20.0;
    Q[0] = SimpleQ; // ramping done in StarParticleRadTransfer.C
    break;

  case NormalStar:
    nbins = 1;
    E[0] = 21.0;  // Good for [Z/H] > -1.3  (Schaerer 2003)
    // Calculate Delta(M_SF) for Cen & Ostriker star particles
#ifdef TRANSFER
    Mform = this->CalculateMassLoss(dtPhoton) / StarMassEjectionFraction;
    // units of Msun/(time in code units)
    L_UV = 4 * pi * StarEnergyToStellarUV * Mform * clight * clight / dtPhoton;
    cgs_convert = SolarMass / TimeUnits;
    Q[0] = cgs_convert * L_UV * eV_erg / E[0]; // ph/s
#else
    Q[0] = 0.0;
#endif
    break;


  case IndividualStarWD:
  case IndividualStarRemnant:
  case IndividualStarUnresolved:

    nbins = 3;
    E[0]  = 0.0; E[1] = 0.0; E[2] = 0.0; E[3] = 0.0;
    Q[0]  = 0.0; Q[1] = 0.0; Q[2] = 0.0; Q[3] = 0.0;
    printf("Star_ComputePhotonRates: WARNING IndividualStarWD and IndividualStarRemnant particles should not be ionizing sources\n");

    break;


  default:
    ENZO_VFAIL("Star type = %"ISYM" not understood.\n", this->type)

  } // ENDSWITCH

  return SUCCESS;
}
