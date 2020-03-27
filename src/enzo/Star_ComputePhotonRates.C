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
  float x, x2, EnergyFractionLW, MeanEnergy, XrayLuminosityFraction;
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



  switch(ABS(this->type)) {

    /* Luminosities from Schaerer (2002) */

  case IndividualStarPopIII:
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
    _mass = max(min((float)(_mass), 500), 5);
    if (_mass > 9 && _mass <= 500) {
      Q[0] = pow(10.0, 43.61 + 4.9*x   - 0.83*x2);
      Q[1] = pow(10.0, 42.51 + 5.69*x  - 1.01*x2);
      Q[2] = pow(10.0, 26.71 + 18.14*x - 3.58*x2);
      Q[3] = pow(10.0, 44.03 + 4.59*x  - 0.77*x2);
    } else if (_mass > 5 && _mass <= 9) {
      Q[0] = pow(10.0, 39.29 + 8.55*x);
      Q[1] = pow(10.0, 29.24 + 18.49*x);
      Q[2] = pow(10.0, 26.71 + 18.14*x - 3.58*x2);
      Q[3] = pow(10.0, 44.03 + 4.59*x  - 0.77*x2);
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

      this->ComputeIonizingRates(Q[0], Q[1]);

      // compute average energy by integrating over the black body spectrum
      ComputeAverageEnergy(E[0], (HI_ionizing_energy / eV_erg), this->Teff);
      ComputeAverageEnergy(E[1], (HeI_ionizing_energy / eV_erg), this->Teff);

      // convert to eV from cgs
      E[0] = E[0] * eV_erg;
      E[1] = E[1] * eV_erg;
      E[2] = HeII_ionizing_energy; // HeII

      // Functions above return the ionizing flux at stellar surface.
      // Convert to ionizing photon rate
      Q[0] = Q[0] * 4.0 * pi * R*R;
      Q[1] = Q[1] * 4.0 * pi * R*R;
      Q[2] = 0.0; // not tracked in this model

    }

    /* compute optically thin rates */
    if( (IndividualStarFUVHeating || IndividualStarLWRadiation) &&
        (M >= IndividualStarOTRadiationMass)){

      nbins = 8 ; // + LW and + FUV

      E[3] = LW_photon_energy; // LW-band radiation (11.2 - 13.6 eV)

      // For now, no IR, Xray, or spectrum types (no need to actually zero)
      // E[4] = 0.0; E[5] = 0.0; E[6] = 0.0;
      // Q[4] = 0.0; Q[5] = 0.0; Q[6] = 0.0;

      // FUV in the 5.6 - 11.2 eV band --- 11.2 to 13.6 is in LW
      //      - set to 8.4 eV for convenience
      E[7] = FUV_photon_energy ;  // FUV radiation - average of 5.6 to 11.2 eV range

      if(IndividualStarLWRadiation){
        float l_lw;
        this->ComputeLWLuminosity(l_lw);
        Q[3] = l_lw / (E[3] / eV_erg);
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
