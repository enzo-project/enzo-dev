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

  const float eV_erg = 6.241509e11;

  int i;
  double L_UV, cgs_convert, _mass;
  float x, x2, EnergyFractionLW, MeanEnergy, XrayLuminosityFraction;
  float Mform, EnergyFractionHeI, EnergyFractionHeII;

  /* for individual star */
  float Teff, g, Z, M, tau;
  float DensityUnits, LengthUnits, TemperatureUnits, tunits, VelocityUnits;
  float H_ionizing_energy   = 13.5984; // eV
  float HeI_ionizing_energy = 24.5874; // eV
  float R;

  int *se_table_position, *rad_table_position;

  if (this->Mass < 0.1 && (ABS(this->type) != IndividualStar))  // Not "born" yet
    _mass = this->FinalMass;
  else
    _mass = this->Mass;
  x = log10((float)(_mass));
  x2 = x*x;



  switch(ABS(this->type)) {

    /* Luminosities from Schaerer (2002) */

  case PopIII:
    nbins = (PopIIIHeliumIonization &&
	     !RadiativeTransferHydrogenOnly) ? 3 : 1;
#ifdef TRANSFER    
    if (!RadiativeTransferOpticallyThinH2) nbins++;
#endif
    E[0] = 28.0;
    E[1] = 30.0;
    E[2] = 58.0;
    E[3] = 12.8;
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

    nbins = 2;

    M   = this->BirthMass;   // interpolate grids on initial ZAMS mass
    Z   = this->Metallicity;
    tau = this->LifeTime;
    tau = tau * (TimeUnits); // convert to cgs

    se_table_position  = this->ReturnSETablePosition();
    rad_table_position = this->ReturnRadTablePosition();

    /* set ionizing radiation photon rates */
    if (M >= IndividualStarIonizingRadiationMinimumMass){
      IndividualStarInterpolateProperties(Teff, R,
                                        se_table_position[0], se_table_position[1],
                                        M, Z);

      g = IndividualStarSurfaceGravity( M, R); // M in solar - R in cgs

      IndividualStarComputeIonizingRates(Q[0], Q[1],
                                         rad_table_position[0], rad_table_position[1], rad_table_position[2],
                                         Teff, g, Z);

      // compute average energy by integrating over the black body spectrum
      H_ionizing_energy   /= eV_erg; // convert to ergs
      HeI_ionizing_energy /= eV_erg; // convert to ergs
      ComputeAverageEnergy(&E[0], &H_ionizing_energy, &Teff);
      ComputeAverageEnergy(&E[1], &HeI_ionizing_energy, &Teff);

      // convert to eV from cgs
      E[0] = E[0] * eV_erg;
      E[1] = E[1] * eV_erg;
      E[2] = 58.0; // HeII

      // Functions above return the ionizing flux at stellar surface.
      // Convert to ionizing photon rate
      Q[0] = Q[0] * 4.0 * pi * R*R;
      Q[1] = Q[1] * 4.0 * pi * R*R;
      Q[2] = 0.0; // do not track
    } else{
      E[0] = 21.0; //irrelevant
      E[1] = 30.0; //irrelevant
      E[2] = 60.0; //irrelevant

      Q[0] = 0.0; Q[1] = 0.0; Q[2] = 0.0;
    } // end ionizing radiation

    /* compute optically thin rates */
    if( (IndividualStarFUVHeating || IndividualStarLWRadiation) &&
       IndividualStarOTRadiationMethod == 1){
      nbins = 5; // + LW and + FUV

      E[3] = 12.8; // LW radiation
      E[4] = 8.8;  // FUV radiation

      if(IndividualStarFUVHeating && M > IndividualStarOTRadiationMass){
          float l_fuv;
          IndividualStarComputeFUVLuminosity(l_fuv, this);
          Q[4] = l_fuv / (E[3] / eV_erg);
      } else{
          Q[4] = 0.0;
      }

      if(IndividualStarLWRadiation && M > IndividualStarOTRadiationMass){
        float l_lw;
        IndividualStarComputeLWLuminosity(l_lw, this);
        Q[3] = l_lw / (E[4] / eV_erg);
      } else{
        Q[3] = 0.0;
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
    E[3] = 12.8;
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
    E[3] = 12.8;
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

    nbins = 2;
    E[0]  = 0.0; E[1] = 0.0;
    Q[0]  = 0.0; Q[1] = 0.0;
    printf("Star_ComputePhotonRates: WARNING IndividualStarWD and IndividualStarRemnant particles should not be ionizing sources\n");

    break;


  default:
    ENZO_VFAIL("Star type = %"ISYM" not understood.\n", this->type)

  } // ENDSWITCH

  return SUCCESS;
}
