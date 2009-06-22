/********************************************************************************
/  Solve cooling using Multispecies chemistry/cooling and Cloudy tables.
/
/  written by: Britton Smith
/  date:       January, 2008
/  modified1:  
/
/  PURPOSE:  Calculate de/dt for a grid and return it to higher level solver.
/
/  RETURNS:
/    SUCCESS or FAIL
/
/*******************************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

#define MH 1.67e-24
#define DEFAULT_MU 0.6
#define H_MASS_FRACTION 0.76

#define MU_TOLERANCE 1e-4

float coolingGridInterpolate1D(float temperature,float *dataField);

float coolingGridInterpolate2D(float parameter1,float temperature,float *dataField);

float coolingGridInterpolate3D(float parameter1,float parameter2,float temperature,
			       float *dataField);

float coolingGridInterpolate4D(float parameter1,float parameter2,float parameter3,
			       float temperature,float *dataField);

float coolingGridInterpolate5D(float parameter1,float parameter2,float parameter3,
			       float parameter4,float temperature,float *dataField);

int SolveCloudyCooling(float *density,float *totalenergy,float *gasenergy,
		       float *velocity1,float *velocity2,float *velocity3,
		       float *De,float *HI,float *HII,
		       float *HeI,float *HeII,float *HeIII,
		       float *HM,float *H2I,float *H2II,
		       float *DI,float *DII,float *HDI,
		       float *metalDensity,
		       int *GridDimension,int GridRank,float dtFixed,
		       float afloat,float TemperatureUnits,float LengthUnits,
		       float aUnits,float DensityUnits,float TimeUnits,
		       int RadiationShield,float HIShieldFactor,
		       float HeIShieldFactor,float HeIIShieldFactor,
		       bool *iterationMask,float *edot)
{

  /* Compute size (in floats) of the current grid. */

  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  // Declare some variables.

  int i;
  int mu_iter;
  float mu,energy,h_mass_fraction,electron_fraction,cloudy_electron_fraction;
  float gamma2,nH2,nOther,x,factor;
  float logHNumberDensity,logMetallicity,cellTemperature,logElectronFraction;
  float mu_metalfree,mu_old;
  float cloudyCooling,cloudyHeating,edotMetals,cloudyCoolingAtTCMB;

  // Variables for interpolation over enzo cooling data

  float logtem0 = log10(CoolData.TemperatureStart);
  float logtem9 = log10(CoolData.TemperatureEnd);
  float dlogtem = (logtem9 - logtem0)/(CoolData.NumberOfTemperatureBins - 1);
  float logCellTemperature,logTIndex;
  int multispeciesCoolIndex;

  // H/He atomic cooling variables

  float ceHI,ceHeI,ceHeII,ciHI,ciHeI,ciHeIS,ciHeII,
    reHII,reHeII1,reHeII2,reHeIII,brem;

  // H2 cooling variables

  float kH_H2,kH2_H2,vibrationH,vibrationL,rotationH,rotationL,Qn;

  // H2 cooling variables for Galli & Palla 1998

  float GPLow,GPHigh,GPHigh1;

  // H2 cooling variables for Simons's fix

  float gaHI, gaH2, gaHe, gaHp, gael, gphdl, galdl, gphdl1;

  // HD cooling variables

  float hdlte,hdlow,hdlte1,hdlow1;

  // Compton cooling variables

  float compton1,compton2;

  /* Get conversion units (ripped from InitializeEquilibriumCoolData.C) */

  double dom    = DensityUnits*pow(afloat,3)/MH;
  double tbase1 = TimeUnits;
  double xbase1 = LengthUnits/(afloat*aUnits);
  double dbase1 = DensityUnits * POW(afloat*aUnits, 3);
  double CoolUnit = (POW(aUnits,5) * POW(xbase1,2) * POW(MH,2)) /
                    (POW(tbase1,3) * dbase1);
  double CurrentRedshift = 1.0/(afloat*aUnits) - 1.0;

  // Set compton cooling

  if (ComovingCoordinates && CloudyCoolingData.CMBTemperatureFloor) {
    compton1 = CoolData.comp * pow((1 + CurrentRedshift),4);
    compton2 = 2.73 * (1 + CurrentRedshift);

    // If constant temperature floor is higher than T_CMB, use that.
    // Emulate a redshift independent CMB floor.

    if (CloudyCoolingData.ConstantTemperatureFloor > compton2) {
      compton2 = CloudyCoolingData.ConstantTemperatureFloor;
      compton1 = CoolData.comp * pow((compton2 / 2.73),4);
    }

  }
  
  // If CMBTemperatureFloor is off, but ConstantTemperatureFloor is used.
  else if (CloudyCoolingData.ConstantTemperatureFloor) {
    compton2 = CloudyCoolingData.ConstantTemperatureFloor;
    compton1 = CoolData.comp * pow((compton2 / 2.73),4);
  }

  else {
    compton1 = tiny_number;
    compton2 = tiny_number;
  }

  // Loop over all grid cells.

  for (i = 0;i < size;i++) {

    if (iterationMask[i]) {

      // Initialize edot.

      edot[i] = 0;

      // calculate H mass fraction

      if (MultiSpecies == 0) {
	h_mass_fraction = H_MASS_FRACTION;
      }

      else{

	h_mass_fraction = HI[i] + HII[i];

	if (MultiSpecies > 1) {
	  h_mass_fraction += HM[i] + H2I[i] + H2II[i];
	}
	if (MultiSpecies > 2) {
	  h_mass_fraction += (HDI[i] / 3.0);
	}

	h_mass_fraction /= density[i];

      }

      // Get metallicity.

      if (CloudyCoolingData.CloudyCoolingGridRank > 2) {
	
	// Calculate metallicity from metal density.
	// Actually, store it as log(metallicity), since that's what we interpolate over.

	logMetallicity = log10(metalDensity[i] / (h_mass_fraction * density[i] * CloudyCoolingData.CloudyMetallicityNormalization));

      }

      // Get electron fraction.

      if (CloudyCoolingData.CloudyCoolingGridRank > 3) {

	if (MultiSpecies == 0) {
	  electron_fraction = 1.0;
	}

	else {
	  electron_fraction = 2 * De[i] / (density[i] * (1 + h_mass_fraction));
	}

	logElectronFraction = log10(electron_fraction);

	// Get extra electrons contributed by metals.

	cloudy_electron_fraction = electron_fraction * 
	  (1 + (2 * CloudyCoolingData.CloudyElectronFractionFactor * pow(10,logMetallicity) * h_mass_fraction) / (1 + h_mass_fraction));

      }

      // calculate mu

      if (MultiSpecies == 0) {
	mu = DEFAULT_MU;
      }

      else {

	mu = De[i] + HI[i] + HII[i] + (HeI[i] + HeII[i] + HeIII[i])/4.0;

	if (MultiSpecies > 1) {
	  mu += HM[i] + (H2I[i] + H2II[i])/2.0;
	}
	if (MultiSpecies > 2) {
	  mu += (DI[i] + DII[i])/2.0 + (HDI[i]/3.0);
	}

	mu = density[i] / mu;

      }

      // Density has been converted into proper by the function calling this, 
      // so we need to multiply by a^3 because DensityUnits converts from 
      // comoving [code units] to proper [cgs].
      logHNumberDensity = log10(density[i]*h_mass_fraction*DensityUnits*pow(afloat,3)/MH);

      // calculate cell temperature

      // Zeus - total energy is really gas energy
      if (HydroMethod == 2) {
	energy = totalenergy[i];
      }

      // PPM
      else {
	// with Dual Energy Formalism
	if (DualEnergyFormalism) {
	  energy = gasenergy[i];
	}
	// without Dual Energy Formalism
	else {
	  energy = totalenergy[i] - 0.5*velocity1[i]*velocity1[i];
	  if (GridRank > 1)
	    energy -= 0.5*velocity2[i]*velocity2[i];
	  if (GridRank > 2)
	    energy -= 0.5*velocity3[i]*velocity3[i];
	}
      }

      cellTemperature = (Gamma - 1.0) * energy * TemperatureUnits * mu;

      // Correct temperature for gamma from H2 (ripped from cool1d_multi.src)

      if (MultiSpecies > 1) {

	nH2 = (H2I[i] + H2II[i])/2.0;
	nOther = De[i] + HI[i] + HII[i] + (HeI[i] + HeII[i] + HeIII[i])/4.0; // why no H-?

	if (nH2/nOther > 1e-3) {
	  x = 6100 / cellTemperature;
	  if (x > 10) {
	    gamma2 = 2.5;
	  }
	  else {
	    gamma2 = 0.5*(5.0 + 2.0*pow(x,2)*exp(x)/pow(exp(x)-1,2));
	  }
	}
	else {
	  gamma2 = 2.5;
	}
	gamma2 = 1.0 + (nH2 + nOther)/(nH2*gamma2 + nOther/(Gamma-1.0));
	cellTemperature *= (gamma2 - 1.0)/(Gamma - 1.0);

      }

#define DONT_ADD_MU_FROM_METALS
#ifdef ADD_MU_FROM_METALS

      // Include additional mu from the metals.

      if (CloudyCoolingData.IncludeCloudyMMW > 0) {

	mu_metalfree = mu;
	mu_old = 0;
	mu_iter = 0;

	while ((fabs(mu-mu_old)/mu > MU_TOLERANCE) && (mu_iter < 10)) {
	  mu_old = mu;

	  switch(CloudyCoolingData.CloudyCoolingGridRank) {

	  // Interpolate over temperature.
	  case 1:
	    mu = mu_metalfree +
	      coolingGridInterpolate1D(cellTemperature,CloudyCoolingData.CloudyMeanMolecularWeight);
	    break;

	  // Interpolate over density and temperature.
	  case 2:
	    mu = mu_metalfree +
	      coolingGridInterpolate2D(logHNumberDensity,cellTemperature,
				       CloudyCoolingData.CloudyMeanMolecularWeight);
	    break;

	  // Interpolate over density, metallicity, and temperature.
	  case 3:
	    mu = mu_metalfree +
	      coolingGridInterpolate3D(logHNumberDensity,logMetallicity,cellTemperature,
				       CloudyCoolingData.CloudyMeanMolecularWeight);
	    break;

	  // Interpolate over density, metallicity, electron fraction, and temperature.
	  case 4:
	    mu = mu_metalfree +
	      coolingGridInterpolate4D(logHNumberDensity,logMetallicity,logElectronFraction,
				       cellTemperature,
				       CloudyCoolingData.CloudyMeanMolecularWeight);
	    break;

	  // Interpolate over density, metallicity, electron fraction, redshift, and temperature.
	  case 5:
	    mu = mu_metalfree +
	      coolingGridInterpolate5D(logHNumberDensity,logMetallicity,logElectronFraction,
				       CurrentRedshift,cellTemperature,
				       CloudyCoolingData.CloudyMeanMolecularWeight);
	    break;

	  default:
	    fprintf(stderr,"Cloudy cooling grid rank must be less than or equal to %"ISYM".\n",
		    CLOUDY_COOLING_MAX_DIMENSION);
	    return FAIL;

	  }

	  cellTemperature *= mu/mu_old;

	  mu_iter++;

	} // while ((fabs(mu-mu_old)/mu > MU_TOLERANCE) && (mu_iter < 10))

      } // if (CloudyCoolingData.IncludeCloudyMMW > 0)

#endif // ADD_MU_FROM_METALS

      // If cloudy is supplying only metals, get all H/He cooling

      if (MultiSpecies > 0) {

	// Get Index for multispecies interpolation

	logCellTemperature = log10(cellTemperature);

	multispeciesCoolIndex = int((logCellTemperature - logtem0)/dlogtem);
	multispeciesCoolIndex = max(multispeciesCoolIndex,0);
	multispeciesCoolIndex = min(multispeciesCoolIndex,CoolData.NumberOfTemperatureBins - 2);
	logTIndex = multispeciesCoolIndex*dlogtem + logtem0;

	ceHI = CoolData.ceHI[multispeciesCoolIndex] + (logCellTemperature - logTIndex)
	  *(CoolData.ceHI[multispeciesCoolIndex+1] - CoolData.ceHI[multispeciesCoolIndex])/dlogtem;
	ceHeI = CoolData.ceHeI[multispeciesCoolIndex] + (logCellTemperature - logTIndex)
	  *(CoolData.ceHeI[multispeciesCoolIndex+1] - CoolData.ceHeI[multispeciesCoolIndex])/dlogtem;
	ceHeII = CoolData.ceHeII[multispeciesCoolIndex] + (logCellTemperature - logTIndex)
	  *(CoolData.ceHeII[multispeciesCoolIndex+1] - CoolData.ceHeII[multispeciesCoolIndex])/dlogtem;
	ciHI = CoolData.ciHI[multispeciesCoolIndex] + (logCellTemperature - logTIndex)
	  *(CoolData.ciHI[multispeciesCoolIndex+1] - CoolData.ciHI[multispeciesCoolIndex])/dlogtem;
	ciHeI = CoolData.ciHeI[multispeciesCoolIndex] + (logCellTemperature - logTIndex)
	  *(CoolData.ciHeI[multispeciesCoolIndex+1] - CoolData.ciHeI[multispeciesCoolIndex])/dlogtem;
	ciHeIS = CoolData.ciHeIS[multispeciesCoolIndex] + (logCellTemperature - logTIndex)
	  *(CoolData.ciHeIS[multispeciesCoolIndex+1] - CoolData.ciHeIS[multispeciesCoolIndex])/dlogtem;
	ciHeII = CoolData.ciHeII[multispeciesCoolIndex] + (logCellTemperature - logTIndex)
	  *(CoolData.ciHeII[multispeciesCoolIndex+1] - CoolData.ciHeII[multispeciesCoolIndex])/dlogtem;
	reHII = CoolData.reHII[multispeciesCoolIndex] + (logCellTemperature - logTIndex)
	  *(CoolData.reHII[multispeciesCoolIndex+1] - CoolData.reHII[multispeciesCoolIndex])/dlogtem;
	reHeII1= CoolData.reHeII1[multispeciesCoolIndex] + (logCellTemperature - logTIndex)
	  *(CoolData.reHeII1[multispeciesCoolIndex+1]- CoolData.reHeII1[multispeciesCoolIndex])/dlogtem;
	reHeII2= CoolData.reHeII2[multispeciesCoolIndex] + (logCellTemperature - logTIndex)
	  *(CoolData.reHeII2[multispeciesCoolIndex+1]- CoolData.reHeII2[multispeciesCoolIndex])/dlogtem;
	reHeIII= CoolData.reHeIII[multispeciesCoolIndex] + (logCellTemperature - logTIndex)
	  *(CoolData.reHeIII[multispeciesCoolIndex+1]- CoolData.reHeIII[multispeciesCoolIndex])/dlogtem;
	brem = CoolData.brem[multispeciesCoolIndex] + (logCellTemperature - logTIndex)
	  *(CoolData.brem[multispeciesCoolIndex+1] - CoolData.brem[multispeciesCoolIndex])/dlogtem;

	edot[i] += (

		    // Collisional excitations

		    - ceHI  *HI  [i]*De[i]             // ce of HI
		    - ceHeI *HeII[i]*De[i]*De[i]*dom/4.0  // ce of HeI
		    - ceHeII*HeII[i]*De[i]/4.0         // ce of HeII
		   
		    // Collisional ionizations

		    - ciHI  *HI  [i]*De[i]             // ci of HI
		    - ciHeI *HeI [i]*De[i]/4.0         // ci of HeI
		    - ciHeII*HeII[i]*De[i]/4.0         // ci of HeII
		    - ciHeIS*HeII[i]*De[i]*De[i]*dom/4.0  // ci of HeIS

		    // Recombinations

		    - reHII  *HII  [i]*De[i]          // re of HII
		    - reHeII1*HeII [i]*De[i]/4.0      // re of HeII
		    - reHeII2*HeII [i]*De[i]/4.0      // re of HeII
		    - reHeIII*HeIII[i]*De[i]/4.0      // re of HeIII

		    // Compton cooling

		    - compton1*(cellTemperature-compton2)*De[i]/dom

		    // X-ray compton heating
		   
		    - CoolData.comp_xray * (cellTemperature-CoolData.temp_xray)*De[i]/dom
		   
		    // Bremsstrahlung

		    - brem*(HII[i]+HeII[i]/4.0+HeIII[i])*De[i]

		    ) / density[i];

      }

      // H2 cooling

      if (MultiSpecies > 1) {

	// Simon's fix to the H_2 cooling.

#define USE_GLOVER_ABEL2008
#ifdef USE_GLOVER_ABEL2008
	gaHI = CoolData.GAHI[multispeciesCoolIndex] + (logCellTemperature - logTIndex)
	  *(CoolData.GAHI[multispeciesCoolIndex+1] - CoolData.GAHI[multispeciesCoolIndex])/dlogtem;
	gaH2 = CoolData.GAH2[multispeciesCoolIndex] + (logCellTemperature - logTIndex)
	  *(CoolData.GAH2[multispeciesCoolIndex+1] - CoolData.GAH2[multispeciesCoolIndex])/dlogtem;
	gaHe = CoolData.GAHe[multispeciesCoolIndex] + (logCellTemperature - logTIndex)
	  *(CoolData.GAHe[multispeciesCoolIndex+1] - CoolData.GAHe[multispeciesCoolIndex])/dlogtem;
	gaHp = CoolData.GAHp[multispeciesCoolIndex] + (logCellTemperature - logTIndex)
	  *(CoolData.GAHp[multispeciesCoolIndex+1] - CoolData.GAHp[multispeciesCoolIndex])/dlogtem;
	gael = CoolData.GAel[multispeciesCoolIndex] + (logCellTemperature - logTIndex)
	  *(CoolData.GAel[multispeciesCoolIndex+1] - CoolData.GAHI[multispeciesCoolIndex])/dlogtem;
	gphdl = CoolData.GP99HighDensityLimit[multispeciesCoolIndex] + (logCellTemperature - logTIndex)
	  *(CoolData.GP99HighDensityLimit[multispeciesCoolIndex+1] - 
	    CoolData.GP99HighDensityLimit[multispeciesCoolIndex])/dlogtem;

	galdl = gaHI * HI[i]  + gaH2 * H2I[i]
	  + gaHe * HeI[i] + gaHp * HII[i]
	  + gael * De[i];
	gphdl1 = gphdl/dom;
	edot[i] -= CoolData.ih2co * H2I[i] *
	  gphdl / (1.0 + gphdl1/galdl) / (2.0*dom*density[i]);

#else

#define USE_GALLI_PALLA1999
#ifdef USE_GALLI_PALLA1999

	// Get H2 cooling from Galli & Palla 1998

	GPLow = (logCellTemperature - logTIndex)*
	  (CoolData.GP99LowDensityLimit[multispeciesCoolIndex+1]-
	   CoolData.GP99LowDensityLimit[multispeciesCoolIndex])/dlogtem
	  + CoolData.GP99LowDensityLimit[multispeciesCoolIndex];

	GPHigh = (logCellTemperature - logTIndex)*
	  (CoolData.GP99HighDensityLimit[multispeciesCoolIndex+1]-
	   CoolData.GP99HighDensityLimit[multispeciesCoolIndex])/dlogtem
	  + CoolData.GP99HighDensityLimit[multispeciesCoolIndex];

	GPHigh1 = GPHigh/(HI[i]*dom);

	edot[i] -= (H2I[i]/2.0/dom/density[i]) * GPHigh/(1 + GPHigh1/GPLow);

#else /* USE_GALLI_PALLA1999 */

        // Get H2 cooling from Lepp & Shull 1983

	kH_H2 = (logCellTemperature - logTIndex)*
	  (CoolData.hyd01k[multispeciesCoolIndex+1]-CoolData.hyd01k[multispeciesCoolIndex])/dlogtem
	  + CoolData.hyd01k[multispeciesCoolIndex];
	kH2_H2 = (logCellTemperature - logTIndex)*
	  (CoolData.h2k01[multispeciesCoolIndex+1]-CoolData.h2k01[multispeciesCoolIndex])/dlogtem
	  + CoolData.h2k01[multispeciesCoolIndex];
	vibrationH = (logCellTemperature - logTIndex)*
	  (CoolData.vibh[multispeciesCoolIndex+1]-CoolData.vibh[multispeciesCoolIndex])/dlogtem
	  + CoolData.vibh[multispeciesCoolIndex];
	rotationH = (logCellTemperature - logTIndex)*
	  (CoolData.roth[multispeciesCoolIndex+1]-CoolData.roth[multispeciesCoolIndex])/dlogtem
	  + CoolData.roth[multispeciesCoolIndex];
	rotationL = (logCellTemperature - logTIndex)*
	  (CoolData.rotl[multispeciesCoolIndex+1]-CoolData.rotl[multispeciesCoolIndex])/dlogtem
	  + CoolData.rotl[multispeciesCoolIndex];

	Qn = 1.2*pow((HI[i]*dom),0.77) + pow((H2I[i]*dom/2.0),0.77);
	rotationL *= Qn;
	vibrationL = (HI[i]*kH_H2 + H2I[i]/2.0*kH2_H2)*dom*8.18e-13;

	edot[i] -= (H2I[i]/2.0/dom/density[i]) * (rotationH/(1 + rotationH/rotationL) + 
						  vibrationH/(1 + vibrationH/vibrationL));

#endif /* USE_GALLI_PALLA1999 */
#endif /* USE_GLOVER_ABEL2008 */

      } // if (MultiSpecies > 1)

      // HD Cooling

      if (MultiSpecies > 2) {

	if (cellTemperature > compton2) {
	  hdlte = CoolData.HDlte[multispeciesCoolIndex] + (logCellTemperature - logTIndex)
	    *(CoolData.HDlte[multispeciesCoolIndex+1] - CoolData.HDlte[multispeciesCoolIndex])/dlogtem;
	  hdlow = CoolData.HDlow[multispeciesCoolIndex] + (logCellTemperature - logTIndex)
	    *(CoolData.HDlow[multispeciesCoolIndex+1] - CoolData.HDlow[multispeciesCoolIndex])/dlogtem;
	}
	else {
	  hdlte = 0;
	  hdlow = 0;
	}

	hdlte1 = hdlte/(HI[i]*dom);
	hdlow1 = max(hdlow, tiny_number);
	edot[i] -= HDI[i] * (hdlte/(1.0 + hdlte1/hdlow1)) / (3.0*dom*density[i]);

      }

      // Photo-ionization heating of H and He
      // taken from cool1d_multi.src.

#define USE_PHOTOIONIZATION
#ifdef USE_PHOTOIONIZATION
      if (RadiationShield == 0) {

	// regular version

	if (RadiationFieldType == 8) {

	  // 1) heating assuming high energy photons produces secondary
	  // electrons which do the heating (Shull & Steenberg, 1985).

	  x = max(HII[i]/(HI[i]+HII[i]), 1.0e-4);
	  factor = 0.9971 * (1.0 - pow((1.0 - pow(x,0.2663)),1.3163));
	  edot[i] += CoolData.ipiht * factor * (
						CoolData.piHI  *HI [i]      // pi of HI
						+ CoolData.piHeI *HeI [i]*0.25 // pi of HeI
						+ CoolData.piHeII*HeII[i]*0.25 // pi of HeII
						) / dom / density[i];

	}

	else {

	  // 2) standard heating

	  edot[i] += CoolData.ipiht * (
				       CoolData.piHI  *HI  [i]      // pi of HI
				       + CoolData.piHeI *HeI [i]*0.25 // pi of HeI
				       + CoolData.piHeII*HeII[i]*0.25 // pi of HeII
				       ) / dom / density[i];

	}

      }

      else {

	// version with approximate self-shielding

	edot[i] += CoolData.ipiht * (
				     (CoolData.piHI  *HI  [i] * 
				      exp(-HIShieldFactor*HI[i]*dom))
				     + (CoolData.piHeI *HeI [i] * 0.25 * 
					exp(-HeIShieldFactor*0.25*HeI[i]*dom))
				     + (CoolData.piHeII*HeII[i] * 0.25 * 
					exp(-HeIIShieldFactor*0.25*HeII[i]*dom))
				     ) / dom / density[i];

      }

#endif /* USE_PHOTOIONIZATION */

      // Metal Cooling

#define USE_CLOUDY_COOLING
#ifdef USE_CLOUDY_COOLING

      switch(CloudyCoolingData.CloudyCoolingGridRank) {

      // Interpolate over temperature.
      case 1:

	// Get Cloudy cooling

	cloudyCooling = coolingGridInterpolate1D(cellTemperature,
						 CloudyCoolingData.CloudyCooling);

	edotMetals = -pow(10,cloudyCooling);

	// If CMBTemperatureFloor is on calculate cooling at T_CMB and subtract from cooling

	if (CloudyCoolingData.CMBTemperatureFloor || CloudyCoolingData.ConstantTemperatureFloor) {
	  cloudyCoolingAtTCMB = coolingGridInterpolate1D(compton2,
							 CloudyCoolingData.CloudyCooling);

	  edotMetals += pow(10,cloudyCoolingAtTCMB);
	}

	// include heating only if requested

	if (CloudyCoolingData.IncludeCloudyHeating) {
	  cloudyHeating = coolingGridInterpolate1D(cellTemperature,
						   CloudyCoolingData.CloudyHeating);

	  edotMetals += pow(10,cloudyHeating);
	}

	break;

      // Interpolate over density and temperature.
      case 2:

	// Get Cloudy cooling

	cloudyCooling = coolingGridInterpolate2D(logHNumberDensity,cellTemperature,
						 CloudyCoolingData.CloudyCooling);

	edotMetals = -pow(10,cloudyCooling);

	// If CMBTemperatureFloor is on calculate cooling at T_CMB and subtract from cooling

	if (CloudyCoolingData.CMBTemperatureFloor || CloudyCoolingData.ConstantTemperatureFloor) {
	  cloudyCoolingAtTCMB = coolingGridInterpolate2D(logHNumberDensity,compton2,
							 CloudyCoolingData.CloudyCooling);

	  edotMetals += pow(10,cloudyCoolingAtTCMB);
	}

	// include heating only if requested

	if (CloudyCoolingData.IncludeCloudyHeating) {
	  cloudyHeating = coolingGridInterpolate2D(logHNumberDensity,cellTemperature,
						   CloudyCoolingData.CloudyHeating);

	  edotMetals += pow(10,cloudyHeating);
	}

	break;

      // Interpolate over density, metallicity, and temperature.

      case 3:

	// Get Cloudy cooling

	cloudyCooling = coolingGridInterpolate3D(logHNumberDensity,logMetallicity,cellTemperature,
						 CloudyCoolingData.CloudyCooling);

	edotMetals = -pow(10,cloudyCooling);
 
	// If CMBTemperatureFloor is on calculate cooling at T_CMB and subtract from cooling

	if (CloudyCoolingData.CMBTemperatureFloor || CloudyCoolingData.ConstantTemperatureFloor) {
	  cloudyCoolingAtTCMB = coolingGridInterpolate3D(logHNumberDensity,logMetallicity,compton2,
							 CloudyCoolingData.CloudyCooling);

	  edotMetals += pow(10,cloudyCoolingAtTCMB);
	}

	// include heating only if requested

	if (CloudyCoolingData.IncludeCloudyHeating) {
	  cloudyHeating = coolingGridInterpolate3D(logHNumberDensity,logMetallicity,cellTemperature,
						   CloudyCoolingData.CloudyHeating);

	  edotMetals += pow(10,cloudyHeating);
	}

      break;

      // Interpolate over density, metallicity, electron fraction, and temperature.

      case 4:

	// Get Cloudy cooling

	cloudyCooling = coolingGridInterpolate4D(logHNumberDensity,logMetallicity,logElectronFraction,
						 cellTemperature,
						 CloudyCoolingData.CloudyCooling);

	edotMetals = -pow(10,cloudyCooling);
 
	// If CMBTemperatureFloor is on calculate cooling at T_CMB and subtract from cooling

	if (CloudyCoolingData.CMBTemperatureFloor || CloudyCoolingData.ConstantTemperatureFloor) {
	  cloudyCoolingAtTCMB = coolingGridInterpolate4D(logHNumberDensity,logMetallicity,logElectronFraction,
							 compton2,
							 CloudyCoolingData.CloudyCooling);

	  edotMetals += pow(10,cloudyCoolingAtTCMB);
	}

	// include heating only if requested

	if (CloudyCoolingData.IncludeCloudyHeating) {
	  cloudyHeating = coolingGridInterpolate4D(logHNumberDensity,logMetallicity,logElectronFraction,
						   cellTemperature,
						   CloudyCoolingData.CloudyHeating);

	  edotMetals += pow(10,cloudyHeating);
	}

      break;

      // Interpolate over density, metallicity, electron fraction, redshift, and temperature.

      case 5:

	// Get Cloudy cooling

	cloudyCooling = coolingGridInterpolate5D(logHNumberDensity,logMetallicity,logElectronFraction,
						 CurrentRedshift,cellTemperature,
						 CloudyCoolingData.CloudyCooling);

	edotMetals = -pow(10,cloudyCooling);
 
	// If CMBTemperatureFloor is on calculate cooling at T_CMB and subtract from cooling

	if (CloudyCoolingData.CMBTemperatureFloor || CloudyCoolingData.ConstantTemperatureFloor) {
	  cloudyCoolingAtTCMB = coolingGridInterpolate5D(logHNumberDensity,logMetallicity,logElectronFraction,
							 CurrentRedshift,compton2,
							 CloudyCoolingData.CloudyCooling);

	  edotMetals += pow(10,cloudyCoolingAtTCMB);
	}

	// include heating only if requested

	if (CloudyCoolingData.IncludeCloudyHeating) {
	  cloudyHeating = coolingGridInterpolate5D(logHNumberDensity,logMetallicity,logElectronFraction,
						   CurrentRedshift,cellTemperature,
						   CloudyCoolingData.CloudyHeating);

	  edotMetals += pow(10,cloudyHeating);
	}

      break;

    default:
	fprintf(stderr,"Cloudy cooling grid rank must be less than or equal to %"ISYM".\n",
		    CLOUDY_COOLING_MAX_DIMENSION);
	return FAIL;

      } // switch(CloudyCoolingData.CloudyCoolingGridRank)

      edotMetals *= density[i] * h_mass_fraction;

      // Multiply by electron fraction if using 4d or 5d cooling.

      if (CloudyCoolingData.CloudyCoolingGridRank > 3) {

	edotMetals *= cloudy_electron_fraction;

      }

      edot[i] += edotMetals;

#endif /* USE_CLOUDY_COOLING */

      // Absolute cooling floor at Tmin of data.
      if ((cellTemperature <= CoolData.TemperatureStart) && (edot[i] < 0)) {
	edot[i] = tiny_number * 1e-20;
      }

      if (fabs(edot[i]) < tiny_number) {
	edot[i] = tiny_number;
      }

    } // if (iterationMask[i])

  } // for (i = 0;i < size;i++)

  return SUCCESS;

}
