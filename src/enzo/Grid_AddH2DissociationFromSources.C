#define MYPROC MyProcessorNumber == ProcessorNumber
#define DEBUG 0

#define JEANS_LENGTH 1
/***********************************************************************
/
/  ADD H2 DISSOCIATION EMISSION FROM SHINING PARTICLES
/
/  written by: John Wise
/  date:       March, 2006
/  modified1: John Regan - added support for IR radiation and active particles
/  date:      May 2018
/
/  PURPOSE:
/
************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "CosmologyParameters.h"
#include "Star.h"
#include "phys_constants.h"

#define THRESHOLD_DENSITY_DB36 1e14
#define THRESHOLD_DENSITY_DB37 5e14
#define sigma_unit  2.08854e-4  //convert Msolar pc^-2 -> cm^-2
#define THRESHOLD_DENSITY_CO_CO   1e13
#define THRESHOLD_DENSITY_OH_OH   1e15
#define THRESHOLD_DENSITY_H2O_H2O 1e16
#define THRESHOLD_DENSITY_CO_H2   1e20
#define THRESHOLD_DENSITY_OH_H2   1e19
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	       float *VelocityUnits, FLOAT Time);
static double CalculateH2IICrossSection(float Energy);
static double CalculateIRCrossSection(float Energy);
static double JeansLength(float T, float dens, float density_units);
int grid::AddH2DissociationFromSources(Star *AllStars)
{

  Star *cstar;
  FLOAT DomainWidth[MAX_DIMENSION];
  FLOAT *ddr2[MAX_DIMENSION];
  FLOAT innerFront, outerFront, innerFront2, outerFront2;
  double Luminosity[MAX_ENERGY_BINS];
  float energies[MAX_ENERGY_BINS], kdiss_r2;
  int ipart, dim, i, j, k, index, indixe, nbins;
  int ActiveDims[MAX_DIMENSION];
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, B1Num, B2Num, B3Num;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
  int HeHIINum, DMNum   , HDIINum
    , CINum   , CIINum  , CONum     , CO2Num   , OINum   , OHNum
    , H2ONum  , O2Num   , SiINum    , SiOINum  , SiO2INum
    , CHNum   , CH2Num  , COIINum   , OIINum   , OHIINum , H2OIINum, H3OIINum, O2IINum
    , MgNum   , AlNum   , SNum      , FeNum
    , SiMNum  , FeMNum  , Mg2SiO4Num, MgSiO3Num, Fe3O4Num
    , ACNum   , SiO2DNum, MgONum    , FeSNum   , Al2O3Num
    , DustNum ;

  double l_char = 0.0, N_H2 = 0.0;
  double H2ISigma = 3.71e-18;

  double N_HDI = 0.0, N_CO = 0.0, N_OH = 0.0, N_H2O = 0.0;
  double HDISigma = 4.22276e-18;
  double COSigma  = 1.74554e-17;
  double OHSigma  = 3.56213e-18;
  double H2OSigma = 1.01407e-17;

  double N_CO_CO   = 9.48921e+14 ;//  +/- 3.906e+13    (4.116%)
  double a_CO_CO   = -0.822589   ;//  +/- 0.003422     (0.416%)
  double b_CO_CO   = -4.31591e-07;//  +/- 1.743e-08    (4.039%)
  double N_CO_H2   = 1.01533e+21 ;//  +/- 9.38e+19     (9.239%)
  double a_CO_H2   = -1.44968    ;//  +/- 0.06207      (4.282%)
  double b_CO_H2   = -0.171213   ;//  +/- 0.01413      (8.252%)
  
  double N_OH_OH   = 1.13614e+17 ;//  +/- 9.017e+15    (7.936%)
  double a_OH_OH   = -1.08514    ;//  +/- 0.01284      (1.183%)
  double b_OH_OH   = -0.00162352 ;//  +/- 0.0001285    (7.914%)
  double N_OH_H2   = 7.18749e+21 ;//  +/- 1.457e+21    (20.27%)
  double a_OH_H2   = -4.22208    ;//  +/- 0.6456       (15.29%)
  double b_OH_H2   = -0.901778   ;//  +/- 0.1155       (12.81%)
  
  double N_H2O_H2O = 1.02787e+17 ;//  +/- 1.365e+16    (13.28%)
  double a_H2O_H2O = -1.37727    ;//  +/- 0.02189      (1.589%)
  /* shielding by H2 is negligible for H2O */

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  this->DebugCheck((char*) "Grid_AddH2Dissociation");


  /* Find fields: density, total energy, velocity1-3. */

  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
        ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }

  /* Find Multi-species fields. */

  if (this->IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, 
				  HeIIINum, HMNum, H2INum, H2IINum, DINum, 
				  DIINum, HDINum) == FAIL) {
    ENZO_FAIL("Error in grid->IdentifySpeciesFields.\n");
  }

  if (this->IdentifySpeciesFieldsMD( HeHIINum, DMNum   , HDIINum
                                   , CINum   , CIINum  , CONum     , CO2Num   , OINum   , OHNum
                                   , H2ONum  , O2Num   , SiINum    , SiOINum  , SiO2INum
                                   , CHNum   , CH2Num  , COIINum   , OIINum   , OHIINum , H2OIINum,  H3OIINum,  O2IINum
                                   , MgNum   , AlNum   , SNum      , FeNum
                                   , SiMNum  , FeMNum  , Mg2SiO4Num, MgSiO3Num, Fe3O4Num
                                   , ACNum   , SiO2DNum, MgONum    , FeSNum   , Al2O3Num
                                   , DustNum ) == FAIL) {
    ENZO_FAIL("Error in grid->IdentifySpeciesFieldsMD.\n");
  }

  /* Get photo-ionization fields */

  int kphHINum, kphHeINum, kphHeIINum, kdissH2INum, kphHMNum, kdissH2IINum;
  int gammaNum;
  IdentifyRadiativeTransferFields(kphHINum, gammaNum, kphHeINum, kphHeIINum, 
				  kdissH2INum, kphHMNum, kdissH2IINum);

  int kdissHDINum, kphCINum, kphOINum, kdissCONum, kdissOHNum, kdissH2ONum;
  IdentifyRadiativeTransferFieldsMD(kdissHDINum, kphCINum, kphOINum, kdissCONum, kdissOHNum, kdissH2ONum);

  /* For now, initialize H2 photo-dissociation field. */

  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  // Now done in Grid_InitializeRadiativeTransferFields.C
//  for (i = 0; i < size; i++)
//    BaryonField[kdissH2INum][i] = 0;
  
  if (AllStars == NULL && ProblemType != 50 && EnabledActiveParticlesCount == 0)
    return SUCCESS;

  /* If using cosmology, get units. */

  float TemperatureUnits, DensityUnits, LengthUnits, VelocityUnits, 
    TimeUnits, aUnits = 1;

  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, PhotonTime);

  // Absorb the unit conversions into the cross-section
  H2ISigma *= (double)TimeUnits / ((double)LengthUnits * (double)LengthUnits);
  HDISigma *= (double)TimeUnits / ((double)LengthUnits * (double)LengthUnits);
  COSigma  *= (double)TimeUnits / ((double)LengthUnits * (double)LengthUnits);
  OHSigma  *= (double)TimeUnits / ((double)LengthUnits * (double)LengthUnits);
  H2OSigma *= (double)TimeUnits / ((double)LengthUnits * (double)LengthUnits);

  for (dim = 0; dim < GridRank; dim++) {
    DomainWidth[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];
    ActiveDims[dim] = GridEndIndex[dim] - GridStartIndex[dim] + 1;
    ddr2[dim] = new FLOAT[ActiveDims[dim]];
  }

  // Dilution factor (prevent breaking in the rate solver near the
  // star).  The value of 1 AU was chosen to be within the radius
  // where the gas is fully molecular.  If the solver still breaks,
  // then this parameter can be safely increased to ~100 AU or turning
  // on H2 self-shielding with RadiationShield = 2.
  float dilutionRadius = 4.848e-6 * pc_cm / (double) LengthUnits;  // 1 AU
  float dilRadius2 = dilutionRadius * dilutionRadius;
  float LightTravelDist = TimeUnits * clight / LengthUnits;

  // Convert from #/s to RT units
  double LConv = (double) TimeUnits / pow(LengthUnits,3);

  /* Loop over radiation sources or star particles in the grid */

  if (ProblemType == 50) {

    RadiationSourceEntry *RS;
    for (RS = GlobalRadiationSources->NextSource; RS; RS = RS->NextSource) {

      if (PhotonTime < RS->CreationTime && 
	  PhotonTime > RS->CreationTime + RS->LifeTime)
	continue;
      for(int ebin = 0; ebin < RS->EnergyBins; ebin++) {
	float IRSED = 0.0, LWSED = 0.0, H2IISED = 0.0;
	double LWLuminosity = 0.0, IRLuminosity = 0.0, H2IILuminosity = 0.0;
	double HMSigma = 0.0, H2IISigma = 0.0;
	if(RS->Energy[ebin] < 11.2) {
	  H2IISED = RS->SED[ebin];
	  IRSED = RS->SED[ebin];
	}
	else if(RS->Energy[ebin] < 13.6) {
	  H2IISED = RS->SED[ebin];
	  LWSED = RS->SED[ebin];
	}
	HMSigma = CalculateIRCrossSection(RS->Energy[ebin]);
	H2IISigma = CalculateH2IICrossSection(RS->Energy[ebin]);
	// Absorb the unit conversions into the cross-section
	H2IISigma *= (double)TimeUnits / ((double)LengthUnits * (double)LengthUnits);
	HMSigma*= (double)TimeUnits / ((double)LengthUnits * (double)LengthUnits);
	IRLuminosity = IRSED * RS->Luminosity / LConv;
	LWLuminosity = LWSED * RS->Luminosity / LConv;
	H2IILuminosity = IRLuminosity + LWLuminosity;
	if(DEBUG) {
	  printf("%s: E = %f eV\t IRSED = %f\t IRLuminosity = %g\n", __FUNCTION__, RS->Energy[ebin],
		 IRSED, IRLuminosity);
	  printf("%s: E = %f eV\t LWSED = %f\t LWLuminosity = %g\n", __FUNCTION__, RS->Energy[ebin],
		 LWSED, LWLuminosity);
	}
	/* Pre-calculate distances from cells to source */

	for (dim = 0; dim < GridRank; dim++)
	  for (i = 0, index = GridStartIndex[dim]; i < ActiveDims[dim]; 
	       i++, index++) {
	    
	    // Calculate dr_i first, then square it
	    ddr2[dim][i] = 
	      fabs(CellLeftEdge[dim][index] + 0.5*CellWidth[dim][0] - 
		   RS->Position[dim]);
	    ddr2[dim][i] = min(ddr2[dim][i], DomainWidth[dim]-ddr2[dim][i]);
	    ddr2[dim][i] = ddr2[dim][i] * ddr2[dim][i];
	  }

      /* Loop over cells */

	double radius2, radius2_yz;
	double colden = 0.0, shield = 1.0, b = 0.0, b5 = 0.0, XN = 0.0;
	double H2mass = mh*2.0, alpha = 1.1;
	double kph_hm = 0.0, kdiss_H2II = 0.0;
	int TemperatureField = 0;

        double kdiss_HDI = 0.0, kdiss_CO = 0.0, kdiss_OH = 0.0, kdiss_H2O = 0.0;
        double shieldHDI = 1.0, shieldCO = 1.0, shieldOH = 1.0, shieldH2O = 1.0;

	/* Pre-compute some quantities to speed things up */
	if (RadiativeTransferH2ShieldType > 0) {
		TemperatureField = this->GetTemperatureFieldNumberForH2Shield();
	}
	kdiss_r2 = (float) (LWLuminosity * H2ISigma / (4.0 * pi));
	kph_hm = (float) (IRLuminosity * HMSigma / (4.0 * pi));
	kdiss_H2II = (float) (H2IILuminosity * H2IISigma / (4.0 * pi));
        if(MultiSpecies > 2)
	  kdiss_HDI = (float) (LWLuminosity * HDISigma / (4.0 * pi));
        if(MetalChemistry) {
	  kdiss_CO  = (float) (LWLuminosity * COSigma  / (4.0 * pi));
	  kdiss_OH  = (float) (LWLuminosity * OHSigma  / (4.0 * pi));
	  kdiss_H2O = (float) (LWLuminosity * H2OSigma / (4.0 * pi));
        }
	for (k = 0; k < ActiveDims[2]; k++) {
	  for (j = 0; j < ActiveDims[1]; j++) {
	    radius2_yz = ddr2[1][j] + ddr2[2][k];
	    index = GRIDINDEX(0, j, k);
	    for (i = 0; i < ActiveDims[0]; i++, index++) {
	      radius2 = radius2_yz + ddr2[0][i];
	      if (radius2 < dilRadius2) {
		BaryonField[kdissH2INum][index] += kdiss_r2 / dilRadius2;
		BaryonField[kdissH2IINum][index] += kdiss_H2II / dilRadius2;
		BaryonField[kphHMNum][index] += kph_hm / dilRadius2;
                if(MultiSpecies > 2)
	          BaryonField[kdissHDINum][index] += kdiss_HDI / dilRadius2;
                if(MetalChemistry) {
	          BaryonField[kdissCONum ][index] += kdiss_CO  / dilRadius2;
	          BaryonField[kdissOHNum ][index] += kdiss_OH  / dilRadius2;
	          BaryonField[kdissH2ONum][index] += kdiss_H2O / dilRadius2;
                }
	      }
	      else {
		BaryonField[kdissH2IINum][index] += kdiss_H2II / radius2;
		BaryonField[kphHMNum][index] += kph_hm / radius2;
		BaryonField[kdissH2INum][index] += kdiss_r2 / radius2;
                if(MultiSpecies > 2)
	          BaryonField[kdissHDINum][index] += kdiss_HDI / radius2;
                if(MetalChemistry) {
	          BaryonField[kdissCONum ][index] += kdiss_CO  / radius2;
	          BaryonField[kdissOHNum ][index] += kdiss_OH  / radius2;
	          BaryonField[kdissH2ONum][index] += kdiss_H2O / radius2;
                }
	      }
	      /* Include Shielding */

	      if (RadiativeTransferH2ShieldType > 0) {
#if JEANS_LENGTH
	        l_char = JeansLength(BaryonField[TemperatureField][index],
		  BaryonField[DensNum][index], DensityUnits)*RadiativeTransferOpticallyThinH2CharLength; //cm
#else
	        l_char = CellWidth[0][0]*LengthUnits/2.0; //cm
#endif
	        N_H2 = BaryonField[H2INum][index]*l_char*DensityUnits/(2.0*mh);
                if(MultiSpecies > 2)
	          N_H2 += BaryonField[HDINum][index]*l_char*DensityUnits/(3.0*mh);
#if DEBUG
	        if(l_char/(CellWidth[0][0]*LengthUnits) != 1.0) {
		  printf("l_char/cellwidth = %g\n", l_char/(CellWidth[0][0]*LengthUnits));
		  printf("N_H2 = %g (cm^-2)\n", N_H2);
		  printf("H2I Density[%d] = %g (code)\n", index, BaryonField[H2INum][index]);
		  printf("Local H2I Density = %g (cgs)\n",
		  BaryonField[H2INum][index]*DensityUnits*CellWidth[0][0]*LengthUnits/mh);
		  printf("Temp = %f\n", BaryonField[TemperatureField][index]);
	        }
#endif
	        shield = 1.0;
                if(MultiSpecies > 2)
	          shieldHDI = 1.0;
	        if(N_H2 >= 0.01 * THRESHOLD_DENSITY_DB37) {
	          b = sqrt(2.0*kboltz*BaryonField[TemperatureField][index]/H2mass);
	          b5 = b/1e5;
	          XN = N_H2/THRESHOLD_DENSITY_DB37;
	          shield = 0.965/pow(1+XN/b5, alpha) + (0.035/sqrt(1+XN))*exp(-8.5e-4*sqrt(1+XN));
                  if(MultiSpecies > 2)
	            shieldHDI = 0.965/pow(1+XN/(b5 * 0.8165), alpha) + (0.035/sqrt(1+XN))*exp(-8.5e-4*sqrt(1+XN));
                         // multiply by sqrt(2/3)
#if DEBUG		
		  printf("%s: shield = %g\t kdiss = %g\t kdiss*shield = %g\n", __FUNCTION__, 
		        shield, BaryonField[kdissH2INum][index], BaryonField[kdissH2INum][index]*shield);
#endif
	        }
	        BaryonField[kdissH2INum][index] *= shield;
	        BaryonField[kdissHDINum][index] *= shieldHDI;
	      } // ENDIF: H2 shielding

              if(MetalChemistry) {
	        N_CO  = BaryonField[CONum ][index]*l_char*DensityUnits/(28.0*mh);
	        N_OH  = BaryonField[OHNum ][index]*l_char*DensityUnits/(17.0*mh);
	        N_H2O = BaryonField[H2ONum][index]*l_char*DensityUnits/(18.0*mh);

	        shieldCO  = 1.0;
	        shieldOH  = 1.0;
	        shieldH2O = 1.0;

                if (N_CO >= THRESHOLD_DENSITY_CO_CO) {
                  XN = N_CO / N_CO_CO;
                  shieldCO = pow(1.0 + XN, a_CO_CO) * exp(b_CO_CO * XN);
                }
                if(MultiSpecies > 1) {
                  if (N_H2 >= THRESHOLD_DENSITY_CO_H2) {
                    XN = N_H2 / N_CO_H2;
                    shieldCO *= pow(1.0 + XN, a_CO_H2) * exp(b_CO_H2 * XN);
                  }
                }

                if (N_OH >= THRESHOLD_DENSITY_OH_OH) {
                  XN = N_OH / N_OH_OH;
                  shieldOH = pow(1.0 + XN, a_OH_OH) * exp(b_OH_OH * XN);
                }
                if(MultiSpecies > 1) {
                  if (N_H2 >= THRESHOLD_DENSITY_OH_H2) {
                    XN = N_H2 / N_OH_H2;
                    shieldOH *= pow(1.0 + XN, a_OH_H2) * exp(b_OH_H2 * XN);
                  }
                }

                if (N_H2O >= THRESHOLD_DENSITY_H2O_H2O) {
                  XN = N_H2O / N_H2O_H2O;
                  shieldH2O = pow(1.0 + XN, a_H2O_H2O);
                }

	        BaryonField[kdissCONum ][index] *= shieldCO;
	        BaryonField[kdissOHNum ][index] *= shieldOH;
	        BaryonField[kdissH2ONum][index] *= shieldH2O;
              }
	    } // END: i-direction
	  } // END: j-direction
	} // END: k-direction
      }
    } // ENDFOR sources

  } // ENDIF ProblemType == 50

  else if (AllStars != NULL) {
    double LWLuminosity = 0.0;
    for (cstar = AllStars; cstar; cstar = cstar->NextStar) {

      // Skip if not 'living'
      if (!(cstar->FeedbackFlag == NO_FEEDBACK ||
	    cstar->FeedbackFlag == CONT_SUPERNOVA)) 
	continue;
      
      /* Determine H2 emission rate */

      if (cstar->ComputePhotonRates(TimeUnits, nbins, energies, Luminosity) == FAIL) {
	ENZO_FAIL("Error in ComputePhotonRates.\n");
      }
      LWLuminosity = Luminosity[3];

      /* Pre-calculate distances from cells to source */

      for (dim = 0; dim < GridRank; dim++)
	for (i = 0, index = GridStartIndex[dim]; i < ActiveDims[dim]; 
	     i++, index++) {
	  
	  // Calculate dr_i first, then square it
	  ddr2[dim][i] = 
	    fabs(CellLeftEdge[dim][index] + 0.5*CellWidth[dim][0] -
		 cstar->pos[dim]);
	  ddr2[dim][i] = min(ddr2[dim][i], DomainWidth[dim]-ddr2[dim][i]);
	  ddr2[dim][i] = ddr2[dim][i] * ddr2[dim][i];
	}

      /* Loop over cells */

      double radius2, radius2_yz;

      kdiss_r2 = (float) (LWLuminosity * H2ISigma / (4.0 * pi));

      for (k = 0; k < ActiveDims[2]; k++) {
	for (j = 0; j < ActiveDims[1]; j++) {
	  radius2_yz = ddr2[1][j] + ddr2[2][k];
	  index = GRIDINDEX(0, j, k);
	  for (i = 0; i < ActiveDims[0]; i++, index++) {
	    radius2 = radius2_yz + ddr2[0][i];
	    if (radius2 < dilRadius2)
	      BaryonField[kdissH2INum][index] += kdiss_r2 / dilRadius2;
	    else
	      BaryonField[kdissH2INum][index] += kdiss_r2 / radius2;
	  } // END: i-direction
	} // END: j-direction
      } // END: k-direction
    } // ENDFOR stars

  } // ENDELSE AllStars
  /* The APs should be hooked into the Global Radiation Sources. */
  else if (EnabledActiveParticlesCount > 0) { 
    RadiationSourceEntry *RS;
    for (RS = GlobalRadiationSources->NextSource; RS; RS = RS->NextSource) {
      
      if (PhotonTime < RS->CreationTime && 
	  PhotonTime > RS->CreationTime + RS->LifeTime)
	continue;
      for(int ebin = 0; ebin < RS->EnergyBins; ebin++) {
	float IRSED = 0.0, LWSED = 0.0, H2IISED = 0.0;
	double LWLuminosity = 0.0, IRLuminosity = 0.0, H2IILuminosity = 0.0;
	double HMSigma = 0.0, H2IISigma = 0.0;
	if(RS->Energy[ebin] < 11.2) {
	  H2IISED = RS->SED[ebin];
	  IRSED = RS->SED[ebin];
	}
	else if(RS->Energy[ebin] < 13.6) {
	  H2IISED = RS->SED[ebin];
	  LWSED = RS->SED[ebin];
	}
	HMSigma = CalculateIRCrossSection(RS->Energy[ebin]);
	H2IISigma = CalculateH2IICrossSection(RS->Energy[ebin]);
	// Absorb the unit conversions into the cross-section
	H2IISigma *= (double)TimeUnits / ((double)LengthUnits * (double)LengthUnits);
	HMSigma*= (double)TimeUnits / ((double)LengthUnits * (double)LengthUnits);
	IRLuminosity = IRSED * RS->Luminosity / LConv;
	LWLuminosity = LWSED * RS->Luminosity / LConv;
	H2IILuminosity = IRLuminosity + LWLuminosity;
#if(DEBUG)
	printf("MyProcessorNumber %d: IRSED = %f\t IRLuminosity = %g\n", MyProcessorNumber, IRSED,
	       IRLuminosity);
	printf("MyProcessorNumber %d: LWSED = %f\t LWLuminosity = %g\n", MyProcessorNumber, LWSED,
	       LWLuminosity);
#endif
	/* Pre-calculate distances from cells to source */

	for (dim = 0; dim < GridRank; dim++)
	  for (i = 0, index = GridStartIndex[dim]; i < ActiveDims[dim]; 
	       i++, index++) {
	    
	    // Calculate dr_i first, then square it
	    ddr2[dim][i] = 
	      fabs(CellLeftEdge[dim][index] + 0.5*CellWidth[dim][0] - 
		   RS->Position[dim]);
	    ddr2[dim][i] = min(ddr2[dim][i], DomainWidth[dim]-ddr2[dim][i]);
	    ddr2[dim][i] = ddr2[dim][i] * ddr2[dim][i];
	  }

	double radius2, radius2_yz;
	double colden = 0.0, shield = 1.0, b = 0.0, b5 = 0.0, XN = 0.0;
	double H2mass = mh*2.0, alpha = 1.1;
	double kph_hm = 0.0, kdiss_H2II = 0.0;
	int TemperatureField = 0;

        double kdiss_HDI = 0.0, kdiss_CO = 0.0, kdiss_OH = 0.0, kdiss_H2O = 0.0;
        double shieldHDI = 1.0, shieldCO = 1.0, shieldOH = 1.0, shieldH2O = 1.0;

	/* Pre-compute some quantities to speed things up */
	if (RadiativeTransferH2ShieldType > 0) {
		TemperatureField = this->GetTemperatureFieldNumberForH2Shield();
	}
	kdiss_r2 = (float) (LWLuminosity * H2ISigma / (4.0 * pi));
	kph_hm = (float) (IRLuminosity * HMSigma / (4.0 * pi));
	kdiss_H2II = (float) (H2IILuminosity * H2IISigma / (4.0 * pi));
        if(MultiSpecies > 2)
	  kdiss_HDI = (float) (LWLuminosity * HDISigma / (4.0 * pi));
        if(MetalChemistry) {
	  kdiss_CO  = (float) (LWLuminosity * COSigma  / (4.0 * pi));
	  kdiss_OH  = (float) (LWLuminosity * OHSigma  / (4.0 * pi));
	  kdiss_H2O = (float) (LWLuminosity * H2OSigma / (4.0 * pi));
        }
	for (k = 0; k < ActiveDims[2]; k++) {
	  for (j = 0; j < ActiveDims[1]; j++) {
	    radius2_yz = ddr2[1][j] + ddr2[2][k];
	    index = GRIDINDEX(0, j, k);
	    for (i = 0; i < ActiveDims[0]; i++, index++) {
	      radius2 = radius2_yz + ddr2[0][i];
	      if (radius2 < dilRadius2) {
		BaryonField[kdissH2INum][index] += kdiss_r2 / dilRadius2;
		BaryonField[kdissH2IINum][index] += kdiss_H2II / dilRadius2;
		BaryonField[kphHMNum][index] += kph_hm / dilRadius2;
                if(MultiSpecies > 2)
	          BaryonField[kdissHDINum][index] += kdiss_HDI / dilRadius2;
                if(MetalChemistry) {
	          BaryonField[kdissCONum ][index] += kdiss_CO  / dilRadius2;
	          BaryonField[kdissOHNum ][index] += kdiss_OH  / dilRadius2;
	          BaryonField[kdissH2ONum][index] += kdiss_H2O / dilRadius2;
                }
	      }
	      else {
		BaryonField[kdissH2INum][index] += kdiss_r2 / radius2;
		BaryonField[kdissH2IINum][index] += kdiss_H2II / radius2;
		BaryonField[kphHMNum][index] += kph_hm / radius2;
                if(MultiSpecies > 2)
	          BaryonField[kdissHDINum][index] += kdiss_HDI / radius2;
                if(MetalChemistry) {
	          BaryonField[kdissCONum ][index] += kdiss_CO  / radius2;
	          BaryonField[kdissOHNum ][index] += kdiss_OH  / radius2;
	          BaryonField[kdissH2ONum][index] += kdiss_H2O / radius2;
                }
	      }
	      /* Include Shielding */
	      //printf("%s: kdissH2I = %e for Grid %p\n", __FUNCTION__, BaryonField[kdissH2INum][index], this);
	      if (RadiativeTransferH2ShieldType > 0) {
#if(JEANS_LENGTH)
	        l_char = JeansLength(BaryonField[TemperatureField][index],
		  BaryonField[DensNum][index], DensityUnits)*RadiativeTransferOpticallyThinH2CharLength; //cm
#else
	        l_char = CellWidth[0][0]*LengthUnits/2.0; //cm
#endif
	        N_H2 = BaryonField[H2INum][index]*l_char*DensityUnits/(2.0*mh);
                if(MultiSpecies > 2)
	          N_H2 += BaryonField[HDINum][index]*l_char*DensityUnits/(3.0*mh);
#if(DEBUG)
		  printf("l_char = %g\n", l_char);
		  printf("CellWidth = %g\n", CellWidth[0][0]*LengthUnits);
		  printf("l_char/cellwidth = %g\n", l_char/(CellWidth[0][0]*LengthUnits));
		  printf("N_H2 = %g (cm^-2)\n", N_H2);
		  printf("H2I Density[%d] = %g (code)\n", index, BaryonField[H2INum][index]);
		  printf("Local H2I Density = %g (cgs)\n",
		  BaryonField[H2INum][index]*DensityUnits*CellWidth[0][0]*LengthUnits/mh);
		  printf("Temp = %f\n", BaryonField[TemperatureField][index]);
#endif
	        shield = 1.0;
                if(MultiSpecies > 2)
	          shieldHDI = 1.0;
	        if(N_H2 >= THRESHOLD_DENSITY_DB37) {
		  b = sqrt(2.0*kboltz*BaryonField[TemperatureField][index]/H2mass);
		  b5 = b/1e5;
		  XN = N_H2/THRESHOLD_DENSITY_DB37;
		  shield = 0.965/pow(1+XN/b5, alpha) + (0.035/sqrt(1+XN))*exp(-8.5e-4*sqrt(1+XN));
                  if(MultiSpecies > 2)
	            shieldHDI = 0.965/pow(1+XN/(b5 * 0.8165), alpha) + (0.035/sqrt(1+XN))*exp(-8.5e-4*sqrt(1+XN));
                         // multiply by sqrt(2/3)
	        }
	        BaryonField[kdissH2INum][index] *= shield;
	        BaryonField[kdissHDINum][index] *= shieldHDI;
	      } // end H2 shield

              if(MetalChemistry) {
	        N_CO  = BaryonField[CONum ][index]*l_char*DensityUnits/(28.0*mh);
	        N_OH  = BaryonField[OHNum ][index]*l_char*DensityUnits/(17.0*mh);
	        N_H2O = BaryonField[H2ONum][index]*l_char*DensityUnits/(18.0*mh);

	        shieldCO  = 1.0;
	        shieldOH  = 1.0;
	        shieldH2O = 1.0;

                if (N_CO >= THRESHOLD_DENSITY_CO_CO) {
                  XN = N_CO / N_CO_CO;
                  shieldCO = pow(1.0 + XN, a_CO_CO) * exp(b_CO_CO * XN);
                }
                if(MultiSpecies > 1) {
                  if (N_H2 >= THRESHOLD_DENSITY_CO_H2) {
                    XN = N_H2 / N_CO_H2;
                    shieldCO *= pow(1.0 + XN, a_CO_H2) * exp(b_CO_H2 * XN);
                  }
                }

                if (N_OH >= THRESHOLD_DENSITY_OH_OH) {
                  XN = N_OH / N_OH_OH;
                  shieldOH = pow(1.0 + XN, a_OH_OH) * exp(b_OH_OH * XN);
                }
                if(MultiSpecies > 1) {
                  if (N_H2 >= THRESHOLD_DENSITY_OH_H2) {
                    XN = N_H2 / N_OH_H2;
                    shieldOH *= pow(1.0 + XN, a_OH_H2) * exp(b_OH_H2 * XN);
                  }
                }

                if (N_H2O >= THRESHOLD_DENSITY_H2O_H2O) {
                  XN = N_H2O / N_H2O_H2O;
                  shieldH2O = pow(1.0 + XN, a_H2O_H2O);
                }

	        BaryonField[kdissCONum ][index] *= shieldCO;
	        BaryonField[kdissOHNum ][index] *= shieldOH;
	        BaryonField[kdissH2ONum][index] *= shieldH2O;
              }
	    } //end i loop
	  } //end j loop
	} //end k loop
      }
    } //end loop over sources
  } //end loop over APs

  for (dim = 0; dim < GridRank; dim++)
    delete [] ddr2[dim];

  return SUCCESS;
  
}

static double CalculateH2IICrossSection(float Energy)
{
  float X = 0.0;  
  float a = 3.35485518, b = 0.93891875, c = -0.01176537;
  /* Fit created from Stancil et al. 1994 */
  X = log(Energy/8.0);
  return a*pow(10.0, -b*X*X)*exp(-c*X) * 1e-18;

}
static double CalculateIRCrossSection(float Energy)
{
  float X = 0.0;
  float A = 3.486e-16;
  /* Fit taken from Tegmark et al. (1997) */
  X = Energy / 0.74;
  return (A*pow((X - 1), 1.5)/pow(X, 3.11));
}

/*
 * Calculate the Jeans Length of a gas cell and return in cgs 
*/

static double JeansLength(float T, float dens, float density_units)
{
  float jeans_length = 0.0;

  jeans_length = 15*kboltz*T/(4.0*pi*GravConst*mh*dens*density_units);
  return sqrt(jeans_length);
}
  
