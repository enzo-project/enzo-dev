#define MYPROC MyProcessorNumber == ProcessorNumber
#define DEBUG 0

#define JEANS_LENGTH 1
#define SHIELD 1
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
#ifdef FUTUREAP
#include "ActiveParticle.h"
#include "ActiveParticle_RadiationParticle.h"
#include "ActiveParticle_SmartStar.h"
#endif
#include "phys_constants.h"

#define THRESHOLD_DENSITY_DB36 1e14
#define THRESHOLD_DENSITY_DB37 5e14
#define sigma_unit  2.08854e-4  //convert Msolar pc^-2 -> cm^-2
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
  double l_char = 0.0, N_H2 = 0.0;
  double H2ISigma = 3.71e-18;

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

  /* Get photo-ionization fields */

  int kphHINum, kphHeINum, kphHeIINum, kdissH2INum, kphHMNum, kdissH2IINum;
  int gammaNum;
  IdentifyRadiativeTransferFields(kphHINum, gammaNum, kphHeINum, kphHeIINum, 
				  kdissH2INum, kphHMNum, kdissH2IINum);

  /* For now, initialize H2 photo-dissociation field. */

  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  // Now done in Grid_InitializeRadiativeTransferFields.C
//  for (i = 0; i < size; i++)
//    BaryonField[kdissH2INum][i] = 0;
  
  if (AllStars == NULL && ProblemType != 50
#if FUTUREAP
    && EnabledActiveParticlesCount == 0)
#else
    )
#endif
      return SUCCESS;

  /* If using cosmology, get units. */

  float TemperatureUnits, DensityUnits, LengthUnits, VelocityUnits, 
    TimeUnits, aUnits = 1;

  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, PhotonTime);

  // Absorb the unit conversions into the cross-section
  H2ISigma *= (double)TimeUnits / ((double)LengthUnits * (double)LengthUnits);

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
  float dilutionRadius = 4.848e-6 * pc / (double) LengthUnits;  // 1 AU
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
	/* Pre-compute some quantities to speed things up */
	TemperatureField = this->GetTemperatureFieldNumberForH2Shield();
	kdiss_r2 = (float) (LWLuminosity * H2ISigma / (4.0 * M_PI));
	kph_hm = (float) (IRLuminosity * HMSigma / (4.0 * M_PI));
	kdiss_H2II = (float) (H2IILuminosity * H2IISigma / (4.0 * M_PI));
	for (k = 0; k < ActiveDims[2]; k++) {
	  for (j = 0; j < ActiveDims[1]; j++) {
	    radius2_yz = ddr2[1][j] + ddr2[2][k];
	    index = GRIDINDEX(0, j, k);
	    for (i = 0; i < ActiveDims[0]; i++, index++) {
	      radius2 = radius2_yz + ddr2[0][i];
	      //if (radius2 < outerFront2 && radius2 > innerFront2) {
	      //radius2 = max(radius2, dilRadius2);
	      if (radius2 < dilRadius2) {
		BaryonField[kdissH2INum][index] += kdiss_r2 / dilRadius2;
		BaryonField[kdissH2IINum][index] += kdiss_H2II / dilRadius2;
		BaryonField[kphHMNum][index] += kph_hm / dilRadius2;
	      }
	      else {
		BaryonField[kdissH2IINum][index] += kdiss_H2II / radius2;
		BaryonField[kphHMNum][index] += kph_hm / radius2;
		BaryonField[kdissH2INum][index] += kdiss_r2 / radius2;
	      }
	      /* Include Shielding */

#if JEANS_LENGTH
	      l_char = JeansLength(BaryonField[TemperatureField][index],
		BaryonField[DensNum][index], DensityUnits)*RadiativeTransferOpticallyThinH2CharLength; //cm
#else
	      l_char = CellWidth[0][0]*LengthUnits/2.0; //cm
#endif
	      N_H2 = BaryonField[H2INum][index]*l_char*DensityUnits/mh;
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
	      if(N_H2 >= THRESHOLD_DENSITY_DB37) {
		b = sqrt(2.0*kboltz*BaryonField[TemperatureField][index]/H2mass);
		b5 = b/1e5;
		XN = N_H2/THRESHOLD_DENSITY_DB37;
		shield = 0.965/pow(1+XN/b5, alpha) + (0.035/sqrt(1+XN))*exp(-8.5e-4*sqrt(1+XN));
#if DEBUG		
		printf("%s: shield = %g\t kdiss = %g\t kdiss*shield = %g\n", __FUNCTION__, 
		      shield, BaryonField[kdissH2INum][index], BaryonField[kdissH2INum][index]*shield);
#endif
	      }
	      BaryonField[kdissH2INum][index] *= shield;
	     
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

      kdiss_r2 = (float) (LWLuminosity * H2ISigma / (4.0 * M_PI));
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
	    //} // ENDIF
	  } // END: i-direction
	} // END: j-direction
      } // END: k-direction
    } // ENDFOR stars

  } // ENDELSE AllStars
#if FUTUREAP
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
	/* Pre-compute some quantities to speed things up */
	TemperatureField = this->GetTemperatureFieldNumberForH2Shield();
	kdiss_r2 = (float) (LWLuminosity * H2ISigma / (4.0 * M_PI));
	kph_hm = (float) (IRLuminosity * HMSigma / (4.0 * M_PI));
	kdiss_H2II = (float) (H2IILuminosity * H2IISigma / (4.0 * M_PI));
	for (k = 0; k < ActiveDims[2]; k++) {
	  for (j = 0; j < ActiveDims[1]; j++) {
	    radius2_yz = ddr2[1][j] + ddr2[2][k];
	    index = GRIDINDEX(0, j, k);
	    for (i = 0; i < ActiveDims[0]; i++, index++) {
	      radius2 = radius2_yz + ddr2[0][i];
	      //if (radius2 < outerFront2 && radius2 > innerFront2) {
	      //radius2 = max(radius2, dilRadius2);
	      if (radius2 < dilRadius2) {
		BaryonField[kdissH2INum][index] += kdiss_r2 / dilRadius2;
		BaryonField[kdissH2IINum][index] += kdiss_H2II / dilRadius2;
		BaryonField[kphHMNum][index] += kph_hm / dilRadius2;
	      }
	      else {
		BaryonField[kdissH2INum][index] += kdiss_r2 / radius2;
		BaryonField[kdissH2IINum][index] += kdiss_H2II / radius2;
		BaryonField[kphHMNum][index] += kph_hm / radius2;
	      }
	      /* Include Shielding */
	      //printf("%s: kdissH2I = %e for Grid %p\n", __FUNCTION__, BaryonField[kdissH2INum][index], this);
#if SHIELD
#if(JEANS_LENGTH)
	      l_char = JeansLength(BaryonField[TemperatureField][index],
		BaryonField[DensNum][index], DensityUnits)*RadiativeTransferOpticallyThinH2CharLength; //cm
#else
	      l_char = CellWidth[0][0]*LengthUnits/2.0; //cm
#endif
	      N_H2 = BaryonField[H2INum][index]*l_char*DensityUnits/mh;
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
	      if(N_H2 >= THRESHOLD_DENSITY_DB37) {
		b = sqrt(2.0*kboltz*BaryonField[TemperatureField][index]/H2mass);
		b5 = b/1e5;
		XN = N_H2/THRESHOLD_DENSITY_DB37;
		shield = 0.965/pow(1+XN/b5, alpha) + (0.035/sqrt(1+XN))*exp(-8.5e-4*sqrt(1+XN));
	      }
	      BaryonField[kdissH2INum][index] *= shield;
#endif
	    } //end i loop
	  } //end j loop
	} //end k loop
      }
    } //end loop over sources
  } //end loop over APs
#endif
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

  jeans_length = 15*kboltz*T/(4.0*M_PI*GravConst*mh*dens*density_units);
  return sqrt(jeans_length);
}
  
