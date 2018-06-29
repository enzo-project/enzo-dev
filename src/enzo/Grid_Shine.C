#define DEBUG 0
#define MYPROC MyProcessorNumber == ProcessorNumber
/***********************************************************************
/
/  GRID CLASS (CREATES PHOTON PACKES AT RADIATION SOURCE POSITION)
/
/  written by: Tom Abel
/  date:       August, 2003
/  modified1:
/
/  PURPOSE:  Creates PhotonPackages for a give radiation source
/
/  INPUTS:
/    RadiationSource - type RadiationSourceEntry defining source
************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "RadiativeTransferHealpixRoutines64.h"
#define MAX_HEALPIX_LEVEL 29



RadiationSourceEntry* DeleteRadiationSource(RadiationSourceEntry *RS);
FLOAT FindCrossSection(int type, float energy);

int grid::Shine(RadiationSourceEntry *RadiationSource)
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  const float EnergyThresholds[] = {13.6, 24.6, 54.4, 100.0}; //Only used for determining HI,HeI,HeII
 
  RadiationSourceEntry *RS = RadiationSource;
  FLOAT min_beam_zvec, dot_prod;
  double vec[3];
  long long BasePackages = 0, NumberOfNewPhotonPackages = 0;
  int dim = 0;
  int64_t ray = 0;
  int count = 0;
  int min_level = RadiativeTransferInitialHEALPixLevel;

  /* base number of rays to star with: for min_level=2 this is 192
     photon packages per source */

  BasePackages = 12*(int)pow(4,min_level);

  /* If using a beamed source, calculate the minimum z-component of
     the ray normal (always beamed in the polar coordinate). */

  if (RS->Type == Beamed)
    min_beam_zvec = cos(M_PI * RadiativeTransferSourceBeamAngle / 180.0);

  int stype = RS->EnergyBins;

  /* At most how many new Photon Packages should be allocated and
     created?  */
  
  NumberOfNewPhotonPackages = BasePackages*stype;
  if (MYPROC && DEBUG) 
    fprintf(stdout, "grid::Shine: Maximum Number of New Photon Packages %"ISYM"\n",
	    NumberOfNewPhotonPackages);

  if (MyProcessorNumber != ProcessorNumber) {
    NumberOfPhotonPackages += NumberOfNewPhotonPackages; 
    return SUCCESS;
  }

  /* Find Multi-species fields. */

  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
    DINum, DIINum, HDINum;
  IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
			HMNum, H2INum, H2IINum, DINum, DIINum, HDINum);

  /* allocate temporary array of pointers */
  FLOAT AllCellWidth[3];
  for (dim=0; dim<GridRank; dim++) 
    AllCellWidth[dim] =(GridRightEdge[dim] - GridLeftEdge[dim])/
      FLOAT(GridEndIndex[dim] - GridStartIndex[dim] + 1);

  if (rand_init == 0) {
    srand( time(NULL) );
    rand_init = 1;
  }
  
  if (MYPROC && DEBUG) fprintf(stdout, "grid::Shine: Loop over sources and packages \n");

  int ebin = 0, this_type = 0, type_count = 0;
  FLOAT FuzzyLength = 0.0;
  FLOAT ShakeSource[3];
  double RampPercent = 1;

  if (RS->Type == Episodic) {
    const float sigma_inv = 4.0;
    float t = PhotonTime - RS->CreationTime + dtPhoton;
    float frac = 2.0 * fabs(t - round(t/RS->RampTime) * RS->RampTime) /
      RS->RampTime;
    RampPercent = exp((frac-1)*sigma_inv);
  } // ENDIF episodic
  else if (PhotonTime < (RS->CreationTime + RS->RampTime)) {   
    float t = PhotonTime-RS->CreationTime+dtPhoton;
    float frac = t / (RS->RampTime+dtPhoton);
    RampPercent = (exp(frac)-1) / (M_E-1);   // M_E = e = 2.71828...
    RampPercent = max(min(RampPercent, 1), 0);
  }

  /* Shake source within the grid cell every time it shines */

  for (dim = 0; dim < MAX_DIMENSION; dim++)
    ShakeSource[dim] = 0.0;
  //ShakeSource[dim] = (-0.01 + 0.02*float(rand())/RAND_MAX) * CellWidth[dim][0];

  if (MYPROC && DEBUG)
    printf("Shine: ramp = %lf, lapsed = %lf/%"FSYM", L = %"GSYM"\n", RampPercent, 
	   PhotonTime-RS->CreationTime+dtPhoton, RS->LifeTime, 
	   RS->Luminosity);


  /* Loop over Energy Bins  - stype is the number of energy bins */ 
  /* The types are:
   * Type0: HI Ionising Radiation
   * Type1: HeI Ionising Radiation
   * Type2: HeII Ionising Radiation
   * Type3: H2I photo-dissociating Radiation
   * Type4: Infrared Radiation (HM Photodetachment)
   * Type5: XRAYS (which heat and ionise HI, HeI and HeII also)
   * Type6: 
   */
  for (ebin=0; ebin < stype; ebin++) {

    if (MYPROC && DEBUG)
      fprintf(stdout, "Shine: Energy = %f eV", RS->Energy[ebin]);
    float photons_per_package;

    if(RS->Energy[ebin] <= 0)
      continue;
    /* If we are doing simple H2I, H2II and HM rates continue here. */
    if(RS->Energy[ebin] <= 13.6 && RadiativeTransferOpticallyThinH2 == 1)
      continue;
    
    photons_per_package = RampPercent * RS->Luminosity * 
      RS->SED[ebin] * dtPhoton / float(BasePackages);

    if (ebin == 0)
      EscapedPhotonCount[0] += photons_per_package * BasePackages;

    /* 
     * Associate Energy Bin with type e.g. IR -> XRAYS 
     */
    if(RS->Energy[ebin] >= 0.01 && RS->Energy[ebin] < 11.2) {
      this_type = IR; /* IR Case */
    }
    else if(RS->Energy[ebin] < 13.6) {
      this_type = LW;  /* LW Case */
    }
    else if (RS->Energy[ebin] < 100.0) { //Set iHI or iHeI or iHeII
      for (type_count = 0; type_count < 3; type_count++)
	{
	  if (RS->Energy[ebin] >= EnergyThresholds[type_count] &&
	      RS->Energy[ebin] <  EnergyThresholds[type_count+1]) {
	    this_type = type_count;
	    break;
	  }
	}
    }
    else if (RS->Energy[ebin] >= 100.0)
      this_type = XRAYS;
    else {
      fprintf(stderr, "Undefined energy bin found\n");
      return FAIL;
    }

    if (MYPROC && DEBUG)
      {
	fprintf(stdout, "Shine: Photons/package[%"ISYM"]: %"GSYM" eV, Luminosity = %"GSYM"\n " \
		"Ramp Luminosity = %"GSYM" \n " \
		"SED = %"GSYM"\n Photons per Package = %"GSYM"\n Type = %"ISYM"\n",
		ebin, RS->Energy[ebin], RS->Luminosity, RampPercent*RS->Luminosity, 
		RS->SED[ebin], photons_per_package, this_type);
      }

    /* Loop over each Ray */
    for (ray=0; ray<BasePackages; ray++) {

      if (RS->Type == Beamed) {
	pix2vec_nest64((int64_t) (1 << min_level), (int64_t) ray, vec);
	// Dot product of the source orientation (already normalized
	// to 1) and ray normal must be greater than cos(beaming angle)
	dot_prod = 0.0;
	for (dim = 0; dim < MAX_DIMENSION; dim++)
	  dot_prod += RS->Orientation[dim] * vec[dim];
	if (fabs(dot_prod) < min_beam_zvec)
	  continue;
      }

      //      for (j=0; j<1; j++) {
      //	if (photons_per_package>tiny_number) { //removed and changed to below by Ji-hoon Kim in Sep.2009
      if (!isnan(photons_per_package) && photons_per_package > 0) { 
	PhotonPackageEntry *NewPack = new PhotonPackageEntry;
	NewPack->NextPackage = PhotonPackages->NextPackage;
	PhotonPackages->NextPackage = NewPack;
	NewPack->PreviousPackage = PhotonPackages;
	if (NewPack->NextPackage != NULL) 
	  NewPack->NextPackage->PreviousPackage  = NewPack;
	NewPack->Photons = photons_per_package;
	NewPack->Type = this_type;

	// Type 6 = tracing spectrum (check Grid_WalkPhotonPackage)
	if (RadiativeTransferTraceSpectrum) NewPack->Type = TRACINGSPECTRUM;

	// Override if we're only doing hydrogen ionization
	if (RadiativeTransferHydrogenOnly) NewPack->Type = 0;

	NewPack->EmissionTimeInterval = dtPhoton;
	NewPack->EmissionTime = PhotonTime;
	NewPack->CurrentTime  = PhotonTime;
	NewPack->ColumnDensity = 0;
	NewPack->Radius = 0.;
	NewPack->ipix = ray;
	NewPack->level = min_level;
	NewPack->Energy = RS->Energy[ebin];
	NewPack->CrossSection = 0.0;
	double dir_vec[3];
	pix2vec_nest64((int64_t) (1 << NewPack->level), NewPack->ipix, dir_vec);
	/* Find the cross section for each radiation type */
	if (NewPack->Type < 4)
	  NewPack->CrossSection = 
	    FindCrossSection(NewPack->Type, NewPack->Energy);
	else
	  NewPack->CrossSection = tiny_number;
	/* Set the photon origin to the source radius (0 = point src) */

	NewPack->SourcePositionDiff = 0.0;

	for (dim = 0; dim < MAX_DIMENSION; dim++) {
	  FuzzyLength = RadiativeTransferSourceRadius * dir_vec[dim] * 
	    AllCellWidth[dim] + ShakeSource[dim];
	  NewPack->SourcePosition[dim] = RS->Position[dim] + FuzzyLength;
	  NewPack->SourcePositionDiff += FuzzyLength * FuzzyLength;
	}
	NewPack->SourcePositionDiff = sqrt(NewPack->SourcePositionDiff);
	NewPack->CurrentSource = RS->SuperSource;

	/* Consider the first super source with a leaf size greater
	   than the cell size. */

#define NO_PRE_MERGE
#ifdef PRE_MERGE
	while (NewPack->CurrentSource != NULL &&
	       RadiativeTransferPhotonMergeRadius * 
	       NewPack->CurrentSource->ClusteringRadius < CellWidth[0][0])
	  NewPack->CurrentSource = NewPack->CurrentSource->ParentSource;
#endif
	if(MYPROC && DEBUG == 2)  //Dumps info on 
	  NewPack->PrintInfo();

//	if (DEBUG) {
//	  printf("Shine: MBH = %d, RS->Type = %d, E=%g, NewPack->Type = %d\n", 
//	         MBH, RS->Type, RS->Energy[ebin], NewPack->Type);  
//	  NewPack->PrintInfo();
//	}

	count++;
      } // if enough photons
    } // Loop over BasePackages (rays)
  } //Loop over energy bins

 
  // Set the new number of photon packages on this grid
  NumberOfPhotonPackages += NumberOfNewPhotonPackages;

  if (MYPROC && DEBUG) {
    printf("Shine: created %"ISYM" packages \n", count);
    PhotonPackageEntry *PP;
    PP = PhotonPackages;
    count = 0;
    while ((PP->NextPackage) != NULL) {
      count++;
      PP = PP->NextPackage;
    }
    if (DEBUG) fprintf(stdout,"Shine: done.\n");
    if (DEBUG) fprintf(stdout,"counted %"ISYM" packages\n", count);
  }

  if (DEBUG) fprintf(stdout, "Shine: PhotonPackages : %p   NextPackage  %p\n", 
		     PhotonPackages, PhotonPackages->NextPackage);
  
  return SUCCESS;
};
