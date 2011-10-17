#define DEBUG 0
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

RadiationSourceEntry* DeleteRadiationSource(RadiationSourceEntry *RS);
FLOAT FindCrossSection(int type, float energy);

int grid::Shine(RadiationSourceEntry *RadiationSource)
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  const float EnergyThresholds[] = {13.6, 24.6, 54.4, 100.0};

  RadiationSourceEntry *RS = RadiationSource;
  FLOAT min_beam_zvec, dot_prod, vec[3];
  int BasePackages, NumberOfNewPhotonPackages;
  int i, j, dim;
  int count=0;
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
  if (DEBUG) 
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
  
  if (DEBUG) fprintf(stdout, "grid::Shine: Loop over sources and packages \n");

  int ebin, this_type, type_count;
  FLOAT FuzzyLength;
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

  switch (RS->Type) {
  case PopII:
    break;
  case Episodic:
  case PopIII:
    if (MyProcessorNumber == ProcessorNumber)
      printf("Shine: ramp = %lf, lapsed = %lf/%"FSYM", L = %"GSYM"\n", RampPercent,
	     PhotonTime-RS->CreationTime+dtPhoton, RS->LifeTime, 
	     RS->Luminosity);
    break;
  case BlackHole:
    break;
  case MBH:
    if (MyProcessorNumber == ProcessorNumber)
      printf("Shine: ramp = %lf, lapsed = %lf/%"FSYM", L = %"GSYM"\n", RampPercent, 
	     PhotonTime-RS->CreationTime+dtPhoton, RS->LifeTime, 
	     RS->Luminosity);
    break;
  } // ENDSWITCH type

  for (i=0; i < stype; i++) {

    float photons_per_package;

    // Type 3 = H2I_LW
    if (!RadiativeTransferOpticallyThinH2 && MultiSpecies > 1 &&
	RS->Energy[i] < 13.6)
      ebin = 3;
    else
      ebin = i;

    photons_per_package = RampPercent * RS->Luminosity * 
      RS->SED[ebin] * dtPhoton / float(BasePackages);

    if (ebin == 0)
      EscapedPhotonCount[0] += photons_per_package * BasePackages;

    if (RS->Energy[ebin] < 100.0) {
      for (type_count = 0; type_count < 3; type_count++)
	if (RS->Energy[ebin] >= EnergyThresholds[type_count] &&
	    RS->Energy[ebin] <  EnergyThresholds[type_count+1]) {
	  this_type = type_count;
	  break;
	}
    }

    if (DEBUG)
      printf("Shine: Photons/package[%"ISYM"]: %"GSYM" eV, %"GSYM", %"GSYM", %"GSYM", %"GSYM"\n", 
	     ebin, RS->Energy[ebin], RS->Luminosity, RampPercent*RS->Luminosity, 
	     RS->SED[ebin], photons_per_package);

    if (RadiativeTransferInterpolateField)
      switch (ebin) {
      case 0:
	FieldsToInterpolate[HINum] = TRUE;
	break;
      case 1:
	FieldsToInterpolate[HeINum] = TRUE;
	break;
      case 2:
	FieldsToInterpolate[HeIINum] = TRUE;
	break;
      case 3:
	FieldsToInterpolate[H2INum] = TRUE;
	break;
      } // ENDSWITCH ebin
      
    // DEBUG fudge
    for (j=0; j<BasePackages; j++) {

      if (RS->Type == Beamed) {
	if (pix2vec_nest((long) (1 << min_level), (long) j, vec) == FAIL)
	  ENZO_FAIL("Error in pix2vec_nested: beamed source");
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

	// Type 4 = X-Ray
	if (((RS->Type == BlackHole || RS->Type == MBH) && i == 0) ||
	    RS->Energy[ebin] > 100)
	  NewPack->Type = 4;
	else
	  NewPack->Type = this_type;

	// Type 5 = tracing spectrum (check Grid_WalkPhotonPackage)
	if (RadiativeTransferTraceSpectrum) NewPack->Type = 5;  //#####

	// Override if we're only doing hydrogen ionization
	if (RadiativeTransferHydrogenOnly) NewPack->Type = 0;

	NewPack->EmissionTimeInterval = dtPhoton;
	NewPack->EmissionTime = PhotonTime;
	NewPack->CurrentTime  = PhotonTime;
	NewPack->ColumnDensity = 0;
	NewPack->Radius = 0.;
	NewPack->ipix = j;
	NewPack->level = min_level;
	NewPack->Energy = RS->Energy[ebin];
	FLOAT dir_vec[3];
	if (pix2vec_nest((long) (1<<NewPack->level), NewPack->ipix, 
			 dir_vec) == FAIL) {
	  NewPack->Photons=-1;
	  ENZO_VFAIL("grid::WalkPhotonPackage:  pix2vec_nest outor %"ISYM" %"ISYM" %"GSYM" %"ISYM"\n",
		  (long) (pow(2,NewPack->level)), NewPack->ipix, 
		  NewPack->Photons, NewPack )
	}

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

	while (NewPack->CurrentSource != NULL &&
	       NewPack->CurrentSource->ClusteringRadius < CellWidth[0][0])
	  NewPack->CurrentSource = NewPack->CurrentSource->ParentSource;

	if (DEBUG) {
	  printf("Shine: MBH = %d, RS->Type = %d, E=%g, NewPack->Type = %d\n", 
	         MBH, RS->Type, RS->Energy[ebin], NewPack->Type);  
	  NewPack->PrintInfo();
	}

	count++;
      } // if enough photons
    } // Loop over BasePackages
  } //Loop over energy bins

  if (DEBUG) printf("grid::Shine: created %"ISYM" packages \n", count);

  // Set the new number of photon packages on this grid
  NumberOfPhotonPackages += NumberOfNewPhotonPackages;

  if (DEBUG) {
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

  if (DEBUG) fprintf(stdout, "Shine: PhotonPackages : %"ISYM"   NextPackage  %"ISYM"\n", 

		     PhotonPackages, PhotonPackages->NextPackage);
  
  return SUCCESS;
};
