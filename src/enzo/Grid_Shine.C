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

#define ONE_ENERGY

RadiationSourceEntry* DeleteRadiationSource(RadiationSourceEntry *RS);
int  pix2vec_nest(long nside,long ipix, FLOAT *v);
FLOAT FindCrossSection(int type, float energy);

int grid::Shine(RadiationSourceEntry *RadiationSource)
{
  int BasePackages, NumberOfNewPhotonPackages;
  int i, j, dim;
  int count=0;
  int min_level = RadiativeTransferInitialHEALPixLevel;

  /* base number of rays to star with: for min_level=2 this is 192
     photon packages per source */

  BasePackages = 12*(int)pow(4,min_level);

  int stype = 3;
#ifdef ONE_ENERGY
  stype = 1;
#endif
  if (MultiSpecies>1 && !RadiativeTransferOpticallyThinH2) stype++;

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
  if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
                      HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL) {
    fprintf(stdout, "Error in grid->IdentifySpeciesFields.\n");
    ENZO_FAIL("");
  }

  /* allocate temporary array of pointers */
  FLOAT AllCellWidth[3];
  for (dim=0; dim<GridRank; dim++) 
    AllCellWidth[dim] =(GridRightEdge[dim] - GridLeftEdge[dim])/
      FLOAT(GridEndIndex[dim] - GridStartIndex[dim] + 1);

  srand(time(NULL));
  
  if (DEBUG) fprintf(stdout, "grid::Shine: Loop over sources and packages \n");

  int ebin;
  RadiationSourceEntry *RS;
  FLOAT FuzzyLength;
  FLOAT ShakeSource[3];
  double RampPercent = 1;
  RS = RadiationSource;
  if (RS != NULL) {

#ifdef ONE_ENERGY
    //    RS->SED[0] = 1.0;
#endif

    if (PhotonTime < (RS->CreationTime + RS->RampTime)) {
      float t = PhotonTime-RS->CreationTime+dtPhoton;
      float frac = t / (RS->RampTime+dtPhoton);
      RampPercent = (exp(frac)-1) / (M_E-1);   // M_E = e = 2.71828...
      RampPercent = max(min(RampPercent, 1), 0);
    }

    /* Shake source within the grid cell every time it shines */

    for (dim = 0; dim < MAX_DIMENSION; dim++)
      //ShakeSource[dim] = 0.0;
      ShakeSource[dim] = (-0.01 + 0.02*float(rand())/RAND_MAX) * CellWidth[dim][0];

    switch (RS->Type) {
    case PopII:
      break;
    case PopIII:
      if (MyProcessorNumber == ProcessorNumber)
	printf("Shine: ramp = %lf, lapsed = %lf/%"FSYM", L = %"GSYM"\n", RampPercent,
	       PhotonTime-RS->CreationTime+dtPhoton, RS->LifeTime, 
	       RS->Luminosity);
      break;
    case BlackHole:
      if (RadiativeTransferInterpolateField) {
	FieldsToInterpolate[HINum] = TRUE;
	FieldsToInterpolate[HeINum] = TRUE;
	FieldsToInterpolate[HeIINum] = TRUE;
      } // ENDIF interpolate fields
//      if (MyProcessorNumber == ProcessorNumber)
//	printf("Shine: ramp = %lf, lapsed = %lf\n", RampPercent,
//	       PhotonTime-RS->CreationTime+dtPhoton);
      break;
    } // ENDSWITCH type

    int ebin;

    for (i=0; i < stype; i++) {

      float photons_per_package;
//      ebin = (i == stype-1 && !RadiativeTransferOpticallyThinH2 && 
//	      MultiSpecies > 1) ? 3 : i;
      ebin = (i == stype-1 && !RadiativeTransferOpticallyThinH2 && 
	      MultiSpecies > 1 && RS->Type == 1) ? 3 : i;

      photons_per_package = RampPercent * RS->Luminosity * 
	RS->SED[ebin] * dtPhoton / float(BasePackages);

      if (ebin == 0)
	EscapedPhotonCount[0] += photons_per_package * BasePackages;

      if (DEBUG)
	printf("Shine: Photons/package[%"ISYM"]: %"GSYM" eV, %"GSYM", %"GSYM", %"GSYM"\n", 
	       ebin, RS->Energy[ebin], RampPercent*RS->Luminosity, 
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
	//      for (j=0; j<1; j++) {
	if (photons_per_package>tiny_number) {
	  PhotonPackageEntry *NewPack = new PhotonPackageEntry;
	  NewPack->NextPackage = PhotonPackages->NextPackage;
	  PhotonPackages->NextPackage = NewPack;
	  NewPack->PreviousPackage = PhotonPackages;
	  if (NewPack->NextPackage != NULL) 
	    NewPack->NextPackage->PreviousPackage  = NewPack;
	  NewPack->Photons = photons_per_package;

	  // Type 4 = X-Ray
	  NewPack->Type = (RS->Type == 2 && i == 0) ? 4 : ebin;

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
	    fprintf(stderr, 
		    "grid::WalkPhotonPackage:  pix2vec_nest outor %"ISYM" %"ISYM" %"GSYM" %"ISYM"\n",
		    (long) (pow(2,NewPack->level)), NewPack->ipix, 
		    NewPack->Photons, NewPack );
	    NewPack->Photons=-1;
	    ENZO_FAIL("");
	  }

	  NewPack->CrossSection = 
	    FindCrossSection(NewPack->Type, NewPack->Energy);

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
	  count++;
	} // if enough photons
      } // Loop over BasePackages
    } //Loop over energy bins
    RS = RS->NextSource;
  }; // if it is valid Sources (not NULL)
  if (DEBUG) printf("grid::Shine: created %"ISYM" packages \n", count);

  // Set the new number of photon packages on this grid
  NumberOfPhotonPackages += NumberOfNewPhotonPackages;

  PhotonPackageEntry *PP ;
  PP = PhotonPackages;
  if (DEBUG) {
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
