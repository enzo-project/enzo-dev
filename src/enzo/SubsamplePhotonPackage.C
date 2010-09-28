#define DEBUG 0
#include <stdlib.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "PhotonPackage.h"

// mode = 0 (SPLIT_RAY) -- splits main ray into subsamples
// mode = 1 (COMBINE_RAY) -- combines subsampling rays to update main ray

int SubsamplePhotonPackage(PhotonPackageEntry *&PP, PhotonPackageEntry **Subrays,
			   int SubrayLevel, int NumberOfSubrays, int mode)
{

  int i, dim, NewLevel;
  long nipix;
  float NewFlux;

  if (PP->PreviousPackage == NULL)
    fprintf(stderr, "SubsamplePhotonPackage: Warning Previous Photon Package"
	    " == NULL\n");
  
  if (DEBUG) 
    fprintf(stdout, "subsample package ipix:%"ISYM" level:%"ISYM"\n",
	    PP->ipix, PP->level);

  if (mode == SPLIT_RAY) {

    NewFlux = PP->Photons / NumberOfSubrays;
    NewLevel = PP->level + SubrayLevel;
    nipix = (PP->ipix)*NumberOfSubrays;  // NumberOfSubrays = 4^SubgridLevel

    for (i = 0; i < NumberOfSubrays; i++) {

      Subrays[i]->Photons         = NewFlux;
      Subrays[i]->Type            = PP->Type;
      Subrays[i]->EmissionTimeInterval = PP->EmissionTimeInterval;
      Subrays[i]->EmissionTime    = PP->EmissionTime;
      Subrays[i]->CurrentTime     = PP->CurrentTime;
      Subrays[i]->ColumnDensity   = PP->ColumnDensity;
      Subrays[i]->Radius          = PP->Radius;
      Subrays[i]->ipix            = nipix++;
      Subrays[i]->level           = NewLevel;
      Subrays[i]->Energy          = PP->Energy;
      Subrays[i]->CrossSection    = PP->CrossSection;
      for (dim = 0; dim < 3; dim++)
	Subrays[i]->SourcePosition[dim] = PP->SourcePosition[dim];
      Subrays[i]->SourcePositionDiff  = PP->SourcePositionDiff;
      Subrays[i]->CurrentSource   = PP->CurrentSource;

    } // ENDFOR (i) childrays

  } // ENDIF SPLIT_RAY

  else if (mode == COMBINE_RAY) {

    /* The main ray gets the sum of the photons left from the
       sub-rays. */

    PP->Photons = 0.0;
    for (i = 0; i < NumberOfSubrays; i++)
      if (Subrays[i]->Photons > 0)
	PP->Photons += Subrays[i]->Photons;

  } // ENDIF COMBINE_RAY

  else
    ENZO_FAIL("Unknown mode in SubsamplePhotonPackage.");
  
  return SUCCESS;
}
