/***********************************************************************
/
/  WRITES RADIATIVE TRANSFER PARAMETERS TO AN OUTPUT FILE
/
/  written by: Tom Abel
/  date:       April, 2004
/  modified1:
/
/  PURPOSE:
/
/  NOTE: 
/
************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

/* function prototypes */

int RadiativeTransferWriteParameters(FILE *fptr)
{

  fprintf(fptr, "dtPhoton                                  = %"GOUTSYM"\n",
	  dtPhoton);
  fprintf(fptr, "RadiativeTransferLoadBalance              = %"ISYM"\n", 
	  RadiativeTransferLoadBalance);
  fprintf(fptr, "RadiativeTransferRadiationPressure        = %"ISYM"\n", 
	  RadiationPressure);
  fprintf(fptr, "RadiativeTransferRadiationPressureScale   = %"FSYM"\n", 
	  RadiationPressureScale);
  fprintf(fptr, "RadiativeTransferSourceRadius             = %"GSYM"\n", 
	  RadiativeTransferSourceRadius);
  fprintf(fptr, "RadiativeTransferPropagationSpeedFraction = %"GSYM"\n", 
	  RadiativeTransferPropagationSpeedFraction);
  fprintf(fptr, "RadiativeTransferPropagationDistance      = %"GOUTSYM"\n", 
	  RadiativeTransferPropagationDistance);
  fprintf(fptr, "RadiativeTransferCoupledRateSolver        = %"ISYM"\n", 
	  RadiativeTransferCoupledRateSolver);
  fprintf(fptr, "RadiativeTransferOpticallyThinH2          = %"ISYM"\n", 
	  RadiativeTransferOpticallyThinH2);
  fprintf(fptr, "RadiativeTransferOpticallyThinH2CharLength = %"GOUTSYM"\n", 
	  RadiativeTransferOpticallyThinH2CharLength);
  fprintf(fptr, "RadiativeTransferFLDCallOnLevel           = %"ISYM"\n", 
	  RadiativeTransferFLDCallOnLevel);
  fprintf(fptr, "RadiativeTransferPeriodicBoundary         = %"ISYM"\n", 
	  RadiativeTransferPeriodicBoundary);
  fprintf(fptr, "RadiativeTransferSplitPhotonRadius        = %"FSYM"\n", 
	  RadiativeTransferSplitPhotonRadius);
  fprintf(fptr, "RadiativeTransferFluxBackgroundLimit      = %"FSYM"\n", 
	  RadiativeTransferFluxBackgroundLimit);
  fprintf(fptr, "RadiativeTransferHubbleTimeFraction       = %"FSYM"\n", 
	  RadiativeTransferHubbleTimeFraction);
  fprintf(fptr, "RadiativeTransferRaysPerCell              = %"FSYM"\n", 
	  RadiativeTransferRaysPerCell);
  fprintf(fptr, "RadiativeTransferTimestepVelocityLimit    = %"FSYM"\n", 
	  RadiativeTransferTimestepVelocityLimit);
  fprintf(fptr, "RadiativeTransferTimestepVelocityLevel    = %"ISYM"\n", 
	  RadiativeTransferTimestepVelocityLevel);
  fprintf(fptr, "RadiativeTransferInitialHEALPixLevel      = %"ISYM"\n", 
	  RadiativeTransferInitialHEALPixLevel);
  fprintf(fptr, "RadiativeTransferPhotonEscapeRadius       = %"FSYM"\n", 
	  RadiativeTransferPhotonEscapeRadius);
  fprintf(fptr, "RadiativeTransferInterpolateField         = %"ISYM"\n", 
	  RadiativeTransferInterpolateField);
  fprintf(fptr, "RadiativeTransferSourceClustering         = %"ISYM"\n", 
	  RadiativeTransferSourceClustering);
  fprintf(fptr, "RadiativeTransferPhotonMergeRadius        = %"FSYM"\n", 
	  RadiativeTransferPhotonMergeRadius);
  fprintf(fptr, "RadiativeTransferSourceBeamAngle          = %"FSYM"\n", 
	  RadiativeTransferSourceBeamAngle);
  fprintf(fptr, "RadiativeTransferHIIRestrictedTimestep    = %"ISYM"\n", 
	  RadiativeTransferHIIRestrictedTimestep);
  fprintf(fptr, "RadiativeTransferAdaptiveTimestep         = %"ISYM"\n",
	  RadiativeTransferAdaptiveTimestep);
  fprintf(fptr, "RadiativeTransferHydrogenOnly             = %"ISYM"\n", 
	  RadiativeTransferHydrogenOnly);
  fprintf(fptr, "RadiativeTransferTraceSpectrum            = %"ISYM"\n", 
	  RadiativeTransferTraceSpectrum);
  fprintf(fptr, "RadiativeTransferRayMaximumLength         = %"FSYM"\n", 
	  RadiativeTransferRayMaximumLength);
  fprintf(fptr, "RadiativeTransferH2ShieldType             = %"ISYM"\n", 
	  RadiativeTransferH2ShieldType);
  fprintf(fptr, "RadiativeTransferUseH2Shielding           = %"ISYM"\n", 
	  RadiativeTransferUseH2Shielding);
  fprintf(fptr, "RadiativeTransferH2IIDiss                 = %"ISYM"\n", 
	  RadiativeTransferH2IIDiss);
  fprintf(fptr, "RadiativeTransferTraceSpectrumTable       = %s\n\n", 
	  RadiativeTransferTraceSpectrumTable);

  return SUCCESS;
}
