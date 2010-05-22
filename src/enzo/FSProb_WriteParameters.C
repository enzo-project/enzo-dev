/*****************************************************************************
 *                                                                           *
 * Copyright 2009 Daniel R. Reynolds                                         *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Free-streaming Radiation Implicit Problem Class
/  Problem parameter output routine
/
/  written by: Daniel Reynolds
/  date:       March, 2009
/
/  PURPOSE: Writes all necessary internal parameters for problem 
/           restart.
/
************************************************************************/
#ifdef TRANSFER
#include "FSProb.h"

int FSProb::WriteParameters(FILE *fptr)
{

//   if (debug)  printf("Entering FSProb::WriteParameters routine\n");
  
  // set restart initial time step to current time step
  fprintf(fptr, "FSRadiationScaling = %22.16e\n", EScale);
  fprintf(fptr, "FSRadiationTheta = %22.16e\n", theta);
  fprintf(fptr, "FSRadiationOpacity = %22.16e\n", kappa);
  fprintf(fptr, "FSRadiationNGammaDot = %22.16e\n", NGammaDot);
  fprintf(fptr, "FSRadiationEtaRadius = %22.16e\n", EtaRadius);
  fprintf(fptr, "FSRadiationEtaCenter = %22.16e %22.16e %22.16e\n", 
	  EtaCenter[0], EtaCenter[1], EtaCenter[2]);
  fprintf(fptr, "FSRadiationLimiterType = %"ISYM"\n", LimType);

  fprintf(fptr, "FSRadiationBoundaryX0Faces = %"ISYM" %"ISYM"\n", 
	  BdryType[0][0], BdryType[0][1]);
  if (rank > 1) {
    fprintf(fptr, "FSRadiationBoundaryX1Faces = %"ISYM" %"ISYM"\n", 
	    BdryType[1][0], BdryType[1][1]);
    if (rank > 2) {
      fprintf(fptr, "FSRadiationBoundaryX2Faces = %"ISYM" %"ISYM"\n", 
	      BdryType[2][0], BdryType[2][1]);
    }
  }

  fprintf(fptr, "FSRadiationMaxDt = %22.16e\n", maxdt);    
  fprintf(fptr, "FSRadiationInitialGuess = %"ISYM"\n", initial_guess);    
  fprintf(fptr, "FSRadiationTolerance = %22.16e\n", sol_tolerance);    
  fprintf(fptr, "FSRadiationMaxMGIters = %i\n", sol_maxit);    
  fprintf(fptr, "FSRadiationMGRelaxType = %i\n", sol_rlxtype);    
  fprintf(fptr, "FSRadiationMGPreRelax = %i\n", sol_npre);    
  fprintf(fptr, "FSRadiationMGPostRelax = %i\n", sol_npost);    

  return SUCCESS;
}
#endif
