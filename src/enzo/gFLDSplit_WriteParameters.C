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
/  Gray Flux-Limited Diffusion Split Implicit Problem Class problem 
/  Parameter output routine
/
/  written by: Daniel Reynolds
/  date:       July 2009
/
/  PURPOSE: Writes all necessary internal parameters for problem 
/           restart.
/
************************************************************************/
#ifdef TRANSFER
#include "gFLDSplit.h"

int gFLDSplit::WriteParameters(FILE *fptr)
{

//   if (debug)  printf("Entering gFLDSplit::WriteParameters routine\n");
  
  fprintf(fptr, "RadHydroESpectrum = %"ISYM"\n", ESpectrum);
  fprintf(fptr, "RadHydroChemistry = %"ISYM"\n", Nchem);
  fprintf(fptr, "RadHydroHFraction = %22.16e\n", HFrac);
  fprintf(fptr, "RadHydroModel = %"ISYM"\n", Model);

  // set restart initial time step to current time step
  fprintf(fptr, "RadHydroInitDt = %22.16e\n", dtrad);
  fprintf(fptr, "RadHydroMaxDt = %22.16e\n", maxdt);
  fprintf(fptr, "RadHydroMinDt = %22.16e\n", mindt);
  fprintf(fptr, "RadHydroMaxSubcycles = %"FSYM"\n", maxsubcycles);
  fprintf(fptr, "RadHydroMaxChemSubcycles = %"FSYM"\n", maxchemsub);
  fprintf(fptr, "RadHydroDtNorm = %22.16e\n", dtnorm);
  fprintf(fptr, "RadHydroDtGrowth = %22.16e\n", dtgrowth);
  fprintf(fptr, "RadHydroDtRadFac = %22.16e\n", dtfac[0]);
  fprintf(fptr, "RadHydroDtGasFac = %22.16e\n", dtfac[1]);
  fprintf(fptr, "RadHydroDtChemFac = %22.16e\n", dtfac[2]);

  fprintf(fptr, "RadiationScaling = %22.16e\n", ErScale);
  fprintf(fptr, "EnergyCorrectionScaling = %22.16e\n", ecScale);
  fprintf(fptr, "ChemistryScaling = %22.16e\n", NiScale);
  if (autoScale) {
    fprintf(fptr, "AutomaticScaling = 1\n");
  } else {
    fprintf(fptr, "AutomaticScaling = 0\n");
  }

  fprintf(fptr, "RadHydroTheta = %22.16e\n", theta);

  fprintf(fptr, "RadiationBoundaryX0Faces = %"ISYM" %"ISYM"\n", 
	  BdryType[0][0], BdryType[0][1]);
  if (rank > 1) {
    fprintf(fptr, "RadiationBoundaryX1Faces = %"ISYM" %"ISYM"\n", 
	    BdryType[1][0], BdryType[1][1]);
    if (rank > 2) {
      fprintf(fptr, "RadiationBoundaryX2Faces = %"ISYM" %"ISYM"\n", 
	      BdryType[2][0], BdryType[2][1]);
    }
  }

  fprintf(fptr, "RadHydroInitialGuess = %"ISYM"\n", initial_guess);    
  fprintf(fptr, "RadHydroKrylovMethod = %"ISYM"\n", Krylov_method);
  fprintf(fptr, "RadHydroSolTolerance = %22.16e\n", sol_tolerance);
  fprintf(fptr, "RadHydroMaxMGIters = %i\n", sol_maxit);    
  fprintf(fptr, "RadHydroMGRelaxType = %i\n", sol_rlxtype);    
  fprintf(fptr, "RadHydroMGPreRelax = %i\n", sol_npre);    
  fprintf(fptr, "RadHydroMGPostRelax = %i\n", sol_npost);    

  fprintf(fptr, "EnergyOpacityC0 = %22.16e\n", EnergyOpacityC0);
  fprintf(fptr, "EnergyOpacityC1 = %22.16e\n", EnergyOpacityC1);
  fprintf(fptr, "EnergyOpacityC2 = %22.16e\n", EnergyOpacityC2);

  // if doing an ionization problem (ProblemTypes 410-415),
  // output additional parameters 
  if ((ProblemType >= 410) && (ProblemType <= 415)) {
    fprintf(fptr, "NGammaDot = %22.16e\n", NGammaDot);
    fprintf(fptr, "EtaRadius = %22.16e\n", EtaRadius);
    fprintf(fptr, "EtaCenter = %22.16e %22.16e %22.16e\n",
	    EtaCenter[0], EtaCenter[1], EtaCenter[2]);
  }

  // output relevant units: although these aren't required for restart, 
  // cosmology runs never output the units (why?), making data analyisis tough
  fprintf(fptr, "DensityUnits = %22.16e\n", DenUnits);
  fprintf(fptr, "LengthUnits = %22.16e\n",  LenUnits);
  fprintf(fptr, "TimeUnits = %22.16e\n",    TimeUnits);

  return SUCCESS;
}
#endif
