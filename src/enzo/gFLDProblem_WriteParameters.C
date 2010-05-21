/*****************************************************************************
 *                                                                           *
 * Copyright 2006 Daniel R. Reynolds                                         *
 * Copyright 2006 Laboratory for Computational Astrophysics                  *
 * Copyright 2006 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Gray Flux-Limited Diffusion Implicit Problem Class problem 
/  parameter output routine
/
/  written by: Daniel Reynolds
/  date:       March, 2008
/
/  PURPOSE: Writes all necessary internal parameters for problem 
/           restart.
/
************************************************************************/
#ifdef TRANSFER
#include "gFLDProblem.h"

int gFLDProblem::WriteParameters(FILE *fptr)
{

//   if (debug)  printf("Entering gFLDProblem::WriteParameters routine\n");
  
  fprintf(fptr, "RadHydroESpectrum = %"ISYM"\n", ESpectrum);
  fprintf(fptr, "RadHydroChemistry = %"ISYM"\n", Nchem);
  fprintf(fptr, "RadHydroHFraction = %22.16e\n", HFrac);
  fprintf(fptr, "RadHydroModel = %"ISYM"\n", Model);

  // set restart initial time step to current time step
  if (dt == 0.0) 
    fprintf(fptr, "RadHydroInitDt = %22.16e\n", initdt);
  else
    fprintf(fptr, "RadHydroInitDt = %22.16e\n", dt);
  fprintf(fptr, "RadHydroMaxDt = %22.16e\n", maxdt);
  fprintf(fptr, "RadHydroMinDt = %22.16e\n", mindt);
  fprintf(fptr, "RadHydroDtNorm = %22.16e\n", dtnorm);
  fprintf(fptr, "RadHydroDtRadFac = %22.16e\n", dtfac[0]);
  fprintf(fptr, "RadHydroDtGasFac = %22.16e\n", dtfac[1]);
  fprintf(fptr, "RadHydroDtChemFac = %22.16e\n", dtfac[2]);

  fprintf(fptr, "RadiationScaling = %22.16e\n", ErScale);
  fprintf(fptr, "EnergyCorrectionScaling = %22.16e\n", ecScale);
  fprintf(fptr, "ChemistryScaling = %22.16e\n", NiScale);

  fprintf(fptr, "RadHydroTheta = %22.16e\n", theta);
  fprintf(fptr, "RadHydroLimiterType = %"ISYM"\n", LimType);

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

  fprintf(fptr, "RadHydroAprxJacobian = %"ISYM"\n", approx_jac);    
  fprintf(fptr, "RadHydroInitialGuess = %"ISYM"\n", initial_guess);    
  fprintf(fptr, "RadHydroAnalyticChem = %"ISYM"\n", AnalyticChem);
  fprintf(fptr, "RadHydroNewtLinesearch = %"ISYM"\n", newt_linesearch);
  fprintf(fptr, "RadHydroNewtIters = %"ISYM"\n", newt_maxit);
  fprintf(fptr, "RadHydroNewtNorm = %"ISYM"\n", newt_norm);    
  fprintf(fptr, "RadHydroINConst = %22.16e\n", newt_INconst);    
  fprintf(fptr, "RadHydroNewtTolerance = %22.16e\n", newt_tol);    
  fprintf(fptr, "RadHydroMinLinesearch = %22.16e\n", 
	  newt_MinLinesearch);    
  fprintf(fptr, "RadHydroMaxMGIters = %i\n", sol_maxit);    
  fprintf(fptr, "RadHydroMGRelaxType = %i\n", sol_rlxtype);    
  fprintf(fptr, "RadHydroMGPreRelax = %i\n", sol_npre);    
  fprintf(fptr, "RadHydroMGPostRelax = %i\n", sol_npost);    

  fprintf(fptr, "PlanckOpacityC0 = %22.16e\n", PlanckOpacityC0);
  fprintf(fptr, "PlanckOpacityC1 = %22.16e\n", PlanckOpacityC1);
  fprintf(fptr, "PlanckOpacityC2 = %22.16e\n", PlanckOpacityC2);
  fprintf(fptr, "PlanckOpacityC3 = %22.16e\n", PlanckOpacityC3);
  fprintf(fptr, "PlanckOpacityC4 = %22.16e\n", PlanckOpacityC4);

  fprintf(fptr, "EnergyOpacityC0 = %22.16e\n", EnergyOpacityC0);
  fprintf(fptr, "EnergyOpacityC1 = %22.16e\n", EnergyOpacityC1);
  fprintf(fptr, "EnergyOpacityC2 = %22.16e\n", EnergyOpacityC2);
  fprintf(fptr, "EnergyOpacityC3 = %22.16e\n", EnergyOpacityC3);
  fprintf(fptr, "EnergyOpacityC4 = %22.16e\n", EnergyOpacityC4);

  // if doing an ionization problem (ProblemTypes 410-415),
  // output additional parameters 
  if ((ProblemType >= 410) && (ProblemType <= 415)) {
    fprintf(fptr, "NGammaDot = %22.16e\n", IonizationParms[0]);
    fprintf(fptr, "EtaRadius = %22.16e\n", IonizationParms[1]);
    fprintf(fptr, "EtaCenter = %22.16e %22.16e %22.16e\n",
	    IonizationParms[2],IonizationParms[3],IonizationParms[4]);
  }
  
  // if doing a Marshak-type problem (20 <= Model < 30), 
  // output additional Marshak parameters 
  if ( Model >= 20 && Model <= 29 ) {
    fprintf(fptr, "SuOlsonGreyEps = %22.16e", &MarshakParms[0]);
  }

  // output relevant units: although these aren't required for restart, 
  // cosmology runs never output the units (why?), making data analyisis tough
  fprintf(fptr, "DensityUnits = %22.16e\n", DenUnits);
  fprintf(fptr, "LengthUnits = %22.16e\n",  LenUnits);
  fprintf(fptr, "TimeUnits = %22.16e\n",    TimeUnits);

  return SUCCESS;
}
#endif
