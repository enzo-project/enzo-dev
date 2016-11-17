/***********************************************************************
/
/  ADD FEEDBACK TO RADIAL PROFILE OVER MULTIPLE GRIDS
/
/  written by: John Wise
/  date:       September, 2005
/  modified1: Ji-hoon Kim
/             October, 2009
/
/ PURPOSE: To apply feedback effects, we must consider multiple grids
/          since sometimes the feedback radius often exceeds the grid
/          boundaries.  This routine makes sure that all of the grids
/          have the same code time to ensure consistency.
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "performance.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "list.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"

#define MAX_TEMPERATURE 1e8

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int RecalibrateMBHFeedbackThermalRadius(FLOAT star_pos[], LevelHierarchyEntry *LevelArray[], 
					int level, float &Radius, 
					double &EjectaDensity, double &EjectaMetalDensity,
					double &EjectaThermalEnergy);
int RemoveParticles(LevelHierarchyEntry *LevelArray[], int level, int ID);
FLOAT FindCrossSection(int type, float energy);

int StarParticleAddFeedback(TopGridData *MetaData, 
			    LevelHierarchyEntry *LevelArray[], int level, 
			    Star* &AllStars, bool* &AddedFeedback)
{

  const double pc = 3.086e18, Msun = 1.989e33, pMass = 1.673e-24, 
    gravConst = 6.673e-8, yr = 3.1557e7, Myr = 3.1557e13;

  Star *cstar;
  bool MarkedSubgrids = false;
  bool SphereCheck;
  int i, l, dim, temp_int, SkipMassRemoval, SphereContained,
      SphereContainedNextLevel, dummy, count;
  float influenceRadius, RootCellWidth, SNe_dt, dtForThisStar, MassLoss;
  double EjectaThermalEnergy, EjectaDensity, EjectaMetalDensity;
  FLOAT Time;
  LevelHierarchyEntry *Temp;

  if (AllStars == NULL)
    return SUCCESS;

  LCAPERF_START("StarParticleAddFeedback");

  /* Get time and SNe timestep */

  Temp = LevelArray[level];
  Time = Temp->GridData->ReturnTime();
  if (LastSupernovaTime < 0)
    SNe_dt = 0.0;
  else
    SNe_dt = Time - LastSupernovaTime;
  LastSupernovaTime = Time;
  RootCellWidth = 1.0 / MetaData->TopGridDims[0];

  /* Set the units. */

  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, 
    VelocityUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, Time);

  count = 0;
  // clear list of Supernovae at each timestep to avoid adding duplicates in Grid_AddFeedbackSphere                                     
  if(UseSupernovaSeedFieldSourceTerms){
    LevelArray[level]->GridData->SuperNovaList.clear();
  }
  for (cstar = AllStars; cstar; cstar = cstar->NextStar, count++) {

    AddedFeedback[count] = false;

    /* Special case for "normal" star particles to account for mass
       loss through supernovae. */

    if (cstar->ReturnType() == NormalStar &&
	cstar->ReturnLevel() == level) {
      MassLoss = cstar->CalculateMassLoss(SNe_dt);
      cstar->SetAccretionMass(-MassLoss);
    }

    if ((cstar->ReturnFeedbackFlag() != MBH_THERMAL) &&
	(cstar->ReturnFeedbackFlag() != MBH_JETS) &&
	!cstar->ApplyFeedbackTrue(SNe_dt))
      continue;

    dtForThisStar = LevelArray[level]->GridData->ReturnTimeStep();
	  
    /* Compute some parameters */

    cstar->CalculateFeedbackParameters
      (influenceRadius, RootCellWidth, SNe_dt, EjectaDensity, 
       EjectaThermalEnergy, EjectaMetalDensity, DensityUnits, LengthUnits, 
       TemperatureUnits, TimeUnits, VelocityUnits, dtForThisStar, 
       Time, SphereCheck);

    if (SphereCheck) {

    /* Recalibrate MBHFeedbackThermalRadius if requested */

    if (cstar->ReturnFeedbackFlag() == MBH_THERMAL)    
      RecalibrateMBHFeedbackThermalRadius(cstar->ReturnPosition(), LevelArray, level, influenceRadius, 
					  EjectaDensity, EjectaMetalDensity, EjectaThermalEnergy);

    /* Determine if a sphere with enough mass (or equivalently radius
       for SNe) is enclosed within grids on this level */

    LCAPERF_START("star_FindFeedbackSphere");
    cstar->FindFeedbackSphere
      (LevelArray, level, influenceRadius, EjectaDensity, EjectaThermalEnergy, 
       SphereContained, SkipMassRemoval, DensityUnits, LengthUnits, 
       TemperatureUnits, TimeUnits, VelocityUnits, Time, MarkedSubgrids);
    LCAPERF_STOP("star_FindFeedbackSphere");

    /* If the particle already had sufficient mass, we still want to
       mark this particle to activate it. */

    if (SkipMassRemoval == TRUE)
      AddedFeedback[count] = true;

    /* If there's no feedback or something weird happens, don't bother. */

    if ( influenceRadius <= tiny_number || 
	 SphereContained == FALSE ||
	((cstar->ReturnFeedbackFlag() == MBH_THERMAL ||
	  cstar->ReturnFeedbackFlag() == MBH_JETS) &&
	 (influenceRadius >= RootCellWidth/2 || 
	  EjectaThermalEnergy <= tiny_number)) )
      continue;

    /* Determine if a sphere is enclosed within the grids on next level
       If that is the case, we perform AddFeedbackSphere not here, 
       but in the EvolveLevel of the next level. */

    SphereContainedNextLevel = FALSE;

    LCAPERF_START("star_FindFeedbackSphere2");
    if ((cstar->ReturnFeedbackFlag() == MBH_THERMAL ||
	 cstar->ReturnFeedbackFlag() == MBH_JETS ||
	 cstar->ReturnFeedbackFlag() == CONT_SUPERNOVA) &&
	LevelArray[level+1] != NULL)
      cstar->FindFeedbackSphere
	(LevelArray, level+1, influenceRadius, EjectaDensity, EjectaThermalEnergy, 
	 SphereContainedNextLevel, dummy, DensityUnits, LengthUnits, 
	 TemperatureUnits, TimeUnits, VelocityUnits, Time, MarkedSubgrids);
    LCAPERF_STOP("star_FindFeedbackSphere2");

//    if (debug) {
//      fprintf(stdout, "EjectaDensity=%g, influenceRadius=%g\n", EjectaDensity, influenceRadius); 
//      fprintf(stdout, "SkipMassRemoval=%d, SphereContained=%d, SphereContainedNextLevel=%d\n", 
//	      SkipMassRemoval, SphereContained, SphereContainedNextLevel); 
//    }

    /* Quit this routine when 
       (1) sphere is not contained, or 
       (2) sphere is contained, but the next level can contain the sphere, too. */ 
    if ((SphereContained == FALSE) ||
	(SphereContained == TRUE && SphereContainedNextLevel == TRUE))
      continue;

    } // ENDIF SphereCheck
    else {
      
      /* When the sphere is completely confined in a grid, only apply
	 feedback at the level at which the star exists. */

      if (level != cstar->ReturnLevel()) 
	continue;

    }
    
    /* Now set cells within the radius to their values after feedback.
       While walking through the hierarchy, look for particle to
       change their properties to post-feedback values. */

    int CellsModified = 0;

    if (SkipMassRemoval == FALSE) {

      /* Determine the H-ionizing photon luminosity to calculate the
	 photo-ionization and heating rate in the initial Stroemgren
	 sphere. */

      int nbins;
      double Q[MAX_ENERGY_BINS], Q_HI, sigma;
      float energies[MAX_ENERGY_BINS], deltaE;
#ifdef TRANSFER
      if (RadiativeTransfer) {
	cstar->ComputePhotonRates(TimeUnits, nbins, energies, Q);
	sigma = (double) FindCrossSection(0, energies[0]);  // HI (cm^2)
	Q_HI = Q[0];
	deltaE = energies[0] - 13.6;  // eV
      } else
#endif /* TRANSFER */
	{
	Q_HI = 0.0;
	sigma = 0.0;
	deltaE = 0.0;
      }

      for (l = level; l < MAX_DEPTH_OF_HIERARCHY; l++)
	for (Temp = LevelArray[l]; Temp; Temp = Temp->NextGridThisLevel) 
	  Temp->GridData->AddFeedbackSphere
	    (cstar, l, influenceRadius, DensityUnits, LengthUnits, 
	     VelocityUnits, TemperatureUnits, TimeUnits, EjectaDensity, 
	     EjectaMetalDensity, EjectaThermalEnergy, Q_HI, sigma, deltaE, 
	     CellsModified);
    } // ENDIF

//    fprintf(stdout, "StarParticleAddFeedback[%"ISYM"][%"ISYM"]: "
//	    "Radius = %e pc, changed %"ISYM" cells.\n", 
//	    cstar->ReturnID(), level, influenceRadius*LengthUnits/pc, CellsModified); 

    /* Remove mass from the star that is added to grids. Also, because EjectaDensity 
       is added with zero net momentum, increase the particle's velocity accordingly. 
       Only for MBH_JETS; currently this is done in Grid_AddFeedbackSphere.C */

    /*
    if (EjectaDensity != 0 && CellsModified > 0)
      if (cstar->ReturnFeedbackFlag() == MBH_THERMAL ||
	  cstar->ReturnFeedbackFlag() == MBH_JETS)
	cstar->RemoveMassFromStarAfterFeedback(influenceRadius, EjectaDensity, 
					       DensityUnits, LengthUnits, CellsModified);
    */

    /* Only kill a Pop III star after it has gone SN */

    if (cstar->ReturnFeedbackFlag() == SUPERNOVA)
      cstar->SetFeedbackFlag(DEATH);

    /* We only color the fields once */

    AddedFeedback[count] = true;

#ifdef UNUSED
    temp_int = CellsModified;
    MPI_Reduce(&temp_int, &CellsModified, 1, MPI_INT, MPI_SUM, ROOT_PROCESSOR,
	       MPI_COMM_WORLD);

    if (debug) {
      if (cstar->ReturnFeedbackFlag() != FORMATION)
	fprintf(stdout, "StarParticleAddFeedback[%"ISYM"][%"ISYM"]: "
		"Radius = %"GSYM" pc\n",
		cstar->ReturnID(), level, influenceRadius*LengthUnits/pc);
      if (cstar->ReturnFeedbackFlag() == DEATH || 

	  cstar->ReturnFeedbackFlag() == CONT_SUPERNOVA ||
	  cstar->ReturnFeedbackFlag() == MBH_THERMAL ||
	  cstar->ReturnFeedbackFlag() == MBH_JETS )
	fprintf(stdout, "StarParticleAddFeedback[%"ISYM"][%"ISYM"]: "
		"Energy = %"GSYM"  , skip = %"ISYM"\n",
		cstar->ReturnID(), level, EjectaThermalEnergy, SkipMassRemoval);
      fprintf(stdout, "StarParticleAddFeedback[%"ISYM"][%"ISYM"]: "
	      "changed %"ISYM" cells.  AddedFeedback[%d] = %d\n", 
	      cstar->ReturnID(), level, CellsModified, 
	      count, AddedFeedback[count]);
    }
#endif
    
  } // ENDFOR stars

  LCAPERF_STOP("StarParticleAddFeedback");
  return SUCCESS;

}
