/***********************************************************************
/
/  EVOLVE LEVEL FUNCTION
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  February, 1995 by GB
/              Overhauled to make sure that all the subgrid's of a grid
/              advance with in lock step (i.e. with the same timestep and
/              in order).  This was done to allow a subgrid to get it's
/              boundary values from another subgrid (with the same parent).
/              Previously, a subgrid' BVs were always interpolated from its
/              parent.
/  modified2:  August, 1995 by GB
/                1) All grids on a level are processed at the same time
/                 (rather than all the subgrids of one parent).
/                2) C routines are called to loop over subgrids
/                 (so parallelizing C compilers can be used).
/                3) Subgrid timesteps are not constant over top grid step.
/              June, 1999 by GB -- Clean up somewhat
/
/  modified3:  August, 2001 by Alexei Kritsuk
/                Added 2nd call of PrepareDensityField() to compute
/                grav. potential (to be written with other baryon fields).
/  modified4:  January, 2004 by Alexei Kritsuk
/                Added support for RandomForcing
/  modified5:  February, 2006 by Daniel Reynolds
/                Added PotentialBdry to EvolveLevel and 
/                PrepareDensityField calls, so that it can be used 
/                within computing isolating BCs for self-gravity.
/  modified6:  January, 2007 by Robert Harkness
/                Group and in-core i/o
/  modified7:  May, 2008 by Robert Harkness
/                Remove Dan Reynolds' isograv mods
/
/  PURPOSE:
/    This routine is the main grid evolution function.  It assumes that the
/    grids of level-1 have already been advanced by dt (passed
/    in the argument) and that their boundary values are properly set.
/    We then perform a complete update on all grids on level, including:
/       - for each grid: set the boundary values from parent/subgrids
/       - for each grid: get a list of its subgrids
/       - determine the timestep based on the minimum timestep for all grids
/       - subcycle over the grid timestep and for each grid:
/           - copy the fields to the old fields
/           - solve the hydro equations (and save fluxes around subgrid)
/           - set the boundary values from parent and/or other grids
/           - update time and check dt(min) for that grid against dt(cycle)
/           - call EvolveLevel(level+1)
/           - accumulate flux around this grid
/       - correct the solution on this grid based on subgrid solutions
/       - correct the solution on this grid based on improved subgrid fluxes
/
/    This routine essentially solves (completely) the grids of this level
/       and all finer levels and then corrects the solution of
/       grids on this level based on the improved subgrid solutions/fluxes.
/
/    Note: as a convenience, we store the current grid's fluxes (i.e. the
/          fluxes around the exterior of this grid) as the last entry in
/          the list of subgrids.
/
************************************************************************/
 
#include <stdlib.h>
#include <stdio.h>
#include "performance.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
 
/* function prototypes */
 
void DeleteFluxes(fluxes *Fluxes);
int  RebuildHierarchy(TopGridData *MetaData,
		      LevelHierarchyEntry *LevelArray[], int level);
int  ReportMemoryUsage(char *header = NULL);
int  UpdateParticlePositions(grid *Grid);
int  CheckEnergyConservation(HierarchyEntry *Grids[], int grid,
			     int NumberOfGrids, int level, float dt);
int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);
int PrepareDensityField(LevelHierarchyEntry *LevelArray[],
			int level, TopGridData *MetaData, FLOAT When);
int SetBoundaryConditions(HierarchyEntry *Grids[], int NumberOfGrids,
			  int level, TopGridData *MetaData,
			  ExternalBoundary *Exterior, LevelHierarchyEntry * Level);
int UpdateFromFinerGrids(int level, HierarchyEntry *Grids[], int NumberOfGrids,
			 int NumberOfSubgrids[],
			 fluxes **SubgridFluxesEstimate[]);
int CommunicationUpdateStarParticleCount(HierarchyEntry *Grids[],
					 TopGridData *MetaData,
					 int NumberOfGrids);
int RadiationFieldUpdate(LevelHierarchyEntry *LevelArray[], int level,
			 TopGridData *MetaData);
int WriteStreamData(HierarchyEntry *Grids[], int NumberOfGrids, 
		    TopGridData *MetaData, int CycleCount, int EndStep = FALSE);
int WriteMovieData(char *basename, int filenumber,
		   LevelHierarchyEntry *LevelArray[], TopGridData *MetaData,
		   FLOAT WriteTime);
int WriteTracerParticleData(char *basename, int filenumber,
		   LevelHierarchyEntry *LevelArray[], TopGridData *MetaData,
		   FLOAT WriteTime);
int WriteAllData(char *basename, int filenumber, HierarchyEntry *TopGrid,
		 TopGridData &MetaData, ExternalBoundary *Exterior,
#ifdef TRANSFER
		 ImplicitProblemABC *ImplicitSolver,
#endif
		 FLOAT WriteTime = -1);
int Group_WriteAllData(char *basename, int filenumber, HierarchyEntry *TopGrid,
		 TopGridData &MetaData, ExternalBoundary *Exterior,
#ifdef TRANSFER
		 ImplicitProblemABC *ImplicitSolver,
#endif
		 FLOAT WriteTime = -1);

void my_exit(int status);
 
int ComputeRandomForcingNormalization(LevelHierarchyEntry *LevelArray[],
                                      int level, TopGridData *MetaData,
                                      float * norm, float * pTopGridTimeStep);
 
 
static int LevelCycleCount[MAX_DEPTH_OF_HIERARCHY];
static float norm = 0.0;            //AK
static float TopGridTimeStep = 0.0; //AK
 
 
/* EvolveGrid function */
 
int EvolveLevel(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[],
		int level, float dtLevelAbove, ExternalBoundary *Exterior
#ifdef TRANSFER
		, ImplicitProblemABC *ImplicitSolver
#endif
		)

#ifdef EMISSIVITY
/* reset Emissivity array here before next step calculate the new values */
  if (StarMakerEmissivityField > 0) {
  /* 
     clear the Emissivity of the level below, after the level below 
     updated the current level (it's parent) and before the next 
     timestep at the current level.
  */
    LevelHierarchyEntry *Temp;
    Temp = LevelArray[level+1];
    while (Temp != NULL) {
      Temp->GridData->UnigridClearEmissivity();
      Temp = Temp->NextGridThisLevel;
    }
  }
#endif

 
      SubgridFluxesEstimate[grid] = new fluxes *[NumberOfSubgrids[grid]];
 
      /* Collect the flux data and store it in the newly minted fluxes.
	 Or rather that's what we should do.  Instead, we create fluxes one
	 by one in this awkward array of pointers to pointers.  This should be
	 changed so that all the routines take arrays of flux rather than
	 arrays of pointers to flux.  Dumb. */
 
      counter = 0;
      NextGrid = Grids[grid]->NextGridNextLevel;
      while (NextGrid != NULL) {
	SubgridFluxesEstimate[grid][counter] = new fluxes;
	Grids[grid]->GridData->ComputeRefinementFactors
	                              (NextGrid->GridData, RefinementFactors);
	NextGrid->GridData->ReturnFluxDims
             (*(SubgridFluxesEstimate[grid][counter++]), RefinementFactors);
	NextGrid = NextGrid->NextGridThisLevel;
      }
 
      /* Add the external boundary of this subgrid to the subgrid list. This
	 makes it easy to keep adding up the fluxes of this grid, but we
	 must keep in mind that the last subgrid should be ignored elsewhere.*/
 
      SubgridFluxesEstimate[grid][counter] = new fluxes;
      Grids[grid]->GridData->ComputeRefinementFactors
                                   (Grids[grid]->GridData, RefinementFactors);
      Grids[grid]->GridData->ReturnFluxDims
               (*(SubgridFluxesEstimate[grid][counter]), RefinementFactors);
 
    } // end loop over grids (create Subgrid list)
 
#ifdef USE_LCAPERF
    LCAPERF_STOP("evolve-level-06"); // create subgrid list
#endif

    /* ------------------------------------------------------- */
    /* Prepare the density field (including particle density). */
 
//  fprintf(stderr, "%"ISYM": EvolveLevel: Enter PrepareDensityField\n", MyProcessorNumber);
 
#ifdef USE_LCAPERF
    LCAPERF_START("evolve-level-07"); // PrepareDensityField()
#endif

    When = 0.5;
 
    if (SelfGravity)
      if (PrepareDensityField(LevelArray, level, MetaData, When) == FAIL) {
	ENZO_FAIL("Error in PrepareDensityField.\n");
      }
 
//  fprintf(stderr, "%"ISYM": EvolveLevel: Exit PrepareDensityField\n", MyProcessorNumber);
 
    /* Prepare normalization for random forcing. Involves top grid only. */
 
    if (RandomForcing && MetaData->CycleNumber > 0 && level == 0)
      if ( ComputeRandomForcingNormalization(LevelArray, 0, MetaData,
                                             &norm, &TopGridTimeStep)
           == FAIL ) {
        ENZO_FAIL("Error in ComputeRandomForcingNormalization.\n");
      }
 
#ifdef USE_LCAPERF
    LCAPERF_STOP("evolve-level-07"); // PrepareDensityField()
#endif

    /* ------------------------------------------------------- */
    /* Evolve all grids by timestep dtThisLevel. */
 
    for (grid = 0; grid < NumberOfGrids; grid++) {
 
      /* Call analysis routines. */
 
#ifdef USE_LCAPERF
      LCAPERF_START_LOW("evolve-level-08"); // Call analysis routines
#endif

      if (ProblemType == 24)
	Grids[grid]->GridData->SphericalInfallGetProfile(level, 1);
      if (ProblemType == 30)
	Grids[grid]->GridData->AnalyzeTrackPeaks(level, 0);
      if (ProblemType == 27)
	if (Grids[grid]->GridData->ReturnProcessorNumber()==MyProcessorNumber){
	  float AM[3], MeanVelocity[3], DMVelocity[3];
	  FLOAT Center[] = {0,0,0}, CenterOfMass[3], DMCofM[3];
	  Grids[grid]->GridData->CalculateAngularMomentum(Center, AM,
			   MeanVelocity, DMVelocity, CenterOfMass, DMCofM);
	  fprintf(stdout, "level = %"ISYM" %"ISYM" %"ISYM"  Vel %"FSYM" %"FSYM" %"FSYM"  DMVel %"FSYM" %"FSYM" %"FSYM"  CofM %"PSYM" %"PSYM" %"PSYM"  DMCofM %"FSYM" %"FSYM" %"FSYM"\n",
		level, LevelCycleCount[level], grid, MeanVelocity[0],
		MeanVelocity[1], MeanVelocity[2],
		DMVelocity[0], DMVelocity[1], DMVelocity[2],
		-CenterOfMass[0], -CenterOfMass[1], -CenterOfMass[2],
		DMCofM[0], DMCofM[1], DMCofM[2]);
	}
 
#ifdef USE_LCAPERF
      LCAPERF_STOP_LOW("evolve-level-08"); // Call analysis routines
#endif

      /* Gravity: compute acceleration field for grid and particles. */
 
#ifdef USE_LCAPERF
      LCAPERF_START("evolve-level-09"); // Compute self-gravity acceleration
#endif
    
      if (SelfGravity) {
	int Dummy;
	if (level <= MaximumGravityRefinementLevel) {
 
	  /* Compute the potential. */
 
	  if (level > 0)
	    if (Grids[grid]->GridData->SolveForPotential(Dummy, level)
		== FAIL) {
	      ENZO_FAIL("Error in grid->SolveForPotential.\n");
	    }
	  if (Grids[grid]->GridData->ComputeAccelerations(level) == FAIL) {
	    ENZO_FAIL("Error in grid->ComputeAccelerations.\n");
	  }
	}
	  /* otherwise, interpolate potential from coarser grid, which is
	     now done in PrepareDensity. */
 
      } // end: if (SelfGravity)
 
#ifdef USE_LCAPERF
      LCAPERF_STOP("evolve-level-09"); // Compute self-gravity acceleration
#endif

      /* Gravity: compute field due to preset sources. */
 
#ifdef USE_LCAPERF
      LCAPERF_START_LOW("evolve-level-10"); // ComputeAccelerationFieldExternal()
#endif

      if (UniformGravity || PointSourceGravity || DiskGravity )
	if (Grids[grid]->GridData->ComputeAccelerationFieldExternal() ==FAIL) {
	  ENZO_FAIL("Error in grid->ComputeAccelerationFieldExternal.\n");
	}
 
#ifdef USE_LCAPERF
      LCAPERF_STOP_LOW("evolve-level-10"); // ComputeAccelerationFieldExternal()
#endif

      /* Check for energy conservation. */
/*
      if (ComputePotential)
	if (CheckEnergyConservation(Grids, grid, NumberOfGrids, level,
				    dtThisLevel) == FAIL) {
	  ENZO_FAIL("Error in CheckEnergyConservation.\n");
	}
*/
      /* Copy current fields (with their boundaries) to the old fields
	  in preparation for the new step. */
 
#ifdef USE_LCAPERF
      LCAPERF_START("evolve-level-11"); // CopyBaryonFieldToOldBaryonField()
#endif

      if (Grids[grid]->GridData->CopyBaryonFieldToOldBaryonField() == FAIL) {
	ENZO_FAIL("Error in grid->CopyBaryonFieldToOldBaryonField.\n");
      }
 
#ifdef USE_LCAPERF
      LCAPERF_STOP("evolve-level-11"); // CopyBaryonFieldToOldBaryonField()
#endif

      /* Add RandomForcing fields to velocities after the copying of current
         fields to old. I also update the total energy accordingly here.
         It makes no sense to force on the very first time step. */
 
#ifdef USE_LCAPERF
      LCAPERF_START_LOW("evolve-level-12"); // AddRandomForcing()
#endif

      if (RandomForcing && MetaData->CycleNumber > 0) //AK
        if(Grids[grid]->GridData->AddRandomForcing(&norm,
                                                   TopGridTimeStep) == FAIL)
          fprintf(stderr, "Error in AddRandomForcing.\n");
 
#ifdef USE_LCAPERF
      LCAPERF_STOP_LOW("evolve-level-12"); // AddRandomForcing()
#endif

      /* Call hydro solver and save fluxes around subgrids. */
 
//      fprintf(stderr, "%"ISYM": Calling Hydro\n", MyProcessorNumber);
 
#ifdef USE_LCAPERF
      LCAPERF_START("evolve-level-13"); // SolveHydroEquations()
#endif

      if (Grids[grid]->GridData->SolveHydroEquations(LevelCycleCount[level],
	 NumberOfSubgrids[grid], SubgridFluxesEstimate[grid], level) == FAIL) {
	ENZO_FAIL("Error in grid->SolveHydroEquations.\n");
      }
 
#ifdef USE_LCAPERF
      LCAPERF_STOP("evolve-level-13"); // SolveHydroEquations()
#endif

//      fprintf(stderr, "%"ISYM": Called Hydro\n", MyProcessorNumber);
 
      /* Solve the cooling and species rate equations. */
 
//      fprintf(stderr, "%"ISYM": Calling SolveCoolAndRateEquations\n", MyProcessorNumber);

      if (MultiSpecies && RadiativeCooling && GadgetEquilibriumCooling == 0) {
 
	LCAPERF_START("evolve-level-14"); // change this?

	int RTCoupledSolverIntermediateStep = FALSE;
	if (Grids[grid1]->GridData->SolveRateAndCoolEquations(RTCoupledSolverIntermediateStep) == FAIL) {
	  ENZO_FAIL("Error in grid->SolveRateEquations.\n");
	}
 
	LCAPERF_STOP("evolve-level-14"); // change this?

//      fprintf(stderr, "%"ISYM": Called SolveCoolAndRateEquations\n", MyProcessorNumber);

      } else {

//      fprintf(stderr, "%"ISYM": Calling MultiSpecies\n", MyProcessorNumber);
 
	LCAPERF_START("evolve-level-14"); // SolveRateEquations()

	if (MultiSpecies)
	  if (Grids[grid1]->GridData->SolveRateEquations() == FAIL) {
	    ENZO_FAIL("Error in grid->SolveRateEquations.\n");
	  }
 
	LCAPERF_STOP("evolve-level-14"); // SolveRateEquations()

//      fprintf(stderr, "%"ISYM": Called MultiSpecies\n", MyProcessorNumber);
 
	/* Include radiative cooling/heating. */
 
//      fprintf(stderr, "%"ISYM": Calling RadiativeCooling\n", MyProcessorNumber);
 
	LCAPERF_START("evolve-level-15"); // SolveRadiativeCooling()

	if (RadiativeCooling)
	  if (Grids[grid1]->GridData->SolveRadiativeCooling() == FAIL) {
	    ENZO_FAIL("Error in grid->SolveRadiativeCooling.\n");
	  }
 
	LCAPERF_STOP("evolve-level-15"); // SolveRadiativeCooling()

//      fprintf(stderr, "%"ISYM": Called RadiativeCooling\n", MyProcessorNumber);

      }
 
      /* Update particle positions (if present). */
 
//      fprintf(stderr, "%"ISYM": Calling UpdatePP\n", MyProcessorNumber);
 
#ifdef USE_LCAPERF
      LCAPERF_START("evolve-level-16"); // UpdateParticlePositions()
#endif

      if (UpdateParticlePositions(Grids[grid]->GridData) == FAIL) {
	ENZO_FAIL("Error in UpdateParticlePositions.\n");
      }
 
#ifdef USE_LCAPERF
      LCAPERF_STOP("evolve-level-16"); // UpdateParticlePositions()
#endif

//      fprintf(stderr, "%"ISYM": Called UpdatePP\n", MyProcessorNumber);
 
      /* Include 'star' particle creation and feedback.
         (first, set the under_subgrid field). */
 
#ifdef USE_LCAPERF
      LCAPERF_START_LOW("evolve-level-17"); // star particle creation/feedback
#endif

      if (StarParticleCreation || StarParticleFeedback) {
	Grids[grid]->GridData->ZeroSolutionUnderSubgrid(NULL,
						 ZERO_UNDER_SUBGRID_FIELD);
	LevelHierarchyEntry *Temp2 = LevelArray[level+1];
	while (Temp2 != NULL) {
	  Grids[grid]->GridData->ZeroSolutionUnderSubgrid(Temp2->GridData,
					 ZERO_UNDER_SUBGRID_FIELD);
	  Temp2 = Temp2->NextGridThisLevel;
	}
      }
      if (StarParticleCreation || StarParticleFeedback) {
	if (Grids[grid]->GridData->StarParticleHandler(level
#ifdef EMISSIVITY
	/* adding the changed StarParticleHandler prototype */
						       ,dtLevelAbove
#endif
           ) == FAIL) {
	  ENZO_FAIL("Error in grid->StarParticleWrapper");
	}
      }
 
#ifdef USE_LCAPERF
      LCAPERF_STOP_LOW("evolve-level-17"); // star particle creation/feedback
#endif

      /* Gravity: clean up AccelerationField. */
 
#ifdef USE_LCAPERF
      LCAPERF_START_LOW("evolve-level-18"); // clean up AccelerationField
#endif

      if (SelfGravity || UniformGravity || PointSourceGravity || DiskGravity ) {
	if (level != MaximumGravityRefinementLevel ||
	    MaximumGravityRefinementLevel == MaximumRefinementLevel)
	  Grids[grid]->GridData->DeleteAccelerationField();
	Grids[grid]->GridData->DeleteParticleAcceleration();
      }
 
#ifdef USE_LCAPERF
      LCAPERF_STOP_LOW("evolve-level-18"); // clean up AccelerationField
#endif

      /* Update current problem time of this subgrid. */
 
#ifdef USE_LCAPERF
      LCAPERF_START_LOW("evolve-level-19"); // SetTimeNextTimestep()
#endif

      Grids[grid]->GridData->SetTimeNextTimestep();
 
#ifdef USE_LCAPERF
      LCAPERF_STOP_LOW("evolve-level-19"); // SetTimeNextTimestep()
#endif

      /* If using comoving co-ordinates, do the expansion terms now. */
 
#ifdef USE_LCAPERF
      LCAPERF_START("evolve-level-20"); // ComovingExpansionTerms()
#endif

      if (ComovingCoordinates)
	Grids[grid]->GridData->ComovingExpansionTerms();
 
#ifdef USE_LCAPERF
      LCAPERF_STOP("evolve-level-20"); // ComovingExpansionTerms()
#endif

    }  // end loop over grids
 
    /* For each grid: a) interpolate boundaries from the parent grid.
                      b) copy any overlapping zones from siblings. */
 
#ifdef USE_LCAPERF
    LCAPERF_START("evolve-level-21"); // SetBoundaryConditions()
#endif

if (SetBoundaryConditions(Grids, NumberOfGrids, level, MetaData,
			  Exterior, LevelArray[level]) == FAIL)
  ENZO_FAIL("Error in SetBoundaryConditions()!\n");
  

#ifdef USE_LCAPERF
    LCAPERF_STOP("evolve-level-21"); // SetBoundaryConditions()
#endif

    /* Update the star particle counters. */
 
#ifdef USE_LCAPERF
    LCAPERF_START("evolve-level-22"); // CommunicationUpdateStarParticleCount()
#endif

    if (StarParticleCreation)
      if (CommunicationUpdateStarParticleCount(Grids, MetaData,
					       NumberOfGrids) == FAIL)
	ENZO_FAIL("Error in CommunicationUpdateStarParticleCount()!\n");
  
#ifdef USE_LCAPERF
    LCAPERF_STOP("evolve-level-22"); // CommunicationUpdateStarParticleCount()
#endif

    /* Check for tracer particle output (only if this bottom of hierarchy). */
 
#ifdef USE_LCAPERF
    LCAPERF_START("evolve-level-24"); // WriteTracerParticleData()
#endif

    if (LevelArray[level+1] == NULL)
      if (LevelArray[level]->GridData->ReturnTime() >=
	  MetaData->TimeLastTracerParticleDump +
	  MetaData->dtTracerParticleDump &&
	  MetaData->dtTracerParticleDump > 0.0) {
	MetaData->TimeLastTracerParticleDump += MetaData->dtTracerParticleDump;
	if (WriteTracerParticleData(MetaData->TracerParticleDumpName,
				    MetaData->TracerParticleDumpNumber++,
				    LevelArray, MetaData,
				    LevelArray[level]->GridData->ReturnTime()) == FAIL) {
	  ENZO_FAIL("Error in WriteTracerParticleData.\n");
	}
      }
 
#ifdef USE_LCAPERF
    LCAPERF_STOP("evolve-level-24"); // WriteTracerParticleData()
#endif

    /* If cosmology, then compute grav. potential for output if needed. */
 
#ifdef USE_LCAPERF
    LCAPERF_START("evolve-level-25"); // PrepareDensityField()
#endif

    if (ComovingCoordinates && SelfGravity && WritePotential) {
      CopyGravPotential = TRUE;
      When = 0.0;
      if (PrepareDensityField(LevelArray, level, MetaData, When) == FAIL) {
        ENZO_FAIL("Error in PrepareDensityField.\n");
      }
      CopyGravPotential = FALSE;
 
      for (grid = 0; grid < NumberOfGrids; grid++) {
        int Dummy;
        if (level <= MaximumGravityRefinementLevel) {
 
          /* Compute the potential. */
 
          if (level > 0)
            if (Grids[grid]->GridData->SolveForPotential(Dummy, level)
                == FAIL) {
              ENZO_FAIL("Error in grid->SolveForPotential.\n");
            }
          // fprintf(stderr, "Call CP from EvolveLevel\n");
          Grids[grid]->GridData->CopyPotentialToBaryonField();
        }
        /* otherwise output empty potential field. */
 
      } //  end loop over grids
    } // if WritePotential
 
#ifdef USE_LCAPERF
    LCAPERF_STOP("evolve-level-25"); // PrepareDensityField()
#endif
 
    /* Check for new level output (only if this is bottom of hierarchy). */
 
#ifdef USE_LCAPERF
    LCAPERF_START("evolve-level-26"); // WriteAllData()
#endif

    if (MetaData->OutputFirstTimeAtLevel > 0 &&
	level >= MetaData->OutputFirstTimeAtLevel &&
	LevelArray[level+1] == NULL) {

      MetaData->OutputFirstTimeAtLevel = level+1;
      LevelHierarchyEntry *Temp2 = LevelArray[0];

      while (Temp2->NextGridThisLevel != NULL)
	Temp2 = Temp2->NextGridThisLevel; /* ugh: find last in linked list */

      //#ifdef USE_HDF5_GROUPS
      if (Group_WriteAllData(MetaData->DataDumpName, MetaData->DataDumpNumber++,
		       Temp2->GridHierarchyEntry, *MetaData, Exterior,
#ifdef TRANSFER
		       ImplicitSolver,
#endif
		       LevelArray[level]->GridData->ReturnTime()) == FAIL) {
	ENZO_FAIL("Error in Group_WriteAllData.\n");
      }
// #else
//       if (WriteAllData(MetaData->DataDumpName, MetaData->DataDumpNumber++,
// 		       Temp2->GridHierarchyEntry, *MetaData, Exterior, 
// #ifdef TRANSFER
// 		       ImplicitSolver,
// #endif
// 		       LevelArray[level]->GridData->ReturnTime()) == FAIL) {
// 	ENZO_FAIL("Error in WriteAllData.\n");
//       }
// #endif
    }
 
#ifdef USE_LCAPERF
    LCAPERF_STOP("evolve-level-26"); // WriteAllData()
#endif

    /* Check for stop (unpleasant to exit from here, but...). */
 
    if (MetaData->StopFirstTimeAtLevel > 0 &&
	level >= MetaData->StopFirstTimeAtLevel &&
	LevelArray[level+1] == NULL) {


      /* Write movie data in all grids if necessary */

      if (MovieSkipTimestep != INT_UNDEFINED) {
	for (int mlevel = 0; mlevel < MAX_DEPTH_OF_HIERARCHY; mlevel++) {
	  if (LevelArray[mlevel] == NULL) break;
	  delete [] Grids;
	  NumberOfGrids = GenerateGridArray(LevelArray, mlevel, &Grids);
	  if (WriteStreamData(Grids, NumberOfGrids, MetaData,
			      LevelCycleCount[mlevel], TRUE) == FAIL) {
	    ENZO_FAIL("Error in WriteStreamData.\n");
	  }
	}
      }

      fprintf(stderr, "Stopping due to request on level %"ISYM"\n", level);
      my_exit(EXIT_SUCCESS);
    }
 
    /* For each grid, delete the GravitatingMassFieldParticles. */
 
#ifdef USE_LCAPERF
    LCAPERF_START("evolve-level-27"); // DeleteGravitatingMassFieldParticles()
#endif

    for (grid = 0; grid < NumberOfGrids; grid++)
      Grids[grid]->GridData->DeleteGravitatingMassFieldParticles();
 
#ifdef USE_LCAPERF
    LCAPERF_STOP("evolve-level-27"); // DeleteGravitatingMassFieldParticles()
#endif

    /* ----------------------------------------- */
    /* Evolve the next level down (recursively). */
 
    if (LevelArray[level+1] != NULL) {
      if (EvolveLevel(MetaData, LevelArray, level+1, dtThisLevel, Exterior) == FAIL) {
	ENZO_VFAIL("Error in EvolveLevel (%"ISYM").\n", level)
      }
    }

    // Streaming movie output (only run if everything is evolved) */

    if (MovieSkipTimestep != INT_UNDEFINED) {
      if (WriteStreamData(Grids, NumberOfGrids, MetaData, 
			  LevelCycleCount[level]) == FAIL) {
	ENZO_FAIL("Error in WriteStreamData.\n");
      }
    }
 

#if defined(USE_LCAPERF) && defined(LCAPERF_LEVELS)
    lcaperf.attribute ("level",&lcaperf_level,LCAPERF_INT);
#endif

    /* ------------------------------------------------------- */
    /* For each grid,
     (a) project the subgrid's solution into this grid (step #18)
     (b) correct for the difference between this grid's fluxes and the
         subgrid's fluxes. (step #19) */
 
#ifdef USE_LCAPERF
    LCAPERF_START("evolve-level-28"); // UpdateFromFinerGrids()
#endif

    if (UpdateFromFinerGrids(level, Grids, NumberOfGrids, NumberOfSubgrids,
			     SubgridFluxesEstimate) == FAIL)
      ENZO_FAIL("Error in UpdateFromFinerGrids()!\n");
 
#ifdef USE_LCAPERF
    LCAPERF_STOP("evolve-level-28"); // UpdateFromFinerGrids()
#endif

    /* ------------------------------------------------------- */
    /* Add the saved fluxes (in the last subsubgrid entry) to the exterior
       fluxes for this subgrid .
       (Note: this must be done after CorrectForRefinedFluxes). */
 
#ifdef USE_LCAPERF
    LCAPERF_START("evolve-level-29"); // AddToBoundaryFluxes()
#endif

    for (grid = 0; grid < NumberOfGrids; grid++) {
 
      if (FluxCorrection)
	if (Grids[grid]->GridData->AddToBoundaryFluxes
	    (SubgridFluxesEstimate[grid][NumberOfSubgrids[grid] - 1])
	    == FAIL) {
	  ENZO_FAIL("Error in grid->AddToBoundaryFluxes.\n");
	}
 
      /* Delete fluxes pointed to by SubgridFluxesEstimate[subgrid]. */
 
      for (subgrid = 0; subgrid < NumberOfSubgrids[grid]; subgrid++) {
	DeleteFluxes(SubgridFluxesEstimate[grid][subgrid]);
	delete       SubgridFluxesEstimate[grid][subgrid];
      }
      delete [] SubgridFluxesEstimate[grid];
 
    } // end of loop over grids
 
#ifdef USE_LCAPERF
    LCAPERF_STOP("evolve-level-29"); // AddToBoundaryFluxes()
#endif

    /* Recompute radiation field, if requested. */
 
#ifdef USE_LCAPERF
    LCAPERF_START("evolve-level-30"); // RadiationFieldUpdate()
#endif

    if (RadiationFieldType >= 10 && RadiationFieldType <= 11 &&
	level <= RadiationFieldLevelRecompute)
      if (RadiationFieldUpdate(LevelArray, level, MetaData) == FAIL) {
	ENZO_FAIL("Error in RecomputeRadiationField.\n");
      }
 
#ifdef USE_LCAPERF
    LCAPERF_STOP("evolve-level-30"); // RadiationFieldUpdate()
#endif

    /* Rebuild the Grids on the next level down.
       Don't bother on the last cycle, as we'll rebuild this grid soon. */
 
#ifdef USE_LCAPERF
    LCAPERF_START("evolve-level-31"); // RebuildHierarchy()
#endif

    if (dtThisLevelSoFar < dtLevelAbove)
      if (RebuildHierarchy(MetaData, LevelArray, level) == FAIL) {
	ENZO_FAIL("Error in RebuildHierarchy.\n");
      }
 
#ifdef USE_LCAPERF
    LCAPERF_STOP("evolve-level-31"); // RebuildHierarchy()
#endif

    cycle++;
    LevelCycleCount[level]++;
 
  } // end of loop over subcycles
 
  if (debug)
    printf("EvolveLevel[%"ISYM"]: NumberOfSubCycles = %"ISYM" (%"ISYM" total)\n", level,
           cycle, LevelCycleCount[level]);
 
  /* If possible & desired, report on memory usage. */
 
  //  if (debug)
  ReportMemoryUsage("Memory usage report: Evolve Level");
 
#if defined(USE_LCAPERF) && defined(LCAPERF_LEVELS)
  lcaperf.attribute ("level",0,LCAPERF_NULL);
#endif

  /* Clean up. */
 
#ifdef UNUSED
  if (level > MaximumGravityRefinementLevel &&

      level == MaximumRefinementLevel)
    ZEUSQuadraticArtificialViscosity /= 1;
#endif
 
  delete [] NumberOfSubgrids;
  delete [] Grids;
  delete [] SubgridFluxesEstimate;
 
  return SUCCESS;
 
}
