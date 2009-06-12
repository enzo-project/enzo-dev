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
/  modified7:  December, 2007 by Robert Harkness
/                Optional StaticSiblingList for root grid
/  modified8:  April, 2009 by John Wise
/                Added star particle class and radiative transfer
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
 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

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
#include "CommunicationUtilities.h"
 
/* function prototypes */
 

int  RebuildHierarchy(TopGridData *MetaData,
		      LevelHierarchyEntry *LevelArray[], int level);
int  ReportMemoryUsage(char *header = NULL);
int  UpdateParticlePositions(grid *Grid);
int  CheckEnergyConservation(HierarchyEntry *Grids[], int grid,
			     int NumberOfGrids, int level, float dt);
int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);
 
#ifdef FAST_SIB
int PrepareDensityField(LevelHierarchyEntry *LevelArray[],
			SiblingGridList SiblingList[],
			int level, TopGridData *MetaData, FLOAT When);
#else  // !FAST_SIB
int PrepareDensityField(LevelHierarchyEntry *LevelArray[],
                        int level, TopGridData *MetaData, FLOAT When);
#endif  // end FAST_SIB
 
#ifdef FAST_SIB
int SetBoundaryConditions(HierarchyEntry *Grids[], int NumberOfGrids,
			  SiblingGridList SiblingList[],
			  int level, TopGridData *MetaData,
			  ExternalBoundary *Exterior, LevelHierarchyEntry * Level);
#else
int SetBoundaryConditions(HierarchyEntry *Grids[], int NumberOfGrids,
                          int level, TopGridData *MetaData,
                          ExternalBoundary *Exterior, LevelHierarchyEntry * Level);
#endif

#ifdef SAB
#ifdef FAST_SIB
int SetAccelerationBoundary(HierarchyEntry *Grids[], int NumberOfGrids,
			    SiblingGridList SiblingList[],
			    int level, TopGridData *MetaData,
			    ExternalBoundary *Exterior,
			    LevelHierarchyEntry * Level,
			    int CycleNumber);
#else
int SetAccelerationBoundary(HierarchyEntry *Grids[], int NumberOfGrids,
			    int level, TopGridData *MetaData, 
			    ExternalBoundary *Exterior,
			    LevelHierarchyEntry * Level,
			    int CycleNumber);
#endif
#endif

#ifdef FLUX_FIX
int UpdateFromFinerGrids(int level, HierarchyEntry *Grids[], int NumberOfGrids,
			 int NumberOfSubgrids[],
			 fluxes **SubgridFluxesEstimate[],
			 LevelHierarchyEntry *SUBlingList[],
			 TopGridData *MetaData);
#else
int UpdateFromFinerGrids(int level, HierarchyEntry *Grids[], int NumberOfGrids,
			 int NumberOfSubgrids[],
			 fluxes **SubgridFluxesEstimate[]);
#endif
int CreateFluxes(HierarchyEntry *Grids[],fluxes **SubgridFluxesEstimate[],
		 int NumberOfGrids,int NumberOfSubgrids[]);		 
int FinalizeFluxes(HierarchyEntry *Grids[],fluxes **SubgridFluxesEstimate[],
		 int NumberOfGrids,int NumberOfSubgrids[]);		 
int RadiationFieldUpdate(LevelHierarchyEntry *LevelArray[], int level,
			 TopGridData *MetaData);


int OutputFromEvolveLevel(LevelHierarchyEntry *LevelArray[],TopGridData *MetaData,
		      int level, ExternalBoundary *Exterior);
 
int ComputeRandomForcingNormalization(LevelHierarchyEntry *LevelArray[],
                                      int level, TopGridData *MetaData,
                                      float * norm, float * pTopGridTimeStep);
int CreateSiblingList(HierarchyEntry ** Grids, int NumberOfGrids, SiblingGridList *SiblingList, 
		      int StaticLevelZero,TopGridData * MetaData,int level);

#ifdef FLUX_FIX
int CreateSUBlingList(TopGridData *MetaData,
		      HierarchyEntry *Grids[],
		      int NumberOfGrids,
		      LevelHierarchyEntry ***SUBlingList);
int DeleteSUBlingList(int NumberOfGrids,
		      LevelHierarchyEntry **SUBlingList);
#endif

int StarParticleInitialize(LevelHierarchyEntry *LevelArray[], int ThisLevel,
			   TopGridData *MetaData, Star *&AllStars);
int StarParticleFinalize(HierarchyEntry *Grids[], TopGridData *MetaData,
			 int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
			 int level, Star *&AllStars);
int AdjustRefineRegion(LevelHierarchyEntry *LevelArray[], 
		       TopGridData *MetaData, int EL_level);

#ifdef TRANSFER
int EvolvePhotons(TopGridData *MetaData,LevelHierarchyEntry *LevelArray[],
		  Star *AllStars, int level);
int RadiativeTransferPrepare(LevelHierarchyEntry *LevelArray[], int level,
			     TopGridData *MetaData, Star *&AllStars,
			     float dtLevelAbove);
#endif

int SetLevelTimeStep(HierarchyEntry *Grids[],
        int NumberOfGrids, int level,
        float *dtThisLevelSoFar, float *dtThisLevel,
        float dtLevelAbove);

void my_exit(int status);
 
 
int LevelCycleCount[MAX_DEPTH_OF_HIERARCHY];
double LevelWallTime[MAX_DEPTH_OF_HIERARCHY];
double LevelZoneCycleCount[MAX_DEPTH_OF_HIERARCHY];
double LevelZoneCycleCountPerProc[MAX_DEPTH_OF_HIERARCHY];
 
static float norm = 0.0;            //AK
static float TopGridTimeStep = 0.0; //AK
#ifdef STATIC_SIBLING_LIST
static int StaticLevelZero = 1;
#else
static int StaticLevelZero = 0;
#endif

int EvolveLevel(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[],
		int level, float dtLevelAbove, ExternalBoundary *Exterior)
{
  /* Declarations */

  int dbx = 0;
 
  FLOAT When;
  //float dtThisLevelSoFar = 0.0, dtThisLevel, dtGrid, dtActual, dtLimit;
  float dtThisLevelSoFar = 0.0, dtThisLevel;
  int cycle = 0, counter = 0, grid1, subgrid, grid2;
  HierarchyEntry *NextGrid;
  int dummy_int;
 
  // Update lcaperf "level" attribute

  Eint32 jb_level = level;
#ifdef USE_JBPERF
  jbPerf.attribute ("level",&jb_level,JB_INT);
#endif

  /* Create an array (Grids) of all the grids. */

  typedef HierarchyEntry* HierarchyEntryPointer;
  HierarchyEntry **Grids;
  int NumberOfGrids = GenerateGridArray(LevelArray, level, &Grids);
  int *NumberOfSubgrids = new int[NumberOfGrids];
  fluxes ***SubgridFluxesEstimate = new fluxes **[NumberOfGrids];

  TIME_MSG("Entered EvolveLevel");

#ifdef FLUX_FIX
  /* Create a SUBling list of the subgrids */
 
  LevelHierarchyEntry **SUBlingList;
#endif


  /* Initialize the chaining mesh used in the FastSiblingLocator. */

  if (dbx) fprintf(stderr, "EL: Initialize FSL \n"); 
  SiblingGridList *SiblingList = new SiblingGridList[NumberOfGrids];
  CreateSiblingList(Grids, NumberOfGrids, SiblingList, StaticLevelZero,MetaData,level);
  
  /* On the top grid, adjust the refine region so that only the finest
     particles are included.  We don't want the more massive particles
     to contaminate the high-resolution region. */

  AdjustRefineRegion(LevelArray, MetaData, level);

  /* ================================================================== */
  /* For each grid: a) interpolate boundaries from its parent.
                    b) copy any overlapping zones.  */
 
#ifdef FAST_SIB
  if (SetBoundaryConditions(Grids, NumberOfGrids, SiblingList,
			    level, MetaData, Exterior, LevelArray[level]) == FAIL)
    ENZO_FAIL("");
#else
  if (SetBoundaryConditions(Grids, NumberOfGrids, level, MetaData,
                            Exterior, LevelArray[level]) == FAIL)
    ENZO_FAIL("");
#endif
 

  /* Clear the boundary fluxes for all Grids (this will be accumulated over
     the subcycles below (i.e. during one current grid step) and used to by the
     current grid to correct the zones surrounding this subgrid (step #18). */
 
  TIME_MSG("after SetBoundaryConditions");

  for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
    Grids[grid1]->GridData->ClearBoundaryFluxes();
 
  /* ================================================================== */
  /* Loop over grid timesteps until the elapsed time equals the timestep
     from the level above (or loop once for the top level). */
 
  while (dtThisLevelSoFar < dtLevelAbove) {
 
    SetLevelTimeStep(Grids, NumberOfGrids, level, 
        &dtThisLevelSoFar, &dtThisLevel, dtLevelAbove);

    /* Initialize the star particles */

    Star *AllStars = NULL;
    StarParticleInitialize(LevelArray, level, MetaData, AllStars);

    /* Initialize the radiative transfer */

#ifdef TRANSFER
    RadiativeTransferPrepare(LevelArray, level, MetaData, AllStars, 
				   dtLevelAbove);
#endif /* TRANSFER */

    /* For each grid, compute the number of it's subgrids. */
 
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
      NextGrid = Grids[grid1]->NextGridNextLevel;
      counter = 0;
      while (NextGrid != NULL) {
	NextGrid = NextGrid->NextGridThisLevel;
	if (++counter > MAX_NUMBER_OF_SUBGRIDS) {
	  fprintf(stderr, "More subgrids than MAX_NUMBER_OF_SUBGRIDS.\n");
	  ENZO_FAIL("");
	}
      }
      NumberOfSubgrids[grid1] = counter + 1;
    }
 
    TIME_MSG("Before subgrid fluxes");

    CreateFluxes(Grids,SubgridFluxesEstimate,NumberOfGrids,NumberOfSubgrids);

    TIME_MSG("After subgrid fluxes");


    /* ------------------------------------------------------- */
    /* Prepare the density field (including particle density). */
 
//  fprintf(stderr, "%"ISYM": EvolveLevel: Enter PrepareDensityField\n", MyProcessorNumber);
 
    When = 0.5;
 
#ifdef FAST_SIB
     PrepareDensityField(LevelArray, SiblingList, level, MetaData, When);
#else   // !FAST_SIB
     PrepareDensityField(LevelArray, level, MetaData, When);
#endif  // end FAST_SIB
 
 
//  fprintf(stderr, "%"ISYM": EvolveLevel: Exit PrepareDensityField\n", MyProcessorNumber);
 
    /* Prepare normalization for random forcing. Involves top grid only. */
 
    ComputeRandomForcingNormalization(LevelArray, 0, MetaData,
                                             &norm, &TopGridTimeStep);
 
    /* ------------------------------------------------------- */
    /* Evolve all grids by timestep dtThisLevel. */
 
    TIME_MSG("EvolveLevel: before main loop");
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
 
      // dcc problem analysis cut  start

      /* Call analysis routines. */

//      if (ProblemType == 24)
//	Grids[grid1]->GridData->SphericalInfallGetProfile(level, 1);
//      if (ProblemType == 30)
//	Grids[grid1]->GridData->AnalyzeTrackPeaks(level, 0);
//      if (ProblemType == 27)
//	if (Grids[grid1]->GridData->ReturnProcessorNumber()==MyProcessorNumber){
//	  float AM[3], MeanVelocity[3], DMVelocity[3];
//	  FLOAT Center[] = {0,0,0}, CenterOfMass[3], DMCofM[3];
//	  Grids[grid1]->GridData->CalculateAngularMomentum(Center, AM,
//			   MeanVelocity, DMVelocity, CenterOfMass, DMCofM);
//	  fprintf(stdout, "level = %"ISYM" %"ISYM" %"ISYM"  Vel %"FSYM" %"FSYM" %"FSYM"  DMVel %"FSYM" %"FSYM" %"FSYM"  CofM %"PSYM" %"PSYM" %"PSYM"  DMCofM %"FSYM" %"FSYM" %"FSYM"\n",
//		level, LevelCycleCount[level], grid1, MeanVelocity[0],
//		MeanVelocity[1], MeanVelocity[2],
//		DMVelocity[0], DMVelocity[1], DMVelocity[2],
//		-CenterOfMass[0], -CenterOfMass[1], -CenterOfMass[2],
//		DMCofM[0], DMCofM[1], DMCofM[2]);
//	}
 
      /* Gravity: compute acceleration field for grid and particles. */
 
      if (SelfGravity) {
	int Dummy;
	if (level <= MaximumGravityRefinementLevel) {
 
	  /* Compute the potential. */
 
	  if (level > 0)
	    Grids[grid1]->GridData->SolveForPotential(Dummy, level);
	  Grids[grid1]->GridData->ComputeAccelerations(level);
	}
	  /* otherwise, interpolate potential from coarser grid, which is
	     now done in PrepareDensity. */
 
      } // end: if (SelfGravity)
 
      /* Gravity: compute field due to preset sources. */
 
      Grids[grid1]->GridData->ComputeAccelerationFieldExternal();
 
      /* Radiation Pressure: add to acceleration field */

#ifdef TRANSFER
      Grids[grid1]->GridData->AddRadiationPressureAcceleration();
#endif /* TRANSFER */

      // dcc gravity cut stop


      /* Check for energy conservation. */
/*
      if (ComputePotential)
	if (CheckEnergyConservation(Grids, grid, NumberOfGrids, level,
				    dtThisLevel) == FAIL) {
	  fprintf(stderr, "Error in CheckEnergyConservation.\n");
	  ENZO_FAIL("");
	}
*/
#ifdef SAB
    } // End of loop over grids

  //dcc cut SAB start 
    //This ensures that all subgrids agree in the boundary.
    //Not a big deal for hydro, but essential for DivB = 0 in MHD runs.
    //Only called on level > 0 because the root grid is dealt with differently than SG's.


    if ( (SelfGravity || UniformGravity || PointSourceGravity) && level > 0) {
#ifdef FAST_SIB
      if( SetAccelerationBoundary(Grids, NumberOfGrids,
				  SiblingList,
				  level, MetaData,
				  Exterior, LevelArray[level], LevelCycleCount[level]) == FAIL ) {
	fprintf(stderr,"Error with AccelerationBoundary.\n");
	ENZO_FAIL("");
      }
#else
      if( SetAccelerationBoundary(Grids, NumberOfGrids,
				  level, MetaData,
				  Exterior, LevelArray[level], LevelCycleCount[level]) == FAIL ) {
	fprintf(stderr,"Error with AccelerationBoundary.\n");
	ENZO_FAIL("");
      }
#endif
    }

    //dcc cut SAB stop
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
#endif //SAB.
      /* Copy current fields (with their boundaries) to the old fields
	  in preparation for the new step. */
 
      Grids[grid1]->GridData->CopyBaryonFieldToOldBaryonField();

      //dcc cut Forcing to problem specific analysis

      /* Add RandomForcing fields to velocities after the copying of current
         fields to old. I also update the total energy accordingly here.
         It makes no sense to force on the very first time step. */
 
      if (MetaData->CycleNumber > 0)
        Grids[grid1]->GridData->AddRandomForcing(&norm, TopGridTimeStep);

      //dcc cut stop Forcing

      /* Call hydro solver and save fluxes around subgrids. */

      Grids[grid1]->GridData->SolveHydroEquations(LevelCycleCount[level],
	    NumberOfSubgrids[grid1], SubgridFluxesEstimate[grid1], level);
      /* Solve the radiative transfer */
	
#ifdef TRANSFER
      EvolvePhotons(MetaData, LevelArray, AllStars, level);
#endif /* TRANSFER */

      /* Solve the cooling and species rate equations. */
 
//      fprintf(stderr, "%"ISYM": Calling SolveCoolAndRateEquations\n", MyProcessorNumber);

      Grids[grid1]->GridData->MultiSpeciesHandler();


    //dcc cut start particles
      /* Update particle positions (if present). */
 
      UpdateParticlePositions(Grids[grid1]->GridData);

      /* Include 'star' particle creation and feedback.
         (first, set the under_subgrid field). */
 
      if (StarParticleCreation || StarParticleFeedback) {
	Grids[grid1]->GridData->ZeroSolutionUnderSubgrid(NULL,
						 ZERO_UNDER_SUBGRID_FIELD);
	LevelHierarchyEntry *Temp2 = LevelArray[level+1];
	while (Temp2 != NULL) {
	  Grids[grid1]->GridData->ZeroSolutionUnderSubgrid(Temp2->GridData,
					 ZERO_UNDER_SUBGRID_FIELD);
	  Temp2 = Temp2->NextGridThisLevel;
	}
      }

      if (StarParticleCreation || StarParticleFeedback) {
	if (Grids[grid1]->GridData->StarParticleHandler(level) == FAIL) {
	  fprintf(stderr, "Error in grid->StarParticleWrapper");
	  ENZO_FAIL("");
	}
      }
 
      /* Gravity: clean up AccelerationField. */

      if (SelfGravity || UniformGravity || PointSourceGravity) {
	if (level != MaximumGravityRefinementLevel ||
	    MaximumGravityRefinementLevel == MaximumRefinementLevel)
	  Grids[grid1]->GridData->DeleteAccelerationField();
	Grids[grid1]->GridData->DeleteParticleAcceleration();
      }
 
      /* Update current problem time of this subgrid. */
 
      Grids[grid1]->GridData->SetTimeNextTimestep();
 
      /* If using comoving co-ordinates, do the expansion terms now. */
 
      if (ComovingCoordinates)
	Grids[grid1]->GridData->ComovingExpansionTerms();
 
    }  // end loop over grids
 
    /* For each grid: a) interpolate boundaries from the parent grid.
                      b) copy any overlapping zones from siblings. */
 
    TIME_MSG("EvolveLevel: after main loop");

#ifdef FAST_SIB
    SetBoundaryConditions(Grids, NumberOfGrids, SiblingList,
			      level, MetaData, Exterior, LevelArray[level]);
#else
    SetBoundaryConditions(Grids, NumberOfGrids, level, MetaData,
                              Exterior, LevelArray[level]);
#endif

    TIME_MSG("after SetBoundaryConditions");

    /* Finalize (accretion, feedback, etc.) star particles */
 
    StarParticleFinalize(Grids, MetaData, NumberOfGrids, LevelArray,
			     level, AllStars);
    /* If cosmology, then compute grav. potential for output if needed. */

    //dcc cut second potential cut: Duplicate?
 
    if (ComovingCoordinates && SelfGravity && WritePotential) {
      CopyGravPotential = TRUE;
      When = 0.0;
 
#ifdef FAST_SIB
      if (PrepareDensityField(LevelArray, SiblingList, level, MetaData, When) == FAIL) {
        fprintf(stderr, "Error in PrepareDensityField.\n");
        ENZO_FAIL("");
      }
#else   // !FAST_SIB
      if (PrepareDensityField(LevelArray, level, MetaData, When) == FAIL) {
        fprintf(stderr, "Error in PrepareDensityField.\n");
        ENZO_FAIL("");
      }
#endif  // end FAST_SIB
 
      CopyGravPotential = FALSE;
 
      for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
        int Dummy;
        if (level <= MaximumGravityRefinementLevel) {
 
          /* Compute the potential. */
 
          if (level > 0)
            if (Grids[grid1]->GridData->SolveForPotential(Dummy, level)
                == FAIL) {
              fprintf(stderr, "Error in grid->SolveForPotential.\n");
              ENZO_FAIL("");
            }
          // fprintf(stderr, "Call CP from EvolveLevel\n");
          Grids[grid1]->GridData->CopyPotentialToBaryonField();
        }
        /* otherwise output empty potential field. */
 
      } //  end loop over grids
    } // if WritePotential
 
    /* For each grid, delete the GravitatingMassFieldParticles. */
 
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
      Grids[grid1]->GridData->DeleteGravitatingMassFieldParticles();
 
    /* ----------------------------------------- */
    /* Evolve the next level down (recursively). */
 
    MetaData->FirstTimestepAfterRestart = FALSE;
    if (dbx) fprintf(stderr, "EL Level %"ISYM" going to Level %"ISYM"\n", level, level+1);
    if (LevelArray[level+1] != NULL) {
      if (EvolveLevel(MetaData, LevelArray, level+1, dtThisLevel, Exterior) == FAIL) {
	fprintf(stderr, "Error in EvolveLevel (%"ISYM").\n", level);
	ENZO_FAIL("");
      }
    }

    // Update lcaperf "level" attribute

#ifdef USE_JBPERF
    jbPerf.attribute ("level",&jb_level,JB_INT);
#endif

    OutputFromEvolveLevel(LevelArray,MetaData,level,Exterior);
    /* Update SubcycleNumber if this is the bottom of the hierarchy --
       Note that this not unique based on which level is the highest,
       it just keeps going */

    if (LevelArray[level+1] == NULL) MetaData->SubcycleNumber += 1;  

    if (dbx) fprintf(stderr, "EL Level %"ISYM" returns from Level %"ISYM"\n", level, level+1);

    /* ------------------------------------------------------- */
    /* For each grid,
     * (a) project the subgrid's solution into this grid (step #18)
     * (b) correct for the difference between this grid's fluxes and the
     *     subgrid's fluxes. (step #19)
     */
 
    TIME_MSG("Before update from finer grids");
    //dcc cut start flux fix
#ifdef FLUX_FIX

    SUBlingList = new LevelHierarchyEntry*[NumberOfGrids];
    for(int list=0; list < NumberOfGrids; list++)
      SUBlingList[list] = NULL;

 
    if (FluxCorrection) {

      /* Fill in the SUBling list */

      if (dbx) fprintf(stderr, "EL: CSL entry \n");
      if (CreateSUBlingList(MetaData, Grids,
                              NumberOfGrids, &SUBlingList) == FAIL) {
        fprintf(stderr, "Error in CreateSUBlingList.\n");
        ENZO_FAIL("");
      }
      if (dbx) fprintf(stderr, "EL: CSL exit \n");
    }
    //dcc cut stop flux fix
/* 
    LevelHierarchyEntry *NextMonkey;
 
    for(grid1 = 0; grid1 < NumberOfGrids; grid1++){
      NextMonkey = SUBlingList[grid1];
      while (NextMonkey != NULL) {
        // fprintf(stderr, "SGcheckEL%"ISYM": SUBling[%"ISYM"]->Grid pointer %p\n",
        //         MyProcessorNumber, grid1, NextMonkey->GridData);
        NextMonkey=NextMonkey->NextGridThisLevel;
      }
    }
*/

#endif

#ifdef FLUX_FIX
    if (UpdateFromFinerGrids(level, Grids, NumberOfGrids, NumberOfSubgrids,
			     SubgridFluxesEstimate,
			     SUBlingList,
			     MetaData) == FAIL)
      ENZO_FAIL("");
#else
    if (UpdateFromFinerGrids(level, Grids, NumberOfGrids, NumberOfSubgrids,
			     SubgridFluxesEstimate) == FAIL)
      ENZO_FAIL("");
#endif

    if (dbx) fprintf(stderr, "OK after UpdateFromFinerGrids \n");

#ifdef FLUX_FIX
    if ( FluxCorrection ) {
      /* Clean up SUBlings */
      if (DeleteSUBlingList( NumberOfGrids, SUBlingList ) == FAIL) {
        fprintf(stderr, "Error in DeleteSUBlingList.\n");
        ENZO_FAIL("");
      }
    }
#endif

    if (dbx) fprintf(stderr, "OK after DeleteSUBlingList \n");

  /* ------------------------------------------------------- */
  /* Add the saved fluxes (in the last subsubgrid entry) to the exterior
     fluxes for this subgrid .
     (Note: this must be done after CorrectForRefinedFluxes). */

    FinalizeFluxes(Grids,SubgridFluxesEstimate,NumberOfGrids,NumberOfSubgrids);

    /* Recompute radiation field, if requested. */
 
    TIME_MSG("Done saving fluxes");

    if (RadiationFieldType >= 10 && RadiationFieldType <= 11 &&
	level <= RadiationFieldLevelRecompute)
      if (RadiationFieldUpdate(LevelArray, level, MetaData) == FAIL) {
	fprintf(stderr, "Error in RecomputeRadiationField.\n");
	ENZO_FAIL("");
      }
 
    /* Rebuild the Grids on the next level down.
       Don't bother on the last cycle, as we'll rebuild this grid soon. */
 
    if (dtThisLevelSoFar < dtLevelAbove) {
      if (RebuildHierarchy(MetaData, LevelArray, level) == FAIL) {
	fprintf(stderr, "Error in RebuildHierarchy.\n");
	ENZO_FAIL("");
      }
    }

    /* Count up number of grids on this level. */
    // ? dcc
    int GridMemory, NumberOfCells, CellsTotal, Particles;
    float AxialRatio, GridVolume;
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
      Grids[grid1]->GridData->CollectGridInformation
        (GridMemory, GridVolume, NumberOfCells, AxialRatio, CellsTotal, Particles);
      LevelZoneCycleCount[level] += NumberOfCells;
      if (MyProcessorNumber == Grids[grid1]->GridData->ReturnProcessorNumber())
	LevelZoneCycleCountPerProc[level] += NumberOfCells;
    }
 
    cycle++;
    LevelCycleCount[level]++;
 
  } // end of loop over subcycles
 
  if (debug)
    fprintf(stdout, "EvolveLevel[%"ISYM"]: NumberOfSubCycles = %"ISYM" (%"ISYM" total)\n", level,
           cycle, LevelCycleCount[level]);
 
  /* If possible & desired, report on memory usage. */
 
  //  if (debug)
  ReportMemoryUsage("Memory usage report: Evolve Level");
 
#ifdef USE_JBPERF
  jbPerf.attribute ("level",0,JB_NULL);
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
 
  /* Clean up the sibling list. */

  if (( StaticLevelZero == 1 && level != 0 ) || StaticLevelZero == 0 ) {
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
      delete [] SiblingList[grid1].GridList;
    delete [] SiblingList;
  }

  if (dbx) fprintf(stderr, "Return from EL Level %"ISYM"\n", level);
 
  return SUCCESS;
 
}
