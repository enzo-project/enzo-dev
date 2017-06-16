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
/  modified9:  June, 2009 by MJT, DC, JHW, TA
/                Cleaned up error handling and created new routines for
/                computing the timestep, output, handling fluxes
/  modified10: July, 2009 by Sam Skillman
/                Added shock analysis
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
#include "preincludes.h"
 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include "EnzoTiming.h"
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
#ifdef TRANSFER
#include "ImplicitProblemABC.h"
#endif
#ifdef NEW_PROBLEM_TYPES
#include "EventHooks.h"
#else
void RunEventHooks(char *, HierarchyEntry *Grid[], TopGridData &MetaData) {}
#endif
 
/* function prototypes */
 
#ifdef TRANSFER
#define IMPLICIT_MACRO , ImplicitSolver
#else
#define IMPLICIT_MACRO 
#endif

#define EXTRA_OUTPUT_MACRO(A,B) ExtraOutput(A,LevelArray,MetaData,level,Exterior IMPLICIT_MACRO,B);
int ExtraOutput(int output_flag, LevelHierarchyEntry *LevelArray[],TopGridData *MetaData, int level, ExternalBoundary *Exterior
#ifdef TRANSFER
			  , ImplicitProblemABC *ImplicitSolver
#endif
        , char * output_string);

int ComputeDednerWaveSpeeds(TopGridData *MetaData,LevelHierarchyEntry *LevelArray[], 
			    int level, FLOAT dt0);
int  RebuildHierarchy(TopGridData *MetaData,
		      LevelHierarchyEntry *LevelArray[], int level);
int  ReportMemoryUsage(char *header = NULL);
int  UpdateParticlePositions(grid *Grid);
int  CheckEnergyConservation(HierarchyEntry *Grids[], int grid,
			     int NumberOfGrids, int level, float dt);
int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);
int WriteStreamData(LevelHierarchyEntry *LevelArray[], int level,
		    TopGridData *MetaData, int *CycleCount, int open=FALSE);
int CallProblemSpecificRoutines(TopGridData * MetaData, HierarchyEntry *ThisGrid,
				int GridNum, float *norm, float TopGridTimeStep, 
				int level, int LevelCycleCount[]);  

#ifdef FAST_SIB
int PrepareDensityField(LevelHierarchyEntry *LevelArray[],
			int level, TopGridData *MetaData, FLOAT When, SiblingGridList **SiblingGridListStorage);
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

int UpdateFromFinerGrids(int level, HierarchyEntry *Grids[], int NumberOfGrids,
			 int NumberOfSubgrids[],
			 fluxes **SubgridFluxesEstimate[],
			 LevelHierarchyEntry *SUBlingList[],
			 TopGridData *MetaData);

int CreateFluxes(HierarchyEntry *Grids[],fluxes **SubgridFluxesEstimate[],
		 int NumberOfGrids,int NumberOfSubgrids[]);		 
int FinalizeFluxes(HierarchyEntry *Grids[],fluxes **SubgridFluxesEstimate[],
		 int NumberOfGrids,int NumberOfSubgrids[]);		 
int RadiationFieldUpdate(LevelHierarchyEntry *LevelArray[], int level,
			 TopGridData *MetaData);


int OutputFromEvolveLevel(LevelHierarchyEntry *LevelArray[],TopGridData *MetaData,
			  int level, ExternalBoundary *Exterior, int OutputNow
#ifdef TRANSFER
			  , ImplicitProblemABC *ImplicitSolver
#endif
			  );
 
int ComputeRandomForcingNormalization(LevelHierarchyEntry *LevelArray[],
                                      int level, TopGridData *MetaData,
                                      float * norm, float * pTopGridTimeStep);

int ComputeStochasticForcing(TopGridData *MetaData,
        HierarchyEntry *Grids[], int NumberOfGrids);

int ClusterSMBHSumGasMass(HierarchyEntry *Grids[], int NumberOfGrids, int level);
int CreateSiblingList(HierarchyEntry ** Grids, int NumberOfGrids, SiblingGridList *SiblingList, 
		      int StaticLevelZero,TopGridData * MetaData,int level);

#ifdef FAST_SIB 
int CreateSUBlingList(TopGridData *MetaData,
		      LevelHierarchyEntry *LevelArray[], int level,
		      SiblingGridList SiblingList[],
		      LevelHierarchyEntry ***SUBlingList);
#else
int CreateSUBlingList(TopGridData *MetaData,
		      LevelHierarchyEntry *LevelArray[], int level,
		      LevelHierarchyEntry ***SUBlingList);
#endif /* FAST_SIB */
int DeleteSUBlingList(int NumberOfGrids,
		      LevelHierarchyEntry **SUBlingList);

int StarParticleInitialize(HierarchyEntry *Grids[], TopGridData *MetaData,
			   int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
			   int ThisLevel, Star *&AllStars,
			   int TotalStarParticleCountPrevious[]);
int StarParticleFinalize(HierarchyEntry *Grids[], TopGridData *MetaData,
			 int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
			 int level, Star *&AllStars,
			 int TotalStarParticleCountPrevious[], int &OutputNow);
int AdjustRefineRegion(LevelHierarchyEntry *LevelArray[], 
		       TopGridData *MetaData, int EL_level);
int AdjustMustRefineParticlesRefineToLevel(TopGridData *MetaData, int EL_level);

#ifdef TRANSFER
int EvolvePhotons(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[],
		  Star *&AllStars, FLOAT GridTime, int level, int LoopTime = TRUE);
int RadiativeTransferPrepare(LevelHierarchyEntry *LevelArray[], int level,
			     TopGridData *MetaData, Star *&AllStars,
			     float dtLevelAbove);
int RadiativeTransferCallFLD(LevelHierarchyEntry *LevelArray[], int level,
			     TopGridData *MetaData, Star *AllStars, 
			     ImplicitProblemABC *ImplicitSolver);
#endif

int SetLevelTimeStep(HierarchyEntry *Grids[],
        int NumberOfGrids, int level,
        float *dtThisLevelSoFar, float *dtThisLevel,
        float dtLevelAbove);

void my_exit(int status);
 
int CallPython(LevelHierarchyEntry *LevelArray[], TopGridData *MetaData,
               int level, int from_topgrid);
int MovieCycleCount[MAX_DEPTH_OF_HIERARCHY];
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

extern int RK2SecondStepBaryonDeposit;


int EvolveLevel(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[],
		int level, float dtLevelAbove, ExternalBoundary *Exterior
#ifdef TRANSFER
		, ImplicitProblemABC *ImplicitSolver
#endif
    , FLOAT dt0, SiblingGridList *SiblingGridListStorage[] 
		)
{
  /* Declarations */

  int dbx = 0;

  FLOAT When, GridTime;
  //float dtThisLevelSoFar = 0.0, dtThisLevel, dtGrid, dtActual, dtLimit;
  //float dtThisLevelSoFar = 0.0, dtThisLevel;
  int cycle = 0, counter = 0, grid1, subgrid, grid2;
  HierarchyEntry *NextGrid;
  int dummy_int, OutputNow = FALSE;

  char level_name[MAX_LINE_LENGTH];
  sprintf(level_name, "Level_%02"ISYM, level);
    
  // Update lcaperf "level" attribute
  Eint32 lcaperf_level = level;
#ifdef USE_LCAPERF
  lcaperf.attribute ("level",&lcaperf_level,LCAPERF_INT);
#endif
  
  /* Create an array (Grids) of all the grids. */

  typedef HierarchyEntry* HierarchyEntryPointer;
  HierarchyEntry **Grids;
  int NumberOfGrids = GenerateGridArray(LevelArray, level, &Grids);
  int *NumberOfSubgrids = new int[NumberOfGrids];
  fluxes ***SubgridFluxesEstimate = new fluxes **[NumberOfGrids];
  int *TotalStarParticleCountPrevious = new int[NumberOfGrids];
  RunEventHooks("EvolveLevelTop", Grids, *MetaData);

  /* Create a SUBling list of the subgrids */
  LevelHierarchyEntry **SUBlingList;

  /* Initialize the chaining mesh used in the FastSiblingLocator. */

  if (dbx) fprintf(stderr, "EL: Initialize FSL \n"); 
  SiblingGridList *SiblingList = new SiblingGridList[NumberOfGrids];
  SiblingGridListStorage[level] = SiblingList;
  CreateSiblingList(Grids, NumberOfGrids, SiblingList, StaticLevelZero,MetaData,level);
  
  /* Adjust the refine region so that only the finest particles 
     are included.  We don't want the more massive particles
     to contaminate the high-resolution region. */

  AdjustRefineRegion(LevelArray, MetaData, level);

  //EMISSIVITY if cleared here will not reach the FLD solver in 2.0, finding better place
  /* Adjust MustRefineParticlesRefineToLevel parameter if requested */
  AdjustMustRefineParticlesRefineToLevel(MetaData, level);

  /* ================================================================== */
  /* For each grid: a) interpolate boundaries from its parent.
                    b) copy any overlapping zones.  */
 
  if (CheckpointRestart == FALSE) {
#ifdef FAST_SIB
    if (SetBoundaryConditions(Grids, NumberOfGrids, SiblingList,
                  level, MetaData, Exterior, LevelArray[level]) == FAIL)
      ENZO_FAIL("Error in SetBoundaryConditions (FastSib)");
#else
    if (SetBoundaryConditions(Grids, NumberOfGrids, level, MetaData,
                              Exterior, LevelArray[level]) == FAIL)
      ENZO_FAIL("Error in SetBoundaryConditions (SlowSib)");
#endif
  }
 
  Grids[0]->GridData->SetNumberOfColours();
  /* Clear the boundary fluxes for all Grids (this will be accumulated over
     the subcycles below (i.e. during one current grid step) and used to by the
     current grid to correct the zones surrounding this subgrid (step #18). 

     If we're just coming in off a CheckpointRestart, instead we take the
     fluxes that were stored in the file and then in the Grid object, and we 
     put them into the SubgridFluxesEstimate array. */
 
  if(CheckpointRestart == TRUE) {
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
      if (Grids[grid1]->GridData->FillFluxesFromStorage(
        &NumberOfSubgrids[grid1],
        &SubgridFluxesEstimate[grid1]) != -1) {
        /*fprintf(stderr, "Level: %"ISYM" Grid: %"ISYM" NS: %"ISYM"\n",
            level, grid1, NumberOfSubgrids[grid1]);*/
      }
    }
  } else {
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++){
      Grids[grid1]->GridData->ClearBoundaryFluxes();
      if ( level > 0 )
      Grids[grid1]->GridData->ClearAvgElectricField();
    }

  }
 
  /* After we calculate the ghost zones, we can initialize streaming
     data files (only on level 0) */

  if (MetaData->FirstTimestepAfterRestart == TRUE && level == 0)
    WriteStreamData(LevelArray, level, MetaData, MovieCycleCount);


  /* ================================================================== */
  /* Loop over grid timesteps until the elapsed time equals the timestep
     from the level above (or loop once for the top level). */
 
  EXTRA_OUTPUT_MACRO(1, "Before Time Loop")

  while ((CheckpointRestart == TRUE)
        || (dtThisLevelSoFar[level] < dtLevelAbove)) {
    if(CheckpointRestart == FALSE) {

    TIMER_START(level_name);
    SetLevelTimeStep(Grids, NumberOfGrids, level, 
        &dtThisLevelSoFar[level], &dtThisLevel[level], dtLevelAbove);

    TimeSinceRebuildHierarchy[level] += dtThisLevel[level];

    /* If StarFormationOncePerRootGridTimeStep, stars are only created
    once per root grid time step and only on MaximumRefinementLevel
    grids. The following sets the MakeStars flag for all
    MaximumRefinementLevel grids when level==0. Post star formation,
    MakeStars is unset in Grid::StarParticleHandler() in order to
    prevent further star formation until the next root grid time
    step. */

    /* Currently (April 2012) this is only implemented for H2REG_STAR,
    and MakeStars is completely ignored in all other star makers. */

    if ( (STARMAKE_METHOD(H2REG_STAR)) && 
	 (level==0) && 
	 (StarFormationOncePerRootGridTimeStep) ) {
      /* At top level, set Grid::MakeStars to 1 for all highest
	 refinement level grids. */
      LevelHierarchyEntry *Temp;
      Temp = LevelArray[MaximumRefinementLevel];
      int count=0;
      while (Temp != NULL) {
	Temp->GridData->SetMakeStars();
	Temp = Temp->NextGridThisLevel;
	count++;
      }
      // if(MyProcessorNumber == ROOT_PROCESSOR) 
      // 	fprintf(stderr,"Set MakeStars=1 for %d MaximumRefinementLevel grids.\n",count);

      TopGridTimeStep = LevelArray[0]->GridData->ReturnTimeStep();

    }

    /* Streaming movie output (write after all parent grids are
       updated) */

    WriteStreamData(LevelArray, level, MetaData, MovieCycleCount);

    /* Initialize the star particles */

    Star *AllStars = NULL;
    StarParticleInitialize(Grids, MetaData, NumberOfGrids, LevelArray,
			   level, AllStars, TotalStarParticleCountPrevious);

    /* Calculate ClusterSMBHColdGasMass */

    ClusterSMBHSumGasMass(Grids, NumberOfGrids, level);

#ifdef TRANSFER
    /* Initialize the radiative transfer */

    TIMER_STOP(level_name);
    RadiativeTransferPrepare(LevelArray, level, MetaData, AllStars, 
			     dtLevelAbove);
    RadiativeTransferCallFLD(LevelArray, level, MetaData, AllStars, 
			     ImplicitSolver);

    /* Solve the radiative transfer */
	
    GridTime = Grids[0]->GridData->ReturnTime() + dtThisLevel[level];
    EvolvePhotons(MetaData, LevelArray, AllStars, GridTime, level);
    TIMER_START(level_name);
 
#endif /* TRANSFER */

    /* trying to clear Emissivity here after FLD uses it, doesn't work */
 
    CreateFluxes(Grids,SubgridFluxesEstimate,NumberOfGrids,NumberOfSubgrids);

    if ((HydroMethod == MHD_RK) && (level == 0))
      ComputeDednerWaveSpeeds(MetaData, LevelArray, level, dt0);
	
    if (debug1 && HydroMethod == MHD_RK && (MyProcessorNumber == ROOT_PROCESSOR)) 
      fprintf(stderr, "wave speeds: timestep: %"GSYM"  C_h: %"GSYM"  C_p: %"GSYM"\n ", 
	       dt0, C_h, C_p);
    /* ------------------------------------------------------- */
    /* Prepare the density field (including particle density). */

    When = 0.5;

#ifdef FAST_SIB
     PrepareDensityField(LevelArray,  level, MetaData, When, SiblingGridListStorage);
#else   // !FAST_SIB
     PrepareDensityField(LevelArray, level, MetaData, When);
#endif  // end FAST_SIB
 
 
    /* Prepare normalization for random forcing. Involves top grid only. */
 
    ComputeRandomForcingNormalization(LevelArray, 0, MetaData,
				      &norm, &TopGridTimeStep);


    /* Compute stochastic force field via FFT from the spectrum. */
    if (DrivenFlowProfile) {
        if (ComputeStochasticForcing(MetaData, Grids, NumberOfGrids) == FAIL) {
            fprintf(stderr, "Error in ComputeStochasticForcing.\n");
            return FAIL;
        }
    }

    /* ------------------------------------------------------- */
    /* Evolve all grids by timestep dtThisLevel. */

    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
 
      CallProblemSpecificRoutines(MetaData, Grids[grid1], grid1, &norm, 
				  TopGridTimeStep, level, LevelCycleCount);

      /* Gravity: compute acceleration field for grid and particles. */
 
      if (SelfGravity) {
	if (level <= MaximumGravityRefinementLevel) {
 
	  /* Compute the potential. */
 
	  if (level > 0)
	    Grids[grid1]->GridData->SolveForPotential(level);
	  Grids[grid1]->GridData->ComputeAccelerations(level);
	  Grids[grid1]->GridData->CopyPotentialToBaryonField();
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

      /* Check for energy conservation. */
/*
      if (ComputePotential)
	if (CheckEnergyConservation(Grids, grid, NumberOfGrids, level,
				    dtThisLevel) == FAIL) {
	  ENZO_FAIL("Error in CheckEnergyConservation.\n");
	}
*/
#ifdef SAB
    } // End of loop over grids
    
    //Ensure the consistency of the AccelerationField
    SetAccelerationBoundary(Grids, NumberOfGrids,SiblingList,level, MetaData,
			    Exterior, LevelArray[level], LevelCycleCount[level]);
    
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
#endif //SAB.
      /* Copy current fields (with their boundaries) to the old fields
	  in preparation for the new step. */
 
      Grids[grid1]->GridData->CopyBaryonFieldToOldBaryonField();

      /* Call hydro solver and save fluxes around subgrids. */

      if( HydroMethod != HD_RK && HydroMethod != MHD_RK ){
	Grids[grid1]->GridData->SolveHydroEquations(LevelCycleCount[level],
	    NumberOfSubgrids[grid1], SubgridFluxesEstimate[grid1], level);
      }else{
          if( UseHydro ) {
        if (HydroMethod == HD_RK)
          Grids[grid1]->GridData->RungeKutta2_1stStep
              (SubgridFluxesEstimate[grid1], NumberOfSubgrids[grid1], level, Exterior);
        else if (HydroMethod == MHD_RK) {
          Grids[grid1]->GridData->MHDRK2_1stStep
              (SubgridFluxesEstimate[grid1], NumberOfSubgrids[grid1], level, Exterior);
        }
        //	dcc notes;  expansion terms probably not here.
        //	SetNextTimestep should be here.
	//	if (ComovingCoordinates)
	//	  Grids[grid1]->GridData->ComovingExpansionTerms();
    //  Grids[grid1]->GridData->SetTimeNextTimestep();

          }//use hydro
      }//hydro method
    }//grids

    if( HydroMethod == HD_RK || HydroMethod == MHD_RK ){
#ifdef FAST_SIB
      SetBoundaryConditions(Grids, NumberOfGrids, SiblingList, level, MetaData, Exterior, LevelArray[level]);
#else
      SetBoundaryConditions(Grids, NumberOfGrids, level, MetaData, Exterior, LevelArray[level]);
#endif


    if (RK2SecondStepBaryonDeposit && SelfGravity && UseHydro) {  
  
      When = 0.5;
#ifdef FAST_SIB
      PrepareDensityField(LevelArray,  level, MetaData, When, SiblingGridListStorage);
#else  
      PrepareDensityField(LevelArray, level, MetaData, When);
#endif  // end FAST_SIB


    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {

      /* Gravity: compute acceleration field for grid and particles. */
      if (RK2SecondStepBaryonDeposit && SelfGravity) {
	int Dummy;
	if (level <= MaximumGravityRefinementLevel) {
	  if (level > 0) 
	    Grids[grid1]->GridData->SolveForPotential(level) ;
	  Grids[grid1]->GridData->ComputeAccelerations(level) ;
	}
      } // end: if (SelfGravity)

      Grids[grid1]->GridData->ComputeAccelerationFieldExternal() ;

    } // End of loop over grids


#ifdef SAB    
    //Ensure the consistency of the AccelerationField
    SetAccelerationBoundary(Grids, NumberOfGrids,SiblingList,level, MetaData,
			    Exterior, LevelArray[level], LevelCycleCount[level]);

#endif //SAB.    

      }
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {

      if (UseHydro) {
        if (HydroMethod == HD_RK)
          Grids[grid1]->GridData->RungeKutta2_2ndStep
              (SubgridFluxesEstimate[grid1], NumberOfSubgrids[grid1], level, Exterior);

        else if (HydroMethod == MHD_RK) {

          Grids[grid1]->GridData->MHDRK2_2ndStep
              (SubgridFluxesEstimate[grid1], NumberOfSubgrids[grid1], level, Exterior);
          if (UseAmbipolarDiffusion) 
            Grids[grid1]->GridData->AddAmbipolarDiffusion();

          if (UseResistivity) 
            Grids[grid1]->GridData->AddResistivity();

        } // ENDIF MHD_RK

        //time1 = ReturnWallTime(); dcc get this, I guess.

        /* Add viscosity */

        if (UseViscosity) 
          Grids[grid1]->GridData->AddViscosity();

        /* If using comoving co-ordinates, do the expansion terms now. */
        if (ComovingCoordinates)
          Grids[grid1]->GridData->ComovingExpansionTerms();

      } // ENDIF UseHydro
    }//grid
    }//RK hydro

      /* Solve the cooling and species rate equations. */
 
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
      Grids[grid1]->GridData->MultiSpeciesHandler();

      /* Update particle positions (if present). */
 
      UpdateParticlePositions(Grids[grid1]->GridData);

    /*Trying after solving for radiative transfer */
#ifdef EMISSIVITY
    /*                                                                                                           
        clear the Emissivity of the level below, after the level below                                            
        updated the current level (it's parent) and before the next
        timestep at the current level.                                                                            
    */
      /*    if (StarMakerEmissivityField > 0) {
    LevelHierarchyEntry *Temp;
    Temp = LevelArray[level];
    while (Temp != NULL) {
      Temp->GridData->ClearEmissivity();
      Temp = Temp->NextGridThisLevel;
      }
      }*/
#endif


      /* Include 'star' particle creation and feedback. */

      Grids[grid1]->GridData->StarParticleHandler
	(Grids[grid1]->NextGridNextLevel, level ,dtLevelAbove, TopGridTimeStep);

      /* Include shock-finding */

      Grids[grid1]->GridData->ShocksHandler();

      /* Compute and apply thermal conduction. */
      if(IsotropicConduction || AnisotropicConduction){
	if(Grids[grid1]->GridData->ConductHeat() == FAIL){
	  ENZO_FAIL("Error in grid->ConductHeat.\n");
	}
      }

      /* Compute and Apply Cosmic Ray Diffusion */
      if(CRModel && CRDiffusion){
        if(Grids[grid1]->GridData->ComputeCRDiffusion() == FAIL){
          fprintf(stderr, "Error in grid->ComputeCRDiffusion.\n");
          return FAIL;
        } // end ComputeCRDiffusion if
      }// end CRDiffusion if

      /* Gravity: clean up AccelerationField. */

#ifndef SAB
      if ((level != MaximumGravityRefinementLevel ||
	   MaximumGravityRefinementLevel == MaximumRefinementLevel) &&
	  !PressureFree)
	Grids[grid1]->GridData->DeleteAccelerationField();
#endif //!SAB

      Grids[grid1]->GridData->DeleteParticleAcceleration();
 
      /* Update current problem time of this subgrid. */
 
      Grids[grid1]->GridData->SetTimeNextTimestep();
 
      /* If using comoving co-ordinates, do the expansion terms now. */
 
      if (ComovingCoordinates)
	Grids[grid1]->GridData->ComovingExpansionTerms();
 
    }  // end loop over grids
 
    /* Finalize (accretion, feedback, etc.) star particles */

    StarParticleFinalize(Grids, MetaData, NumberOfGrids, LevelArray,
			 level, AllStars, TotalStarParticleCountPrevious, OutputNow);

    /* For each grid: a) interpolate boundaries from the parent grid.
                      b) copy any overlapping zones from siblings. */
 
    EXTRA_OUTPUT_MACRO(2,"After SolveHydroEquations grid loop")

#ifdef FAST_SIB
    SetBoundaryConditions(Grids, NumberOfGrids, SiblingList,
			  level, MetaData, Exterior, LevelArray[level]);
#else
    SetBoundaryConditions(Grids, NumberOfGrids, level, MetaData,
			  Exterior, LevelArray[level]);
#endif
    EXTRA_OUTPUT_MACRO(25,"After SBC")

    /* If cosmology, then compute grav. potential for output if needed. */


    /* For each grid, delete the GravitatingMassFieldParticles. */
 
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
      Grids[grid1]->GridData->DeleteGravitatingMassFieldParticles();

    TIMER_STOP(level_name);
    /* ----------------------------------------- */
    /* Evolve the next level down (recursively). */
 
    MetaData->FirstTimestepAfterRestart = FALSE;

    } else { // CheckpointRestart
        // dtThisLevelSoFar set during restart
        // dtThisLevel set during restart
        // Set dtFixed on each grid to dtThisLevel
        for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
          Grids[grid1]->GridData->SetTimeStep(dtThisLevel[level]);
    }

    if (LevelArray[level+1] != NULL) {
      if (EvolveLevel(MetaData, LevelArray, level+1, dtThisLevel[level], Exterior
#ifdef TRANSFER
		      , ImplicitSolver
#endif
          ,dt0,SiblingGridListStorage
		      ) == FAIL) {
	ENZO_VFAIL("Error in EvolveLevel (%"ISYM").\n", level)
      }
    }

#ifdef USE_LCAPERF
    // Update lcaperf "level" attribute

    lcaperf.attribute ("level",&lcaperf_level,LCAPERF_INT);
#endif

    OutputFromEvolveLevel(LevelArray, MetaData, level, Exterior, OutputNow
#ifdef TRANSFER
			  , ImplicitSolver
#endif
			  );
#ifdef USE_PYTHON
    LCAPERF_START("CallPython");
    CallPython(LevelArray, MetaData, level, 0);
    LCAPERF_STOP("CallPython");
#endif

    /* Update SubcycleNumber and the timestep counter for the
       streaming data if this is the bottom of the hierarchy -- Note
       that this not unique based on which level is the highest, it
       just keeps going */

    if (LevelArray[level+1] == NULL) {
      MetaData->SubcycleNumber++;
      MetaData->MovieTimestepCounter++;
    }

    /* Once MBH particles are inserted throughout the whole grid hierarchy,
       turn off MBH creation (at the bottom of the hierarchy) */

    if (STARMAKE_METHOD(MBH_PARTICLE) && (LevelArray[level+1] == NULL)) { 
      StarParticleCreation -= pow(2, MBH_PARTICLE);  
    }

    /* ------------------------------------------------------- */
    /* For each grid,
     * (a) project the subgrid's solution into this grid (step #18)
     * (b) correct for the difference between this grid's fluxes and the
     *     subgrid's fluxes. (step #19)
     */
 
    SUBlingList = new LevelHierarchyEntry*[NumberOfGrids];
#ifdef FAST_SIB
    CreateSUBlingList(MetaData, LevelArray, level, SiblingList,
		      &SUBlingList);
#else
    CreateSUBlingList(MetaData, LevelArray, level, &SUBlingList);
#endif /* FAST_SIB */


    EXTRA_OUTPUT_MACRO(3,"Before UFG")

    UpdateFromFinerGrids(level, Grids, NumberOfGrids, NumberOfSubgrids,
			     SubgridFluxesEstimate,SUBlingList,MetaData);

    DeleteSUBlingList( NumberOfGrids, SUBlingList );

    EXTRA_OUTPUT_MACRO(4,"After UFG")


    if(UseMHDCT == TRUE && MHD_ProjectE == TRUE){
      for(grid1=0;grid1<NumberOfGrids; grid1++){
        Grids[grid1]->GridData->MHD_UpdateMagneticField(level, LevelArray[level+1], FALSE);
        }
    }//MHD True

    EXTRA_OUTPUT_MACRO(5,"After UMF")

  /* ------------------------------------------------------- */
  /* Add the saved fluxes (in the last subsubgrid entry) to the exterior
     fluxes for this subgrid .
     (Note: this must be done after CorrectForRefinedFluxes). */

    if( UseMHDCT ){
#ifdef FAST_SIB
    SetBoundaryConditions(Grids, NumberOfGrids, SiblingList,
			  level, MetaData, Exterior, LevelArray[level]);
#else
    SetBoundaryConditions(Grids, NumberOfGrids, level, MetaData,
			  Exterior, LevelArray[level]);
#endif
    }

    EXTRA_OUTPUT_MACRO(51, "After SBC")

    FinalizeFluxes(Grids,SubgridFluxesEstimate,NumberOfGrids,NumberOfSubgrids);

    /* Recompute radiation field, if requested. */
    RadiationFieldUpdate(LevelArray, level, MetaData);
 
//     //dcc cut second potential cut: Duplicate?
 
//     if (SelfGravity && WritePotential) {
//       CopyGravPotential = TRUE;
//       When = 0.0;
 
// #ifdef FAST_SIB
//       PrepareDensityField(LevelArray, SiblingList, level, MetaData, When);
// #else   // !FAST_SIB
//       PrepareDensityField(LevelArray, level, MetaData, When);
// #endif  // end FAST_SIB
 
 
//       for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
//         if (level <= MaximumGravityRefinementLevel) {
 
//           /* Compute the potential. */
 
//           if (level > 0)
//             Grids[grid1]->GridData->SolveForPotential(level);
//           Grids[grid1]->GridData->CopyPotentialToBaryonField();
//         }
//       } //  end loop over grids
//        CopyGravPotential = FALSE;

//     } // if WritePotential
 
    /* Count up number of grids on this level. */

    int GridMemory, NumberOfCells, CellsTotal, Particles;
    float AxialRatio, GridVolume;
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
      Grids[grid1]->GridData->CollectGridInformation
        (GridMemory, GridVolume, NumberOfCells, AxialRatio, CellsTotal, Particles);
      LevelZoneCycleCount[level] += NumberOfCells;
      TIMER_ADD_CELLS(level, NumberOfCells);
      if (MyProcessorNumber == Grids[grid1]->GridData->ReturnProcessorNumber())
	LevelZoneCycleCountPerProc[level] += NumberOfCells;
    }
    TIMER_SET_NGRIDS(level, NumberOfGrids);

    /* Rebuild the Grids on the next level down.
       Don't bother on the last cycle, as we'll rebuild this grid soon. */
 
    if (dtThisLevelSoFar[level] < dtLevelAbove)
      RebuildHierarchy(MetaData, LevelArray, level);

    cycle++;
    LevelCycleCount[level]++;
    LevelSubCycleCount[level]++;
    if ((MetaData->StaticHierarchy == 0) && (level < MaximumRefinementLevel)) {
      LevelSubCycleCount[level+1] = 0;
    }
 
  } // end of loop over subcycles
 
    EXTRA_OUTPUT_MACRO(6, "After Subcycle Loop")
  if (debug)
    fprintf(stdout, "EvolveLevel[%"ISYM"]: NumberOfSubCycles = %"ISYM" (%"ISYM" total, %"ISYM" sub)\n", 
            level, cycle, LevelCycleCount[level], LevelSubCycleCount[level]);
 
  /* If possible & desired, report on memory usage. */
 
  ReportMemoryUsage("Memory usage report: Evolve Level");
 
#ifdef USE_LCAPERF
  lcaperf.attribute ("level",0,LCAPERF_NULL);
#endif

  
  /* Clean up. */
 
  delete [] NumberOfSubgrids;
  delete [] Grids;
  delete [] SubgridFluxesEstimate;
  delete [] TotalStarParticleCountPrevious;

  dtThisLevel[level] = dtThisLevelSoFar[level] = 0.0;
 
  /* Clean up the sibling list. */


  if ((NumberOfGrids >1) || ( StaticLevelZero == 1 && level != 0 ) || StaticLevelZero == 0 ) {

    for (grid1 = 0; grid1 < NumberOfGrids; grid1++){
      // Need this check to match that from CreateSiblingList.C and
      // Grid_FastSiblingLocatorFindSiblings.C
      if (NumberOfGrids == 1)
        delete SiblingList[grid1].GridList;
      else
        delete [] SiblingList[grid1].GridList;
    }
    delete [] SiblingList;
    SiblingGridListStorage[level] = NULL;
  }

  return SUCCESS;
 
}
