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
 
/* function prototypes */
 
void DeleteFluxes(fluxes *Fluxes);
int  RebuildHierarchy(TopGridData *MetaData,
		      LevelHierarchyEntry *LevelArray[], int level);
int  ReportMemoryUsage(char *header = NULL);
int  UpdateParticlePositions(grid *Grid);
int  CheckEnergyConservation(HierarchyEntry *Grids[], int grid,
			     int NumberOfGrids, int level, float dt);
float CommunicationMinValue(float Value);
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

#ifdef USE_HDF5_GROUPS
int Group_WriteAllData(char *basename, int filenumber, HierarchyEntry *TopGrid,
		 TopGridData &MetaData, ExternalBoundary *Exterior,
		 FLOAT WriteTime = -1);
#else
int WriteAllData(char *basename, int filenumber, HierarchyEntry *TopGrid,
		 TopGridData &MetaData, ExternalBoundary *Exterior,
		 FLOAT WriteTime = -1);
#endif
 
int ComputeRandomForcingNormalization(LevelHierarchyEntry *LevelArray[],
                                      int level, TopGridData *MetaData,
                                      float * norm, float * pTopGridTimeStep);

int FastSiblingLocatorInitializeStaticChainingMesh(ChainingMeshStructure *Mesh, int Rank,
						   int TopGridDims[]); 
int FastSiblingLocatorInitialize(ChainingMeshStructure *Mesh, int Rank,
				 int TopGridDims[]);
int FastSiblingLocatorFinalize(ChainingMeshStructure *Mesh);
 
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
		       TopGridData *MetaData);

#ifdef TRANSFER
int EvolvePhotons(TopGridData *MetaData,LevelHierarchyEntry *LevelArray[],
		  Star *AllStars);
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

static int StaticSiblingListInitialized = 0;

#ifdef STATIC_SIBLING_LIST
static SiblingGridList StaticSiblingList[MAX_NUMBER_OF_SUBGRIDS];
static int StaticLevelZero = 1;
#else
static int StaticLevelZero = 0;
#endif

#define TIME_MESSAGING 


int EvolveLevel(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[],
		int level, float dtLevelAbove, ExternalBoundary *Exterior)
{
  /* Declarations */

  int dbx = 0;
 
  FLOAT When;
  //float dtThisLevelSoFar = 0.0, dtThisLevel, dtGrid, dtActual, dtLimit;
  float dtThisLevelSoFar = 0.0, dtThisLevel;
  int RefinementFactors[MAX_DIMENSION];
  int cycle = 0, counter = 0, grid1, subgrid, grid2;
  HierarchyEntry *NextGrid;
  int dummy_int;
 
#if defined(USE_JBPERF) && defined(JB_PERF_LEVELS)
  Eint32 jb_level = level;
  jbPerf.attribute ("level",&jb_level,JB_INT);
#endif

  /* Create an array (Grids) of all the grids. */

  JBPERF_START("evolve-level-01"); // GenerateGridArray ()

  typedef HierarchyEntry* HierarchyEntryPointer;
  HierarchyEntry **Grids;
  int NumberOfGrids = GenerateGridArray(LevelArray, level, &Grids);
  int *NumberOfSubgrids = new int[NumberOfGrids];
  fluxes ***SubgridFluxesEstimate = new fluxes **[NumberOfGrids];

  JBPERF_STOP("evolve-level-01"); // GenerateGridArray ()

  JBPERF_START("evolve-level-02"); // SetBoundaryConditions()
  TIME_MSG("Entered EvolveLevel");

#ifdef FLUX_FIX
  /* Create a SUBling list of the subgrids */
 
  LevelHierarchyEntry **SUBlingList;
#endif

  /* Initialize the chaining mesh used in the FastSiblingLocator. */

  if (dbx) fprintf(stderr, "EL: Initialize FSL \n"); 

  // If this is level 0 the SiblingList does not change and can be static

#ifdef STATIC_SIBLING_LIST
  if ( StaticLevelZero == 1 && level == 0 ) {

    if (!StaticSiblingListInitialized) {

      if (debug) fprintf(stderr, "INITIALIZE Level 0 StaticSiblingList\n");

      ChainingMeshStructure StaticChainingMesh;

      FastSiblingLocatorInitializeStaticChainingMesh
	(&StaticChainingMesh, MetaData->TopGridRank, MetaData->TopGridDims);

      for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
        Grids[grid1]->GridData->FastSiblingLocatorAddGrid(&StaticChainingMesh);

      for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
        if (Grids[grid1]->GridData->FastSiblingLocatorFindSiblings(
                              &StaticChainingMesh, &StaticSiblingList[grid1],
                              MetaData->LeftFaceBoundaryCondition,
                              MetaData->RightFaceBoundaryCondition) == FAIL) {
          fprintf(stderr, "Error in grid->FastSiblingLocatorFindSiblings.\n");
          ENZO_FAIL("");
        }

      /* Clean up the chaining mesh. */

      FastSiblingLocatorFinalize(&StaticChainingMesh);

      StaticSiblingListInitialized = 1;

    }

  } // if StaticLevelZero && level == 0
#endif

  SiblingGridList *SiblingList = new SiblingGridList[NumberOfGrids];

  ChainingMeshStructure ChainingMesh;

#ifdef STATIC_SIBLING_LIST
  if (StaticLevelZero == 1 && level == 0 ) {

    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
      SiblingList[grid1].NumberOfSiblings = StaticSiblingList[grid1].NumberOfSiblings;
      SiblingList[grid1].GridList = StaticSiblingList[grid1].GridList;
    }

  }
#endif

  if (( StaticLevelZero == 1 && level != 0 ) || StaticLevelZero == 0 ) {

  FastSiblingLocatorInitialize(&ChainingMesh, MetaData->TopGridRank,
			       MetaData->TopGridDims);
 
  /* Add all the grids to the chaining mesh. */

  if (dbx) fprintf(stderr, "EL: FSL AddGrid entry \n");
  for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
    Grids[grid1]->GridData->FastSiblingLocatorAddGrid(&ChainingMesh);

  if (dbx) fprintf(stderr, "EL: FSL AddGrid exit \n");
 
  /* For each grid, get a list of possible siblings from the chaining mesh. */
 
  for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
    if (Grids[grid1]->GridData->FastSiblingLocatorFindSiblings(
                              &ChainingMesh, &SiblingList[grid1],
			      MetaData->LeftFaceBoundaryCondition,
			      MetaData->RightFaceBoundaryCondition) == FAIL) {
      fprintf(stderr, "Error in grid->FastSiblingLocatorFindSiblings.\n");
      ENZO_FAIL("");
    }
 
  /* Clean up the chaining mesh. */
 
  FastSiblingLocatorFinalize(&ChainingMesh);

  }

  /* On the top grid, adjust the refine region so that only the finest
     particles are included.  We don't want the more massive particles
     to contaminate the high-resolution region. */

  if (RefineRegionAutoAdjust && level == 0)
    if (AdjustRefineRegion(LevelArray, MetaData) == FAIL) {
      fprintf(stderr, "Error in AdjustRefineRegion.\n");
      ENZO_FAIL("");
    }

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
 
  JBPERF_STOP("evolve-level-02"); // SetBoundaryConditions()

  /* Clear the boundary fluxes for all Grids (this will be accumulated over
     the subcycles below (i.e. during one current grid step) and used to by the
     current grid to correct the zones surrounding this subgrid (step #18). */
 
  TIME_MSG("after SetBoundaryConditions");
  JBPERF_START("evolve-level-03"); // ClearBoundaryFluxes()

  for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
    Grids[grid1]->GridData->ClearBoundaryFluxes();
 
  JBPERF_STOP("evolve-level-03"); // ClearBoundaryFluxes()

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
 
    JBPERF_START("evolve-level-05"); // compute number of subgrids

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
 
    JBPERF_STOP("evolve-level-05"); // compute number of subgrids
    TIME_MSG("Before subgrid fluxes");


    /* For each grid, create the subgrid list. */
 
    JBPERF_START("evolve-level-06"); // create subgrid list

    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
 
      /* Allocate the subgrid fluxes for this grid. */
 
      SubgridFluxesEstimate[grid1] = new fluxes *[NumberOfSubgrids[grid1]];
 
      for (subgrid = 0; subgrid < NumberOfSubgrids[grid1]; subgrid++)
	SubgridFluxesEstimate[grid1][subgrid] = NULL;
 
      /* Collect the flux data and store it in the newly minted fluxes.
	 Or rather that's what we should do.  Instead, we create fluxes one
	 by one in this awkward array of pointers to pointers.  This should be
	 changed so that all the routines take arrays of flux rather than
	 arrays of pointers to flux.  Dumb. */
 
      counter = 0;

      // Only allocate fluxes for local grids: saves a *lot* of storage

      if (MyProcessorNumber ==
          Grids[grid1]->GridData->ReturnProcessorNumber()) {
 
	NextGrid = Grids[grid1]->NextGridNextLevel;
	while (NextGrid != NULL) {
	  SubgridFluxesEstimate[grid1][counter] = new fluxes;
	  Grids[grid1]->GridData->ComputeRefinementFactors
	                              (NextGrid->GridData, RefinementFactors);
	  NextGrid->GridData->ReturnFluxDims
             (*(SubgridFluxesEstimate[grid1][counter++]), RefinementFactors);
	  NextGrid = NextGrid->NextGridThisLevel;
	}
 
	/* Add the external boundary of this subgrid to the subgrid list. This
	   makes it easy to keep adding up the fluxes of this grid, but we must
	   keep in mind that the last subgrid should be ignored elsewhere. */
 
	SubgridFluxesEstimate[grid1][counter] = new fluxes;
	Grids[grid1]->GridData->ComputeRefinementFactors
                                   (Grids[grid1]->GridData, RefinementFactors);
	Grids[grid1]->GridData->ReturnFluxDims
               (*(SubgridFluxesEstimate[grid1][counter]), RefinementFactors);

      }
 
    } // end loop over grids (create Subgrid list)
 
    JBPERF_STOP("evolve-level-06"); // create subgrid list
    TIME_MSG("After subgrid fluxes");


    /* ------------------------------------------------------- */
    /* Prepare the density field (including particle density). */
 
//  fprintf(stderr, "%"ISYM": EvolveLevel: Enter PrepareDensityField\n", MyProcessorNumber);
 
    JBPERF_START("evolve-level-07"); // PrepareDensityField()

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
 
    JBPERF_STOP("evolve-level-07"); // PrepareDensityField()

    /* ------------------------------------------------------- */
    /* Evolve all grids by timestep dtThisLevel. */
 
    TIME_MSG("EvolveLevel: before main loop");
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
 
      /* Call analysis routines. */
 
      JBPERF_START_LOW("evolve-level-08"); // Call analysis routines

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
 
      JBPERF_STOP_LOW("evolve-level-08"); // Call analysis routines

      /* Gravity: compute acceleration field for grid and particles. */
 
      JBPERF_START("evolve-level-09"); // Compute self-gravity acceleration


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
 
      JBPERF_STOP("evolve-level-09"); // Compute self-gravity acceleration


      /* Gravity: compute field due to preset sources. */
 
      JBPERF_START_LOW("evolve-level-10"); // ComputeAccelerationFieldExternal()
	  Grids[grid1]->GridData->ComputeAccelerationFieldExternal();
 
      JBPERF_STOP_LOW("evolve-level-10"); // ComputeAccelerationFieldExternal()

      /* Radiation Pressure: add to acceleration field */

#ifdef TRANSFER
	  Grids[grid1]->GridData->AddRadiationPressureAcceleration() == FAIL);
#endif /* TRANSFER */

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


    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
#endif //SAB.
      /* Copy current fields (with their boundaries) to the old fields
	  in preparation for the new step. */
 
      JBPERF_START("evolve-level-11"); // CopyBaryonFieldToOldBaryonField()

	if (Grids[grid1]->GridData->CopyBaryonFieldToOldBaryonField() == FAIL) {
	  fprintf(stderr, "Error in grid->CopyBaryonFieldToOldBaryonField.\n");
	  ENZO_FAIL("");
	}
 
      JBPERF_STOP("evolve-level-11"); // CopyBaryonFieldToOldBaryonField()

      /* Add RandomForcing fields to velocities after the copying of current
         fields to old. I also update the total energy accordingly here.
         It makes no sense to force on the very first time step. */
 
      JBPERF_START_LOW("evolve-level-12"); // AddRandomForcing()

      if (RandomForcing && MetaData->CycleNumber > 0) //AK
        if(Grids[grid1]->GridData->AddRandomForcing(&norm,
                                                   TopGridTimeStep) == FAIL)
          fprintf(stderr, "Error in AddRandomForcing.\n");
 
      JBPERF_STOP_LOW("evolve-level-12"); // AddRandomForcing()

      /* Call hydro solver and save fluxes around subgrids. */
 
      JBPERF_START("evolve-level-13"); // SolveHydroEquations()

//      fprintf(stderr, "%"ISYM": Calling Hydro\n", MyProcessorNumber);
 
      if (Grids[grid1]->GridData->SolveHydroEquations(LevelCycleCount[level],
	 NumberOfSubgrids[grid1], SubgridFluxesEstimate[grid1], level) == FAIL) {
	fprintf(stderr, "Error in grid->SolveHydroEquations.\n");
	ENZO_FAIL("");
      }
 
      JBPERF_STOP("evolve-level-13"); // SolveHydroEquations()

//      fprintf(stderr, "%"ISYM": Called Hydro\n", MyProcessorNumber);
 
      /* Solve the radiative transfer */
	
#ifdef TRANSFER
      while ((dtPhoton > 0.) && RadiativeTransfer &&
	     (Grids[grid1]->GridData->ReturnTime() >= PhotonTime))  {
	if (debug) 
	  printf("EvolvePhotons[%"ISYM"]: dt = %"GSYM", Time = %"FSYM", ", 
		 level, dtPhoton, PhotonTime);
	if (EvolvePhotons(MetaData, LevelArray, AllStars) == FAIL) {
	  fprintf(stderr, "Error in EvolvePhotons.\n");
	  ENZO_FAIL("");
	}
      } /* ENDWHILE evolve photon */
#endif /* TRANSFER */

      /* Solve the cooling and species rate equations. */
 
//      fprintf(stderr, "%"ISYM": Calling SolveCoolAndRateEquations\n", MyProcessorNumber);

      if (MultiSpecies && RadiativeCooling) {
 
	JBPERF_START("evolve-level-14"); // change this?

	if (Grids[grid1]->GridData->SolveRateAndCoolEquations() == FAIL) {
	  fprintf(stderr, "Error in grid->SolveRateEquations.\n");
	  ENZO_FAIL("");
	}
 
	JBPERF_STOP("evolve-level-14"); // change this?

//      fprintf(stderr, "%"ISYM": Called SolveCoolAndRateEquations\n", MyProcessorNumber);

      } else {

//      fprintf(stderr, "%"ISYM": Calling MultiSpecies\n", MyProcessorNumber);
 
	JBPERF_START("evolve-level-14"); // SolveRateEquations()

	if (MultiSpecies)
	  if (Grids[grid1]->GridData->SolveRateEquations() == FAIL) {
	    fprintf(stderr, "Error in grid->SolveRateEquations.\n");
	    ENZO_FAIL("");
	  }
 
	JBPERF_STOP("evolve-level-14"); // SolveRateEquations()

//      fprintf(stderr, "%"ISYM": Called MultiSpecies\n", MyProcessorNumber);
 
	/* Include radiative cooling/heating. */
 
//      fprintf(stderr, "%"ISYM": Calling RadiativeCooling\n", MyProcessorNumber);
 
	JBPERF_START("evolve-level-15"); // SolveRadiativeCooling()

	if (RadiativeCooling)
	  if (Grids[grid1]->GridData->SolveRadiativeCooling() == FAIL) {
	    fprintf(stderr, "Error in grid->SolveRadiativeCooling.\n");
	    ENZO_FAIL("");
	  }
 
	JBPERF_STOP("evolve-level-15"); // SolveRadiativeCooling()

//      fprintf(stderr, "%"ISYM": Called RadiativeCooling\n", MyProcessorNumber);

      }

      /* Update particle positions (if present). */
 
//      fprintf(stderr, "%"ISYM": Calling UpdatePP\n", MyProcessorNumber);
 
      JBPERF_START("evolve-level-16"); // UpdateParticlePositions()

      if (UpdateParticlePositions(Grids[grid1]->GridData) == FAIL) {
	fprintf(stderr, "Error in UpdateParticlePositions.\n");
	ENZO_FAIL("");
      }
 
      JBPERF_STOP("evolve-level-16"); // UpdateParticlePositions()

//      fprintf(stderr, "%"ISYM": Called UpdatePP\n", MyProcessorNumber);
 
      /* Include 'star' particle creation and feedback.
         (first, set the under_subgrid field). */
 
      JBPERF_START_LOW("evolve-level-17"); // star particle creation/feedback

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
 
      JBPERF_STOP_LOW("evolve-level-17"); // star particle creation/feedback

      /* Gravity: clean up AccelerationField. */

      JBPERF_START_LOW("evolve-level-18"); // clean up AccelerationField

      if (SelfGravity || UniformGravity || PointSourceGravity) {
	if (level != MaximumGravityRefinementLevel ||
	    MaximumGravityRefinementLevel == MaximumRefinementLevel)
	  Grids[grid1]->GridData->DeleteAccelerationField();
	Grids[grid1]->GridData->DeleteParticleAcceleration();
      }
 
      JBPERF_STOP_LOW("evolve-level-18"); // clean up AccelerationField

      /* Update current problem time of this subgrid. */
 
      JBPERF_START_LOW("evolve-level-19"); // SetTimeNextTimestep()

      Grids[grid1]->GridData->SetTimeNextTimestep();
 
      JBPERF_STOP_LOW("evolve-level-19"); // SetTimeNextTimestep()

      /* If using comoving co-ordinates, do the expansion terms now. */
 
      JBPERF_START("evolve-level-20"); // ComovingExpansionTerms()

      if (ComovingCoordinates)
	Grids[grid1]->GridData->ComovingExpansionTerms();
 
      JBPERF_STOP("evolve-level-20"); // ComovingExpansionTerms()

    }  // end loop over grids
 
    /* For each grid: a) interpolate boundaries from the parent grid.
                      b) copy any overlapping zones from siblings. */
 
    JBPERF_START("evolve-level-21"); // SetBoundaryConditions()
    TIME_MSG("EvolveLevel: after main loop");

#ifdef FAST_SIB
    if (SetBoundaryConditions(Grids, NumberOfGrids, SiblingList,
			      level, MetaData, Exterior, LevelArray[level]) == FAIL)
      ENZO_FAIL("");
#else
    if (SetBoundaryConditions(Grids, NumberOfGrids, level, MetaData,
                              Exterior, LevelArray[level]) == FAIL)
      ENZO_FAIL("");
#endif

    JBPERF_STOP("evolve-level-21"); // SetBoundaryConditions()
    TIME_MSG("after SetBoundaryConditions");

    /* Finalize (accretion, feedback, etc.) star particles */
 
    JBPERF_START("evolve-level-22"); // StarParticleFinalize()

    if (StarParticleFinalize(Grids, MetaData, NumberOfGrids, LevelArray,
			     level, AllStars) == FAIL) {
      fprintf(stderr, "Error in StarParticleFinalize.\n");
      ENZO_FAIL("");
    }

    JBPERF_STOP("evolve-level-22"); // StarParticleFinalize()

    /* Check for movie output (only check if this is bottom of hierarchy). */
 
    JBPERF_START("evolve-level-23"); // WriteMovieData()

    if (LevelArray[level+1] == NULL)
      if (LevelArray[level]->GridData->ReturnTime() >=
	  MetaData->TimeLastMovieDump + MetaData->dtMovieDump &&
	  MetaData->dtMovieDump > 0.0) {
	MetaData->TimeLastMovieDump += MetaData->dtMovieDump;
	if (WriteMovieData(MetaData->MovieDumpName,
			  MetaData->MovieDumpNumber++, LevelArray, MetaData,
			  LevelArray[level]->GridData->ReturnTime()) == FAIL) {
	  fprintf(stderr, "Error in WriteMovieData.\n");
	  ENZO_FAIL("");
	}
      }
 
    JBPERF_STOP("evolve-level-23"); // WriteMovieData()

    /* Check for tracer particle output (only if this bottom of hierarchy). */
 
    JBPERF_START("evolve-level-24"); // WriteTracerParticleData()

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
	  fprintf(stderr, "Error in WriteTracerParticleData.\n");
	  ENZO_FAIL("");
	}
      }
 
    JBPERF_STOP("evolve-level-24"); // WriteTracerParticleData()

    /* If cosmology, then compute grav. potential for output if needed. */
 
    JBPERF_START("evolve-level-25"); // PrepareDensityField()

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
 
    JBPERF_STOP("evolve-level-25"); // PrepareDensityField()

    /* Check for new level output (only if this is bottom of hierarchy). */
 
    JBPERF_START("evolve-level-26"); // WriteAllData()

    if (MetaData->OutputFirstTimeAtLevel > 0 &&
	level >= MetaData->OutputFirstTimeAtLevel &&
	LevelArray[level+1] == NULL) {

      MetaData->OutputFirstTimeAtLevel = level+1;
      LevelHierarchyEntry *Temp2 = LevelArray[0];

      while (Temp2->NextGridThisLevel != NULL)
	Temp2 = Temp2->NextGridThisLevel; /* ugh: find last in linked list */

#ifdef USE_HDF5_GROUPS
      if (Group_WriteAllData(MetaData->DataDumpName, MetaData->DataDumpNumber++,
		       Temp2->GridHierarchyEntry, *MetaData, Exterior,
		       LevelArray[level]->GridData->ReturnTime()) == FAIL) {
	fprintf(stderr, "Error in Group_WriteAllData.\n");
	ENZO_FAIL("");
      }
#else
      if (WriteAllData(MetaData->DataDumpName, MetaData->DataDumpNumber++,
		       Temp2->GridHierarchyEntry, *MetaData, Exterior, 
		       LevelArray[level]->GridData->ReturnTime()) == FAIL) {
	fprintf(stderr, "Error in WriteAllData.\n");
	ENZO_FAIL("");
      }
#endif

    }
 
    /* Check for new output from new file in file system, but only
       when at bottom of hierarchy */
    
#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if ( (access("outputNow", F_OK) != -1 ) &&
    LevelArray[level+1] == NULL) {
#ifdef USE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      LevelHierarchyEntry *Temp2 = LevelArray[0];
      while (Temp2->NextGridThisLevel != NULL)
    Temp2 = Temp2->NextGridThisLevel; /* ugh: find last in linked list */
      printf("Detected outputNow.\n");
#ifdef USE_HDF5_GROUPS
      if (Group_WriteAllData(MetaData->DataDumpName, MetaData->DataDumpNumber++,
		       Temp2->GridHierarchyEntry, *MetaData, Exterior,
		       LevelArray[level]->GridData->ReturnTime()) == FAIL) {
	fprintf(stderr, "Error in Group_WriteAllData.\n");
	ENZO_FAIL("");
      }
#else
      if (WriteAllData(MetaData->DataDumpName, MetaData->DataDumpNumber++,
		       Temp2->GridHierarchyEntry, *MetaData, Exterior, 
		       LevelArray[level]->GridData->ReturnTime()) == FAIL) {
	fprintf(stderr, "Error in WriteAllData.\n");
	ENZO_FAIL("");
      }
#endif
      if (MyProcessorNumber == ROOT_PROCESSOR)
      if (unlink("outputNow")) {
    fprintf(stderr, "Error deleting 'outputNow'\n");
    ENZO_FAIL("");
      }
    } 
      
    /* Check to see if new subcycle information has been given to us */
    
#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if ( (access("subcycleCount", F_OK) != -1 ) &&
        LevelArray[level+1] == NULL) {
      printf("Detected subcycleCount\n");
#ifdef USE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      FILE *fptr; 
      if ((fptr = fopen("subcycleCount", "r")) == NULL) {
        fprintf(stderr, "Error opening subcycle file subcycleCount.  Continuing.\n");
      }
      else {
        /* Grab the number of cycles to dump on */
        char line[MAX_LINE_LENGTH];
        if (fgets(line, MAX_LINE_LENGTH, fptr) == NULL) {
          fprintf(stderr, "Error reading subcycle file subcycleCount.  Skipping.\n");
        } else {
          sscanf(line, "%"ISYM, &MetaData->SubcycleSkipDataDump);
          MetaData->SubcycleLastDataDump = MetaData->SubcycleNumber;
        }
        fclose(fptr);
#ifdef USE_MPI
        MPI_Barrier(MPI_COMM_WORLD);
#endif
        if (MyProcessorNumber == ROOT_PROCESSOR)
          if (unlink("subcycleCount")) {
            fprintf(stderr, "Error deleting subcycleCount.\n");
          }
      }
    }

    /* Check to see if we should start outputting interpolated data based on
       the cycles of the highest level */

    if (MetaData->SubcycleNumber >= MetaData->SubcycleLastDataDump +
        MetaData->SubcycleSkipDataDump   &&
        MetaData->SubcycleSkipDataDump > 0) {
      printf("Writing data based on SubcycleDumpSkipping (%"ISYM" %"ISYM" %"ISYM")\n",
          MetaData->SubcycleNumber, MetaData->SubcycleLastDataDump,
          MetaData->SubcycleSkipDataDump);
      MetaData->SubcycleLastDataDump += MetaData->SubcycleSkipDataDump;
#ifdef USE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      LevelHierarchyEntry *Temp2 = LevelArray[0];
      while (Temp2->NextGridThisLevel != NULL)
        Temp2 = Temp2->NextGridThisLevel; /* ugh: find last in linked list */

#ifdef USE_HDF5_GROUPS
      if (Group_WriteAllData(MetaData->DataDumpName, MetaData->DataDumpNumber++,
            Temp2->GridHierarchyEntry, *MetaData, Exterior,
            LevelArray[level]->GridData->ReturnTime()) == FAIL) {
        fprintf(stderr, "Error in Group_WriteAllData.\n");
        ENZO_FAIL("");
      }
#else
      if (WriteAllData(MetaData->DataDumpName, MetaData->DataDumpNumber++,
            Temp2->GridHierarchyEntry, *MetaData, Exterior,
            LevelArray[level]->GridData->ReturnTime()) == FAIL) {
        fprintf(stderr, "Error in WriteAllData.\n");
        ENZO_FAIL("");
      }
#endif

    } 

    /* Check for stop by file-touching */

#ifdef USE_MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
    if ( (access("stopNow", F_OK) != -1 ) &&
    LevelArray[level+1] == NULL) {
#ifdef USE_MPI
      MPI_Barrier(MPI_COMM_WORLD);
#endif
      LevelHierarchyEntry *Temp2 = LevelArray[0];
      while (Temp2->NextGridThisLevel != NULL)
    Temp2 = Temp2->NextGridThisLevel; /* ugh: find last in linked list */
      printf("Detected stopNow\n");
#ifdef USE_HDF5_GROUPS
      if (Group_WriteAllData(MetaData->DataDumpName, MetaData->DataDumpNumber++,
		       Temp2->GridHierarchyEntry, *MetaData, Exterior,
		       LevelArray[level]->GridData->ReturnTime()) == FAIL) {
	fprintf(stderr, "Error in Group_WriteAllData.\n");
	ENZO_FAIL("");
      }
#else
      if (WriteAllData(MetaData->DataDumpName, MetaData->DataDumpNumber++,
		       Temp2->GridHierarchyEntry, *MetaData, Exterior, 
		       LevelArray[level]->GridData->ReturnTime()) == FAIL) {
	fprintf(stderr, "Error in WriteAllData.\n");
	ENZO_FAIL("");
      }
#endif
      if (MyProcessorNumber == ROOT_PROCESSOR)
    if (unlink("stopNow")) {
      fprintf(stderr, "Error deleting stopNow\n");
      ENZO_FAIL("");
    } 
      fprintf(stderr, "Stopping due to request on level %"ISYM"\n", level);
      my_exit(EXIT_SUCCESS);
    }


    JBPERF_STOP("evolve-level-26"); // WriteAllData()

    /* Check for stop (unpleasant to exit from here, but...). */
 
    if (MetaData->StopFirstTimeAtLevel > 0 &&
	level >= MetaData->StopFirstTimeAtLevel &&
	LevelArray[level+1] == NULL) {

      // Write movie data in all grids if necessary

      if (MovieSkipTimestep != INT_UNDEFINED)
	for (int mlevel = 0; mlevel < MAX_DEPTH_OF_HIERARCHY; mlevel++) {
	  if (LevelArray[mlevel] == NULL) break;
	  delete [] Grids;
	  NumberOfGrids = GenerateGridArray(LevelArray, mlevel, &Grids);
	  if (WriteStreamData(Grids, NumberOfGrids, MetaData,
			      LevelCycleCount[mlevel], TRUE) == FAIL) {
	    fprintf(stderr, "Error in WriteStreamData.\n");
	    ENZO_FAIL("");
	  }
	}

      fprintf(stderr, "Stopping due to request on level %"ISYM"\n", level);
      my_exit(EXIT_SUCCESS);
    }
 
    /* For each grid, delete the GravitatingMassFieldParticles. */
 
    JBPERF_START("evolve-level-27"); // DeleteGravitatingMassFieldParticles()

    for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
      Grids[grid1]->GridData->DeleteGravitatingMassFieldParticles();
 
    JBPERF_STOP("evolve-level-27"); // DeleteGravitatingMassFieldParticles()

    /* Update SubcycleNumber if this is the bottom of the hierarchy --
       Note that this not unique based on which level is the highest,
       it just keeps going */

    if (LevelArray[level+1] == NULL) MetaData->SubcycleNumber += 1;  

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

#if defined(USE_JBPERF) && defined(JB_PERF_LEVELS)
    jbPerf.attribute ("level",&jb_level,JB_INT);
#endif

    // Streaming movie output (only run if everything is evolved)

    if (MovieSkipTimestep != INT_UNDEFINED) {
      if (WriteStreamData(Grids, NumberOfGrids, MetaData, 
			  LevelCycleCount[level]) == FAIL) {
	fprintf(stderr, "Error in WriteStreamData.\n");
	ENZO_FAIL("");
      }
    }

    if (dbx) fprintf(stderr, "EL Level %"ISYM" returns from Level %"ISYM"\n", level, level+1);

    /* ------------------------------------------------------- */
    /* For each grid,
     * (a) project the subgrid's solution into this grid (step #18)
     * (b) correct for the difference between this grid's fluxes and the
     *     subgrid's fluxes. (step #19)
     */
 
    JBPERF_START("evolve-level-28"); // UpdateFromFinerGrids()
    TIME_MSG("Before update from finer grids");

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

    JBPERF_STOP("evolve-level-28"); // UpdateFromFinerGrids()

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
 
    JBPERF_START("evolve-level-29"); // AddToBoundaryFluxes()
    TIME_MSG("Before saving fluxes");

    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {

      // Only deallocate fluxes for local grids

      if (MyProcessorNumber ==
          Grids[grid1]->GridData->ReturnProcessorNumber()) {
 
      if (FluxCorrection)
	if (Grids[grid1]->GridData->AddToBoundaryFluxes
	    (SubgridFluxesEstimate[grid1][NumberOfSubgrids[grid1] - 1])
	    == FAIL) {
	  fprintf(stderr, "Error in grid->AddToBoundaryFluxes.\n");
	  ENZO_FAIL("");
	}
 
      /* Delete fluxes pointed to by SubgridFluxesEstimate[subgrid]. */
 
      for (subgrid = 0; subgrid < NumberOfSubgrids[grid1]; subgrid++) {
	DeleteFluxes(SubgridFluxesEstimate[grid1][subgrid]);
	delete       SubgridFluxesEstimate[grid1][subgrid];
      }
      delete [] SubgridFluxesEstimate[grid1];

      }
 
    } // end of loop over grids
 
    JBPERF_STOP("evolve-level-29"); // AddToBoundaryFluxes()

    /* Recompute radiation field, if requested. */
 
    JBPERF_START("evolve-level-30"); // RadiationFieldUpdate()
    TIME_MSG("Done saving fluxes");

    if (RadiationFieldType >= 10 && RadiationFieldType <= 11 &&
	level <= RadiationFieldLevelRecompute)
      if (RadiationFieldUpdate(LevelArray, level, MetaData) == FAIL) {
	fprintf(stderr, "Error in RecomputeRadiationField.\n");
	ENZO_FAIL("");
      }
 
    JBPERF_STOP("evolve-level-30"); // RadiationFieldUpdate()

    /* Rebuild the Grids on the next level down.
       Don't bother on the last cycle, as we'll rebuild this grid soon. */
 
    JBPERF_START("evolve-level-31"); // RebuildHierarchy()

    if (dtThisLevelSoFar < dtLevelAbove) {
      if (RebuildHierarchy(MetaData, LevelArray, level) == FAIL) {
	fprintf(stderr, "Error in RebuildHierarchy.\n");
	ENZO_FAIL("");
      }
    }

    /* Count up number of grids on this level. */
 
    int GridMemory, NumberOfCells, CellsTotal, Particles;
    float AxialRatio, GridVolume;
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
      Grids[grid1]->GridData->CollectGridInformation
        (GridMemory, GridVolume, NumberOfCells, AxialRatio, CellsTotal, Particles);
      LevelZoneCycleCount[level] += NumberOfCells;
      if (MyProcessorNumber == Grids[grid1]->GridData->ReturnProcessorNumber())
	LevelZoneCycleCountPerProc[level] += NumberOfCells;
    }
 
    JBPERF_STOP("evolve-level-31"); // RebuildHierarchy()

    cycle++;
    LevelCycleCount[level]++;
 
  } // end of loop over subcycles
 
  if (debug)
    fprintf(stdout, "EvolveLevel[%"ISYM"]: NumberOfSubCycles = %"ISYM" (%"ISYM" total)\n", level,
           cycle, LevelCycleCount[level]);
 
  /* If possible & desired, report on memory usage. */
 
  //  if (debug)
  ReportMemoryUsage("Memory usage report: Evolve Level");
 
#if defined(USE_JBPERF) && defined(JB_PERF_LEVELS)
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
