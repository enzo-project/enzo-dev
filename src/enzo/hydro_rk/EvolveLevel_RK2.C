/***********************************************************************
/
/  EVOLVE LEVEL USING RUNGE-KUTTA2 METHOD
/
/  written by: Peng Wang
/  date:       May, 2007
/  modified1:
/
/  PURPOSE:
/         Evolve Hydro & MHD using 2nd order Runge-Kutta method.
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#include "performance.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "Grid.h"
#include "LevelHierarchy.h"
#include "../hydro_rk/tools.h"

/* function prototypes */

int StarParticleInitialize(LevelHierarchyEntry *LevelArray[], int ThisLevel,
			   TopGridData *MetaData, Star *&AllStars);
int StarParticleFinalize(HierarchyEntry *Grids[], TopGridData *MetaData,
			 int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
			 int level, Star *&AllStars);

#ifdef TRANSFER
int EvolvePhotons(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[],
		  Star *AllStars, FLOAT GridTime, int level, int LoopTime = TRUE);
int RadiativeTransferPrepare(LevelHierarchyEntry *LevelArray[], int level,
			     TopGridData *MetaData, Star *&AllStars,
			     float dtLevelAbove);
#endif

int CreateSiblingList(HierarchyEntry ** Grids, int NumberOfGrids, 
		      SiblingGridList *SiblingList, int StaticLevelZero, 
		      TopGridData* MetaData, int level);
void DeleteFluxes(fluxes *Fluxes);
int  RebuildHierarchy(TopGridData *MetaData,
		      LevelHierarchyEntry *LevelArray[], int level);
int  ReportMemoryUsage(char *header = NULL);
int  UpdateParticlePositions(grid *Grid);
int  CheckEnergyConservation(HierarchyEntry *Grids[], int grid, 
			     int NumberOfGrids, int level, float dt);
int CommunicationMergeStarParticle(HierarchyEntry *Grids[], int NumberOfGrids);
#ifdef USE_MPI
int CommunicationReduceValues(float *Values, int Number, MPI_Op ReduceOperation);
#endif
int CommunicationAllSumValues(float *Values, int Number);
float CommunicationMinValue(float Value);
float CommunicationMaxValue(float Value);
int CommunicationBarrier();
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

int SetBoundaryConditions(HierarchyEntry *Grids[], int NumberOfGrids,
			  SiblingGridList SiblingList[],
			  int level, TopGridData *MetaData, 
			  ExternalBoundary *Exterior);
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




int OutputFromEvolveLevel(LevelHierarchyEntry *LevelArray[],TopGridData *MetaData,
		      int level, ExternalBoundary *Exterior);

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

#ifdef FLUX_FIX
int CreateSUBlingList(TopGridData *MetaData,
		      HierarchyEntry *Grids[],
		      int NumberOfGrids,
		      LevelHierarchyEntry ***SUBlingList);
int DeleteSUBlingList(int NumberOfGrids,
		      LevelHierarchyEntry **SUBlingList);
#endif

// int UpdateFromFinerGrids(HierarchyEntry *Grids[], int NumberOfGrids,
// 			 int NumberOfSubgrids[], 
// 			 fluxes **SubgridFluxesEstimate[]);
int CommunicationUpdateStarParticleCount(HierarchyEntry *Grids[],
					 TopGridData *MetaData,
					 int NumberOfGrids);
int FinalizeFluxes(HierarchyEntry *Grids[],fluxes **SubgridFluxesEstimate[],
		 int NumberOfGrids,int NumberOfSubgrids[]);		 
int RadiationFieldUpdate(LevelHierarchyEntry *LevelArray[], int level,
			 TopGridData *MetaData);
//int WriteStreamData(LevelHierarchyEntry *LevelArray[], int level,
//                    TopGridData *MetaData, int *CycleCount, int open=FALSE);
//int WriteMovieData(char *basename, int filenumber, 
//		   LevelHierarchyEntry *LevelArray[], TopGridData *MetaData, 
//		   FLOAT WriteTime);
int WriteTracerParticleData(char *basename, int filenumber, 
		   LevelHierarchyEntry *LevelArray[], TopGridData *MetaData, 
		   FLOAT WriteTime);
int FastSiblingLocatorInitialize(ChainingMeshStructure *Mesh, int Rank,
                                 int TopGridDims[]);
int FastSiblingLocatorFinalize(ChainingMeshStructure *Mesh);
int FastSiblingLocatorInitializeStaticChainingMesh(ChainingMeshStructure *Mesh, int Rank,
						   int TopGridDims[]); 

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
double ReturnWallTime();
int CallPython(LevelHierarchyEntry *LevelArray[], TopGridData *MetaData,
               int level);
int SetLevelTimeStep(HierarchyEntry *Grids[], int NumberOfGrids, int level, 
		     float *dtThisLevelSoFar, float *dtThisLevel, 
		     float dtLevelAbove);

void my_exit(int status);
int CallPython(LevelHierarchyEntry *LevelArray[], TopGridData *MetaData,
               int level);

/* Counters for performance and cycle counting. */

static int LevelCycleCount[MAX_DEPTH_OF_HIERARCHY];
static int MovieCycleCount[MAX_DEPTH_OF_HIERARCHY];
//double LevelWallTime[MAX_DEPTH_OF_HIERARCHY];
//double LevelZoneCycleCount[MAX_DEPTH_OF_HIERARCHY];
//double LevelZoneCycleCountPerProc[MAX_DEPTH_OF_HIERARCHY];
 
static float norm = 0.0;            //AK
static float TopGridTimeStep = 0.0; //AK

static int StaticSiblingListInitialized = 0;

#ifdef STATIC_SIBLING_LIST
static SiblingGridList StaticSiblingList[MAX_NUMBER_OF_SUBGRIDS];
static int StaticLevelZero = 1;
#else
static int StaticLevelZero = 0;
#endif

/* EvolveGrid function */

int EvolveLevel_RK2(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[],
		    int level, float dtLevelAbove, ExternalBoundary *Exterior, FLOAT dt0)
{

  float dtThisLevelSoFar = 0.0, dtThisLevel, dtGrid;
  int RefinementFactors[MAX_DIMENSION];
  int cycle = 0, counter = 0, grid1, subgrid, iLevel;
  HierarchyEntry *NextGrid;
  double time1 = ReturnWallTime();

  // Update lcaperf "level" attribute

  Eint32 jb_level = level;
#ifdef USE_JBPERF
  jbPerf.attribute ("level",&jb_level,JB_INT);
#endif

#ifdef FLUX_FIX
  /* Create a SUBling list of the subgrids */
  LevelHierarchyEntry **SUBlingList;
#endif

  FLOAT When;
  float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1, TimeUnits, 
    VelocityUnits, CriticalDensity = 1, BoxLength = 1, MagneticUnits;
  double MassUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, 1.0);
  MassUnits = DensityUnits*pow(LengthUnits,3);
  
  /* Create an array (Grids) of all the grids. */

  typedef HierarchyEntry* HierarchyEntryPointer;
  HierarchyEntry **Grids;
  int NumberOfGrids = GenerateGridArray(LevelArray, level, &Grids);
  int *NumberOfSubgrids = new int[NumberOfGrids];
  fluxes ***SubgridFluxesEstimate = new fluxes **[NumberOfGrids];

  SiblingGridList *SiblingList = new SiblingGridList[NumberOfGrids];
  CreateSiblingList(Grids, NumberOfGrids, SiblingList, StaticLevelZero, 
		    MetaData, level);

  /*if (SetBoundaryConditions(Grids, NumberOfGrids,SiblingList,level, MetaData,
			      Exterior) == FAIL) {
    ENZO_FAIL("");
    }*/
#ifdef FAST_SIB
  if (SetBoundaryConditions(Grids, NumberOfGrids, SiblingList,
			    level, MetaData, Exterior, LevelArray[level]) == FAIL)
    ENZO_FAIL("");
#else
  if (SetBoundaryConditions(Grids, NumberOfGrids, level, MetaData,
                            Exterior, LevelArray[level]) == FAIL)
    ENZO_FAIL("");
#endif


  
  /* Count the number of colours in the first grid (to define Ncolor) */

  Grids[0]->GridData->SetNumberOfColours();

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

 
#ifdef TRANSFER
    RadiativeTransferPrepare(LevelArray, level, MetaData, AllStars, dtLevelAbove);
#endif /* TRANSFER */


    /* For each grid, compute the number of it's subgrids. */

    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
      NextGrid = Grids[grid1]->NextGridNextLevel;
      counter = 0;
      while (NextGrid != NULL) {
	NextGrid = NextGrid->NextGridThisLevel;
	if (++counter > MAX_NUMBER_OF_SUBGRIDS)
	  ENZO_FAIL("More subgrids than MAX_NUMBER_OF_SUBGRIDS.");
      }
      NumberOfSubgrids[grid1] = counter + 1;
    }


    /* For each grid, create the subgrid list. */

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

      } // end: (is this grid on this processor)

    } // end loop over grids (create Subgrid list)

    /* compute wave speed
       Reference: Matsumoto, PASJ, 2007, 59, 905 */

    if (HydroMethod == MHD_RK) {
      int lmax;
      LevelHierarchyEntry *Temp;
      for (lmax = MAX_DEPTH_OF_HIERARCHY-1; lmax >= 0; lmax--) {
	Temp = LevelArray[lmax];
	if (Temp != NULL) {
	  break;
	} 
      }
      //      lmax = 0; // <- Pengs version had lmax = 6
      //      lmax = 1;
      FLOAT dx0, dy0, dz0, h_min, DivBDampingLength = 1.0;

      dx0 = (DomainRightEdge[0] - DomainLeftEdge[0]) / MetaData->TopGridDims[0];
      dy0 = (MetaData->TopGridRank > 1) ? 
	(DomainRightEdge[1] - DomainLeftEdge[1]) / MetaData->TopGridDims[1] : 1e8;
      dz0 = (MetaData->TopGridRank > 2) ? 
	(DomainRightEdge[2] - DomainLeftEdge[2]) / MetaData->TopGridDims[2] : 1e8;
      h_min = my_MIN(dx0, dy0, dz0);
      h_min /= pow(RefineBy, lmax);
      C_h = 0.5*MetaData->CourantSafetyNumber*(h_min/dt0);
      C_h = min( C_h, 1e6/VelocityUnits); // never faster than __ cm/s (for very small dt0 a problems)
      if (EOSType == 3)  // for isothermal runs just use the constant sound speed
	C_h = EOSSoundSpeed;

      C_p = sqrt(0.18*DivBDampingLength*C_h);
      //      C_p = sqrt(0.18*DivBDampingLength)*C_h;
      if (debug) 
	fprintf(stderr, "lengthscale %g timestep: %g  C_h: %g  C_p: %g\n ", 
		h_min, dt0, C_h, C_p);
    }

//     if (SelfGravity && MetaData->TopGridRank == 3) {
//       if (PrepareDensityField(LevelArray, SiblingList, level, MetaData) == FAIL) {
//       //      if (PrepareDensityField(LevelArray, level, MetaData) == FAIL) {
// 	fprintf(stderr, "Error in PrepareDensityField.\n");
// 	ENZO_FAIL("");
//       }
//     }

    When = 0.5;
    if (SelfGravity) {
#ifdef FAST_SIB
      PrepareDensityField(LevelArray, SiblingList, level, MetaData, When);
#else   // !FAST_SIB
      PrepareDensityField(LevelArray, level, MetaData, When);
#endif  // end FAST_SIB
    }

    /* Solve the radiative transfer */

#ifdef TRANSFER
    FLOAT GridTime = Grids[0]->GridData->ReturnTime();
    EvolvePhotons(MetaData, LevelArray, AllStars, GridTime, level);
#endif /* TRANSFER */

    /* Compute particle-particle acceleration */

    if (NBodyDirectSummation == TRUE) 
      for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
	Grids[grid1]->GridData->ComputeParticleParticleAcceleration(level);

    /* ------------------------------------------------------- */
    /* Evolve all grids by timestep dtThisLevel. Predictor Step*/

    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
      
      CallProblemSpecificRoutines(MetaData, Grids[grid1], grid1, &norm, 
				  TopGridTimeStep, level, LevelCycleCount);

      /* Gravity: compute acceleration field for grid and particles. */
      if (SelfGravity) {
	int Dummy;
	if (level <= MaximumGravityRefinementLevel) {
	  if (level > 0) 
	    Grids[grid1]->GridData->SolveForPotential(level) ;
	  Grids[grid1]->GridData->ComputeAccelerations(level) ;
	}
	// otherwise, use interpolated potential from coarser grid, which is
	// done in PrepareDensity.
      } // end: if (SelfGravity)

      Grids[grid1]->GridData->ComputeAccelerationFieldExternal() ;

#ifdef TRANSFER
      /* Radiation Pressure: add to acceleration field */
      Grids[grid1]->GridData->AddRadiationPressureAcceleration();
#endif /* TRANSFER */

      Grids[grid1]->GridData->CopyBaryonFieldToOldBaryonField();

      if (UseHydro) {
	if (HydroMethod == HD_RK)
	  Grids[grid1]->GridData->RungeKutta2_1stStep
	    (SubgridFluxesEstimate[grid1], NumberOfSubgrids[grid1], level, Exterior);
	else if (HydroMethod == MHD_RK) {
	  Grids[grid1]->GridData->MHDRK2_1stStep
	    (SubgridFluxesEstimate[grid1], NumberOfSubgrids[grid1], level, Exterior);
	}
      } // ENDIF UseHydro
	
      /* Do this here so that we can get the correct
	 time interpolated boundary condition */
      Grids[grid1]->GridData->SetTimeNextTimestep();
      
    }  // end loop over grids
      
    /*if (SetBoundaryConditions(Grids, NumberOfGrids,SiblingList,level, MetaData,
      Exterior) == FAIL) {
      ENZO_FAIL("");
      }*/
#ifdef FAST_SIB
    SetBoundaryConditions(Grids, NumberOfGrids, SiblingList, level, MetaData, Exterior, LevelArray[level]);
#else
    SetBoundaryConditions(Grids, NumberOfGrids, level, MetaData, Exterior, LevelArray[level]);
#endif



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
	
	 
	  time1 = ReturnWallTime();
	  
	  Grids[grid1]->GridData->PoissonSolver(level);
	
	} // ENDIF MHD_RK
      } // ENDIF UseHydro

      /* Add viscosity */

      if (UseViscosity) {
	Grids[grid1]->GridData->AddViscosity();
      }

      /* Solve the cooling and species rate equations. */
 
      Grids[grid1]->GridData->MultiSpeciesHandler();

      /* Update particle positions (if present). */

      UpdateParticlePositions(Grids[grid1]->GridData);

      /* Include 'star' particle creation and feedback. */

      Grids[grid1]->GridData->StarParticleHandler
	(Grids[grid1]->NextGridNextLevel, level);
 
      /* Gravity: clean up AccelerationField. */

	 if (level != MaximumGravityRefinementLevel ||
	     MaximumGravityRefinementLevel == MaximumRefinementLevel)
	     Grids[grid1]->GridData->DeleteAccelerationField();

      Grids[grid1]->GridData->DeleteParticleAcceleration();
 
      if (UseFloor) 
	Grids[grid1]->GridData->SetFloor();
      

      /* If using comoving co-ordinates, do the expansion terms now. */

      if (ComovingCoordinates)
	Grids[grid1]->GridData->ComovingExpansionTerms();

    }  // end loop over grids


#ifdef FAST_SIB
    SetBoundaryConditions(Grids, NumberOfGrids, SiblingList, level, 
			  MetaData, Exterior, LevelArray[level]);
#else
    SetBoundaryConditions(Grids, NumberOfGrids, level, MetaData, 
			  Exterior, LevelArray[level]);
#endif



    /* Finalize (accretion, feedback, etc.) star particles */
 
    StarParticleFinalize(Grids, MetaData, NumberOfGrids, LevelArray,
			 level, AllStars);

    /* Reompute potential if output is requested */
    if (SelfGravity && WritePotential) {
      CopyGravPotential = TRUE;
      When = 0.0;
 
#ifdef FAST_SIB
      PrepareDensityField(LevelArray, SiblingList, level, MetaData, When);
#else   // !FAST_SIB
      PrepareDensityField(LevelArray, level, MetaData, When);
#endif  // end FAST_SIB
 
      //      CopyGravPotential = FALSE;
 
      for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
        if (level <= MaximumGravityRefinementLevel) {
 
          /* Compute the potential. */
 
          if (level > 0)
            Grids[grid1]->GridData->SolveForPotential(level);
          Grids[grid1]->GridData->CopyPotentialToBaryonField();
        }
      } //  end loop over grids
    } // if WritePotential


    OutputFromEvolveLevel(LevelArray,MetaData,level,Exterior);
    CallPython(LevelArray, MetaData, level);

    /* For each grid, delete the GravitatingMassFieldParticles. */

    for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
      Grids[grid1]->GridData->DeleteGravitatingMassFieldParticles();


    /* ----------------------------------------- */
    /* Evolve the next level down (recursively). */

    MetaData->FirstTimestepAfterRestart = FALSE;

    //    LevelWallTime[level] += ReturnWallTime() - time1;
    if (LevelArray[level+1] != NULL) {
      if (EvolveLevel_RK2(MetaData, LevelArray, level+1, dtThisLevel, Exterior, dt0) 
	  == FAIL) {
	fprintf(stderr, "Error in EvolveLevel_RK2 (%d).\n", level);
	ENZO_FAIL("");
      }
    }
    time1 = ReturnWallTime();

#ifdef USE_JBPERF
    // Update lcaperf "level" attribute

    jbPerf.attribute ("level",&jb_level,JB_INT);
#endif

    OutputFromEvolveLevel(LevelArray,MetaData,level,Exterior);
    CallPython(LevelArray, MetaData, level);

    /* Update SubcycleNumber and the timestep counter for the
       streaming data if this is the bottom of the hierarchy -- Note
       that this not unique based on which level is the highest, it
       just keeps going */

    if (LevelArray[level+1] == NULL) {
      MetaData->SubcycleNumber++;
      MetaData->TimestepCounter++;
    }

    /* ------------------------------------------------------- */
    /* For each grid,
     (a) project the subgrid's solution into this grid (step #18)
     (b) correct for the difference between this grid's fluxes and the
         subgrid's fluxes. (step #19) */


#ifdef FLUX_FIX
    SUBlingList = new LevelHierarchyEntry*[NumberOfGrids];
    CreateSUBlingList(MetaData, Grids,NumberOfGrids, &SUBlingList);
#endif

#ifdef FLUX_FIX
    UpdateFromFinerGrids(level, Grids, NumberOfGrids, NumberOfSubgrids,
			 SubgridFluxesEstimate,
			 SUBlingList,
			 MetaData);
#else
    UpdateFromFinerGrids(level, Grids, NumberOfGrids, NumberOfSubgrids, SubgridFluxesEstimate);
#endif
    
#ifdef FLUX_FIX        /* Clean up SUBlings */
    DeleteSUBlingList( NumberOfGrids, SUBlingList );
#endif

     
    /* ------------------------------------------------------- */
    /* Add the saved fluxes (in the last subsubgrid entry) to the exterior
       fluxes for this subgrid .
       (Note: this must be done after CorrectForRefinedFluxes). */

    FinalizeFluxes(Grids,SubgridFluxesEstimate,NumberOfGrids,NumberOfSubgrids);
    
    /* Recompute radiation field, if requested. */

    RadiationFieldUpdate(LevelArray, level, MetaData);

    /* Rebuild the Grids on the next level down.
       Don't bother on the last cycle, as we'll rebuild this grid soon. */

    //    LevelWallTime[level] += ReturnWallTime() - time1;
    if (dtThisLevelSoFar < dtLevelAbove) 
      RebuildHierarchy(MetaData, LevelArray, level);

    time1 = ReturnWallTime();

    /* Count up number of grids on this level. */

    int GridMemory, NumberOfCells, CellsTotal, Particles;
    float AxialRatio, GridVolume;
#ifdef UNUSED
    for (grid1 = 0; grid1 < NumberOfGrids; grid1++) {
      Grids[grid1]->GridData->CollectGridInformation
        (GridMemory, GridVolume, NumberOfCells, AxialRatio, CellsTotal, Particles);
      //      LevelZoneCycleCount[level] += NumberOfCells;
      //if (MyProcessorNumber == Grids[grid1]->GridData->ReturnProcessorNumber())
      //	LevelZoneCycleCountPerProc[level] += NumberOfCells;
    }
#endif

    cycle++;
    LevelCycleCount[level]++;

  } // while (dtThisLevelSoFar < dtLevelAbove)


  if (debug)
    printf("EvolveLevelRK2[%d]: NumberOfSubCycles = %d (%d total)\n", level, 
           cycle, LevelCycleCount[level]);

  /* If possible & desired, report on memory usage. */

  ReportMemoryUsage("Memory usage report: Evolve Level");

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
    for (int grid1 = 0; grid1 < NumberOfGrids; grid1++)
      delete [] SiblingList[grid1].GridList;
    delete [] SiblingList;
  }

  return SUCCESS;

}
