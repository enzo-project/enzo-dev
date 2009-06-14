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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <unistd.h>
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
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
int AdjustRefineRegion(LevelHierarchyEntry *LevelArray[], 
		       TopGridData *MetaData);

#ifdef TRANSFER
int EvolvePhotons(TopGridData *MetaData,LevelHierarchyEntry *LevelArray[],
		  Star *AllStars);
int RadiativeTransferPrepare(LevelHierarchyEntry *LevelArray[], int level,
			     TopGridData *MetaData, Star *&AllStars,
			     float dtLevelAbove);
#endif

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
int PrepareDensityField(LevelHierarchyEntry *LevelArray[],
                        SiblingGridList SiblingList[],
                        int level, TopGridData *MetaData);
int PrepareDensityField(LevelHierarchyEntry *LevelArray[],
                        int level, TopGridData *MetaData, double When);
int SetBoundaryConditions(HierarchyEntry *Grids[], int NumberOfGrids,
			  SiblingGridList SiblingList[],
			  int level, TopGridData *MetaData, 
			  ExternalBoundary *Exterior);
#ifdef SIB2
int SetBoundaryConditions(HierarchyEntry *Grids[], int NumberOfGrids,
			  SiblingGridList SiblingList[],
			  int level, TopGridData *MetaData,
			  ExternalBoundary *Exterior, LevelHierarchyEntry * Level);
#else
int SetBoundaryConditions(HierarchyEntry *Grids[], int NumberOfGrids,
                          int level, TopGridData *MetaData,
                          ExternalBoundary *Exterior, LevelHierarchyEntry * Level);
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
#ifdef USE_HDF5_GROUPS
int Group_WriteAllData(char *basename, int filenumber, HierarchyEntry *TopGrid,
		 TopGridData &MetaData, ExternalBoundary *Exterior,
		 FLOAT WriteTime = -1);
#else
int WriteAllData(char *basename, int filenumber, HierarchyEntry *TopGrid,
		 TopGridData &MetaData, ExternalBoundary *Exterior,
		 FLOAT WriteTime = -1);
#endif
int FastSiblingLocatorInitialize(ChainingMeshStructure *Mesh, int Rank,
                                 int TopGridDims[]);
int FastSiblingLocatorFinalize(ChainingMeshStructure *Mesh);
int FastSiblingLocatorInitializeStaticChainingMesh(ChainingMeshStructure *Mesh, int Rank,
						   int TopGridDims[]); 

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
double ReturnWallTime();
void my_exit(int status);

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
  int cycle = 0, counter = 0, grid, subgrid, iLevel, ErrorSignal = 0;
  HierarchyEntry *NextGrid;
  double time1 = ReturnWallTime();
  Star *AllStars = NULL;

#ifdef FLUX_FIX
  /* Create a SUBling list of the subgrids */
  LevelHierarchyEntry **SUBlingList;
#endif

  if (StarParticleCreation || StarParticleFeedback)
    if (StarParticleInitialize(LevelArray, level, MetaData, 
			       AllStars) == FAIL) {
      fprintf(stderr, "Error in StarParticleInitalize.\n");
      return FAIL;
    }

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


#ifdef STATIC_SIBLING_LIST
  if ( StaticLevelZero == 1 && level == 0 ) {

    if (!StaticSiblingListInitialized) {

      if (debug) fprintf(stderr, "INITIALIZE Level 0 StaticSiblingList\n");

      ChainingMeshStructure StaticChainingMesh;

      FastSiblingLocatorInitializeStaticChainingMesh
	(&StaticChainingMesh, MetaData->TopGridRank, MetaData->TopGridDims);

      for (grid = 0; grid < NumberOfGrids; grid++)
        Grids[grid]->GridData->FastSiblingLocatorAddGrid(&StaticChainingMesh);

      for (grid = 0; grid < NumberOfGrids; grid++)
        if (Grids[grid]->GridData->FastSiblingLocatorFindSiblings(
                              &StaticChainingMesh, &StaticSiblingList[grid],
                              MetaData->LeftFaceBoundaryCondition,
                              MetaData->RightFaceBoundaryCondition) == FAIL) {
          fprintf(stderr, "Error in grid->FastSiblingLocatorFindSiblings.\n");
          return FAIL;
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

    for (grid = 0; grid < NumberOfGrids; grid++) {
      SiblingList[grid].NumberOfSiblings = StaticSiblingList[grid].NumberOfSiblings;
      SiblingList[grid].GridList = StaticSiblingList[grid].GridList;
    }

  }
#endif

  if (( StaticLevelZero == 1 && level != 0 ) || StaticLevelZero == 0 ) {

  FastSiblingLocatorInitialize(&ChainingMesh, MetaData->TopGridRank,
			       MetaData->TopGridDims);
 
  /* Add all the grids to the chaining mesh. */

  for (grid = 0; grid < NumberOfGrids; grid++)
    Grids[grid]->GridData->FastSiblingLocatorAddGrid(&ChainingMesh);

 
  /* For each grid, get a list of possible siblings from the chaining mesh. */
 
  for (grid = 0; grid < NumberOfGrids; grid++)
    if (Grids[grid]->GridData->FastSiblingLocatorFindSiblings(
                              &ChainingMesh, &SiblingList[grid],
			      MetaData->LeftFaceBoundaryCondition,
			      MetaData->RightFaceBoundaryCondition) == FAIL) {
      fprintf(stderr, "Error in grid->FastSiblingLocatorFindSiblings.\n");
      return FAIL;
    }
 
  /* Clean up the chaining mesh. */
 
  FastSiblingLocatorFinalize(&ChainingMesh);

  }

  //  PerformanceTimers[31] += ReturnWallTime() - time1;


  /*if (SetBoundaryConditions(Grids, NumberOfGrids,SiblingList,level, MetaData,
			      Exterior) == FAIL) {
    return FAIL;
    }*/
#ifdef SIB2
  if (SetBoundaryConditions(Grids, NumberOfGrids, SiblingList,
			    level, MetaData, Exterior, LevelArray[level]) == FAIL)
    return FAIL;
#else
  if (SetBoundaryConditions(Grids, NumberOfGrids, level, MetaData,
                            Exterior, LevelArray[level]) == FAIL)
    return FAIL;
#endif



  for (grid = 0; grid < NumberOfGrids; grid++)
    Grids[grid]->GridData->ClearBoundaryFluxes();

  /* Create a list of shining (no radiative transfer, only a 1/r^2
     radiation profile) particles from all grids in all levels.  For
     now, this only applies to Pop III star particles. Then check if
     the stellar feedback is contained in grids on this level and
     finer grids.  If so, apply changes to the grid(s). */

#ifdef TRANSFER
    if (RadiativeTransfer)
      if (RadiativeTransferPrepare(LevelArray, level, MetaData, AllStars, 
				   dtLevelAbove) == FAIL) {
	fprintf(stderr, "Error in RadiativeTransferPrepare.\n");
	return FAIL;
      }
#endif /* TRANSFER */


  /* ================================================================== */
  /* Loop over grid timesteps until the elapsed time equals the timestep 
     from the level above (or loop once for the top level). */

  while (dtThisLevelSoFar < dtLevelAbove) {

    if (level == 0) {

      dtThisLevel      = dtLevelAbove;
      dtThisLevelSoFar = dtLevelAbove;

    } else {

      dtThisLevel = huge_number;
      for (grid = 0; grid < NumberOfGrids; grid++) {
	dtGrid      = Grids[grid]->GridData->ComputeTimeStep();
	dtThisLevel = min(dtThisLevel, dtGrid);
      }
      dtThisLevel = CommunicationMinValue(dtThisLevel);

      /* Advance dtThisLevelSoFar (don't go over dtLevelAbove). */
      
      if (dtThisLevelSoFar+dtThisLevel*1.05 >= dtLevelAbove) {
	dtThisLevel      = dtLevelAbove - dtThisLevelSoFar;
	dtThisLevelSoFar = dtLevelAbove;
      }
      else
	dtThisLevelSoFar += dtThisLevel;

    }

    if (debug) {
      float utime = 1.0;
      if (UsePhysicalUnit) {
	utime = TimeUnits/3.1558e7;
      }
      printf("Level[%d]: dt = %g(%g/%g)\n", level, dtThisLevel*utime,
	     dtThisLevelSoFar*utime, dtLevelAbove*utime);
    }

    for (grid = 0; grid < NumberOfGrids; grid++)
      Grids[grid]->GridData->SetTimeStep(dtThisLevel);

//     /* Streaming movie output (write before everything is evolved) */

//     if (MovieSkipTimestep != INT_UNDEFINED) {
//       if (WriteStreamData(LevelArray, level, MetaData, MovieCycleCount) == FAIL) {
//         fprintf(stderr, "Error in WriteStreamData.\n");
//         return FAIL;
//       }
//       if (MovieCycleCount[level] == MovieSkipTimestep)
//         for (iLevel = level; iLevel < MAX_DEPTH_OF_HIERARCHY; iLevel++)
//           MovieCycleCount[iLevel] = 0;
//       MovieCycleCount[level]++;
//     }


    /* For each grid, compute the number of it's subgrids. */

    for (grid = 0; grid < NumberOfGrids; grid++) {
      NextGrid = Grids[grid]->NextGridNextLevel;
      counter = 0;
      while (NextGrid != NULL) {
	NextGrid = NextGrid->NextGridThisLevel;
	if (++counter > MAX_NUMBER_OF_SUBGRIDS) {
	  fprintf(stderr, "More subgrids than MAX_NUMBER_OF_SUBGRIDS.\n");
	  return FAIL;
	}
      }
      NumberOfSubgrids[grid] = counter + 1;
    }



    /* For each grid, create the subgrid list. */

    for (grid = 0; grid < NumberOfGrids; grid++) {

      /* Allocate the subgrid fluxes for this grid. */

      SubgridFluxesEstimate[grid] = new fluxes *[NumberOfSubgrids[grid]];
      for (subgrid = 0; subgrid < NumberOfSubgrids[grid]; subgrid++)
	SubgridFluxesEstimate[grid][subgrid] = NULL;

      /* Collect the flux data and store it in the newly minted fluxes.
	 Or rather that's what we should do.  Instead, we create fluxes one
	 by one in this awkward array of pointers to pointers.  This should be
	 changed so that all the routines take arrays of flux rather than
	 arrays of pointers to flux.  Dumb. */

      counter = 0;
      if (MyProcessorNumber == 
	  Grids[grid]->GridData->ReturnProcessorNumber()) {
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
	   makes it easy to keep adding up the fluxes of this grid, but we must
	   keep in mind that the last subgrid should be ignored elsewhere. */

	SubgridFluxesEstimate[grid][counter] = new fluxes;
	Grids[grid]->GridData->ComputeRefinementFactors
                                   (Grids[grid]->GridData, RefinementFactors);
	Grids[grid]->GridData->ReturnFluxDims
               (*(SubgridFluxesEstimate[grid][counter]), RefinementFactors);

      } // end: (is this grid on this processor)

    } // end loop over grids (create Subgrid list)

    /* compute wave speed
       Reference: Matsumoto, PASJ, 2006 */

    if (HydroMethod == MHD_RK) {
      int lmax;
      LevelHierarchyEntry *Temp;
      for (lmax = MAX_DEPTH_OF_HIERARCHY-1; lmax >= 0; lmax--) {
	Temp = LevelArray[lmax];
	if (Temp != NULL) {
	  break;
	} 
      }
      lmax = 6;
      FLOAT dx0 = (DomainRightEdge[0] - DomainLeftEdge[0]) / MetaData->TopGridDims[0];
      FLOAT dy0 = (MetaData->TopGridRank > 1) ? 
	(DomainRightEdge[1] - DomainLeftEdge[1]) / MetaData->TopGridDims[1] : 1e8;
      FLOAT dz0 = (MetaData->TopGridRank > 2) ? 
	(DomainRightEdge[2] - DomainLeftEdge[2]) / MetaData->TopGridDims[2] : 1e8;
      FLOAT h_min = my_MIN(dx0, dy0, dz0);
      h_min /= pow(RefineBy, lmax);
      C_h = MetaData->CourantSafetyNumber*h_min/dt0;
      C_p = sqrt(0.18*C_h);

    }

//     if (SelfGravity && MetaData->TopGridRank == 3) {
//       if (PrepareDensityField(LevelArray, SiblingList, level, MetaData) == FAIL) {
//       //      if (PrepareDensityField(LevelArray, level, MetaData) == FAIL) {
// 	fprintf(stderr, "Error in PrepareDensityField.\n");
// 	return FAIL;
//       }
//     }

    When = 0.5;

#ifdef FAST_SIB
    if (SelfGravity)
      if (PrepareDensityField(LevelArray, SiblingList,
			      level, MetaData, When) == FAIL) {
	fprintf(stderr, "Error in PrepareDensityField.\n");
	return FAIL;
      }
#else   // !SIB3
    if (SelfGravity)
      if (PrepareDensityField(LevelArray, level, MetaData, When) == FAIL) {
        fprintf(stderr, "Error in PrepareDensityField.\n");
        return FAIL;
      }
#endif  // end SIB3


    /* Compute particle-particle acceleration */

    //if (NBodyDirectSummation == TRUE) 
    //for (grid = 0; grid < NumberOfGrids; grid++)
    //Grids[grid]->GridData->ComputeParticleParticleAcceleration(level);

    /* ------------------------------------------------------- */
    /* Evolve all grids by timestep dtThisLevel. */

    for (grid = 0; grid < NumberOfGrids; grid++) {

      /* Gravity: compute acceleration field for grid and particles. */

      if (SelfGravity && MetaData->TopGridRank == 3) {
	int Dummy;
	if (level <= MaximumGravityRefinementLevel) {
	  if (level > 0) {
	    if (Grids[grid]->GridData->SolveForPotential(Dummy, level) 
		== FAIL) {
	      fprintf(stderr, "Error in grid->SolveForPotential.\n");
	      return FAIL;
	    }
	  }
	  if (Grids[grid]->GridData->ComputeAccelerations(level) == FAIL) {
	    fprintf(stderr, "Error in grid->ComputeAccelerations.\n");
	    return FAIL;
	  }
	}
	// otherwise, interpolate potential from coarser grid, which is
	//   now done in PrepareDensity.
      } // end: if (SelfGravity)

      if (UniformGravity || PointSourceGravity || ExternalGravity) {
	if (Grids[grid]->GridData->ComputeAccelerationFieldExternal() ==FAIL) {
	  fprintf(stderr,"Error in grid->ComputeAccelerationFieldExternal.\n");
	  return FAIL;
	}
      }

#ifdef TRANSFER

      /* Radiation Pressure: add to acceleration field */

      if (RadiativeTransfer && RadiationPressure)
	if (Grids[grid]->GridData->AddRadiationPressureAcceleration() == FAIL) {
	  fprintf(stderr,"Error in grid->AddRadiationPressureAcceleration.\n");
	  return FAIL;
	}

#endif /* TRANSFER */


      if (Grids[grid]->GridData->CopyBaryonFieldToOldBaryonField() == FAIL) {
	fprintf(stderr, "Error in grid->CopyBaryonFieldToOldBaryonField.\n");
	return FAIL;
      }

      if (UseHydro) {
	if (HydroMethod == HD_RK) {
	  if (Grids[grid]->GridData->RungeKutta2_1stStep(LevelCycleCount[level], 
							 SubgridFluxesEstimate[grid],
							 NumberOfSubgrids[grid], level,
							 Exterior) == FAIL) {
	    fprintf(stderr, "Error in grid->RungeKutta2_1stStep.\n");
	    return FAIL;
	  }
	} 
	else if (HydroMethod == MHD_RK) {
	  if (Grids[grid]->GridData->MHDRK2_1stStep(LevelCycleCount[level], 
						    SubgridFluxesEstimate[grid],
						    NumberOfSubgrids[grid], level,
						    Exterior) == FAIL) {
	    fprintf(stderr, "Error in grid->MHDRK2_1stStep.\n");
	    return FAIL;
	  }
	}
      }
	
      /* Do this here so that we can get the correct
	 time interpolated boundary condition */
      Grids[grid]->GridData->SetTimeNextTimestep();
      
    }  // end loop over grids
      
    /*if (SetBoundaryConditions(Grids, NumberOfGrids,SiblingList,level, MetaData,
      Exterior) == FAIL) {
      return FAIL;
      }*/
#ifdef SIB2
  if (SetBoundaryConditions(Grids, NumberOfGrids, SiblingList,
			    level, MetaData, Exterior, LevelArray[level]) == FAIL)
    return FAIL;
#else
  if (SetBoundaryConditions(Grids, NumberOfGrids, level, MetaData,
                            Exterior, LevelArray[level]) == FAIL)
    return FAIL;
#endif
    
    for (grid = 0; grid < NumberOfGrids; grid++) {

      if (UseHydro) {
	if (HydroMethod == HD_RK) {
	  if (Grids[grid]->GridData->RungeKutta2_2ndStep(LevelCycleCount[level], 
							 SubgridFluxesEstimate[grid],
							 NumberOfSubgrids[grid], level,
							 Exterior) == FAIL) {
	    fprintf(stderr, "Error in grid->RungeKutta2_2ndStep.\n");
	    ErrorSignal = 1;
	    continue;
	  }
	}
	else if (HydroMethod == MHD_RK) {
	  if (Grids[grid]->GridData->MHDRK2_2ndStep(LevelCycleCount[level], 
						    SubgridFluxesEstimate[grid],
						    NumberOfSubgrids[grid], level,
						    Exterior) == FAIL) {
	    fprintf(stderr, "Error in grid->MHDRK2_2ndStep.\n");
	    return FAIL;
	  }

	  if (UseAmbipolarDiffusion) {
	    Grids[grid]->GridData->AddAmbipolarDiffusion();
	  }

	  if (UseResistivity) {
	    Grids[grid]->GridData->AddResistivity();
	  }
	
	  if(UseDivergenceCleaning){
	    
	    time1 = ReturnWallTime();

	    if (Grids[grid]->GridData->PoissonSolver(UseDivergenceCleaning, level)==FAIL){
	      fprintf(stderr, "Error in grid->PoissonSolver.\n");
	      return FAIL;
	    }

	    if (Grids[grid]->GridData->PoissonCleanStep(level)==FAIL){
	      fprintf(stderr, "Error in grid->PoissonCleaning.\n");
	      return FAIL;
	    }
	    
	    //	    PerformanceTimers[32] += ReturnWallTime() - time1;

	  }

	}
      }

      /* Add viscosity */

      if (UseViscosity) {
	Grids[grid]->GridData->AddViscosity();
      }

      /* Solve the chemical and species rate equations. */
      
      if (MultiSpecies && RadiativeCooling) {
	
	if (Grids[grid]->GridData->SolveRateAndCoolEquations() == FAIL) {
	  fprintf(stderr, "Error in grid->SolveRateAndCoolEquations\n");
	  ErrorSignal = 1;
	  continue;
	}
	
	/* 3: Cooling power modulated by atomic hydrogen density
	   4: Cooling function modulated by atomic hydrogen density */
	
	if (RadiativeCooling == 3 || RadiativeCooling == 4) {
	  if (Grids[grid]->GridData->SolveRadiativeCooling() == FAIL) {
	    fprintf(stderr, "Error in grid->SolveRadiativeCooling.\n");
	    ErrorSignal = 1;
	    continue;
	  }
	}
	
      } 
      else {
	
	if (MultiSpecies) {
	  if (Grids[grid]->GridData->SolveRateEquations() == FAIL) {
	    fprintf(stderr, "Error in grid->SolveRateEquations.\n");
	    ErrorSignal = 1;
	    continue;
	  }
	}
	
	/* Equilibrium radiative cooling/heating. */
	
	if (RadiativeCooling) {
	  if (Grids[grid]->GridData->SolveRadiativeCooling() == FAIL) {
	    fprintf(stderr, "Error in grid->SolveRadiativeCooling.\n");
	    ErrorSignal = 1;
	    continue;
	  }
	}
	
      } /* ENDELSE MultiSpecies && RadiativeCooling */


      /* Update particle positions (if present). */

      if (UpdateParticlePositions(Grids[grid]->GridData) == FAIL) {
	fprintf(stderr, "Error in UpdateParticlePositions.\n");
	ErrorSignal = 1;
	continue;
      }

      /* Include 'star' particle creation and feedback. */

      if (StarParticleCreation || StarParticleFeedback) {

	/* First, set the under_subgrid field. */

	Grids[grid]->GridData->ZeroSolutionUnderSubgrid(NULL, 
						 ZERO_UNDER_SUBGRID_FIELD);
	LevelHierarchyEntry *Temp2 = LevelArray[level+1];
	while (Temp2 != NULL) {
	  Grids[grid]->GridData->ZeroSolutionUnderSubgrid(Temp2->GridData, 
					 ZERO_UNDER_SUBGRID_FIELD);
	  Temp2 = Temp2->NextGridThisLevel;
	}

	/* Do star particle creation and feedback */


      Grids[grid]->GridData->StarParticleHandler
	(Grids[grid]->NextGridNextLevel, level);
      }

      if ((SelfGravity || UniformGravity || PointSourceGravity || ExternalGravity) 
	  && MetaData->TopGridRank == 3) {
	if (level != MaximumGravityRefinementLevel ||
	    MaximumGravityRefinementLevel == MaximumRefinementLevel) {
	  Grids[grid]->GridData->DeleteAccelerationField();
	}
	Grids[grid]->GridData->DeleteParticleAcceleration();
      }

      if (UseFloor) {
	Grids[grid]->GridData->SetFloor();
      }

      /* If using comoving co-ordinates, do the expansion terms now. */

      if (ComovingCoordinates)
	Grids[grid]->GridData->ComovingExpansionTerms();

    }  // end loop over grids


    /* Solve the radiative transfer */

#ifdef TRANSFER
    Grids[0]->GridData->SetTimePreviousTimestep();
      while ((dtPhoton > 0.) && RadiativeTransfer &&
	     (Grids[grid]->GridData->ReturnTime() >= PhotonTime))  {
	if (debug) 
	  printf("EvolvePhotons[%"ISYM"]: dt = %"GSYM", Time = %"FSYM", ", 
		 level, dtPhoton, PhotonTime);
	if (EvolvePhotons(MetaData, LevelArray, AllStars) == FAIL) {
	  fprintf(stderr, "Error in EvolvePhotons.\n");
	  return FAIL;
	}
      } /* ENDWHILE evolve photon */
    Grids[0]->GridData->SetTimeNextTimestep();
#endif /* TRANSFER */

    /*if (SetBoundaryConditions(Grids, NumberOfGrids,SiblingList,level, MetaData,
				Exterior) == FAIL) {
      return FAIL;
      }*/
#ifdef SIB2
  if (SetBoundaryConditions(Grids, NumberOfGrids, SiblingList,
			    level, MetaData, Exterior, LevelArray[level]) == FAIL)
    return FAIL;
#else
  if (SetBoundaryConditions(Grids, NumberOfGrids, level, MetaData,
                            Exterior, LevelArray[level]) == FAIL)
    return FAIL;
#endif

    
    CommunicationBarrier();

//     if (StarParticleCreation >> SINK_PARTICLE & 1 && level == MaximumRefinementLevel) {
//       if (CommunicationMergeStarParticle(Grids, NumberOfGrids) == FAIL) {
// 	printf("CommunicationMergeStarParticle failed.\n");
// 	return FAIL;
//       }
//     }

    CommunicationBarrier();

    if (StarParticleCreation) {
      if (CommunicationUpdateStarParticleCount(Grids, MetaData,
					       NumberOfGrids) == FAIL) {
	return FAIL;
      }
    }
    
    /* Collect all sink particle masses and report them to STDOUT */

    CommunicationBarrier();

    if (StarParticleCreation >> SINK_PARTICLE & 1 &&
	level == MaximumRefinementLevel) {
      float totalMass = 0;
      for (int ilevel = 0; ilevel <= MaximumRefinementLevel; ilevel++) {
	LevelHierarchyEntry *Temp2 = LevelArray[ilevel];
	while (Temp2 != NULL) {
	  totalMass += Temp2->GridData->ReturnTotalSinkMass();
	  Temp2 = Temp2->NextGridThisLevel;
	}
      }
      CommunicationAllSumValues(&totalMass, 1);

// #ifdef UNUSED
//       /* If crossed the critical mass, make a new shining particle */
//       float ShiningCriticalMass = 100.0;
//       float massdiff = 500.0 - ShiningCriticalMass;
//       if (RadiativeTransfer == 1) { 
// 	if ((int)((totalMass*MassUnits/1.989e33+massdiff)/500.0) != 
// 	    (int)((TotalSinkMass*MassUnits/1.989e33+massdiff)/500.0)) {
	  
// 	  printf("Finding new shining particle...\n");
	  
// 	  /* Find out the current maximum non-shining particle */
	  
// 	  float maxMass1 = -0.1, maxMass2;
// 	  LevelHierarchyEntry *maxGrid;
// 	  for (int ilevel = 0; ilevel <= MaximumRefinementLevel; ilevel++) {
// 	    LevelHierarchyEntry *Temp2 = LevelArray[ilevel];
// 	    while (Temp2 != NULL) {
// 	      maxMass2 = Temp2->GridData->ReturnMaximumNonRadiatingSinkMass();
// 	      if (maxMass2 > maxMass1) {
// 		maxMass1 = maxMass2;
// 		maxGrid = Temp2;
// 	      }
// 	      Temp2 = Temp2->NextGridThisLevel;
// 	    }
// 	  }

// 	  float maxMass = CommunicationMaxValue(maxMass1);
	  
// 	  /* Set the current maximum mass non-shining particle to shining */
	  
// 	  if (maxMass == maxMass1) {
// 	    if (maxGrid->GridData->SetRadiatingSinkParticle() == FAIL) {
// 	      printf("EvolveLevel_RK2: SetRadiatingSinkParticle failed.\n");
// 	      return FAIL;
// 	    }
// 	  }
	  
// 	}
//       } // if (RadiativeTransfer
//       TotalSinkMass = totalMass;
// #endif /* UNUSED */

      FLOAT ThisTime = LevelArray[MaximumRefinementLevel]->GridData->ReturnTime();
      if (debug) {
	fprintf(stdout, "SinkParticle: Time, Total Mass = %"GOUTSYM" %g\n",
		ThisTime, totalMass*MassUnits/1.989e33);
      }
    }

    /* Check for movie output (only check if this is bottom of hierarchy). */

//     if (LevelArray[level+1] == NULL)
//       if (LevelArray[level]->GridData->ReturnTime() >=
// 	  MetaData->TimeLastMovieDump + MetaData->dtMovieDump &&
// 	  MetaData->dtMovieDump > 0.0) {
// 	MetaData->TimeLastMovieDump += MetaData->dtMovieDump;
// 	if (WriteMovieData(MetaData->MovieDumpName, 
// 			  MetaData->MovieDumpNumber++, LevelArray, MetaData, 
// 			  LevelArray[level]->GridData->ReturnTime()) == FAIL) {
// 	  fprintf(stderr, "Error in WriteMovieData.\n");
// 	  return FAIL;
// 	}
//       }

    /* Check for tracer particle output (only if this bottom of hierarchy). */

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
	  return FAIL;
	}
      }

    /* Communicate ErrorSignal to stop all processes and output data */

#ifdef USE_MPI
    int value = ErrorSignal;
    MPI_Allreduce(&value, &ErrorSignal, 1, MPI_INT, MPI_MAX, 
		  MPI_COMM_WORLD);
#endif /* USE_MPI */    

    /* Check for error signal and output if there is an error */
    CommunicationBarrier();
    if (ErrorSignal) {
      LevelHierarchyEntry *Temp2 = LevelArray[0];
      while (Temp2->NextGridThisLevel != NULL)
	Temp2 = Temp2->NextGridThisLevel; /* ugh: find last in linked list */
      
      fprintf(stderr, "Error in EvolveLevel.\n");
      fprintf(stderr, "--> Dumping data (output number %d).\n",
	      MetaData->DataDumpNumber);

#ifdef USE_HDF5_GROUPS
      if (Group_WriteAllData(MetaData->DataDumpName, MetaData->DataDumpNumber++,
		       Temp2->GridHierarchyEntry, *MetaData, Exterior,
		       LevelArray[level]->GridData->ReturnTime()) == FAIL) {
	fprintf(stderr, "Error in Group_WriteAllData.\n");
	return FAIL;
      }
#else
      if (WriteAllData(MetaData->DataDumpName, MetaData->DataDumpNumber++,
		       Temp2->GridHierarchyEntry, *MetaData, Exterior, 
		       LevelArray[level]->GridData->ReturnTime()) == FAIL) {
	fprintf(stderr, "Error in WriteAllData.\n");
	return FAIL;
      }
#endif
    }

    /* Check for new level output (only if this is bottom of hierarchy). */

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
	return FAIL;
      }
#else
      if (WriteAllData(MetaData->DataDumpName, MetaData->DataDumpNumber++,
		       Temp2->GridHierarchyEntry, *MetaData, Exterior, 
		       LevelArray[level]->GridData->ReturnTime()) == FAIL) {
	fprintf(stderr, "Error in WriteAllData.\n");
	return FAIL;
      }
#endif
    }

    /* Check for stop (unpleasant to exit from here, but...). */

    if (MetaData->StopFirstTimeAtLevel > 0 &&
	level >= MetaData->StopFirstTimeAtLevel &&
	LevelArray[level+1] == NULL) {
      fprintf(stderr, "Stopping due to request on level %d\n", level);

      my_exit(EXIT_SUCCESS);
    }

    /* For each grid, delete the GravitatingMassFieldParticles. */

    for (grid = 0; grid < NumberOfGrids; grid++)
      Grids[grid]->GridData->DeleteGravitatingMassFieldParticles();

    /* ----------------------------------------- */
    /* Evolve the next level down (recursively). */

    //    LevelWallTime[level] += ReturnWallTime() - time1;
    if (LevelArray[level+1] != NULL) {
      if (EvolveLevel_RK2(MetaData, LevelArray, level+1, dtThisLevel, Exterior, dt0) 
	  == FAIL) {
	fprintf(stderr, "Error in EvolveLevel (%d).\n", level);
	return FAIL;
      }
    }
    time1 = ReturnWallTime();


    /* ------------------------------------------------------- */
    /* For each grid,
     (a) project the subgrid's solution into this grid (step #18)
     (b) correct for the difference between this grid's fluxes and the
         subgrid's fluxes. (step #19) */


#ifdef FLUX_FIX
    SUBlingList = new LevelHierarchyEntry*[NumberOfGrids];
    for(int list=0; list < NumberOfGrids; list++)
      SUBlingList[list] = NULL;

    if (FluxCorrection) {
      /* Fill in the SUBling list */
      if (CreateSUBlingList(MetaData, Grids,
                              NumberOfGrids, &SUBlingList) == FAIL) {
        fprintf(stderr, "Error in CreateSUBlingList.\n");
        return FAIL;
      }
    }

#endif

#ifdef FLUX_FIX
    if (UpdateFromFinerGrids(level, Grids, NumberOfGrids, NumberOfSubgrids,
			     SubgridFluxesEstimate,
			     SUBlingList,
			     MetaData) == FAIL)
      return FAIL;
#else
    if (UpdateFromFinerGrids(level, Grids, NumberOfGrids, NumberOfSubgrids,
			     SubgridFluxesEstimate) == FAIL)
      return FAIL;
#endif

#ifdef FLUX_FIX
    if ( FluxCorrection ) {
      /* Clean up SUBlings */
      if (DeleteSUBlingList( NumberOfGrids, SUBlingList ) == FAIL) {
        fprintf(stderr, "Error in DeleteSUBlingList.\n");
        return FAIL;
      }
    }
#endif

     
    /* ------------------------------------------------------- */
    /* Add the saved fluxes (in the last subsubgrid entry) to the exterior
       fluxes for this subgrid .
       (Note: this must be done after CorrectForRefinedFluxes). */

    for (grid = 0; grid < NumberOfGrids; grid++) {

      if (FluxCorrection)
	if (Grids[grid]->GridData->AddToBoundaryFluxes(SubgridFluxesEstimate[grid][NumberOfSubgrids[grid] - 1])
	    == FAIL) {
	  fprintf(stderr, "Error in grid->AddToBoundaryFluxes.\n");
	  return FAIL;
	}
      /* Delete fluxes pointed to by SubgridFluxesEstimate[subgrid]. */

      for (subgrid = 0; subgrid < NumberOfSubgrids[grid]; subgrid++) {
	DeleteFluxes(SubgridFluxesEstimate[grid][subgrid]);
	delete SubgridFluxesEstimate[grid][subgrid];
      }

      delete [] SubgridFluxesEstimate[grid];

    } // end of loop over grids

    /* Recompute radiation field, if requested. */

    if (RadiationFieldType >= 10 && RadiationFieldType <= 11 && 
	level <= RadiationFieldLevelRecompute)
      if (RadiationFieldUpdate(LevelArray, level, MetaData) == FAIL) {
	fprintf(stderr, "Error in RecomputeRadiationField.\n");
	return FAIL;
      }

    /* Rebuild the Grids on the next level down.
       Don't bother on the last cycle, as we'll rebuild this grid soon. */

    //    LevelWallTime[level] += ReturnWallTime() - time1;
    if (dtThisLevelSoFar < dtLevelAbove) {
      if (RebuildHierarchy(MetaData, LevelArray, level) == FAIL) {
	fprintf(stderr, "Error in RebuildHierarchy.\n");
	return FAIL;
      }
    }
    time1 = ReturnWallTime();

    /* Count up number of grids on this level. */

    int GridMemory, NumberOfCells, CellsTotal, Particles;
    float AxialRatio, GridVolume;
    for (grid = 0; grid < NumberOfGrids; grid++) {
      Grids[grid]->GridData->CollectGridInformation
        (GridMemory, GridVolume, NumberOfCells, AxialRatio, CellsTotal, Particles);
      //      LevelZoneCycleCount[level] += NumberOfCells;
      //      if (MyProcessorNumber == Grids[grid]->GridData->ReturnProcessorNumber())
      //	LevelZoneCycleCountPerProc[level] += NumberOfCells;
    }


    cycle++;
    //    LevelCycleCount[level]++;

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
