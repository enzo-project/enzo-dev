/***********************************************************************
/
/  INITIALIZE A NEW SIMULATION
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  Robert Harkness
/  date:       September 2004
/  modified2:  Stephen Skory
/  date:       May, 2008
/  modified3:  Alexei Kritsuk
/  date:       May, 2008
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
// This routine intializes a new simulation based on the parameter file.
 
 
#include <string.h>
#include <stdio.h>
 
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
 
// Function prototypes
 
int InitializeMovieFile(TopGridData &MetaData, HierarchyEntry &TopGrid);
int WriteHierarchyStuff(FILE *fptr, HierarchyEntry *Grid,
                        char* base_name, int &GridID, FLOAT WriteTime);
int ReadParameterFile(FILE *fptr, TopGridData &MetaData, float *Initialdt);
int WriteParameterFile(FILE *fptr, TopGridData &MetaData);
void ConvertTotalEnergyToGasEnergy(HierarchyEntry *Grid);
int SetDefaultGlobalValues(TopGridData &MetaData);
int CommunicationPartitionGrid(HierarchyEntry *Grid, int gridnum);
int CommunicationBroadcastValue(PINT *Value, int BroadcastProcessor);
 
// Initialization function prototypes
 
int ShockTubeInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid);
int WavePoolInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
		       TopGridData &MetaData);
int ShockPoolInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
			TopGridData &MetaData);
int DoubleMachInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
			TopGridData &MetaData, ExternalBoundary &Exterior);
int ShockInABoxInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
			  TopGridData &MetaData, ExternalBoundary &Exterior);
int ImplosionInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
                        TopGridData &MetaData);
int RotatingCylinderInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
			       TopGridData &MetaData);
int KHInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
                          TopGridData &MetaData);
int NohInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
                          TopGridData &MetaData);
int SedovBlastInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
                          TopGridData &MetaData);
int RadiatingShockInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
			     TopGridData &MetaData);
int ZeldovichPancakeInitialize(FILE *fptr, FILE *Outfptr,
			       HierarchyEntry &TopGrid);
int PressurelessCollapseInitialize(FILE *fptr, FILE *Outfptr,
			       HierarchyEntry &TopGrid, TopGridData &MetaData);
int AdiabaticExpansionInitialize(FILE *fptr, FILE *Outfptr,
				 HierarchyEntry &TopGrid);
int TestGravityInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
			  TopGridData &MetaData);
int TestOrbitInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
                        TopGridData &MetaData);
int GalaxySimulationInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
                        TopGridData &MetaData);
int TestGravitySphereInitialize(FILE *fptr, FILE *Outfptr,
			       HierarchyEntry &TopGrid, TopGridData &MetaData);
int SphericalInfallInitialize(FILE *fptr, FILE *Outfptr,
			      HierarchyEntry &TopGrid, TopGridData &MetaData);
int GravityEquilibriumTestInitialize(FILE *fptr, FILE *Outfptr,
			      HierarchyEntry &TopGrid, TopGridData &MetaData);
int CollapseTestInitialize(FILE *fptr, FILE *Outfptr,
			  HierarchyEntry &TopGrid, TopGridData &MetaData);
int TestGravityMotion(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
			  TopGridData &MetaData);
int SupernovaRestartInitialize(FILE *fptr, FILE *Outfptr,
			       HierarchyEntry &TopGrid, TopGridData &MetaData,
			       ExternalBoundary &Exterior);
int PutSinkRestartInitialize(FILE *fptr, FILE *Outfptr,
			       HierarchyEntry &TopGrid, TopGridData &MetaData,
			       ExternalBoundary &Exterior);
int ProtostellarCollapseInitialize(FILE *fptr, FILE *Outfptr,
				   HierarchyEntry &TopGrid,
				   TopGridData &MetaData);
int CoolingTestInitialize(FILE *fptr, FILE *Outfptr, 
			  HierarchyEntry &TopGrid, TopGridData &MetaData);
int CoolingTestInitialize_BDS(FILE *fptr, FILE *Outfptr, 
			  HierarchyEntry &TopGrid, TopGridData &MetaData);
 
int CosmologySimulationInitialize(FILE *fptr, FILE *Outfptr,
                                  HierarchyEntry &TopGrid,
                                  TopGridData &MetaData);
int CosmologySimulationReInitialize(HierarchyEntry *TopGrid,
                                    TopGridData &MetaData);
 
int NestedCosmologySimulationInitialize(FILE *fptr, FILE *Outfptr,
                                        HierarchyEntry &TopGrid,
                                        TopGridData &MetaData);
int NestedCosmologySimulationReInitialize(HierarchyEntry *TopGrid,
                                          TopGridData &MetaData);
 
int TurbulenceSimulationInitialize(FILE *fptr, FILE *Outfptr,
                                  HierarchyEntry &TopGrid,
                                  TopGridData &MetaData);
int TurbulenceSimulationReInitialize(HierarchyEntry *TopGrid,
                                    TopGridData &MetaData);
 
int TracerParticleCreation(FILE *fptr, HierarchyEntry &TopGrid,
                           TopGridData &MetaData);

int ShearingBoxInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
                        TopGridData &MetaData);
int ShearingBox2DInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
                        TopGridData &MetaData);
int ShearingBoxStratifiedInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
                        TopGridData &MetaData);
#ifdef TRANSFER
int PhotonTestInitialize(FILE *fptr, FILE *Outfptr, 
			 HierarchyEntry &TopGrid, TopGridData &MetaData);
#endif /* TRANSFER */


int Hydro1DTestInitialize(FILE *fptr, FILE *Outfptr,
			  HierarchyEntry &TopGrid, TopGridData &MetaData);
int TurbulenceInitialize(FILE *fptr, FILE *Outfptr, 
			 HierarchyEntry &TopGrid, TopGridData &MetaData, int SetBaryonFields);
int Collapse3DInitialize(FILE *fptr, FILE *Outfptr,
			 HierarchyEntry &TopGrid, TopGridData &MetaData);
int Collapse1DInitialize(FILE *fptr, FILE *Outfptr,
			 HierarchyEntry &TopGrid, TopGridData &MetaData);
int MHD1DTestInitialize(FILE *fptr, FILE *Outfptr,
                        HierarchyEntry &TopGrid, TopGridData &MetaData);
int MHD2DTestInitialize(FILE *fptr, FILE *Outfptr,
                        HierarchyEntry &TopGrid, TopGridData &MetaData);
int MHD3DTestInitialize(FILE *fptr, FILE *Outfptr, 
			HierarchyEntry &TopGrid, TopGridData &MetaData);
int CollapseMHD3DInitialize(FILE *fptr, FILE *Outfptr, 
			    HierarchyEntry &TopGrid, TopGridData &MetaData);
int MHDTurbulenceInitialize(FILE *fptr, FILE *Outfptr, 
			    HierarchyEntry &TopGrid, TopGridData &MetaData, int SetBaryonFields);
int GalaxyDiskInitialize(FILE *fptr, FILE *Outfptr, 
			 HierarchyEntry &TopGrid, TopGridData &MetaData);
int AGNDiskInitialize(FILE *fptr, FILE *Outfptr, 
		      HierarchyEntry &TopGrid, TopGridData &MetaData);
int FreeExpansionInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
			    TopGridData &MetaData);

int PoissonSolverTestInitialize(FILE *fptr, FILE *Outfptr, 
				HierarchyEntry &TopGrid, TopGridData &MetaData);

void PrintMemoryUsage(char *str);



 
// Character strings
 
char outfilename[] = "amr.out";
 
 
 
 
int InitializeNew(char *filename, HierarchyEntry &TopGrid,
		  TopGridData &MetaData, ExternalBoundary &Exterior,
		  float *Initialdt)
{

  
 
  // Declarations
 
  FILE *fptr, *BCfptr, *Outfptr;
  float Dummy[MAX_DIMENSION];
  int dim, i;

 
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    Dummy[dim] = 0.0;
 
  // Open parameter file
 
  if ((fptr = fopen(filename, "r")) == NULL) {
        ENZO_FAIL("Error opening parameter file.");
  }
 
  // Open output file
 
  if (MyProcessorNumber == ROOT_PROCESSOR)
    if ((Outfptr = fopen(outfilename, "w")) == NULL) {
      fprintf(stderr, "Error opening parameter output file %s\n", outfilename);
      ENZO_FAIL("");
    }
 
  // set the default MetaData values
 
  SetDefaultGlobalValues(MetaData);
 
  // Read the MetaData/global values from the Parameter file
 
  if (ReadParameterFile(fptr, MetaData, Initialdt) == FAIL) {
        ENZO_FAIL("Error in ReadParameterFile.");
  }

  // Set the number of particle attributes, if left unset
 
  if (NumberOfParticleAttributes == INT_UNDEFINED)
    if (StarParticleCreation || StarParticleFeedback)
      NumberOfParticleAttributes = 3;
    else
      NumberOfParticleAttributes = 0;
 
  // Give unset parameters their default values
 
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    if (RefineRegionLeftEdge[dim] == FLOAT_UNDEFINED)
      RefineRegionLeftEdge[dim]   = DomainLeftEdge[dim];
    if (RefineRegionRightEdge[dim] == FLOAT_UNDEFINED)
      RefineRegionRightEdge[dim]  = DomainRightEdge[dim];
  }
 
  // If the problem reads in a restart dump, then skip over the following
 
  if (ProblemType != 40) {
 
  // Error check the rank
 
  if (MetaData.TopGridRank < 0 || MetaData.TopGridRank > 3) {
    fprintf(stderr, "TopGridRank = %"ISYM" ill defined.\n", MetaData.TopGridRank);
    ENZO_FAIL("");
  }
 
  // Error check the dimensions and at the same time add ghost zones
 
  for (dim = 0; dim < MetaData.TopGridRank; dim++) {
    if (MetaData.TopGridDims[dim] < 1 || MetaData.TopGridDims[dim] > 8192) {
      fprintf(stderr, "TopGridDims[%"ISYM"] = %"ISYM" ill defined.\n", dim,
	      MetaData.TopGridDims[dim]);
      ENZO_FAIL("");
    }
    MetaData.TopGridDims[dim] = (MetaData.TopGridDims[dim] > 1) ?
                     MetaData.TopGridDims[dim] + 2*DEFAULT_GHOST_ZONES : 1;
  }
 
  // Create the top grid, prepare it, set the time and parameters
 
  TopGrid.GridData = new grid;
 
  TopGrid.GridData->PrepareGrid(MetaData.TopGridRank, MetaData.TopGridDims,
				DomainLeftEdge, DomainRightEdge,
				MetaData.NumberOfParticles);
  TopGrid.GridData->SetTime(MetaData.Time);
  TopGrid.GridData->SetHydroParameters(MetaData.CourantSafetyNumber,
				       MetaData.PPMFlatteningParameter,
				       MetaData.PPMDiffusionParameter,
				       MetaData.PPMSteepeningParameter);
  TopGrid.GridData->SetGravityParameters(MetaData.GravityBoundary);
 
  // Repair TopGridDims (subtract ghost zones added earlier)
 
  for (dim = 0; dim < MetaData.TopGridRank; dim++)
    MetaData.TopGridDims[dim] = max(MetaData.TopGridDims[dim] -
				    2*DEFAULT_GHOST_ZONES, 1);
 
  // Set TopGrid Hierarchy Entry
 
  TopGrid.NextGridThisLevel = NULL;  // always true
  TopGrid.ParentGrid        = NULL;  // always true
  TopGrid.NextGridNextLevel = NULL;  // can be reset by initializer
 
  } // end: if (ProblemType != 40)
 
  // Call problem initializer

  PrintMemoryUsage("Call problem init");
 
  if (ProblemType == 0) {
        ENZO_FAIL("No problem specified.");
  }
 
  int ret = INT_UNDEFINED;
 
  if (debug)
    printf("InitializeNew: Starting problem initialization.\n");
 
  // 1) Shocktube problem
 
  if (ProblemType == 1)
    ret = ShockTubeInitialize(fptr, Outfptr, TopGrid);
 
  // 2) Wave pool
 
  if (ProblemType == 2)
    ret = WavePoolInitialize(fptr, Outfptr, TopGrid, MetaData);
 
  // 3) Shock pool
 
  if (ProblemType == 3)
    ret = ShockPoolInitialize(fptr, Outfptr, TopGrid, MetaData);
 
  // 4) Double Mach reflection
 
  if (ProblemType == 4)
    ret = DoubleMachInitialize(fptr, Outfptr, TopGrid, MetaData, Exterior);
 
  // 5) ShockInABox
 
  if (ProblemType == 5)
    ret = ShockInABoxInitialize(fptr, Outfptr, TopGrid, MetaData, Exterior);
 
  // 6) Implosion
 
  if (ProblemType == 6)
    ret = ImplosionInitialize(fptr, Outfptr, TopGrid, MetaData);
 
  // 7) SedovBlast
 
  if (ProblemType == 7)
    ret = SedovBlastInitialize(fptr, Outfptr, TopGrid, MetaData);

  // 8) KH Instability

  if (ProblemType == 8)
    ret = KHInitialize(fptr, Outfptr, TopGrid, MetaData);

  // 9) 2D/3D Noh Problem

  if (ProblemType == 9)
    ret = NohInitialize(fptr, Outfptr, TopGrid, MetaData);

  // 10) RotatingCylinder
 
  if (ProblemType == 10)
    ret = RotatingCylinderInitialize(fptr, Outfptr, TopGrid, MetaData);

  // 11) RadiatingShock
 
  if (ProblemType == 11)
    ret = RadiatingShockInitialize(fptr, Outfptr, TopGrid, MetaData);

  // 12) Free expansion blast wave
  if (ProblemType == 12)
    ret = FreeExpansionInitialize(fptr, Outfptr, TopGrid, MetaData);
 
  // 20) Zeldovich Pancake
 
  if (ProblemType == 20)
    ret = ZeldovichPancakeInitialize(fptr, Outfptr, TopGrid);
 
  // 21) 1D Pressureless collapse
 
  if (ProblemType == 21)
    ret = PressurelessCollapseInitialize(fptr, Outfptr, TopGrid, MetaData);
 
  // 22) Adiabatic expansion
 
  if (ProblemType == 22)
    ret = AdiabaticExpansionInitialize(fptr, Outfptr, TopGrid);
 
  // 23) GravityTest
 
  if (ProblemType == 23)
    ret = TestGravityInitialize(fptr, Outfptr, TopGrid, MetaData);
 
  // 24) Spherical Infall
 
  if (ProblemType == 24)
    ret = SphericalInfallInitialize(fptr, Outfptr, TopGrid, MetaData);
 
  // 25) TestGravitySphere
 
  if (ProblemType == 25)
    ret = TestGravitySphereInitialize(fptr, Outfptr, TopGrid, MetaData);
 
  // 26) GravityEquilibriumTest
 
  if (ProblemType == 26)
    ret = GravityEquilibriumTestInitialize(fptr, Outfptr, TopGrid, MetaData);
 
  // 27) CollapseTest
 
  if (ProblemType == 27)
    ret = CollapseTestInitialize(fptr, Outfptr, TopGrid, MetaData);
 
  // 28) TestGravityMotion
 
  if (ProblemType == 28)
    ret = TestGravityMotion(fptr, Outfptr, TopGrid, MetaData);

  // 29) TestOrbit
  if (ProblemType == 29)
    ret = TestOrbitInitialize(fptr, Outfptr, TopGrid, MetaData);
 
  // 30) Cosmology Simulation
 
  if (ProblemType == 30) {
    if (PartitionNestedGrids) {
      ret = NestedCosmologySimulationInitialize(fptr, Outfptr, TopGrid, MetaData);
    } else {
      ret = CosmologySimulationInitialize(fptr, Outfptr, TopGrid, MetaData);
    }
  }

  
 
  // 31) GalaxySimulation
  if (ProblemType == 31)
    ret = GalaxySimulationInitialize(fptr, Outfptr, TopGrid, MetaData);


// 35) Shearing Box Simulation
  if (ProblemType == 35) 
      ret = ShearingBoxInitialize(fptr, Outfptr, TopGrid, MetaData);
  if (ProblemType == 36) 
      ret = ShearingBox2DInitialize(fptr, Outfptr, TopGrid, MetaData);
  if (ProblemType == 37) 
      ret = ShearingBoxStratifiedInitialize(fptr, Outfptr, TopGrid, MetaData);
 
   
  // 40) Supernova Explosion from restart
 
  if (ProblemType == 40)
    ret = SupernovaRestartInitialize(fptr, Outfptr, TopGrid, MetaData,
				     Exterior);
 
  // 50) Photon Test

#ifdef TRANSFER
  if (ProblemType == 50)
    ret = PhotonTestInitialize(fptr, Outfptr, TopGrid, MetaData);
#endif /* TRANSFER */

  // 60) Turbulence Simulation.
 
  if (ProblemType == 60)
    ret = TurbulenceSimulationInitialize(fptr, Outfptr, TopGrid, MetaData);

  // 61) Protostellar Collapse
  if (ProblemType == 61)
    ret = ProtostellarCollapseInitialize(fptr, Outfptr, TopGrid, MetaData);

  // 62) Cooling test problem
  if (ProblemType == 62)
    ret = CoolingTestInitialize(fptr, Outfptr, TopGrid, MetaData);

  // Insert new problem intializer here...

    if (ProblemType ==300) {
    ret = PoissonSolverTestInitialize(fptr, Outfptr, TopGrid, MetaData);
  }




  /* 100) 1D HD Test */
  if (ProblemType == 100) {
    ret = Hydro1DTestInitialize(fptr, Outfptr, TopGrid, MetaData);
  }

  /* 101) 3D Collapse */
  if (ProblemType == 101) {
    ret = Collapse3DInitialize(fptr, Outfptr, TopGrid, MetaData);
  }

  /* 102) 1D Spherical Collapse */
  if (ProblemType == 102) {
    ret = Collapse1DInitialize(fptr, Outfptr, TopGrid, MetaData);
  }

  /* 106) Hydro and MHD Turbulence problems/Star Formation */
  if (ProblemType == 106) {
    ret = TurbulenceInitialize(fptr, Outfptr, TopGrid, MetaData, 0);
  }

  // 107) Put Sink from restart
 
  if (ProblemType == 107)
    ret = PutSinkRestartInitialize(fptr, Outfptr, TopGrid, MetaData,
				     Exterior);
 

  /* 200) 1D MHD Test */
  if (ProblemType == 200) {
    ret = MHD1DTestInitialize(fptr, Outfptr, TopGrid, MetaData);
  }

  /* 201) 2D MHD Test */
  if (ProblemType == 201) {
    ret = MHD2DTestInitialize(fptr, Outfptr, TopGrid, MetaData);
  }

  /* 202) 3D MHD Collapse */
  if (ProblemType == 202) {
    ret = CollapseMHD3DInitialize(fptr, Outfptr, TopGrid, MetaData);
  }

  /* 203) MHD Turbulence Collapse */
  if (ProblemType == 203) {
    ret = MHDTurbulenceInitialize(fptr, Outfptr, TopGrid, MetaData, 0);
  }

  /* 204) 3D MHD Test */
  if (ProblemType == 204) {
    ret = MHD3DTestInitialize(fptr, Outfptr, TopGrid, MetaData);
  }

  /* 207) Galaxy Disk */
  if (ProblemType == 207) {
    ret = GalaxyDiskInitialize(fptr, Outfptr, TopGrid, MetaData);
  }

  /* 208) AGN Disk */
  if (ProblemType == 208) {
    ret = AGNDiskInitialize(fptr, Outfptr, TopGrid, MetaData);
  }


 
  
  if (ret == INT_UNDEFINED) {
    fprintf(stderr, "Problem Type %"ISYM" undefined.\n", ProblemType);
    ENZO_FAIL("");
  }
 
  if (ret == FAIL) {
        ENZO_FAIL("Error in problem initialization.");
  }
 
  if (debug)
    printf("InitializeNew: Finished problem initialization.\n");

  /* If requested, initialize streaming data files. */

  InitializeMovieFile(MetaData, TopGrid);
 
  /* Do some error checking */
 
  if (MetaData.StopTime == FLOAT_UNDEFINED)
    ENZO_FAIL("StopTime never set.");

  int nFields = TopGrid.GridData->ReturnNumberOfBaryonFields();
  if (nFields >= MAX_NUMBER_OF_BARYON_FIELDS) {
    printf("NumberOfBaryonFields (%"ISYM") exceeds "
	   "MAX_NUMBER_OF_BARYON_FIELDS (%"ISYM").\n", 
	   nFields, MAX_NUMBER_OF_BARYON_FIELDS);
    ENZO_FAIL("");
  }

  PrintMemoryUsage("1st Initialization done");

 
  if (debug)
    printf("Initialize Exterior\n");
 
  // Initialize the exterior (unless it was set in the problem initializer)
 
  if (Exterior.AmIPrepared() == FALSE) {

    Exterior.Prepare(TopGrid.GridData);   // set rank and dims

    if (MetaData.BoundaryConditionName != NULL) {

      if ((BCfptr = fopen(MetaData.BoundaryConditionName, "r")) == NULL) {
	fprintf(stderr, "Error opening BC file: %s\n",
		MetaData.BoundaryConditionName);
	ENZO_FAIL("");
      }

      fprintf(stderr, "Opened BC file mode r\n");

      if (Exterior.ReadExternalBoundary(BCfptr) == FAIL) {
		ENZO_FAIL("Error in ReadExternalBoundary.");
      }
      fclose(BCfptr);
    } else 
      {
      if (debug) 
        fprintf(stderr, "InitializeExternalBoundaryFace\n");

      SimpleConstantBoundary = TRUE;

      for (dim = 0; dim < MetaData.TopGridRank; dim++) {
        if (MetaData.LeftFaceBoundaryCondition[dim] != periodic ||
            MetaData.RightFaceBoundaryCondition[dim] != periodic) {
          SimpleConstantBoundary = FALSE;
        }
      }

      if (debug) {
        if (SimpleConstantBoundary) {
          fprintf(stderr, "SimpleConstantBoundary TRUE\n");
        } else {
          fprintf(stderr, "SimpleConstantBoundary FALSE\n");
        }
      }
        
      for (dim = 0; dim < MetaData.TopGridRank; dim++)
	if (Exterior.InitializeExternalBoundaryFace(dim,
				    MetaData.LeftFaceBoundaryCondition[dim],
				    MetaData.RightFaceBoundaryCondition[dim],
				    Dummy, Dummy)
	    == FAIL) {
	  	  ENZO_FAIL("Error in InitializeExternalBoundaryFace.");
	}
 
      // Initialize particle boundary conditions
 
      Exterior.InitializeExternalBoundaryParticles(
					  MetaData.ParticleBoundaryType);
 
    }  // end: if (MetaData.BoundaryConditionName != NULL)
 
  }  // end of set Exterior


  PrintMemoryUsage("Exterior set");

  if (debug) {
    fprintf(stderr, "End of set exterior\n");
  }
 
  // Set values that were left undefined (above)
 
  if (MetaData.TimeLastDataDump == FLOAT_UNDEFINED)
    MetaData.TimeLastDataDump = MetaData.Time - MetaData.dtDataDump*1.00001;
  if (MetaData.TimeLastHistoryDump == FLOAT_UNDEFINED)
    MetaData.TimeLastHistoryDump = MetaData.Time - MetaData.dtHistoryDump;
 
  if (MetaData.TimeLastTracerParticleDump == FLOAT_UNDEFINED)
    MetaData.TimeLastTracerParticleDump =
                           MetaData.Time - MetaData.dtTracerParticleDump;
 
  if (MetaData.CycleLastDataDump == INT_UNDEFINED)
    MetaData.CycleLastDataDump = MetaData.CycleNumber -
                                 MetaData.CycleSkipDataDump;
  if (MetaData.CycleLastHistoryDump == INT_UNDEFINED)
    MetaData.CycleLastHistoryDump = MetaData.CycleNumber -
                                    MetaData.CycleSkipHistoryDump;
 
  // Make changes required for Zeus solver, and turn the TotalEnergy
  // variable (should be renamed just Energy) into GasEnergy
 
  if (HydroMethod == Zeus_Hydro &&
      ProblemType != 10 &&  // BWO (Rotating cylinder)
      ProblemType != 11 &&  // BWO (radiating shock)
      ProblemType != 20 &&
      ProblemType != 27 &&
      ProblemType != 30 &&
      ProblemType != 31 &&  // BWO (isolated galaxies)
      ProblemType != 60) //AK
    ConvertTotalEnergyToGasEnergy(&TopGrid);
 
  // If using StarParticles, set the number to zero 
  // (assuming it hasn't already been set)
  if (NumberOfStarParticles == NULL)
    if (StarParticleCreation || StarParticleFeedback)
      NumberOfStarParticles = 0;
 
  // Convert minimum initial overdensity for refinement to mass
  // (unless MinimumMass itself was actually set)


  for (i = 0; i < MAX_FLAGGING_METHODS; i++)
    if (MinimumMassForRefinement[i] == FLOAT_UNDEFINED) {
      MinimumMassForRefinement[i] = MinimumOverDensityForRefinement[i];
      for (dim = 0; dim < MetaData.TopGridRank; dim++)
	MinimumMassForRefinement[i] *=
	  (DomainRightEdge[dim]-DomainLeftEdge[dim])/
	  float(MetaData.TopGridDims[dim]);
    }
 
  // Check for the creation of tracer particles
  // Tracer particles will not be created at this point if ||rgio in ON
 
  if (TracerParticleCreation(fptr, TopGrid, MetaData) == FAIL) {
        ENZO_FAIL("Error in TracerParticleCreation");
  }
 
  // Write the MetaData/global values to the Parameter file
 
  if (MyProcessorNumber == ROOT_PROCESSOR)
    if (WriteParameterFile(Outfptr, MetaData) == FAIL) {
            ENZO_FAIL("Error in WriteParameterFile.");
    }
 
  if (debug)
    printf("InitializeNew: Initial grid hierarchy set\n");
 
  // Walk the grids
 
  HierarchyEntry *CurrentGrid;
  FLOAT WT = -1.0;
  int GP = 1;
  int gridcounter = 0;
 
  CurrentGrid = &TopGrid;


 
  while (CurrentGrid != NULL) {
 
    // WriteHierarchyStuff(stderr, CurrentGrid, "UUUU", GP, WT);
 
    if (debug)
      printf("InitializeNew: Partition Initial Grid %"ISYM"\n", gridcounter);
 
    CommunicationPartitionGrid(CurrentGrid, gridcounter);
 
    gridcounter++;

    if (PartitionNestedGrids)
      CurrentGrid = CurrentGrid->NextGridNextLevel;
    else
      CurrentGrid = NULL;
 
  }
 
  // For problem 30, using ParallelGridIO,
  // read in data only after partitioning the grid
 
  if (debug)
    if (ParallelRootGridIO == TRUE && ProblemType == 30) {
      if (PartitionNestedGrids) {
        printf("InitializeNew: Re-initialize NestedCosmologySimulation\n");
      } else {
        printf("InitializeNew: Re-initialize CosmologySimulation\n");
      }
    }

  PrintMemoryUsage("Before 2nd pass");
 
  if (ParallelRootGridIO == TRUE && ProblemType == 30) {
    if (PartitionNestedGrids) {
      if (NestedCosmologySimulationReInitialize(&TopGrid, MetaData) == FAIL) {
                ENZO_FAIL("Error in NestedCosmologySimulationReInitialize.");
      }
    } else {
      if (CosmologySimulationReInitialize(&TopGrid, MetaData) == FAIL) {
                ENZO_FAIL("Error in CosmologySimulationReInitialize.");
      }
    }
  }

  PrintMemoryUsage("After 2nd pass");

  // For problem 60, using ParallelGridIO, read in data only after
  // partitioning grid.
 
  if (ParallelRootGridIO == TRUE && ProblemType == 60)
    if (TurbulenceSimulationReInitialize(&TopGrid, MetaData) == FAIL) {
            ENZO_FAIL("Error in TurbulenceSimulationReInitialize.");
    }
 
 if (ProblemType == 106)
   if (TurbulenceInitialize(fptr, Outfptr, TopGrid, MetaData, 1)
       == FAIL) {
     fprintf(stderr, "Error in TurbulenceReInitialize.\n");
     ENZO_FAIL("");
   }

  // For ProblemType 203 (Turbulence Simulation we only initialize the data
  // once the topgrid has been split.
 if (ProblemType == 203)
   if (MHDTurbulenceInitialize(fptr, Outfptr, TopGrid, MetaData, 1)
       == FAIL) {
     fprintf(stderr, "Error in MHDTurbulenceReInitialize.\n");
     ENZO_FAIL("");
   }
   
 
  // Close parameter files
 
  fclose(fptr);
 
  if (MyProcessorNumber == ROOT_PROCESSOR)
    fclose(Outfptr);

  PrintMemoryUsage("Exit X_Init");

  // 2006-12-11 Skory bug fix for star particle miscounts
  // Added the following line:

  CommunicationBroadcastValue(&MetaData.NumberOfParticles, ROOT_PROCESSOR);
 
  MetaData.FirstTimestepAfterRestart = FALSE;
  


  return SUCCESS;
 
}
 
 
 
 
void ConvertTotalEnergyToGasEnergy(HierarchyEntry *Grid)
{
  if (Grid != NULL) {
    Grid->GridData->ConvertTotalEnergyToGasEnergy();
    ConvertTotalEnergyToGasEnergy(Grid->NextGridThisLevel);
    ConvertTotalEnergyToGasEnergy(Grid->NextGridNextLevel);
  }
}
