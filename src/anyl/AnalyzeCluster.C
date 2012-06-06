/***********************************************************************
/
/  COMPUTE A CLUSTER PROFILE
/
/  written by: Greg Bryan
/  date:       August, 1997
/
/  modified1: Ji-hoon Kim
/             July, 2009 
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#ifdef USE_HDF4
#include <df.h>
#endif /* USE_HDF4 */
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
#define DEFINE_STORAGE
#include "../enzo/EnzoTiming.h"
#include "../enzo/ErrorExceptions.h"
#include "../enzo/macros_and_parameters.h"
#include "../enzo/units.h"
#include "../enzo/typedefs.h"
#include "../enzo/global_data.h"
#include "../enzo/Fluxes.h"
#include "../enzo/GridList.h"
#include "../enzo/ExternalBoundary.h"
#include "../enzo/Grid.h"
#include "../enzo/Hierarchy.h"
#include "../enzo/LevelHierarchy.h"
#include "../enzo/TopGridData.h"
#include "../enzo/CosmologyParameters.h"  
#include "../enzo/communication.h"
#include "../enzo/flowdefs.h"
#include "../enzo/PhotonCommunication.h"

#undef DEFINE_STORAGE

#define MAX_BINS 200

/* function prototypes */
int ReadAllData(char *filename, HierarchyEntry *TopGrid, TopGridData &tgd,
		ExternalBoundary *Exterior, float *Inititaldt);
int Group_ReadAllData(char *filename, HierarchyEntry *TopGrid, TopGridData &tgd,
              ExternalBoundary *Exterior, float *Initialdt,
              bool ReadParticlesOnly=false);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int SetDefaultGlobalValues(TopGridData &MetaData);
int CopyOverlappingZones(grid* CurrentGrid, TopGridData *MetaData, 
			 LevelHierarchyEntry *LevelArray[], int level);
int AnalyzeClusterReadParameterFile(char *filename, int &NumberOfCenters,
				    FLOAT *CenterList[],
				    AnalyzeClusterParameters *parameters);
int AnalyzeClusterComputeClumpingFactor(LevelHierarchyEntry *LevelArray[],
					TopGridData *MetaData,
					int NumberOfGridPoints,
					int NumberOfParticles,
					FLOAT SphereCenter[], 
					float SphereRadius, char *Name);
int GenerateGridArray(LevelHierarchyEntry *LevelArray[], int level,
		      HierarchyEntry **Grids[]);
int SetBoundaryConditions(HierarchyEntry *Grids[], int NumberOfGrids,
			  SiblingGridList SiblingList[],
			  int level, TopGridData *MetaData, 
			  ExternalBoundary *Exterior, LevelHierarchyEntry *Level);
int CommunicationInitialize(Eint32 *argc, char **argv[]);
int CommunicationFinalize();
int CommunicationSumValues(FLOAT *values, int number);
int CommunicationAllSumValues(FLOAT *values, int number);
void FindVirialRadius(LevelHierarchyEntry *LevelArray[], float &rvir, float &r500,
		      float critical_density, float BoxSize,
		      FLOAT *Center, FLOAT OuterEdge, float MeanVelocity[MAX_DIMENSION][3], 
		      int NumberOfPoints, float *ProfileRadius, float ProfileValue[][MAX_PROFILES],
		      float ProfileWeight[][MAX_PROFILES], char *ProfileName[MAX_PROFILES],
		      AnalyzeClusterParameters *parameters);
void my_exit(int status);
int FastSiblingLocatorInitialize(ChainingMeshStructure *Mesh, int Rank,
				 int TopGridDims[]);
int FastSiblingLocatorFinalize(ChainingMeshStructure *Mesh);
int CosmologyGetUnits(float *DensityUnits, float *LengthUnits, 
              float *TemperatureUnits, float *TimeUnits,
              float *VelocityUnits, FLOAT Time);
int GetUnits(float *DensityUnits, float *LengthUnits, 
              float *TemperatureUnits, float *TimeUnits,
              float *VelocityUnits, FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

Eint32 main(Eint32 argc, char *argv[])
{
  CommunicationInitialize(&argc, &argv);

  /* Main declarations */

  TopGridData MetaData;
  HierarchyEntry TopGrid;
  ExternalBoundary Exterior;
  LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];

  /* Initialize */

  debug                = (MyProcessorNumber == ROOT_PROCESSOR) ? TRUE : FALSE;
  char *myname         = argv[0];
  AnalyzeClusterParameters parameters;

  /* general declarations. */

  int i, j, profile, dim, ret, level, NumberOfGrids, grid, grid1;
  float pi = 3.14159;
  const double Mpc = 3.0856e24, SolarMass = 1.989e33, GravConst = 6.67e-8;
  char Name[MAX_LINE_LENGTH], DiskImageName[MAX_LINE_LENGTH];
  LevelHierarchyEntry *Temp2, *Temp;
  HierarchyEntry **Grids;
  float DiskVector[3];

  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    LevelArray[level] = NULL;

  /* Error check */

  if (argc != 3) {
    fprintf(stderr, "usage: %s amr_file anyl_parameter_file\n", myname);
    my_exit(EXIT_FAILURE);
  }

  /* Read the saved file. */

  SetDefaultGlobalValues(MetaData); 

  // First expect to read in packed-HDF5
  float dummy;
#ifdef USE_HDF5_GROUPS
  if (Group_ReadAllData(argv[1], &TopGrid, MetaData, &Exterior, &dummy) == FAIL) {
      if (MyProcessorNumber == ROOT_PROCESSOR) {
	fprintf(stderr, "Error in Group_ReadAllData %s\n", argv[1]);
	fprintf(stderr, "Probably not in a packed-HDF5 format. Trying other read routines.\n");
      }
#endif
      // If not packed-HDF5, then try usual HDF5 or HDF4
      if (ReadAllData(argv[1], &TopGrid, MetaData, &Exterior, &dummy) == FAIL) {
	if (MyProcessorNumber == ROOT_PROCESSOR) {
	  fprintf(stderr, "Error in ReadAllData %s.\n", argv[1]);
	}
	my_exit(EXIT_FAILURE);
      }
#ifdef USE_HDF5_GROUPS
    }
#endif

  // recursively add levels

  AddLevel(LevelArray, &TopGrid, 0);    

  /* Read the anyl parameter file. */

  int NumberOfCenters;
  FLOAT *CenterList[MAX_DIMENSION]; 

  AnalyzeClusterReadParameterFile(argv[2], NumberOfCenters, CenterList,  
				  &parameters);

  int NumberOfPoints = parameters.npoints;

  /* Set default constants for Cosmology Simulation(Type=30) if ComovingCoordinates = 0 */

  if (ComovingCoordinates != 1) {
    InitialRedshift = 0; 
    //    FinalRedshift = 0;
    HubbleConstantNow = 0.7; 
    OmegaMatterNow = 0.3;
    OmegaLambdaNow = 0.7;
    //    float ComovingBoxSize = 1;
    //    float MaxExpansionRate = 1;
  }  

  /* From the time, compute the current redshift. */

  FLOAT a=1, dadt, CurrentRedshift = 0.0;
  float OmegaCurvatureNow = 0.0, Esquared, average_dens, critical_density = 1.0;
  if (ComovingCoordinates) {
    if (CosmologyComputeExpansionFactor(LevelArray[0]->GridData->ReturnTime(),
					&a, &dadt) == FAIL) {
      fprintf(stderr, "Error in ComputeExpansionFactor.\n");
      my_exit(EXIT_FAILURE);
    }

    CurrentRedshift = (1.0+InitialRedshift)/a - 1.0;

    /* If unset, then set virial_dens according to spherical top-hot collapse
       (see Bryan & Norman 1998). */

    OmegaCurvatureNow = 1.0 - OmegaMatterNow - OmegaLambdaNow;
    Esquared = OmegaMatterNow    * pow(1.0+CurrentRedshift, 3) +
               OmegaCurvatureNow * pow(1.0+CurrentRedshift, 2) +
               OmegaLambdaNow;
    float x = OmegaMatterNow * pow(1.0+CurrentRedshift, 3) / Esquared - 1.0;
    if (parameters.virial_dens < 0) {
      if (OmegaCurvatureNow < 1.0e-4)
	parameters.virial_dens = 18.0*pi*pi + 82.0*x - 39.0*x*x;
      else if (OmegaLambdaNow < 1.0e-4)
	parameters.virial_dens = 18.0*pi*pi + 60.0*x - 32.0*x*x;
      else {
	printf("This cosmology not supported.\n");
	my_exit(EXIT_FAILURE);
      }
    }

    average_dens = 2.78e11*OmegaMatterNow*pow(HubbleConstantNow, 2) *
      pow((1+InitialRedshift)/a, 3);
    critical_density = 2.78e11*pow(HubbleConstantNow, 2) * Esquared;

  } // end: if (ComovingCoordinates)

  else { 

    CurrentRedshift = (1.0+InitialRedshift)/a - 1.0;

    OmegaCurvatureNow = 1.0 - OmegaMatterNow - OmegaLambdaNow;
    Esquared = OmegaMatterNow    * pow(1.0+CurrentRedshift, 3) +
               OmegaCurvatureNow * pow(1.0+CurrentRedshift, 2) +
               OmegaLambdaNow;
  }

  /* Set BoxSize. */

  float BoxSize = 1; 
  float DensityConversion = 1, VelocityConversion = 1;  

  float DensityUnits = 1, LengthUnits = 1, VelocityUnits = 1, TimeUnits = 1,
        TemperatureUnits = 1; 

  if (ComovingCoordinates)  {
    BoxSize = ComovingBoxSize/HubbleConstantNow;  
    // Added by Matt Turk, to get TimeUnits.
    if (CosmologyGetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
			  &TimeUnits, &VelocityUnits, 
			  LevelArray[0]->GridData->ReturnTime()) == FAIL) {
      fprintf(stderr, "Error in CosmologyGetUnits.\n");
      return FAIL;
    }
  }
  else {
    
    if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,  
              &TimeUnits, &VelocityUnits, LevelArray[0]->GridData->ReturnTime()) == FAIL) {
      fprintf(stderr, "Error in GetUnits.\n");
      return FAIL;
    }

    BoxSize = LengthUnits/Mpc;  

  }

  FLOAT InnerEdge = parameters.rinner/BoxSize, OuterEdge = parameters.router/BoxSize;  

  BoxSize *= a/(1+InitialRedshift); 

  /* ------------------------------------------------------------ */
  /* Set the boundary conditions.                                 */

  if (debug) printf("->setting BCs\n");
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
    if (debug) printf("level %d\n", level);
    NumberOfGrids = GenerateGridArray(LevelArray, level, &Grids);
    if (level > 0)
      for (grid = 0; grid < NumberOfGrids; grid++)
	Grids[grid]->ParentGrid->GridData->SetTime(
				       Grids[grid]->GridData->ReturnTime());

    /*----------------------------------------------*/
    /* Create chaining mesh for FastSiblingLocator. */
    /*----------------------------------------------*/
    
    /* Initialize the chaining mesh used in the FastSiblingLocator. */

    ChainingMeshStructure ChainingMesh;
    FastSiblingLocatorInitialize(&ChainingMesh, MetaData.TopGridRank,
				 MetaData.TopGridDims);
    SiblingGridList *SiblingList = new SiblingGridList[NumberOfGrids];
    if (SiblingList == NULL) {
      fprintf(stderr, "Error allocating SiblingList\n");
      return FAIL;
    }

    /* Add all the grids to the chaining mesh. */

    for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
      Grids[grid1]->GridData->FastSiblingLocatorAddGrid(&ChainingMesh);

    /* For each grid, get a list of possible siblings from the chaining mesh. */

    for (grid1 = 0; grid1 < NumberOfGrids; grid1++)
      if (Grids[grid1]->GridData->
	  FastSiblingLocatorFindSiblings(&ChainingMesh, &SiblingList[grid1],
				     MetaData.LeftFaceBoundaryCondition,
				     MetaData.RightFaceBoundaryCondition) == FAIL) {
	fprintf(stderr, "Error in grid->FastSiblingLocatorFindSiblings.\n");
	return FAIL;
      }

    /* Clean up the chaining mesh. */

    FastSiblingLocatorFinalize(&ChainingMesh);	

    SetBoundaryConditions(Grids, NumberOfGrids, SiblingList, level, &MetaData, 
			  &Exterior, LevelArray[level]);

  } // ENDFOR level

  /* ------------------------------------------------------------ */
  /* Zero the solution (on this grid) which is underneath any subgrid
     (so we get only the high resolution solution from the subgrid). */

  if (debug) printf("->zeroing redundant solution\n");
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY-1; level++) {
    Temp = LevelArray[level];
    while (Temp != NULL) {
      Temp->GridData->ZeroSolutionUnderSubgrid(NULL, ZERO_UNDER_SUBGRID_FIELD);
      Temp2 = LevelArray[level+1];
      while (Temp2 != NULL) {
	Temp->GridData->ZeroSolutionUnderSubgrid(Temp2->GridData, 
						 ZERO_UNDER_SUBGRID_FIELD);
	Temp2 = Temp2->NextGridThisLevel;
      }
      Temp = Temp->NextGridThisLevel;
    }
  }

  /* ------------------------------------------------------------ */
  /* Loop over centers. */

  FLOAT Center[MAX_DIMENSION];
  FLOAT CTCenter[MAX_DIMENSION];
  for (int center = 0; center < NumberOfCenters; center++) {
    if (NumberOfCenters > 1) {
      debug = FALSE;
      printf("Computing center %d\n", center);
    }
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      Center[dim] = CenterList[dim][center];
      //fprintf(stderr, "Center[dim] = %g, CenterList[dim][0] = %g", Center[dim], CenterList[dim][center]);
    }

    /* Set base name. */

    char AnalyzeBaseName[MAX_LINE_LENGTH];
    strcpy(AnalyzeBaseName, argv[1]);
    strcat(AnalyzeBaseName, "_analyze");

    if (NumberOfCenters == 1)
      sprintf(Name, "%s", AnalyzeBaseName);
    else
      sprintf(Name, "%s%.3d", AnalyzeBaseName, center);

    /* ------------------------------------------------------------ */
    /* Find the highest density spot if Center not specified. */

    float MaxDensity = -huge_number;
    //float MinCoolingTime = huge_number;
    int gridID = 0;
    if (Center[0] == FLOAT_UNDEFINED) {
      if (debug) printf("->finding maximum\n");
      for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
    //float MinCoolingTime = huge_number;
	Temp = LevelArray[level];
     printf("LEVEL = %i\n", level+1); 
	while (Temp != NULL) {
      //printf("\tlevel = %i grid = %i\n", level+1, gridID++);
	  Temp->GridData->FindMaximumBaryonDensity(Center, &MaxDensity);
      //Temp->GridData->FindMinimumCoolingTime(CTCenter, &MinCoolingTime);
	  Temp = Temp->NextGridThisLevel;
	}
     //printf("-> -> MinCoolingTime = %g (TU = %0.4e)\n", MinCoolingTime,TimeUnits);
     //printf("-> -> Center = %g %g %g\n", CTCenter[0],CTCenter[1],CTCenter[2]);
      }

#ifdef USE_MPI

      if (NumberOfProcessors > 1) {

	/* Copy center and maxdensity into a buffer and send to root. */

	FLOAT Buffer1[MAX_DIMENSION+1], *Buffer2;
	Buffer2 = new FLOAT[NumberOfProcessors*(MAX_DIMENSION+1)];
	for (dim = 0; dim < MAX_DIMENSION; dim++)
	  Buffer1[dim+1] = Center[dim];
	Buffer1[0] = FLOAT(MaxDensity);
	MPI_Gather(Buffer1, MAX_DIMENSION+1, MY_MPIFLOAT, 
		   Buffer2, NumberOfProcessors*(MAX_DIMENSION+1), MY_MPIFLOAT,
		   ROOT_PROCESSOR, MPI_COMM_WORLD);

	/* Find max on root. */

	int MaxIndex = 0;
	if (MyProcessorNumber == ROOT_PROCESSOR)
	  for (i = 0; i < NumberOfProcessors; i++)
	    if (Buffer2[i*(MAX_DIMENSION+1)] > MaxDensity) {
	      MaxIndex = i*(MAX_DIMENSION+1);
	      MaxDensity = Buffer2[i*(MAX_DIMENSION+1)];
	    }

	/* Broadcast back center and maxdensity. */

	MPI_Bcast(Buffer2+MaxIndex, MAX_DIMENSION+1, MY_MPIFLOAT, 
		  ROOT_PROCESSOR, MPI_COMM_WORLD);
	MaxDensity = float(Buffer2[MaxIndex]);
	for (dim = 0; dim < MAX_DIMENSION; dim++)
	  Center[dim] = Buffer2[MaxIndex+dim+1];
	delete [] Buffer2;

      }

#endif /* USE_MPI */      

    } // end: if (Center[0] == FLOAT_UNDEFINED)

    /* ------------------------------------------------------------ */
    /* Set up radial grid. */

    FLOAT ProfileValue[MAX_BINS][MAX_PROFILES], ProfileRadius[MAX_BINS+1],
          ProfileWeight[MAX_BINS][MAX_PROFILES];
    char  *ProfileName[MAX_PROFILES];

    for (i = 0; i < MAX_PROFILES; i++) {
      for (j = 0; j < MAX_BINS; j++) {
	ProfileValue[j][i] = 0.0;
	ProfileWeight[j][i] = 0.0;
      }
      ProfileName[i] = NULL;
    }

    /* Compute radii (ProfileRadius is inner edge of bin). */ 

    ProfileRadius[0] = 0;
    float dlogRadius = (log10(OuterEdge) - log10(InnerEdge)) /
      float(NumberOfPoints-1);
    for (i = 1; i < NumberOfPoints+1; i++)
      ProfileRadius[i] = POW(10, log10(InnerEdge) + float(i-1)*dlogRadius); 

    /* Linearly binned "ProfileRadius2" for .Vertical profile only 
       becuase logarithmic interval may not make sense for vertical profile. - Ji-hoon Kim */

    FLOAT ProfileRadius2[MAX_BINS+1];  
    ProfileRadius2[0] = 0;
    float dRadius = (OuterEdge - InnerEdge) /
      float(NumberOfPoints-1);
    for (i = 1; i < NumberOfPoints+1; i++)
      ProfileRadius2[i] = InnerEdge + float(i-1)*dRadius;  

    /* ------------------------------------------------------------ */  
    /* Calculate rvirial. */
    /* Compute total density and find rvir (dvir is overdensity of virial_dens 
       M(solar)/Mpc^3). */

    /* First, compute profiles. */

    if (debug) printf("->computing r_virial (now with critical_dens)\n");
    FLOAT MeanVelocity[MAX_DIMENSION][3], MeanVelocityWeight[MAX_DIMENSION][3];
    for (i = 0; i < MAX_DIMENSION; i++) {
      MeanVelocity[i][0] = 0;
      MeanVelocity[i][1] = 0;
    }

    float rvir = 0, mvir = 0, mvir_gas = 0, spin_gas = 0, spin_dm = 0, mvir_star = 0, r500 = 0;

    FindVirialRadius(LevelArray, rvir, r500,
                     critical_density, BoxSize, Center, OuterEdge, MeanVelocity,
                     NumberOfPoints, ProfileRadius, ProfileValue, ProfileWeight, ProfileName,
                     &parameters);

    if (rvir > OuterEdge || rvir == 0) {
      printf("warning: rvir (%0.9e) > OuterEdge (%0.9e)\n", rvir, OuterEdge); 
      rvir = OuterEdge;
    }

    FLOAT NewCenter[MAX_DIMENSION], NewCenterWeight;
    float rnew = rvir;

// Loop 100 times to find the centroid of the halo.
#define MAX_ITERATIONS 100
    for (i = 0; i < MAX_ITERATIONS; i++) {
      rnew = 0.90*rnew;
      printf("NewCenter(%d) = %f %f %f  r = %g\n", i, Center[0], Center[1], Center[2],rnew);
      for (dim = 0; dim < MAX_DIMENSION; dim++)
	    NewCenter[dim] = 0;
      NewCenterWeight = 0;
      for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
	    LevelHierarchyEntry *Temp = LevelArray[level];
	    while (Temp != NULL) {
	      Temp->GridData->FindMeanVelocityAndCenter(Center, rnew, NewCenter, NewCenterWeight,
						    MeanVelocity, MeanVelocityWeight);
	      Temp = Temp->NextGridThisLevel;
	    }
      }
      CommunicationAllSumValues(NewCenter, MAX_DIMENSION);
      CommunicationAllSumValues(&NewCenterWeight, 1);

      for (dim = 0; dim < MAX_DIMENSION; dim++)
	if (NewCenterWeight > 0) {
	  NewCenter[dim] /= NewCenterWeight;
	  Center[dim] += NewCenter[dim];
	}
//      printf("shift = %g %g %g  weight = %g\n", NewCenter[0], NewCenter[1], NewCenter[2], NewCenterWeight);
    }

    FindVirialRadius(LevelArray, rvir, r500, 
		     critical_density, BoxSize, Center, OuterEdge, MeanVelocity,
		     NumberOfPoints, ProfileRadius, ProfileValue, ProfileWeight, ProfileName,
		     &parameters);
    
    mvir = pow(rvir*BoxSize, 3)*4.0/3.0*pi*critical_density*
           parameters.virial_dens;

    /* ------------------------------------------------------------ */
    /* Calculate the mean velocity. */

    for (i = 0; i < MAX_DIMENSION; i++) {
      MeanVelocity[i][0] = 0;
      MeanVelocity[i][1] = 0;
      MeanVelocityWeight[i][0] = 0;
      MeanVelocityWeight[i][1] = 0;
    }
    FLOAT MeanVelocityOuterEdge = parameters.MeanVelocityVirialFraction*rvir;

    if (debug) printf("->finding mean velocity\n");
    for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
      LevelHierarchyEntry *Temp = LevelArray[level];
      while (Temp != NULL) {
	Temp->GridData->FindMeanVelocityAndCenter(Center, MeanVelocityOuterEdge,
						  NewCenter, NewCenterWeight,
						  MeanVelocity, MeanVelocityWeight);
	Temp = Temp->NextGridThisLevel;
      }
    }
    CommunicationAllSumValues(&MeanVelocity[0][0], MAX_DIMENSION*2);
    CommunicationAllSumValues(&MeanVelocityWeight[0][0], MAX_DIMENSION*2);

    /* Compute mean velocities within MeanVelocityOuterEdge for gas [0], 
       dm [1], and total [2]. */

    for (i = 0; i < MAX_DIMENSION; i++) {
      if ((MeanVelocityWeight[i][0] + MeanVelocityWeight[i][1]) > 0)
	MeanVelocity[i][2] = (MeanVelocity[i][0] + MeanVelocity[i][1])/
	  (MeanVelocityWeight[i][0] + MeanVelocityWeight[i][1]);
      if (MeanVelocityWeight[i][0] > 0)
	MeanVelocity[i][0] /= MeanVelocityWeight[i][0];
      if (MeanVelocityWeight[i][1] > 0)
	MeanVelocity[i][1] /= MeanVelocityWeight[i][1];
    }

    /* Set Tvirial (in K) and (if required) the cold temperature cutoff
       (from Bryan & Norman, 1998). */

    float mu = 1.22;    // Set the mean mass per particle either to fully
    if (mvir > 1.0e13)  //   ionized or fully neutral depending on (arbitrary)
      mu = 0.59;        //    cutoff.
    float Tvir = parameters.VirialTemperatureNormalization * 2.73e7 * mu *
                 pow(mvir/1.0e15, 2.0/3.0) * 
                 pow(HubbleConstantNow*HubbleConstantNow*
		     parameters.virial_dens * Esquared, float(1.0/3.0));
    if (parameters.ColdTemperatureCutoffVirialFraction > 0)
      parameters.ColdTemperatureCutoff = 
	parameters.ColdTemperatureCutoffVirialFraction * Tvir;

    /* ------------------------------------------------------------ */
    /* Calculate profiles. */

    for (i = 0; i < MAX_PROFILES; i++) {
      for (j = 0; j < MAX_BINS; j++) {
	ProfileValue[j][i] = 0;
	ProfileWeight[j][i] = 0;
      }
      ProfileName[i] = NULL;
    }

    if (debug) printf("->finding profiles\n");
    for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
      LevelHierarchyEntry *Temp = LevelArray[level];
      while (Temp != NULL) {
	Temp->GridData->AddToRadialProfile(Center, OuterEdge, MeanVelocity, 
					   NumberOfPoints, ProfileRadius,
					   ProfileValue, ProfileWeight,
					   ProfileName, &parameters);
	Temp = Temp->NextGridThisLevel;
      }
    }
    CommunicationAllSumValues(&ProfileValue[0][0], MAX_BINS*MAX_PROFILES);
    CommunicationAllSumValues(&ProfileWeight[0][0], MAX_BINS*MAX_PROFILES);

    /* Compute the cummulative mass and overdensity. */    

    for (j = 0; j < NumberOfPoints; j++) {

      /* Compute the mass (gas+dm+star), in M(solar), within this annalus. */

      ProfileValue[j][38] = ProfileValue[j][0] + ProfileValue[j][30] +
	                    ProfileValue[j][60];

      /* Keep a running sum of the total mass within radius j+1. */  

      ProfileValue[j][29] = ProfileValue[max(j-1,0)][29] + ProfileValue[j][0];
      ProfileValue[j][57] = ProfileValue[max(j-1,0)][57] + ProfileValue[j][30];
      ProfileValue[j][68] = ProfileValue[max(j-1,0)][68] + ProfileValue[j][60];
      ProfileValue[j][90] = ProfileValue[max(j-1,0)][90] + ProfileValue[j][83]; 

      if (j > 0)
	ProfileValue[j][38] += ProfileValue[j-1][38];

      /* Compute the average overdensity within the sphere with radius [j+1].*/
      
      ProfileValue[j][39] = ProfileValue[j][38] / 
	(pow(ProfileRadius[j+1]*BoxSize, 3) * 4.0*pi/3.0) / critical_density;   

      /* Compute the dynamical time within the sphere with radius j+1. */

      ProfileValue[j][22] = sqrt(3*pi/(16*GravConst*
				       max(ProfileValue[j][39],1e-20)*
				       critical_density*SolarMass/(Mpc*Mpc*Mpc)));
    }

    ProfileName[22] = "T_dyn (s)";
    ProfileName[39] = "od_total";   
    ProfileName[29] = "m_gas (Ms)";   
    ProfileName[38] = "m_total (Ms)";  
    ProfileName[57] = "m_dm (Ms)";   
    if (ProfileName[60] != NULL)
      ProfileName[68] = "m_star (Ms)"; 
    ProfileName[90] = "m_gas_metal (Ms)";   
    
    /* If the ProfileWeight is < 0, it is to be replaced by the effective 
       volume in the annulus. */  

    for (i = 0; i < MAX_PROFILES; i++)
      for (j = 0; j < NumberOfPoints; j++)
	if (ProfileWeight[j][i] < 0)
	  ProfileWeight[j][i] = pow(BoxSize, 3) * 4.0*pi/3.0 *
	    (pow(ProfileRadius[j+1], 3) - pow(ProfileRadius[j], 3));

    /* Normalize profile values by dividing each by it's cumulative weight. */

    for (i = 0; i < MAX_PROFILES; i++)
      for (j = 0; j < NumberOfPoints; j++) {

	if (ProfileWeight[j][i] > 0.0){
	  ProfileValue[j][i] /= ProfileWeight[j][i];  
	  /*        fprintf(stderr,"profile is greater than zero! %i\n ",i);*/

	}

	/* If this is a radial velocity dispersion, subtract the mean velocity   
	   which MUST be the profile immediately before. */

	if (ProfileName[i] != NULL)
	  if (strstr(ProfileName[i], "vr_rms") != NULL)
	    ProfileValue[j][i] -= pow(ProfileValue[j][i-1],2);

	/* If this is a velocity rms profile, convert from 3D ms to 3D rms. */  
 
	if (ProfileName[i] != NULL)
	  if (strstr(ProfileName[i], "_rms") != NULL) {
	    ProfileValue[j][i] = sqrt(max(fabs(ProfileValue[j][i]),0.0)/1.0);
	    //	    fprintf(stderr, "%"GOUTSYM"\n ", ProfileValue[j][i]);
	  }
      }

    /* Compute the dark matter relaxation time due to the finite
       number of particles (from (8-71) in Binney & Tremaine). */
    
    for (j = 0; j < NumberOfPoints; j++)
      if (ProfileValue[j][32] > 1) {
	ProfileValue[j][56] = 0.34 * pow(ProfileValue[j][31]*1e5, 3) /
	  (pow(GravConst, 2) * 
	   (ProfileValue[j][30]*ProfileWeight[j][30] / 
	    ProfileValue[j][32]) * SolarMass *
	   ProfileValue[j][30] * SolarMass / (Mpc*Mpc*Mpc) *
	   log(0.4*ProfileValue[j][32]) );
	ProfileName[56] = "T_relax (s)";
      }

    /* Compute the cumulative luminosity and luminosity-weighted
       temperatures. */

    float CumLum = 0, CumTemp = 0, CumClump = 0;
    for (j = 0; j < NumberOfPoints; j++) {
      CumLum += ProfileValue[j][3];
      CumTemp += ProfileValue[j][84]*ProfileValue[j][3];
      ProfileValue[j][85] = CumLum;
      ProfileValue[j][86] = CumTemp/CumLum;
    }
    ProfileName[85] = "xray_cumulative_luminosity (erg/s)";
    ProfileName[86] = "temp_gas_cumulative_xray_weighted (K)";

    /* Compute the clumping factor and the cumulative clumping factor. */

    for (j = 0; j < NumberOfPoints; j++) {
      if (ProfileValue[j][0] > 0)
	ProfileValue[j][87] /= ProfileValue[j][0]*ProfileValue[j][0];
      CumClump += ProfileValue[j][87]*(ProfileValue[j][29] - 
	       ( (j==0) ? 0 : ProfileValue[j-1][29]));
      ProfileValue[j][88] = CumClump/ProfileValue[j][29];
    }
    ProfileName[88] = "cumulative clumping factor";

    /* ------------------------------------------------------------ */
    /* Compute properties averaged over rvir. */  

    FLOAT RvirRadius[2], RvirValue[1][MAX_PROFILES],  
          RvirWeight[1][MAX_PROFILES];

    /* Allocate and clear rvir data. */

    RvirRadius[0] = 0;
    RvirRadius[1] = rvir;
    for (i = 0; i < MAX_PROFILES; i++) {
      RvirValue[0][i] = 0;    
      RvirWeight[0][i] = 0;
    }

    if (rvir > 0) {

      /* Find values within rvir. */

      if (debug) printf("->finding rvir average\n");  
      for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
	LevelHierarchyEntry *Temp = LevelArray[level];
	while (Temp != NULL) {
	  Temp->GridData->AddToRadialProfile(Center, rvir, MeanVelocity, 1,
					     RvirRadius, RvirValue, RvirWeight,
					     ProfileName, &parameters);
	  Temp = Temp->NextGridThisLevel;
	}
      }
      CommunicationAllSumValues(&RvirValue[0][0], MAX_BINS);
      CommunicationAllSumValues(&RvirWeight[0][0], MAX_BINS);

      /* Set gas, dm and star masses within rvir. */

      mvir_gas = RvirValue[0][29] = RvirValue[0][0];
      RvirValue[0][57] = RvirValue[0][30];
      mvir_star = RvirValue[0][60] = RvirValue[0][60];

      /* If the RvirWeight is < 0, it is replaced by the annulus volume. */

      for (i = 0; i < MAX_PROFILES; i++)
	if (RvirWeight[0][i] < 0) 
	  RvirWeight[0][i] = pow(BoxSize, 3) * 4.0*pi/3.0 *
	    (pow(RvirRadius[1], 3) - pow(RvirRadius[0], 3));

      /* Normalize profile values by dividing each by cumulative weight. */

      for (i = 0; i < MAX_PROFILES; i++) {

	if (RvirWeight[0][i] > 0) 
	  RvirValue[0][i] /= RvirWeight[0][i];

	/* If this is a radial velocity dispersion, subtract the mean velocity 
	   which MUST be the profile immediately before. */

	if (ProfileName[i] != NULL)
	  if (strstr(ProfileName[i], "vr_rms") != NULL) {
	    RvirValue[0][i] -= pow(RvirValue[0][i-1], 2);
	  }

	/* If this is an rms profile, convert from 3D rms^2 to 3D rms. */   

	if (ProfileName[i] != NULL)
	  if (strstr(ProfileName[i], "_rms") != NULL)
	    RvirValue[0][i] = sqrt(max(RvirValue[0][i],0.0)/1.0);

      }

      /* Compute the spin parameter.  Note: this is usually defined as    
	 L E^1/2 / GM^5/2, but here is l e^1/2 / GM  where l (e) is the
	 specific angular momentum (specific energy).*/

      //Changed to sqrt(0.5) * sqrt(mv^2) by Matt Turk, 2008-02-12"
      float ang_mom, SpinUnits = Mpc * 1.0e5 * 1.0e5 /
	                         (GravConst * SolarMass);
      ang_mom = sqrt(RvirValue[0][17]*RvirValue[0][17] +
		     RvirValue[0][18]*RvirValue[0][18] +
		     RvirValue[0][19]*RvirValue[0][19] );
      spin_gas = SpinUnits * ang_mom * sqrt(0.5) *RvirValue[0][1] / mvir;
      //spin_gas = SpinUnits * ang_mom *RvirValue[0][1] / mvir;

      ang_mom = sqrt(RvirValue[0][35]*RvirValue[0][35] +
		     RvirValue[0][36]*RvirValue[0][36] +
		     RvirValue[0][37]*RvirValue[0][37] );
      spin_dm = SpinUnits * ang_mom * sqrt(0.5) *RvirValue[0][31] / mvir;
      //spin_dm = SpinUnits * ang_mom *RvirValue[0][31] / mvir;

      /* Compute clumping factor. */

      RvirValue[0][87] /= RvirValue[0][0]*RvirValue[0][0];
    } // end: if (rvir > 0)

    /* ------------------------------------------------------------ */
    /* If requested, compute special clumping factor. */

    if (parameters.ComputeClumpingFactor)
      if (AnalyzeClusterComputeClumpingFactor(LevelArray, &MetaData,
		 nint(RvirValue[0][4]), nint(RvirValue[0][32]), 
                 Center, rvir, Name) == FAIL) {
	fprintf(stderr, "Error in AnalyzeClusterComputeClumpingFactor\n");
	my_exit(EXIT_FAILURE);
      }

    /* ------------------------------------------------------------ */
    /* If requested, compute disk properties. */

    if (rvir > 0 && parameters.ComputeDiskInformation) {  

      /* Allocate space. */

      int const nimages = 6;
      int image_size = parameters.DiskImageSize*parameters.DiskImageSize;
      FLOAT *DiskImage[nimages];
      for (j = 0; j < nimages; j++) {
	DiskImage[j] = new FLOAT[image_size];
	for (i = 0; i < image_size; i++)
	  DiskImage[j][i] = 0;
      }

      /* Compute unit disk vector. */

     float length = sqrt(RvirValue[0][80]*RvirValue[0][80] +
		    RvirValue[0][81]*RvirValue[0][81] +
		    RvirValue[0][82]*RvirValue[0][82]);
      for (dim = 0; dim < 3; dim++)
	DiskVector[dim] = RvirValue[0][80+dim]/length;  

      /* Get info from grids. */

      if (debug) printf("->finding disk profile\n");

      for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
	LevelHierarchyEntry *Temp = LevelArray[level];
	while (Temp != NULL) {
	  Temp->GridData->AddToDiskProfile(Center, rvir, MeanVelocity,   
				       NumberOfPoints,
				       ProfileRadius,
				       ProfileValue, ProfileWeight,
				       ProfileName, &parameters,
				       DiskVector, DiskImage,
				       parameters.DiskImageSize, 
					   parameters.DiskRadius);
	  Temp = Temp->NextGridThisLevel;
	}
      }
      CommunicationSumValues(&ProfileValue[0][100], MAX_BINS*6);
      CommunicationSumValues(&ProfileWeight[0][100], MAX_BINS*6);
      for (j = 0; j < nimages; j++)
	CommunicationSumValues(DiskImage[j], image_size);

      if (MyProcessorNumber == ROOT_PROCESSOR) {

	/* If the RvirWeight is < 0, it is replaced by the annulus volume. */  

      for (i = 100; i < 120; i++)
	for (j = 0; j < NumberOfPoints; j++)
	  if (ProfileWeight[j][i] < 0)             
	    ProfileWeight[j][i] = pow(BoxSize, 2) * pi *
	      (pow(ProfileRadius[j+1], 2) - pow(ProfileRadius[j], 2));

      /* Normalize profile values by dividing each by cumulative weight. */  

      for (i = 100; i < 120; i++)
	for (j = 0; j < NumberOfPoints; j++)
	  if (ProfileWeight[j][i] > 0) 
	    ProfileValue[j][i] /= ProfileWeight[j][i];  

#ifdef USE_HDF4
      /* Save disk image. */

      Eint32 OutDims[2];
      OutDims[0] = OutDims[1] = parameters.DiskImageSize;
      if (DFSDsetdims(2, OutDims) == HDF_FAIL) {
	fprintf(stderr, "Error in DFSDsetdims.\n");
	return FAIL;
      }

      strcpy(DiskImageName, Name);
      strcat(DiskImageName, ".DiskImage");

      float32 *float_temp = new float32[image_size];
      for (j = 0; j < nimages; j++) {
	
	for (i = 0; i < image_size; i++)
	  float_temp[i] = float32(DiskImage[j][i]);

	if (j == 0)
	  ret = DFSDputdata(DiskImageName, 2, OutDims, (VOIDP) float_temp);
	if (j != 0)
	  ret = DFSDadddata(DiskImageName, 2, OutDims, (VOIDP) float_temp);

	if (ret == HDF_FAIL) {
	  fprintf(stderr, "Error in DFSDput/adddata.\n");
	  return FAIL;
	}

      }

      /* Clean up. */

      delete float_temp;
#endif  //USE_HDF4

      } // end: if (MyProcessorNumber == ROOT_PROCESSOR)

      for (j = 0; j < nimages; j++)
	delete DiskImage[j];

      /* ------------------------------------------------------------ */
      /* While calculating disk profile, ompute vertical profile as well
	 using ProfileRadius2 defined above.  - Ji-hoon Kim */

      if (debug) printf("->finding vertical profile\n");

      for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
	LevelHierarchyEntry *Temp = LevelArray[level];
	while (Temp != NULL) {
	  if(parameters.LinearProfileRadiusForVertical){
	    Temp->GridData->AddToVerticalProfile(Center, rvir, MeanVelocity,   
				       NumberOfPoints,
				       ProfileRadius2,   
				       ProfileValue, ProfileWeight,
				       ProfileName, &parameters,
				       DiskVector);  
	  } else {
	    Temp->GridData->AddToVerticalProfile(Center, rvir, MeanVelocity,   
				       NumberOfPoints,
				       ProfileRadius,   
				       ProfileValue, ProfileWeight,
				       ProfileName, &parameters,
				       DiskVector);  
	  }
	  Temp = Temp->NextGridThisLevel;
	}
      }
      CommunicationSumValues(&ProfileValue[0][100], MAX_BINS*6);
      CommunicationSumValues(&ProfileWeight[0][100], MAX_BINS*6);

      if (MyProcessorNumber == ROOT_PROCESSOR) {

      /* If the RvirWeight is < 0, it is replaced by the slab edge-on area */  

      for (i = 150; i < 190; i++)
	for (j = 0; j < NumberOfPoints; j++)
	  if (ProfileWeight[j][i] < 0)             
	    if(parameters.LinearProfileRadiusForVertical){
	      ProfileWeight[j][i] = parameters.DiskRadiusCutoff *
		BoxSize * (ProfileRadius2[j+1] - ProfileRadius2[j]);   
	    } else {
	      ProfileWeight[j][i] = parameters.DiskRadiusCutoff *
		BoxSize * (ProfileRadius[j+1] - ProfileRadius[j]);   
	    }

      /* Normalize profile values by dividing each by cumulative weight. */  

      for (i = 150; i < 190; i++)
	for (j = 0; j < NumberOfPoints; j++)
	  if (ProfileWeight[j][i] > 0) 
	    ProfileValue[j][i] /= ProfileWeight[j][i];  

      } // end: if (MyProcessorNumber == ROOT_PROCESSOR)
      
    } //end: if (rvir > 0 && parameters.ComputeDiskInformation) 
  

    /* ------------------------------------------------------------ */ 
    /* Output the data: naming convention:
       file 0: Name
       file 1: Name.Inertial
       file 2: Name.Species
       file 3: Name.DarkMatter
       file 4: Name.StarParticles
       file 5: Name.Disk            
       file 6: Name.Vertical */

    if (MyProcessorNumber == ROOT_PROCESSOR) {
    
#define NUMBER_OF_FILES 7 

    FILE *fptrs[NUMBER_OF_FILES];
    char *FileName[NUMBER_OF_FILES];
    int   ColumnNumber[NUMBER_OF_FILES];

    /* Generate output file names. */

    for (i = 0; i < NUMBER_OF_FILES; i++) {
      ColumnNumber[i] = 3;
      FileName[i] = new char[MAX_LINE_LENGTH];
      strcpy(FileName[i], Name);
    }
    strcat(FileName[1], ".Inertial");  
    strcat(FileName[2], ".Species");
    strcat(FileName[3], ".DarkMatter");
    strcat(FileName[4], ".StarParticles");
    strcat(FileName[5], ".Disk");
    strcat(FileName[6], ".Vertical"); 

    /* Open output files. */

    for (i = 0; i < NUMBER_OF_FILES; i++) {
      fptrs[i] = NULL;
      if ((i != 4 || StarParticleCreation > 0) &&
	  (i != 5 || parameters.ComputeDiskInformation == TRUE))
	fptrs[i] = fopen(FileName[i], "w");
    }

    /* Generate an array of file pointers so each profile knows which   
       file it should be sent to. */

    int ProfileFile[MAX_PROFILES];
    for (profile = 0; profile < 8; profile++)
      ProfileFile[profile] = 0;
    for (profile = 8; profile < 17; profile++)
      ProfileFile[profile] = 2;
    for (profile = 17; profile < 30; profile++)
      ProfileFile[profile] = 0;
    for (profile = 25; profile < 28; profile++)
      ProfileFile[profile] = 2;
    for (profile = 40; profile < 46; profile++)
      ProfileFile[profile] = 1;
    ProfileFile[46] = 0;
    for (profile = 30; profile < 40; profile++)
      ProfileFile[profile] = 3;
    for (profile = 50; profile < 56; profile++)
      ProfileFile[profile] = 1;
    ProfileFile[56] = ProfileFile[57] = 3;
    if (StarParticleCreation > 0)
      for (profile = 60; profile < 70; profile++)
	ProfileFile[profile] = 4;
    for (profile = 80; profile < 89; profile++)
      ProfileFile[profile] = 0;
    ProfileFile[83] = 2; 
    ProfileFile[90] = 2; 
    for (profile = 100; profile < 120; profile++)  
      ProfileFile[profile] = 5;
    for (profile = 120; profile < 130; profile++)
      ProfileFile[profile] = 2;
    for (profile = 150; profile < 190; profile++)  
      ProfileFile[profile] = 6;
    for (profile = 191; profile < 194; profile++)  
      ProfileFile[profile] = 0;

    /* Output global values. */

    fprintf(fptrs[0], "# Center             = %"GOUTSYM" %"GOUTSYM" %"GOUTSYM"\n", Center[0], 
	  Center[1], Center[2]);
    fprintf(fptrs[0], "# MaxBaryonValue     = %g\n", MaxDensity);
    //fprintf(fptrs[0], "# MinCoolingTime     = %g\n", MinCoolingTime);
    //fprintf(fptrs[0], "# CoolingCenter      = %"GOUTSYM" %"GOUTSYM" %"GOUTSYM"\n",
          //CTCenter[0], CTCenter[1], CTCenter[2]);
    fprintf(fptrs[0], "# MeanVelocity (gas) = %"GOUTSYM" %"GOUTSYM" %"GOUTSYM" (km/s)\n", 
	  MeanVelocity[0][0], MeanVelocity[1][0], MeanVelocity[2][0]);
    fprintf(fptrs[0], "# MeanVelocity (dm ) = %"GOUTSYM" %"GOUTSYM" %"GOUTSYM" (km/s)\n", 
	  MeanVelocity[0][1], MeanVelocity[1][1], MeanVelocity[2][1]);
    fprintf(fptrs[0], "# MeanVelocity (tot) = %"GOUTSYM" %"GOUTSYM" %"GOUTSYM" (km/s)\n",
	  MeanVelocity[0][2], MeanVelocity[1][2], MeanVelocity[2][2]);
    fprintf(fptrs[0], "#  (Within r = %"GOUTSYM" Mpc)\n#\n",  
	    MeanVelocityOuterEdge*BoxSize);
    fprintf(fptrs[0], "# L (gas)            = %"GOUTSYM" %"GOUTSYM" %"GOUTSYM" (Mpc km/s)\n",
	    RvirValue[0][17], RvirValue[0][18], RvirValue[0][19]);
    fprintf(fptrs[0], "# L (dm)             = %"GOUTSYM" %"GOUTSYM" %"GOUTSYM" (Mpc km/s)\n",
	    RvirValue[0][35], RvirValue[0][36], RvirValue[0][37]);
    fprintf(fptrs[0], "# L (star)           = %"GOUTSYM" %"GOUTSYM" %"GOUTSYM" (Mpc km/s)\n",
	    RvirValue[0][65], RvirValue[0][66], RvirValue[0][67]);
    fprintf(fptrs[0], "# InnerEdge = %g Mpc  OuterEdge = %g Mpc\n\n",    
	  InnerEdge*BoxSize, OuterEdge*BoxSize);

    /* Output general information about simulation. */

    fprintf(fptrs[0], "# OmegaMatterNow     = %g\n", OmegaMatterNow);
    fprintf(fptrs[0], "# OmegaLambdaNow     = %g\n", OmegaLambdaNow);
    fprintf(fptrs[0], "# HubbleConstantNow  = %g\n", HubbleConstantNow);
    fprintf(fptrs[0], "# CurrentRedshift    = %"GOUTSYM"\n", CurrentRedshift);
    fprintf(fptrs[0], "# AMR file name      = %s\n\n", argv[1]);
    if (parameters.MetaData != NULL) {
      fprintf(fptrs[0], "# MetaData           = %s\n\n", parameters.MetaData);
    }

    /* Output rvir information. */

    fprintf(fptrs[0], "# rvir               = %g (Mpc)\n", rvir*BoxSize);
    fprintf(fptrs[0], "# r500               = %g (Mpc)\n", r500*BoxSize);
    fprintf(fptrs[0], "# mvir               = %g (M_solar)\n", mvir);
    fprintf(fptrs[0], "# VirialDensity      = %g (to critical density)\n", 
	    parameters.virial_dens);
    fprintf(fptrs[0], "# mvir (gas,dm,star) = %g %g %g (M_solar)\n", 
	    mvir_gas, mvir - mvir_gas - mvir_star, mvir_star);
    fprintf(fptrs[0], "# spin (gas, dm)     = (%g, %g)\n", spin_gas, spin_dm); 
    fprintf(fptrs[0], "# Tvir predicted (K) = %g\n", Tvir);
    if (parameters.XrayTableFileName != NULL) {
      fprintf(fptrs[0], "# XrayLowerCutoffkeV = %g\n", 
	      parameters.XrayLowerCutoffkeV);
      fprintf(fptrs[0], "# XrayUpperCutoffkeV = %g\n", 
	      parameters.XrayUpperCutoffkeV);
      fprintf(fptrs[0], "# XrayTableFileName  = %s\n", 
	      parameters.XrayTableFileName);
    }
    fprintf(fptrs[0], "# ColdTempCutoff (K) = %g\n\n", 
	    parameters.ColdTemperatureCutoff);
    if (parameters.ComputeDiskInformation)
        fprintf(fptrs[0], "# Disk Vector = %f %f %f\n", DiskVector[0], DiskVector[1], DiskVector[2]);

    for (i = 0; i < NUMBER_OF_FILES; i++)
      if (fptrs[i] != NULL)
	fprintf(fptrs[i], "#");

    for (profile = 0; profile < MAX_PROFILES; profile++)
      if (ProfileName[profile] != NULL) {
	    if (debug) {
          fprintf(stderr,"point %i Rvir: %g \n ", profile, RvirValue[0]);
        }
      	fprintf(fptrs[ProfileFile[profile]], "%"GOUTSYM" ", RvirValue[0][profile]); 
      }
    
    for (i = 0; i < NUMBER_OF_FILES; i++)
      if (fptrs[i] != NULL)
	fprintf(fptrs[i], "\n#\n");

    /* Print column info. */

    for (i = 0; i < NUMBER_OF_FILES; i++)
      if (fptrs[i] != NULL) {
	fprintf(fptrs[i], "#\n# COLUMN   DESCRIPTION\n");
	fprintf(fptrs[i], "#   %3d    %s\n", 1, "bin central radius (Mpc)");
	fprintf(fptrs[i], "#   %3d    %s\n", 2, "bin outer radius (Mpc)");
      }

    for (profile = 0; profile < MAX_PROFILES; profile++)     
      if (ProfileName[profile] != NULL) {
	fprintf(fptrs[ProfileFile[profile]], "#   %3d    %s\n", 
		ColumnNumber[ProfileFile[profile]]++,
		ProfileName[profile]);
      }

    /* Print header. */
  
    for (i = 0; i < NUMBER_OF_FILES; i++)
      if (fptrs[i] != NULL)
	fprintf(fptrs[i], "#\n#%12.12s %12.12s ", "r_cent (Mpc)", 
		"r_edge (Mpc)");

    for (profile = 0; profile < MAX_PROFILES; profile++)
      if (ProfileName[profile] != NULL)
	fprintf(fptrs[ProfileFile[profile]], "%12.12s ", ProfileName[profile]);
        fprintf(fptrs[0], "%12.12s %12.12s %12.12s ", "vol_gas (real)", 
		"vol (calc)", "vol_dm (real)");    

    for (i = 0; i < NUMBER_OF_FILES; i++)
      if (fptrs[i] != NULL)
	fprintf(fptrs[i], "\n");

    /* Loop over all radial bins, printing each. */

    for (j = 0; j < NumberOfPoints; j++) {

      FLOAT rmid = 0.5*(ProfileRadius[j]+ProfileRadius[j+1])*BoxSize;            

      for (i = 0; i < NUMBER_OF_FILES; i++)
	if (fptrs[i] != NULL) 
	  fprintf(fptrs[i], "%"GOUTSYM" %"GOUTSYM" ", 
	          rmid, ProfileRadius[j+1]*BoxSize);

      for (profile = 0; profile < MAX_PROFILES; profile++)  
	if (ProfileName[profile] != NULL)
	  fprintf(fptrs[ProfileFile[profile]], "%"GOUTSYM" ", 
		  (ProfileValue[j][profile] == 0) ? 
		  tiny_number : ProfileValue[j][profile]);  

      fprintf(fptrs[0], "%"GOUTSYM" %"GOUTSYM" %"GOUTSYM" ", ProfileWeight[j][0],   
	      ((POW(ProfileRadius[j+1], 3) - POW(ProfileRadius[j], 3) ) * 
	      POW(BoxSize, 3) * 4.0*pi/3.0),
	      ProfileWeight[j][30]);

      /* Print vertical profile using ProfileRadius2.  JHK in Dec.2007 */

      if ((fptrs[6] != NULL) && (parameters.LinearProfileRadiusForVertical)) {
	FLOAT rmid = 0.5*(ProfileRadius2[j]+ProfileRadius2[j+1])*BoxSize;                   
	fprintf(fptrs[6], "%"GOUTSYM" %"GOUTSYM" ", 
		rmid, ProfileRadius2[j+1]*BoxSize); 
      }

      for (i = 0; i < NUMBER_OF_FILES; i++)
	if (fptrs[i] != NULL)
	  fprintf(fptrs[i], "\n");

    }//end of j

    /* close files */

    for (i = 0; i < NUMBER_OF_FILES; i++)
      if (fptrs[i] != NULL)
	fclose(fptrs[i]);


    /* Print additional info for Kennicutt-Schmidt relations, JHK in Nov.2007 */

    if (parameters.PrintGlobalProfileValues) { 

      FILE *fpLKS, *fpGKS, *fpSSD;
      double AverageGasSurfaceDensity = 0.0;  
      double AverageSFRSurfaceDensity = 0.0;  
      FLOAT ObservableDiskRadiusForKS = ProfileRadius[NumberOfPoints-1];

      fpLKS = fopen("local_KS.dat","w");      
      fpGKS = fopen("global_KS.dat","a");
      fpSSD = fopen("stellar_surface_density.dat", "a");
//      fprintf(fpGKS, "# time  Sigma_gas  Sigma_SFR  R_vir  M_vir  M_vir_star  M_vir_gas  M_vir_DM\n");

      /* Print gas surface density and SFR surface density for local K-S law.  
	 Msun/pc^2 vs. Msun/yr/kpc^2 */

      for (j = 1; j < NumberOfPoints; j++) {

	fprintf(fpLKS, "%"GOUTSYM"   %"GOUTSYM"   %"GOUTSYM" \n", ProfileRadius[j],
		log10( (ProfileValue[j][106] == 0) ? tiny_number : ProfileValue[j][106]/POW(1e6, 2) ),
		log10( (ProfileValue[j][116] == 0) ? tiny_number : ProfileValue[j][116]/POW(1e3, 2) )); 

	AverageGasSurfaceDensity += ProfileValue[j][106]*ProfileWeight[j][106];//weights are annuli
	AverageSFRSurfaceDensity += ProfileValue[j][116]*ProfileWeight[j][116];

      }

      /* Print disk-averaged gas surface density and disk-averaged SFR surface
	 density for global K-S law.  Msun/pc^2 vs. Msun/yr/kpc^2  */

      fprintf(fpGKS, "%g   %"GOUTSYM"   %"GOUTSYM"    %g    %g    %g    %g    %g\n", 
	      MetaData.Time,
	      log10( AverageGasSurfaceDensity/pi/POW((ObservableDiskRadiusForKS*BoxSize),2)/POW(1e6, 2) ),
	      log10( AverageSFRSurfaceDensity/pi/POW((ObservableDiskRadiusForKS*BoxSize),2)/POW(1e3, 2) ),
	      rvir*BoxSize, mvir, mvir_star, mvir_gas, mvir - mvir_gas - mvir_star); 

      /* Print stellar/gas/SFR surface density w.r.t. time and radius (in Msun/pc^2, Msun/yr/kpc^2) */

      for (j = 1; j < NumberOfPoints; j++) {

	if (ProfileValue[j][114] != 0)
	  fprintf(fpSSD, "%g   %"GOUTSYM"   %"GOUTSYM"   %"GOUTSYM"   %"GOUTSYM"\n", 
		  MetaData.Time, ProfileRadius[j],
		  log10( (ProfileValue[j][114] == 0) ? tiny_number : ProfileValue[j][114]/POW(1e6, 2) ),
		  log10( (ProfileValue[j][106] == 0) ? tiny_number : ProfileValue[j][106]/POW(1e6, 2) ),
		  log10( (ProfileValue[j][116] == 0) ? tiny_number : ProfileValue[j][116]/POW(1e3, 2) ));

      }

      fclose(fpLKS);
      fclose(fpGKS);
      fclose(fpSSD);
      
    } // end: if (parameters.PrintGlobalProfileValues)

    } // end: if (MyProcessorNumber == ROOT_PROCESSOR)

  } // end: loop over centers

  my_exit(EXIT_SUCCESS);

}


void my_exit(int status)
{
  CommunicationFinalize();
  exit(status);
}

void FindVirialRadius(LevelHierarchyEntry *LevelArray[], float &rvir, float &r500,
		      float critical_density, float BoxSize,
		      FLOAT *Center, FLOAT OuterEdge, float MeanVelocity[MAX_DIMENSION][3], 
		      int NumberOfPoints, float *ProfileRadius, float ProfileValue[][MAX_PROFILES],
		      float ProfileWeight[][MAX_PROFILES], char *ProfileName[MAX_PROFILES],
		      AnalyzeClusterParameters *parameters)
{
  int i, j, level;
  const float pi = 3.14159;


  for (i = 0; i < MAX_PROFILES; i++) {
    for (j = 0; j < MAX_BINS; j++) {
      ProfileValue[j][i] = 0.0;
      ProfileWeight[j][i] = 0.0;
    }
    ProfileName[i] = NULL;
  }


  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {
    printf("level %d\n", level);
    LevelHierarchyEntry *Temp = LevelArray[level];
    while (Temp != NULL) {
      if (Temp->GridData->AddToRadialProfile(Center, OuterEdge, MeanVelocity,
					     NumberOfPoints, ProfileRadius,
					     ProfileValue, ProfileWeight,
					     ProfileName, parameters) == FAIL) {
	fprintf(stderr, "Error in grid->AddToRadialProfile.\n");
	exit(EXIT_FAILURE);
      }
      Temp = Temp->NextGridThisLevel;
    }
  }
  CommunicationAllSumValues(&ProfileValue[0][0], MAX_BINS*MAX_PROFILES);
  CommunicationAllSumValues(&ProfileWeight[0][0], MAX_BINS*MAX_PROFILES);

  /* Compute rvir based on critical density (as of Feb 21/2000). */

  for (j = 0; j < NumberOfPoints; j++) {

    /* Compute the mass (gas+dm+star), in M(solar), within this annulus. */
    
    ProfileValue[j][38] = ProfileValue[j][0] + ProfileValue[j][30] +
	                  ProfileValue[j][60];

    /* Keep a running sum of the total mass within radius j+1. */

    if (j > 0)
      ProfileValue[j][38] += ProfileValue[j-1][38];

    /* Compute the average overdensity within the sphere with radius [j+1].*/

    ProfileValue[j][39] = ProfileValue[j][38] / 
	(pow(ProfileRadius[j+1]*BoxSize, 3) * 4.0*pi/3.0) / critical_density;
  }


  /* Find the radius at which the cumulative overdensity goes through 
     virial_dens (typically 200); also compute r500 while we're at it.
     Update: 06.01.2010 CBH - Start from outside and move into center 
     to calculate rvir.  This assure you have sufficient mass in center
     (even if you aren't *right* on a halo). */

    if (ProfileValue[NumberOfPoints-1][39] >= parameters->virial_dens) {
        fprintf(stderr, "Error: Choose a larger OuterEdge; rvir > OuterEdge (%g)\n", OuterEdge);
        my_exit(EXIT_FAILURE);
    }
    for (j = NumberOfPoints-1; j >= 0; j--) {

    /* Check to see if the overdensity has increased beyond virial_dens. */
    if (ProfileValue[j][39] >= parameters->virial_dens && rvir == 0) {

        /* interpolate in log-space since overdensity ~ r^alpha */

        rvir = log(ProfileRadius[j+1]) +
            (log(ProfileRadius[j+2])     - log(ProfileRadius[j+1])     ) *
            (log(parameters->virial_dens) - log(ProfileValue[j][39])) /
            (log(ProfileValue[j+1][39])    - log(ProfileValue[j][39]) +
            tiny_number);

        rvir = exp(rvir);
    }

    /* Check to see if the overdensity has increased beyond 500. */

    if (ProfileValue[j][39] >= 500 && r500 == 0) {
        r500 = log(ProfileRadius[j+1]) +
               (log(ProfileRadius[j+2])     - log(ProfileRadius[j+1])     ) *
               (log(500.0)                  - log(ProfileValue[j][39])) /
               (log(ProfileValue[j+1][39])    - log(ProfileValue[j][39]) +
                tiny_number);
        r500 = exp(r500);
        break;
    }
  }
  return;
}
