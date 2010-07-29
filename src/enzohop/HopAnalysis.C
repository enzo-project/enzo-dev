/***********************************************************************
/
/  PERFORMS A "HOP" CLUSTERING ANALYSIS (astro-ph/9712200)
/
/  written by: Greg Bryan
/  date:       February, 1999
/  modified1: chummels 5.24.2010 - added additional output of HopParticles.out
/     for use in tracking halos from timestep to timestep. must define: 
/     PARTICLE_OUTPUT.
/
/  PURPOSE:
/
************************************************************************/

//
//
// Define PARTICLE_OUTPUT for output of HopParticles.out file.
// #define PARTICLE_OUTPUT
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef USE_HDF4
#include <df.h>
#endif /* USE_HDF4 */
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
#include "../enzo/macros_and_parameters.h"
#include "../enzo/typedefs.h"
#define DEFINE_STORAGE
#include "../enzo/ErrorExceptions.h"
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
#include "../enzo/units.h"
#include "../enzo/flowdefs.h"
#include "../enzo/PhotonCommunication.h"
#undef DEFINE_STORAGE
extern "C" {
#include "kd.h"
}

/* function prototypes */

int ReadAllData(char *filename, HierarchyEntry *TopGrid, TopGridData &tgd,
		ExternalBoundary *Exterior, float *Inititaldt);
int Group_ReadAllData(char *filename, HierarchyEntry *TopGrid, TopGridData &tgd,
		      ExternalBoundary *Exterior, float *Initialdt,
		      bool ReadParticlesOnly=false);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int SetDefaultGlobalValues(TopGridData &MetaData);
int  DepositParticleMassField(HierarchyEntry *Grid, float RequestTime = -1.0);
int  CopyOverlappingZones(grid* CurrentGrid, TopGridData *MetaData, 
			 LevelHierarchyEntry *LevelArray[], int level);
int  CopyOverlappingParticleMassFields(grid* CurrentGrid, 
				      TopGridData *MetaData, 
				      LevelHierarchyEntry *LevelArray[], 
				      int level);
int CommunicationInitialize(Eint32 *argc, char **argv[]);
int CommunicationFinalize();
void my_exit(int status);
void hop_main(KD kd);
void regroup_main(float dens_outer);
int kdInit(KD *kd, int nBucket);
Eint32 hide_isdigit(Eint32 c);

main(int argc, char *argv[])
{
  CommunicationInitialize(&argc, &argv);

  /* Main declarations */

  TopGridData MetaData;
  HierarchyEntry TopGrid;
  ExternalBoundary Exterior;
  LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];

  /* Initialize */

  float HopDensityThreshold    = 160;   // reproduces FOF b=0.2
  debug                        = FALSE;
  char *ParameterFile          = NULL;
  char *myname                 = argv[0];
  FLOAT RegionLeft[MAX_DIMENSION], RegionRight[MAX_DIMENSION];
  int level, dim, i, j, part, UseParticleType[NUM_PARTICLE_TYPES];
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    LevelArray[level] = NULL;
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    RegionLeft[dim] = RegionRight[dim] = FLOAT_UNDEFINED;
  for (i = 0; i < NUM_PARTICLE_TYPES; i++)
    UseParticleType[i] = 0;
#ifdef USE_MPI
  if (NumberOfProcessors > 1)
    ENZO_FAIL("Can only use one processor.");
#endif /* USE_MPI */

  /* --------------------------------------------------------------- */
  /* Interpret command-line arguments. */

  char c;
  while (--argc > 0 && (*++argv)[0] == '-')
    while (c = *++argv[0])
      switch (c) {

	/* get beginning of region selection (float coordinate). */

      case 'b':
	dim = 0;
	while (dim < MAX_DIMENSION && argc > 1 && hide_isdigit(*argv[1])) {
	  argc--;
	  if (sscanf((*++argv), "%"FSYM, &RegionLeft[dim++]) != 1) {
	    fprintf(stderr, "%s: error reading Begin coordinates.\n", myname);
	    my_exit(EXIT_FAILURE);
	  }
	  while (*(argv[0]+1))
	    ++argv[0];
	}
	break;

	/* debug */

      case 'd':
	debug = TRUE;
	break;

	/* get finish of region selection (float coordinate). */

      case 'f':
	dim = 0;
	while (dim < MAX_DIMENSION && argc > 1 && hide_isdigit(*argv[1])) {
	  argc--;
	  if (sscanf((*++argv), "%"FSYM, &RegionRight[dim++]) != 1) {
	    fprintf(stderr, "%s: error reading Finish coordinates.\n", myname);
	    my_exit(EXIT_FAILURE);
	  }
	  while (*(argv[0]+1))
	    ++argv[0];
	}
	break;

	/* Use gas particles */

      case 'g':
	UseParticleType[0] = TRUE;
	break;

	/* Use dark matter particles */

      case 'm':
	UseParticleType[1] = TRUE;
	break;

	/* Star particles */

      case 's':
	UseParticleType[2] = TRUE;
	break;

	/* get density threshold. */

      case 't':
	if (--argc > 0) {
	  sscanf((*++argv), "%"FSYM, &HopDensityThreshold);
	  while (*(argv[0]+1))
	    ++argv[0];
	}
	break;

	/* Unknown */

      default:
	fprintf(stderr, "%s: unknown command-line option: -%s.\n", myname, &c);
	my_exit(EXIT_FAILURE);
	
      } // end of switch(c)

  /* Error check for number of parameters, and set parameter file. */

  if (argc != 1) {
    fprintf(stderr, "enzohop [-b #] [-f #] [-t #] [-g] [-s] [-m] [-d] amr_file\n");
    fprintf(stderr, "  -b)egin region\n");
    fprintf(stderr, "  -f)inish region\n");
    fprintf(stderr, "  -t)hreshold for hop (default 160)\n");
    fprintf(stderr, "  -g)as particles used in hop analysis\n");
    fprintf(stderr, "  -s)tar particles used in hop analysis\n");
    fprintf(stderr, "  -m) dark matter particles used in hop analysis (default)\n");
    fprintf(stderr, "  -d)ebug\n");
    my_exit(FAIL);
  }
  ParameterFile = argv[0];

  /* If none set, use dark matter only. */

  int count = 0;
  for (i = 0; i < NUM_PARTICLE_TYPES; i++)
    count += UseParticleType[i];
  if (count == 0)
    UseParticleType[1] = 1;

  /* If using dark matter, then set must-refine particles also. */

  if (UseParticleType[1])
    UseParticleType[4] = 1;

  /* Read the saved file. */

  SetDefaultGlobalValues(MetaData);


  // First expect to read in packed-HDF5
  float dummy;
#ifdef USE_HDF5_GROUPS
  if (Group_ReadAllData(ParameterFile, &TopGrid, MetaData, &Exterior, &dummy) == FAIL) {
      if (MyProcessorNumber == ROOT_PROCESSOR) {
	fprintf(stderr, "Error in Group_ReadAllData %s\n", argv[1]);
	fprintf(stderr, "Probably not in a packed-HDF5 format. Trying other read routines.\n");
      }
#endif
      // If not packed-HDF5, then try usual HDF5 or HDF4
      if (ReadAllData(ParameterFile, &TopGrid, MetaData, &Exterior, &dummy) == FAIL) {
	if (MyProcessorNumber == ROOT_PROCESSOR) {
	  fprintf(stderr, "Error in ReadAllData %s.\n", argv[1]);
	}
	my_exit(EXIT_FAILURE);
      }
#ifdef USE_HDF5_GROUPS
    }
#endif

  AddLevel(LevelArray, &TopGrid, 0);    // recursively add levels

  /* Set the Cell width of the root grid. */

  float BaseRadius = (DomainRightEdge[0] - DomainLeftEdge[0])/
      float(MetaData.TopGridDims[0]);

  /* If undefined, set parameters. */

  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    if (RegionLeft[dim] == FLOAT_UNDEFINED)
      RegionLeft[dim] = DomainLeftEdge[dim];
    if (RegionRight[dim] == FLOAT_UNDEFINED)
      RegionRight[dim] = DomainRightEdge[dim];
  }

  if (debug)
    printf("HopAnalysis region: Left = %g %g %g   Right = %g %g %g\n",
	   RegionLeft[0], RegionLeft[1], RegionLeft[2],
	   RegionRight[0], RegionRight[1], RegionRight[2]);

  /* Initialize Particle List info */

  ListOfParticles *ListOfParticlesHead[NUM_PARTICLE_TYPES];
  int TotalNumberOfParticles[NUM_PARTICLE_TYPES];
  for (i = 0; i < NUM_PARTICLE_TYPES; i++) {
    ListOfParticlesHead[i] = NULL;
    TotalNumberOfParticles[i] = 0;
  }

  /* --------------------------------------------------------------- */
  /* Loop over all the levels, and collect particles */

  printf("Collecting particles...\n");
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++) {

    /* If SelfGravity, set all the particle mass fields. */

    LevelHierarchyEntry *Temp = LevelArray[level];
#ifdef UNUSED
    if (SelfGravity)
      while (Temp != NULL) {
	DepositParticleMassField(Temp->GridHierarchyEntry);
	Temp = Temp->NextGridThisLevel;
      }
#endif /* UNUSED */

    /* Loop over all the grids. */

    Temp = LevelArray[level];
    while (Temp != NULL) {

      /* Allocate a new ListOfParticles for this grid. */

      for (i = 0; i < NUM_PARTICLE_TYPES; i++) {
	ListOfParticles *Temp = ListOfParticlesHead[i];
	ListOfParticlesHead[i] = new ListOfParticles;
	ListOfParticlesHead[i]->NextList = Temp;
      }

      /* Set particle density. */

#ifdef UNUSED
      if (SelfGravity) {
	CopyOverlappingParticleMassFields(Temp->GridData, &MetaData, 
					  LevelArray, level);
	if (Temp->GridHierarchyEntry->ParentGrid != NULL)
	  Temp->GridHierarchyEntry->ParentGrid->GridData->DepositParticlePositions(Temp->GridData, Temp->GridHierarchyEntry->ParentGrid->GridData->ReturnTime(), GRAVITATING_MASS_FIELD_PARTICLES);
      }
#endif /* UNUSED */

      /* Initialize the UNDER_SUBGRID_FIELD for this grid. */

      Temp->GridData->ZeroSolutionUnderSubgrid(NULL, ZERO_UNDER_SUBGRID_FIELD);

      /* Zero the solution (on this grid) which is underneath any subgrid
	 (so we get only the high resolution solution from the subgrid). */

      LevelHierarchyEntry *Temp2 = LevelArray[level+1];
      while (Temp2 != NULL) {
	Temp->GridData->ZeroSolutionUnderSubgrid(Temp2->GridData, 
        					 ZERO_UNDER_SUBGRID_FIELD);
	Temp2 = Temp2->NextGridThisLevel;
      }

      /* Generate particle list for this grid. */

      Temp->GridData->OutputAsParticleData(RegionLeft, RegionRight,
					   ListOfParticlesHead, BaseRadius);

      /* Delete if not required. */

      for (i = 0; i < NUM_PARTICLE_TYPES; i++)
	if (UseParticleType[i])
	  TotalNumberOfParticles[i] += 
	    ListOfParticlesHead[i]->NumberOfParticles;
	else {
	  if (ListOfParticlesHead[i]->NumberOfParticles > 0) {
	    for (dim = 0; dim < MetaData.TopGridRank; dim++) {
	      delete [] ListOfParticlesHead[i]->ParticlePosition[dim];
	      delete [] ListOfParticlesHead[i]->ParticleVelocity[dim];
	    }
	    delete [] ListOfParticlesHead[i]->ParticleRadius;
	    for (j = 0; j < ListOfParticlesHead[i]->NumberOfValues; j++)
	      delete [] ListOfParticlesHead[i]->ParticleValue[j];
	  }
	}

      /* Delete grid. */

      delete Temp->GridData;

      /* Next grid on this level. */

      Temp = Temp->NextGridThisLevel;

    } // end loop over grids

  } // end loop over levels

  /* --------------------------------------------------------------- */
  /* Convert particles into kd structure. */

  if (debug)
    printf("TotalNumberOfParticles = %d %d %d %d %d\n", 
	   TotalNumberOfParticles[0], TotalNumberOfParticles[1], 
	   TotalNumberOfParticles[2], TotalNumberOfParticles[3],
	   TotalNumberOfParticles[4]);

  /* Initialize the kd structure. */

  if (sizeof(FLOAT) != 4) {
    fprintf(stderr, "error: hop is hard-coded for 4-byte floats.\n");
    my_exit(EXIT_FAILURE);
  }

  KD kd;
  int nBucket = 16, kdcount = 0;
  kdInit(&kd, nBucket);
  kd->nActive = 0;
  for (i = 0; i < NUM_PARTICLE_TYPES; i++)
    kd->nActive += TotalNumberOfParticles[i];
  kd->p = new PARTICLE[kd->nActive];
  if (kd->p == NULL) {
    fprintf(stderr, "failed allocating particles.\n");
    my_exit(EXIT_FAILURE);
  }

  float DensityFactor = 1.0/(2.78e11*OmegaMatterNow*pow(HubbleConstantNow, 2) *
			     pow(ComovingBoxSize/HubbleConstantNow, 3));
  ListOfParticles FullList[NUM_PARTICLE_TYPES];

  for (i = 0; i < NUM_PARTICLE_TYPES; i++) {

    if (TotalNumberOfParticles[i] > 0) {

      /* Allocate a field for these particles. */

      FullList[i].NumberOfValues = ListOfParticlesHead[i]->NumberOfValues;
      if (debug) 
	printf("NumberOfValues[%d] = %d\n", i, FullList[i].NumberOfValues);
      for (dim = 0; dim < MetaData.TopGridRank; dim++) {
	FullList[i].ParticlePosition[dim] = new 
	  float[TotalNumberOfParticles[i]];
	FullList[i].ParticleVelocity[dim] = new 
	  float[TotalNumberOfParticles[i]];
      }
      FullList[i].ParticleRadius = new float[TotalNumberOfParticles[i]];
      for (j = 0; j < FullList[i].NumberOfValues; j++)
	FullList[i].ParticleValue[j] = new float[TotalNumberOfParticles[i]];
      if (i >= 1)
	FullList[i].ParticleIndex = new int[TotalNumberOfParticles[i]];

      /* Copy grid lists into the full list. */

      ListOfParticles *Temp = ListOfParticlesHead[i];
      count = 0;
      while (Temp != NULL) {
	//printf("i=%d part=%d count=%d\n", i, Temp->NumberOfParticles, count);
	for (dim = 0; dim < MetaData.TopGridRank; dim++)
	  for (part = 0; part < Temp->NumberOfParticles; part++) {
	    FullList[i].ParticlePosition[dim][count+part] = 
	      Temp->ParticlePosition[dim][part];
	    FullList[i].ParticleVelocity[dim][count+part] = 
	      Temp->ParticleVelocity[dim][part];
	  }
	for (part = 0; part < Temp->NumberOfParticles; part++)
	  FullList[i].ParticleRadius[count+part] = Temp->ParticleRadius[part];
	for (j = 0; j < FullList[i].NumberOfValues; j++)
	  for (part = 0; part < Temp->NumberOfParticles; part++)
	    FullList[i].ParticleValue[j][count+part] = 
	      Temp->ParticleValue[j][part];
	if (i >= 1)
	  for (part = 0; part < Temp->NumberOfParticles; part++)
	    FullList[i].ParticleIndex[count+part] = Temp->ParticleIndex[part];

	/* Copy positions into kd structure. */

	int izero = SelfGravity - 1;  // always zero
	for (part = 0; part < Temp->NumberOfParticles; part++) {
	  kd->p[kdcount].r[0] = Temp->ParticlePosition[0][part];
	  kd->p[kdcount].r[1] = Temp->ParticlePosition[1][part];
	  kd->p[kdcount].r[2] = Temp->ParticlePosition[2][part];
	  for (dim = 0; dim < MetaData.TopGridRank; dim++)
	    if (kd->p[kdcount].r[dim] != 
		kd->p[kdcount].r[dim+izero]) {
	      printf("warn: %d %d %d %d %g\n", kdcount, part, i, dim,
		     kd->p[kdcount].r[dim]);
	      kd->p[kdcount].r[dim] = 0;
	    }
	  kd->p[kdcount].fMass = Temp->ParticleValue[0][part]*DensityFactor;
	  kdcount++;
	}

	/* Delete this grid's particles. */

	if (Temp->NumberOfParticles > 0) {
	  for (dim = 0; dim < MetaData.TopGridRank; dim++) {
	    delete [] Temp->ParticlePosition[dim];
	    delete [] Temp->ParticleVelocity[dim];
	  }
	  delete [] Temp->ParticleRadius;
	  if (i >= 1)
	    delete [] Temp->ParticleIndex;
	  for (j = 0; j < FullList[i].NumberOfValues; j++)
	    delete [] Temp->ParticleValue[j];
	}

	count += Temp->NumberOfParticles;
	Temp = Temp->NextList;
      }

    } // end: if (TotalNumberOfParticles[i] > 0)

  } // end: loop over particle types

  /* --------------------------------------------------------------- */
  /* Call hop. */

  fprintf(stderr, "Calling hop...\n");
  hop_main(kd);

  fprintf(stderr, "Calling regroup...\n");
  regroup_main(HopDensityThreshold);


  /* --------------------------------------------------------------- */
  /* Read the group membership and compute group properties. */

  FILE *fptr;
#ifdef PARTICLE_OUTPUT
  FILE *fparticles;
#endif /* PARTICLE_OUTPUT */

  if ((fptr = fopen("zregroup.tag", "r")) == NULL) {
    fprintf(stderr, "Error opening regroup output zregroup.hop\n");
    my_exit(EXIT_FAILURE);
  }

  int nActive, nGroups;
  fread(&nActive, 4, 1, fptr);
  fread(&nGroups, 4, 1, fptr);
  printf("nActive = %d(=%d)   nGroups = %d\n", nActive, kd->nActive, nGroups);

  /* Allocate space and read group memberships for the particles. */

  int *GroupID = new int[nActive];
  if (fread(GroupID, 4, nActive, fptr) != nActive) {
    fprintf(stderr, "Error reading GroupID file zregroup.hop\n");
    my_exit(EXIT_FAILURE);
  }
  fclose(fptr);

  /* Allocate space and read group memberships for the particles. */

  float *Density = new float[nActive];
  if ((fptr = fopen("output_hop.den", "r")) == NULL) {
    fprintf(stderr, "Error opening regroup output output_hop.den\n");
    my_exit(EXIT_FAILURE);
  }
  fread(&nActive, 4, 1, fptr);
  if (fread(Density, 4, nActive, fptr) != nActive) {
    fprintf(stderr, "Error reading GroupID file output_hop.den\n");
    my_exit(EXIT_FAILURE);
  }
  fclose(fptr);

#ifdef PARTICLE_OUTPUT
  /* Output particle file properties. */
  if ((fparticles = fopen("HopParticles.out", "w")) == NULL) {
    fprintf(stderr, "Error opening regroup output HopParticles.out\n");
    my_exit(EXIT_FAILURE);
  }

  fprintf(fparticles, "#ParticleID GroupID Density\n");
#endif /* PARTICLE_OUTPUT */

  /* Allocate and initialize group information. */

  int const NumberOfGroupProperties = 14;
  float *GroupProperties[NumberOfGroupProperties];
  for (i = 0; i < NumberOfGroupProperties; i++) {
    GroupProperties[i] = new float[nGroups];
    for (j = 0; j < nGroups; j++)
      GroupProperties[i][j] = 0;
  }

  /* Loop over particles, adding to group properties. */

  int itype, index;
#ifdef PARTICLE_OUTPUT
  int ParticleIndex, GroupIndex, DensityValue;
#endif /* PARTICLE_OUTPUT */
  float Luminosity, Mass;
  for (i = 0; i < nActive; i++)
    if ((j = GroupID[i]) >= 0) {

      /* Convert GroupID into twopart index for FullList. */

      itype = 0;  /* particle type 0 */
      index = i;
      while (index >= TotalNumberOfParticles[itype])
	index -= TotalNumberOfParticles[itype++];

#ifdef PARTICLE_OUTPUT
      ParticleIndex = FullList[itype].ParticleIndex[index];
      GroupIndex = j;
      DensityValue = Density[i];
      fprintf(fparticles, "%d     %d     %d\n",ParticleIndex, GroupIndex, DensityValue);
#endif /* PARTICLE_OUTPUT */

      /* Total mass. */

      Mass = FullList[itype].ParticleValue[0][index];
      GroupProperties[0][j] += Mass;

      /* Free-free emission (only for baryons). */

      if (itype == 0)
	GroupProperties[1][j] += (Luminosity =
				  FullList[itype].ParticleValue[1][index] *
				  FullList[itype].ParticleValue[0][index] *
			      sqrt(FullList[itype].ParticleValue[2][index]));

      /* Luminosity-weighted Temperature (only for baryons). */

      if (itype == 0)
	GroupProperties[2][j] += Luminosity*FullList[itype].ParticleValue[2][index];

      /* Number of points. */

      GroupProperties[3][j] += 1.0;

      /* Max density (and position). */

      if (Density[i] > GroupProperties[4][j]) {
	GroupProperties[4][j] = Density[i];
	GroupProperties[5][j] = FullList[itype].ParticlePosition[0][index];
	GroupProperties[6][j] = FullList[itype].ParticlePosition[1][index];
	GroupProperties[7][j] = FullList[itype].ParticlePosition[2][index];
      }

      /* center-of-mass position */

      GroupProperties[8][j] += FullList[itype].ParticlePosition[0][index]*Mass;
      GroupProperties[9][j] += FullList[itype].ParticlePosition[1][index]*Mass;
      GroupProperties[10][j] +=FullList[itype].ParticlePosition[2][index]*Mass;

      /* mass-weighted velocity. */

      GroupProperties[11][j] +=FullList[itype].ParticleVelocity[0][index]*Mass;
      GroupProperties[12][j] +=FullList[itype].ParticleVelocity[1][index]*Mass;
      GroupProperties[13][j] +=FullList[itype].ParticleVelocity[2][index]*Mass;

    }

  for (j = 0; j < nGroups; j++) {

    /* Normalize temperature. */

    if (GroupProperties[1][j] > 0)
      GroupProperties[2][j] /= GroupProperties[1][j];

    /* Normalize center-of-mass position and velocity. */

    if (GroupProperties[0][j] > 0)
      for (i = 8; i < 14; i++)
	GroupProperties[i][j] /= GroupProperties[0][j];

  }

  /* Output group properties. */

  if ((fptr = fopen("HopAnalysis.out", "w")) == NULL) {
    fprintf(stderr, "Error opening regroup output HopAnalysis.out\n");
    my_exit(EXIT_FAILURE);
  }

  fprintf(fptr, "#Group     Mass    # part    max dens      x          y         z     center-of-mass x      y           z         vx         vy       vz\n");
  for (j = 0; j < nGroups; j++) {
    fprintf(fptr, "%d      ", j);
    for (i = 0; i < 14; i++)
      if (i < 1 || i > 2)
	fprintf(fptr, " %.9g ", GroupProperties[i][j]);
    fprintf(fptr, "\n");
  }

  fclose(fptr);

#ifdef PARTICLE_OUTPUT
  fclose(fparticles);
#endif /* PARTICLE_OUTPUT */

  my_exit(EXIT_SUCCESS);
}

void my_exit(int status)
{
  CommunicationFinalize();
  exit(status);
}
