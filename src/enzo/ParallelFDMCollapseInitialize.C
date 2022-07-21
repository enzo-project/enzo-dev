
/***********************************************************************
/
/  INITIALIZE A COLLAPSE TEST
/
/  written by: Greg Bryan
/  date:       May, 1998
/  modified1:
/
/  PURPOSE:
/    Set up a number of spherical objects
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

// This routine intializes a new simulation based on the parameter file.
//
#ifdef USE_MPI
#include "mpi.h"
#ifdef USE_MPE
#include "mpe.h"
#endif /* USE_MPE */
#endif /* USE_MPI */

#include <string.h>
#include <stdio.h>
#include <math.h>
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
#include "fortran.def"
#include "error.h"
#include "message.h"

/* function prototypes */
 
void WriteListOfFloats(FILE *fptr, int N, float floats[]);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void WriteListOfInts(FILE *fptr, int N, int nums[]);
int CommunicationAllSumIntegerValues(int *Values, int Number);

static char *FDMCollapseRePsiName          = NULL;
static char *FDMCollapseImPsiName          = NULL;
static char *FDMCollapseAbsBdName          = NULL;
static int   FDMCollapseSubgridsAreStatic    = TRUE;
static int   FDMCollapseNumberOfInitialGrids = 1;
static int   FDMUseParticles               = FALSE;
static float FDMParticleMeanDensity        = FLOAT_UNDEFINED;
#define MAX_INITIAL_GRIDS 10

int ParallelFDMCollapseInitialize(FILE *fptr, FILE *Outfptr, 
			  HierarchyEntry &TopGrid, TopGridData &MetaData)
{
  const char *DensName = "Density";
  const char *TEName   = "TotalEnergy";
  const char *GEName   = "GasEnergy";
  const char *Vel1Name = "x-velocity";
  const char *Vel2Name = "y-velocity";
  const char *Vel3Name = "z-velocity";
  const char *ColourName = "SN_Colour";
  const char *ElectronName = "Electron_Density";
  const char *HIName    = "HI_Density";
  const char *HIIName   = "HII_Density";
  const char *HeIName   = "HeI_Density";
  const char *HeIIName  = "HeII_Density";
  const char *HeIIIName = "HeIII_Density";
  const char *HMName    = "HM_Density";
  const char *H2IName   = "H2I_Density";
  const char *H2IIName  = "H2II_Density";
  const char *DIName    = "DI_Density";
  const char *DIIName   = "DII_Density";
  const char *HDIName   = "HDI_Density";
  const char *MetalName = "Metal_Density";
  const char *RePsiName = "Re_Psi"; 
  const char *ImPsiName = "Im_Psi"; 
  const char *FDMDensName = "FDMDensity"; 
  const char *GravPotName = "GravPotential";

  /* declarations */

  char line[MAX_LINE_LENGTH];
  int i, j, dim, gridnum, ret, SubgridsAreStatic, region;
  HierarchyEntry *Subgrid;
  gridnum = 0;
 
  char *RealPsiName = NULL, *ImagPsiName = NULL, *AbsBdName = NULL;

  /* Set default parameters and names */
 
  int   FDMCollapseGridDimension[MAX_INITIAL_GRIDS][MAX_DIMENSION];
  int   FDMCollapseGridLevel[MAX_INITIAL_GRIDS];
  FLOAT FDMCollapseGridLeftEdge[MAX_INITIAL_GRIDS][MAX_DIMENSION];
  FLOAT FDMCollapseGridRightEdge[MAX_INITIAL_GRIDS][MAX_DIMENSION];
  for (i = 0; i < MAX_INITIAL_GRIDS; i++)
    FDMCollapseGridLevel[i] = 1;
  for (dim = 0; dim < MetaData.TopGridRank; dim++) {
    FDMCollapseGridLeftEdge[0][dim] = DomainLeftEdge[dim];
    FDMCollapseGridRightEdge[0][dim] = DomainRightEdge[dim];
    FDMCollapseGridDimension[0][dim] = MetaData.TopGridDims[dim];
  }

  FDMCollapseGridLevel[0] = 0;

  /* read parameters */
  /* Read input from file. */
 
  char *dummy = new char[MAX_LINE_LENGTH];
  dummy[0] = 0;
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
 
    ret = 0;
 
    /* Read parameters */
 
    if (sscanf(line, "FDMCollapseRePsiName = %s", dummy) == 1)
      FDMCollapseRePsiName = dummy;
    if (sscanf(line, "FDMCollapseImPsiName = %s", dummy) == 1)
      FDMCollapseImPsiName = dummy;
    if (sscanf(line, "FDMCollapseAbsBdName = %s", dummy) == 1)
      FDMCollapseAbsBdName = dummy;

    ret += sscanf(line, "FDMUseParticles = %"ISYM, &FDMUseParticles);
    ret += sscanf(line, "FDMParticleMeanDensity = %"FSYM, &FDMParticleMeanDensity);
    ret += sscanf(line, "FDMCollapseAbsorbingBoundary = %"ISYM, &FDMCollapseAbsorbingBoundary);

        /* If the dummy char space was used, then make another. */
 
    if (dummy[0] != 0) {
      dummy = new char[MAX_LINE_LENGTH];
      ret++;
    }
      
    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && strstr(line, "FDMCollapse") && line[0] != '#')
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);
  } 

   /* -------------------------------------------------------------------- */
  /* Generate the root grid and set-up the hierarchy. */
 
  HierarchyEntry *GridsList;
  GridsList = &TopGrid;
 
  /* Initialize the root grid. */
 
  RealPsiName            = FDMCollapseRePsiName;
  ImagPsiName            = FDMCollapseImPsiName;
  AbsBdName              = FDMCollapseAbsBdName;

  int TotalRefinement = nint(POW(FLOAT(RefineBy), FDMCollapseGridLevel[gridnum]));

  if (GridsList->GridData->ParallelFDMCollapseInitializeGrid(
			   RealPsiName,
			   ImagPsiName,
			   AbsBdName,
         FDMUseParticles,
         FDMParticleMeanDensity,
         FDMCollapseSubgridsAreStatic,
			   TotalRefinement) == FAIL) {
      ENZO_FAIL("Error in grid->ParallelFDMCollapseInitializeGrid.\n");
  }

  /* set up field names and units */

  int count = 0;
  DataLabel[count++] = (char*) DensName;
  DataLabel[count++] = (char*) TEName;
  DataLabel[count++] = (char*) Vel1Name;
  DataLabel[count++] = (char*) Vel2Name;
  DataLabel[count++] = (char*) Vel3Name;

  DataLabel[count++] = (char*) RePsiName;
  DataLabel[count++] = (char*) ImPsiName;
  DataLabel[count++] = (char*) FDMDensName;
  
  if(WritePotential)
	DataLabel[count++] = (char*) GravPotName;

  for (i = 0; i < count; i++)
    DataUnits[i] = NULL;

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    if (FDMCollapseRePsiName)
    fprintf(Outfptr, "FDMCollapseRePsiName          = %s\n",
	    FDMCollapseRePsiName);
    if (FDMCollapseImPsiName)
    fprintf(Outfptr, "FDMCollapseImPsiName      = %s\n",
	    FDMCollapseImPsiName);
    if (FDMCollapseAbsBdName)
    fprintf(Outfptr, "FDMCollapseAbsBdName      = %s\n",
	    FDMCollapseAbsBdName);
  }

    delete dummy;

  return SUCCESS;

}

int ParallelFDMCollapseReInitialize(HierarchyEntry *TopGrid, TopGridData &MetaData)
{
  int dim, gridnum = 0;
  char *RealPsiName = NULL, *ImagPsiName = NULL, *AbsBdName = NULL;

  if (MyProcessorNumber == ROOT_PROCESSOR)
  printf("FDMCollapse: ReInitializing grid %"ISYM"\n", gridnum);

  if (FDMCollapseNumberOfInitialGrids > 1) {
    if (FDMCollapseRePsiName)
      sprintf(RealPsiName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
	      FDMCollapseRePsiName, gridnum);
    if (FDMCollapseImPsiName)
      sprintf(ImagPsiName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
	      FDMCollapseRePsiName, gridnum);
    if (FDMCollapseAbsBdName)
      sprintf(AbsBdName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
	      FDMCollapseAbsBdName, gridnum);
  } else {
    RealPsiName            = FDMCollapseRePsiName;
    ImagPsiName            = FDMCollapseImPsiName;
    AbsBdName              = FDMCollapseAbsBdName;
  }

  /* Call grid initializer.  Use TotalRefinement = -1 to flag real read. */
 
  int TotalRefinement = -1;
  HierarchyEntry *Temp = TopGrid;
  while (Temp != NULL) {
 
    if (Temp->GridData->ParallelFDMCollapseInitializeGrid(
			   RealPsiName,
			   ImagPsiName,
			   AbsBdName,
         FDMUseParticles,
         FDMParticleMeanDensity,
         FDMCollapseSubgridsAreStatic,
			   TotalRefinement) == FAIL) {
      ENZO_FAIL("Error in grid->ParallelFDMCollapseInitializeGrid.\n");

    }
 
    Temp = Temp->NextGridThisLevel;
  }
 
  return SUCCESS;
}

