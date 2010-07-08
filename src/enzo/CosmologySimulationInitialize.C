/***********************************************************************
/
/  INITIALIZE A COSMOLOGY SIMULATION
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified:   Robert Harkness
/  date:       September 2004
/
/  PURPOSE:  Initialize for cosmology simulations.  Reads in a number
/      of initial grids.  If more than one, then they are all numbered.
/      We currently assume that all are subgrids except the first.
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
// This routine intializes a new simulation based on the parameter file.
//
 
#ifdef USE_MPI
#include "mpi.h"
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
#include "CosmologyParameters.h"
#include "fortran.def"
#include "CommunicationUtilities.h"
 
// Function prototypes
 
void WriteListOfFloats(FILE *fptr, int N, float floats[]);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void WriteListOfInts(FILE *fptr, int N, int nums[]);
void PrintMemoryUsage(char *str);
int InitializeRateData(FLOAT Time);


// Cosmology Parameters (that need to be shared)
 
static float CosmologySimulationOmegaBaryonNow       = 1.0;  // standard
static float CosmologySimulationOmegaCDMNow          = 0.0;  // no dark matter
static float CosmologySimulationInitialTemperature   = FLOAT_UNDEFINED;
 
static char *CosmologySimulationDensityName          = NULL;
static char *CosmologySimulationTotalEnergyName      = NULL;
static char *CosmologySimulationGasEnergyName        = NULL;
static char *CosmologySimulationParticlePositionName = NULL;
static char *CosmologySimulationParticleVelocityName = NULL;
static char *CosmologySimulationParticleMassName     = NULL;
static char *CosmologySimulationParticleTypeName     = NULL;
static char *CosmologySimulationVelocityNames[MAX_DIMENSION];
static char *CosmologySimulationParticleVelocityNames[MAX_DIMENSION];

static int   CosmologySimulationSubgridsAreStatic    = TRUE;
 
static float CosmologySimulationInitialFractionHII   = 1.2e-5;
static float CosmologySimulationInitialFractionHeII  = 1.0e-14;
static float CosmologySimulationInitialFractionHeIII = 1.0e-17;
static float CosmologySimulationInitialFractionHM    = 2.0e-9;
static float CosmologySimulationInitialFractionH2I   = 2.0e-20;
static float CosmologySimulationInitialFractionH2II  = 3.0e-14;
static int   CosmologySimulationUseMetallicityField  = FALSE;
 
static int CosmologySimulationManuallySetParticleMassRatio = FALSE;
static float CosmologySimulationManualParticleMassRatio = 1.0;

static int   CosmologySimulationCalculatePositions   = FALSE; 

#ifdef TRANSFER
static float RadHydroInitialRadiationEnergy = 1.0e-32;
#endif

#define MAX_INITIAL_GRIDS 10
 
 
 
 
int CosmologySimulationInitialize(FILE *fptr, FILE *Outfptr,
			       HierarchyEntry &TopGrid, TopGridData &MetaData)
{
  char *DensName = "Density";
  char *TEName   = "TotalEnergy";
  char *GEName   = "GasEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  char *ElectronName = "Electron_Density";
  char *HIName    = "HI_Density";
  char *HIIName   = "HII_Density";
  char *HeIName   = "HeI_Density";
  char *HeIIName  = "HeII_Density";
  char *HeIIIName = "HeIII_Density";
  char *HMName    = "HM_Density";
  char *H2IName   = "H2I_Density";
  char *H2IIName  = "H2II_Density";
  char *DIName    = "DI_Density";
  char *DIIName   = "DII_Density";
  char *HDIName   = "HDI_Density";
  char *MetalName = "Metal_Density";
  char *GPotName  = "Grav_Potential";
  char *ForbidName  = "ForbiddenRefinement";
  char *MachName   = "Mach";
  char *CRName     = "CR_Density";
  char *PSTempName = "PreShock_Temperature";
  char *PSDenName  = "PreShock_Density";
  char *ExtraNames[2] = {"Z_Field1", "Z_Field2"};

#ifdef TRANSFER
  char *RadName    = "Grey_Radiation_Energy";
#endif
#ifdef EMISSIVITY
  char *EtaName    = "Emissivity";
#endif
 

  char line[MAX_LINE_LENGTH];
  int i, j, dim, gridnum, ret, SubgridsAreStatic, region;
  HierarchyEntry *Subgrid;
 
  char *DensityName = NULL, *TotalEnergyName = NULL, *GasEnergyName = NULL,
       *ParticlePositionName = NULL, *ParticleVelocityName = NULL,
       *ParticleMassName = NULL, *VelocityNames[MAX_DIMENSION],
       *ParticleTypeName = NULL, *ParticleVelocityNames[MAX_DIMENSION];
 
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    VelocityNames[dim] = NULL;
 
  // Set default parameters: parameters, names and subgrid info
 
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    VelocityNames[dim] = NULL;
    ParticleVelocityNames[dim] = NULL;
    CosmologySimulationParticleVelocityNames[dim] = NULL;
    CosmologySimulationVelocityNames[dim] = NULL;
  }
 
  int   CosmologySimulationGridDimension[MAX_INITIAL_GRIDS][MAX_DIMENSION];
  int   CosmologySimulationGridLevel[MAX_INITIAL_GRIDS];
  FLOAT CosmologySimulationGridLeftEdge[MAX_INITIAL_GRIDS][MAX_DIMENSION];
  FLOAT CosmologySimulationGridRightEdge[MAX_INITIAL_GRIDS][MAX_DIMENSION];
  for (i = 0; i < MAX_INITIAL_GRIDS; i++)
    CosmologySimulationGridLevel[i] = 1;
  for (dim = 0; dim < MetaData.TopGridRank; dim++) {
    CosmologySimulationGridLeftEdge[0][dim] = DomainLeftEdge[dim];
    CosmologySimulationGridRightEdge[0][dim] = DomainRightEdge[dim];
    CosmologySimulationGridDimension[0][dim] = MetaData.TopGridDims[dim];
  }
  CosmologySimulationGridLevel[0] = 0;
 
  if (!ComovingCoordinates) {
    ENZO_FAIL("ComovingCoordinates must be TRUE!\n");
  }
 
  if (DualEnergyFormalism == FALSE && HydroMethod != Zeus_Hydro)
    fprintf(stderr, "CosmologySimulation: DualEnergyFormalism is off!\n");
  if (!SelfGravity)
    fprintf(stderr, "CosmologySimulation: gravity is off!?!\n");
 
  // Read input from file
 
  char *dummy = new char[MAX_LINE_LENGTH];
  dummy[0] = 0;
  int comment_count = 0;
 
  while ((fgets(line, MAX_LINE_LENGTH, fptr) != NULL) 
      && (comment_count < 2)) {
 
    ret = 0;
 
    // Read parameters
 
    ret += sscanf(line, "CosmologySimulationOmegaBaryonNow = %"FSYM,
		  &CosmologySimulationOmegaBaryonNow);
    ret += sscanf(line, "CosmologySimulationOmegaCDMNow = %"FSYM,
		  &CosmologySimulationOmegaCDMNow);
    ret += sscanf(line, "CosmologySimulationInitialTemperature = %"FSYM,
		  &CosmologySimulationInitialTemperature);
 
    if (sscanf(line, "CosmologySimulationDensityName = %s", dummy) == 1)
      CosmologySimulationDensityName = dummy;
    if (sscanf(line, "CosmologySimulationTotalEnergyName = %s", dummy) == 1)
      CosmologySimulationTotalEnergyName = dummy;
    if (sscanf(line, "CosmologySimulationGasEnergyName = %s", dummy) == 1)
      CosmologySimulationGasEnergyName = dummy;
    if (sscanf(line, "CosmologySimulationVelocity1Name = %s", dummy) == 1)
      CosmologySimulationVelocityNames[0] = dummy;
    if (sscanf(line, "CosmologySimulationVelocity2Name = %s", dummy) == 1)
      CosmologySimulationVelocityNames[1] = dummy;
    if (sscanf(line, "CosmologySimulationVelocity3Name = %s", dummy) == 1)
      CosmologySimulationVelocityNames[2] = dummy;
    if (sscanf(line, "CosmologySimulationParticlePositionName = %s", dummy)
	== 1) CosmologySimulationParticlePositionName = dummy;
    if (sscanf(line, "CosmologySimulationParticleVelocityName = %s", dummy)
	== 1) CosmologySimulationParticleVelocityName = dummy;
    if (sscanf(line, "CosmologySimulationParticleMassName = %s", dummy) == 1)
      CosmologySimulationParticleMassName = dummy;
    if (sscanf(line, "CosmologySimulationParticleTypeName = %s", dummy) == 1)
      CosmologySimulationParticleTypeName = dummy;
    if (sscanf(line, "CosmologySimulationParticleVelocity1Name = %s", dummy) == 1)
      CosmologySimulationParticleVelocityNames[0] = dummy;
    if (sscanf(line, "CosmologySimulationParticleVelocity2Name = %s", dummy) == 1)
      CosmologySimulationParticleVelocityNames[1] = dummy;
    if (sscanf(line, "CosmologySimulationParticleVelocity3Name = %s", dummy) == 1)
      CosmologySimulationParticleVelocityNames[2] = dummy;    

    ret += sscanf(line, "CosmologySimulationNumberOfInitialGrids = %"ISYM,
		  &CosmologySimulationNumberOfInitialGrids);
    ret += sscanf(line, "CosmologySimulationSubgridsAreStatic = %"ISYM,
		  &CosmologySimulationSubgridsAreStatic);
 
    if (sscanf(line, "CosmologySimulationGridLeftEdge[%"ISYM"]", &gridnum) > 0)
      ret += sscanf(line, "CosmologySimulationGridLeftEdge[%"ISYM"] = %"PSYM" %"PSYM" %"PSYM,
		    &gridnum, &CosmologySimulationGridLeftEdge[gridnum][0],
		    &CosmologySimulationGridLeftEdge[gridnum][1],
		    &CosmologySimulationGridLeftEdge[gridnum][2]);
    if (sscanf(line, "CosmologySimulationGridRightEdge[%"ISYM"]", &gridnum) > 0)
      ret += sscanf(line, "CosmologySimulationGridRightEdge[%"ISYM"] = %"PSYM" %"PSYM" %"PSYM,
		    &gridnum, &CosmologySimulationGridRightEdge[gridnum][0],
		    &CosmologySimulationGridRightEdge[gridnum][1],
		    &CosmologySimulationGridRightEdge[gridnum][2]);
    if (sscanf(line, "CosmologySimulationGridDimension[%"ISYM"]", &gridnum) > 0)
      ret += sscanf(line, "CosmologySimulationGridDimension[%"ISYM"] = %"ISYM" %"ISYM" %"ISYM,
		    &gridnum, &CosmologySimulationGridDimension[gridnum][0],
		    &CosmologySimulationGridDimension[gridnum][1],
		    &CosmologySimulationGridDimension[gridnum][2]);
    if (sscanf(line, "CosmologySimulationGridLevel[%"ISYM"]", &gridnum) > 0)
      ret += sscanf(line, "CosmologySimulationGridLevel[%"ISYM"] = %"ISYM,
		    &gridnum, &CosmologySimulationGridLevel[gridnum]);
 
    ret += sscanf(line, "CosmologySimulationInitialFractionHII = %"FSYM,
		  &CosmologySimulationInitialFractionHII);
    ret += sscanf(line, "CosmologySimulationInitialFractionHeII = %"FSYM,
		  &CosmologySimulationInitialFractionHeII);
    ret += sscanf(line, "CosmologySimulationInitialFractionHeIII = %"FSYM,
		  &CosmologySimulationInitialFractionHeIII);
    ret += sscanf(line, "CosmologySimulationInitialFractionHM = %"FSYM,
		  &CosmologySimulationInitialFractionHM);
    ret += sscanf(line, "CosmologySimulationInitialFractionH2I = %"FSYM,
		  &CosmologySimulationInitialFractionH2I);
    ret += sscanf(line, "CosmologySimulationInitialFractionH2II = %"FSYM,
		  &CosmologySimulationInitialFractionH2II);
    ret += sscanf(line, "CosmologySimulationUseMetallicityField = %"ISYM,
		  &CosmologySimulationUseMetallicityField);
 
    ret += sscanf(line, "CosmologySimulationManuallySetParticleMassRatio = %"ISYM,
		  &CosmologySimulationManuallySetParticleMassRatio);
    ret += sscanf(line, "CosmologySimulationManualParticleMassRatio = %"FSYM,
		  &CosmologySimulationManualParticleMassRatio);

    ret += sscanf(line, "CosmologySimulationCalculatePositions = %"ISYM,
		  &CosmologySimulationCalculatePositions);

    // If the dummy char space was used, then make another
 
    if (dummy[0] != 0) {
      dummy = new char[MAX_LINE_LENGTH];
      ret++;
    }
 
    // If the line is suspicious, issue a warning

    if (strstr(line, "\"\"\"")) comment_count++;
 
    if (ret == 0 && strstr(line, "=") && strstr(line, "CosmologySimulation") &&
	line[0] != '#' && MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);
 
  }

#ifdef TRANSFER
  // if using FLD-based solver for radiation, initialize the relevant fields
  if (RadiativeTransferFLD > 1) {
    // Read RadHydro input from secondary input file
    if (MetaData.RadHydroParameterFname != NULL) {
      FILE *RHfptr;
      if ((RHfptr = fopen(MetaData.RadHydroParameterFname, "r")) != NULL) {
	while (fgets(line, MAX_LINE_LENGTH, RHfptr) != NULL) {
	  ret = 0;
	  // read relevant problem parameters
	  ret += sscanf(line, "RadHydroRadiationEnergy = %"FSYM, 
			&RadHydroInitialRadiationEnergy);
	} // end input from parameter file
	fclose(RHfptr);
      }
    }
    
    // set up CoolData object if not already set up
    if (CoolData.ceHI == NULL) 
      if (InitializeRateData(MetaData.Time) == FAIL) {
	fprintf(stderr,"Error in InitializeRateData.\n");
	return FAIL;
      }
  }

#endif
 
  // More error checking
 
  if (CosmologySimulationDensityName == NULL &&
      (CosmologySimulationParticlePositionName == NULL &&
       !CosmologySimulationCalculatePositions)) {
    ENZO_FAIL("Missing initial data.\n");
  }
 
  if (CosmologySimulationDensityName != NULL && CellFlaggingMethod[0] != 2)
    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, "CosmologySimulation: check CellFlaggingMethod.\n");
 
  if (CosmologySimulationDensityName == NULL && CellFlaggingMethod[0] != 4)
    if (MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, "CosmologySimulation: check CellFlaggingMethod.\n");
 
  if (CosmologySimulationNumberOfInitialGrids > MAX_INITIAL_GRIDS) {
    ENZO_FAIL("Too many InitialGrids! increase MAX_INITIAL_GRIDS\n");
  }
 
  if (CosmologySimulationDensityName == NULL &&
      MultiSpecies+RadiativeCooling > 0) {
    fprintf(stderr, "warning: no density field; setting MultiSpecies/RadiativeCooling = 0\n");
    MultiSpecies = RadiativeCooling = 0;
  }

  if (CosmologySimulationParticleVelocityNames[0] != NULL &&
      !CosmologySimulationCalculatePositions) {
    ENZO_FAIL("CosmologySimulation: 1-component files only valid for use with "
	    "CosmologySimualtionCalculatePositions.\n");
  }
 
  // If temperature is left unset, set it assuming that T=550 K at z=200
 
  if (CosmologySimulationInitialTemperature == FLOAT_UNDEFINED)
    CosmologySimulationInitialTemperature = 550.0 *
      POW((1.0 + InitialRedshift)/(1.0 + 200), 2);
 
  // Generate the grids and set-up the hierarchy
 
  HierarchyEntry *GridsList[MAX_INITIAL_GRIDS];
  GridsList[0] = &TopGrid;
  for (gridnum = 1; gridnum < CosmologySimulationNumberOfInitialGrids;
       gridnum++) {
 
    // Create a spot in the hierarchy
 
    Subgrid    = new HierarchyEntry;
 
    // Find where to put this new grid
 
    int ParentGrid = INT_UNDEFINED;
    for (i = 0; i < gridnum; i++)
      if (CosmologySimulationGridLevel[i] ==
	  CosmologySimulationGridLevel[gridnum]-1)
	for (dim = 0; dim < MetaData.TopGridRank; dim++) {
	  if (CosmologySimulationGridLeftEdge[gridnum][dim] <
	      CosmologySimulationGridLeftEdge[i][dim]       ||
	      CosmologySimulationGridRightEdge[gridnum][dim] >
	      CosmologySimulationGridRightEdge[i][dim]       )
	    break;
	  ParentGrid = i;
	}
 
    if (ParentGrid == INT_UNDEFINED) {
      ENZO_VFAIL("Grid %"ISYM" has no valid parent.\n", gridnum)
    }
 
    // Insert this grid at the appropriate position in the subgrid chain
 
    GridsList[gridnum] = Subgrid;
    Subgrid->NextGridNextLevel = NULL;
    Subgrid->NextGridThisLevel = GridsList[ParentGrid]->NextGridNextLevel;
    Subgrid->ParentGrid        = GridsList[ParentGrid];
    GridsList[ParentGrid]->NextGridNextLevel = Subgrid;
 
    // Error check for consistency and add ghost zones to dimension
 
    for (dim = 0; dim < MetaData.TopGridRank; dim++) {
      FLOAT SubgridCellSize = (DomainRightEdge[dim] - DomainLeftEdge[dim])/
	FLOAT(MetaData.TopGridDims[dim]*
	      POW(FLOAT(RefineBy), CosmologySimulationGridLevel[gridnum]));
      if (debug1) {
	printf("%"GSYM"\n", SubgridCellSize);
	printf("%"GSYM"\n", POW(FLOAT(RefineBy), CosmologySimulationGridLevel[gridnum]));
	printf("%"ISYM" %"ISYM"\n", MetaData.TopGridDims[dim], dim);
	printf("%"GSYM" %"GSYM"\n", CosmologySimulationGridRightEdge[gridnum][dim],
	       CosmologySimulationGridLeftEdge[gridnum][dim]);
	printf("%"ISYM"\n", nint((CosmologySimulationGridRightEdge[gridnum][dim] -
				  CosmologySimulationGridLeftEdge[gridnum][dim]   )
				 /SubgridCellSize));
      } // ENDIF debug1
 
      // Check if declared size matches left/right edges
 
      if (nint((CosmologySimulationGridRightEdge[gridnum][dim] -
		CosmologySimulationGridLeftEdge[gridnum][dim]   )
	       /SubgridCellSize) !=
	  CosmologySimulationGridDimension[gridnum][dim]) {
	fprintf(stderr, "Subgrid inconsistency: grid %"ISYM", dim %"ISYM"\n",
		gridnum, dim);
	fprintf(stderr, " subgrid: %"GOUTSYM" -> %"GOUTSYM", CellSize = %"GOUTSYM"\n",
	      CosmologySimulationGridLeftEdge[gridnum][dim],
	      CosmologySimulationGridRightEdge[gridnum][dim], SubgridCellSize);
	ENZO_FAIL("Subgrid Inconsistency!");
      }
 
      // Check if left/right edge fall on Parent cell boundary
 
      if (nint((CosmologySimulationGridLeftEdge[gridnum][dim] -
		CosmologySimulationGridLeftEdge[ParentGrid][dim])/
	       SubgridCellSize) % RefineBy != 0 ||
	  nint((CosmologySimulationGridRightEdge[gridnum][dim] -
		CosmologySimulationGridLeftEdge[ParentGrid][dim])/
	       SubgridCellSize) % RefineBy != 0 ) {
	fprintf(stderr, "Subgrid inconsistency: grid %"ISYM", dim %"ISYM"\n",
		gridnum, dim);
	fprintf(stderr, "left or right edges are not on parent cell edge.\n");
	ENZO_FAIL("Subgrid Inconsistency!");
      }
 
      // Add ghost zones
 
      CosmologySimulationGridDimension[gridnum][dim] += 2*DEFAULT_GHOST_ZONES;
    }
 
    // Create a new subgrid and initialize it
 
    Subgrid->GridData = new grid;
    Subgrid->GridData->InheritProperties(Subgrid->ParentGrid->GridData);
    Subgrid->GridData->PrepareGrid(MetaData.TopGridRank,
				   CosmologySimulationGridDimension[gridnum],
				   CosmologySimulationGridLeftEdge[gridnum],
				   CosmologySimulationGridRightEdge[gridnum],
				   0);
 
    // If subgrids are static, convert to static regions
 
    if (CosmologySimulationSubgridsAreStatic == TRUE) {
      for (region = 0; region < MAX_STATIC_REGIONS; region++)
	if (StaticRefineRegionLevel[region] == INT_UNDEFINED) {
	  StaticRefineRegionLevel[region] =
	    CosmologySimulationGridLevel[gridnum] - 1;
	  for (dim = 0; dim < MetaData.TopGridRank; dim++) {
	    StaticRefineRegionLeftEdge[region][dim] =
	      CosmologySimulationGridLeftEdge[gridnum][dim];
	    StaticRefineRegionRightEdge[region][dim] =
	      CosmologySimulationGridRightEdge[gridnum][dim];
	  }
	  for (dim = MetaData.TopGridRank; dim < MAX_DIMENSION; dim++) {
	    StaticRefineRegionLeftEdge[region][dim] = DomainLeftEdge[dim];
	    StaticRefineRegionRightEdge[region][dim] = DomainRightEdge[dim];
	  }
	  break;
	}
      if (region == MAX_STATIC_REGIONS) {
	ENZO_FAIL("Out of static refine regions\n");
      }
    }
 
    // Remove ghost zones from dim
 
    for (dim = 0; dim < MetaData.TopGridRank; dim++)
      CosmologySimulationGridDimension[gridnum][dim] -=
	2*DEFAULT_GHOST_ZONES;
 
  } // end: loop over gridnums
 
 
  // Initialize the previously-generated grids
 
  for (gridnum = 0; gridnum < CosmologySimulationNumberOfInitialGrids;
       gridnum++) {
 
    // If there is more than one grid, add the grid number to the name
 
    if (CosmologySimulationNumberOfInitialGrids > 1) {
 
      if (MyProcessorNumber == ROOT_PROCESSOR)
	printf("CosmologySimulation: Initializing grid %"ISYM"\n", gridnum);
 
      if (CosmologySimulationDensityName)
	sprintf(DensityName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
		CosmologySimulationDensityName, gridnum);
      if (CosmologySimulationTotalEnergyName)
	sprintf(TotalEnergyName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
		CosmologySimulationTotalEnergyName, gridnum);
      if (CosmologySimulationGasEnergyName)
	sprintf(GasEnergyName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
		CosmologySimulationGasEnergyName, gridnum);
      for (dim = 0; dim < MetaData.TopGridRank; dim++)
	if (CosmologySimulationVelocityNames[dim])
	  sprintf(VelocityNames[dim] = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
		  CosmologySimulationVelocityNames[dim], gridnum);
      if (CosmologySimulationParticlePositionName)
	sprintf(ParticlePositionName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
		CosmologySimulationParticlePositionName, gridnum);
      if (CosmologySimulationParticleVelocityName)
	sprintf(ParticleVelocityName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
		CosmologySimulationParticleVelocityName, gridnum);
      if (CosmologySimulationParticleMassName)
	sprintf(ParticleMassName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
		CosmologySimulationParticleMassName, gridnum);
      if (CosmologySimulationParticleTypeName)
        sprintf(ParticleTypeName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
                CosmologySimulationParticleTypeName, gridnum);
      for (dim = 0; dim < MetaData.TopGridRank; dim++)
	if (CosmologySimulationParticleVelocityNames[dim])
	  sprintf(ParticleVelocityNames[dim] = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
		  CosmologySimulationParticleVelocityNames[dim], gridnum);
 
    } else {
      DensityName            = CosmologySimulationDensityName;
      TotalEnergyName        = CosmologySimulationTotalEnergyName;
      GasEnergyName          = CosmologySimulationGasEnergyName;
      for (dim = 0; dim < MetaData.TopGridRank; dim++)
	VelocityNames[dim]   = CosmologySimulationVelocityNames[dim];
      ParticlePositionName   = CosmologySimulationParticlePositionName;
      ParticleVelocityName   = CosmologySimulationParticleVelocityName;
      ParticleMassName       = CosmologySimulationParticleMassName;
      ParticleTypeName       = CosmologySimulationParticleTypeName;
      for (dim = 0; dim < MetaData.TopGridRank; dim++)
	ParticleVelocityNames[dim]   = CosmologySimulationParticleVelocityNames[dim];
    }
 
    // If there is a subgrid, use CosmologySimulationSubgridsAreStatic,
    // otherwise just set to false
 
    SubgridsAreStatic = (GridsList[gridnum]->NextGridNextLevel == NULL) ?
      FALSE : CosmologySimulationSubgridsAreStatic;
 
    // Initialize the grid by reading in data
 
    int TotalRefinement = nint(POW(FLOAT(RefineBy),
				   CosmologySimulationGridLevel[gridnum]));
    if (GridsList[gridnum]->GridData->CosmologySimulationInitializeGrid(
			     gridnum,
			     CosmologySimulationOmegaBaryonNow,
			       CosmologySimulationOmegaCDMNow,
			       CosmologySimulationInitialTemperature,
			     DensityName, TotalEnergyName,
			       GasEnergyName, VelocityNames,
			       ParticlePositionName, ParticleVelocityName,
			       ParticleMassName, ParticleTypeName,
			     ParticleVelocityNames,
			     SubgridsAreStatic, TotalRefinement,
			     CosmologySimulationInitialFractionHII,
			     CosmologySimulationInitialFractionHeII,
			     CosmologySimulationInitialFractionHeIII,
			     CosmologySimulationInitialFractionHM,
			     CosmologySimulationInitialFractionH2I,
			     CosmologySimulationInitialFractionH2II,
#ifdef TRANSFER
			     RadHydroInitialRadiationEnergy,
#endif
			     CosmologySimulationUseMetallicityField,
			     MetaData.NumberOfParticles,
			     CosmologySimulationManuallySetParticleMassRatio,
			     CosmologySimulationManualParticleMassRatio,
			     CosmologySimulationCalculatePositions
						       ) == FAIL) {
      ENZO_FAIL("Error in grid->CosmologySimulationInitializeGrid.\n");
    }
 
    // Set boundary conditions if necessary
 
  } // end loop over initial grids
 
 
  /* Convert minimum initial overdensity for refinement to mass
     (unless MinimumMass itself was actually set).
     Note: multiply MinimumMassForRefinement by the OmegaBaryonNow since the
     routine that uses this parameter only counts baryonic mass. */
 
  for (i = 0; i < MAX_FLAGGING_METHODS; i++) {
    if (MinimumMassForRefinement[i] == FLOAT_UNDEFINED) {
 
      MinimumMassForRefinement[i] = CosmologySimulationOmegaBaryonNow/
	                            OmegaMatterNow;
      if (CellFlaggingMethod[i] == 4)
	MinimumMassForRefinement[i] = CosmologySimulationOmegaCDMNow/
	                              OmegaMatterNow;
 
      MinimumMassForRefinement[i] *= MinimumOverDensityForRefinement[i];
      for (dim = 0; dim < MetaData.TopGridRank; dim++)
	MinimumMassForRefinement[i] *=
	  (DomainRightEdge[dim]-DomainLeftEdge[dim])/
	  float(MetaData.TopGridDims[dim]);
    }
  }
 
  // Set up field names and units
 
  i = 0;
  DataLabel[i++] = DensName;
  DataLabel[i++] = TEName;
  if (DualEnergyFormalism)
    DataLabel[i++] = GEName;
  DataLabel[i++] = Vel1Name;
  DataLabel[i++] = Vel2Name;
  DataLabel[i++] = Vel3Name;
#ifdef TRANSFER
  if (RadiativeTransferFLD > 1)
    DataLabel[i++] = RadName;
#endif
  if (MultiSpecies) {
    DataLabel[i++] = ElectronName;
    DataLabel[i++] = HIName;
    DataLabel[i++] = HIIName;
    DataLabel[i++] = HeIName;
    DataLabel[i++] = HeIIName;
    DataLabel[i++] = HeIIIName;
    if (MultiSpecies > 1) {
      DataLabel[i++] = HMName;
      DataLabel[i++] = H2IName;
      DataLabel[i++] = H2IIName;
    }
    if (MultiSpecies > 2) {
      DataLabel[i++] = DIName;
      DataLabel[i++] = DIIName;
      DataLabel[i++] = HDIName;
    }
  }
  if (CosmologySimulationUseMetallicityField) {
    DataLabel[i++] = MetalName;
    if(MultiMetals){
      DataLabel[i++] = ExtraNames[0];
      DataLabel[i++] = ExtraNames[1];
    }
  }
  if(STARMAKE_METHOD(COLORED_POP3_STAR)){
    DataLabel[i++] = ForbidName;
  }
 
  if (WritePotential)
    DataLabel[i++] = GPotName;  

  if (CRModel) {
    DataLabel[i++] = MachName;
    if(StorePreShockFields){
      DataLabel[i++] = PSTempName;
      DataLabel[i++] = PSDenName;
    }
    DataLabel[i++] = CRName;
  } 

#ifdef EMISSIVITY
  if (StarMakerEmissivityField > 0)
    DataLabel[i++] = EtaName;
#endif
 
  for (j = 0; j < i; j++)
    DataUnits[j] = NULL;
 
  // Write parameters to parameter output file
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "CosmologySimulationOmegaBaryonNow       = %"FSYM"\n",
	    CosmologySimulationOmegaBaryonNow);
    fprintf(Outfptr, "CosmologySimulationOmegaCDMNow          = %"FSYM"\n",
	    CosmologySimulationOmegaCDMNow);
    fprintf(Outfptr, "CosmologySimulationInitialTemperature   = %"FSYM"\n\n",
	    CosmologySimulationInitialTemperature);
 
    fprintf(Outfptr, "CosmologySimulationDensityName          = %s\n",
	    CosmologySimulationDensityName);
    if (CosmologySimulationTotalEnergyName)
    fprintf(Outfptr, "CosmologySimulationTotalEnergyName      = %s\n",
	    CosmologySimulationTotalEnergyName);
    if (CosmologySimulationGasEnergyName)
    fprintf(Outfptr, "CosmologySimulationGasEnergyName        = %s\n",
	    CosmologySimulationGasEnergyName);
    fprintf(Outfptr, "CosmologySimulationVelocity1Name        = %s\n",
	    CosmologySimulationVelocityNames[0]);
    fprintf(Outfptr, "CosmologySimulationVelocity2Name        = %s\n",
	    CosmologySimulationVelocityNames[1]);
    fprintf(Outfptr, "CosmologySimulationVelocity3Name        = %s\n",
	    CosmologySimulationVelocityNames[2]);
    if (CosmologySimulationParticlePositionName)
    fprintf(Outfptr, "CosmologySimulationParticlePositionName = %s\n",
	    CosmologySimulationParticlePositionName);
    if (CosmologySimulationParticleVelocityName)
    fprintf(Outfptr, "CosmologySimulationParticleVelocityName = %s\n",
	    CosmologySimulationParticleVelocityName);
    if (CosmologySimulationParticleMassName)
    fprintf(Outfptr, "CosmologySimulationParticleMassName     = %s\n\n",
	    CosmologySimulationParticleMassName);
    if (sscanf(line, "CosmologySimulationParticleVelocity1Name = %s", dummy) == 1)
      CosmologySimulationParticleVelocityNames[0] = dummy;
    if (sscanf(line, "CosmologySimulationParticleVelocity2Name = %s", dummy) == 1)
      CosmologySimulationParticleVelocityNames[1] = dummy;
    if (sscanf(line, "CosmologySimulationParticleVelocity3Name = %s", dummy) == 1)
      CosmologySimulationParticleVelocityNames[2] = dummy;    
    
    fprintf(Outfptr, "CosmologySimulationNumberOfInitialGrids = %"ISYM"\n",
	    CosmologySimulationNumberOfInitialGrids);
    fprintf(Outfptr, "CosmologySimulationSubgridsAreStatic    = %"ISYM"\n",
	    CosmologySimulationSubgridsAreStatic);
    fprintf(Outfptr, "CosmologySimulationCalculatePositions   = %"ISYM"\n",
	    CosmologySimulationCalculatePositions);

    for (gridnum = 1; gridnum < CosmologySimulationNumberOfInitialGrids;
	 gridnum++) {
      fprintf(Outfptr, "CosmologySimulationGridLeftEdge[%"ISYM"]     = ", gridnum);
      WriteListOfFloats(Outfptr, MetaData.TopGridRank,
			CosmologySimulationGridLeftEdge[gridnum]);
      fprintf(Outfptr, "CosmologySimulationGridRightEdge[%"ISYM"]    = ", gridnum);
      WriteListOfFloats(Outfptr, MetaData.TopGridRank,
			CosmologySimulationGridRightEdge[gridnum]);
      fprintf(Outfptr, "CosmologySimulationGridDimension[%"ISYM"]    = ", gridnum);
      WriteListOfInts(Outfptr, MetaData.TopGridRank,
		      CosmologySimulationGridDimension[gridnum]);
    }
    fprintf(Outfptr, "\n");
 
    fprintf(Outfptr, "CosmologySimulationInitialFractionHII   = %"GSYM"\n",
	    CosmologySimulationInitialFractionHII);
    fprintf(Outfptr, "CosmologySimulationInitialFractionHeII  = %"GSYM"\n",
	    CosmologySimulationInitialFractionHeII);
    fprintf(Outfptr, "CosmologySimulationInitialFractionHeIII = %"GSYM"\n",
	    CosmologySimulationInitialFractionHeIII);
    fprintf(Outfptr, "CosmologySimulationInitialFractionHM    = %"GSYM"\n",
	    CosmologySimulationInitialFractionHM);
    fprintf(Outfptr, "CosmologySimulationInitialFractionH2I   = %"GSYM"\n",
	    CosmologySimulationInitialFractionH2I);
    fprintf(Outfptr, "CosmologySimulationInitialFractionH2II  = %"GSYM"\n",
	    CosmologySimulationInitialFractionH2II);
    fprintf(Outfptr, "CosmologySimulationUseMetallicityField  = %"ISYM"\n\n",
	    CosmologySimulationUseMetallicityField);
  }
 
  delete [] dummy;
 
  return SUCCESS;
}
 
 
void RecursivelySetParticleCount(HierarchyEntry *GridPoint, PINT *Count);
 
// Re-call the initializer on level zero grids.
// Used in case of ParallelRootGridIO.
 
int CosmologySimulationReInitialize(HierarchyEntry *TopGrid,
				    TopGridData &MetaData)
{
 
 
  int dim, gridnum = 0;
 
  char *DensityName = NULL, *TotalEnergyName = NULL, *GasEnergyName = NULL,
    *ParticlePositionName = NULL, *ParticleVelocityName = NULL,
    *ParticleMassName = NULL, *VelocityNames[MAX_DIMENSION],
    *ParticleTypeName = NULL, *ParticleVelocityNames[MAX_DIMENSION];
  
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    ParticleVelocityNames[dim] = NULL;
    VelocityNames[dim] = NULL;
  }
 
  if (MyProcessorNumber == ROOT_PROCESSOR)
    printf("CosmologySimulation: ReInitializing grid %"ISYM"\n", gridnum);
 
  // If there is more than one grid, add the grid number to the name
 
  if (CosmologySimulationNumberOfInitialGrids > 1) {
 
    if (CosmologySimulationDensityName)
      sprintf(DensityName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
	      CosmologySimulationDensityName, gridnum);
    if (CosmologySimulationTotalEnergyName)
      sprintf(TotalEnergyName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
	      CosmologySimulationTotalEnergyName, gridnum);
    if (CosmologySimulationGasEnergyName)
      sprintf(GasEnergyName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
	      CosmologySimulationGasEnergyName, gridnum);
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      if (CosmologySimulationVelocityNames[dim])
	sprintf(VelocityNames[dim] = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
		CosmologySimulationVelocityNames[dim], gridnum);
    if (CosmologySimulationParticlePositionName)
      sprintf(ParticlePositionName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
	      CosmologySimulationParticlePositionName, gridnum);
    if (CosmologySimulationParticleVelocityName)
      sprintf(ParticleVelocityName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
	      CosmologySimulationParticleVelocityName, gridnum);
    if (CosmologySimulationParticleMassName)
      sprintf(ParticleMassName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
	      CosmologySimulationParticleMassName, gridnum);
    if (CosmologySimulationParticleTypeName)
      sprintf(ParticleTypeName = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
              CosmologySimulationParticleTypeName, gridnum);
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      if (CosmologySimulationParticleVelocityNames[dim])
	sprintf(ParticleVelocityNames[dim] = new char[MAX_LINE_LENGTH], "%s.%1"ISYM,
		CosmologySimulationParticleVelocityNames[dim], gridnum);
 
  } else {
    DensityName            = CosmologySimulationDensityName;
    TotalEnergyName        = CosmologySimulationTotalEnergyName;
    GasEnergyName          = CosmologySimulationGasEnergyName;
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      VelocityNames[dim]   = CosmologySimulationVelocityNames[dim];
    ParticlePositionName   = CosmologySimulationParticlePositionName;
    ParticleVelocityName   = CosmologySimulationParticleVelocityName;
    ParticleMassName       = CosmologySimulationParticleMassName;
    ParticleTypeName       = CosmologySimulationParticleTypeName;
    for (dim = 0; dim < MAX_DIMENSION; dim++)
      ParticleVelocityNames[dim]   = CosmologySimulationParticleVelocityNames[dim];
  }
 
  // If there is a subgrid, use CosmologySimulationSubgridsAreStatic,
  // otherwise just set to false
 
  int SubgridsAreStatic = (TopGrid->NextGridNextLevel == NULL) ?
    FALSE : CosmologySimulationSubgridsAreStatic;
 
  // Call grid initializer.  Use TotalRefinement = -1 to flag real read
 
  int TotalRefinement = -1;
 
  // Loop over level zero grid
 
  HierarchyEntry *Temp = TopGrid;

  PrintMemoryUsage("Call G_CSIG");
 
  while (Temp != NULL) {
 
    if (Temp->GridData->CosmologySimulationInitializeGrid(gridnum,
			     CosmologySimulationOmegaBaryonNow,
			       CosmologySimulationOmegaCDMNow,
			       CosmologySimulationInitialTemperature,
			     DensityName, TotalEnergyName,
			       GasEnergyName, VelocityNames,
			       ParticlePositionName, ParticleVelocityName,
			       ParticleMassName, ParticleTypeName,
			     ParticleVelocityNames,
			     SubgridsAreStatic, TotalRefinement,
			     CosmologySimulationInitialFractionHII,
			     CosmologySimulationInitialFractionHeII,
			     CosmologySimulationInitialFractionHeIII,
			     CosmologySimulationInitialFractionHM,
			     CosmologySimulationInitialFractionH2I,
			     CosmologySimulationInitialFractionH2II,
#ifdef TRANSFER
			     RadHydroInitialRadiationEnergy,
#endif
			     CosmologySimulationUseMetallicityField,
			     MetaData.NumberOfParticles,
			     CosmologySimulationManuallySetParticleMassRatio,
			     CosmologySimulationManualParticleMassRatio,
			     CosmologySimulationCalculatePositions
						       ) == FAIL) {
      ENZO_FAIL("Error in grid->CosmologySimulationInitializeGrid.\n");
    }
 
    Temp = Temp->NextGridThisLevel;
  }

  PrintMemoryUsage("Called G_CSIG");
 
    //  Create tracer particles

    PINT DummyNumberOfParticles = 0;
 
    Temp = TopGrid;
 
    while (Temp != NULL) {
      if (Temp->GridData->TracerParticleCreateParticles(
		TracerParticleCreationLeftEdge,
		TracerParticleCreationRightEdge,
		TracerParticleCreationSpacing,
		DummyNumberOfParticles) == FAIL) {
	ENZO_FAIL("Error in grid->TracerParticleCreateParticles\n");
      }
      
      Temp = Temp->NextGridThisLevel;
    }
 
  // Get the global particle count
 
  int LocalNumberOfParticles;
 
  Temp = TopGrid;
  while (Temp != NULL) {
 
    LocalNumberOfParticles = 0;
    LocalNumberOfParticles = Temp->GridData->ReturnNumberOfParticles();
    // printf("OldLocalParticleCount: %"ISYM"\n", LocalNumberOfParticles );

#ifdef USE_MPI 
    CommunicationAllReduceValues(&LocalNumberOfParticles, 1, MPI_SUM);
#endif /* USE_MPI */
    Temp->GridData->SetNumberOfParticles(LocalNumberOfParticles);
 
    LocalNumberOfParticles = Temp->GridData->ReturnNumberOfParticles();
    // printf("NewLocalParticleCount: %"ISYM"\n", LocalNumberOfParticles );
 
    Temp = Temp->NextGridThisLevel;
  }
 
  PrintMemoryUsage("Local pc set");
 
 
  // Loop over grids and set particle ID number
 
  Temp = TopGrid;
  PINT ParticleCount = 0;
  RecursivelySetParticleCount(Temp, &ParticleCount);
  if (debug)
    printf("FinalParticleCount = %"ISYM"\n", ParticleCount);

  // 2006-12-11 Skory bug fix for star particle miscounts
  // Removed the following line:
  // MetaData.NumberOfParticles = 0;
  // Added the following line:
  MetaData.NumberOfParticles = ParticleCount;

  PrintMemoryUsage("Recursive pc set");

 
  return SUCCESS;
}
 
void RecursivelySetParticleCount(HierarchyEntry *GridPoint, PINT *Count)
{
  // Add Count to the particle id's on this grid (which start from zero
  // since we are doing a parallel root grid i/o)
 
  GridPoint->GridData->AddToParticleNumber(Count);
 
  // Recursively apply this to siblings and children
 
  if (GridPoint->NextGridThisLevel != NULL)
    RecursivelySetParticleCount(GridPoint->NextGridThisLevel, Count);
 
  if (GridPoint->NextGridNextLevel != NULL)

    RecursivelySetParticleCount(GridPoint->NextGridNextLevel, Count);
 
  return;
}
