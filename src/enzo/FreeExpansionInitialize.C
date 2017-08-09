/***********************************************************************
/
/  INITIALIZE FREE EXPANSION STAGE OF A BLAST WAVE
/
/  written by: John Wise
/  date:       August, 2009
/  modified1:  
/
/  PURPOSE:
/
/   REFERENCE: Draine & Woods, 1991, ApJ
/              Truelove & McKee, 1999, ApJS, 120, 299
/
/   Two dimensional parameters: explosion energy E and ambient density rho_1
/   Two independent variables: radius r, time t
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/

// This routine intializes a new simulation based on the parameter file.

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

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int FreeExpansionInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
			    TopGridData &MetaData)
{
  char	*DensName  = "Density";
  char	*TEName	   = "TotalEnergy";
  char	*GEName	   = "GasEnergy";
  char	*Vel1Name  = "x-velocity";
  char	*Vel2Name  = "y-velocity";
  char	*Vel3Name  = "z-velocity";
  char	*B1Name	   = "Bx";
  char	*B2Name	   = "By";
  char	*B3Name	   = "Bz";
  char	*PhiName   = "Phi";
  char	*Phi_pName = "Phip";

  /* parameter declarations */

  float FreeExpansionTotalEnergy, FreeExpansionGasEnergy, B2, V2;
  FLOAT FreeExpansionSubgridLeft, FreeExpansionSubgridRight;
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
  
  /* local declarations */

  char line[MAX_LINE_LENGTH];
  int  dim, ret, NumberOfSubgridZones[MAX_DIMENSION],
                          SubgridDims[MAX_DIMENSION];

  int   FreeExpansionFullBox     = FALSE;
  float FreeExpansionVelocity[3]   = {0.0, 0.0, 0.0};   // gas initally (t=0) at rest
  float FreeExpansionBField[3]   = {0.0, 0.0, 0.0};     // Gauss
  float FreeExpansionDensity    = 1.0;
  float FreeExpansionMaxVelocity = FLOAT_UNDEFINED;  // km/s
  float FreeExpansionRadius    = 0.1;
  float FreeExpansionMass    = 1.0;   // Msun
  double FreeExpansionEnergy        = 1e51;   // ergs
  float FreeExpansionTemperature = 100;  // K
  float dx = (DomainRightEdge[0] - DomainLeftEdge[0]) / MetaData.TopGridDims[0];

  /* Use 3.5 zones on the finest level to resolve the initial explosion at t=0. */

  float dr = 3.5*dx*max(POW(RefineBy,-MaximumRefinementLevel), 0.25);

  /* set no subgrids by default. */

  FreeExpansionSubgridLeft         = 0.0;    // start of subgrid(s)
  FreeExpansionSubgridRight        = 0.0;    // end of subgrid(s)

  /* read input from file */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */

    ret += sscanf(line, "FreeExpansionFullBox  = %"ISYM, &FreeExpansionFullBox);
    ret += sscanf(line, "FreeExpansionMass  = %"FSYM, &FreeExpansionMass);
    ret += sscanf(line, "FreeExpansionRadius  = %"FSYM, &FreeExpansionRadius);
    ret += sscanf(line, "FreeExpansionDensity  = %"FSYM, &FreeExpansionDensity);
    ret += sscanf(line, "FreeExpansionEnergy   = %"FSYM, &FreeExpansionEnergy);
    ret += sscanf(line, "FreeExpansionMaxVelocity   = %"FSYM, &FreeExpansionMaxVelocity);
    ret += sscanf(line, "FreeExpansionTemperature   = %"FSYM, &FreeExpansionTemperature);
    ret += sscanf(line, "FreeExpansionVelocity = %"FSYM" %"FSYM" %"FSYM, 
		  FreeExpansionVelocity, FreeExpansionVelocity+1, FreeExpansionVelocity+2);
    ret += sscanf(line, "FreeExpansionBField = %"FSYM" %"FSYM" %"FSYM, 
		  FreeExpansionBField, FreeExpansionBField+1, FreeExpansionBField+2);
    ret += sscanf(line, "FreeExpansionSubgridLeft = %"FSYM, 
		  &FreeExpansionSubgridLeft);
    ret += sscanf(line, "FreeExpansionSubgridRight = %"FSYM, 
		  &FreeExpansionSubgridRight);

    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") && strstr(line, "FreeExpansion") && 
	line[0] != '#' && MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr, 
	 "warning: the following parameter line was not interpreted:\n%s\n", 
	      line);

  } // end input from parameter file


  /* get units */

  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, 
    VelocityUnits, PressureUnits, MagneticUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, 
	   &VelocityUnits, 0.0);

  PressureUnits = DensityUnits * (LengthUnits/TimeUnits)*(LengthUnits/TimeUnits);
  MagneticUnits = sqrt(PressureUnits*4.0*M_PI);
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    FreeExpansionBField[dim] /= MagneticUnits;

  if (debug) 
    printf("Bunits = %"GSYM" G, Bfield(code) = %"GSYM" %"GSYM" %"GSYM"\n",
	   MagneticUnits, FreeExpansionBField[0], FreeExpansionBField[1],
	   FreeExpansionBField[2]);

  /* Set up current problem time, ambient total energy. */

  MetaData.Time         = 0.0;
  FreeExpansionGasEnergy = FreeExpansionTemperature / TemperatureUnits / 
    ((Gamma-1.0)*Mu);

  V2 = 0;
  for (dim = 0; dim < MetaData.TopGridRank; dim++)
    V2 += FreeExpansionVelocity[dim] * FreeExpansionVelocity[dim];
  FreeExpansionTotalEnergy = FreeExpansionGasEnergy + 0.5 * V2;

  if (HydroMethod == MHD_RK) {
    B2 = 0.0;
    for (dim = 0; dim < MetaData.TopGridRank; dim++)
      B2 += FreeExpansionBField[dim] * FreeExpansionBField[dim];
    FreeExpansionTotalEnergy += 0.5 * B2 / FreeExpansionDensity;
  }

  
  /* set periodic boundaries (if FullBox=1),
     otherwise keep reflecting (the default) */

  if (FreeExpansionFullBox)
    for (dim = 0; dim < MetaData.TopGridRank; dim++) {
      MetaData.LeftFaceBoundaryCondition[dim]  = periodic;
      MetaData.RightFaceBoundaryCondition[dim] = periodic;
    }

  /* set up uniform grid as of before explosion */

  TopGrid.GridData->InitializeUniformGrid(FreeExpansionDensity, 
					  FreeExpansionTotalEnergy,
					  FreeExpansionGasEnergy,
					  FreeExpansionVelocity,
					  FreeExpansionBField);

  /* Create as many subgrids as there are refinement levels 
     needed to resolve the initial explosion region upon the start-up. */

  HierarchyEntry ** Subgrid;
  if (MaximumRefinementLevel > 0) 
    Subgrid   = new HierarchyEntry*[MaximumRefinementLevel];

  /* Create new HierarchyEntries. */

  int lev;
  for (lev = 0; lev < MaximumRefinementLevel; lev++) 
    Subgrid[lev] = new HierarchyEntry;

  for (lev = 0; lev < MaximumRefinementLevel; lev++) {

    for (dim = 0; dim < MetaData.TopGridRank; dim++)
      NumberOfSubgridZones[dim] =
	nint((FreeExpansionSubgridRight - FreeExpansionSubgridLeft)/
	     ((DomainRightEdge[dim] - DomainLeftEdge[dim] )/
	      float(MetaData.TopGridDims[dim])))
        *POW(RefineBy, lev + 1);

    if (debug)
      printf("FreeExpansion:: Level[%"ISYM"]: NumberOfSubgridZones[0] = %"ISYM"\n", lev+1, 
	     NumberOfSubgridZones[0]);

    if (NumberOfSubgridZones[0] > 0) {

      /* fill them out */

      if (lev == 0)
	TopGrid.NextGridNextLevel  = Subgrid[0];
      Subgrid[lev]->NextGridThisLevel = NULL;
      if (lev == MaximumRefinementLevel-1)
	Subgrid[lev]->NextGridNextLevel = NULL;
      else
	Subgrid[lev]->NextGridNextLevel = Subgrid[lev+1];
      if (lev == 0)
	Subgrid[lev]->ParentGrid        = &TopGrid;
      else
	Subgrid[lev]->ParentGrid        = Subgrid[lev-1];

      /* compute the dimensions and left/right edges for the subgrid */

      for (dim = 0; dim < MetaData.TopGridRank; dim++) {
	SubgridDims[dim] = NumberOfSubgridZones[dim] + 2*NumberOfGhostZones;
	LeftEdge[dim]    = FreeExpansionSubgridLeft;
	RightEdge[dim]   = FreeExpansionSubgridRight;
      }

      /* create a new subgrid and initialize it */

      Subgrid[lev]->GridData = new grid;
      Subgrid[lev]->GridData->InheritProperties(TopGrid.GridData);
      Subgrid[lev]->GridData->PrepareGrid(MetaData.TopGridRank, SubgridDims,
					  LeftEdge, RightEdge, 0);
      Subgrid[lev]->GridData->InitializeUniformGrid(FreeExpansionDensity,
						    FreeExpansionTotalEnergy,
						    FreeExpansionGasEnergy,
						    FreeExpansionVelocity,
						    FreeExpansionBField);

      /* set up the initial explosion area on the finest resolution subgrid */

      if (lev == MaximumRefinementLevel - 1)
	Subgrid[lev]->GridData->
	  FreeExpansionInitializeGrid(FreeExpansionFullBox, FreeExpansionDensity,
				      FreeExpansionEnergy, FreeExpansionMaxVelocity,
				      FreeExpansionMass, FreeExpansionRadius, 
				      DensityUnits, VelocityUnits, LengthUnits,
				      TimeUnits);

    } else
      printf("FreeExpansion: single grid start-up.\n");
  } // ENDFOR levels

  /* set up subgrids from level 1 to max refinement level -1 */

  for (lev = MaximumRefinementLevel - 1; lev > 0; lev--)
    
    Subgrid[lev]->GridData->ProjectSolutionToParentGrid(*(Subgrid[lev-1]->GridData));
  
  /* set up the root grid */

  if (MaximumRefinementLevel > 0)
    Subgrid[0]->GridData->ProjectSolutionToParentGrid(*(TopGrid.GridData));
  else
    TopGrid.GridData->
      FreeExpansionInitializeGrid(FreeExpansionFullBox, FreeExpansionDensity,
				  FreeExpansionEnergy, FreeExpansionMaxVelocity,
				  FreeExpansionMass, FreeExpansionRadius, 
				  DensityUnits, VelocityUnits, LengthUnits, 
				  TimeUnits);

  /* set up field names and units */
  int i = 0;
  DataLabel[i++] = DensName;
  DataLabel[i++] = TEName;
  if (DualEnergyFormalism)
    DataLabel[i++] = GEName;
  DataLabel[i++] = Vel1Name;
  DataLabel[i++] = Vel2Name;
  DataLabel[i++] = Vel3Name;
  if (HydroMethod == MHD_RK) {
    DataLabel[i++] = B1Name;
    DataLabel[i++] = B2Name;
    DataLabel[i++] = B3Name;
    DataLabel[i++] = PhiName;
    if (UseDivergenceCleaning) {
      DataLabel[i++] = Phi_pName;
    }
  }

  for (int j=0; j< i; j++) 
    DataUnits[j] = NULL;

  /* Write parameters to parameter output file */

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fprintf(Outfptr, "FreeExpansionFullBox         = %"ISYM"\n", FreeExpansionFullBox);
    fprintf(Outfptr, "FreeExpansionDensity         = %"FSYM"\n", FreeExpansionDensity);
    fprintf(Outfptr, "FreeExpansionMass            = %"FSYM"\n", FreeExpansionMass);
    fprintf(Outfptr, "FreeExpansionRadius          = %"FSYM"\n", FreeExpansionRadius);
    fprintf(Outfptr, "FreeExpansionTemperature     = %"GSYM"\n", FreeExpansionTemperature);
    fprintf(Outfptr, "FreeExpansionEnergy          = %lg\n"    , FreeExpansionEnergy);
    fprintf(Outfptr, "FreeExpansionMaxVelocity     = %"FSYM"\n", FreeExpansionMaxVelocity);
  }

  return SUCCESS;

}
