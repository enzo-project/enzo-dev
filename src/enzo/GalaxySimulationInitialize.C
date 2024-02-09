/***********************************************************************
/
/  INITIALIZE A GALAXY SIMULATION
/
/  written by: Greg Bryan
/  date:       May, 1998
/  modified1:  Elizabeth Tasker, March 2004
/
/  PURPOSE:
/
/    Set up a number of spherical objects
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
#include "LevelHierarchy.h"
#include "TopGridData.h"

void WriteListOfFloats(FILE *fptr, int N, float floats[]);
void WriteListOfFloats(FILE *fptr, int N, double floats[]);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);

int GetUnits(float *DensityUnits, float *LengthUnits,
       float *TemperatureUnits, float *TimeUnits,
       float *VelocityUnits, double *MassUnits, FLOAT Time);

int ReadEquilibriumTable(char * name, FLOAT Time);

int GalaxySimulationInitialize(FILE *fptr, FILE *Outfptr, 
			  HierarchyEntry &TopGrid, TopGridData &MetaData, ExternalBoundary &Exterior)
{
  char *DensName    = "Density";
  char *TEName      = "TotalEnergy";
  char *GEName      = "GasEnergy";
  char *Vel1Name    = "x-velocity";
  char *Vel2Name    = "y-velocity";
  char *Vel3Name    = "z-velocity";
  char *CRName      = "CREnergyDensity";
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
  char *MetalName   = "Metal_Density";
  char *MetalIaName = "MetalSNIa_Density";
  char *BxName      = "Bx";
  char *ByName      = "By";
  char *BzName      = "Bz";
  char *PhiName     = "Phi";

  /* declarations */

  char  line[MAX_LINE_LENGTH];
  int   dim, ret, level, disk, i;

  /* make sure it is 3D */
  
  if (MetaData.TopGridRank != 3) {
    ENZO_VFAIL("Cannot do GalaxySimulation in %"ISYM" dimension(s)\n", MetaData.TopGridRank)
  }

  /* set default parameters */

  double GalaxySimulationGasMass,
    GalaxySimulationGalaxyMass,
    GalaxySimulationCR,
    GalaxySimulationDiskTemperature,
    GalaxySimulationAngularMomentum[MAX_DIMENSION],
    GalaxySimulationUniformVelocity[MAX_DIMENSION],
    GalaxySimulationUniformDensity,
    GalaxySimulationEquilibrateChem;

  FLOAT GalaxySimulationDiskPosition[MAX_DIMENSION];

  double GalaxySimulationDiskRadius,
    GalaxySimulationDiskScaleHeightz,
    GalaxySimulationDiskScaleHeightR,
    GalaxySimulationTruncationRadius,
    GalaxySimulationDiskDensityCap;

  int GalaxySimulationDiskPressureBalance;

  double GalaxySimulationInitialTemperature,
        GalaxySimulationDarkMatterConcentrationParameter,
        GalaxySimulationInflowDensity;

  FLOAT GalaxySimulationInflowTime;

  int GalaxySimulationGasHalo;
  double GalaxySimulationGasHaloScaleRadius,
        GalaxySimulationGasHaloDensity,
        GalaxySimulationGasHaloDensity2,
        GalaxySimulationGasHaloTemperature,
        GalaxySimulationGasHaloAlpha,
        GalaxySimulationGasHaloZeta,
        GalaxySimulationGasHaloZeta2,
        GalaxySimulationGasHaloCoreEntropy,
        GalaxySimulationGasHaloRatio,
        GalaxySimulationGasHaloMetallicity,
        GalaxySimulationDiskMetallicityEnhancementFactor;
  char GalaxySimulationEquilibriumFile[MAX_LINE_LENGTH] = "equilibrium_table_50_027-Zsun.h5"; 

  
  int GalaxySimulationGasHaloRotation;
  double GalaxySimulationGasHaloRotationScaleVelocity,
         GalaxySimulationGasHaloRotationScaleRadius,
         GalaxySimulationGasHaloRotationIndex;


  int   GalaxySimulationRefineAtStart,
    GalaxySimulationUseMetallicityField;
 
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];

  /* Default Values */

  GalaxySimulationRefineAtStart      = TRUE;
  GalaxySimulationUseMetallicityField  = FALSE;
  GalaxySimulationInitialTemperature = 1000.0;
  GalaxySimulationDiskRadius         = 0.2;      // CODE UNITS
  GalaxySimulationDiskTemperature    = 1.e4;     // [K]
  GalaxySimulationDiskPressureBalance = 0; // off
  GalaxySimulationDiskScaleHeightz   = 325e-6;   // Mpc
  GalaxySimulationDiskScaleHeightR   = 3500e-6;  // Mpc
  GalaxySimulationTruncationRadius   = .026; // [ Mpc ]
  GalaxySimulationDiskDensityCap     = 0.0;
  GalaxySimulationDarkMatterConcentrationParameter = 12;
  GalaxySimulationGasMass            = 4.0e10;
  GalaxySimulationGalaxyMass         = 1.0e12;
  GalaxySimulationDiskTemperature    = 1000.0;
  GalaxySimulationEquilibrateChem    = 0; // bool
  GalaxySimulationGasHalo            = 0; // uniform halo w/ densicm and UniformTemperature
  GalaxySimulationGasHaloScaleRadius = .001; // Mpc
  GalaxySimulationGasHaloDensity     = 1.8e-27; // cgs
  GalaxySimulationGasHaloDensity2    = 0.0; // cgs
  GalaxySimulationGasHaloTemperature = 1.0e+6;  // Kelvin
  GalaxySimulationGasHaloAlpha       = 2.0/3.0;  // unitless
  GalaxySimulationGasHaloZeta        = 0;
  GalaxySimulationGasHaloZeta2       = 0;
  GalaxySimulationGasHaloCoreEntropy = 5.0;  // keV cm^2
  GalaxySimulationGasHaloRatio       = 10; // ratio of cooling time to freefall time
  GalaxySimulationGasHaloMetallicity = 0.1; // Zsun
  GalaxySimulationGasHaloRotation    = 0; // off
  GalaxySimulationGasHaloRotationScaleVelocity = 180.0; // km/s
  GalaxySimulationGasHaloRotationScaleRadius   = 10.0; // kpc
  GalaxySimulationGasHaloRotationIndex         = 0.0; // unitless
  GalaxySimulationDiskMetallicityEnhancementFactor = 3.0; // w.r.t to halo metallicity
  GalaxySimulationInflowTime         = -1;
  GalaxySimulationInflowDensity      = 0;
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    GalaxySimulationDiskPosition[dim] = 0.5*(DomainLeftEdge[dim] +
					     DomainRightEdge[dim]);
    GalaxySimulationAngularMomentum[dim] = 0.0;
    GalaxySimulationUniformVelocity[dim] = 0.0;
  }
  GalaxySimulationUniformDensity   = 1.0E-28;
  GalaxySimulationCR = .01;

  /* read input from file */
  char *filename_holder = new char[MAX_LINE_LENGTH];
  filename_holder[0] = 0;
      
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
    
    ret = 0;
   
    ret += sscanf(line, "GalaxySimulationRefineAtStart = %"ISYM,
		  &GalaxySimulationRefineAtStart);
    ret += sscanf(line, "GalaxySimulationUseMetallicityField = %"ISYM,
		  &GalaxySimulationUseMetallicityField);
    ret += sscanf(line, "GalaxySimulationInitialTemperature = %"FSYM,
		  &GalaxySimulationInitialTemperature);
    ret += sscanf(line, "GalaxySimulationUniformDensity = %"FSYM,
		  &GalaxySimulationUniformDensity);
    ret += sscanf(line, "GalaxySimulationUniformVelocity = %"FSYM" %"FSYM" %"FSYM,
                  &GalaxySimulationUniformVelocity[0], &GalaxySimulationUniformVelocity[1],
                  &GalaxySimulationUniformVelocity[2]);
    ret += sscanf(line, "GalaxySimulationDiskRadius = %"PSYM,
		  &GalaxySimulationDiskRadius);
    ret += sscanf(line, "GalaxySimulationGalaxyMass = %"FSYM,
		  &GalaxySimulationGalaxyMass);
    ret += sscanf(line, "GalaxySimulationGasMass = %"FSYM,
		  &GalaxySimulationGasMass);
    ret += sscanf(line, "GalaxySimulationCR = %"FSYM,
		  &GalaxySimulationCR);
    ret += sscanf(line, "GalaxySimulationDiskPosition = %"PSYM" %"PSYM" %"PSYM, 
		  &GalaxySimulationDiskPosition[0],
		  &GalaxySimulationDiskPosition[1],
		  &GalaxySimulationDiskPosition[2]);
    ret += sscanf(line, "GalaxySimulationDiskScaleHeightz = %"PSYM,
		  &GalaxySimulationDiskScaleHeightz);
    ret += sscanf(line, "GalaxySimulationDiskScaleHeightR = %"PSYM,
		  &GalaxySimulationDiskScaleHeightR);
    ret += sscanf(line, "GalaxySimulationTruncationRadius = %"PSYM,
		  &GalaxySimulationTruncationRadius);
    ret += sscanf(line, "GalaxySimulationDiskDensityCap = %"FSYM,
		  &GalaxySimulationDiskDensityCap);    
    ret += sscanf(line, "GalaxySimulationDarkMatterConcentrationParameter = %"FSYM,
		  &GalaxySimulationDarkMatterConcentrationParameter);
    ret += sscanf(line, "GalaxySimulationDiskTemperature = %"FSYM,
		  &GalaxySimulationDiskTemperature);
    ret += sscanf(line, "GalaxySimulationDiskPressureBalance = %"ISYM,
      &GalaxySimulationDiskPressureBalance);
    ret += sscanf(line, "GalaxySimulationEquilibrateChem = %"FSYM,
		  &GalaxySimulationEquilibrateChem);
    if (sscanf(line, "GalaxySimulationEquilibriumFile = %s", filename_holder) == 1) {
      strcpy(GalaxySimulationEquilibriumFile, filename_holder);
      ret++;
    }
    ret += sscanf(line, "GalaxySimulationGasHalo = %"ISYM,
		  &GalaxySimulationGasHalo);
    ret += sscanf(line, "GalaxySimulationGasHaloScaleRadius = %"FSYM,
		  &GalaxySimulationGasHaloScaleRadius);
    ret += sscanf(line, "GalaxySimulationGasHaloDensity = %"FSYM,
		  &GalaxySimulationGasHaloDensity);
    ret += sscanf(line, "GalaxySimulationGasHaloDensity2 = %"FSYM,
		  &GalaxySimulationGasHaloDensity2);
    ret += sscanf(line, "GalaxySimulationGasHaloTemperature = %"FSYM,
		  &GalaxySimulationGasHaloTemperature);
    ret += sscanf(line, "GalaxySimulationGasHaloAlpha = %"FSYM,
		  &GalaxySimulationGasHaloAlpha);
    ret += sscanf(line, "GalaxySimulationGasHaloZeta = %"FSYM,
		  &GalaxySimulationGasHaloZeta);
    ret += sscanf(line, "GalaxySimulationGasHaloZeta2 = %"FSYM,
		  &GalaxySimulationGasHaloZeta2);
    ret += sscanf(line, "GalaxySimulationGasHaloCoreEntropy = %"FSYM,
		  &GalaxySimulationGasHaloCoreEntropy);
    ret += sscanf(line, "GalaxySimulationGasHaloRatio = %"FSYM,
		  &GalaxySimulationGasHaloRatio);
    ret += sscanf(line, "GalaxySimulationGasHaloMetallicity = %"FSYM,
		  &GalaxySimulationGasHaloMetallicity);
    ret += sscanf(line, "GalaxySimulationGasHaloRotation = %"ISYM,
		  &GalaxySimulationGasHaloRotation);
    ret += sscanf(line, "GalaxySimulationGasHaloRotationScaleVelocity = %"FSYM,
		  &GalaxySimulationGasHaloRotationScaleVelocity);
    ret += sscanf(line, "GalaxySimulationGasHaloRotationScaleRadius = %"FSYM,
		  &GalaxySimulationGasHaloRotationScaleRadius);
    ret += sscanf(line, "GalaxySimulationGasHaloRotationIndex = %"FSYM,
		  &GalaxySimulationGasHaloRotationIndex);
    ret += sscanf(line, "GalaxySimulationDiskMetallicityEnhancementFactor = %"FSYM,
		  &GalaxySimulationDiskMetallicityEnhancementFactor);
    ret += sscanf(line, "GalaxySimulationInflowTime = %"FSYM,
		  &GalaxySimulationInflowTime);
    ret += sscanf(line, "GalaxySimulationInflowDensity = %"FSYM,
		  &GalaxySimulationInflowDensity);
    ret += sscanf(line, "GalaxySimulationAngularMomentum = %"FSYM" %"FSYM" %"FSYM,
		  &GalaxySimulationAngularMomentum[0],
		  &GalaxySimulationAngularMomentum[1],
		  &GalaxySimulationAngularMomentum[2]);
    
    /* if the line is suspicious, issue a warning */
    if (ret == 0 && strstr(line, "=") && line[0] != '#' 
        && strstr(line, "GalaxySimulation") && !strstr(line, "DiskGravityDarkMatter")
        && !strstr(line,"RPSWind") && !strstr(line,"PreWind") 
        )
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);

  } // end input from parameter file

  // If using DiskGravity, make two GalaxySimulation parameters consistent
  if (DiskGravity > 0 && DiskGravityDarkMatterUseNFW) {
    GalaxySimulationGalaxyMass = DiskGravityDarkMatterMass;
    GalaxySimulationDarkMatterConcentrationParameter = DiskGravityDarkMatterConcentration;
  }

  if (GalaxySimulationEquilibrateChem)
    ReadEquilibriumTable(GalaxySimulationEquilibriumFile, MetaData.Time);

  delete [] filename_holder;

  /* fix wind values wrt units */
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits;
  double MassUnits;
  if (GetUnits(&DensityUnits, &LengthUnits,&TemperatureUnits, &TimeUnits,
               &VelocityUnits, &MassUnits, MetaData.Time) == FAIL){
    fprintf(stderr, "Error in GetUnits.\n");
    return FAIL;
  }
      
  /* Convert RPS parameters to code units */
  GalaxySimulationRPSWindDensity = GalaxySimulationRPSWindDensity/DensityUnits;
  GalaxySimulationRPSWindPressure = GalaxySimulationRPSWindPressure/DensityUnits/LengthUnits/LengthUnits*TimeUnits*TimeUnits;
  GalaxySimulationRPSWindVelocity[0] = GalaxySimulationRPSWindVelocity[0]/LengthUnits*TimeUnits;
  GalaxySimulationRPSWindVelocity[1] = GalaxySimulationRPSWindVelocity[1]/LengthUnits*TimeUnits;
  GalaxySimulationRPSWindVelocity[2] = GalaxySimulationRPSWindVelocity[2]/LengthUnits*TimeUnits;
  GalaxySimulationRPSWindShockSpeed = GalaxySimulationRPSWindShockSpeed/LengthUnits*TimeUnits;
  GalaxySimulationRPSWindDelay = GalaxySimulationRPSWindDelay/TimeUnits;

  /* Align gaseous and stellar disks */
  if( DiskGravity > 0 ){
    for( i = 0 ; i < MAX_DIMENSION ; i++ )
      DiskGravityAngularMomentum[i] = GalaxySimulationAngularMomentum[i];
  } // end DiskGravity if

  /* set up grid */

  if (TopGrid.GridData->GalaxySimulationInitializeGrid(GalaxySimulationDiskRadius,
    GalaxySimulationGalaxyMass, 
    GalaxySimulationGasMass,
    GalaxySimulationDiskPosition, 
    GalaxySimulationDiskScaleHeightz,
    GalaxySimulationDiskScaleHeightR,
    GalaxySimulationTruncationRadius, 
    GalaxySimulationDiskDensityCap,
    GalaxySimulationDarkMatterConcentrationParameter,
    GalaxySimulationDiskTemperature,
    GalaxySimulationDiskPressureBalance, 
    GalaxySimulationInitialTemperature,
    GalaxySimulationUniformDensity,
    GalaxySimulationEquilibrateChem,
    GalaxySimulationGasHalo,
    GalaxySimulationGasHaloScaleRadius,
    GalaxySimulationGasHaloDensity,
    GalaxySimulationGasHaloDensity2,
    GalaxySimulationGasHaloTemperature,
    GalaxySimulationGasHaloAlpha,
    GalaxySimulationGasHaloZeta,
    GalaxySimulationGasHaloZeta2,
    GalaxySimulationGasHaloCoreEntropy,
    GalaxySimulationGasHaloRatio,
    GalaxySimulationGasHaloMetallicity,
    GalaxySimulationGasHaloRotation,
    GalaxySimulationGasHaloRotationScaleVelocity,
    GalaxySimulationGasHaloRotationScaleRadius,
    GalaxySimulationGasHaloRotationIndex,
    GalaxySimulationDiskMetallicityEnhancementFactor,
    GalaxySimulationAngularMomentum,
    GalaxySimulationUniformVelocity,
    GalaxySimulationUseMetallicityField,
    GalaxySimulationInflowTime,
    GalaxySimulationInflowDensity,0,
    GalaxySimulationCR
    ) == FAIL) {
      ENZO_FAIL("Error in GalaxySimulationInitialize[Sub]Grid.");
    } else {
      TopGrid.GridData->_GalaxySimulationInitialization = 1;
    }

  /* Convert minimum initial overdensity for refinement to mass
     (unless MinimumMass itself was actually set). */

  for (i = 0; i < MAX_FLAGGING_METHODS; i++)
    if (MinimumMassForRefinement[i] == FLOAT_UNDEFINED) {
      MinimumMassForRefinement[i] = MinimumOverDensityForRefinement[i];
      for (dim = 0; dim < MetaData.TopGridRank; dim++)
        MinimumMassForRefinement[i] *=
        (DomainRightEdge[dim]-DomainLeftEdge[dim])/
        float(MetaData.TopGridDims[dim]);
    }

  /* If requested, refine the grid to the desired level. */

  if (GalaxySimulationRefineAtStart) {

    /* Declare, initialize and fill out the LevelArray. */

    LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];
    for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
      LevelArray[level] = NULL;
    AddLevel(LevelArray, &TopGrid, 0);

    /* Add levels to the maximum depth or until no new levels are created,
       and re-initialize the level after it is created. */

    for (level = 0; level < MaximumRefinementLevel; level++) {
      if (RebuildHierarchy(&MetaData, LevelArray, level) == FAIL) {
        fprintf(stderr, "Error in RebuildHierarchy.\n");
        return FAIL;
      }
      if (LevelArray[level+1] == NULL)
      	break;
      LevelHierarchyEntry *Temp = LevelArray[level+1];
      while (Temp != NULL) {

	      if (Temp->GridData->GalaxySimulationInitializeGrid(
          GalaxySimulationDiskRadius,
          GalaxySimulationGalaxyMass, 
          GalaxySimulationGasMass,
          GalaxySimulationDiskPosition, 
          GalaxySimulationDiskScaleHeightz,
          GalaxySimulationDiskScaleHeightR,
          GalaxySimulationTruncationRadius, 
          GalaxySimulationDiskDensityCap,
          GalaxySimulationDarkMatterConcentrationParameter,
          GalaxySimulationDiskTemperature, 
          GalaxySimulationDiskPressureBalance,
          GalaxySimulationInitialTemperature,
          GalaxySimulationUniformDensity,
          GalaxySimulationEquilibrateChem,
          GalaxySimulationGasHalo,
          GalaxySimulationGasHaloScaleRadius,
          GalaxySimulationGasHaloDensity,
          GalaxySimulationGasHaloDensity2,
          GalaxySimulationGasHaloTemperature,
          GalaxySimulationGasHaloAlpha,
          GalaxySimulationGasHaloZeta,
          GalaxySimulationGasHaloZeta2,
          GalaxySimulationGasHaloCoreEntropy,
          GalaxySimulationGasHaloRatio,
          GalaxySimulationGasHaloMetallicity,
          GalaxySimulationGasHaloRotation,
          GalaxySimulationGasHaloRotationScaleVelocity,
          GalaxySimulationGasHaloRotationScaleRadius,
          GalaxySimulationGasHaloRotationIndex,
          GalaxySimulationDiskMetallicityEnhancementFactor,
          GalaxySimulationAngularMomentum,
          GalaxySimulationUniformVelocity,
          GalaxySimulationUseMetallicityField,
          GalaxySimulationInflowTime,
          GalaxySimulationInflowDensity,level,
          GalaxySimulationCR
          ) == FAIL) {
	          ENZO_FAIL("Error in GalaxySimulationInitialize[Sub]Grid.");
        } else { // if initialization didn't fail, set flag
          Temp->GridData->_GalaxySimulationInitialization = 1;
        }
	    Temp = Temp->NextGridThisLevel;
      }
    } // end: loop over levels

    /* Loop back from the bottom, restoring the consistency among levels. */

    for (level = MaximumRefinementLevel; level > 0; level--) {
      LevelHierarchyEntry *Temp = LevelArray[level];
      while (Temp != NULL) {
	      if (Temp->GridData->ProjectSolutionToParentGrid(*LevelArray[level-1]->GridData) == FAIL) {
          fprintf(stderr, "Error in grid->ProjectSolutionToParentGrid.\n");
          return FAIL;
        }
	    Temp = Temp->NextGridThisLevel;
      }
    }
  } // end: if (GalaxySimulationRefineAtStart)

  /* If Galaxy is Subject to ICM Wind, Initialize the exterior */

  if ( GalaxySimulationRPSWind > 0 ) {
    Exterior.Prepare(TopGrid.GridData);
	
    const int MAX_BNDRY_VARS = 6;
    float InflowValue[MAX_BNDRY_VARS], Dummy[MAX_BNDRY_VARS];
    InflowValue[0] = GalaxySimulationRPSWindDensity;
    InflowValue[1] = GalaxySimulationRPSWindPressure/(Gamma-1.0)/GalaxySimulationRPSWindDensity;
    if (HydroMethod != 2) {
      InflowValue[1] = InflowValue[1] + 0.5*(   pow(GalaxySimulationRPSWindVelocity[0],2)
	                                            + pow(GalaxySimulationRPSWindVelocity[1],2)
	                                            + pow(GalaxySimulationRPSWindVelocity[2],2));
    }
    InflowValue[2] = GalaxySimulationRPSWindVelocity[0];
    InflowValue[3] = GalaxySimulationRPSWindVelocity[1];
    InflowValue[4] = GalaxySimulationRPSWindVelocity[2];
    if (GalaxySimulationUseMetallicityField)
      InflowValue[5] = 1.0e-10;
  
    if (Exterior.InitializeExternalBoundaryFace(0, inflow, outflow, InflowValue,
						Dummy) == FAIL) {
      fprintf(stderr, "Error in InitializeExternalBoundaryFace.\n");
      return FAIL;
    }
    if (MetaData.TopGridRank > 1)
      Exterior.InitializeExternalBoundaryFace(1, inflow, outflow,
					      InflowValue, Dummy);
    if (MetaData.TopGridRank > 2)
      Exterior.InitializeExternalBoundaryFace(2, inflow, outflow,
					      InflowValue, Dummy);
    Exterior.InitializeExternalBoundaryParticles(MetaData.ParticleBoundaryType);
	
    /* Set Global Variables for RPS Wind (see ExternalBoundary_SetGalaxySimulationBoundary.C)*/

    GalaxySimulationRPSWindDelay += TopGrid.GridData->ReturnTime();
    GalaxySimulationRPSWindTotalEnergy = InflowValue[1]; 	
    GalaxySimulationPreWindDensity     = GalaxySimulationUniformDensity/DensityUnits;
    GalaxySimulationPreWindTotalEnergy = GalaxySimulationInitialTemperature/TemperatureUnits/((Gamma-1.0)*0.6); 
    GalaxySimulationPreWindVelocity[0] = 0.0;
    GalaxySimulationPreWindVelocity[1] = 0.0;
    GalaxySimulationPreWindVelocity[2] = 0.0;
  }

  // If we used the Equilibrium Table, delete it
  if (GalaxySimulationEquilibrateChem){
    if (MultiSpecies) {
      delete [] EquilibriumTable.HI;
      delete [] EquilibriumTable.HII;
      delete [] EquilibriumTable.HeI;
      delete [] EquilibriumTable.HeII;
      delete [] EquilibriumTable.HeIII;
      delete [] EquilibriumTable.de;
      if (MultiSpecies > 1) {
        delete [] EquilibriumTable.HM;
        delete [] EquilibriumTable.H2I;
        delete [] EquilibriumTable.H2II;
      }
      if (MultiSpecies > 2) {
        delete [] EquilibriumTable.DI;
        delete [] EquilibriumTable.DII;
        delete [] EquilibriumTable.HDI;
      }
    }
  }

 /* set up field names and units */

 int count = 0;
 DataLabel[count++] = DensName;
 DataLabel[count++] = TEName;
 if (DualEnergyFormalism)
   DataLabel[count++] = GEName;
 DataLabel[count++] = Vel1Name;
 if(MetaData.TopGridRank > 1)
   DataLabel[count++] = Vel2Name;
 if(MetaData.TopGridRank > 2)
   DataLabel[count++] = Vel3Name;
  if( UseMHD ){
      DataLabel[count++] = BxName;
      DataLabel[count++] = ByName;
      DataLabel[count++] = BzName;
  }
  if (HydroMethod == MHD_RK){
      DataLabel[count++] = PhiName;
  }
 if(CRModel)
   DataLabel[count++] = CRName;
 if (MultiSpecies) {
   DataLabel[count++] = ElectronName;
   DataLabel[count++] = HIName;
   DataLabel[count++] = HIIName;
   DataLabel[count++] = HeIName;
   DataLabel[count++] = HeIIName;
   DataLabel[count++] = HeIIIName;
   if (MultiSpecies > 1) {
     DataLabel[count++] = HMName;
     DataLabel[count++] = H2IName;
     DataLabel[count++] = H2IIName;
   }
   if (MultiSpecies > 2) {
     DataLabel[count++] = DIName;
     DataLabel[count++] = DIIName;
     DataLabel[count++] = HDIName;
   }
 }
 if (GalaxySimulationUseMetallicityField)
   DataLabel[count++] = MetalName;
 if (StarMakerTypeIaSNe)
   DataLabel[count++] = MetalIaName;

 for (i = 0; i < count; i++)
   DataUnits[i] = NULL;

 /* Write parameters to parameter output file */

 if (MyProcessorNumber == ROOT_PROCESSOR) {

   fprintf(Outfptr, "GalaxySimulationRefineAtStart      = %"ISYM"\n",
	   GalaxySimulationRefineAtStart);
   fprintf(Outfptr, "GalaxySimulationUseMetallicityField          = %"ISYM"\n",
	   GalaxySimulationUseMetallicityField);
   fprintf(Outfptr, "GalaxySimulationInitialTemperature = %"GOUTSYM"\n",
	   GalaxySimulationInitialTemperature);
   fprintf(Outfptr, "GalaxySimulationUniformDensity = %"GOUTSYM"\n",
     GalaxySimulationUniformDensity);
   fprintf(Outfptr, "GalaxySimulationUniformVelocity    = %"GOUTSYM" %"GOUTSYM" %"GOUTSYM"\n",
	   GalaxySimulationUniformVelocity[0], GalaxySimulationUniformVelocity[1],
	   GalaxySimulationUniformVelocity[2]);
   fprintf(Outfptr, "GalaxySimulationDiskRadius = %"GOUTSYM"\n",
	   GalaxySimulationDiskRadius);
   fprintf(Outfptr, "GalaxySimulationGalaxyMass = %"GOUTSYM"\n",
	   GalaxySimulationGalaxyMass);
   fprintf(Outfptr, "GalaxySimulationGasMass = %"GOUTSYM"\n",
	   GalaxySimulationGasMass);
   fprintf(Outfptr, "GalaxySimulationCR = %"GOUTSYM"\n",
     GalaxySimulationCR);
   fprintf(Outfptr, "GalaxySimulationDiskScaleHeightz = %"GOUTSYM"\n",
	   GalaxySimulationDiskScaleHeightz);
   fprintf(Outfptr, "GalaxySimulationDiskScaleHeightR = %"GOUTSYM"\n",
	   GalaxySimulationDiskScaleHeightR);
   fprintf(Outfptr, "GalaxySimulationTruncationRadius = %"GOUTSYM"\n",
     GalaxySimulationTruncationRadius);
   fprintf(Outfptr, "GalaxySimulationDarkMatterConcentrationParameter = %"GOUTSYM"\n",
	   GalaxySimulationDarkMatterConcentrationParameter);
   fprintf(Outfptr, "GalaxySimulationDiskTemperature = %"GOUTSYM"\n",
	   GalaxySimulationDiskTemperature);
   fprintf(Outfptr, "GalaxySimulationGasHalo = %"ISYM"\n",
     GalaxySimulationGasHalo);
   fprintf(Outfptr, "GalaxySimulationGasHaloScaleRadius = %"GOUTSYM"\n",
     GalaxySimulationGasHaloScaleRadius);
   fprintf(Outfptr, "GalaxySimulationGasHaloDensity = %"GOUTSYM"\n",
     GalaxySimulationGasHaloDensity);
   fprintf(Outfptr, "GalaxySimulationGasHaloDensity2 = %"GOUTSYM"\n",
     GalaxySimulationGasHaloDensity2);
   fprintf(Outfptr, "GalaxySimulationGasHaloTemperature = %"GOUTSYM"\n",
     GalaxySimulationGasHaloTemperature);
   fprintf(Outfptr, "GalaxySimulationGasHaloAlpha = %"GOUTSYM"\n",
     GalaxySimulationGasHaloAlpha);
   fprintf(Outfptr, "GalaxySimulationGasHaloZeta = %"GOUTSYM"\n",
     GalaxySimulationGasHaloZeta);
   fprintf(Outfptr, "GalaxySimulationGasHaloZeta2 = %"GOUTSYM"\n",
     GalaxySimulationGasHaloZeta2);
   fprintf(Outfptr, "GalaxySimulationGasHaloCoreEntropy = %"GOUTSYM"\n",
     GalaxySimulationGasHaloCoreEntropy);
   fprintf(Outfptr, "GalaxySimulationGasHaloMetallicity = %"GOUTSYM"\n",
     GalaxySimulationGasHaloMetallicity);
   fprintf(Outfptr, "GalaxySimulationDiskMetallicityEnhancementFactor = %"GOUTSYM"\n",
     GalaxySimulationDiskMetallicityEnhancementFactor);
   fprintf(Outfptr, "GalaxySimulationInflowTime = %"GOUTSYM"\n",
	   GalaxySimulationInflowTime);
   fprintf(Outfptr, "GalaxySimulationInflowDensity = %"GOUTSYM"\n",
	   GalaxySimulationInflowDensity);
   fprintf(Outfptr, "GalaxySimulationDiskPosition = ");
   WriteListOfFloats(Outfptr, MetaData.TopGridRank, GalaxySimulationDiskPosition);
   fprintf(Outfptr, "GalaxySimulationAngularMomentum = ");
   WriteListOfFloats(Outfptr, MetaData.TopGridRank, GalaxySimulationAngularMomentum);
 }

#ifdef USE_MPI

 // BWO: this forces the synchronization of the various point source gravity
 // parameters between processors.  If this is not done, things go to pieces!

 MPI_Barrier(MPI_COMM_WORLD);
 MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
 MPI_Bcast(&PointSourceGravityConstant,1,DataType,ROOT_PROCESSOR, MPI_COMM_WORLD);
 MPI_Bcast(&PointSourceGravityCoreRadius,1,DataType,ROOT_PROCESSOR, MPI_COMM_WORLD);

#endif

 return SUCCESS;

}
