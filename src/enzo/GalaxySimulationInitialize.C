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
#include "CommunicationUtilities.h"
#include "phys_constants.h"

void WriteListOfFloats(FILE *fptr, int N, float floats[]);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);

int GetUnits(float *DensityUnits, float *LengthUnits,
       float *TemperatureUnits, float *TimeUnits,
       float *VelocityUnits, double *MassUnits, FLOAT Time);

char* ChemicalSpeciesBaryonFieldLabel(const int &atomic_number, int element_set=1);

void RecursivelySetParticleCount(HierarchyEntry *GridPoint, PINT *Count);

int IndividualStarProperties_Initialize(TopGridData &MetaData);
int IndividualStarRadiationProperties_Initialize(void);
int InitializeStellarYields(const float &time);


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
  char *GravPotentialName = "GravPotential";


  char *ElectronName = "Electron_Density";
  char *HIName       = "HI_Density";
  char *HIIName      = "HII_Density";
  char *HeIName      = "HeI_Density";
  char *HeIIName     = "HeII_Density";
  char *HeIIIName    = "HeIII_Density";
  char *HMName       = "HM_Density";
  char *H2IName      = "H2I_Density";
  char *H2IIName     = "H2II_Density";
  char *DIName       = "DI_Density";
  char *DIIName      = "DII_Density";
  char *HDIName      = "HDI_Density";

  char *MetalName       = "Metal_Density";
//  char *MetallicityName = "Metallicity";
  char *MetalIaName     = "MetalSNIa_Density";

  char *AGBMetalName    = "AGB_Metal_Density";
  char *PopIIIMetalName = "PopIII_Metal_Density";
  char *PopIIIPISNeMetalName = "PopIII_PISNe_Metal_Density";
  char *SNIIMetalName   = "SNII_Metal_Density";
  char *SNIaMetalName   = "SNIa_Metal_Density";
  char *RProcMetalName  = "RProcess_Metal_Density";

  char *ExtraMetalName0    = "SNIa_sCh_Metal_Density";
  char *ExtraMetalName1    = "SNIa_SDS_Metal_Density";
  char *ExtraMetalName2    = "SNIa_HeRS_Metal_Density";


  /* Chemical Tracers */
  /* handled with lookup table - see below */

  /* declarations */

  char  line[MAX_LINE_LENGTH];
  int   dim, ret, level, disk, i;

  /* make sure it is 3D */

  if (MetaData.TopGridRank != 3) {
    ENZO_VFAIL("Cannot do GalaxySimulation in %"ISYM" dimension(s)\n", MetaData.TopGridRank)
  }

  /* set default parameters */

  float GalaxySimulationGasMass,
    GalaxySimulationGalaxyMass,
    GalaxySimulationCR,
    GalaxySimulationDiskTemperature,
    GalaxySimulationAngularMomentum[MAX_DIMENSION],
    GalaxySimulationUniformVelocity[MAX_DIMENSION],
    GalaxySimulationUniformDensity,
    GalaxySimulationUniformCR,
    GalaxySimulationUniformEnergy;

  FLOAT GalaxySimulationDiskRadius,
    GalaxySimulationDiskPosition[MAX_DIMENSION],
    GalaxySimulationDiskScaleHeightz,
    GalaxySimulationDiskScaleHeightR,
    GalaxySimulationTruncationRadius;


  float GalaxySimulationInitialTemperature,
    GalaxySimulationDarkMatterConcentrationParameter,
    GalaxySimulationInflowTime,
    GalaxySimulationInflowDensity;
	
	int GalaxySimulationGasHalo;
	float GalaxySimulationGasHaloScaleRadius,
		GalaxySimulationGasHaloDensity;

  int   GalaxySimulationRefineAtStart,
    GalaxySimulationUseMetallicityField,
    GalaxySimulationUseDensityPerturbation,
    GalaxySimulationDarkMatterParticles,
    GalaxySimulationSMAUGIC;

  float GalaxySimulationPerturbationFraction;
 
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
  float ZeroBField[3] = {0.0, 0.0, 0.0};

  float *DMParticleMass;
  float *DMParticlePosition[MAX_DIMENSION], *DMParticleVelocity[MAX_DIMENSION];

  /* Chemical tracers */
  float GalaxySimulationInitialDiskMetallicity, GalaxySimulationInitialHaloMetallicity;

  /* Default Values */

  GalaxySimulationRefineAtStart      = TRUE;
  GalaxySimulationUseMetallicityField  = FALSE;
  GalaxySimulationInitialTemperature = 1000.0;
  GalaxySimulationDiskRadius         = 0.2;      // CODE UNITS
  GalaxySimulationDiskTemperature    = 1.e4;     // [K]
  GalaxySimulationDiskScaleHeightz   = 325e-6;   // Mpc_cm
  GalaxySimulationDiskScaleHeightR   = 3500e-6;  // Mpc_cm
  GalaxySimulationTruncationRadius   = .026; // [ Mpc_cm ]
  GalaxySimulationDarkMatterConcentrationParameter = 12;
  GalaxySimulationGasMass            = 4.0e10;
  GalaxySimulationGalaxyMass         = 1.0e12;
  GalaxySimulationDiskTemperature    = 1000.0;   // AJE whhy is this repeated??
  GalaxySimulationGasHalo            = 0; // uniform halo w/ densicm and UniformTemperature
  GalaxySimulationGasHaloScaleRadius = .001; // Mpc_cm
  GalaxySimulationGasHaloDensity     = 1.8e-27; // cgs
  GalaxySimulationInflowTime         = -1;
  GalaxySimulationInflowDensity      = 0;
  GalaxySimulationUseDensityPerturbation = 0;
  GalaxySimulationPerturbationFraction   = 0.5;
  GalaxySimulationDarkMatterParticles = FALSE;
  GalaxySimulationSMAUGIC             = FALSE;

  DMParticleMass = NULL;
  for (int dim = 0; dim < MAX_DIMENSION; dim++){
    DMParticlePosition[dim] = NULL;
    DMParticleVelocity[dim] = NULL;
  }

  /* Chemical tracer defaults in SetDefaultGlobalValues.C; set to tiny_number */

  /* Set default abundances for halo and galaxy disk. 
     Default is primordial with 100% neutral galaxy ISM and 100% ionized halo */

  // these are the disk abundances
  TestProblemData.InnerHydrogenFractionByMass         = 0.75;
  TestProblemData.HII_Fraction_Inner                  = tiny_number;
  TestProblemData.HeII_Fraction_Inner                 = tiny_number;
  TestProblemData.HeIII_Fraction_Inner                = tiny_number;
  TestProblemData.HM_Fraction_Inner                   = tiny_number;
  TestProblemData.H2I_Fraction_Inner                  = tiny_number;
  TestProblemData.H2II_Fraction_Inner                 = tiny_number;
  TestProblemData.InnerDeuteriumToHydrogenRatio        = tiny_number;

  // these are the halo abundances
  TestProblemData.HydrogenFractionByMass              = 0.75;
  TestProblemData.HII_Fraction                        = 1.00;
  TestProblemData.HeII_Fraction                       = tiny_number;
  TestProblemData.HeIII_Fraction                      = 1.00;
  TestProblemData.HM_Fraction                         = tiny_number;
  TestProblemData.H2I_Fraction                        = tiny_number;
  TestProblemData.H2II_Fraction                       = tiny_number;
  TestProblemData.DeuteriumToHydrogenRatio             = tiny_number;

  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    GalaxySimulationDiskPosition[dim] = 0.5*(DomainLeftEdge[dim] +
					     DomainRightEdge[dim]);
    GalaxySimulationAngularMomentum[dim] = 0.0;
    GalaxySimulationUniformVelocity[dim] = 0.0;
  }
  GalaxySimulationUniformDensity = 1.0E-28;
  GalaxySimulationUniformEnergy = 1.0;
  GalaxySimulationCR = .01;
  GalaxySimulationUniformCR = .01;


  GalaxySimulationInitialDiskMetallicity = tiny_number;
  GalaxySimulationInitialHaloMetallicity = tiny_number;
  /* read input from file */

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
    ret += sscanf(line, "GalaxySimulationUniformCR = %"FSYM,
		  &GalaxySimulationUniformCR);
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
    ret += sscanf(line, "GalaxySimulationDarkMatterConcentrationParameter = %"FSYM,
		  &GalaxySimulationDarkMatterConcentrationParameter);
    ret += sscanf(line, "GalaxySimulationDiskTemperature = %"FSYM,
		  &GalaxySimulationDiskTemperature);
    ret += sscanf(line, "GalaxySimulationGasHalo = %"ISYM,
		  &GalaxySimulationGasHalo);
    ret += sscanf(line, "GalaxySimulationGasHaloScaleRadius = %"FSYM,
		  &GalaxySimulationGasHaloScaleRadius);
    ret += sscanf(line, "GalaxySimulationGasHaloDensity = %"FSYM,
		  &GalaxySimulationGasHaloDensity);
    ret += sscanf(line, "GalaxySimulationInflowTime = %"FSYM,
		  &GalaxySimulationInflowTime);
    ret += sscanf(line, "GalaxySimulationInflowDensity = %"FSYM,
		  &GalaxySimulationInflowDensity);
    ret += sscanf(line, "GalaxySimulationAngularMomentum = %"FSYM" %"FSYM" %"FSYM,
		  &GalaxySimulationAngularMomentum[0],
		  &GalaxySimulationAngularMomentum[1],
		  &GalaxySimulationAngularMomentum[2]);

    ret += sscanf(line, "GalaxySimulationUseDensityPerturbation = %"ISYM,
                        &GalaxySimulationUseDensityPerturbation);
    ret += sscanf(line, "GalaxySimulationPerturbationFraction = %"FSYM,
                        &GalaxySimulationPerturbationFraction);

    /* Initial abundances for the galaxy disk */
    ret += sscanf(line, "GalaxySimulationHydrogenFractionByMass = %"FSYM,
                        &TestProblemData.InnerHydrogenFractionByMass);
    ret += sscanf(line, "GalaxySimulationHIIFraction = %"FSYM,
                        &TestProblemData.HII_Fraction_Inner);
    ret += sscanf(line, "GalaxySimulationHeIIFraction = %"FSYM,
                        &TestProblemData.HeII_Fraction_Inner);
    ret += sscanf(line, "GalaxySimulationHeIIIFraction = %"FSYM,
                        &TestProblemData.HeIII_Fraction_Inner);
    ret += sscanf(line, "GalaxySimulationHMFraction = %"FSYM,
                        &TestProblemData.HM_Fraction_Inner);
    ret += sscanf(line, "GalaxySimulationH2IFraction = %"FSYM,
                        &TestProblemData.H2I_Fraction_Inner);
    ret += sscanf(line, "GalaxySimulationH2IIFraction = %"FSYM,
                        &TestProblemData.H2II_Fraction_Inner);
    ret += sscanf(line, "GalaxySimulationDeuteriumToHydrogenRatio = %"FSYM,
                        &TestProblemData.InnerDeuteriumToHydrogenRatio);

    /* Initial abundances for the galaxy halo */
    ret += sscanf(line, "GalaxySimulationHydrogenFractionByMassHalo = %"FSYM,
                        &TestProblemData.HydrogenFractionByMass);
    ret += sscanf(line, "GalaxySimulationHIIFractionHalo = %"FSYM,
                        &TestProblemData.HII_Fraction);
    ret += sscanf(line, "GalaxySimulationHeIIFractionHalo = %"FSYM,
                        &TestProblemData.HeII_Fraction);
    ret += sscanf(line, "GalaxySimulationHeIIIFractionHalo = %"FSYM,
                        &TestProblemData.HeIII_Fraction);
    ret += sscanf(line, "GalaxySimulationHMFractionHalo = %"FSYM,
                        &TestProblemData.HM_Fraction);
    ret += sscanf(line, "GalaxySimulationH2IFractionHalo = %"FSYM,
                        &TestProblemData.H2I_Fraction);
    ret += sscanf(line, "GalaxySimulationH2IIFractionHalo = %"FSYM,
                        &TestProblemData.H2II_Fraction);
    ret += sscanf(line, "GalaxySimulationDeuteriumToHydrogenRatioHalo = %"FSYM,
                        &TestProblemData.DeuteriumToHydrogenRatio);


    /* Read in chemical tracer IC's */
    ret += sscanf(line, "GalaxySimulationInitialCIFraction = %"FSYM,
                        &TestProblemData.CI_Fraction);
    ret += sscanf(line, "GalaxySimulationInitialNIFraction = %"FSYM,
                        &TestProblemData.NI_Fraction);
    ret += sscanf(line, "GalaxySimulationInitialOIFraction = %"FSYM,
                        &TestProblemData.OI_Fraction);
    ret += sscanf(line, "GalaxySimulationInitialMgIFraction = %"FSYM,
                        &TestProblemData.MgI_Fraction);
    ret += sscanf(line, "GalaxySimulationInitialSiIFraction = %"FSYM,
                        &TestProblemData.SiI_Fraction);
    ret += sscanf(line, "GalaxySimulationInitialFeIFraction = %"FSYM,
                        &TestProblemData.FeI_Fraction);
    ret += sscanf(line, "GalaxySimulationInitialYIFraction = %"FSYM,
                        &TestProblemData.YI_Fraction);
    ret += sscanf(line, "GalaxySimulationInitialBaIFraction = %"FSYM,
                        &TestProblemData.BaI_Fraction);
    ret += sscanf(line, "GalaxySimulationInitialLaIFraction = %"FSYM,
                        &TestProblemData.LaI_Fraction);
    ret += sscanf(line, "GalaxySimulationInitialEuIFraction = %"FSYM,
                        &TestProblemData.EuI_Fraction);


    /* Initial Chemical tracer values for halo */
    ret += sscanf(line, "GalaxySimulationInitialCIFractionHalo = %"FSYM,
                        &TestProblemData.CI_Fraction_2);
    ret += sscanf(line, "GalaxySimulationInitialNIFractionHalo = %"FSYM,
                        &TestProblemData.NI_Fraction_2);
    ret += sscanf(line, "GalaxySimulationInitialOIFractionHalo = %"FSYM,
                        &TestProblemData.OI_Fraction_2);
    ret += sscanf(line, "GalaxySimulationInitialMgIFractionHalo = %"FSYM,
                        &TestProblemData.MgI_Fraction_2);
    ret += sscanf(line, "GalaxySimulationInitialSiIFractionHalo = %"FSYM,
                        &TestProblemData.SiI_Fraction_2);
    ret += sscanf(line, "GalaxySimulationInitialFeIFractionHalo = %"FSYM,
                        &TestProblemData.FeI_Fraction_2);
    ret += sscanf(line, "GalaxySimulationInitialYIFractionHalo = %"FSYM,
                        &TestProblemData.YI_Fraction_2);
    ret += sscanf(line, "GalaxySimulationInitialBaIFractionHalo = %"FSYM,
                        &TestProblemData.BaI_Fraction_2);
    ret += sscanf(line, "GalaxySimulationInitialLaIFractionHalo = %"FSYM,
                        &TestProblemData.LaI_Fraction_2);
    ret += sscanf(line, "GalaxySimulationInitialEuIFractionHalo = %"FSYM,
                        &TestProblemData.EuI_Fraction_2);

   ret += sscanf(line, "GalaxySimulationInitialSpeciesFractionsDisk = %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 0,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 1,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 2,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 3,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 4,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 5,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 6,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 7,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 8,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 9,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 10,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 11,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 12,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 13,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 14,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 15,
                        TestProblemData.ChemicalTracerSpecies_Fractions + 16 );

   ret += sscanf(line, "GalaxySimulationInitialSpeciesFractionsHalo = %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM,
                        TestProblemData.ChemicalTracerSpecies_Fractions_2 + 0,
                        TestProblemData.ChemicalTracerSpecies_Fractions_2 + 1,
                        TestProblemData.ChemicalTracerSpecies_Fractions_2 + 2,
                        TestProblemData.ChemicalTracerSpecies_Fractions_2 + 3,
                        TestProblemData.ChemicalTracerSpecies_Fractions_2 + 4,
                        TestProblemData.ChemicalTracerSpecies_Fractions_2 + 5,
                        TestProblemData.ChemicalTracerSpecies_Fractions_2 + 6,
                        TestProblemData.ChemicalTracerSpecies_Fractions_2 + 7,
                        TestProblemData.ChemicalTracerSpecies_Fractions_2 + 8,
                        TestProblemData.ChemicalTracerSpecies_Fractions_2 + 9,
                        TestProblemData.ChemicalTracerSpecies_Fractions_2 + 10,
                        TestProblemData.ChemicalTracerSpecies_Fractions_2 + 11,
                        TestProblemData.ChemicalTracerSpecies_Fractions_2 + 12,
                        TestProblemData.ChemicalTracerSpecies_Fractions_2 + 13,
                        TestProblemData.ChemicalTracerSpecies_Fractions_2 + 14,
                        TestProblemData.ChemicalTracerSpecies_Fractions_2 + 15,
                        TestProblemData.ChemicalTracerSpecies_Fractions_2 + 16 );



    ret += sscanf(line, "GalaxySimulationDarkMatterParticles = %"ISYM,
                        &GalaxySimulationDarkMatterParticles);

    ret += sscanf(line, "GalaxySimulationSMAUGIC = %"ISYM,
                        &GalaxySimulationSMAUGIC);

    ret += sscanf(line, "TestProblemUseMetallicityField = %"ISYM,
                        &TestProblemData.UseMetallicityField);

    ret += sscanf(line, "GalaxySimulationInitialDiskMetallicity = %"FSYM,
                        &GalaxySimulationInitialDiskMetallicity);

    ret += sscanf(line, "GalaxySimulationInitialHaloMetallicity = %"FSYM,
                        &GalaxySimulationInitialHaloMetallicity);

    /* if the line is suspicious, issue a warning */
    if (ret == 0 && strstr(line, "=") && strstr(line, "GalaxySimulation") 
	&& line[0] != '#' && !strstr(line,"RPSWind") && !strstr(line,"PreWind"))
      fprintf(stderr, "warning: the following parameter line was not interpreted:\n%s\n", line);

  } // end input from parameter file


  /* fix wind values wrt units */
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits;
  double MassUnits;
  if (GetUnits(&DensityUnits, &LengthUnits,&TemperatureUnits, &TimeUnits,
               &VelocityUnits, &MassUnits, MetaData.Time) == FAIL){
    fprintf(stderr, "Error in GetUnits.\n");
    return FAIL;
  }
  GalaxySimulationRPSWindDensity = GalaxySimulationRPSWindDensity/DensityUnits;
  GalaxySimulationRPSWindPressure = GalaxySimulationRPSWindPressure/DensityUnits/LengthUnits/LengthUnits*TimeUnits*TimeUnits;
  GalaxySimulationRPSWindVelocity[0] = GalaxySimulationRPSWindVelocity[0]/LengthUnits*TimeUnits;
  GalaxySimulationRPSWindVelocity[1] = GalaxySimulationRPSWindVelocity[1]/LengthUnits*TimeUnits;
  GalaxySimulationRPSWindVelocity[2] = GalaxySimulationRPSWindVelocity[2]/LengthUnits*TimeUnits;
  GalaxySimulationRPSWindShockSpeed = GalaxySimulationRPSWindShockSpeed/LengthUnits*TimeUnits;
  GalaxySimulationRPSWindDelay = GalaxySimulationRPSWindDelay/TimeUnits;

  TestProblemData.MultiSpecies = MultiSpecies;

  /* initialize dark matter */
  /* Align gaseous and stellar disks */
  if( DiskGravity > 0 ){
    for( i = 0 ; i < MAX_DIMENSION ; i++ )
      DiskGravityAngularMomentum[i] = GalaxySimulationAngularMomentum[i];


    // set central density if DM mass given
    if( DiskGravityDarkMatterDensity < 0){
      float xtemp = DiskGravityDarkMatterMassInteriorR / DiskGravityDarkMatterR; // convenience for below
      DiskGravityDarkMatterDensity = DiskGravityDarkMatterMassInterior*SolarMass /
                                     ( (2.0 * pi * POW(DiskGravityDarkMatterR*Mpc_cm,3.0) ) *
                                       (0.5 * log(1.0 + xtemp*xtemp) + log(1.0 + xtemp) - atan(xtemp)));

    }


  } // end DiskGravity if

  /* set up grid */

  /* try something here */
  HierarchyEntry *Temp = &TopGrid;

  while (Temp != NULL){

    if (TopGrid.GridData->GalaxySimulationInitializeGrid(GalaxySimulationDiskRadius,
						       GalaxySimulationGalaxyMass, 
						       GalaxySimulationGasMass,
						       GalaxySimulationDiskPosition, 
						       GalaxySimulationDiskScaleHeightz,
						       GalaxySimulationDiskScaleHeightR,
						       GalaxySimulationTruncationRadius, 
						       GalaxySimulationDarkMatterConcentrationParameter,
						       GalaxySimulationDiskTemperature,
                                                       GalaxySimulationInitialTemperature,
                                                       GalaxySimulationInitialDiskMetallicity,
                                                       GalaxySimulationInitialHaloMetallicity,
						       GalaxySimulationUniformDensity,
						       GalaxySimulationGasHalo,
						       GalaxySimulationGasHaloScaleRadius,
						       GalaxySimulationGasHaloDensity,
						       GalaxySimulationAngularMomentum,
						       GalaxySimulationUniformVelocity,
						       GalaxySimulationUseMetallicityField,
						       GalaxySimulationInflowTime,
						       GalaxySimulationInflowDensity,0,
						       GalaxySimulationCR,
                                                       GalaxySimulationUseDensityPerturbation,
                                                       GalaxySimulationPerturbationFraction,
                                                       GalaxySimulationSMAUGIC)
	      == FAIL) {
      ENZO_FAIL("Error in GalaxySimulationInitialize[Sub]Grid.");
    }// end subgrid if
    Temp = Temp->NextGridThisLevel;

  } // end top grid loop

  /* If we are using dark matter particles, loop through grids depositing the particles */
  if (GalaxySimulationDarkMatterParticles){
    const int MAXIMUM_NUMBER_OF_INITIAL_PARTICLES = 2000000;

    /* Read in the particles */
    DMParticleMass = new float[MAXIMUM_NUMBER_OF_INITIAL_PARTICLES];
    for(int dim = 0; dim < MAX_DIMENSION; dim++){
      DMParticlePosition[dim] = new float[MAXIMUM_NUMBER_OF_INITIAL_PARTICLES];
      DMParticleVelocity[dim] = new float[MAXIMUM_NUMBER_OF_INITIAL_PARTICLES];
    }

    char line[MAX_LINE_LENGTH];
    int err;

    FILE *fptr = fopen("GalaxySimulationParticleIC.in", "r");
    if (fptr == NULL){
      ENZO_FAIL("Error opening galaxy simulation dark matter particle positions\n");
    }

    i = 0;
    while(fgets(line, MAX_LINE_LENGTH, fptr) != NULL){
      if (line[0] != '#'){
        err = sscanf(line, "%"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM,
                          &DMParticleMass[i], 
                          &DMParticlePosition[0][i], &DMParticlePosition[1][i], &DMParticlePosition[2][i],
                          &DMParticleVelocity[0][i], &DMParticleVelocity[1][i], &DMParticleVelocity[2][i]);
        i++;
      }
    }
    fclose(fptr);
    int NumberOfDMParticles = i;

    Temp = &TopGrid;
    while (Temp != NULL){
      if(Temp->GridData->GalaxySimulationInitializeParticles(NumberOfDMParticles,
                                                               DMParticleMass, DMParticlePosition,
                                                               DMParticleVelocity) == FAIL){
        fprintf(stderr, "Error in grid->GalaxySimulationInitializeParticles.\n");
        return FAIL;
      } // end if particle initialize

      Temp = Temp->NextGridThisLevel;
    }

    /* clean up */
    delete [] DMParticleMass;

    for(int dim = 0; dim <MAX_DIMENSION; dim++){
      delete [] DMParticlePosition[dim];
      delete [] DMParticleVelocity[dim];
    }

    int LocalNumberOfParticles = 0;
    Temp = &TopGrid;
    /* figure out particle count across grids and processors */
    while (Temp != NULL){
      LocalNumberOfParticles = 0;
      LocalNumberOfParticles = Temp->GridData->ReturnNumberOfParticles();

#ifdef USE_MPI
      CommunicationAllReduceValues(&LocalNumberOfParticles, 1, MPI_SUM);
#endif /* USE_MPI */

      Temp->GridData->SetNumberOfParticles(LocalNumberOfParticles);

      LocalNumberOfParticles = Temp->GridData->ReturnNumberOfParticles();

      Temp = Temp->NextGridThisLevel;

    } // done communicating particle count

    Temp = &TopGrid;
    PINT ParticleCount = 0;
    RecursivelySetParticleCount(Temp, &ParticleCount);

    MetaData.NumberOfParticles = ParticleCount;

  } // end if using particles check

//  if (GalaxySimulationInitialStellarDist){
    if(FALSE){ // NOT WORKING NOT SURE WHY - causes seg fault
      // initialize star properties
      IndividualStarProperties_Initialize(MetaData);
      IndividualStarRadiationProperties_Initialize();
      InitializeStellarYields(MetaData.Time);

      Temp = &TopGrid;
      int TotalNumberOfNewParticles = 0;
      while (Temp != NULL){
        int MaximumNumberOfNewParticles = 10;
        int NumberOfNewParticles = 0;
        //TopGrid.GridData->AllocateNewParticles(MaximumNumberOfNewParticles);

        if(Temp->GridData->GalaxySimulationInitialStars(&MaximumNumberOfNewParticles, &NumberOfNewParticles) == FAIL){
          fprintf(stderr, "Error in grid->GalaxySimulationInitialStars.\n");
          return FAIL;
        } // end if particle initialize

        Temp = Temp->NextGridThisLevel;
        TotalNumberOfNewParticles += NumberOfNewParticles;
      }

//      if (TotalNumberOfNewParticles > 0) {
//        TogGrid->NumberOfParticles = NumberOfNewParticles;
 //     } else{
//        ENZO_FAIL("Was not able to deposit stars in chemical evolution test\n");
//      } // end: if (NumberOfNewParticles > 0)

    int LocalNumberOfParticles = 0;
    Temp = &TopGrid;
    /* figure out particle count across grids and processors */
    while (Temp != NULL){
      LocalNumberOfParticles = 0;
      LocalNumberOfParticles = Temp->GridData->ReturnNumberOfParticles();

#ifdef USE_MPI
      CommunicationAllReduceValues(&LocalNumberOfParticles, 1, MPI_SUM);
#endif /* USE_MPI */

      Temp->GridData->SetNumberOfParticles(LocalNumberOfParticles);

      LocalNumberOfParticles = Temp->GridData->ReturnNumberOfParticles();

      Temp = Temp->NextGridThisLevel;

    } // done communicating particle count

    Temp = &TopGrid;
    PINT ParticleCount = 0;
    RecursivelySetParticleCount(Temp, &ParticleCount);

    MetaData.NumberOfParticles = ParticleCount;


  } // end initial stellar dist

    /* Convert minimum initial overdensity for refinement to mass
       (unless MinimumMass itself was actually set). */

  if (MinimumMassForRefinement[0] == FLOAT_UNDEFINED) {
    MinimumMassForRefinement[0] = MinimumOverDensityForRefinement[0];
    for (int dim = 0; dim < MetaData.TopGridRank; dim++)
      MinimumMassForRefinement[0] *=(DomainRightEdge[dim]-DomainLeftEdge[dim])/
        float(MetaData.TopGridDims[dim]);
  }

  /* If requested, refine the grid to the desired level. */
  LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];
  if (GalaxySimulationRefineAtStart) {

    /* Declare, initialize and fill out the LevelArray. */

  //  LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];
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

	if (Temp->GridData->GalaxySimulationInitializeGrid(GalaxySimulationDiskRadius,
						       GalaxySimulationGalaxyMass, 
						       GalaxySimulationGasMass,
						       GalaxySimulationDiskPosition, 
						       GalaxySimulationDiskScaleHeightz,
						       GalaxySimulationDiskScaleHeightR,
						       GalaxySimulationTruncationRadius, 
						       GalaxySimulationDarkMatterConcentrationParameter,
						       GalaxySimulationDiskTemperature, 
						       GalaxySimulationInitialTemperature,
                                                       GalaxySimulationInitialDiskMetallicity,
                                                       GalaxySimulationInitialHaloMetallicity,
						       GalaxySimulationUniformDensity,
						       GalaxySimulationGasHalo,
						       GalaxySimulationGasHaloScaleRadius,
						       GalaxySimulationGasHaloDensity,
						       GalaxySimulationAngularMomentum,
						       GalaxySimulationUniformVelocity,
						       GalaxySimulationUseMetallicityField,
						       GalaxySimulationInflowTime,
						       GalaxySimulationInflowDensity,level,
						       GalaxySimulationCR,
                                                       GalaxySimulationUseDensityPerturbation,
                                                       GalaxySimulationPerturbationFraction,
                                                       GalaxySimulationSMAUGIC)
	      == FAIL) {
	    ENZO_FAIL("Error in GalaxySimulationInitialize[Sub]Grid.");
	}// end subgrid if

	Temp = Temp->NextGridThisLevel;
      }
    } // end: loop over levels

    /* Loop back from the bottom, restoring the consistency among levels. */

    for (level = MaximumRefinementLevel; level > 0; level--) {
      LevelHierarchyEntry *Temp = LevelArray[level];
      while (Temp != NULL) {
	if (Temp->GridData->ProjectSolutionToParentGrid(
				   *LevelArray[level-1]->GridData) == FAIL) {
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
	
    /* Set Global Variables for RPS Wind (see ExternalBoundary_SetGalaxySimulationBoundary.C)*/

    GalaxySimulationRPSWindDelay += TopGrid.GridData->ReturnTime();
    GalaxySimulationRPSWindTotalEnergy = InflowValue[1]; 	
    GalaxySimulationPreWindDensity     = GalaxySimulationUniformDensity/DensityUnits;
    GalaxySimulationPreWindTotalEnergy = GalaxySimulationInitialTemperature/TemperatureUnits/((Gamma-1.0)*0.6); 
    GalaxySimulationPreWindVelocity[0] = 0.0;
    GalaxySimulationPreWindVelocity[1] = 0.0;
    GalaxySimulationPreWindVelocity[2] = 0.0;
  }

  if (IndividualStarICSupernovaRate > 0){ // save galaxy properties for SN driving if used
    if (IndividualStarICSupernovaR < 0) 
      IndividualStarICSupernovaR = GalaxySimulationDiskScaleHeightR; // in Mpc_cm
    if (IndividualStarICSupernovaZ < 0)
      IndividualStarICSupernovaZ = GalaxySimulationDiskScaleHeightz; // in Mpc_cm
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
 if(CRModel)
   DataLabel[count++] = CRName;

 if(WritePotential)
   DataLabel[count++] = GravPotentialName;

 if (MultiSpecies) {
   DataLabel[count++] = ElectronName;
   DataLabel[count++] = HIName;
   DataLabel[count++] = HIIName;
   DataLabel[count++] = HeIName;
   DataLabel[count++] = HeIIName;
   DataLabel[count++] = HeIIIName;

   if (MultiSpecies > 1){
     DataLabel[count++] = HMName;
     DataLabel[count++] = H2IName;
     DataLabel[count++] = H2IIName;
   }

   if (MultiSpecies > 2){
     DataLabel[count++] = DIName;
     DataLabel[count++] = DIIName;
     DataLabel[count++] = HDIName;
   }
 }

 if (GalaxySimulationUseMetallicityField)
   DataLabel[count++] = MetalName;

 if (TestProblemData.UseMetallicityField){
   DataLabel[count++] = MetalName;
 }

 if (StarMakerTypeIaSNe)
   DataLabel[count++] = MetalIaName;

 if (MultiMetals == 2){

   for(int i =0; i < StellarYieldsNumberOfSpecies; i ++){
     if(StellarYieldsAtomicNumbers[i] > 2){
       DataLabel[count++] = ChemicalSpeciesBaryonFieldLabel(StellarYieldsAtomicNumbers[i]);
     }
   } // yields loop

   if (IndividualStarTrackAGBMetalDensity){
     DataLabel[i++] = AGBMetalName;
   }

   if (IndividualStarPopIIIFormation){
     DataLabel[i++] = PopIIIMetalName;
     DataLabel[i++] = PopIIIPISNeMetalName;
   }

   if (IndividualStarTrackSNMetalDensity){
     DataLabel[i++] = SNIaMetalName;

     if (IndividualStarSNIaModel == 2){
      DataLabel[i++] = ExtraMetalName0;
      DataLabel[i++] = ExtraMetalName1;
      DataLabel[i++] = ExtraMetalName2;
     }

     DataLabel[i++] = SNIIMetalName;
   }

   if (IndividualStarRProcessModel){
     DataLabel[i++] = RProcMetalName;
   }

 }

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
   fprintf(Outfptr, "GalaxySimulationUniformCR = %"GOUTSYM"\n",
     GalaxySimulationUniformCR);
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
   fprintf(Outfptr, "GalaxySimulationInflowTime = %"GOUTSYM"\n",
	   GalaxySimulationInflowTime);
   fprintf(Outfptr, "GalaxySimulationInflowDensity = %"GOUTSYM"\n",
	   GalaxySimulationInflowDensity);

   /* Output chemical abundance IC's */

   // galaxy chemistry
   fprintf(Outfptr, "GalaxySimulationHydrogenFractionByMass = %"GOUTSYM"\n",
           TestProblemData.InnerHydrogenFractionByMass);
   fprintf(Outfptr, "GalaxySimulationDeuteriumToHydrogenRatio = %"GOUTSYM"\n",
           TestProblemData.InnerDeuteriumToHydrogenRatio);
   fprintf(Outfptr, "GalaxySimulationHIIFraction = %"GOUTSYM"\n",
           TestProblemData.HII_Fraction_Inner);
   fprintf(Outfptr, "GalaxySimulationHeIIFraction = %"GOUTSYM"\n",
           TestProblemData.HeII_Fraction_Inner);
   fprintf(Outfptr, "GalaxySimulationHeIIIFraction = %"GOUTSYM"\n",
           TestProblemData.HeIII_Fraction_Inner);
   fprintf(Outfptr, "GalaxySimulationHMFraction = %"GOUTSYM"\n",
           TestProblemData.HM_Fraction_Inner);
   fprintf(Outfptr, "GalaxySimulationH2IFraction = %"GOUTSYM"\n",
           TestProblemData.H2I_Fraction_Inner);
   fprintf(Outfptr, "GalaxySimulationH2IIFraction = %"GOUTSYM"\n",
           TestProblemData.H2II_Fraction_Inner);

   // halo chemistry
   fprintf(Outfptr, "GalaxySimulationHydrogenFractionByMassHalo = %"GOUTSYM"\n",
              TestProblemData.HydrogenFractionByMass);
   fprintf(Outfptr, "GalaxySimulationDeuteriumToHydrogenRatioHalo = %"GOUTSYM"\n",
              TestProblemData.DeuteriumToHydrogenRatio);
   fprintf(Outfptr, "GalaxySimulationHIIFractionHalo = %"GOUTSYM"\n",
              TestProblemData.HII_Fraction);
   fprintf(Outfptr, "GalaxySimulationHeIIFractionHalo = %"GOUTSYM"\n",
              TestProblemData.HeII_Fraction);
   fprintf(Outfptr, "GalaxySimulationHeIIIFractionHalo = %"GOUTSYM"\n",
              TestProblemData.HeIII_Fraction);
   fprintf(Outfptr, "GalaxySimulationHMFraction = %"GOUTSYM"\n",
              TestProblemData.HM_Fraction);
   fprintf(Outfptr, "GalaxySimulationH2IFractionHalo = %"GOUTSYM"\n",
              TestProblemData.H2I_Fraction);
   fprintf(Outfptr, "GalaxySimulationH2IIFractionHalo = %"GOUTSYM"\n",
              TestProblemData.H2II_Fraction);

   // galaxy chemical tracers
   fprintf(Outfptr, "GalaxySimulationInitialCIFraction = %"GOUTSYM"\n",
           TestProblemData.CI_Fraction);
   fprintf(Outfptr, "GalaxySimulationInitialNIFraction = %"GOUTSYM"\n",
           TestProblemData.NI_Fraction);
   fprintf(Outfptr, "GalaxySimulationInitialOIFraction = %"GOUTSYM"\n",
           TestProblemData.OI_Fraction);
   fprintf(Outfptr, "GalaxySimulationInitialMgIFraction = %"GOUTSYM"\n",
           TestProblemData.MgI_Fraction);
   fprintf(Outfptr, "GalaxySimulationInitialSiIFraction = %"GOUTSYM"\n",
           TestProblemData.SiI_Fraction);
   fprintf(Outfptr, "GalaxySimulationInitialFeIFraction = %"GOUTSYM"\n",
           TestProblemData.FeI_Fraction);
   fprintf(Outfptr, "GalaxySimulationInitialYIFraction = %"GOUTSYM"\n",
           TestProblemData.YI_Fraction);
   fprintf(Outfptr, "GalaxySimulationInitialBaIFraction = %"GOUTSYM"\n",
           TestProblemData.BaI_Fraction);
   fprintf(Outfptr, "GalaxySimulationInitialLaIFraction = %"GOUTSYM"\n",
           TestProblemData.LaI_Fraction);
   fprintf(Outfptr, "GalaxySimulationInitialEuIFraction = %"GOUTSYM"\n",
           TestProblemData.EuI_Fraction);

   // halo chemical tracers
   fprintf(Outfptr, "GalaxySimulationInitialCIFractionHalo = %"GOUTSYM"\n",
           TestProblemData.CI_Fraction_2);
   fprintf(Outfptr, "GalaxySimulationInitialNIFractionHalo = %"GOUTSYM"\n",
           TestProblemData.NI_Fraction_2);
   fprintf(Outfptr, "GalaxySimulationInitialOIFractionHalo = %"GOUTSYM"\n",
           TestProblemData.OI_Fraction_2);
   fprintf(Outfptr, "GalaxySimulationInitialMgIFractionHalo = %"GOUTSYM"\n",
           TestProblemData.MgI_Fraction_2);
   fprintf(Outfptr, "GalaxySimulationInitialSiIFractionHalo = %"GOUTSYM"\n",
           TestProblemData.SiI_Fraction_2);
   fprintf(Outfptr, "GalaxySimulationInitialFeIFractionHalo = %"GOUTSYM"\n",
           TestProblemData.FeI_Fraction_2);
   fprintf(Outfptr, "GalaxySimulationInitialYIFractionHalo = %"GOUTSYM"\n",
           TestProblemData.YI_Fraction_2);
   fprintf(Outfptr, "GalaxySimulationInitialBaIFractionHalo = %"GOUTSYM"\n",
           TestProblemData.BaI_Fraction_2);
   fprintf(Outfptr, "GalaxySimulationInitialLaIFractionHalo = %"GOUTSYM"\n",
           TestProblemData.LaI_Fraction_2);
   fprintf(Outfptr, "GalaxySimulationInitialEuIFractionHalo = %"GOUTSYM"\n",
           TestProblemData.EuI_Fraction_2);

   fprintf(Outfptr, "GalaxySimulationInitialSpeciesFractionsDisk = ");
   WriteListOfFloats(Outfptr, MAX_STELLAR_YIELDS, TestProblemData.ChemicalTracerSpecies_Fractions);

   fprintf(Outfptr, "GalaxySimulationInitialSpeciesFractionsHalo = ");
   WriteListOfFloats(Outfptr, MAX_STELLAR_YIELDS, TestProblemData.ChemicalTracerSpecies_Fractions_2);

   fprintf(Outfptr, "TestProblemUseMetallicityField = %"ISYM"\n",
           TestProblemData.UseMetallicityField);

   fprintf(Outfptr, "GalaxySimulationInitialDiskMetallicity = %"GOUTSYM"\n",
           GalaxySimulationInitialDiskMetallicity);

   fprintf(Outfptr, "GalaxySimulationInitialHaloMetallicity = %"GOUTSYM"\n",
           GalaxySimulationInitialHaloMetallicity);

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


