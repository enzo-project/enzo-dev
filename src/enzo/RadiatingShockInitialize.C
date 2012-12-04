/***********************************************************************
/
/  INITIALIZE RADIATING SEDOV BLAST WAVE
/
/  written by: Brian O'Shea
/  date:       August 2007
/  modified1:  
/
/  PURPOSE:
/
/   REFERENCE: Self-similar solution: L.I. Sedov (1946);
/              see also: Sedov (1959), Similarity and Dimensional Methods
/              in Mechanics, pp. 210, 219, 228;
/              see also: Landau & Lifshitz, Fluid Dynamics, Sect. 99
/              "The Propagation of Strong Shock Waves" (1959).
/              Experiments, terrestrial/numerical: Taylor (1941, 1949).
/
/   Two dimensional parameters: explosion energy E and ambient density rho_1
/   Two independent variables: radius r, time t
/   One dimensionless combination: r*(rho_1/E/t^2)^(1/5)
/
/
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
// This routine intializes a new simulation based on the parameter file.
//
 
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
#define DEFINE_STORAGE
#include "RadiatingShockGlobalData.h"
#undef DEFINE_STORAGE

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);
 
int RadiatingShockInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
			 TopGridData &MetaData)
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
  char *COIName   = "COI_Density";
  char *CIName    = "CI_Density";
  char *CIIName   = "CII_Density";
  char *OIName    = "OI_Density";
  char *OIIName   = "OII_Density";
  char *SiIName   = "SiI_Density";
  char *SiIIName  = "SiII_Density";
  char *SiIIIName = "SiIII_Density";
  char *CHIName   = "CHI_Density";
  char *CH2IName  = "CH2I_Density";
  char *CH3IIName = "CH3II_Density";
  char *C2IName   = "C2I_Density";
  char *HCOIIName = "HCOII_Density";
  char *OHIName   = "OHI_Density";
  char *H2OIName  = "H2OI_Density";
  char *O2IName   = "O2I_Density";
  char *MetalName = "Metal_Density";
  char *ExtraNames[2] = {"Z_Field1", "Z_Field2"};

  /* parameter declarations */
 
  FLOAT RadiatingShockSubgridLeft, RadiatingShockSubgridRight;
  FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
  FLOAT RadiatingShockCenterPosition[MAX_DIMENSION];
  float ZeroBField[3] = {0.0, 0.0, 0.0};

  /* local declarations */
 
  char line[MAX_LINE_LENGTH];
  int  i, j, k, dim, ret, NumberOfSubgridZones[MAX_DIMENSION],
    SubgridDims[MAX_DIMENSION];
 
  /* make sure it is 2D or 3D */
 
  if (MetaData.TopGridRank < 2 || MetaData.TopGridRank > 3) {
    ENZO_VFAIL("Cannot do RadiatingShock in %"ISYM" dimension(s)\n", MetaData.TopGridRank)
  }
 
  /* There are many parameters:  geometry (cylindrical or spherical symmetry),
                                 gamma,
				 E,
				 rho_1.
     Set their default values */
 
  int RadiatingShockRandomSeedInitialize = 0;

  for(i=0; i<MAX_DIMENSION; i++)
    RadiatingShockCenterPosition[i] = 0.5;  // right in the middle of the box

  float Pi                      = 3.14159;
  float RadiatingShockVelocity[3]   = {0.0, 0.0, 0.0};   // gas initally at rest
  float RadiatingShockPressure      = 1e-5;
  float RadiatingShockInnerDensity             = 1.0;
  float RadiatingShockOuterDensity             = 1.0;
  float RadiatingShockEnergy        = 1.0;
  int RadiatingShockUseDensityFluctuations = 0;
  int RadiatingShockRandomSeed = 123456789;
  float RadiatingShockDensityFluctuationLevel = 0.1;
  int RadiatingShockInitializeWithKE = 0;
  int RadiatingShockUseSedovProfile = 0;
  FLOAT RadiatingShockSedovBlastRadius = 0.05;
  float RadiatingShockKineticEnergyFraction = 0.0;

  FLOAT RadiatingShockSpreadOverNumZones = 3.5;
  float dx = (DomainRightEdge[0] - DomainLeftEdge[0])/
    MetaData.TopGridDims[0];
 
  /* set no subgrids by default. */
 
  RadiatingShockSubgridLeft         = 0.0;    // start of subgrid(s)
  RadiatingShockSubgridRight        = 0.0;    // end of subgrid(s)


  TestProblemData.MultiSpecies = MultiSpecies;  // set this from global data (kind of a hack, but necessary)
  TestProblemData.GloverChemistryModel = GloverChemistryModel; // set this from global data (kind of a hack, but necessary)
 
  /* read input from file */
 
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
 
    ret = 0;
 
    /* read parameters specifically for radiating shock problem*/
 
    ret += sscanf(line, "RadiatingShockInnerDensity  = %"FSYM, &RadiatingShockInnerDensity); // density inside heated region
    ret += sscanf(line, "RadiatingShockOuterDensity  = %"FSYM, &RadiatingShockOuterDensity); // ambient density
    ret += sscanf(line, "RadiatingShockPressure = %"FSYM, &RadiatingShockPressure);  // ambient pressure
    ret += sscanf(line, "RadiatingShockEnergy   = %"FSYM, &RadiatingShockEnergy);  // supernova explosion energy
    ret += sscanf(line, "RadiatingShockSubgridLeft = %"PSYM,
		        &RadiatingShockSubgridLeft);
    ret += sscanf(line, "RadiatingShockSubgridRight = %"PSYM,
		        &RadiatingShockSubgridRight);
    ret += sscanf(line, "RadiatingShockUseDensityFluctuations   = %"ISYM, &RadiatingShockUseDensityFluctuations);
    ret += sscanf(line, "RadiatingShockRandomSeed   = %"ISYM, &RadiatingShockRandomSeed);
    ret += sscanf(line, "RadiatingShockDensityFluctuationLevel   = %"FSYM, &RadiatingShockDensityFluctuationLevel);
    ret += sscanf(line, "RadiatingShockInitializeWithKE = %"ISYM, &RadiatingShockInitializeWithKE);
    ret += sscanf(line, "RadiatingShockUseSedovProfile = %"ISYM, &RadiatingShockUseSedovProfile);
    ret += sscanf(line, "RadiatingShockSedovBlastRadius = %"PSYM, &RadiatingShockSedovBlastRadius);
    ret += sscanf(line, "RadiatingShockUseSedovProfile = %"ISYM, &RadiatingShockUseSedovProfile);

    ret += sscanf(line, "RadiatingShockKineticEnergyFraction = %"FSYM, &RadiatingShockKineticEnergyFraction);

    ret += sscanf(line, "RadiatingShockCenterPosition = %"PSYM" %"PSYM" %"PSYM,
		  RadiatingShockCenterPosition, RadiatingShockCenterPosition+1,
		  RadiatingShockCenterPosition+2);

    ret += sscanf(line, "RadiatingShockSpreadOverNumZones  = %"PSYM, &RadiatingShockSpreadOverNumZones);


    /* read in more general test parameters to set species, turn on color fields, etc. */
    ret += sscanf(line, "TestProblemHydrogenFractionByMass = %"FSYM, &TestProblemData.HydrogenFractionByMass);
    ret += sscanf(line, "TestProblemDeuteriumToHydrogenRatio = %"FSYM, &TestProblemData.DeuteriumToHydrogenRatio);

    ret += sscanf(line, "TestProblemInitialHIFractionInner  = %"FSYM, &TestProblemData.HI_Fraction_Inner);
    ret += sscanf(line, "TestProblemInitialHIIFractionInner  = %"FSYM, &TestProblemData.HII_Fraction_Inner);
    ret += sscanf(line, "TestProblemInitialHeIFractionInner  = %"FSYM, &TestProblemData.HeI_Fraction_Inner);
    ret += sscanf(line, "TestProblemInitialHeIIFractionInner  = %"FSYM, &TestProblemData.HeII_Fraction_Inner);
    ret += sscanf(line, "TestProblemInitialHeIIIFractionInner  = %"FSYM, &TestProblemData.HeIII_Fraction_Inner);
    ret += sscanf(line, "TestProblemInitialHMFractionInner  = %"FSYM, &TestProblemData.HM_Fraction_Inner);
    ret += sscanf(line, "TestProblemInitialH2IFractionInner  = %"FSYM, &TestProblemData.H2I_Fraction_Inner);
    ret += sscanf(line, "TestProblemInitialH2IIFractionInner  = %"FSYM, &TestProblemData.H2II_Fraction_Inner);

    ret += sscanf(line, "TestProblemInitialDIFractionInner  = %"FSYM, &TestProblemData.DI_Fraction_Inner);
    ret += sscanf(line, "TestProblemInitialDIIFractionInner  = %"FSYM, &TestProblemData.DII_Fraction_Inner);
    ret += sscanf(line, "TestProblemInitialHDIFractionInner  = %"FSYM, &TestProblemData.HDI_Fraction_Inner);

    ret += sscanf(line, "TestProblemInitialCOIFractionInner  = %"FSYM, &TestProblemData.COI_Fraction_Inner);
    ret += sscanf(line, "TestProblemInitialCIFractionInner  = %"FSYM, &TestProblemData.CI_Fraction_Inner);
    ret += sscanf(line, "TestProblemInitialCIIFractionInner  = %"FSYM, &TestProblemData.CII_Fraction_Inner);
    ret += sscanf(line, "TestProblemInitialOIFractionInner  = %"FSYM, &TestProblemData.OI_Fraction_Inner);
    ret += sscanf(line, "TestProblemInitialOIIFractionInner  = %"FSYM, &TestProblemData.OII_Fraction_Inner);
    ret += sscanf(line, "TestProblemInitialSiIFractionInner  = %"FSYM, &TestProblemData.SiI_Fraction_Inner);
    ret += sscanf(line, "TestProblemInitialSiIIFractionInner  = %"FSYM, &TestProblemData.SiII_Fraction_Inner);
    ret += sscanf(line, "TestProblemInitialSiIIIFractionInner  = %"FSYM, &TestProblemData.SiIII_Fraction_Inner);
    ret += sscanf(line, "TestProblemInitialCHIFractionInner  = %"FSYM, &TestProblemData.CHI_Fraction_Inner);
    ret += sscanf(line, "TestProblemInitialCH2IFractionInner  = %"FSYM, &TestProblemData.CH2I_Fraction_Inner);
    ret += sscanf(line, "TestProblemInitialCH3IIFractionInner  = %"FSYM, &TestProblemData.CH3II_Fraction_Inner);
    ret += sscanf(line, "TestProblemInitialC2IFractionInner  = %"FSYM, &TestProblemData.C2I_Fraction_Inner);
    ret += sscanf(line, "TestProblemInitialHCOIIFractionInner  = %"FSYM, &TestProblemData.HCOII_Fraction_Inner);
    ret += sscanf(line, "TestProblemInitialOHIFractionInner  = %"FSYM, &TestProblemData.OHI_Fraction_Inner);
    ret += sscanf(line, "TestProblemInitialH2OIFractionInner  = %"FSYM, &TestProblemData.H2OI_Fraction_Inner);
    ret += sscanf(line, "TestProblemInitialO2IFractionInner  = %"FSYM, &TestProblemData.O2I_Fraction_Inner);

    ret += sscanf(line, "TestProblemInitialHIFraction  = %"FSYM, &TestProblemData.HI_Fraction);
    ret += sscanf(line, "TestProblemInitialHIIFraction  = %"FSYM, &TestProblemData.HII_Fraction);
    ret += sscanf(line, "TestProblemInitialHeIFraction  = %"FSYM, &TestProblemData.HeI_Fraction);
    ret += sscanf(line, "TestProblemInitialHeIIFraction  = %"FSYM, &TestProblemData.HeII_Fraction);
    ret += sscanf(line, "TestProblemInitialHeIIIFraction  = %"FSYM, &TestProblemData.HeIII_Fraction);
    ret += sscanf(line, "TestProblemInitialHMFraction  = %"FSYM, &TestProblemData.HM_Fraction);
    ret += sscanf(line, "TestProblemInitialH2IFraction  = %"FSYM, &TestProblemData.H2I_Fraction);
    ret += sscanf(line, "TestProblemInitialH2IIFraction  = %"FSYM, &TestProblemData.H2II_Fraction);

    ret += sscanf(line, "TestProblemInitialDIFraction  = %"FSYM, &TestProblemData.DI_Fraction);
    ret += sscanf(line, "TestProblemInitialDIIFraction  = %"FSYM, &TestProblemData.DII_Fraction);
    ret += sscanf(line, "TestProblemInitialHDIFraction  = %"FSYM, &TestProblemData.HDI_Fraction);

    ret += sscanf(line, "TestProblemInitialCOIFraction  = %"FSYM, &TestProblemData.COI_Fraction);
    ret += sscanf(line, "TestProblemInitialCIFraction  = %"FSYM, &TestProblemData.CI_Fraction);
    ret += sscanf(line, "TestProblemInitialCIIFraction  = %"FSYM, &TestProblemData.CII_Fraction);
    ret += sscanf(line, "TestProblemInitialOIFraction  = %"FSYM, &TestProblemData.OI_Fraction);
    ret += sscanf(line, "TestProblemInitialOIIFraction  = %"FSYM, &TestProblemData.OII_Fraction);
    ret += sscanf(line, "TestProblemInitialSiIFraction  = %"FSYM, &TestProblemData.SiI_Fraction);
    ret += sscanf(line, "TestProblemInitialSiIIFraction  = %"FSYM, &TestProblemData.SiII_Fraction);
    ret += sscanf(line, "TestProblemInitialSiIIIFraction  = %"FSYM, &TestProblemData.SiIII_Fraction);
    ret += sscanf(line, "TestProblemInitialCHIFraction  = %"FSYM, &TestProblemData.CHI_Fraction);
    ret += sscanf(line, "TestProblemInitialCH2IFraction  = %"FSYM, &TestProblemData.CH2I_Fraction);
    ret += sscanf(line, "TestProblemInitialCH3IIFraction  = %"FSYM, &TestProblemData.CH3II_Fraction);
    ret += sscanf(line, "TestProblemInitialC2IFraction  = %"FSYM, &TestProblemData.C2I_Fraction);
    ret += sscanf(line, "TestProblemInitialHCOIIFraction  = %"FSYM, &TestProblemData.HCOII_Fraction);
    ret += sscanf(line, "TestProblemInitialOHIFraction  = %"FSYM, &TestProblemData.OHI_Fraction);
    ret += sscanf(line, "TestProblemInitialH2OIFraction  = %"FSYM, &TestProblemData.H2OI_Fraction);
    ret += sscanf(line, "TestProblemInitialO2IFraction  = %"FSYM, &TestProblemData.O2I_Fraction);

    ret += sscanf(line, "TestProblemUseMetallicityField  = %"ISYM, &TestProblemData.UseMetallicityField);
    ret += sscanf(line, "TestProblemInitialMetallicityFraction  = %"FSYM, &TestProblemData.MetallicityField_Fraction);

    ret += sscanf(line, "TestProblemUseMassInjection  = %"ISYM, &TestProblemData.UseMassInjection);
    ret += sscanf(line, "TestProblemInitialHydrogenMass  = %"FSYM, &TestProblemData.InitialHydrogenMass);
    ret += sscanf(line, "TestProblemInitialDeuteriumMass  = %"FSYM, &TestProblemData.InitialDeuteriumMass);
    ret += sscanf(line, "TestProblemInitialHeliumMass  = %"FSYM, &TestProblemData.InitialHeliumMass);
    ret += sscanf(line, "TestProblemInitialMetalMass  = %"FSYM, &TestProblemData.InitialMetalMass);

    ret += sscanf(line, "TestProblemMultiMetals  = %"ISYM, &TestProblemData.MultiMetals);
    ret += sscanf(line, "TestProblemInitialMultiMetalsField1Fraction  = %"FSYM, &TestProblemData.MultiMetalsField1_Fraction);
    ret += sscanf(line, "TestProblemInitialMultiMetalsField2Fraction  = %"FSYM, &TestProblemData.MultiMetalsField2_Fraction);

    /* if the line is suspicious, issue a warning */
 
    if (ret == 0 && strstr(line, "=") && (strstr(line, "RadiatingShock") || strstr(line, "TestProblem")) &&
	line[0] != '#' && MyProcessorNumber == ROOT_PROCESSOR)
      fprintf(stderr,
	 "*** warning: the following parameter line was not interpreted:\n%s\n",
	      line);
 
  } // end input from parameter file
 
  /* set number of zones on the finest level to resolve the initial explosion */
  FLOAT dr = RadiatingShockSpreadOverNumZones*dx*POW(RefineBy,-MaximumRefinementLevel);
 
  /* compute p_2 as a function of explosion energy E, initial explosion
     radius dr, and gamma.
 
     2D:  p_2 = (gamma-1)*E/(pi*r^2)*rho_in
     3D:  p_2 = (gamma-1)*E/(4/3*pi*r^3)*rho_in
          rho_2 = 1 = rho_1
 
     If p_2 (inner) is too high, the code crashes in euler.src (dnu<0) after
     a few time steps.
     In 2D it is stable for RadiatingShockEnergy <= 0.025 (tested with uniform grid 200^2).
  */
 
  float DensityUnits=1, LengthUnits=1, TemperatureUnits=1, TimeUnits=1,
    VelocityUnits=1;
  double MassUnits=1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, &MassUnits, 0.0) == FAIL) {
        ENZO_FAIL("Error in GetUnits.");
  }

  double RadiatingShockEnergyDouble;


  // calculate some values if we aren't using the Sedov profile.  If we are, this is all
  // calculated inside Grid::RadiatingShockInitialize
  if(!RadiatingShockUseSedovProfile){
    // input is in units of 1 FOE, this converts to CGS
    RadiatingShockEnergyDouble = double(RadiatingShockEnergy) * 1.0e+51;

    // converts this to ergs per gram in code units
    RadiatingShockEnergyDouble /= (double(DensityUnits) * POW(double(LengthUnits),5.0) * POW(double(TimeUnits),-2.0) );
  
    RadiatingShockEnergy = float(RadiatingShockEnergyDouble);

  } else {
    printf("using sedov profile: using RadiatingShockEnergy alone!\n");
 
  }

  if(debug)
    printf("RadiatingShockInitialize:  RadiatingShockEnergy is %e in code units\n",RadiatingShockEnergy);

  float numberInjectionCells = 0;
  double InjectionMass2Density_scaleFactor;

  // If injecting a mass of gas into the center, 
  // effective number of cells is simply area/volume of circle/sphere
  // with r = RadiatingShockSpreadOverNumZones.
  if (TestProblemData.UseMassInjection) {
    // Make sure total mass is not zero.
    if ((TestProblemData.InitialHydrogenMass <= 0.0) && (TestProblemData.InitialHeliumMass <= 0.0)) {
      ENZO_FAIL("Hydrogen and helium mass cannot both be zero.  That would be zero mass in the center.\n");
    }

    // 2D
    if(MetaData.TopGridRank==2){
      numberInjectionCells = Pi * pow(RadiatingShockSpreadOverNumZones,2);
    }
    // 3D
    else {
      numberInjectionCells = (4./3.) * Pi * pow(RadiatingShockSpreadOverNumZones,3);
    }

    InjectionMass2Density_scaleFactor = 1.989e33 / 
      (MassUnits * numberInjectionCells * 
       pow((dx*POW(RefineBy,-MaximumRefinementLevel)),3));

    // ignore D mass
    RadiatingShockInnerDensity = (TestProblemData.InitialHydrogenMass + 
				  TestProblemData.InitialHeliumMass) * 
      float(InjectionMass2Density_scaleFactor);
    fprintf(stderr,"Setting inner density to %e.\n",RadiatingShockInnerDensity);

    // Adjust H mass fraction
    TestProblemData.InnerHydrogenFractionByMass = TestProblemData.InitialHydrogenMass /
      (TestProblemData.InitialHydrogenMass + TestProblemData.InitialHeliumMass);
    fprintf(stderr,"Inner H mass fraction is %.2f.\n",TestProblemData.InnerHydrogenFractionByMass);

    // Adjust D/H ratio
    TestProblemData.InnerDeuteriumToHydrogenRatio = TestProblemData.InitialDeuteriumMass /
      TestProblemData.InitialHydrogenMass;
    fprintf(stderr,"Inner D/H ratio is %.2e.\n",TestProblemData.InnerDeuteriumToHydrogenRatio);

    // Set metal fraction
    TestProblemData.MetallicityField_Fraction = TestProblemData.InitialMetalMass /
      (TestProblemData.InitialHydrogenMass + TestProblemData.InitialHeliumMass);
    fprintf(stderr,"Inner metal fraction is %.2e.\n",TestProblemData.MetallicityField_Fraction);
  }
  else {
    TestProblemData.InnerHydrogenFractionByMass = TestProblemData.HydrogenFractionByMass;
    TestProblemData.InnerDeuteriumToHydrogenRatio = TestProblemData.DeuteriumToHydrogenRatio;
  }

  // BWO: modified to include the idea that the inner region might have more gas
  float RadiatingShockInnerPressure = 3.0*(Gamma-1.0)*RadiatingShockEnergy*RadiatingShockInnerDensity/
                                  (MetaData.TopGridRank + 1.0)/
                                  POW(dr,MetaData.TopGridRank)/Pi;
 
  /* Check the self-similarity condition: p2/p1 >> (gamma+1)/(gamma-1). */
 
  float pjump = RadiatingShockInnerPressure/RadiatingShockPressure;
  if ( pjump < 10.*(Gamma+1)/(Gamma-1) )
    printf("SBI: WARNING! No self-similarity. Pressure jump %"FSYM".\n", pjump);
 
  /* Compute energies */

  // ambient gas internal energy
  RadiatingShockTotalEnergy = RadiatingShockPressure/((Gamma - 1.0)*RadiatingShockOuterDensity);

  if(debug)
    printf("ambient gas energy should be %e\n",RadiatingShockTotalEnergy);

  // heated gas internal energy
  RadiatingShockInnerTotalEnergy= RadiatingShockInnerPressure/((Gamma - 1.0)*
						   RadiatingShockInnerDensity);

  /* calculate kinetic energy quantities (if we're initializing with a sawtooth) 
     this is similar to what we do above, for 2D (cylindrical) and 3D(spherical) */
  float RadiatingShockRhoZero, RadiatingShockVelocityZero,  MassZero;

  RadiatingShockRhoZero = RadiatingShockVelocityZero =  MassZero = 0.0;

  /* we calculate all of this stuff when not using a Sedov blast profile.  If we ARE using
     Sedov values, ignore it because it doesn't matter. */
  if(RadiatingShockInitializeWithKE && !RadiatingShockUseSedovProfile){

    if(MetaData.TopGridRank==2){  // cylindrical supernova (2D)

      MassZero = Pi * POW(dr,2.0) * RadiatingShockInnerDensity; // Mzero in code units (actually mass per unit length)

      RadiatingShockRhoZero = MassZero / (2.0*Pi/3.0* POW(dr,2.0) );  // this is really density

      // units actually work out to velocity (assuming input energy is really energy per unit length)
      RadiatingShockVelocityZero = POW( 10.0/3.0 * (RadiatingShockEnergy*RadiatingShockKineticEnergyFraction) / MassZero, 0.5 );

    } else {  // spherical supernova (3D)

      MassZero = 1.3333 * Pi * POW(dr,3.0) * RadiatingShockInnerDensity;  // Mzero in code units

      RadiatingShockRhoZero = MassZero / Pi / POW(dr,3.0);

      fprintf(stderr,"SBI: RadiatingShockKineticEnergyFraction is %e\n",RadiatingShockKineticEnergyFraction);

      RadiatingShockVelocityZero = POW( 3.0 * (RadiatingShockEnergy*RadiatingShockKineticEnergyFraction) / MassZero, 0.5);

    }

  }

  /* set the periodic boundaries */

  for (dim = 0; dim < MetaData.TopGridRank; dim++) {
    MetaData.LeftFaceBoundaryCondition[dim]  = periodic;
    MetaData.RightFaceBoundaryCondition[dim] = periodic;
  }
 
  /* set up uniform top grid before setting up explosion */
 
  if (TopGrid.GridData->InitializeUniformGrid(RadiatingShockOuterDensity,
					      RadiatingShockTotalEnergy,
					      RadiatingShockTotalEnergy,
					      RadiatingShockVelocity,
                          ZeroBField) == FAIL) {
        ENZO_FAIL("Error in InitializeUniformGrid.");
  }
 
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
	nint((RadiatingShockSubgridRight - RadiatingShockSubgridLeft)/
	     ((DomainRightEdge[dim] - DomainLeftEdge[dim] )/
	      float(MetaData.TopGridDims[dim])))
        *int(POW(RefineBy, lev + 1));
 
    if (debug)
      printf("RadiatingShock:: Level[%"ISYM"]: NumberOfSubgridZones[0] = %"ISYM"\n", lev+1,
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
	LeftEdge[dim]    = RadiatingShockSubgridLeft;
	RightEdge[dim]   = RadiatingShockSubgridRight;
      }
 
      /* create a new subgrid and initialize it */
 
      Subgrid[lev]->GridData = new grid;
      Subgrid[lev]->GridData->InheritProperties(TopGrid.GridData);
      Subgrid[lev]->GridData->PrepareGrid(MetaData.TopGridRank, SubgridDims,
				     LeftEdge, RightEdge, 0);
      if (Subgrid[lev]->GridData->InitializeUniformGrid(RadiatingShockOuterDensity,
						   RadiatingShockTotalEnergy,
						   RadiatingShockTotalEnergy,
						   RadiatingShockVelocity,
                           ZeroBField) == FAIL) {
		ENZO_FAIL("Error in InitializeUniformGrid (subgrid).");
      }
 
      /* set up the initial explosion area on the finest resolution subgrid */
 
      if (lev == MaximumRefinementLevel - 1)
	if (Subgrid[lev]->GridData->RadiatingShockInitializeGrid(dr,
				    RadiatingShockInnerDensity,
				    RadiatingShockInnerTotalEnergy,
				    RadiatingShockUseDensityFluctuations,
				    RadiatingShockRandomSeed,
				    RadiatingShockDensityFluctuationLevel,
				    RadiatingShockInitializeWithKE,
				    RadiatingShockUseSedovProfile,
				    RadiatingShockSedovBlastRadius,
				    RadiatingShockEnergy,
				    RadiatingShockPressure,
				    RadiatingShockKineticEnergyFraction,
				    RadiatingShockRhoZero,
				    RadiatingShockVelocityZero,
				    RadiatingShockRandomSeedInitialize,
 				    RadiatingShockCenterPosition) 
	    == FAIL) {
	  	  ENZO_FAIL("Error in RadiatingShockInitialize[Sub]Grid.");
	}

      RadiatingShockRandomSeedInitialize = 1;  // random number generator is now seeded - don't do it for topgrid
    }
    else{
      printf("RadiatingShock: single grid start-up.\n");
    }
  }

  /* set up subgrids from level 1 to max refinement level -1 */
 
  for (lev = MaximumRefinementLevel - 1; lev > 0; lev--)
    if (Subgrid[lev]->GridData->ProjectSolutionToParentGrid(
				       *(Subgrid[lev-1]->GridData))
	== FAIL) {
            ENZO_FAIL("Error in ProjectSolutionToParentGrid.");
    }

  /* set up the root grid */
 
  if (MaximumRefinementLevel > 0) {
    if (Subgrid[0]->GridData->ProjectSolutionToParentGrid(*(TopGrid.GridData))
	== FAIL) {
            ENZO_FAIL("Error in ProjectSolutionToParentGrid.");
    }
  }
  else
    if (TopGrid.GridData->RadiatingShockInitializeGrid(dr,
				    RadiatingShockInnerDensity,
			            RadiatingShockInnerTotalEnergy,
				    RadiatingShockUseDensityFluctuations,
				    RadiatingShockRandomSeed,
				    RadiatingShockDensityFluctuationLevel,
				    RadiatingShockInitializeWithKE,
				    RadiatingShockUseSedovProfile,
				    RadiatingShockSedovBlastRadius,
				    RadiatingShockEnergy,
				    RadiatingShockPressure,
				    RadiatingShockKineticEnergyFraction,
				    RadiatingShockRhoZero,
				    RadiatingShockVelocityZero,
				    RadiatingShockRandomSeedInitialize,
 				    RadiatingShockCenterPosition ) == FAIL) {
            ENZO_FAIL("Error in RadiatingShockInitializeGrid.");
    }

  /* set up field names and units -- NOTE: these absolutely MUST be in 
     the same order that they are in Grid_InitializeUniformGrids.C, or 
     else you'll find out that data gets written into incorrectly-named
     fields.  Just FYI. */

  i = 0;
  DataLabel[i++] = DensName;
  DataLabel[i++] = TEName;
  if(DualEnergyFormalism)
    DataLabel[i++] = GEName;
  DataLabel[i++] = Vel1Name;

  if(MetaData.TopGridRank > 1)
    DataLabel[i++] = Vel2Name;

  if(MetaData.TopGridRank > 2)
    DataLabel[i++] = Vel3Name;

  if (TestProblemData.MultiSpecies) {
    DataLabel[i++] = ElectronName;
    DataLabel[i++] = HIName;
    DataLabel[i++] = HIIName;
    DataLabel[i++] = HeIName;
    DataLabel[i++] = HeIIName;
    DataLabel[i++] = HeIIIName;
    if (TestProblemData.MultiSpecies > 1) {
      DataLabel[i++] = HMName;
      DataLabel[i++] = H2IName;
      DataLabel[i++] = H2IIName;
    }
    if (TestProblemData.MultiSpecies > 2) {
      DataLabel[i++] = DIName;
      DataLabel[i++] = DIIName;
      DataLabel[i++] = HDIName;
    }
  }
 
  if (TestProblemData.UseMetallicityField) {
    DataLabel[i++] = MetalName;

    if(TestProblemData.MultiMetals){
      DataLabel[i++] = ExtraNames[0];
      DataLabel[i++] = ExtraNames[1];
    }
  }

  if(TestProblemData.GloverChemistryModel){

    int GCM = TestProblemData.GloverChemistryModel;  // purely for convenience

    DataLabel[i++] = HIIName;
    DataLabel[i++] = HIName;
    DataLabel[i++] = H2IName;

    if( (GCM==1) || (GCM==2) || (GCM==3) || (GCM==7) ){
      DataLabel[i++] = DIName;
      DataLabel[i++] = DIIName;
      DataLabel[i++] = HDIName;
      DataLabel[i++] = HeIName;
      DataLabel[i++] = HeIIName;
      DataLabel[i++] = HeIIIName;
    }

    if( (GCM==3) || (GCM==5) || (GCM==7) ){
      DataLabel[i++] = COIName;
    }

    if( (GCM==2) || (GCM==3) || (GCM==7) ){
      DataLabel[i++] = CIName;
      DataLabel[i++] = CIIName;
      DataLabel[i++] = OIName;
      DataLabel[i++] = OIIName;
    }

    if( (GCM==2) || (GCM==3) ){
      DataLabel[i++] = SiIName;
      DataLabel[i++] = SiIIName;
      DataLabel[i++] = SiIIIName;
    }

    if( (GCM==3) || (GCM==7) ){
      DataLabel[i++] = CHIName;
      DataLabel[i++] = CH2IName;
      DataLabel[i++] = CH3IIName;
      DataLabel[i++] = C2IName;
      DataLabel[i++] = HCOIIName;
      DataLabel[i++] = OHIName;
      DataLabel[i++] = H2OIName;
      DataLabel[i++] = O2IName;
    }

  } //   if(TestProblemData.GloverChemistryModel)


  for(j=0; j < i; j++)
    DataUnits[j] = NULL;
 
  /* Write parameters to parameter output file */
 
  if (MyProcessorNumber == ROOT_PROCESSOR) {

    fprintf(Outfptr, "RadiatingShockInnerDensity         = %"FSYM"\n"  , RadiatingShockInnerDensity);
    fprintf(Outfptr, "RadiatingShockOuterDensity         = %"FSYM"\n"  , RadiatingShockOuterDensity);
    fprintf(Outfptr, "RadiatingShockPressure        = %"FSYM"\n"  , RadiatingShockPressure);
    fprintf(Outfptr, "RadiatingShockEnergy          = %"FSYM"\n"  , RadiatingShockEnergy);
    fprintf(Outfptr, "RadiatingShockInnerPressure   = %"FSYM"\n"  ,
	    RadiatingShockInnerPressure);
 
    fprintf(Outfptr, "RadiatingShockUseDensityFluctuations   = %"ISYM"\n", RadiatingShockUseDensityFluctuations);
    fprintf(Outfptr, "RadiatingShockRandomSeed   = %"ISYM"\n", RadiatingShockRandomSeed);
    fprintf(Outfptr, "RadiatingShockDensityFluctuationLevel   = %"FSYM"\n", RadiatingShockDensityFluctuationLevel);
    fprintf(Outfptr, "RadiatingShockInitializeWithKE = %"ISYM"\n", RadiatingShockInitializeWithKE);
    fprintf(Outfptr, "RadiatingShockUseSedovProfile = %"ISYM"\n", RadiatingShockUseSedovProfile);

    fprintf(Outfptr, "RadiatingShockSedovBlastRadius = %"PSYM"\n", RadiatingShockSedovBlastRadius);

    fprintf(Outfptr,  "RadiatingShockKineticEnergyFraction = %"FSYM"\n", RadiatingShockKineticEnergyFraction);

    fprintf(Outfptr, "RadiatingShockCenterPosition = %"PSYM" %"PSYM" %"PSYM"\n",
		  RadiatingShockCenterPosition, RadiatingShockCenterPosition+1,
		  RadiatingShockCenterPosition+2);

    fprintf(Outfptr, "RadiatingShockSpreadOverNumZones  = %"PSYM"\n", RadiatingShockSpreadOverNumZones);

    fprintf(Outfptr, "TestProblemHydrogenFractionByMass = %"FSYM"\n",   TestProblemData.HydrogenFractionByMass);
    fprintf(Outfptr, "TestProblemDeuteriumToHydrogenRatio = %"FSYM"\n", TestProblemData.DeuteriumToHydrogenRatio);

    fprintf(Outfptr, "TestProblemInitialHIFractionInner  = %"FSYM"\n", TestProblemData.HI_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialHIIFractionInner  = %"FSYM"\n", TestProblemData.HII_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialHeIFractionInner  = %"FSYM"\n", TestProblemData.HeI_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialHeIIFractionInner  = %"FSYM"\n", TestProblemData.HeII_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialHeIIIFractionInner  = %"FSYM"\n", TestProblemData.HeIII_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialHMFractionInner  = %"FSYM"\n", TestProblemData.HM_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialH2IFractionInner  = %"FSYM"\n", TestProblemData.H2I_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialH2IIFractionInner  = %"FSYM"\n", TestProblemData.H2II_Fraction_Inner);

    fprintf(Outfptr, "TestProblemInitialDIFractionInner  = %"FSYM"\n", TestProblemData.DI_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialDIIFractionInner  = %"FSYM"\n", TestProblemData.DII_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialHDIFractionInner  = %"FSYM"\n", TestProblemData.HDI_Fraction_Inner);

    fprintf(Outfptr, "TestProblemInitialCOIFractionInner  = %"FSYM"\n", TestProblemData.COI_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialCIFractionInner  = %"FSYM"\n", TestProblemData.CI_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialCIIFractionInner  = %"FSYM"\n", TestProblemData.CII_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialOIFractionInner  = %"FSYM"\n", TestProblemData.OI_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialOIIFractionInner  = %"FSYM"\n", TestProblemData.OII_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialSiIFractionInner  = %"FSYM"\n", TestProblemData.SiI_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialSiIIFractionInner  = %"FSYM"\n", TestProblemData.SiII_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialSiIIIFractionInner  = %"FSYM"\n", TestProblemData.SiIII_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialCHIFractionInner  = %"FSYM"\n", TestProblemData.CHI_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialCH2IFractionInner  = %"FSYM"\n", TestProblemData.CH2I_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialCH3IIFractionInner  = %"FSYM"\n", TestProblemData.CH3II_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialC2IFractionInner  = %"FSYM"\n", TestProblemData.C2I_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialHCOIIFractionInner  = %"FSYM"\n", TestProblemData.HCOII_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialOHIFractionInner  = %"FSYM"\n", TestProblemData.OHI_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialH2OIFractionInner  = %"FSYM"\n", TestProblemData.H2OI_Fraction_Inner);
    fprintf(Outfptr, "TestProblemInitialO2IFractionInner  = %"FSYM"\n", TestProblemData.O2I_Fraction_Inner);

    fprintf(Outfptr, "TestProblemInitialHIFraction  = %"FSYM"\n", TestProblemData.HI_Fraction);
    fprintf(Outfptr, "TestProblemInitialHIIFraction  = %"FSYM"\n", TestProblemData.HII_Fraction);
    fprintf(Outfptr, "TestProblemInitialHeIFraction  = %"FSYM"\n", TestProblemData.HeI_Fraction);
    fprintf(Outfptr, "TestProblemInitialHeIIFraction  = %"FSYM"\n", TestProblemData.HeII_Fraction);
    fprintf(Outfptr, "TestProblemInitialHeIIIFraction  = %"FSYM"\n", TestProblemData.HeIII_Fraction);
    fprintf(Outfptr, "TestProblemInitialHMFraction  = %"FSYM"\n", TestProblemData.HM_Fraction);
    fprintf(Outfptr, "TestProblemInitialH2IFraction  = %"FSYM"\n", TestProblemData.H2I_Fraction);
    fprintf(Outfptr, "TestProblemInitialH2IIFraction  = %"FSYM"\n", TestProblemData.H2II_Fraction);

    fprintf(Outfptr, "TestProblemInitialDIFraction  = %"FSYM"\n", TestProblemData.DI_Fraction);
    fprintf(Outfptr, "TestProblemInitialDIIFraction  = %"FSYM"\n", TestProblemData.DII_Fraction);
    fprintf(Outfptr, "TestProblemInitialHDIFraction  = %"FSYM"\n", TestProblemData.HDI_Fraction);

    fprintf(Outfptr, "TestProblemInitialCOIFraction  = %"FSYM"\n", TestProblemData.COI_Fraction);
    fprintf(Outfptr, "TestProblemInitialCIFraction  = %"FSYM"\n", TestProblemData.CI_Fraction);
    fprintf(Outfptr, "TestProblemInitialCIIFraction  = %"FSYM"\n", TestProblemData.CII_Fraction);
    fprintf(Outfptr, "TestProblemInitialOIFraction  = %"FSYM"\n", TestProblemData.OI_Fraction);
    fprintf(Outfptr, "TestProblemInitialOIIFraction  = %"FSYM"\n", TestProblemData.OII_Fraction);
    fprintf(Outfptr, "TestProblemInitialSiIFraction  = %"FSYM"\n", TestProblemData.SiI_Fraction);
    fprintf(Outfptr, "TestProblemInitialSiIIFraction  = %"FSYM"\n", TestProblemData.SiII_Fraction);
    fprintf(Outfptr, "TestProblemInitialSiIIIFraction  = %"FSYM"\n", TestProblemData.SiIII_Fraction);
    fprintf(Outfptr, "TestProblemInitialCHIFraction  = %"FSYM"\n", TestProblemData.CHI_Fraction);
    fprintf(Outfptr, "TestProblemInitialCH2IFraction  = %"FSYM"\n", TestProblemData.CH2I_Fraction);
    fprintf(Outfptr, "TestProblemInitialCH3IIFraction  = %"FSYM"\n", TestProblemData.CH3II_Fraction);
    fprintf(Outfptr, "TestProblemInitialC2IFraction  = %"FSYM"\n", TestProblemData.C2I_Fraction);
    fprintf(Outfptr, "TestProblemInitialHCOIIFraction  = %"FSYM"\n", TestProblemData.HCOII_Fraction);
    fprintf(Outfptr, "TestProblemInitialOHIFraction  = %"FSYM"\n", TestProblemData.OHI_Fraction);
    fprintf(Outfptr, "TestProblemInitialH2OIFraction  = %"FSYM"\n", TestProblemData.H2OI_Fraction);
    fprintf(Outfptr, "TestProblemInitialO2IFraction  = %"FSYM"\n", TestProblemData.O2I_Fraction);

    fprintf(Outfptr, "TestProblemUseMetallicityField  = %"ISYM"\n", TestProblemData.UseMetallicityField);
    fprintf(Outfptr, "TestProblemInitialMetallicityFraction  = %"FSYM"\n", TestProblemData.MetallicityField_Fraction);

    fprintf(Outfptr, "TestProblemUseMassInjection  = %"ISYM"\n", TestProblemData.UseMassInjection);
    fprintf(Outfptr, "TestProblemInitialHydrogenMass  = %"ESYM"\n", TestProblemData.InitialHydrogenMass);
    fprintf(Outfptr, "TestProblemInitialDeuteriumMass  = %"ESYM"\n", TestProblemData.InitialDeuteriumMass);
    fprintf(Outfptr, "TestProblemInitialHeliumMass  = %"ESYM"\n", TestProblemData.InitialHeliumMass);
    fprintf(Outfptr, "TestProblemInitialMetalMass  = %"ESYM"\n", TestProblemData.InitialMetalMass);

    fprintf(Outfptr, "TestProblemMultiMetals  = %"ISYM"\n", TestProblemData.MultiMetals);
    fprintf(Outfptr, "TestProblemInitialMultiMetalsField1Fraction  = %"FSYM"\n", TestProblemData.MultiMetalsField1_Fraction);
    fprintf(Outfptr, "TestProblemInitialMultiMetalsField2Fraction  = %"FSYM"\n", TestProblemData.MultiMetalsField2_Fraction);

  } //   if (MyProcessorNumber == ROOT_PROCESSOR) 

 
  return SUCCESS;
 
}
