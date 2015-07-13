/***********************************************************************
/
/  ROTATING CYLINDER PROBLEM TYPE
/
/  written by: Matthew Turk, Brian O'Shea
/  date:       July, 2010
/
/  PURPOSE:
/
************************************************************************/

#ifdef NEW_PROBLEM_TYPES
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
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
#include "ProblemType.h"
#include "EventHooks.h"

class ProblemType_CollapsingCoolingCloud;

class CollapsingCoolingCloudGrid : private grid {
    friend class ProblemType_CollapsingCoolingCloud;
};

int FindField(int field, int farray[], int numfields);

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);

/* stuff from BWO for this problem type */
void mt_init(unsigned_int seed);

unsigned_long_int mt_random();

void calculate_radial_profiles(float central_density, float central_temperature, FLOAT outer_radius);
float dTdr(float r, float T);
float Mass_of_r(float r);
float n_of_r(float r);
float dn_dr(float r);
float g_of_r(float r);

#define PC_CGS 3.0857e+18
#define RADIUS_BINS 2048

double mu = 1.22, mp=1.67e-24, kb=1.38e-16, gravconst=6.67e-8;
float n_core, r_core, n0, r0, r_outer, T_center, this_radius;
float numdens_of_r[RADIUS_BINS],radius_bins[RADIUS_BINS],T_of_r[RADIUS_BINS];

/* end of BWO stuff */

class ProblemType_CollapsingCoolingCloud : public EnzoProblemType
{
    private:
        FLOAT CollapsingCoolingCloudSubgridLeft, CollapsingCoolingCloudSubgridRight;
        FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
        FLOAT CollapsingCoolingCloudCenterPosition[MAX_DIMENSION];
        float CollapsingCoolingCloudVelocity[3];   // gas initally at rest
        float CollapsingCoolingCloudBField[3];   // gas initally at rest
        FLOAT CollapsingCoolingCloudRadius;
        float CollapsingCoolingCloudLambda;
        float CollapsingCoolingCloudExternalDensity;
        float CollapsingCoolingCloudCentralDensity;
        float CollapsingCoolingCloudExternalTemperature;
        float CollapsingCoolingCloudCentralTemperature;
        float CollapsingCoolingCloudTotalEnergy;
        int CollapsingCoolingCloudUseDensityFluctuations;
        int CollapsingCoolingCloudRandomSeedInitialize;
        int CollapsingCoolingCloudRandomSeed;
        float CollapsingCoolingCloudFluctuationLevel;

    public:
    ProblemType_CollapsingCoolingCloud() : EnzoProblemType()
    { 
        std::cout << "Creating problem type Collapsing Cloud" << std::endl;
    }

    ~ProblemType_CollapsingCoolingCloud()
    {
    }

    virtual int InitializeFromRestart(
            HierarchyEntry &TopGrid, TopGridData &MetaData)
    {
        return SUCCESS;
    }

    virtual int InitializeSimulation(FILE *fptr, FILE *Outfptr,
            HierarchyEntry &TopGrid, TopGridData &MetaData)
    {
      if(debug){
        printf("Entering CollapsingCoolingCloudInitialize\n");
        fflush(stdout);
      }

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
      char *ExtraNames[2] = {"Z_Field1", "Z_Field2"};

      /* local declarations */

      char line[MAX_LINE_LENGTH];
      int  i, j, dim, ret, NumberOfSubgridZones[MAX_DIMENSION],
           SubgridDims[MAX_DIMENSION];

      /* make sure it is 3D */

      if (MetaData.TopGridRank != 3) {
        printf("Cannot do CollapsingCoolingCloud in %"ISYM" dimension(s)\n", MetaData.TopGridRank);
        ENZO_FAIL("");
      }

      for(i=0; i<MAX_DIMENSION; i++)
        CollapsingCoolingCloudCenterPosition[i] = 0.5;  // right in the middle of the box

      this->CollapsingCoolingCloudVelocity[0] = 
        this->CollapsingCoolingCloudVelocity[1] = 
        this->CollapsingCoolingCloudVelocity[2] = 0.0; // gas initally at rest
      this->CollapsingCoolingCloudBField[0] =
        this->CollapsingCoolingCloudBField[1] =
        this->CollapsingCoolingCloudBField[2] = 0.0; // gas initally at rest
      this->CollapsingCoolingCloudRadius = 0.3;
      this->CollapsingCoolingCloudLambda = 0.05;
      this->CollapsingCoolingCloudCentralDensity = 100.0;
      this->CollapsingCoolingCloudExternalDensity = 1.0;
      this->CollapsingCoolingCloudCentralTemperature = 1000.0;
      this->CollapsingCoolingCloudExternalTemperature = 100.0;
      this->CollapsingCoolingCloudTotalEnergy = 1.0;
      float Pi                      = 3.14159;

      this->CollapsingCoolingCloudUseDensityFluctuations=0;
      this->CollapsingCoolingCloudRandomSeed=123456789;
      this->CollapsingCoolingCloudFluctuationLevel=0.1;

      this->CollapsingCoolingCloudRandomSeedInitialize=0;

      /* set no subgrids by default. */

      this->CollapsingCoolingCloudSubgridLeft         = 0.0;    // start of subgrid(s)
      this->CollapsingCoolingCloudSubgridRight        = 0.0;    // end of subgrid(s)

      TestProblemData.MultiSpecies = MultiSpecies;  // set this from global data (kind of a hack, but necessary)

      /* read input from file */

      while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

        ret = 0;

        /* read parameters specifically for radiating shock problem*/

        ret += sscanf(line, "CollapsingCoolingCloudCentralDensity  = %"FSYM, &CollapsingCoolingCloudCentralDensity);
        ret += sscanf(line, "CollapsingCoolingCloudExternalDensity  = %"FSYM, &CollapsingCoolingCloudExternalDensity);
        ret += sscanf(line, "CollapsingCoolingCloudCentralTemperature  = %"FSYM, &CollapsingCoolingCloudCentralTemperature);
        ret += sscanf(line, "CollapsingCoolingCloudExternalTemperature  = %"FSYM, &CollapsingCoolingCloudExternalTemperature);

        ret += sscanf(line, "CollapsingCoolingCloudSubgridLeft = %"PSYM,
            &CollapsingCoolingCloudSubgridLeft);
        ret += sscanf(line, "CollapsingCoolingCloudSubgridRight = %"PSYM,
            &CollapsingCoolingCloudSubgridRight);
        ret += sscanf(line, "CollapsingCoolingCloudLambda = %"FSYM,
            &CollapsingCoolingCloudLambda);

        ret += sscanf(line, "CollapsingCoolingCloudTotalEnergy = %"FSYM,
            &CollapsingCoolingCloudTotalEnergy);

        ret += sscanf(line, "CollapsingCoolingCloudRadius = %"PSYM,
            &CollapsingCoolingCloudRadius);
        ret += sscanf(line, "CollapsingCoolingCloudCenterPosition = %"PSYM" %"PSYM" %"PSYM,
            CollapsingCoolingCloudCenterPosition, CollapsingCoolingCloudCenterPosition+1,
            CollapsingCoolingCloudCenterPosition+2);

        ret += sscanf(line, "CollapsingCoolingCloudUseDensityFluctuations = %"ISYM,
            &CollapsingCoolingCloudUseDensityFluctuations);
        ret += sscanf(line, "CollapsingCoolingCloudRandomSeed = %"ISYM,
            &CollapsingCoolingCloudRandomSeed);
        ret += sscanf(line, "CollapsingCoolingCloudFluctuationLevel = %"FSYM,
            &CollapsingCoolingCloudFluctuationLevel);


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

        ret += sscanf(line, "TestProblemUseMetallicityField  = %"ISYM, &TestProblemData.UseMetallicityField);
        ret += sscanf(line, "TestProblemInitialMetallicityFraction  = %"FSYM, &TestProblemData.MetallicityField_Fraction);

	ret += sscanf(line, "TestProblemMultiMetals  = %"ISYM, &TestProblemData.MultiMetals);
	ret += sscanf(line, "TestProblemInitialMultiMetalsField1Fraction  = %"FSYM, &TestProblemData.MultiMetalsField1_Fraction);
	ret += sscanf(line, "TestProblemInitialMultiMetalsField2Fraction  = %"FSYM, &TestProblemData.MultiMetalsField2_Fraction);

        /* if the line is suspicious, issue a warning */

        if (ret == 0 && strstr(line, "=") && (strstr(line, "CollapsingCoolingCloud") || strstr(line, "TestProblem")) &&
            line[0] != '#' && MyProcessorNumber == ROOT_PROCESSOR)
          fprintf(stderr,
              "*** warning: the following parameter line was not interpreted:\n%s\n",
              line);

      } // end input from parameter file

      /* Set the units. */
      
      float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1.0,
	TimeUnits = 1.0, VelocityUnits = 1.0;
      double MassUnits=1.0;
      FLOAT Time=0.0;
      if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
		   &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
        ENZO_FAIL("Error in GetUnits.");
      }
      float EnergyUnits;
      float TempToEnergyConversion;
      EnergyUnits = POW(LengthUnits, 2.0) / POW(TimeUnits, 2.0);
      TempToEnergyConversion =  kb/((Gamma - 1.0)*mu*mp); 
      TempToEnergyConversion /= EnergyUnits;  // this times temperature gives you energy units in ENZO UNITS (K -> Enzo)

      // returns three arrays:  n(r), T(r), r, in cm^-3, Kelvin, pc respectively.
      calculate_radial_profiles(CollapsingCoolingCloudCentralDensity,CollapsingCoolingCloudCentralTemperature,CollapsingCoolingCloudRadius);

      // convert from number density, Kelvin, parsecs to Enzo internal units.
      for(int i=0; i<RADIUS_BINS;i++){
	numdens_of_r[i] *= mu*mp/DensityUnits;  // now in Enzo internal density units 
	T_of_r[i] *= TempToEnergyConversion;  // now in internal energy units
	radius_bins[i] /= LengthUnits;
      }

      float ExternalDensity, ExternalEnergy;

      ExternalDensity = numdens_of_r[RADIUS_BINS-1];
      ExternalEnergy = T_of_r[RADIUS_BINS-1];

      printf("internal/external density/energy are:  %e/%e   %e/%e\n",numdens_of_r[0], numdens_of_r[RADIUS_BINS-1],
	     T_of_r[0],T_of_r[RADIUS_BINS-1]);

      CollapsingCoolingCloudRadius *= PC_CGS/LengthUnits;  // now in internal length units

      this->InitializeUniformGrid(TopGrid.GridData,
				  ExternalDensity,
				  ExternalEnergy,
				  ExternalEnergy,
				  CollapsingCoolingCloudVelocity,
				  CollapsingCoolingCloudBField);

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
            nint((CollapsingCoolingCloudSubgridRight - CollapsingCoolingCloudSubgridLeft)/
                ((DomainRightEdge[dim] - DomainLeftEdge[dim] )/
                 float(MetaData.TopGridDims[dim])))
            *int(POW(RefineBy, lev + 1));

        if (debug)
          printf("CollapsingCoolingCloud:: Level[%"ISYM"]: NumberOfSubgridZones[0] = %"ISYM"\n", lev+1,
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
            LeftEdge[dim]    = CollapsingCoolingCloudSubgridLeft;
            RightEdge[dim]   = CollapsingCoolingCloudSubgridRight;
          }

          /* create a new subgrid and initialize it */

          Subgrid[lev]->GridData = this->CreateNewUniformGrid(
                                        TopGrid.GridData,
                                        MetaData.TopGridRank, SubgridDims,
                                        LeftEdge, RightEdge, 0,
                                        ExternalDensity,
					ExternalEnergy,
					ExternalEnergy,
                                        CollapsingCoolingCloudVelocity,
                                        CollapsingCoolingCloudBField);

          /* set up the initial explosion area on the finest resolution subgrid */

          if (lev == MaximumRefinementLevel - 1)
            if (this->InitializeGrid(Subgrid[lev]->GridData, TopGrid, MetaData)
                == FAIL) {
              ENZO_FAIL("Error in CollapsingCoolingCloudInitialize[Sub]Grid.");
            }

        }
        else{
          printf("CollapsingCoolingCloud: single grid start-up.\n");
        }
      }

      this->FinalizeGrids(Subgrid, TopGrid, MetaData);

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

      for(j=0; j < i; j++)
        DataUnits[j] = NULL;

      /* Write parameters to parameter output file */

      if (MyProcessorNumber == ROOT_PROCESSOR) {
        fprintf(Outfptr, "CollapsingCoolingCloudCentralDensity         = %"FSYM"\n"  , CollapsingCoolingCloudCentralDensity);
        fprintf(Outfptr, "CollapsingCoolingCloudExternalDensity         = %"FSYM"\n"  , CollapsingCoolingCloudExternalDensity);
        fprintf(Outfptr, "CollapsingCoolingCloudCentralTemperature         = %"FSYM"\n"  , CollapsingCoolingCloudCentralTemperature);
        fprintf(Outfptr, "CollapsingCoolingCloudExternalTemperature         = %"FSYM"\n"  , CollapsingCoolingCloudExternalTemperature);

        fprintf(Outfptr, "CollapsingCoolingCloudLambda         = %"FSYM"\n"  , CollapsingCoolingCloudLambda);
        fprintf(Outfptr, "CollapsingCoolingCloudTotalEnergy         = %"FSYM"\n"  , CollapsingCoolingCloudTotalEnergy);
        fprintf(Outfptr, "CollapsingCoolingCloudRadius         = %"PSYM"\n"  , CollapsingCoolingCloudRadius);
        fprintf(Outfptr, "CollapsingCoolingCloudCenterPosition = %"PSYM" %"PSYM" %"PSYM"\n",
            CollapsingCoolingCloudCenterPosition, CollapsingCoolingCloudCenterPosition+1,
            CollapsingCoolingCloudCenterPosition+2);

        fprintf(Outfptr, "CollapsingCoolingCloudUseDensityFluctuations    = %"ISYM"\n"  , CollapsingCoolingCloudUseDensityFluctuations);
        fprintf(Outfptr, "CollapsingCoolingCloudRandomSeed             = %"ISYM"\n"  , CollapsingCoolingCloudRandomSeed);
        fprintf(Outfptr, "CollapsingCoolingCloudFluctuationLevel       = %"FSYM"\n"  , CollapsingCoolingCloudFluctuationLevel);

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

        fprintf(Outfptr, "TestProblemUseMetallicityField  = %"ISYM"\n", TestProblemData.UseMetallicityField);
        fprintf(Outfptr, "TestProblemInitialMetallicityFraction  = %"FSYM"\n", TestProblemData.MetallicityField_Fraction);

	fprintf(Outfptr, "TestProblemMultiMetals  = %"ISYM"\n", TestProblemData.MultiMetals);
	fprintf(Outfptr, "TestProblemInitialMultiMetalsField1Fraction  = %"FSYM"\n", TestProblemData.MultiMetalsField1_Fraction);
	fprintf(Outfptr, "TestProblemInitialMultiMetalsField2Fraction  = %"FSYM"\n", TestProblemData.MultiMetalsField2_Fraction);

      } //   if (MyProcessorNumber == ROOT_PROCESSOR) 


      if(debug){
        printf("Exiting CollapsingCoolingCloudInitialize\n");
        fflush(stdout);
      }

      return SUCCESS;

    }

/*

This is the grid-by-grid initializer.

*/
    int InitializeGrid(grid *thisgrid_orig,
            HierarchyEntry &TopGrid, TopGridData &MetaData)
    {

      CollapsingCoolingCloudGrid *thisgrid =
        static_cast<CollapsingCoolingCloudGrid *>(thisgrid_orig);

      if (thisgrid->ProcessorNumber != MyProcessorNumber)
        return SUCCESS;

      if(debug){
        printf("Entering CollapsingCoolingCloudInitializeGrid\n");
        fflush(stdout);
      }

      printf("CollapsingCoolingCloudRadius = %e\n", this->CollapsingCoolingCloudRadius);
      printf("CollapsingCoolingCloudCenterPosition = %e %e %e\n", 
          this->CollapsingCoolingCloudCenterPosition[0],
          this->CollapsingCoolingCloudCenterPosition[1],
          this->CollapsingCoolingCloudCenterPosition[2]);
      printf("CollapsingCoolingCloudLambda = %e\n",this->CollapsingCoolingCloudLambda);
      printf("CollapsingCoolingCloudCentralDensity = %e\n",this->CollapsingCoolingCloudCentralDensity);


      /* declarations */

      int size = 1, dim, cellindex;

      for (dim = 0; dim < thisgrid->GridRank; dim++)
        size *= thisgrid->GridDimension[dim];

      FLOAT r,x,y,z, radius, xyradius, zdist;

      float sintheta, costheta, omega;

      int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, MetalNum;

      int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
	DINum, DIINum, HDINum;

      if (thisgrid->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
            Vel3Num, TENum) == FAIL) {
        fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
        ENZO_FAIL("");
      }

      if (MultiSpecies)
	if (thisgrid->IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
				  HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL) {
	  ENZO_FAIL("Error in grid->IdentifySpeciesFields.");
	}

      int MetallicityField = FALSE;
      if ((MetalNum = FindField(Metallicity, thisgrid->FieldType, thisgrid->NumberOfBaryonFields))
          != -1)
        MetallicityField = TRUE;
      else
        MetalNum = 0;

      float therandomfraction;

      if( CollapsingCoolingCloudRandomSeedInitialize == 0)
	mt_init(((unsigned_int) CollapsingCoolingCloudRandomSeed));
      
      unsigned_long_int therandominteger;

      /* set fields in the cylinder region */

      int index, jndex, i, j, k;
      float outside_rho, outside_TE, outside_GE;

      outside_rho =  thisgrid->BaryonField[DensNum][0];

      omega = CollapsingCoolingCloudLambda * sqrt(GravitationalConstant * numdens_of_r[0]/numdens_of_r[RADIUS_BINS-1]) / 0.117;

      if(HydroMethod==2){  // ZEUS

        outside_TE = thisgrid->BaryonField[TENum][0];

      } else { // PPM

        outside_TE = thisgrid->BaryonField[TENum][0];

        if(DualEnergyFormalism){
          outside_GE = thisgrid->BaryonField[GENum][0];
        }

      }  // if(HydroMethod==2)

      int cell_radial_index;


      for (k = 0; k < thisgrid->GridDimension[2]; k++)
        for (j = 0; j < thisgrid->GridDimension[1]; j++)
          for (i = 0; i < thisgrid->GridDimension[0]; i++){

            /* Compute position */
            x=y=z=0.0;

            cellindex = i + j*thisgrid->GridDimension[0]
                          + k*thisgrid->GridDimension[0]
                             *thisgrid->GridDimension[1];

            x = thisgrid->CellLeftEdge[0][i] + 0.5*thisgrid->CellWidth[0][i];
            y = thisgrid->CellLeftEdge[1][j] + 0.5*thisgrid->CellWidth[1][j];
            z = thisgrid->CellLeftEdge[2][k] + 0.5*thisgrid->CellWidth[2][k];

            /* Find distance from center. */

            // it's REALLY r^2 right now
            xyradius = POW(x-CollapsingCoolingCloudCenterPosition[0], 2.0) +
              POW(y-CollapsingCoolingCloudCenterPosition[1], 2.0);

	    radius = xyradius + POW(z-CollapsingCoolingCloudCenterPosition[2], 2.0);

            xyradius = sqrt(xyradius);  // ok, now it's just xy-radius
            radius = sqrt(radius);  // ok, now it's just radius

            if ( radius <= CollapsingCoolingCloudRadius ){

	      for(int thisindex=1;thisindex<RADIUS_BINS;thisindex++)
		if(radius_bins[thisindex-1] <= radius &&
		   radius < radius_bins[thisindex])
		  cell_radial_index=thisindex;
	      
	      therandomfraction=1.0;

	      if(CollapsingCoolingCloudUseDensityFluctuations){
		therandominteger = mt_random();

		therandomfraction =  float((therandominteger%32768)) /  32768.0;

		if(therandomfraction < 0.0 || therandomfraction > 1.0){
		  fprintf(stderr,"yarr!  random number generator went bad!  %e\n",therandomfraction);
		  return FAIL;
		}

		therandomfraction =  1.0 +  CollapsingCoolingCloudFluctuationLevel*(therandomfraction - 0.5);

		//printf(stderr,"therandomfraction = %e\n"FSYM,therandomfraction);

	      }

              thisgrid->BaryonField[DensNum][cellindex] = therandomfraction * numdens_of_r[cell_radial_index];

              sintheta = (y-CollapsingCoolingCloudCenterPosition[1])/xyradius;
              costheta = (x-CollapsingCoolingCloudCenterPosition[0])/xyradius;


              // x,y, and maybe z velocity.  
              thisgrid->BaryonField[Vel1Num][cellindex] = -1.0*sintheta*omega*xyradius;

              thisgrid->BaryonField[Vel2Num][cellindex] = costheta*omega*xyradius;

              thisgrid->BaryonField[Vel3Num][cellindex] = 0.0;

              if(HydroMethod == 2){

                // ZEUS
                thisgrid->BaryonField[TENum][cellindex] = T_of_r[cell_radial_index] ;

              } else {

                // PPM
                thisgrid->BaryonField[TENum][cellindex] =  T_of_r[cell_radial_index] 
                  + 0.5 * thisgrid->BaryonField[Vel1Num][cellindex] *
                  thisgrid->BaryonField[Vel1Num][cellindex]
                  + 0.5 * thisgrid->BaryonField[Vel2Num][cellindex] *
                  thisgrid->BaryonField[Vel2Num][cellindex]
                  + 0.5 * thisgrid->BaryonField[Vel3Num][cellindex] *
                  thisgrid->BaryonField[Vel3Num][cellindex];

                // gas energy (PPM dual energy formalims)
                if(DualEnergyFormalism)
                  thisgrid->BaryonField[GENum][cellindex] =  T_of_r[cell_radial_index]; 

              } // if(HydroMethod == 2)

	  // Set multispecies fields!
	  // this attempts to set them such that species conservation is maintained,
	  // using the method in CosmologySimulationInitializeGrid.C
	  if(TestProblemData.MultiSpecies) {

	    thisgrid->BaryonField[HIINum][cellindex] = TestProblemData.HII_Fraction_Inner * thisgrid->BaryonField[DensNum][cellindex];
	      
	    thisgrid->BaryonField[HeIINum][cellindex] = TestProblemData.HeII_Fraction_Inner * thisgrid->BaryonField[DensNum][cellindex];
	      
	    thisgrid->BaryonField[HeIIINum][cellindex] = TestProblemData.HeIII_Fraction_Inner * thisgrid->BaryonField[DensNum][cellindex];

	    thisgrid->BaryonField[HeINum][cellindex] = (1.0-TestProblemData.HydrogenFractionByMass)*thisgrid->BaryonField[DensNum][cellindex] -
	      thisgrid->BaryonField[HeIINum][cellindex] - thisgrid->BaryonField[HeIIINum][cellindex];
	      
	    if(TestProblemData.MultiSpecies > 1){
	      thisgrid->BaryonField[HMNum][cellindex] = TestProblemData.HM_Fraction_Inner * thisgrid->BaryonField[HIINum][cellindex];
		
	      thisgrid->BaryonField[H2INum][cellindex] = TestProblemData.H2I_Fraction_Inner * thisgrid->BaryonField[0][cellindex];
		
	      thisgrid->BaryonField[H2IINum][cellindex] = TestProblemData.H2II_Fraction_Inner * 2.0 * thisgrid->BaryonField[HIINum][cellindex];
	    }

	    // HI density is calculated by subtracting off the various ionized fractions
	    // from the total
	    thisgrid->BaryonField[HINum][cellindex] = TestProblemData.HydrogenFractionByMass*thisgrid->BaryonField[0][cellindex]
	      - thisgrid->BaryonField[HIINum][cellindex];

	    if (MultiSpecies > 1)
	      thisgrid->BaryonField[HINum][cellindex] -= (thisgrid->BaryonField[HMNum][cellindex] + thisgrid->BaryonField[H2IINum][cellindex]
						+ thisgrid->BaryonField[H2INum][cellindex]);

	    // Electron "density" (remember, this is a factor of m_p/m_e scaled from the 'normal'
	    // density for convenience) is calculated by summing up all of the ionized species.
	    // The factors of 0.25 and 0.5 in front of HeII and HeIII are to fix the fact that we're
	    // calculating mass density, not number density (because the BaryonField values are 4x as
	    // heavy for helium for a single electron)
	    thisgrid->BaryonField[DeNum][cellindex] = thisgrid->BaryonField[HIINum][cellindex] +
	      0.25*thisgrid->BaryonField[HeIINum][cellindex] + 0.5*thisgrid->BaryonField[HeIIINum][cellindex];

	    if (MultiSpecies > 1)
	      thisgrid->BaryonField[DeNum][cellindex] += 0.5*thisgrid->BaryonField[H2IINum][cellindex] -
		thisgrid->BaryonField[HMNum][cellindex];
	      
	    // Set deuterium species (assumed to be a negligible fraction of the total, so not
	    // counted in the conservation)
	    if(TestProblemData.MultiSpecies > 2){
	      thisgrid->BaryonField[DINum ][cellindex] = TestProblemData.DeuteriumToHydrogenRatio * thisgrid->BaryonField[HINum][cellindex];
	      thisgrid->BaryonField[DIINum][cellindex] = TestProblemData.DeuteriumToHydrogenRatio * thisgrid->BaryonField[HIINum][cellindex];
	      thisgrid->BaryonField[HDINum][cellindex] = 0.75 * TestProblemData.DeuteriumToHydrogenRatio * thisgrid->BaryonField[H2INum][cellindex];
	    }

	  } // if(TestProblemData.MultiSpecies)

              if(TestProblemData.UseMetallicityField > 0 && MetalNum != FALSE)
                thisgrid->BaryonField[MetalNum][cellindex] = thisgrid->BaryonField[DensNum][cellindex]*TestProblemData.MetallicityField_Fraction;

	      if(TestProblemData.UseMetallicityField){
		thisgrid->BaryonField[MetalNum][cellindex] = TestProblemData.MetallicityField_Fraction * thisgrid->BaryonField[DensNum][cellindex];

		if(TestProblemData.MultiMetals){
		  thisgrid->BaryonField[MetalNum+1][cellindex] = TestProblemData.MultiMetalsField1_Fraction * thisgrid->BaryonField[DensNum][cellindex];
		  thisgrid->BaryonField[MetalNum+2][cellindex] = TestProblemData.MultiMetalsField2_Fraction * thisgrid->BaryonField[DensNum][cellindex];
		}
	      } // if(TestProblemData.UseMetallicityField)
	      

            } // if (r <= CollapsingCoolingCloudRadius)
	    else {  // set external density, temperature so it all doesn't fall apart.
              thisgrid->BaryonField[DensNum][cellindex] = numdens_of_r[RADIUS_BINS-1];
	      thisgrid->BaryonField[TENum][cellindex] = T_of_r[RADIUS_BINS-1] ;
	      if(DualEnergyFormalism)
		thisgrid->BaryonField[GENum][cellindex] = T_of_r[RADIUS_BINS-1] ;

	    } // else if  (r <= CollapsingCoolingCloudRadius)


          } // for (i = 0; i < GridDimension[0]; i++)

      if(debug){
        printf("Exiting CollapsingCoolingCloudInitialize\n");
        fflush(stdout);
      }

      return SUCCESS;

    } 

};

//.. register:
namespace{
    EnzoProblemType_creator_concrete<ProblemType_CollapsingCoolingCloud>
        collapsing_cooling_cloud("CollapsingCoolingCloud");
}



void calculate_radial_profiles(float central_density, float central_temperature, FLOAT outer_radius){

  // user sets these
  n_core = central_density;  // n_H: particles/cc

  r_outer = outer_radius; // pc
  T_center = central_temperature;  // Kelvin

  double dr, k1, k2, k3, k4, this_temperature;

  n0 = 1.0e+3; // particles/cc
  r0 = 1.0 * PC_CGS;  // in parsecs

  r_core = r0 * pow( n_core/n0, -1.0/2.2);

  printf("n0,n_c, r0, r_c = %e %e    %e %e\n",n0,n_core,r0,r_core);

  dr = (r_outer - 0.0) * PC_CGS / RADIUS_BINS;  // dr in cm 

  this_temperature = T_center;

  this_radius = dr/10.0;  // small starting radius (in cm)

  for(int i=0; i<RADIUS_BINS; i++){
    numdens_of_r[i]=radius_bins[i]=T_of_r[i]= -1.0;
  }

  printf("%e    %e    %e\n", this_radius/PC_CGS, this_temperature, n_of_r(this_radius) );

  int counter=0;

  T_of_r[counter]=this_temperature;  // Kelvin
  numdens_of_r[counter]=n_of_r(this_radius);  // particles/CC
  radius_bins[counter]=this_radius;  // in CGS
  counter++;

  while(this_radius <= r_outer*PC_CGS && counter < RADIUS_BINS){

    k1 = dTdr(this_radius,          this_temperature);
    k2 = dTdr(this_radius + 0.5*dr, this_temperature + 0.5*dr*k1);
    k3 = dTdr(this_radius + 0.5*dr, this_temperature + 0.5*dr*k2);
    k4 = dTdr(this_radius + dr,     this_temperature + dr*k3);

    this_temperature += (1.0/6.0) * dr * (k1 + 2.0*k2 + 2.0*k3 + k4);
    this_radius += dr;  // new radius

    T_of_r[counter]=this_temperature;  // Kelvin
    numdens_of_r[counter]=n_of_r(this_radius);  // particles/CC
    radius_bins[counter]=this_radius;  // in CGS

    printf("%e    %e    %e   %e  %e  %d\n", this_radius/PC_CGS, this_temperature, n_of_r(this_radius), dr, r_outer*PC_CGS, counter );

    counter++;

  } // while()

}

float dTdr(float r, float T){
  return (-1.0*dn_dr(r) * (T / n_of_r(r)) - g_of_r(r) * mu*mp/kb);
}

float n_of_r(float r){
  if(r <= r_core)
    return n_core;
  else
    return n0*pow(r/r0, -2.2);
}

float dn_dr(float r){
  if(r <= r_core)
    return 0.0;
  else
    return -2.2*(n0/r0)*pow(r/r0, -3.2);
}

float g_of_r(float r){
  return -gravconst*Mass_of_r(r)/(r*r);
}

float Mass_of_r(float r){
  if(r <= r_core)
    return ( (4.0/3.0) * 3.14159 * pow(r, 3.0) * mu * mp * n_core);
   else 
    return ( (4.0/3.0) * 3.14159 * pow(r_core, 3.0) * mu * mp * n_core
	     + 5.0 * 3.14159 * mu * mp * n0 * pow(r0, 2.2) * (pow(r,0.8) - pow(r_core, 0.8) ) );
}


#endif
