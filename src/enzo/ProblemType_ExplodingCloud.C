/***********************************************************************
/
/  EXPLODING CLOUD PROBLEM TYPE
/
/  written by: Elizabeth Tasker
/  date:       September, 2011
/
/  PURPOSE: Sedov-like explosion inside a dense cloud of gas
/
************************************************************************/

#ifdef NEW_PROBLEM_TYPES
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <iostream>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "phys_constants.h"
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

class ProblemType_ExplodingCloud;

class ExplodingCloudGrid : private grid 
{
    friend class ProblemType_ExplodingCloud;
};

int FindField(int field, int farray[], int numfields);

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);

float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1.0, TimeUnits = 1.0, VelocityUnits = 1.0, EnergyUnits = 1.0, TempToEnergyConversion = 1.0;
double MassUnits=1.0;

class ProblemType_ExplodingCloud : public EnzoProblemType
{
    private:
        FLOAT ExplodingCloudSubgridLeft, ExplodingCloudSubgridRight;
        FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];

        FLOAT ExplodingCloudEnergyDumpPosition[MAX_DIMENSION];
        FLOAT ExplodingCloudCenterPosition[MAX_DIMENSION];

        FLOAT ExplodingCloudRadius;
   
        float ExplodingCloudExternalDensity;
        float ExplodingCloudInnerDensity;
        float ExplodingCloudExternalTemperature;
        float ExplodingCloudInnerTemperature;
        float ExplodingCloudTotalInputEnergy;

  public:
    ProblemType_ExplodingCloud() : EnzoProblemType()
    { 
        std::cout << "Creating problem type Exploding Cloud" << std::endl;
    }

    ~ProblemType_ExplodingCloud()
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
        printf("Entering ExplodingCloudInitialize\n");
        fflush(stdout);
      }

      char *DensName = "Density";
      char *TEName   = "TotalEnergy";
      char *GEName   = "GasEnergy";
      char *Vel1Name = "x-velocity";
      char *Vel2Name = "y-velocity";
      char *Vel3Name = "z-velocity";
      
      /* local declarations */

      char line[MAX_LINE_LENGTH];
      int  i, j, dim, ret, NumberOfSubgridZones[MAX_DIMENSION],
	SubgridDims[MAX_DIMENSION];

      for (dim = 0; dim < MAX_DIMENSION; dim++) {
        ExplodingCloudCenterPosition[dim] = 0.5*(DomainRightEdge[dim]-DomainLeftEdge[dim]);  // right in the middle of the box
	ExplodingCloudEnergyDumpPosition[dim] = 0.5*(DomainRightEdge[dim]-DomainLeftEdge[dim]);  
      }
      this->ExplodingCloudRadius = 0.045;                                             // [kpc]
      this->ExplodingCloudInnerDensity = 2.47e9;                                  // [Ms / kpc3] ~ 100 cm-3
      this->ExplodingCloudExternalDensity = 8.83e6;                             // [Ms / kpc3] ~ 0.3 cm-3
      this->ExplodingCloudInnerTemperature = 18.0;                              // [K]
      this->ExplodingCloudExternalTemperature = 6000.0;                     // [K]
      this->ExplodingCloudTotalInputEnergy = 5e52;                              // [ergs]
  
      /* set no subgrids by default. */

      this->ExplodingCloudSubgridLeft         = 0.0;    
      this->ExplodingCloudSubgridRight        = 0.0;    

      /* read input from file */

      while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

	ret = 0;

	/* read parameters */
	
	ret += sscanf(line, "ExplodingCloudExternalDensity  = %"FSYM, &ExplodingCloudExternalDensity);
	ret += sscanf(line, "ExplodingCloudInnerDensity = %"FSYM, &ExplodingCloudInnerDensity);
	ret += sscanf(line, "ExplodingCloudExternalTemperature  = %"FSYM, &ExplodingCloudExternalTemperature);
	ret += sscanf(line, "ExplodingCloudInnerTemperature = %"FSYM, &ExplodingCloudInnerTemperature);
	ret += sscanf(line, "ExplodingCloudTotalInputEnergy   = %"FSYM, &ExplodingCloudTotalInputEnergy);
	ret += sscanf(line, "ExplodingCloudRadius = %"FSYM, &ExplodingCloudRadius);
	ret += sscanf(line, "ExplodingCloudSubgridLeft = %"FSYM, &ExplodingCloudSubgridLeft);
	ret += sscanf(line, "ExplodingCloudSubgridRight = %"FSYM, &ExplodingCloudSubgridRight);
	ret += sscanf(line, "ExplodingCloudEnergyDumpPosition = %"PSYM" %"PSYM" %"PSYM, 
		      ExplodingCloudEnergyDumpPosition, ExplodingCloudEnergyDumpPosition+1, ExplodingCloudEnergyDumpPosition+2);
	ret += sscanf(line, "ExplodingCloudCenterPosition = %"PSYM" %"PSYM" %"PSYM, 
		      ExplodingCloudCenterPosition, ExplodingCloudCenterPosition+1, ExplodingCloudCenterPosition+2);

	/* if the line is suspicious, issue a warning */

	if (ret == 0 && strstr(line, "=") && strstr(line, "ExplodingCloud") && 
	    line[0] != '#' && MyProcessorNumber == ROOT_PROCESSOR)
	  fprintf(stderr, 
		  "warning: the following parameter line was not interpreted:\n%s\n",
		  line);

      } // end input from parameter file
      
      /* Set the units. */
      
      FLOAT Time=0.0;
      if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
        ENZO_FAIL("Error in GetUnits.");
      }
      EnergyUnits = POW(LengthUnits, 2.0) / POW(TimeUnits, 2.0);
      TempToEnergyConversion =  kboltz/((Gamma - 1.0)*Mu*mh); 
      TempToEnergyConversion /= EnergyUnits;  // this times temperature gives you energy units in ENZO UNITS (K -> Enzo)

      float ExternalEnergy = ExplodingCloudExternalTemperature*TempToEnergyConversion;
      float ExternalVelocity[] = {0.0, 0.0, 0.0};
      float ExternalBField[] = {0.0, 0.0, 0.0};

      this->InitializeUniformGrid(TopGrid.GridData,
				  ExplodingCloudExternalDensity,
				  ExternalEnergy,
				  ExternalEnergy,
				  ExternalVelocity,
				  ExternalBField);
      

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
            nint((ExplodingCloudSubgridRight - ExplodingCloudSubgridLeft)/
		 ((DomainRightEdge[dim] - DomainLeftEdge[dim] )/
                 float(MetaData.TopGridDims[dim])))
            *int(POW(RefineBy, lev + 1));

        if (debug)
          printf("ExplodingCloud:: Level[%"ISYM"]: NumberOfSubgridZones[0] = %"ISYM"\n", lev+1,
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
            LeftEdge[dim]    = ExplodingCloudSubgridLeft;
            RightEdge[dim]   = ExplodingCloudSubgridRight;
          }

          /* create a new subgrid and initialize it */

          Subgrid[lev]->GridData = this->CreateNewUniformGrid(
							      TopGrid.GridData,
							      MetaData.TopGridRank, SubgridDims,
							      LeftEdge, RightEdge, 0,
							      ExplodingCloudExternalDensity,
							      ExternalEnergy,
							      ExternalEnergy,
							      ExternalVelocity,
							      ExternalBField);

	  /* set up the initial explosion area on the finest resolution subgrid */

          if (lev == MaximumRefinementLevel - 1)
            if (this->InitializeGrid(Subgrid[lev]->GridData, TopGrid, MetaData)
                == FAIL) {
              ENZO_FAIL("Error in ExplodingCloudInitialize[Sub]Grid.");
            }

        }
        else {
          printf("ExplodingCloud: single grid start-up.\n");
        }
      }

      this->FinalizeGrids(Subgrid, TopGrid, MetaData);

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

      for(j=0; j < i; j++)
        DataUnits[j] = NULL;

      /* Write parameters to parameter output file */

      if (MyProcessorNumber == ROOT_PROCESSOR) {
	
        fprintf(Outfptr, "ExplodingCloudInnerDensity         = %"FSYM"\n"  , ExplodingCloudInnerDensity);
        fprintf(Outfptr, "ExplodingCloudExternalDensity         = %"FSYM"\n"  , ExplodingCloudExternalDensity);
        fprintf(Outfptr, "ExplodingCloudInnerTemperature         = %"FSYM"\n"  , ExplodingCloudInnerTemperature);
        fprintf(Outfptr, "ExplodingCloudExternalTemperature         = %"FSYM"\n"  , ExplodingCloudExternalTemperature);
	fprintf(Outfptr, "ExplodingCloudRadius         = %"FSYM"\n"  , ExplodingCloudRadius);
	fprintf(Outfptr, "ExplodingCloudTotalInputEnergy         = %"FSYM"\n"  , ExplodingCloudTotalInputEnergy);
	fprintf(Outfptr, "ExplodingCloudCenterPosition = %"PSYM" %"PSYM" %"PSYM"\n",
		ExplodingCloudCenterPosition, ExplodingCloudCenterPosition+1, ExplodingCloudCenterPosition+2);
	fprintf(Outfptr, "ExplodingCloudEnergyDumpPosition = %"PSYM" %"PSYM" %"PSYM"\n",
		ExplodingCloudEnergyDumpPosition, ExplodingCloudEnergyDumpPosition+1, ExplodingCloudEnergyDumpPosition+2);
	
      } //   if (MyProcessorNumber == ROOT_PROCESSOR) 
      
      if ( debug ) 
	{	
	  printf("Exiting ExplodingCloudInitialize. \n");
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

    ExplodingCloudGrid *thisgrid =
      static_cast<ExplodingCloudGrid *>(thisgrid_orig);

    if (thisgrid->ProcessorNumber != MyProcessorNumber)
      return SUCCESS;

    if (debug)
      {
	printf("Entering ExplodingCloudInitializeGrid\n");
	fflush(stdout);
      }

    int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
    if (thisgrid->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
					     Vel3Num, TENum) == FAIL) 
      ENZO_FAIL("Error in IdentifyPhysicalQuantities.");

     /* declarations */

    int i, j, k, n = 0, size = 1, dim,  cellindex;
    float x, y, z, r, cellvolume = 1.0;

    // Calculate size of baryon field
    for (dim = 0; dim < thisgrid->GridRank; dim++) {
      size *= thisgrid->GridDimension[dim];
      cellvolume*= thisgrid->CellWidth[dim][0];
    }


    // Iterative over grid cells putting in cloud
    for (k = 0; k < thisgrid->GridDimension[2]; k++)
      for (j = 0; j < thisgrid->GridDimension[1]; j++)
	for (i = 0; i < thisgrid->GridDimension[0]; i++, n++ ) {
	    
	    // Physical position from cloud center
	    x = thisgrid->CellLeftEdge[0][i] + 0.5*thisgrid->CellWidth[0][i] 
	      - ExplodingCloudCenterPosition[0];
            y = thisgrid->CellLeftEdge[1][j] + 0.5*thisgrid->CellWidth[1][j] 
	      - ExplodingCloudCenterPosition[1];
            z = thisgrid->CellLeftEdge[2][k] + 0.5*thisgrid->CellWidth[2][k]
	      - ExplodingCloudCenterPosition[2];

	    r = sqrt(x*x + y*y +z*z);

	    if (r <= ExplodingCloudRadius) {
	      // Inside cloud
	      
	      thisgrid->BaryonField[DensNum][n] = ExplodingCloudInnerDensity;
	      thisgrid->BaryonField[TENum][n] = ExplodingCloudInnerTemperature*TempToEnergyConversion;
	      if (DualEnergyFormalism)
		thisgrid->BaryonField[GENum][n] = ExplodingCloudInnerTemperature*TempToEnergyConversion;
	      
	    } // end if inside cloud

	    for (dim = 0; dim < thisgrid->GridRank; dim++)
	      thisgrid->BaryonField[Vel1Num+dim][n] = 0.0;
	    
	} // end loop over grid

    // Add in the energy dump

    i = int(thisgrid->GridDimension[0]*(ExplodingCloudEnergyDumpPosition[0]-thisgrid->GridLeftEdge[0])/(thisgrid->GridRightEdge[0]-thisgrid->GridLeftEdge[0]));
    j = int(thisgrid->GridDimension[1]*(ExplodingCloudEnergyDumpPosition[1]-thisgrid->GridLeftEdge[1])/(thisgrid->GridRightEdge[1]-thisgrid->GridLeftEdge[1]));
    k = int(thisgrid->GridDimension[2]*(ExplodingCloudEnergyDumpPosition[2]-thisgrid->GridLeftEdge[2])/(thisgrid->GridRightEdge[2]-thisgrid->GridLeftEdge[2]));
    n =  i + j*thisgrid->GridDimension[0] + k*thisgrid->GridDimension[0]*thisgrid->GridDimension[1];

    thisgrid->BaryonField[TENum][n] = ExplodingCloudTotalInputEnergy/EnergyUnits/MassUnits/(thisgrid->BaryonField[DensNum][n]*cellvolume);
    if (DualEnergyFormalism)
          thisgrid->BaryonField[GENum][n] = ExplodingCloudTotalInputEnergy/EnergyUnits/MassUnits/(thisgrid->BaryonField[DensNum][n]*cellvolume);


    if ( debug ) {
      printf("Exiting CollapsingCoolingCloudInitialize\n");
      fflush(stdout);
    }
    
    return SUCCESS;
      
  }

};

//.. register:                                                                                                 
namespace{                                                                                                     
  EnzoProblemType_creator_concrete<ProblemType_ExplodingCloud>                                       
    exploding_cloud("ExplodingCloud");                                                    
}                                                                                                              
                                                                                                               

#endif


    
      
	  
	      
