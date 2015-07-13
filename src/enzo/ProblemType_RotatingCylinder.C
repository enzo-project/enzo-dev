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

class ProblemType_RotatingCylinder;

class RotatingCylinderGrid : private grid {
    friend class ProblemType_RotatingCylinder;
};

int FindField(int field, int farray[], int numfields);

void JustPrintSomething(HierarchyEntry *Grids[], TopGridData &MetaData)
{
    std::cout << "I am being called!" << std::endl;
}

class ProblemType_RotatingCylinder : public EnzoProblemType
{
    private:
        FLOAT RotatingCylinderSubgridLeft, RotatingCylinderSubgridRight;
        FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
        FLOAT RotatingCylinderCenterPosition[MAX_DIMENSION];
        float RotatingCylinderVelocity[3];   // gas initally at rest
        float RotatingCylinderBField[3];   // gas initally at rest
        FLOAT RotatingCylinderRadius;
        float RotatingCylinderLambda;
        float RotatingCylinderOverdensity;
        float RotatingCylinderDensity;
        float RotatingCylinderTotalEnergy;

    public:
    ProblemType_RotatingCylinder() : EnzoProblemType()
    { 
        std::cout << "Creating problem type Rotating Cylinder" << std::endl;
        RegisterEventPlugin("Printing", &JustPrintSomething);
        RegisterEventHook("EvolveLevelTop", "Printing");
    }

    ~ProblemType_RotatingCylinder()
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
        printf("Entering RotatingCylinderInitialize\n");
        fflush(stdout);
      }

      char *DensName = "Density";
      char *TEName   = "TotalEnergy";
      char *GEName   = "GasEnergy";
      char *Vel1Name = "x-velocity";
      char *Vel2Name = "y-velocity";
      char *Vel3Name = "z-velocity";
      char *MetalName = "Metal_Density";

      /* local declarations */

      char line[MAX_LINE_LENGTH];
      int  i, j, dim, ret, NumberOfSubgridZones[MAX_DIMENSION],
           SubgridDims[MAX_DIMENSION];

      /* make sure it is 3D */

      if (MetaData.TopGridRank != 3) {
        printf("Cannot do RotatingCylinder in %"ISYM" dimension(s)\n", MetaData.TopGridRank);
        ENZO_FAIL("");
      }

      for(i=0; i<MAX_DIMENSION; i++)
        RotatingCylinderCenterPosition[i] = 0.5;  // right in the middle of the box

      this->RotatingCylinderVelocity[0] = 
        this->RotatingCylinderVelocity[1] = 
        this->RotatingCylinderVelocity[2] = 0.0; // gas initally at rest
      this->RotatingCylinderBField[0] =
        this->RotatingCylinderBField[1] =
        this->RotatingCylinderBField[2] = 0.0; // gas initally at rest
      this->RotatingCylinderRadius = 0.3;
      this->RotatingCylinderLambda = 0.05;
      this->RotatingCylinderOverdensity = 20.0;
      this->RotatingCylinderDensity = 1.0;
      this->RotatingCylinderTotalEnergy = 1.0;
      float Pi                      = 3.14159;

      /* set no subgrids by default. */

      RotatingCylinderSubgridLeft         = 0.0;    // start of subgrid(s)
      RotatingCylinderSubgridRight        = 0.0;    // end of subgrid(s)

      /* read input from file */

      while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

        ret = 0;

        /* read parameters specifically for radiating shock problem*/

        ret += sscanf(line, "RotatingCylinderOverdensity  = %"FSYM, &RotatingCylinderOverdensity);
        ret += sscanf(line, "RotatingCylinderSubgridLeft = %"PSYM,
            &RotatingCylinderSubgridLeft);
        ret += sscanf(line, "RotatingCylinderSubgridRight = %"PSYM,
            &RotatingCylinderSubgridRight);
        ret += sscanf(line, "RotatingCylinderLambda = %"FSYM,
            &RotatingCylinderLambda);

        ret += sscanf(line, "RotatingCylinderTotalEnergy = %"FSYM,
            &RotatingCylinderTotalEnergy);

        ret += sscanf(line, "RotatingCylinderRadius = %"PSYM,
            &RotatingCylinderRadius);
        ret += sscanf(line, "RotatingCylinderCenterPosition = %"PSYM" %"PSYM" %"PSYM,
            RotatingCylinderCenterPosition, RotatingCylinderCenterPosition+1,
            RotatingCylinderCenterPosition+2);

        ret += sscanf(line, "TestProblemUseMetallicityField  = %"ISYM, &TestProblemData.UseMetallicityField);
        ret += sscanf(line, "TestProblemInitialMetallicityFraction  = %"FSYM, &TestProblemData.MetallicityField_Fraction);

        /* if the line is suspicious, issue a warning */

        if (ret == 0 && strstr(line, "=") && (strstr(line, "RotatingCylinder") || strstr(line, "TestProblem")) &&
            line[0] != '#' && MyProcessorNumber == ROOT_PROCESSOR)
          fprintf(stderr,
              "*** warning: the following parameter line was not interpreted:\n%s\n",
              line);

      } // end input from parameter file


      this->InitializeUniformGrid(TopGrid.GridData,
            RotatingCylinderDensity,
            RotatingCylinderTotalEnergy,
            RotatingCylinderTotalEnergy,
            RotatingCylinderVelocity,
            RotatingCylinderBField);

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
            nint((RotatingCylinderSubgridRight - RotatingCylinderSubgridLeft)/
                ((DomainRightEdge[dim] - DomainLeftEdge[dim] )/
                 float(MetaData.TopGridDims[dim])))
            *int(POW(RefineBy, lev + 1));

        if (debug)
          printf("RotatingCylinder:: Level[%"ISYM"]: NumberOfSubgridZones[0] = %"ISYM"\n", lev+1,
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
            LeftEdge[dim]    = RotatingCylinderSubgridLeft;
            RightEdge[dim]   = RotatingCylinderSubgridRight;
          }

          /* create a new subgrid and initialize it */

          Subgrid[lev]->GridData = this->CreateNewUniformGrid(
                                        TopGrid.GridData,
                                        MetaData.TopGridRank, SubgridDims,
                                        LeftEdge, RightEdge, 0,
                                        RotatingCylinderDensity,
                                        RotatingCylinderTotalEnergy,
                                        RotatingCylinderTotalEnergy,
                                        RotatingCylinderVelocity,
                                        RotatingCylinderBField);

          /* set up the initial explosion area on the finest resolution subgrid */

          if (lev == MaximumRefinementLevel - 1)
            if (this->InitializeGrid(Subgrid[lev]->GridData, TopGrid, MetaData)
                == FAIL) {
              ENZO_FAIL("Error in RotatingCylinderInitialize[Sub]Grid.");
            }

        }
        else{
          printf("RotatingCylinder: single grid start-up.\n");
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

      if (TestProblemData.UseMetallicityField)
        DataLabel[i++] = MetalName;

      for(j=0; j < i; j++)
        DataUnits[j] = NULL;

      /* Write parameters to parameter output file */

      if (MyProcessorNumber == ROOT_PROCESSOR) {
        fprintf(Outfptr, "RotatingCylinderOverdensity         = %"FSYM"\n"  , RotatingCylinderOverdensity);
        fprintf(Outfptr, "RotatingCylinderLambda         = %"FSYM"\n"  , RotatingCylinderLambda);
        fprintf(Outfptr, "RotatingCylinderTotalEnergy         = %"FSYM"\n"  , RotatingCylinderTotalEnergy);
        fprintf(Outfptr, "RotatingCylinderRadius         = %"PSYM"\n"  , RotatingCylinderRadius);
        fprintf(Outfptr, "RotatingCylinderCenterPosition = %"PSYM" %"PSYM" %"PSYM"\n",
            RotatingCylinderCenterPosition, RotatingCylinderCenterPosition+1,
            RotatingCylinderCenterPosition+2);
        fprintf(Outfptr, "TestProblemUseMetallicityField  = %"ISYM"\n", TestProblemData.UseMetallicityField);
        fprintf(Outfptr, "TestProblemInitialMetallicityFraction  = %"FSYM"\n", TestProblemData.MetallicityField_Fraction);

      } //   if (MyProcessorNumber == ROOT_PROCESSOR) 


      if(debug){
        printf("Exiting RotatingCylinderInitialize\n");
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

      RotatingCylinderGrid *thisgrid =
        static_cast<RotatingCylinderGrid *>(thisgrid_orig);

      if (thisgrid->ProcessorNumber != MyProcessorNumber)
        return SUCCESS;

      if(debug){
        printf("Entering RotatingCylinderInitializeGrid\n");
        fflush(stdout);
      }

      printf("RotatingCylinderRadius = %e\n", this->RotatingCylinderRadius);
      printf("RotatingCylinderCenterPosition = %e %e %e\n", 
          this->RotatingCylinderCenterPosition[0],
          this->RotatingCylinderCenterPosition[1],
          this->RotatingCylinderCenterPosition[2]);
      printf("RotatingCylinderLambda = %e\n",this->RotatingCylinderLambda);
      printf("RotatingCylinderOverdensity = %e\n",this->RotatingCylinderOverdensity);


      /* declarations */

      int size = 1, dim, cellindex;

      for (dim = 0; dim < thisgrid->GridRank; dim++)
        size *= thisgrid->GridDimension[dim];

      FLOAT r,x,y,z, radius, zdist;

      float sintheta, costheta, omega;

      int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, MetalNum;

      if (thisgrid->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
            Vel3Num, TENum) == FAIL) {
        fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
        ENZO_FAIL("");
      }

      int MetallicityField = FALSE;
      if ((MetalNum = FindField(Metallicity, thisgrid->FieldType, thisgrid->NumberOfBaryonFields))
          != -1)
        MetallicityField = TRUE;
      else
        MetalNum = 0;

      /* set fields in the cylinder region */

      int index, jndex, i, j, k;
      float outside_rho, outside_TE, outside_GE;

      outside_rho =  thisgrid->BaryonField[DensNum][0];

      omega = RotatingCylinderLambda * sqrt(GravitationalConstant * RotatingCylinderOverdensity * outside_rho) / 0.117;

      if(HydroMethod==2){  // ZEUS

        outside_TE = thisgrid->BaryonField[TENum][0];

      } else { // PPM

        outside_TE = thisgrid->BaryonField[TENum][0];

        if(DualEnergyFormalism){
          outside_GE = thisgrid->BaryonField[GENum][0];
        }

      }  // if(HydroMethod==2)

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
            radius = POW(x-RotatingCylinderCenterPosition[0], 2.0) +
              POW(y-RotatingCylinderCenterPosition[1], 2.0);

            radius = sqrt(radius);  // ok, now it's just radius

            zdist = fabs(z-RotatingCylinderCenterPosition[2]);

            if ( (radius <= RotatingCylinderRadius) && (zdist <= RotatingCylinderRadius) ){

              thisgrid->BaryonField[DensNum][cellindex] = outside_rho * RotatingCylinderOverdensity;

              if(TestProblemData.UseMetallicityField>0 && MetalNum != FALSE)
                thisgrid->BaryonField[MetalNum][cellindex] = thisgrid->BaryonField[DensNum][cellindex]*TestProblemData.MetallicityField_Fraction;

              sintheta = (y-RotatingCylinderCenterPosition[1])/radius;
              costheta = (x-RotatingCylinderCenterPosition[0])/radius;


              // x,y, and maybe z velocity.  
              thisgrid->BaryonField[Vel1Num][cellindex] = -1.0*sintheta*omega*radius;

              thisgrid->BaryonField[Vel2Num][cellindex] = costheta*omega*radius;

              thisgrid->BaryonField[Vel3Num][cellindex] = 0.0;

              if(HydroMethod == 2){

                // ZEUS
                thisgrid->BaryonField[TENum][cellindex] = outside_TE / RotatingCylinderOverdensity;

              } else {

                // PPM
                thisgrid->BaryonField[TENum][cellindex] = outside_TE / RotatingCylinderOverdensity
                  + 0.5 * thisgrid->BaryonField[Vel1Num][cellindex] *
                  thisgrid->BaryonField[Vel1Num][cellindex]
                  + 0.5 * thisgrid->BaryonField[Vel2Num][cellindex] *
                  thisgrid->BaryonField[Vel2Num][cellindex]
                  + 0.5 * thisgrid->BaryonField[Vel3Num][cellindex] *
                  thisgrid->BaryonField[Vel3Num][cellindex];

                // gas energy (PPM dual energy formalims)
                if(DualEnergyFormalism)
                  thisgrid->BaryonField[GENum][cellindex] = outside_GE / RotatingCylinderOverdensity;

              } // if(HydroMethod == 2)

            } // if (r <= RotatingCylinderRadius)

          } // for (i = 0; i < GridDimension[0]; i++)

      if(debug){
        printf("Exiting RotatingCylinderInitialize\n");
        fflush(stdout);
      }

      return SUCCESS;

    } 


};

//.. register:
namespace{
    EnzoProblemType_creator_concrete<ProblemType_RotatingCylinder>
        rotating_cylinder("RotatingCylinder");
}

#endif
