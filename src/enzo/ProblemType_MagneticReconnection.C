/***********************************************************************
/
/  MAGNETIC RECONNECTION PROBLEM TYPE
/
/  written by: J. S. Oishi, Matthew Turk, Brian O'Shea
/  date:       March, 2013
/
/  PURPOSE: to reconnect field lines
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

class ProblemType_MagneticReconnection;

class MagneticReconnectionGrid : private grid {
    friend class ProblemType_MagneticReconnection;
};

void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);
int GetUnits (float *DensityUnits, float *LengthUnits,
	      float *TemperatureUnits, float *TimeUnits,
	      float *VelocityUnits, double *MassUnits, FLOAT Time);

int FindField(int field, int farray[], int numfields);

class ProblemType_MagneticReconnection : public EnzoProblemType
{
    private:
        FLOAT MagneticReconnectionSubgridLeft, MagneticReconnectionSubgridRight;
        FLOAT LeftEdge[MAX_DIMENSION], RightEdge[MAX_DIMENSION];
        FLOAT MagneticReconnectionCenterPosition[MAX_DIMENSION];
        float MagneticReconnectionVelocity[3];   // gas initally at rest
        float MagneticReconnectionBField[3];   
        float MagneticReconnectionBperturbk[MAX_DIMENSION];   // perturbation wavenumbers
        float MagneticReconnectionLambda;
        float MagneticReconnectionBperturbation;
        float MagneticReconnectionOverdensity;
        float MagneticReconnectionDensity;
        float MagneticReconnectionTotalEnergy;
        int MagneticReconnectionRefineAtStart;

    public:
    ProblemType_MagneticReconnection() : EnzoProblemType()
    { 
        std::cout << "Creating problem type Magnetic Reconnection" << std::endl;
        //RegisterEventPlugin("Printing", &JustPrintSomething);
        //RegisterEventHook("EvolveLevelTop", "Printing");
    }

    ~ProblemType_MagneticReconnection()
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
      char *DensName = "Density";
      char *TEName   = "TotalEnergy";
      char *GEName   = "GasEnergy";
      char *Vel1Name = "x-velocity";
      char *Vel2Name = "y-velocity";
      char *Vel3Name = "z-velocity";
      char *MetalName = "Metal_Density";
      char *BxName = "Bx";
      char *ByName = "By";
      char *BzName = "Bz";
      char *PhiName = "Phi";
      if ( UseMHDCT ){
         MHDLabel[0] = "BxF";
         MHDLabel[1] = "ByF";
         MHDLabel[2] = "BzF";
         
         MHDeLabel[0] = "Ex";
         MHDeLabel[1] = "Ey";
         MHDeLabel[2] = "Ez";
         
         MHDUnits[0] = "None";
         MHDUnits[1] = "None";
         MHDUnits[2] = "None";
         
         MHDeUnits[0] = "None";
         MHDeUnits[1] = "None";
         MHDeUnits[2] = "None";
      }
      

      /* local declarations */

      char line[MAX_LINE_LENGTH];
      int  i, j, dim, ret, level, NumberOfSubgridZones[MAX_DIMENSION],
           SubgridDims[MAX_DIMENSION];
      float Pi                      = 3.14159;

      for(i=0; i<MAX_DIMENSION; i++) 
        MagneticReconnectionCenterPosition[i] = 0.5;  // right in the middle of the box

      MagneticReconnectionBperturbk[0] = 2 * Pi;  
      MagneticReconnectionBperturbk[1] = Pi;  
      MagneticReconnectionBperturbk[2] = 0;


      this->MagneticReconnectionVelocity[0] = 
        this->MagneticReconnectionVelocity[1] = 
        this->MagneticReconnectionVelocity[2] = 0.0; // gas initally at rest
      this->MagneticReconnectionBField[0] = 1.0;
      this->MagneticReconnectionBField[1] =
        this->MagneticReconnectionBField[2] = 0.0; // gas initally at rest
      this->MagneticReconnectionLambda = 0.5; // in problem units!
      this->MagneticReconnectionBperturbation = 0.1; // in problem units!
      this->MagneticReconnectionOverdensity = 1.;
      this->MagneticReconnectionDensity = 0.2;
      this->MagneticReconnectionTotalEnergy = 0.5;
      this->MagneticReconnectionRefineAtStart = TRUE;

      /* set no subgrids by default. */

      MagneticReconnectionSubgridLeft         = 0.0;    // start of subgrid(s)
      MagneticReconnectionSubgridRight        = 0.0;    // end of subgrid(s)

      /* read input from file */

      while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

        ret = 0;

        /* read parameters specifically for radiating shock problem*/

        ret += sscanf(line, "MagneticReconnectionOverdensity  = %"FSYM, &MagneticReconnectionOverdensity);
        ret += sscanf(line, "MagneticReconnectionSubgridLeft = %"PSYM,
            &MagneticReconnectionSubgridLeft);
        ret += sscanf(line, "MagneticReconnectionSubgridRight = %"PSYM,
            &MagneticReconnectionSubgridRight);
        ret += sscanf(line, "MagneticReconnectionLambda = %"FSYM,
            &MagneticReconnectionLambda);
        ret += sscanf(line, "MagneticReconnectionBperturbation = %"FSYM,
            &MagneticReconnectionBperturbation);

        ret += sscanf(line, "MagneticReconnectionTotalEnergy = %"FSYM,
            &MagneticReconnectionTotalEnergy);

        ret += sscanf(line, "MagneticReconnectionCenterPosition = %"PSYM" %"PSYM" %"PSYM,
            MagneticReconnectionCenterPosition, MagneticReconnectionCenterPosition+1,
            MagneticReconnectionCenterPosition+2);
        ret += sscanf(line, "MagneticReconnectionBperturbk = %"PSYM" %"PSYM" %"PSYM,
            MagneticReconnectionBperturbk, MagneticReconnectionBperturbk+1,
            MagneticReconnectionBperturbk+2);
        ret += sscanf(line, "MagneticReconnectionBField = %"PSYM" %"PSYM" %"PSYM,
            MagneticReconnectionBField, MagneticReconnectionBField+1,
            MagneticReconnectionBField+2);
        ret += sscanf(line, "MagneticReconnectionRefineAtStart = %"ISYM, &MagneticReconnectionRefineAtStart);

        ret += sscanf(line, "TestProblemUseMetallicityField  = %"ISYM, &TestProblemData.UseMetallicityField);
        ret += sscanf(line, "TestProblemInitialMetallicityFraction  = %"FSYM, &TestProblemData.MetallicityField_Fraction);

        /* if the line is suspicious, issue a warning */

        if (ret == 0 && strstr(line, "=") && (strstr(line, "MagneticReconnection") || strstr(line, "TestProblem")) &&
            line[0] != '#' && MyProcessorNumber == ROOT_PROCESSOR)
          fprintf(stderr,
              "*** warning: the following parameter line was not interpreted:\n%s\n",
              line);

      } // end input from parameter file



      MHD_ProjectE = FALSE;
      MHD_ProjectB = TRUE;
      this->InitializeUniformGrid(TopGrid.GridData,
            MagneticReconnectionDensity,
            MagneticReconnectionTotalEnergy,
            MagneticReconnectionTotalEnergy,
            MagneticReconnectionVelocity,
            MagneticReconnectionBField);
      this->InitializeGrid(TopGrid.GridData, TopGrid, MetaData);
      if (MagneticReconnectionRefineAtStart)
         {
            /* Declare, initialize, and fill out the LevelArray. */
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
                  if (this->InitializeGrid(Temp->GridData, TopGrid, MetaData) == FAIL)
                     {
                        ENZO_FAIL("Error in MagneticReconnection->InitializeGrid");
                     }
                  Temp = Temp->NextGridThisLevel;
               } // end: loop over grids on this level
            } // end: loop over levels
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

      if(MetaData.TopGridRank > 1 || UseMHD )
        DataLabel[i++] = Vel2Name;

      if(MetaData.TopGridRank > 2 || UseMHD )
        DataLabel[i++] = Vel3Name;

      if( UseMHD ){
          DataLabel[i++] = BxName;
          DataLabel[i++] = ByName;
          DataLabel[i++] = BzName;
      }
      if (HydroMethod == MHD_RK){
          DataLabel[i++] = PhiName;
      }

      if (TestProblemData.UseMetallicityField)
        DataLabel[i++] = MetalName;

      for(j=0; j < i; j++)
        DataUnits[j] = NULL;

      /* Write parameters to parameter output file */

      if (MyProcessorNumber == ROOT_PROCESSOR) {
        fprintf(Outfptr, "MagneticReconnectionOverdensity         = %"FSYM"\n"  , MagneticReconnectionOverdensity);
        fprintf(Outfptr, "MagneticReconnectionLambda         = %"FSYM"\n"  , MagneticReconnectionLambda);
        fprintf(Outfptr, "MagneticReconnectionBperturbation  = %"FSYM"\n"  , MagneticReconnectionBperturbation);
        fprintf(Outfptr, "MagneticReconnectionTotalEnergy         = %"FSYM"\n"  , MagneticReconnectionTotalEnergy);
        fprintf(Outfptr, "MagneticReconnectionCenterPosition = %"PSYM" %"PSYM" %"PSYM"\n",
            MagneticReconnectionCenterPosition, MagneticReconnectionCenterPosition+1,
            MagneticReconnectionCenterPosition+2);
        fprintf(Outfptr, "MagneticReconnectionBperturbk = %"PSYM" %"PSYM" %"PSYM"\n",
            MagneticReconnectionBperturbk, MagneticReconnectionBperturbk+1,
            MagneticReconnectionBperturbk+2);
        fprintf(Outfptr, "MagneticReconnectionBField = %"PSYM" %"PSYM" %"PSYM"\n",
            MagneticReconnectionBField, MagneticReconnectionBField+1,
            MagneticReconnectionBField+2);
        fprintf(Outfptr, "MagneticReconnectionRefineAtStart           = %"ISYM"\n", MagneticReconnectionRefineAtStart);
        fprintf(Outfptr, "TestProblemUseMetallicityField  = %"ISYM"\n", TestProblemData.UseMetallicityField);
        fprintf(Outfptr, "TestProblemInitialMetallicityFraction  = %"FSYM"\n", TestProblemData.MetallicityField_Fraction);

      } //   if (MyProcessorNumber == ROOT_PROCESSOR) 

      MHD_ProjectE = TRUE;
      MHD_ProjectB = FALSE;

      if(debug){
        printf("Exiting MagneticReconnectionInitialize\n");
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

      MagneticReconnectionGrid *thisgrid =
        static_cast<MagneticReconnectionGrid *>(thisgrid_orig);

      if (thisgrid->ProcessorNumber != MyProcessorNumber)
        return SUCCESS;

      /* declarations */

      int size = 1, dim, cellindex;

      for (dim = 0; dim < thisgrid->GridRank; dim++)
        size *= thisgrid->GridDimension[dim];

      FLOAT r, x, y, z, xy, yx;


      float TemperatureUnits = 1.0, DensityUnits = 1.0, LengthUnits = 1.0;
      float VelocityUnits = 1.0, TimeUnits = 1.0;
      double MassUnits = 1.0;
      
      // Get system of units
      if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, 
                   &TimeUnits, &VelocityUnits, &MassUnits, thisgrid->Time) == FAIL) {
         ENZO_FAIL("Error in GetUnits.");
      }

      float lambda_scaled = MagneticReconnectionLambda/LengthUnits;

      int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum, B1Num, B2Num, B3Num, MetalNum;
      thisgrid->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                           Vel3Num, TENum, B1Num, B2Num, B3Num);

      int MetallicityField = FALSE;
      if ((MetalNum = FindField(Metallicity, thisgrid->FieldType, thisgrid->NumberOfBaryonFields))
          != -1)
        MetallicityField = TRUE;
      else
        MetalNum = 0;

      /* set fields in the cylinder region */

      int index, jndex, i, j, k;
      float outside_rho, outside_TE, outside_GE, kx, ky;
      kx = MagneticReconnectionBperturbk[0];
      ky = MagneticReconnectionBperturbk[1];

      outside_rho =  thisgrid->BaryonField[DensNum][0];

      if(HydroMethod==2){  // ZEUS

        outside_TE = thisgrid->BaryonField[TENum][0];

      } else { // PPM

        outside_TE = thisgrid->BaryonField[TENum][0];

        if(DualEnergyFormalism){
          outside_GE = thisgrid->BaryonField[GENum][0];
        }

      }  // if(HydroMethod==2)

  
      // initialize B fields
      int field;
      for (field = 0; field < 3; field++)
         for (k = 0; k < thisgrid->MagneticDims[field][2]; k++)
            for (j = 0; j < thisgrid->MagneticDims[field][1]; j++)
               for (i = 0; i < thisgrid->MagneticDims[field][0]; i++){

                  /* Compute position */
                  x=y=z=0.0;
                  
                  cellindex = i + j*thisgrid->MagneticDims[field][0]
                     + k*thisgrid->MagneticDims[field][0]
                     *thisgrid->MagneticDims[field][1];
                  
                  // for face centered fields, use functions f,g of position like so:
                  // Bx = f(x,xy)
                  // By = g(yx,y)
                  x = thisgrid->CellLeftEdge[0][i];
                  y = thisgrid->CellLeftEdge[1][j];
                  xy = thisgrid->CellLeftEdge[0][i] + 0.5*thisgrid->CellWidth[0][i];
                  yx = thisgrid->CellLeftEdge[1][j] + 0.5*thisgrid->CellWidth[1][j];

                  if (thisgrid->GridRank == 3)
                     z = thisgrid->CellLeftEdge[2][k] + k*thisgrid->CellWidth[2][k];
                  
                  switch (field) {
                  case 0:
                    thisgrid->MagneticField[field][cellindex] = MagneticReconnectionBField[0] * tanh((yx-MagneticReconnectionCenterPosition[1])/lambda_scaled) + MagneticReconnectionBperturbation * ky * cos(kx*LengthUnits * x) * sin(ky*LengthUnits * (yx-MagneticReconnectionCenterPosition[1]));
                     break;
                  case 1:
                    thisgrid->MagneticField[field][cellindex] = MagneticReconnectionBField[1] -MagneticReconnectionBperturbation * kx * sin(kx*LengthUnits * xy) * cos(ky*LengthUnits * (y-MagneticReconnectionCenterPosition[1]));
                     break;
                  case 2:
                     thisgrid->MagneticField[field][cellindex] = MagneticReconnectionBField[2];
                     break;
                  }
               }


      if( thisgrid->CenterMagneticField() == FAIL ) {
         fprintf(stderr," error with CenterMagneticField\n");
         return FAIL;
      }

      // Thermal Physics compensates for the field profile with a density bump
      // and an assumed isothermal initial condition

      // Use the MagneticReconnectionTotalEnergy as an isothermal sound speed...
      MagneticReconnectionTotalEnergy = POW(MagneticReconnectionBField[0],2)/(2.* (MagneticReconnectionOverdensity + MagneticReconnectionDensity));
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
            if (thisgrid->GridRank == 3)
               z = thisgrid->CellLeftEdge[2][k] + 0.5*thisgrid->CellWidth[2][k];

            thisgrid->BaryonField[DensNum][cellindex] = 
                MagneticReconnectionOverdensity/POW(cosh((y-MagneticReconnectionCenterPosition[1])/lambda_scaled),2) 
                + MagneticReconnectionDensity;

            thisgrid->BaryonField[TENum][cellindex] = 
                MagneticReconnectionTotalEnergy/(Gamma - 1.0) 
                        + 0.5 * (POW(thisgrid->BaryonField[B1Num][cellindex], 2) 
                               + POW(thisgrid->BaryonField[B2Num][cellindex], 2) 
                               + POW(thisgrid->BaryonField[B3Num][cellindex], 2))/
                                  thisgrid->BaryonField[DensNum][cellindex];
          } // for (i = 0; i < GridDimension[0]; i++)

      return SUCCESS;

    } 


};

//.. register:
namespace{
    EnzoProblemType_creator_concrete<ProblemType_MagneticReconnection>
        magnetic_reconnection("MagneticReconnection");
}

#endif
