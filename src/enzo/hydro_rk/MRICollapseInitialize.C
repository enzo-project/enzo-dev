#include <string.h>
#include <stdio.h>
#include <math.h>
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
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);
int GetUnits(float *DensityUnits, float *LengthUnits,
		      float *TemperatureUnits, float *TimeUnits,
		      float *VelocityUnits, FLOAT Time);

int MRICollapseInitialize(FILE *fptr, FILE *Outfptr, 
			    HierarchyEntry &TopGrid, TopGridData &MetaData)
{
  char *DensName = "Density";
  char *TEName   = "TotalEnergy";
  char *GEName   = "GasEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  char *BxName = "Bx";
  char *ByName = "By";
  char *BzName = "Bz";
  char *PhiName = "Phi";
  char *DebugName = "Debug";
  char *Phi_pName = "Phip";

  /* declarations */

  char  line[MAX_LINE_LENGTH];
  int   dim, ret, level, sphere, i;

  /* set default parameters */

 
  /* read input from file */

  float AngularVelocity=1e-3;
  float VelocityGradient=1.5;
  float ThermalMagneticRatio=400; 
  float FluctuationAmplitudeFraction=0.1;
  float Radius=0.25;
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    /* read parameters */

    ret += sscanf(line, "AngularVelocity= %f", &AngularVelocity);
    ret += sscanf(line, "VelocityGradient= %f", &VelocityGradient);
    ret += sscanf(line, "ThermalMagneticRatio= %f", &ThermalMagneticRatio);
    ret += sscanf(line, "FluctuationAmplitudeFraction = %f", &FluctuationAmplitudeFraction);
    ret += sscanf(line, "Radius = %f", &Radius);
   

  } // end input from parameter file
  
  //  float rhou = 1.0, lenu = 1.0, tempu = 1.0, tu = 1.0, velu = 1.0, presu = 1.0, bfieldu = 1.0;
 //  if (UsePhysicalUnit) {
//     GetUnits(&rhou, &lenu, &tempu, &tu, &velu, MetaData.Time);
//     presu = rhou*lenu*lenu/tu/tu;
//     bfieldu = sqrt(presu*4.0*Pi);
//   }
  
//  printf("rhou=%g,velu=%g,lenu=%g,tu=%g (%g yr),tempu=%g,presu=%g, bfieldu=%g\n", 
//	 rhou, velu,lenu,tu,tu/3.1558e7,tempu,presu,bfieldu);

  
  // B0 /= bfieldu;

  //printf("t=%g\n", MediumPressure/MediumDensity*tempu);



  if (TopGrid.GridData->MRICollapseInitializeGrid(AngularVelocity, VelocityGradient, ThermalMagneticRatio, FluctuationAmplitudeFraction, Radius)
      == FAIL) {
    fprintf(stderr, "Error in ShearingBoxInitializeGrid.\n");
    return FAIL;
  }


  /* set up field names and units */

  int count = 0;
  DataLabel[count++] = DensName;
  DataLabel[count++] = Vel1Name;
  DataLabel[count++] = Vel2Name;
  DataLabel[count++] = Vel3Name;
  DataLabel[count++] = TEName;
  if (DualEnergyFormalism) {
    DataLabel[count++] = GEName;
  }
  DataLabel[count++] = BxName;
  DataLabel[count++] = ByName;
  DataLabel[count++] = BzName;
  DataLabel[count++] = PhiName;
  if(UseDivergenceCleaning){
    DataLabel[count++] = Phi_pName;
    DataLabel[count++] = DebugName;
  }

  for (i = 0; i < count; i++) {
    DataUnits[i] = NULL;
  }

  return SUCCESS;

}
