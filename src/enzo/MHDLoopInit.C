#include <math.h>
#include <string.h>
#include <stdio.h>

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
int MHDLoopInit(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
                          TopGridData &MetaData, ExternalBoundary &Exterior)
{
  float Pi = 3.14159265;
  float Density = 25.0/(36*Pi); 
  float Pressure = 5.0/(12*Pi);
  float Vx=0.0, Vy = 0.0, Vz = 0.0;
  float B0=1e-3;
  FLOAT R0 = 0.3;
  FLOAT Center[3];
  int CurrentAxis = 2;
  for(int dim =0;dim<3;dim++)
    Center[dim] = 0.5*(DomainLeftEdge[dim]+DomainRightEdge[dim]);

  char *DensName = "Density";
  char *TEName = "TotalEnergy";
  char *GEName   = "GasEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  char *BxName = "Bx";
  char *ByName = "By";
  char *BzName = "Bz";
  char *PhiName = "Phi";
  int i=0;
  
    DataLabel[i++] = DensName;
  if( EquationOfState == 0 ){
    DataLabel[i++] = TEName;
  }
    DataLabel[i++] = Vel1Name;
    DataLabel[i++] = Vel2Name;
    DataLabel[i++] = Vel3Name;
  if( UseMHD ){
    DataLabel[i++] = BxName;
    DataLabel[i++] = ByName;
    DataLabel[i++] = BzName;
    DataLabel[i++] = PhiName;
  }
   if (DualEnergyFormalism)
    DataLabel[i++] = GEName;

  for ( int j=0;j<i;j++)
    DataUnits[j] = NULL;

  MHDcLabel[0] = "Bx";
  MHDcLabel[1] = "By";
  MHDcLabel[2] = "Bz";

  MHDLabel[0] = "BxF";
  MHDLabel[1] = "ByF";
  MHDLabel[2] = "BzF";

  MHDeLabel[0] = "Ex";
  MHDeLabel[1] = "Ey";
  MHDeLabel[2] = "Ez";
  

  MHDUnits[0] = "FourPiGauss";
  MHDUnits[1] = "FourPiGauss";
  MHDUnits[2] = "FourPiGauss";

  MHDeUnits[0] = "FourPiGauss";
  MHDeUnits[1] = "FourPiGauss";
  MHDeUnits[2] = "FourPiGauss";

  MHDUnits[0] = "None";
  MHDUnits[1] = "None";
  MHDUnits[2] = "None";
  
  MHDeUnits[0] = "None";
  MHDeUnits[1] = "None";
  MHDeUnits[2] = "None";

  
  char line[MAX_LINE_LENGTH];
  int ret=0;

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
    ret = 0;
    ret += sscanf(line, "MHDLoopDensity = %"PSYM, &Density);
    ret += sscanf(line, "MHDLoopPressure = %"PSYM, &Pressure);
    ret += sscanf(line, "MHDLoopVx = %"PSYM, &Vx);
    ret += sscanf(line, "MHDLoopVy = %"PSYM, &Vy);
    ret += sscanf(line, "MHDLoopVz = %"PSYM, &Vz);
    ret += sscanf(line, "MHDLoopB0 = %"PSYM, &B0);
    ret += sscanf(line, "MHDLoopR0 = %"PSYM, &R0);
    ret += sscanf(line, "MHDLoopCurrentAxis = %d", &CurrentAxis);
    ret += sscanf(line, "MHDLoopCenter = %"PSYM" %"PSYM" %"PSYM, 
		  Center, Center+1, Center+2);
    
  }  
  if( TopGrid.GridData->MHDLoopInitGrid(Density,Pressure,Vx,Vy, Vz, B0, R0,Center, CurrentAxis) == FAIL ){
    ENZO_FAIL("Error in MHDLoopInitGrid\n");
  }

  return SUCCESS;


}




