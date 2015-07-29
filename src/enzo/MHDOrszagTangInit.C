//
/***********************************************************************
/
/  GRID CLASS ()
/
/  written by: David Collins
/  date:       2005?
/  modified1:
/
/  PURPOSE: The Orszag Tang Vortex.
/           ICs from Ryu, Dongsu; Miniati, Francesco; Jones, T. W.; Frank, Adam, 1998 
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/


#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"

int MHDOrszagTangInit(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
		      TopGridData &MetaData, ExternalBoundary &Exterior){

  float Pi = 3.14159265;
  float Density = 25.0/(36*Pi); 
  float Pressure = 5.0/(12*Pi);
  float V0=1.0;
  float B0=1/sqrt(4*Pi);
  //The Isothermal initial conditions come from Mignone, JCP, 2007
  if( EquationOfState == 1 ){
    Density = 1;
    V0= IsothermalSoundSpeed ;
    B0 = IsothermalSoundSpeed * sqrt(3./5.);
  }



  int count=0;
  char *DensName = "Density";
  
  char *TEName = "TotalEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  char *GEName   = "GasEnergy";
  char *BxName = "Bx";
  char *ByName = "By";
  char *BzName = "Bz";
  char *PhiName = "Phi";

  DataLabel[count++] = DensName;
  if( EquationOfState == 0 )  
    DataLabel[count++] = TEName;
  if (DualEnergyFormalism){
    DataLabel[count++] = GEName;
    //<dbg>
    fprintf(stderr,"MHDORszagTangInit: Dual Energy Formailism not completely installed.\n");
    return FAIL;
  }  
  DataLabel[count++] = Vel1Name;
  DataLabel[count++] = Vel2Name;
  DataLabel[count++] = Vel3Name;
  if( UseMHD ){
      DataLabel[count++] = BxName;
      DataLabel[count++] = ByName;
      DataLabel[count++] = BzName;
  }
  if (HydroMethod == MHD_RK){
      DataLabel[count++] = PhiName;
  }

  for (int i = 0; i < count; i++)
    DataUnits[i] = NULL;

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


  if( TopGrid.GridData->MHDOrszagTangInitGrid(Density,Pressure,V0,B0) == FAIL ){
    fprintf(stderr, " Error in Orszag Tang\n");
    return FAIL;
  }

  return SUCCESS;

}
