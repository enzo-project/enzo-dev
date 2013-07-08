/**************************
MHDCT_ParameterJuggle.

If the parameter file is from the old vesion of MHDCT, there is some parameter 
collisions.  In order to not have to change every old MHDCT parameter file,
this routine automatically changes the parameters if MHD_Used (exclusively and old flag)
is set to true.  Parameters are written with the new numbering.


**************************/

 
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
#include "StarParticleData.h"

int MHDCT_ParameterJuggle(){

  if( UseMHDCT == FALSE)
    return SUCCESS;

  //
  //Problem Types
  //

  //MHDBlast
  if( ProblemType == 100)
    ProblemType = 500;

  if( HydroMethod == HD_RK || HydroMethod == MHD_RK )
    ENZO_VFAIL("Error: HydroMethod %"ISYM" from Old MHD incompatable with New MHD.\n Either set UseMHDCT = 0, or use HydroMethod 6.",HydroMethod);
  
  return SUCCESS;
}
