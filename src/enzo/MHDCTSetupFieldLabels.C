/***********************************************************************
/
/  MHDCTSetupFieldLabels
/
/  written by: Sam Skillman
/  date:       September, 2012
/  modified1: 
/
/  PURPOSE: Helper function to set MHD Labels when reading data.
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/


#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

void MHDCTSetupFieldLabels(void){
  if (!UseMHDCT)
    return;

  if(MHDLabel[0]==NULL)
    MHDLabel[0] = "BxF";
  if(MHDLabel[1]==NULL)
    MHDLabel[1] = "ByF";
  if(MHDLabel[2]==NULL)
    MHDLabel[2] = "BzF";

  if(MHDUnits[0]==NULL)
    MHDUnits[0] = "None";
  if(MHDUnits[1]==NULL)
    MHDUnits[1] = "None";
  if(MHDUnits[2]==NULL)
    MHDUnits[2] = "None";
  if(MHDeLabel[0] == NULL)

    MHDeLabel[0] = "Ex";
  if(MHDeLabel[1] == NULL)
    MHDeLabel[1] = "Ey";
  if(MHDeLabel[2] == NULL)
    MHDeLabel[2] = "Ez";
  if(MHDeUnits[0] == NULL)

    MHDeUnits[0] = "None";
  if(MHDeUnits[1] == NULL)
    MHDeUnits[1] = "None";
  if(MHDeUnits[2] == NULL)
    MHDeUnits[2] = "None";

  return;
}




