/***********************************************************************
/
/  GRID CLASS (BEFORE REBUILDING, REMOVED UNNEEDED ARRAYS)
/
/  written by: Greg Bryan
/  date:       June, 1995
/  modified1:
/
/  PURPOSE:
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
 
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
/* function prototypes */
 
 
void grid::CleanUp()
{
 
  int i;
 
  for (i = 0; i < MAX_DIMENSION; i++) {
    delete [] ParticleAcceleration[i];
//    delete [] AccelerationField[i];
 
    ParticleAcceleration[i]      = NULL;
//    AccelerationField[i]         = NULL;
  }
  delete [] ParticleAcceleration[MAX_DIMENSION];
  ParticleAcceleration[MAX_DIMENSION] = NULL;
 
  for (i = 0; i < MAX_NUMBER_OF_BARYON_FIELDS; i++) {
    delete [] OldBaryonField[i];
    OldBaryonField[i] = NULL;
  }
 
  delete [] GravitatingMassField;
  delete [] GravitatingMassFieldParticles;
 
  GravitatingMassField          = NULL;
  GravitatingMassFieldParticles = NULL;

#ifdef SAB
  for (i = 0; i < MAX_DIMENSION; i++)
    if (OldAccelerationField[i] != NULL) {
      delete [] OldAccelerationField[i];
      OldAccelerationField[i] = NULL;
    }
#endif

  if( UseMHDCT ){
    for(i=0;i<3;i++){
      if( OldMagneticField[i] != NULL ){
	delete [] OldMagneticField[i];
	OldMagneticField[i] = NULL;
      }

      
    }
  }

}
