
/***********************************************************************
/
/  Determine the types of active particles we need to create to 
/  replace the star objects. 
/
/  written by: John Regan
/  date:       May, 2018
/  modified1:  
/
/  PURPOSE:
/
************************************************************************/
#include "preincludes.h"
#include "ErrorExceptions.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "CosmologyParameters.h"
#include "ActiveParticle.h"

int grid::DetermineActiveParticleTypes(char **ActiveParticleType)
{


  bool present = false;
  int type = 0, numtypes = 0;
  
  for(int ptype = 0; ptype <=  MAX_ACTIVE_PARTICLE_TYPES; ptype++)
    {
      present = false;
      switch(ptype) {

      case PARTICLE_TYPE_SINGLE_STAR:
	for(int part = 0; part < NumberOfParticles; part++) {
	  if(ParticleType == NULL || ParticleType[part] != ptype)
	    continue;
	  if(ParticleType[part] == ptype){
	    present = true;
	    numtypes++;
	  }
	}
	if(present == true){
	  strcpy(ActiveParticleType[ptype], "PopIII");
	  SetActiveParticleTypeCounts(ptype, numtypes);
	}
	break;

      case PARTICLE_TYPE_CLUSTER:
	for(int part = 0; part < NumberOfParticles; part++) {
	   if(ParticleType == NULL || ParticleType[part] != ptype)
	    continue;
	  if(ParticleType[part] == ptype) {
	    present = true;
	    numtypes++;
	  }
	}
	if(present == true) {
	  strcpy(ActiveParticleType[ptype], "CenOstriker");
	  SetActiveParticleTypeCounts(ptype, numtypes);
	}
	break;
      }
      
    }
  
  return numtypes;
}
