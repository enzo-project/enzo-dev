/***********************************************************************
/
/  READS THE STAR PARTICLE DATA
/
/  written by: Greg Bryan
/  date:       March, 1997
/  modified1: Michael Kuhlen, October 2010, HDF5 hierarchy
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
#include <hdf5.h>
#include <stdio.h>
#include <string.h>


#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "StarParticleData.h"

// defined in Grid_ReadHierarchyInformationHDF5.C 
int HDF5_ReadAttribute(hid_t group_id, const char *AttributeName, int &Attribute, FILE *log_fptr);

int ReadStarParticleData(FILE *fptr, hid_t Hfile_id, FILE *log_fptr)
{
 
  if (StarParticleCreation == FALSE)
    return SUCCESS;
 
  /* read in number data. */

  if (!AddParticleAttributes) {

    if (HierarchyFileInputFormat % 2 == 0) {
      HDF5_ReadAttribute(Hfile_id, "NumberOfStarParticles", NumberOfStarParticles, log_fptr);
    }

    if (HierarchyFileInputFormat == 1) {
      if (fscanf(fptr, "NumberOfStarParticles = %"ISYM"\n",
		 &NumberOfStarParticles) != 1) {
	//      ENZO_FAIL("Error reading NumberOfStarParticles.\n");
	
      }
    }

  } else 
    NumberOfStarParticles = 0;
  
  NumberOfOtherParticles = 0; 

  return SUCCESS;
}
