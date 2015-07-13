/***********************************************************************
/
/  SET THE PARAMETER DEFAULTS
/
/  written by: Greg Bryan
/  date:       June, 1997
/  modified1:  Robert Harkness
/  date:       November, 2003
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
// This routine reads the parameter file in the argument
// and sets parameters based upon it.
 
#include <stdio.h>
#include <string.h>
#include <math.h>
 
#include "macros_and_parameters.h"
#include "global_data.h"
#include "Parameters.h"
 
// function prototypes
 
 
 
 
int SetParameterDefaults(parmstruct *Parameters)
{
 
  int dim;
 
  char ppos_name[] = "ParticlePositions";
  char pvel_name[] = "ParticleVelocities";
  char gden_name[] = "GridDensity";
  char gvel_name[] = "GridVelocity";
 
  // Set Defaults
 
  Parameters->Rank                = 3;
  Parameters->GridRefinement      = 1;
  Parameters->ParticleRefinement  = 1;
  Parameters->WaveNumberCutoff    = INT_UNDEFINED;
  Parameters->InitializeParticles = TRUE;
  Parameters->InitializeGrids     = TRUE;
  Parameters->RandomNumberGenerator = 0;
 
  Parameters->ParticlePositionName = ppos_name;
  Parameters->ParticleVelocityName = pvel_name;
  Parameters->ParticleMassName     = NULL;
  Parameters->ParticleTypeName     = NULL;
  Parameters->GridDensityName      = gden_name;
  Parameters->GridVelocityName     = gvel_name;
 
  for (dim = 0; dim < Parameters->Rank; dim++) {
    Parameters->GridDims[dim]     = INT_UNDEFINED;
    Parameters->ParticleDims[dim] = INT_UNDEFINED;
    Parameters->MaxDims[dim]      = INT_UNDEFINED;
    Parameters->NewCenter[dim]    = INT_UNDEFINED;
    Parameters->NewCenterFloat[dim] = FLOAT_UNDEFINED;
    Parameters->StartIndex[dim]   = 0;
    Parameters->StartIndexInNewCenterTopGridSystem[dim] = INT_UNDEFINED;
    Parameters->EndIndexInNewCenterTopGridSystem[dim] = INT_UNDEFINED;
    Parameters->TopGridStart[dim] = INT_UNDEFINED;
    Parameters->TopGridEnd[dim] = INT_UNDEFINED;
    Parameters->RootGridDims[dim] = INT_UNDEFINED;
    Parameters->RefineRegionLeftEdge[dim] = FLOAT_UNDEFINED;
    Parameters->RefineRegionRightEdge[dim] = FLOAT_UNDEFINED;
  }

  Parameters->RefineBy = 2;
  Parameters->MaximumInitialRefinementLevel = INT_UNDEFINED;
  Parameters->AutomaticSubgridBuffer = 1;
 
  return SUCCESS;
}
