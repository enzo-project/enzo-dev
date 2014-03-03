/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE PARTICLES, GIVEN A VELOCITY FIELD)
/
/  written by: John H. Wise
/  date:       January, 2009
/  modified1:  
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include <hdf5.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "h5utilities.h"

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int ReadFile(char *name, int Rank, int Dim[], int StartIndex[],
	     int EndIndex[], int BufferOffset[], float *buffer,
	     inits_type **tempbuffer, int Part, int Npart);

int ReadIntFile(char *name, int Rank, int Dims[], int StartIndex[],
		int EndIndex[], int BufferOffset[], int *buffer,
		int **tempbuffer, int Part, int Npart);

//#define ICPART_SHIFT8

int grid::CosmologyInitializeParticles(
		   char *CosmologySimulationParticleVelocityName,
		   char *CosmologySimulationParticleDisplacementName,
		   char *CosmologySimulationParticleMassName,
		   char *CosmologySimulationParticleTypeName,
		   char *CosmologySimulationParticleVelocityNames[],
		   char *CosmologySimulationParticleDisplacementNames[],
		   float CosmologySimulationOmegaBaryonNow,
		   int *Offset, int level)
{

  int dim, i, j, k, index, index1, index2;
  bool OneComponentPerFile;

  float DensityUnits=1, LengthUnits=1, TemperatureUnits=1, TimeUnits=1,
    VelocityUnits=1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits,
	       InitialTimeInCodeUnits) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }

  /* Check if we have one component per file */

  if (CosmologySimulationParticleVelocityName != NULL &&
      CosmologySimulationParticleVelocityNames[0] != NULL) {
    ENZO_FAIL("grid::CosmologyInitializeParticles: Both 3-component and 1-component "
	    "particle velocity files are defined.  Choose one or the other!\n");
  }

  if (CosmologySimulationParticleVelocityNames[0] != NULL)
    OneComponentPerFile = true;
  else
    OneComponentPerFile = false;

  /* First create a mask in which we will omit particles because of
     static subgrids. */

  int count, region, size = 1;
  int StartRegion[MAX_DIMENSION], EndRegion[MAX_DIMENSION];
  int ActiveDim[MAX_DIMENSION];

  for (dim = 0; dim < GridRank; dim++) {
    ActiveDim[dim] = GridEndIndex[dim] - GridStartIndex[dim] + 1;
    size *= ActiveDim[dim];
  }
  bool *mask = new bool[size];
  bool skip;

  NumberOfParticles = size;
  for (i = 0; i < size; i++)
    mask[i] = true;

  for (region = 0; region < MAX_STATIC_REGIONS; region++) {
    skip = false;
    if (StaticRefineRegionLevel[region] == level) {
      for (dim = 0; dim < MAX_DIMENSION; dim++)
	if (max(StaticRefineRegionLeftEdge[region][dim],  GridLeftEdge[dim]) >=
	    min(StaticRefineRegionRightEdge[region][dim], GridRightEdge[dim]))
	  skip = true;
      if (skip) break;
      for (dim = 0; dim < MAX_DIMENSION; dim++) {
	StartRegion[dim] = (int) nint((StaticRefineRegionLeftEdge[region][dim] - 
				       GridLeftEdge[dim]) / CellWidth[dim][0]);
	EndRegion[dim] = (int) nint((StaticRefineRegionRightEdge[region][dim] - 
				     GridLeftEdge[dim]) / CellWidth[dim][0]);
	StartRegion[dim] = max(StartRegion[dim], 0);
	EndRegion[dim] = min(EndRegion[dim], ActiveDim[dim]);
      }

      for (k = StartRegion[2]; k < EndRegion[2]; k++)
	for (j = StartRegion[1]; j < EndRegion[1]; j++) {
	  index = (k*ActiveDim[1] + j)*ActiveDim[0] + StartRegion[0];
	  for (i = StartRegion[0]; i < EndRegion[0]; i++, index++) {
	    mask[index] = false;
	    NumberOfParticles--;
	  } // ENDFOR i
	} // ENDFOR j
    } // ENDIF level
  } // ENDFOR region

  hid_t file_id;
  herr_t err = 0, h5_error = -1;

  if (OneComponentPerFile)
    file_id = openH5File(CosmologySimulationParticleVelocityNames[0]);
  else
    file_id = openH5File(CosmologySimulationParticleVelocityName);
  assert(file_id != h5_error);

  /* Calculate conversion factor for displacement factor (partly
     stored in the file from mpgrafic)

     * 1/vfact/(ComovingBoxSize/h0) converts velocity in km/s
     to a particle displacement.

  */

  float32 vfact;
  float disp_factor;
  hid_t grp_id = H5Gopen(file_id,"/");
  err |= readAttribute(grp_id, H5T_NATIVE_FLOAT, "vfact", &vfact);
  H5Gclose(grp_id);
  assert(err == 0);
  H5Fclose(file_id);

  disp_factor = 1e-5*VelocityUnits / vfact / (ComovingBoxSize / HubbleConstantNow);

  /* Read in velocities, compute particle positions (with
     displacements), store velocities, and then clean up. */

  float *temp_vel[MAX_DIMENSION];
  float *mass = NULL;
  int  *types = NULL;
  inits_type *tempbuffer = NULL;
  int *int_tempbuffer = NULL;
  for (dim = 0; dim < GridRank; dim++) {
    if (OneComponentPerFile) {
      if (ReadFile(CosmologySimulationParticleVelocityNames[dim], GridRank,
		   GridDimension, GridStartIndex, GridEndIndex, Offset,
		   NULL, &tempbuffer, 0, 1) == FAIL) {
	ENZO_VFAIL("Error reading particle velocity field %"ISYM".\n", dim)
      }
    } else {
      if (ReadFile(CosmologySimulationParticleVelocityName, GridRank,
		   GridDimension, GridStartIndex, GridEndIndex, Offset,
		   NULL, &tempbuffer, dim, 3) == FAIL) {
      ENZO_VFAIL("Error reading particle velocity field %"ISYM".\n", dim)
      }
    } // ENDELSE OneComponentPerFile
    temp_vel[dim] = new float[size];
    for (i = 0; i < size; i++)
      temp_vel[dim][i] = (float) tempbuffer[i];
  }

  if (CosmologySimulationParticleMassName != NULL) {
    if (ReadFile(CosmologySimulationParticleMassName, GridRank,
		 GridDimension, GridStartIndex, GridEndIndex, Offset,
		 NULL, &tempbuffer, 0, 1) == FAIL) {
      ENZO_FAIL("Error reading particle mass.\n");
    }
    mass = new float[size];
    for (i = 0; i < size; i++)
      mass[i] = (float) tempbuffer[i];
  } // ENDIF read masses

  
  if (CosmologySimulationParticleTypeName != NULL) {
    int_tempbuffer = new int[size];
    if (ReadIntFile(CosmologySimulationParticleTypeName, GridRank,
		 GridDimension, GridStartIndex, GridEndIndex, Offset,
		 NULL, &int_tempbuffer, 0, 1) == FAIL) {
      ENZO_FAIL("Error reading particle types.\n");
    }
    
    types = new int[size];
    for (i = 0; i < size; i++)
      types[i] = (int) int_tempbuffer[i];
  } // ENDIF read masses
  

  // Cleanup
  delete [] tempbuffer;
  if (int_tempbuffer != NULL)
    delete [] int_tempbuffer;

  /* Note that this block of code uses the Zel'Dovich approximation to
     extrapolate particle displacements (and thus particle positions) from
     the particle velocity file(s).  This is doable only because it's a linear
     perturbation.  The piece of code right after this will over-write these 
     particle positions if the user has specified a file or files for particle
     displacements.  --BWO, 12 April 2011 */

  // Store positions and velocities
  FLOAT this_pos[3];
  count = 0;
  this->AllocateNewParticles(NumberOfParticles);
  for (k = 0; k < ActiveDim[2]; k++) {
    this_pos[2] = GridLeftEdge[2] + (k+0.5)*CellWidth[2][0];
    for (j = 0; j < ActiveDim[1]; j++) {
      this_pos[1] = GridLeftEdge[1] + (j+0.5)*CellWidth[1][0];
      index = (k*ActiveDim[1] + j)*ActiveDim[0];
      for (i = 0; i < ActiveDim[0]; i++, index++) {
	if (mask[index]) {
	  this_pos[0] = GridLeftEdge[0] + (i+0.5)*CellWidth[0][0];
	  for (dim = 0; dim < MAX_DIMENSION; dim++) {
	    ParticleVelocity[dim][count] = temp_vel[dim][index];
	    ParticlePosition[dim][count] = this_pos[dim] + 
	      disp_factor * ParticleVelocity[dim][count];
	  }
	  count++;
	} // ENDIF mask
      } // ENDFOR i
    } // ENDFOR j
  } // ENDFOR k
	

  /* If the user has specified a particle displacement file (or files), 
     then calculate the particle positions directly from this.  Note that if
     this is done, it writes over the positions calculated previously (in the 
     code directly above this comment).  --BWO, 12 April 2011 */

  if( CosmologySimulationParticleDisplacementName != NULL ||
      CosmologySimulationParticleDisplacementNames[0] != NULL )
    {
      if (CosmologySimulationParticleDisplacementName != NULL &&
	  CosmologySimulationParticleDisplacementNames[0] != NULL) {
	ENZO_FAIL("grid::CosmologyInitializeParticles: Both 3-component and 1-component "
		  "particle displacement files are defined.  Choose one or the other!\n");
      }
		
		
      inits_type *tempbuffer = new inits_type[size];
		
      for (dim = 0; dim < GridRank; dim++) {
	if (CosmologySimulationParticleDisplacementNames[dim] != NULL) {
	  if (ReadFile(CosmologySimulationParticleDisplacementNames[dim], GridRank,
		       GridDimension, GridStartIndex, GridEndIndex, Offset,
		       NULL, &tempbuffer, 0, 1) == FAIL) {
	    ENZO_VFAIL("Error reading particle displacement field %"ISYM".\n", dim)
	      }
	} else {
	  if (ReadFile(CosmologySimulationParticleDisplacementName, GridRank,
		       GridDimension, GridStartIndex, GridEndIndex, Offset,
		       NULL, &tempbuffer, dim, 3) == FAIL) {
	    ENZO_VFAIL("Error reading particle displacement field %"ISYM".\n", dim)
	      }
	} // ENDELSE OneComponentPerFile

	// Write over the particle positions that were calculated previously!
	count = 0;
	for (k = 0; k < ActiveDim[2]; k++) {
	  this_pos[2] = GridLeftEdge[2] + (k+0.5)*CellWidth[2][0];
	  for (j = 0; j < ActiveDim[1]; j++) {
	    this_pos[1] = GridLeftEdge[1] + (j+0.5)*CellWidth[1][0];
	    index = (k*ActiveDim[1] + j)*ActiveDim[0];
	    for (i = 0; i < ActiveDim[0]; i++, index++) {
	      if (mask[index]) {
		this_pos[0] = GridLeftEdge[0] + (i+0.5)*CellWidth[0][0];
		ParticlePosition[dim][count] = this_pos[dim] + tempbuffer[index];
		count++;
	      } // ENDIF mask
	    } // ENDFOR i
	  } // ENDFOR j
	} // ENDFOR k
			
			
      } //ENDFOR dim
		
      delete[] tempbuffer;
    }

  // If provided, store masses
  count = 0;
  if (mass != NULL) {
    for (k = 0; k < ActiveDim[2]; k++)
      for (j = 0; j < ActiveDim[1]; j++) {
	index = (k*ActiveDim[1] + j)*ActiveDim[0];
	for (i = 0; i < ActiveDim[0]; i++, index++)
	  if (mask[index]) {
	    ParticleMass[count] = mass[index];
	    count++;
	  } // ENDIF mask
      } // ENDFOR j
  } // ENDIF mass

  // If provided, store types
  count = 0;
  if (types != NULL) {
    for (k = 0; k < ActiveDim[2]; k++)
      for (j = 0; j < ActiveDim[1]; j++) {
	index = (k*ActiveDim[1] + j)*ActiveDim[0];
	for (i = 0; i < ActiveDim[0]; i++, index++)
	  if (mask[index]) {

	    if (MustRefineParticlesCreateParticles == 2){
	      if (types[index] > 0){
		ParticleType[count] = PARTICLE_TYPE_MUST_REFINE;
	      } else {
		ParticleType[count] = PARTICLE_TYPE_DARK_MATTER;
	      }
	      count++;
	    } else {
	      ParticleType[count] = types[index];
	      count++;
	    }
	  } // ENDIF mask
      } // ENDFOR j
  } // ENDIF types

#ifdef ICPART_SHIFT8
  /* Check to see if the particle is adjacent to static grid boundary
     (level 1 only).  If so, shift it by 1/8th of a cell inwards
     (additional shift accounts for baryons).  This shift fixes an
     error in Enzo's technique to solve the potential with a CIC.  See
     Grid_SolveForPotential.C and GenerateRealization.C (inits) for
     more details. */

  FLOAT shift;
  int dim2, dim3;
  bool inside1[MAX_DIMENSION], inside2[MAX_DIMENSION];

  /* Since we've already converted the 3D particle array to 1D,
     construct a lookup table from (i,j,k) to the particle array
     index */

  count = 0;
  int *index_arr = new int[size];
  for (k = 0; k < ActiveDim[2]; k++)
    for (j = 0; j < ActiveDim[1]; j++) {
      index = (k*ActiveDim[1] + j)*ActiveDim[0];
      for (i = 0; i < ActiveDim[0]; i++, index++)
	if (mask[index]) {
	  index_arr[index] = count;
	  count++;
	} else
	  index_arr[index] = INT_UNDEFINED;
    } // ENDFOR j

  shift = CellWidth[0][0] / 
    (8.0 * (1.0 - CosmologySimulationOmegaBaryonNow/OmegaMatterNow));

  if (level == 0) {
    for (region = 0; region < MAX_STATIC_REGIONS; region++) {
      if (StaticRefineRegionLevel[region] == level) {
	for (dim = 0; dim < MAX_DIMENSION; dim++) {

	  StartRegion[dim] = (int) nint((StaticRefineRegionLeftEdge[region][dim] - 
					 GridLeftEdge[dim]) / CellWidth[dim][0]);
	  EndRegion[dim] = (int) nint((StaticRefineRegionRightEdge[region][dim] - 
				       GridLeftEdge[dim]) / CellWidth[dim][0]);

	  // Check if the faces (dim=StartRegion-1 and dim=EndRegion)
	  // are on this processor.
	  inside1[dim] = StartRegion[dim]-1 >= 0 &&
	    StartRegion[dim]-1 < ActiveDim[dim];
	  inside2[dim] = EndRegion[dim]+1 > 0 && EndRegion[dim]+1 < ActiveDim[dim];
	  StartRegion[dim] = max(StartRegion[dim], 0);
	  EndRegion[dim] = min(EndRegion[dim], ActiveDim[dim]);
	} // ENDFOR dim

	for (dim = 0; dim < MAX_DIMENSION; dim++) {

	  dim2 = (dim+1) % MAX_DIMENSION;
	  dim3 = (dim+2) % MAX_DIMENSION;

	  for (j = StartRegion[dim3]; j < EndRegion[dim3]; j++) {
	    for (i = StartRegion[dim2]; i < EndRegion[dim2]; i++) {
	      switch (dim) {
	      case 0:
		index1 = (j*ActiveDim[1] + i)*ActiveDim[0] + StartRegion[0] - 1;
		index2 = (j*ActiveDim[1] + i)*ActiveDim[0] + EndRegion[0] + 1;
		break;
	      case 1:
		index1 = (i*ActiveDim[1] + StartRegion[1] - 1)*ActiveDim[0] + j;
		index2 = (i*ActiveDim[1] + EndRegion[1] + 1)*ActiveDim[0] + j;
		break;
	      case 2:
		index1 = ((StartRegion[2]-1)*ActiveDim[1] + j)*ActiveDim[0] + i;
		index2 = ((EndRegion[2]+1)*ActiveDim[1] + j)*ActiveDim[0] + i;
		break;
	      } // ENDSWITCH dim

	      if (inside1[dim])
		ParticlePosition[dim][index_arr[index1]] += shift;
	      if (inside2[dim])
		ParticlePosition[dim][index_arr[index2]] -= shift;
	      
	    } // ENDFOR i
	  } // ENDFOR j
	} // ENDFOR dims
      } // ENDIF right level
    } // ENDFOR region
  } // ENDIF level 0
  delete [] index_arr;
#endif /* ICPART_SHIFT8 */

  // Clean up
  for (dim = 0; dim < GridRank; dim++)
    delete [] temp_vel[dim];
  delete [] mask;
  if (mass != NULL)
    delete [] mass;
  if (types != NULL)

    delete [] types;

  return SUCCESS;

}
