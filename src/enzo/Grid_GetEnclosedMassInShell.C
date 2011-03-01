/***********************************************************************
/
/  GET ENCLOSED MASS IN A SHELL
/
/  written by: John Wise
/  date:       September, 2006
/  modified1:  March, 2010 by JHW -- modified from sphere to shell
/
/  PURPOSE:
/
************************************************************************/
#include <stdlib.h>
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
#include "LevelHierarchy.h"
#include "CosmologyParameters.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int grid::GetEnclosedMassInShell(Star *star, float radius0, float radius1, 
				 float &mass, float &metallicity2, 
				 float &metallicity3,
				 float &coldgas_mass, float AvgVelocity[])
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
 
  /* Quick check to see if shell overlaps this grid.  Grid must
     overlap with a sphere of radius=radius1 but not be completely
     contained in a sphere with radius=radius0 */

  int i, j, k, index, dim;
  bool inside_outer = true, contained_inner = true;

  for (dim = 0; dim < GridRank; dim++) {
    contained_inner &= (GridLeftEdge[dim]  >= star->pos[dim] - radius0 &&
			GridRightEdge[dim] <= star->pos[dim] + radius0);
    inside_outer &= !(star->pos[dim] - radius1 > GridRightEdge[dim] ||
		      star->pos[dim] + radius1 < GridLeftEdge[dim]);
  }

  if (!inside_outer || contained_inner)
    return SUCCESS;

  this->DebugCheck((char*) "Grid_GetEnclosedMass");

  FLOAT DomainWidth[MAX_DIMENSION];
  int size = 1;
  for (dim = 0; dim < GridRank; dim++) {
    DomainWidth[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];
    size *= GridDimension[dim];
  }

  /* Set the units. */

  const double Msun = 1.989e33;
  float DensityUnits, LengthUnits, TemperatureUnits, 
    TimeUnits, VelocityUnits;

  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, 
	   &VelocityUnits, Time);

  float CellWidthTemp = float(CellWidth[0][0]);

  float MassConversion = (float) (double(LengthUnits*CellWidthTemp) *
				  double(LengthUnits*CellWidthTemp) *
				  double(LengthUnits*CellWidthTemp) *
				  double(DensityUnits) / Msun);

  int DensNum = FindField(Density, FieldType, NumberOfBaryonFields);
  int Vel1Num = FindField(Velocity1, FieldType, NumberOfBaryonFields);
  int GENum = FindField(InternalEnergy, FieldType, NumberOfBaryonFields);
  int ColorField = FindField(ForbiddenRefinement, FieldType, NumberOfBaryonFields); 
  if ((ColorField < 0) && (star->ReturnFeedbackFlag() == COLOR_FIELD))
      ENZO_FAIL("Couldn't Find Color Field!");

  /* Find metallicity field and set flag. */

  int MetallicityField = FALSE, MetalNum;
  if ((MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields)) 
      != -1)
    MetallicityField = TRUE;
  else
    MetalNum = 0;

  /* Find SN colour field */

  int UseColour = FALSE, SNColourNum;
  if ((SNColourNum = FindField(SNColour, FieldType, NumberOfBaryonFields)) 
      != -1)
    UseColour = TRUE;
  else
    SNColourNum = 0;

  /* Calculate temperature for cold gas accretion (use
     ComputeTemperatureField for non-DualEnergyFormalism).  Divide by
     TemperatureUnits for DualEnergyFormalism.
     ComputeTemperatureField is slow if used repeatedly as in the loop
     in finding a sphere with an enclosed gas mass. */

  float *ThresholdField = NULL;
  float ColdThreshold;
  float MinimumTemperature = (MultiSpecies > 1) ? 1e3 : 1e4;
  if (DualEnergyFormalism) {
    ColdThreshold = MinimumTemperature / (TemperatureUnits * (Gamma-1.0) * 0.6);
    ThresholdField = BaryonField[GENum];
  } else {
    ColdThreshold = MinimumTemperature;
    ThresholdField = new float[size];  // i.e. temperature
    this->ComputeTemperatureField(ThresholdField);
  }

  /* Set under subgrid field to zero */


  FLOAT delx, dely, delz;
  float gasmass, dr2, radius0_2, radius1_2;

  radius0_2 = radius0 * radius0;
  radius1_2 = radius1 * radius1;

  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    delz = CellLeftEdge[2][k] + 0.5*CellWidth[2][k] - star->pos[2];
    delz = min(delz, DomainWidth[2]-delz);
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      dely = CellLeftEdge[1][j] + 0.5*CellWidth[1][j] - star->pos[1];
      dely = min(dely, DomainWidth[1]-dely);
      index = (k*GridDimension[1] + j)*GridDimension[0] + GridStartIndex[0];
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) { 

	if (BaryonField[NumberOfBaryonFields][index] != 0.0)
	  continue;

	delx = CellLeftEdge[0][i] + 0.5*CellWidth[0][i] - star->pos[0];
	delx = min(delx, DomainWidth[0]-delx);

	dr2 = delx*delx + dely*dely + delz*delz;

	if (dr2 >= radius0_2 && dr2 < radius1_2) {
	  gasmass = BaryonField[DensNum][index] * MassConversion;
	  mass += gasmass;
	  if (star->ReturnFeedbackFlag() == COLOR_FIELD)
	    BaryonField[ColorField][index] = BaryonField[DensNum][index];
	  for (dim = 0; dim < GridRank; dim++)
	    AvgVelocity[dim] += BaryonField[Vel1Num+dim][index] * gasmass;
	  if (ThresholdField[index] < ColdThreshold)
	    coldgas_mass += gasmass;
	  if (MetallicityField)
	    metallicity2 += BaryonField[MetalNum][index] * MassConversion;
	  if (UseColour)
	    metallicity3 += BaryonField[SNColourNum][index] * MassConversion;
	}
	
      }
    }
  }

  if (!DualEnergyFormalism)
    delete [] ThresholdField;

  return SUCCESS;

}
