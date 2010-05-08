/***********************************************************************
/
/  GET ENCLOSED MASS WITHIN SOME RADIUS
/
/  written by: John Wise
/  date:       September, 2006
/  modified1:
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

int grid::GetEnclosedMass(Star *star, float radius, float &mass,
			  float &metallicity, float &coldgas_mass, 
			  float AvgVelocity[])
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  /* Quick check to see if sphere overlaps this grid. */

  int i, j, k, index, dim;

  for (dim = 0; dim < GridRank; dim++)
    if (star->pos[dim] - radius > GridRightEdge[dim] ||
	star->pos[dim] + radius < GridLeftEdge[dim]   )
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

  // Cold gas temperature threshold = 10^4 for no H2 cooling and 10^3
  // with H2 cooling
  float ColdTemperature = (MultiSpecies > 1) ? 1e3 : 1e4;
 
  float *ThresholdField = NULL;
  float ColdThreshold;
  if (DualEnergyFormalism) {
    ColdThreshold = ColdTemperature / (TemperatureUnits * (Gamma-1.0) * 0.6);
    ThresholdField = BaryonField[GENum];
  } else {
    ColdThreshold = ColdTemperature;
    ThresholdField = new float[size];  // i.e. temperature
    this->ComputeTemperatureField(ThresholdField);
  }

  /* Set under subgrid field to zero */


  FLOAT delx, dely, delz;
  float gasmass, dr2, radius2 = radius*radius;

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

	if (dr2 < radius2) {
	  gasmass = BaryonField[DensNum][index] * MassConversion;
	  mass += gasmass;
	  if (star->ReturnFeedbackFlag() == COLOR_FIELD)
	    BaryonField[ColorField][index] = BaryonField[DensNum][index];
	  for (dim = 0; dim < GridRank; dim++)
	    AvgVelocity[dim] += BaryonField[Vel1Num+dim][index] * gasmass;
	  if (ThresholdField[index] < ColdThreshold)
	    coldgas_mass += gasmass;
	  if (MetallicityField)
	    metallicity += BaryonField[MetalNum][index] * MassConversion;
	  if (UseColour)
	    metallicity += BaryonField[SNColourNum][index] * MassConversion;
	}
	
      }
    }
  }

  if (!DualEnergyFormalism)
    delete [] ThresholdField;

  return SUCCESS;

}







// GetEnclosedMass when not using Star object - Ji-hoon Kim, Jan, 2010

int grid::GetEnclosedMass(FLOAT star_pos[], float radius, float &mass,
			  float &metallicity, float &coldgas_mass, 
			  float AvgVelocity[], float &OneOverRSquaredSum)
{

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
 
  /* Quick check to see if sphere overlaps this grid. */

  int i, j, k, index, dim;

  for (dim = 0; dim < GridRank; dim++)
    if (star_pos[dim] - radius > GridRightEdge[dim] ||
	star_pos[dim] + radius < GridLeftEdge[dim]   )
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

  // Cold gas temperature threshold = 10^4 for no H2 cooling and 10^3
  // with H2 cooling
  float ColdTemperature = (MultiSpecies > 1) ? 1e3 : 1e4;
 
  float *ThresholdField = NULL;
  float ColdThreshold;
  if (DualEnergyFormalism) {
    ColdThreshold = ColdTemperature / (TemperatureUnits * (Gamma-1.0) * 0.6);
    ThresholdField = BaryonField[GENum];
  } else {
    ColdThreshold = ColdTemperature;
    ThresholdField = new float[size];  // i.e. temperature
    this->ComputeTemperatureField(ThresholdField);
  }

  /* Set under subgrid field to zero */


  FLOAT delx, dely, delz;
  float gasmass, dr2, radius2 = radius*radius;

  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    delz = CellLeftEdge[2][k] + 0.5*CellWidth[2][k] - star_pos[2];
    delz = min(delz, DomainWidth[2]-delz);
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      dely = CellLeftEdge[1][j] + 0.5*CellWidth[1][j] - star_pos[1];
      dely = min(dely, DomainWidth[1]-dely);
      index = (k*GridDimension[1] + j)*GridDimension[0] + GridStartIndex[0];
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) { 

	if (BaryonField[NumberOfBaryonFields][index] != 0.0)
	  continue;

	delx = CellLeftEdge[0][i] + 0.5*CellWidth[0][i] - star_pos[0];
	delx = min(delx, DomainWidth[0]-delx);

	dr2 = delx*delx + dely*dely + delz*delz;

	if (dr2 < radius2) {
	  gasmass = BaryonField[DensNum][index] * MassConversion;
	  mass += gasmass;
	  for (dim = 0; dim < GridRank; dim++)
	    AvgVelocity[dim] += BaryonField[Vel1Num+dim][index] * gasmass;
	  if (ThresholdField[index] < ColdThreshold)
	    coldgas_mass += gasmass;
	  if (MetallicityField)
	    metallicity += BaryonField[MetalNum][index] * MassConversion;
	  if (UseColour)
	    metallicity += BaryonField[SNColourNum][index] * MassConversion;
	  //OneOverRSqauredSum used in Grid_AddFeedbackSphere for MBHFeedback=1
	  //imposed upperbound of 1/(2*CellWidth)^2
	  OneOverRSquaredSum += min(1.0/dr2, 1.0/(4.0*CellWidthTemp*CellWidthTemp));  
	}
	
      }
    }
  }

  if (!DualEnergyFormalism)
    delete [] ThresholdField;

  return SUCCESS;

}
