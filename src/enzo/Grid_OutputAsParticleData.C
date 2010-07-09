/***********************************************************************
/
/  GRID CLASS (CONVERT GRID DATA TO PARTICLE DATA FOR OUTPUT)
/
/  written by: Greg Bryan
/  date:       August, 1996
/  modified1:
/
/  PURPOSE:
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "CosmologyParameters.h"
#include "Grid.h"
 
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
 
 
int grid::OutputAsParticleData(FLOAT RegionLeftEdge[], FLOAT RegionRightEdge[],
			  ListOfParticles *ParticleList[NUM_PARTICLE_TYPES],
			       float BaseRadius)
{
 
  if (BaryonField[NumberOfBaryonFields] == NULL) {
    ENZO_FAIL("UNDER_SUBGRID_FLAG field not set.\n");
  }
 
  if (SelfGravity && GravityResolution != 1) {
    ENZO_FAIL("OutputAsParticleData assumes GravityResolution == 1.\n");
  }
 
  /* Declarations */
 
  int i,j,k, n, dim, field, findex, start[] = {0,0,0}, stop[] = {0,0,0};
  const double SolarMass = 1.989e33, Mpc = 3.086e24;
 
  /* Set the Conversion factor for velocity. */
 
  float DensityConversion = 1, VelocityConversion = 1;
  FLOAT a = 1, dadt;
  float TemperatureUnits=1, DensityUnits=1, LengthUnits=1, VelocityUnits=1, 
    TimeUnits=1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }

  if (ComovingCoordinates) {
 
    CosmologyComputeExpansionFactor(Time, &a, &dadt);
 
    /* Convert cgs units to more reasonable values.
       density:   M(solar)/Mpc^3
       velocity:  km/s */
 
    VelocityConversion = VelocityUnits / 1.0e5;
    DensityConversion = float(double(DensityUnits) / SolarMass * POW(Mpc, 3));
  }
 
  /* Set number of values for particle lists. */
 
  ParticleList[0]->NumberOfValues = 3;  /* gas */
  ParticleList[1]->NumberOfValues = 2;  /* dm */
  ParticleList[2]->NumberOfValues = 2;  /* stars */
 
  for (i=0; i < NUM_PARTICLE_TYPES; i++)
    ParticleList[i]->NumberOfParticles = 0;
 
  /* Check To see if grid overlaps the projected field. */
 
  float BoxSize = 1, CellVolume = 1;
  if (ComovingCoordinates)
    BoxSize = ComovingBoxSize/HubbleConstantNow*a/(1+InitialRedshift);
  int size = 1;
  for (dim = 0; dim < GridRank; dim++) {
    if (GridLeftEdge[dim] > RegionRightEdge[dim] ||
	GridRightEdge[dim] < RegionLeftEdge[dim])
      return SUCCESS;
    CellVolume *= CellWidth[dim][0]*BoxSize;
    size *= GridDimension[dim];
  }
 
  /*-------------------------------------------------------------------*/
  /* Set Baryon particle list */
 
  if (NumberOfBaryonFields > 0) {
 
    /* Find fields: density, total energy, velocity1-3. */
 
    int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
    if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
					 Vel3Num, TENum) == FAIL) {
      ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
    }
 
    /* Compute the temperature. */
 
    float *temperature = new float[size];
    if (this->ComputeTemperatureField(temperature) == FAIL) {
      ENZO_FAIL("Error in grid->ComputeTemperatureField.\n");
    }
 
    /* Find the start and stop indicies in the ProjectionDimension of this
       grid for the projected region. */
 
    for (dim = 0; dim < GridRank; dim++) {
      start[dim] = max(int((RegionLeftEdge[dim] - GridLeftEdge[dim]) /
			   CellWidth[dim][0]), 0) + GridStartIndex[dim];
      stop[dim]  = min(int((RegionRightEdge[dim] - GridLeftEdge[dim]) /
			   CellWidth[dim][0]),
		       GridEndIndex[dim] - GridStartIndex[dim]) +
	GridStartIndex[dim];
    }
 
    /* Count number of particles. */
 
    size = 0;
    for (k = start[2]; k <= stop[2]; k++)
      for (j = start[1]; j <= stop[1]; j++) {
	findex = (j + k*GridDimension[1])*GridDimension[0] + start[0];
	for (i = start[0]; i <= stop[0]; i++, findex++)
	  if (BaryonField[NumberOfBaryonFields][findex] == 0)
	    size++;
      }
 
    /* Allocate space for particles. */
 
    for (dim = 0; dim < GridRank; dim++) {
      ParticleList[0]->ParticlePosition[dim] = new float[size];
      ParticleList[0]->ParticleVelocity[dim] = new float[size];
    }
    ParticleList[0]->ParticleRadius = new float[size];
    for (field = 0; field < ParticleList[0]->NumberOfValues; field++)
      ParticleList[0]->ParticleValue[field] = new float[size];
 
    /* Create particles. */
 
    int index = 0;
    for (k = start[2]; k <= stop[2]; k++)
      for (j = start[1]; j <= stop[1]; j++) {
	findex = (j + k*GridDimension[1])*GridDimension[0] + start[0];
	for (i = start[0]; i <= stop[0]; i++, findex++)
	  if (BaryonField[NumberOfBaryonFields][findex] == 0) {
	    ParticleList[0]->ParticlePosition[0][index] =
	      CellLeftEdge[0][i] + 0.5*CellWidth[0][0];
	    if (GridRank > 1)
	      ParticleList[0]->ParticlePosition[1][index] =
		CellLeftEdge[1][j] + 0.5*CellWidth[1][0];
	    if (GridRank > 2)
	      ParticleList[0]->ParticlePosition[2][index] =
		CellLeftEdge[2][k] + 0.5*CellWidth[2][0];
	    //	    ParticleList[0]->ParticleRadius[index] = CellWidth[0][0];
 
	    for (dim = 0; dim < GridRank; dim++)
	      ParticleList[0]->ParticleVelocity[dim][index] =
		BaryonField[Vel1Num+dim][findex]*VelocityConversion;
 
	    ParticleList[0]->ParticleRadius[index] = BaseRadius *
	      POW(1.0/BaryonField[DensNum][findex], 0.3333);
 
	    ParticleList[0]->ParticleValue[0][index] =
	            BaryonField[DensNum][findex]*CellVolume*DensityConversion;
	    ParticleList[0]->ParticleValue[1][index] =
	      BaryonField[DensNum][findex];
	    ParticleList[0]->ParticleValue[2][index] = temperature[findex];
 
/*	    for (field = 0; field < NumberOfBaryonFields; field++)
	      ParticleList[0]->ParticleValue[field+2][index] =
		BaryonField[field][findex];
*/
	    index++;
	  }
      }
 
    /* Set number of Particles. */
 
    ParticleList[0]->NumberOfParticles = index;
 
    /* Clean up. */
 
    delete temperature;
 
  }
 
  /*-------------------------------------------------------------------*/
  /* Set dark matter particle list */
 
  /* Count number of particles in region. */
 
  size = 0;
  for (n = 0; n < NumberOfParticles; n++) {
 
    /* Check to see if particle is inside region. */
 
    int ParticleInVolume = TRUE;
    for (dim = 0; dim < GridRank; dim++)
      if (ParticlePosition[dim][n] > RegionRightEdge[dim] ||
	  ParticlePosition[dim][n] < RegionLeftEdge[dim])
	ParticleInVolume = FALSE;
 
    if (ParticleInVolume)
      ParticleList[ParticleType[n]]->NumberOfParticles++;
 
  }
 
  /* Allocate space for particles. */
 
  int ntype, ParticleCount[NUM_PARTICLE_TYPES];
  for (ntype = 1; ntype < NUM_PARTICLE_TYPES; ntype++) {
 
    size = ParticleList[ntype]->NumberOfParticles;
    ParticleCount[ntype] = 0;
 
    if (size > 0) {
      for (dim = 0; dim < GridRank; dim++) {
	ParticleList[ntype]->ParticlePosition[dim] = new float[size];
	ParticleList[ntype]->ParticleVelocity[dim] = new float[size];
      }
      ParticleList[ntype]->ParticleRadius    = new float[size];
      ParticleList[ntype]->ParticleValue[0]  = new float[size];
      ParticleList[ntype]->ParticleValue[1]  = new float[size];
      ParticleList[ntype]->ParticleIndex     = new PINT[size];
      ParticleList[ntype]->NumberOfParticles = size;
    } // end: if (size > 0)
 
  } // end: loop over particle types
 
  /* Create Particles. */
 
  int FieldPosition[MAX_DIMENSION], index, count, itype;
  float density;
  for (n = 0; n < NumberOfParticles; n++) {
 
    /* Check to see if particle is inside region. */
 
    int ParticleInVolume = TRUE;
    for (dim = 0; dim < GridRank; dim++)
      if (ParticlePosition[dim][n] > RegionRightEdge[dim] ||
	  ParticlePosition[dim][n] < RegionLeftEdge[dim])
	ParticleInVolume = FALSE;
 
    if (ParticleInVolume) {
 
      /* Look up particle density. */
 
      for (dim = 0; dim < GridRank; dim++)
	FieldPosition[dim] =
	  int((ParticlePosition[dim][n] -
	       GravitatingMassFieldParticlesLeftEdge[dim])
	      /GravitatingMassFieldParticlesCellSize);
 
      index = (FieldPosition[1] +
	       FieldPosition[2]*GravitatingMassFieldParticlesDimension[1])*
	GravitatingMassFieldParticlesDimension[0] + FieldPosition[0];
      if (GravitatingMassFieldParticles != NULL)
	density = GravitatingMassFieldParticles[index];
      else
	density = 0;
      if (density == 0) density = ParticleMass[n];
 
      /* set particle type (for convenience). */
 
      itype = ParticleType[n];
 
      /* Set position, radius and mass. */
 
      count = ParticleCount[itype];
      for (dim = 0; dim < GridRank; dim++) {
	ParticleList[itype]->ParticlePosition[dim][count] =
	  float(ParticlePosition[dim][n]);
	ParticleList[itype]->ParticleVelocity[dim][count] =
	  ParticleVelocity[dim][n] * VelocityConversion;
      }
      //      ParticleList[itype]->ParticleRadius[count] = CellWidth[0][0]*
      //	max(POW(ParticleMass[n]/density, 0.3333), 3);
      ParticleList[itype]->ParticleRadius[count] = BaseRadius *
	POW(1.0/density, 0.3333);
      ParticleList[itype]->ParticleValue[0][count] =
	ParticleMass[n]*CellVolume*DensityConversion;
      ParticleList[itype]->ParticleValue[1][count] = density;
      ParticleList[itype]->ParticleIndex[count] = ParticleNumber[n];
 
      ParticleCount[itype]++;
 
    } // end: if (ParticleInVolume)
 
  } // end: loop over particles
 
  if (debug)

    printf("Grid %"ISYM" %"ISYM" %"ISYM": NumberOfParticles = %"ISYM" %"ISYM"\n", GridDimension[0],
	   GridDimension[1], GridDimension[2],
	   ParticleList[0]->NumberOfParticles,
	   ParticleList[1]->NumberOfParticles);
 
  return SUCCESS;
}
