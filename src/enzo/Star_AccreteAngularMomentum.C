/***********************************************************************
/
/  ADD ACCRETED ANGULAR MOMENTUM ONTO THE STAR PARTICLE (for MBH_JETS)
/
/  written by: Ji-hoon Kim
/  date:       November, 2009
/  modified1: 
/
/  purpose:  For MBHFeedback = 2 to 5 (FeedbackFlag=MBH_JETS) calculate 
/            angular momentum accreted onto MBH; this will not affect 
/            Star_Accrete.  
/
************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
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

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int Star::AccreteAngularMomentum(void)
{

  if (CurrentGrid == NULL || 
      type != MBH || 
      (MBHFeedback < 2 || MBHFeedback > 5))
    return SUCCESS;

  if (CurrentGrid->GridRank != MAX_DIMENSION) {
    ENZO_FAIL("star:AccreteAngularMomentum: 1 or 2 dimension is not implemented! \n");
  }

  int dim, i, j, k, index, ibuff = NumberOfGhostZones;
  double gas_angmom[] = {0.0, 0.0, 0.0}, total_gas_mass = 0.0, gas_mass = 0.0;
  FLOAT CellVolume = 1, BoxSize = 1, DensityConversion = 1, VelocityConversion = 1;
  FLOAT a = 1, dadt;
  FLOAT delx, dely, delz, velx, vely, velz, time = CurrentGrid->Time;
  const double Msun = 1.989e33, Mpc = 3.086e24;

  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, time);

  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  CurrentGrid->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
					  Vel3Num, TENum);

  if (ComovingCoordinates) {
    CosmologyComputeExpansionFactor(time, &a, &dadt);
    BoxSize = ComovingBoxSize/HubbleConstantNow*a/(1+InitialRedshift);
  } else
    BoxSize = LengthUnits/Mpc; // to Mpc

  DensityConversion = FLOAT(double(DensityUnits) / Msun * pow(Mpc, 3)); // to Msun/Mpc^3
  VelocityConversion = FLOAT(double(VelocityUnits) / 1.0e5); // to km/s

  float CellWidthTemp = float(CurrentGrid->CellWidth[0][0]);

  for (dim = 0; dim < MAX_DIMENSION; dim++) 
    CellVolume *= CellWidthTemp*BoxSize; // in Mpc^3

  /* indices for Star particle */

  i = (int)((pos[0] - CurrentGrid->CellLeftEdge[0][0]) / CellWidthTemp);
  j = (int)((pos[1] - CurrentGrid->CellLeftEdge[1][0]) / CellWidthTemp);
  k = (int)((pos[2] - CurrentGrid->CellLeftEdge[2][0]) / CellWidthTemp);

  /* Check whether the current grid contains the whole 27 cells */

//   if (i < ibuff+1 || i > CurrentGrid->GridDimension[0]-ibuff-2 || 
//       j < ibuff+1 || j > CurrentGrid->GridDimension[1]-ibuff-2 ||
//       k < ibuff+1 || k > CurrentGrid->GridDimension[2]-ibuff-2) {
//     fprintf(stdout, "star::AAM: 27 cells around MBH not contained; moving on.\n"); 
//     return SUCCESS;
//   }
  
  /* Find angular momentum in 27 cells */

  for (int kk = -1; kk <= 1; kk++) {
    // relative position
    delz = (CurrentGrid->CellLeftEdge[2][k+kk] + 0.5*CurrentGrid->CellWidth[2][0] 
	    - pos[2]) * BoxSize; // in Mpc

    for (int jj = -1; jj <= 1; jj++) {
      dely = (CurrentGrid->CellLeftEdge[1][j+jj] + 0.5*CurrentGrid->CellWidth[1][0] 
	      - pos[1]) * BoxSize;

      for (int ii = -1; ii <= 1; ii++) {
	delx = (CurrentGrid->CellLeftEdge[0][i+ii] + 0.5*CurrentGrid->CellWidth[0][0] 
		- pos[0]) * BoxSize;

	index = i+ii+(j+jj+(k+kk)*CurrentGrid->GridDimension[1])*CurrentGrid->GridDimension[0];
	gas_mass = CurrentGrid->BaryonField[DensNum][index] * 
	  DensityConversion * CellVolume; // in Msun

	// relative velocity
	velx = (CurrentGrid->BaryonField[Vel1Num][index] - vel[0]) * VelocityConversion; // in km/s
	vely = (CurrentGrid->BaryonField[Vel2Num][index] - vel[1]) * VelocityConversion;
	velz = (CurrentGrid->BaryonField[Vel3Num][index] - vel[2]) * VelocityConversion;
	
	// store gas angular momentum in: Msun * Mpc * km/s
	gas_angmom[0] += gas_mass * ( vely*delz - velz*dely); 
	gas_angmom[1] += gas_mass * (-velx*delz + velz*delx);
	gas_angmom[2] += gas_mass * ( velx*dely - vely*delx);
	total_gas_mass += gas_mass;
      }
    }
  }

//  fprintf(stdout, "star::AAM: gas_angmom = (%g, %g, %g), total_gas_mass = %g, DeltaMass = %g\n", 
//	  gas_angmom[0], gas_angmom[1], gas_angmom[2], total_gas_mass, DeltaMass); 

  // specific gas angular momentum in: Mpc * km/s
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    gas_angmom[dim] /= total_gas_mass;

  // finally store angular momentum onto MBH in: Msun * Mpc * km/s
  for (dim = 0; dim < MAX_DIMENSION; dim++) 
    this->accreted_angmom[dim] += (float)(this->DeltaMass * gas_angmom[dim]);

//  fprintf(stdout, "star::AAM: this->accreted_angmom = (%g, %g, %g)\n", 
//	  this->accreted_angmom[0], this->accreted_angmom[1], this->accreted_angmom[2]); 

  return SUCCESS;
}
