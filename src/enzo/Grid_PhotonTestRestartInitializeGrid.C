/***********************************************************************
/
/  GRID CLASS (RE-INITIALIZE FOR PUTTING IN SINK PARTICLES)
/
/  written by: Elizabeth Harper-Clark
/  date:       March, 2010
/  modified1:
/
/  PURPOSE: based on SupernovaeRestartInitialize
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/
 
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
 int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int grid::PhotonTestRestartInitialize(int level, int *NumberOfCellsSet)
{
  /* declarations */
 
  int dim, i, j, k,l, n;
  FLOAT delx, dely, delz, radius2, DomainWidth[MAX_DIMENSION];
 int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
    DINum, DIINum, HDINum,  kphHINum, gammaNum, kphHeINum,
    kphHeIINum, kdissH2INum, RPresNum1, RPresNum2, RPresNum3; 

  float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1.0, TimeUnits = 1.0,
    VelocityUnits = 1.0;
  if (UsePhysicalUnit)
    GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, Time);
  double MassUnits = DensityUnits*pow(LengthUnits,3);
  //printf("Mass Units = %g \n",MassUnits);
  //printf("Time Units = %g \n",TimeUnits);

  for (dim = 0; dim < GridRank; dim++)
    DomainWidth[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];
 
  /* Find fields: density, total energy, velocity1-3. */
 
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }


  if (RadiativeTransfer && (MultiSpecies < 1)) {
    ENZO_FAIL("Grid_PhotonTestInitialize: Radiative Transfer but not MultiSpecies set");
  }

  //   Allocate fields for photo ionization and heating rates
//   if (RadiativeTransfer)
//     if (MultiSpecies) {
//       FieldType[kphHINum    = NumberOfBaryonFields++] = kphHI;
//       FieldType[gammaNum    = NumberOfBaryonFields++] = PhotoGamma;
//       if (RadiativeTransferHydrogenOnly == FALSE) {
// 	FieldType[kphHeINum   = NumberOfBaryonFields++] = kphHeI;
// 	FieldType[kphHeIINum  = NumberOfBaryonFields++] = kphHeII;
//       }
//       if (MultiSpecies > 1) 
// 	FieldType[kdissH2INum    = NumberOfBaryonFields++] = kdissH2I;
//     } 

//   if (RadiationPressure && RadiativeTransfer) {
//     FieldType[RPresNum1 = NumberOfBaryonFields++] = RadPressure0;
//     FieldType[RPresNum2 = NumberOfBaryonFields++] = RadPressure1;
//     FieldType[RPresNum3 = NumberOfBaryonFields++] = RadPressure2;
//   }

  NumberOfPhotonPackages = 0;
  PhotonPackages-> NextPackage= NULL;


  const double Mpc = 3.0856e24, SolarMass = 1.989e33, GravConst = 6.67e-8,
               pi = 3.14159, mh = 1.67e-24, kboltz = 1.381e-16;

 
  /* Initialize radiation fields - not needed in restart?? */

//   if (this->InitializeRadiativeTransferFields() == FAIL) {
//     ENZO_FAIL("\nError in InitializeRadiativeTransferFields.\n");
//   }



  return SUCCESS;
}
