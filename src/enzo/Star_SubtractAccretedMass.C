
/***********************************************************************
/
/  SUBTRACT ACCRETED MASS FROM CELLS
/
/  written by: Ji-hoon Kim
/  date:       September, 2009
/  modified1: 
/
/  PURPOSE: This routine subtracts the accreted gas mass out 
/           for Star particle type BlackHole and MBH.
/           Note that accretion_rate is calculted in Star_CalculateMassAccretion.C
/           but DeltaMass is calculated in Star_Accrete.C.
/           At the moment, the mass is subtracted 
/           only from the cell the Star particle resides in.
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
#include "CosmologyParameters.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "StarParticleData.h"

int FindField(int field, int farray[], int numfields);

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int Star::SubtractAccretedMass(void)
{

  /* Check if the star type is correct */

  if ((this->type != BlackHole && abs(this->type) != MBH) || (this->CurrentGrid == NULL))
    return SUCCESS;

  int dim, igrid[MAX_DIMENSION], index, size;
  float OldDensity, NewDensity, factor;
  float densgrid, ugrid, vgrid, wgrid, denssink, usink, vsink, wsink, drho;
  double Msun = 1.989e33;
  FLOAT time = CurrentGrid->OldTime;

  if (time <= 0)
    time = CurrentGrid->Time - CurrentGrid->dtFixed;

  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, time);

  /* Find metallicity or SN Color field and set flag. */

  /*
  int ZNum, ZField;
  int MetallicityField = FALSE, MetalNum;
  if ((MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields)) 
      != -1)
    MetallicityField = TRUE;
  else
    MetalNum = 0;

  int UseColour = FALSE, SNColourNum;
  if ((SNColourNum = FindField(SNColour, FieldType, NumberOfBaryonFields)) 
      != -1)
    UseColour = TRUE;
  else
    SNColourNum = 0;

  ZNum = max(MetalNum, SNColourNum);
  ZField = max(MetallicityField, UseColour);
  */

  /* Find fields: density, total energy, velocity1-3. */

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (CurrentGrid->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				       Vel3Num, TENum) == FAIL) {
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
    ENZO_FAIL("");
  }
  
  /* Find Multi-species fields. */

  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
  if (MultiSpecies) 
    if (CurrentGrid->IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, 
				    HeIIINum, HMNum, H2INum, H2IINum, DINum, 
				    DIINum, HDINum) == FAIL) {
      fprintf(stderr, "Error in grid->IdentifySpeciesFields.\n");
      ENZO_FAIL("");
    }

  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    size *= CurrentGrid->GridDimension[dim];
    igrid[dim] = (int) (pos[dim] - CurrentGrid->GridLeftEdge[dim]) /
      CurrentGrid->CellWidth[0][0];
  }

  index = 
    ((igrid[2] + CurrentGrid->GridStartIndex[2]) * CurrentGrid->GridDimension[1] + 
     igrid[1] + CurrentGrid->GridStartIndex[1]) * CurrentGrid->GridDimension[0] + 
     igrid[0] + CurrentGrid->GridStartIndex[0];

  float MassConversion = (float) (pow(LengthUnits * CurrentGrid->CellWidth[0][0], 3.0)
				  * double(DensityUnits) / Msun);

  /* Subtract accreted mass from the grids, and calculate new densities */
  
  OldDensity = CurrentGrid->BaryonField[DensNum][index];
  NewDensity = OldDensity - this->DeltaMass * Msun / 
    pow(CurrentGrid->CellWidth[0][0]*LengthUnits, 3.0) / DensityUnits;    
  factor = NewDensity / OldDensity;

  densgrid  = OldDensity;
  ugrid     = CurrentGrid->BaryonField[Vel1Num][index];
  vgrid     = CurrentGrid->BaryonField[Vel2Num][index];
  wgrid     = CurrentGrid->BaryonField[Vel3Num][index];
  denssink  = (this->Mass - this->DeltaMass) / MassConversion;
  usink     = vel[1];
  vsink     = vel[2];
  wsink     = vel[3];
  drho      = this->DeltaMass / MassConversion;


  CurrentGrid->BaryonField[DensNum][index] *= factor;
  CurrentGrid->BaryonField[Vel1Num][index] = 
    (densgrid*ugrid - drho*ugrid) / (densgrid - drho);
  CurrentGrid->BaryonField[Vel2Num][index] = 
    (densgrid*vgrid - drho*vgrid) / (densgrid - drho);
  CurrentGrid->BaryonField[Vel3Num][index] = 
    (densgrid*wgrid - drho*wgrid) / (densgrid - drho);

  //this->Mass += this->DeltaMass;  // this is already done in Star_Accrete.C
  vel[1] = (denssink*usink + drho*ugrid) / (denssink + drho);
  vel[2] = (denssink*vsink + drho*vgrid) / (denssink + drho);
  vel[3] = (denssink*wsink + drho*wgrid) / (denssink + drho);

  /*
  fprintf(stderr, "star::SubtractAccretedMass:  OldDensity =%g, NewDensity =%g, factor =%g\n", 
	  OldDensity, NewDensity, factor);
  fprintf(stderr, "star::SubtractAccretedMass:  vel_g[1] = %g -> %g, vel_p[1] = %g -> %g\n", 
	  ugrid, CurrentGrid->BaryonField[Vel1Num][index], usink, vel[1]);
  */

  if (MultiSpecies) {
    CurrentGrid->BaryonField[DeNum][index] *= factor;
    CurrentGrid->BaryonField[HINum][index] *= factor;
    CurrentGrid->BaryonField[HIINum][index] *= factor;
    CurrentGrid->BaryonField[HeINum][index] *= factor;
    CurrentGrid->BaryonField[HeIINum][index] *= factor;
    CurrentGrid->BaryonField[HeIIINum][index] *= factor;
  }
  if (MultiSpecies > 1) {
    CurrentGrid->BaryonField[HMNum][index] *= factor;
    CurrentGrid->BaryonField[H2INum][index] *= factor;
    CurrentGrid->BaryonField[H2IINum][index] *= factor;
  }
  if (MultiSpecies > 2) {
    CurrentGrid->BaryonField[DINum][index] *= factor;
    CurrentGrid->BaryonField[DIINum][index] *= factor;
    CurrentGrid->BaryonField[HIINum][index] *= factor;
    CurrentGrid->BaryonField[HDINum][index] *= factor;
  }

  /*  
  if (ZField == TRUE)
    CurrentGrid->BaryonField[ZNum][index] *= factor;
  */

  return SUCCESS;

}
