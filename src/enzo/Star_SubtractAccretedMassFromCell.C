/***********************************************************************
/
/  SUBTRACT ACCRETED MASS FROM CELLS
/
/  written by: Ji-hoon Kim
/  date:       September, 2009
/  modified1: 
/
/  PURPOSE: This routine subtracts the accreted gas mass out of the grid
/           for Star particle type BlackHole..
/           Note that accretion_rate is calculted in Star_CalculateMassAccretion.C
/           but DeltaMass is calculated in Star_Accrete.C.
/           At the moment, this method is used only for BlackHole;
/           for MBH, the job is done in Grid_SubtractAccretedMassFromSphere.C
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

int FindField(int field, int farray[], int numfields);

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int Star::SubtractAccretedMassFromCell(void)
{

  /* Check if the star type is correct */

  if ((this->type != BlackHole) || 
      (this->CurrentGrid == NULL))
    return SUCCESS;

  int dim, igrid[MAX_DIMENSION], index, size;
  double Msun = 1.989e33;
  FLOAT time = CurrentGrid->Time;

  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, time);

  /* Find fields: density, total energy, velocity1-3. */

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (CurrentGrid->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }
  
  /* Find Multi-species fields. */

  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
  if (MultiSpecies) 
    if (CurrentGrid->IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, 
				    HeIIINum, HMNum, H2INum, H2IINum, DINum, 
				    DIINum, HDINum) == FAIL) {
      ENZO_FAIL("Error in grid->IdentifySpeciesFields.\n");
    }

  /* Find Metallicity or SNColour field and set flag. */

  int SNColourNum, MetalNum, MBHColourNum, Galaxy1ColourNum, Galaxy2ColourNum,
    MetalIaNum, MetalIINum;
  int MetallicityField = FALSE;

  if (CurrentGrid->IdentifyColourFields(SNColourNum, MetalNum, MetalIaNum, MetalIINum, 
              MBHColourNum, Galaxy1ColourNum, Galaxy2ColourNum) == FAIL)
    ENZO_FAIL("Error in grid->IdentifyColourFields.\n");

  /* Now let's start working! */

  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    size *= CurrentGrid->GridDimension[dim];
    igrid[dim] = (int) ((pos[dim] - CurrentGrid->GridLeftEdge[dim]) /
			CurrentGrid->CellWidth[0][0]);
  }

  index = 
    ((igrid[2] + CurrentGrid->GridStartIndex[2]) * CurrentGrid->GridDimension[1] + 
     igrid[1] + CurrentGrid->GridStartIndex[1]) * CurrentGrid->GridDimension[0] + 
     igrid[0] + CurrentGrid->GridStartIndex[0];

  float MassConversion = (float) (pow(LengthUnits * CurrentGrid->CellWidth[0][0], 3.0)
				  * double(DensityUnits) / Msun);
  float densgrid, ugrid, vgrid, wgrid, denssink, usink, vsink, wsink, drho;
  double OldDensity, NewDensity, factor;


  /* Subtract accreted mass from the grids, and calculate new densities */
  
  OldDensity = CurrentGrid->BaryonField[DensNum][index];
  NewDensity = OldDensity - this->DeltaMass / MassConversion;  
  factor = NewDensity / OldDensity;

  denssink  = (float(this->Mass) - this->DeltaMass) / MassConversion; //check below
  usink     = this->vel[0];
  vsink     = this->vel[1];
  wsink     = this->vel[2];
  densgrid  = OldDensity;
  ugrid     = CurrentGrid->BaryonField[Vel1Num][index];
  vgrid     = CurrentGrid->BaryonField[Vel2Num][index];
  wgrid     = CurrentGrid->BaryonField[Vel3Num][index];
  drho      = this->DeltaMass / MassConversion;


  /* Modify density and velocity fields, for both the particle and the grid */

  //this->Mass += this->DeltaMass; //this is already done in Star_Accrete.C
  this->vel[0] = (denssink*usink + drho*ugrid) / (denssink + drho);
  this->vel[1] = (denssink*vsink + drho*vgrid) / (denssink + drho);
  this->vel[2] = (denssink*wsink + drho*wgrid) / (denssink + drho);

  CurrentGrid->BaryonField[DensNum][index] *= factor;
  //CurrentGrid->BaryonField[Vel1Num][index] = (densgrid*ugrid - drho*ugrid) / (densgrid - drho);
  //                                         = ugrid;  //velocity of the grids will be unchanged! 


//  fprintf(stdout, "star::SubtractAccretedMass[%d]:  DeltaMass = %e, OldDensity =%e, NewDensity =%e, factor =%e\n", 
//	  this->Identifier, this->DeltaMass, OldDensity, NewDensity, factor); 
//  fprintf(stdout, "star::SubtractAccretedMass[%d]:  vel_p[1] = %g -> %g\n", 
//	  this->Identifier, vsink, vel[1]); 


  /* Update species and colour fields */

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

  if (MetalNum >= 0)
    CurrentGrid->BaryonField[MetalNum][index] *= factor;
  if (MetalIaNum >= 0)
    CurrentGrid->BaryonField[MetalIaNum][index] *= factor;
  if (SNColourNum >= 0)
    CurrentGrid->BaryonField[SNColourNum][index] *= factor;
  if (MBHColourNum > 0)
    CurrentGrid->BaryonField[MBHColourNum][index] *= factor;    

  return SUCCESS;

}
