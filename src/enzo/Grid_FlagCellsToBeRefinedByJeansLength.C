/*****************************************************************************
 *                                                                           *
 * Copyright 2004 Greg Bryan                                                 *
 * Copyright 2004 Laboratory for Computational Astrophysics                  *
 * Copyright 2004 Board of Trustees of the University of Illinois            *
 * Copyright 2004 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  GRID CLASS (FLAG CELLS TO BE REFINED BY THE JEAN'S CRITERION)
/
/  written by: Greg Bryan
/  date:       February, 1998
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    number of flagged cells, or -1 on failure
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

/* function prototypes */

int GetUnits(float *DensityUnits, float *LengthUnits,
		      float *TemperatureUnits, float *TimeUnits,
		      float *VelocityUnits, FLOAT Time);


int grid::FlagCellsToBeRefinedByJeansLength()
{

  /* declarations */

  int i, dim;
#ifndef MHDCT
  float IsothermalSoundSpeed=0;
#endif //MHDCT 

  /* error check */

  if (FlaggingField == NULL) {
    fprintf(stderr, "Flagging Field is undefined.\n");
    return -1;
  }

  /* compute size */

  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
  
  /* Find fields: density, total energy, velocity1-3. */
  
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				       Vel3Num, TENum) == FAIL) {
    fprintf(stderr, "RefineByJeansLength: Error in IdentifyPhysicalQuantities.\n");
    return FAIL;
  }
  
  
  //dcc expanded to account for different unit systems.
  //All three do the same thing: CellWidth^2 > (JeansLength/SafetyFactor)^2
  //0-Cosmology: JeansLength^2 = Pi * BoltzmanConst / (G*ProtonMass) * Temp/Density
  //1-Isothermal: JeansLength^2 = IsothermalSoundSpeed^2*Pi/(G*Rho)
  //                              IsothermalSoundSpeed is a run time parameter.
  //                              Note that GravitationalConstant = 4 Pi G
  //2-Adiabatic: JeansLentgh^2 = Gamma Pi/G * Pressure/Density^2
  
  float *temperature = NULL;
  float *pressure = NULL;
  FLOAT CellWidthSquared = CellWidth[0][0]*CellWidth[0][0];
  FLOAT JLSquared;

  //<dbg>
  float DensityZone;
  //</dbg>
  float DensityUnits = 1, LengthUnits = 1, VelocityUnits, TimeUnits,
    TemperatureUnits;
#ifdef MHDCT
  switch( RefineByJeansLengthUnits ){
    
  case 0: //cosmology
#endif //MHDCT    
    /* Compute the temperature field. */
    
    temperature = new float[size];
    if (this->ComputeTemperatureField(temperature) == FAIL) {
      fprintf(stderr, "Error in grid->ComputeTemperature.\n");
      return -1;
    }
    
    /* Get density units. */

    if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
		 &TimeUnits, &VelocityUnits, Time) == FAIL) {
      fprintf(stderr, "Error in CosmologyGetUnits.\n");
      return FAIL;
    }
    
    /* Compute constant for Jean's length computation.
       l_j = sqrt((pi*k*T) / (G \rho m_p))  . */
    
    JLSquared = (double(3.14159*1.38e-16/6.67e-8)/
		       (double(DensityUnits)*double(1.67e-24))) /
      (double(LengthUnits)*double(LengthUnits));

    if (ProblemType == 60 || ProblemType == 61)
      JLSquared = double(4.0*3.14159*3.14159)/GravitationalConstant; //AK

    
    /* This is the safety factor to decrease the Jean's length by. */
    
    JLSquared /= POW(RefineByJeansLengthSafetyFactor, 2);
    
    /* Loop over grid. */
    
    for (i = 0; i < size; i++)
      if (CellWidthSquared > JLSquared*temperature[i]/BaryonField[DensNum][i])
	FlaggingField[i]++;
    
    /* clean up */
    
    delete temperature;
#ifdef MHDCT    
    break;
  case 1: //isothermal.  pretty easy.  Remember that GravitationalConstant = 4 Pi G.
    JLSquared = IsothermalSoundSpeed*IsothermalSoundSpeed*
      4*3.14159*3.14159/(GravitationalConstant);
    JLSquared /= RefineByJeansLengthSafetyFactor*RefineByJeansLengthSafetyFactor;
    for(i=0;i<size;i++){
      //fprintf(stderr, "dx2 %f jl2 %f d %f", CellWidthSquared, JLSquared, BaryonField[DensNum][i]);
      //      if( CellWidthSquared > JLSquared/BaryonField[DensNum][i] ){

      //<dbg>
      DensityZone = BaryonField[DensNum][i];
      if( DensityZone == 0.0 ) {
	fprintf(stderr,"RefingByJeansLength: Density == 0.  (That's not good.)\n");
      }
      //</dbg>
      if( CellWidthSquared > JLSquared/DensityZone ){
	//fprintf(stderr," ok \n");
	FlaggingField[i]++;
      }else{
	//fprintf(stderr," no \n");
      }
    }

    break;
  case 2: //adiabatic, code units.
    pressure = new float[size];
    int result;
    if (DualEnergyFormalism)
      result = this->ComputePressureDualEnergyFormalism(Time, pressure);
    else
      result = this->ComputePressure(Time, pressure);
    
    if (result == FAIL) 
      {fprintf(stderr, "Error in grid->ComputePressure, called in MHD_Athena\n");return FAIL;}

    JLSquared = Gamma*4*3.14159*3.14159/(GravitationalConstant);
    JLSquared /= RefineByJeansLengthSafetyFactor*RefineByJeansLengthSafetyFactor;
    for(i=0;i<size;i++)
      if( CellWidthSquared > 
	  JLSquared*pressure[i]/(BaryonField[DensNum][i]*BaryonField[DensNum][i] ) )
	FlaggingField[i]++;
    
    delete pressure;
    break;
  case 3: //5
    for(i=0;i<size;i++)
      if( BaryonField[DensNum][i] > 5/(64.*CellWidth[0][0]) ){
	FlaggingField[i]++;
      }
    break;
  }
#endif //MHDCT  
  /* Count number of flagged Cells. */

  int NumberOfFlaggedCells = 0;
  for (i = 0; i < size; i++) {
    FlaggingField[i] = (FlaggingField[i] >= 1)? 1 : 0;
    NumberOfFlaggedCells += FlaggingField[i];
  }

  return NumberOfFlaggedCells;

}
