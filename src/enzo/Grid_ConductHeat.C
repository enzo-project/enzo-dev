/***********************************************************************
/
/  GRID CLASS (Compute and apply thermal conduction)
/
/  written by:  David A. Ventimiglia and Brian O'Shea
/  date:        November, 2009
/
/  PURPOSE:  Calculates and applies heat fluxes due to thermal 
/  conduction
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#include <math.h> 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "phys_constants.h"
#include "CosmologyParameters.h"

int grid::ConductHeat(){

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0)
    return SUCCESS;


  this->DebugCheck("ConductHeat");

  // Some locals
  int size = 1, idx, i,j,k, Nsub=0; 
  int DensNum, TENum, GENum, Vel1Num, Vel2Num, Vel3Num;
  float *e, eold;
  float dtSubcycle, dtSoFar;

  for (int dim = 0; dim < GridRank; dim++) 
    size *= GridDimension[dim];

  float *dedt = new float[size];

  // Get internal energy field
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum) == FAIL) 
    {
      ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
    }

  if (UseMHD){
    iBx=FindField(Bfield1, FieldType, NumberOfBaryonFields);
    iBy=FindField(Bfield2, FieldType, NumberOfBaryonFields);
    iBz=FindField(Bfield3, FieldType, NumberOfBaryonFields);  
  }

  // figure out what is actually internal energy and set it.  If we
  // use PPM, convert total energy into internal energy regardless, 
  // because it has to be updated to keep everything consistent.

  if (HydroMethod==Zeus_Hydro){
    e = BaryonField[TENum];  // 'total energy' is really gas internal energy in zeus
  } else if (HydroMethod==PPM_DirectEuler) {
    if(DualEnergyFormalism){
      e = BaryonField[GENum];  // if DEF, gas energy is really gas energy.
    }
    else {
      // energy will be calculated below inside time loop.
      e = new float[size];
    } 
  } else if (HydroMethod == MHD_RK) {
    e = new float[size];
  } else {  // fails for PPM_LR, HD_RK
    ENZO_FAIL("Error in Grid::ConductHeat - your Hydro/MHD method is not supported!\n");
  }

  // Sub-cycle, computing and applying heat

  // dtSubcycle = timestep of this subcycle
  // dtSoFar = overall timestep taken in heat equation
  // dtFixed = fixed timestep for the entire level.

  dtSoFar = 0.0;

  Nsub=0;  // number of subcycles

  while(dtSoFar < dtFixed){

    // convert total energy into internal energy
    if(HydroMethod==PPM_DirectEuler && DualEnergyFormalism==0) {
      for (i = 0; i < size; i++) {
	e[i] = BaryonField[TENum][i] - 0.5*POW(BaryonField[Vel1Num][i], 2.0);
	if(GridRank > 1)
	  e[i] -= 0.5*POW(BaryonField[Vel2Num][i], 2.0);
	if(GridRank > 2)
	  e[i] -= 0.5*POW(BaryonField[Vel3Num][i], 2.0);
      }
    }

    if(HydroMethod==MHD_RK){  // convert total energy into internal energy: subract off kinetic, magnetic energy
      
      for (i = 0; i < size; i++) {
	e[i] = BaryonField[TENum][i] - 0.5*POW(BaryonField[Vel1Num][i], 2.0);
	if(GridRank > 1)
	  e[i] -= 0.5*POW(BaryonField[Vel2Num][i], 2.0);
	if(GridRank > 2)
	  e[i] -= 0.5*POW(BaryonField[Vel3Num][i], 2.0);
	
	e[i] -= 0.5*(POW(BaryonField[iBx][i],2.0) + POW(BaryonField[iBy][i],2.0) + POW(BaryonField[iBz][i],2.0))/BaryonField[DensNum][i];
      }
    }

    // compute this subcycle timestep
    if (this->ComputeConductionTimeStep(dtSubcycle) == FAIL)
      ENZO_FAIL("Error in ComputeConductionTimeStep.\n");

    // make sure we don't extend past dtFixed
    dtSubcycle = min(dtSubcycle, dtFixed-dtSoFar);

    // compute de/dt for each cell.
    if (this->ComputeHeat(dedt) == FAIL) {
      ENZO_FAIL("Error in ComputeHeat.");
    }

    /* We loop over the whole grid, except for the first and last 
       cell in any given 1-D sweep.*/
    int GridStart[] = {0, 0, 0}, GridEnd[] = {0, 0, 0};

    for (int dim = 0; dim<GridRank; dim++) {
      GridStart[dim] = 0;
      GridEnd[dim] = GridDimension[dim] - 1;
    }

    for (k = GridStart[2]; k <= GridEnd[2]; k++) 
      for (j = GridStart[1]; j <= GridEnd[1]; j++) 
  	for (i = GridStart[0]; i <= GridEnd[0]; i++) {

	  idx = ELT(i,j,k);

	  eold=e[idx];

  	  e[idx] += dedt[idx]*dtSubcycle;

	  if (e[idx]<0) {
	    ENZO_VFAIL("Grid_ConductHeat: e=%g dedt=%g (eold=%g) i,j,k=%d,%d,%d  dtFixed,dtSub: %e, %e\n", 
		       e[ELT(i,j,k)],dedt[idx], eold, i,j,k,
		       dtFixed, dtSubcycle);
	  }

	  // energy needs to be updated here if hydro is PPM or MHD
	  if(HydroMethod != Zeus_Hydro){

	    BaryonField[TENum][idx] = e[idx] + 0.5*POW(BaryonField[Vel1Num][idx], 2.0);
	    if(GridRank > 1)
	      BaryonField[TENum][idx] += 0.5*POW(BaryonField[Vel2Num][idx], 2.0);
	    if(GridRank > 2)
	      BaryonField[TENum][idx] += 0.5*POW(BaryonField[Vel3Num][idx], 2.0);

	    if(UseMHD)
	      BaryonField[TENum][idx] += 0.5*(POW(BaryonField[iBx][idx],2.0) + 
					      POW(BaryonField[iBy][idx],2.0) + 
					      POW(BaryonField[iBz][idx],2.0))/BaryonField[DensNum][idx];

	  } // if(HydroMethod != Zeus_Hydro)

	} // triple for loop

    // increment timestep
    dtSoFar += dtSubcycle;
    Nsub++;

  } // while(dtSoFar < dtFixed)

  if(debug1) printf("Grid::ConductHeat:  Nsubcycles = %"ISYM"\n",Nsub);

  if((HydroMethod==PPM_DirectEuler && DualEnergyFormalism==0) || (HydroMethod==MHD_RK))
    delete [] e;

  delete [] dedt;

  return SUCCESS;
  
}

