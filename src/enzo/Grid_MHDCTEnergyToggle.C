/***********************************************************************
/
/  GRID CLASS ()
/
/  written by: David Collins
/  date:       2004-2013
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
/  Originally, MHDCT was written to use conserved energy throughout, while the rest of the
/  sovlers in Enzo dealt with specific energy.  Converting from one to the other unfortunately
/  yeilds minor differences in the solution.  This is due to the time interpolation in InterpolateBoundaryFromParent. 
/  In order to maintain contact with older solutions, I have provided two sets of conversion routines.  
/  
/  In all cases, use MHDCTUseSpecificEnergy == TRUE (the default).  
/ 
/  For all cases, specific energy is used for initialization and IO, just like the rest of enzo.
/ 
/  MHDCT_ConvertEnergyToConservedS
/  MHDCT_ConvertEnergyToSpecificS
/   When MHDCTUseSpecificEnergy == TRUE, these two routines convert from specific to conserved and back
/   inside Grid_SolveMHDEquations in order to prepare the energy for the MHD solver.  This should be the
/   default usage.
/  
/  MHDCT_ConvertEnergyToSpecificC
/  MHDCT_ConvertEnergyToConservedC
/   When MHDCTUseSpecificEnergy == FALSE, these routines convert from conserved to specific.  This is done 
/   with an extra copy of the energy field, in order to reduce noise.  This is done for pressure computation
/   and IO.
/ 
/  MHDCT_EnergyToggle
/   When MHDCTUseSpecificEnergy == FALSE, this converts the specific energy used in initialization to 
/   conserved.
/ 
************************************************************************/

#include <stdio.h>
#include "ErrorExceptions.h"
#include "EnzoTiming.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"

int dbg = FALSE;

int  MHDCT_EnergyToggle(HierarchyEntry &TopGrid, TopGridData &MetaData, ExternalBoundary *Exterior, LevelHierarchyEntry *LevelArray[]){
    //Converts the TotalEnergy field from specific to conserved.
    //This is a stopgap to ensure the code can be used while dcollins tracks down a problem
    //in a solver.

    if( MHDCTUseSpecificEnergy )
        return SUCCESS;
    LevelHierarchyEntry *Temp;
    MHDCTUseSpecificEnergy = TRUE; //just this once...
    for (int i = 0; i < MAX_DEPTH_OF_HIERARCHY; i++){
        Temp = LevelArray[i];
        while (Temp != NULL) {
          Temp->GridData->MHDCT_ConvertEnergyToConservedS();
          Temp = Temp->NextGridThisLevel;
        }
    }
    MHDCTUseSpecificEnergy = FALSE;
    return SUCCESS;
}

int grid::MHDCT_ConvertEnergyToConservedS()
{

    if ( HydroMethod != MHD_Li ||  EquationOfState != 0 || ProcessorNumber != MyProcessorNumber || MHDCTUseSpecificEnergy == FALSE)
        return SUCCESS;

    int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum, B1Num, B2Num, B3Num, size=1;
    if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                         Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
      ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
    }

    for (int dim = 0; dim < GridRank; dim++)
        size *= GridDimension[dim];
    for (int i=0;i<size;i++){
        BaryonField[TENum][i] *= BaryonField[DensNum][i];
    }
    
    return SUCCESS;
}
    

int grid::MHDCT_ConvertEnergyToConservedC()
{

    if ( HydroMethod != MHD_Li ||  EquationOfState != 0 || ProcessorNumber != MyProcessorNumber || MHDCTUseSpecificEnergy == TRUE )
        return SUCCESS;

    int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum, B1Num, B2Num, B3Num, size=1;
    if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                         Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
      ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
    }

    //Conversion to specific uses a copied temporary variable.
    delete [] BaryonField[TENum];
    BaryonField[TENum] = MHDCT_temp_conserved_energy;
    MHDCT_temp_conserved_energy= NULL;

    return SUCCESS;
}

int grid::MHDCT_ConvertEnergyToSpecificC()
{
    //Stores BaryonField[TENum] (conserved). in MHDCT_temp_energy
    //Creates new BaryonField[TENum] (specific) for output and computation of pressure
    
    if ( HydroMethod != MHD_Li ||  EquationOfState != 0 || ProcessorNumber != MyProcessorNumber || MHDCTUseSpecificEnergy == TRUE )
        return SUCCESS;

    int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum, B1Num, B2Num, B3Num, size=1;
    if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                     Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) 
        ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
    
    for (int dim = 0; dim < GridRank; dim++)
        size *= GridDimension[dim];

    MHDCT_temp_conserved_energy = BaryonField[TENum];
    BaryonField[TENum] = new float[size];
    for (int i=0; i<size; i++)
        BaryonField[TENum][i] = MHDCT_temp_conserved_energy[i]/BaryonField[DensNum][i];

    return SUCCESS;
}

int grid::MHDCT_ConvertEnergyToSpecificS()
{
    //Stores BaryonField[TENum] (conserved). in MHDCT_temp_energy
    //Creates new BaryonField[TENum] (specific) for output and computation of pressure
    
    if ( HydroMethod != MHD_Li ||  EquationOfState != 0 || ProcessorNumber != MyProcessorNumber || MHDCTUseSpecificEnergy == FALSE )
        return SUCCESS;

    int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum, B1Num, B2Num, B3Num, size=1;
    if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                     Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) 
        ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
    
    for (int dim = 0; dim < GridRank; dim++)
        size *= GridDimension[dim];

    for (int i=0;i<size;i++)
        BaryonField[TENum][i]  /= BaryonField[DensNum][i];

    return SUCCESS;
}
