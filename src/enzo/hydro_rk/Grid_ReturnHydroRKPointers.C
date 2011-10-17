/***********************************************************************
/
/  GRID CLASS (RETURNS AN ARRAY OF POINTERS THAT ARE COMPATIBLE WITH 
/              THE HYDRO_RK SOLVERS)
/
/  written by: John Wise
/  date:       July, 2009
/  modified1:
/
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
#include "TopGridData.h"
#include "Grid.h"

int FindField(int field, int farray[], int numfields);

int grid::ReturnHydroRKPointers(float **Prim, bool ReturnMassFractions)
{

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if (NumberOfBaryonFields == 0)
    return SUCCESS;

  int i, n, dim, size, nfield, n0;
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  int B1Num, B2Num, B3Num, PhiNum;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;

  /* Add the physical quantities */

  if (HydroMethod == HD_RK) {
    this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				     Vel3Num, TENum);
    nfield = n0 = NEQ_HYDRO;
  }

  else if (HydroMethod == MHD_RK) {
    this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				     Vel3Num, TENum, B1Num, B2Num, B3Num, 
				     PhiNum);
    nfield = n0 = NEQ_MHD;
  };

  Prim[iden] = BaryonField[DensNum];
  Prim[ivx]  = BaryonField[Vel1Num];
  Prim[ivy]  = BaryonField[Vel2Num];
  Prim[ivz]  = BaryonField[Vel3Num];
  Prim[ietot]= BaryonField[TENum];
  if (DualEnergyFormalism)
    Prim[ieint] = BaryonField[GENum];

  if (HydroMethod == MHD_RK) {
    Prim[iBx] = BaryonField[B1Num];
    Prim[iBy] = BaryonField[B2Num];
    Prim[iBz] = BaryonField[B3Num];
    Prim[iPhi]= BaryonField[PhiNum];
  }
  /*
  printf("Physical Quantities: %"ISYM" %"ISYM"  %"ISYM" %"ISYM" %"ISYM"  %"ISYM"  %"ISYM" %"ISYM" %"ISYM" %"ISYM"\n", 
	 DensNum, GENum, Vel1Num, Vel2Num, 
	 Vel3Num, TENum, B1Num, B2Num, B3Num, 
	 PhiNum);
  */
  /* Add the species */

  if (MultiSpecies) {
    this->IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, 
				HMNum, H2INum, H2IINum, DINum, DIINum, HDINum);

    //Prim[nfield++] = BaryonField[DeNum];
    Prim[nfield++] = BaryonField[HINum];
    Prim[nfield++] = BaryonField[HIINum];
    Prim[nfield++] = BaryonField[HeINum];
    Prim[nfield++] = BaryonField[HeIINum];
    Prim[nfield++] = BaryonField[HeIIINum];

    if (MultiSpecies > 1) {
      Prim[nfield++] = BaryonField[HMNum];
      Prim[nfield++] = BaryonField[H2INum];
      Prim[nfield++] = BaryonField[H2IINum];
    }

    if (MultiSpecies > 2) {
      Prim[nfield++] = BaryonField[DINum];
      Prim[nfield++] = BaryonField[DIINum];
      Prim[nfield++] = BaryonField[HDINum];
    }

  } // ENDIF MultiSpecies

  /* Add the colours (NColor is determined in EvolveLevel_RK2) */  

  int SNColourNum, MetalNum, MetalIaNum, MBHColourNum, Galaxy1ColourNum, 
    Galaxy2ColourNum; 

  if (this->IdentifyColourFields(SNColourNum, MetalNum, MetalIaNum, MBHColourNum, 
				 Galaxy1ColourNum, Galaxy2ColourNum) == FAIL) {
    fprintf(stderr, "Error in grid->IdentifyColourFields.\n");
    return FAIL;
  }
  
  if (MetalNum != -1) {
    Prim[nfield++] = BaryonField[MetalNum];
    if (StarMakerTypeIaSNe)
      Prim[nfield++] = OldBaryonField[MetalIaNum];
    if (MultiMetals || TestProblemData.MultiMetals) {
      Prim[nfield++] = BaryonField[MetalNum+1];
      Prim[nfield++] = BaryonField[MetalNum+2];
    }
  }

  if (SNColourNum      != -1) Prim[nfield++] = BaryonField[SNColourNum];  
  /*   //##### These fields are currently not being used and only causing interpolation problems
  if (MBHColourNum     != -1) Prim[nfield++] = BaryonField[MBHColourNum];
  if (Galaxy1ColourNum != -1) Prim[nfield++] = BaryonField[Galaxy1ColourNum];
  if (Galaxy2ColourNum != -1) Prim[nfield++] = BaryonField[Galaxy2ColourNum];
  */

  /* Convert the species and color fields into mass fractions */

  for (dim = 0, size = 1; dim < GridRank; dim++)
    size *= GridDimension[dim];

  if (ReturnMassFractions)  
    for (n = n0; n < nfield; n++)
      for (i = 0; i < size; i++) 
	Prim[n][i] /= Prim[iden][i];

  //  fprintf(stdout, "grid::ReturnHydroRKPointers: nfield = %"ISYM"\n", nfield);  

  return SUCCESS;

}
