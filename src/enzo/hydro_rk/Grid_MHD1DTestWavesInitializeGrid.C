/***********************************************************************
/
/  GRID CLASS (INITIALIZE MHD 1D WAVE TEST)
/
/  written by: J. S. Oishi
/  date:       April 2011
/  modified1:
/
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
#include "Grid.h"
#include "CosmologyParameters.h"

int grid::MHD1DTestWavesInitializeGrid(float rhol, 
				  float vxl,
				  float vyl,
				  float vzl,
				  float etotl,
				  float Bxl,
				  float Byl,
				  float Bzl)
{  

  /* create fields */
  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  FieldType[NumberOfBaryonFields++] = Velocity2;
  FieldType[NumberOfBaryonFields++] = Velocity3;
  FieldType[NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism) {
    FieldType[NumberOfBaryonFields++] = InternalEnergy;
  }
  if (HydroMethod == MHD_RK) {
    FieldType[NumberOfBaryonFields++] = Bfield1;
    FieldType[NumberOfBaryonFields++] = Bfield2;
    FieldType[NumberOfBaryonFields++] = Bfield3;
    FieldType[NumberOfBaryonFields++] = PhiField;
  }

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }

  int activesize = 1, dim;

  for (dim = 0; dim < GridRank; dim++)
    activesize *= (GridDimension[dim] - 2*NumberOfGhostZones);
  
  this->AllocateGrids();

  
  FLOAT x, rhobg, bxbg, bybg, bzbg, etotbg, pgasbg, B2 = 0., ampl = 1.e-6;
  int i;

  rhobg = 1.; 
  if (HydroMethod == MHD_RK) {
    bxbg = 1.;
    bybg = sqrt(2.);
    bzbg = 0.5;
    B2 = bxbg * bxbg + bybg * bybg + bzbg * bzbg;
  }
  pgasbg = 1/Gamma;
  etotbg = rhobg / ((Gamma - 1.0)*rhobg) + 0.5 * B2/rhobg;

  for (i = 0; i < GridDimension[0]; i++) {

    /* Compute position */
    x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];

    // the following is for a sine wave
    BaryonField[iden ][i] = rhobg + rhol * ampl * sin(2*PI*x);
    BaryonField[ivx  ][i] = vxl * ampl * sin(2*PI*x);
    BaryonField[ivy  ][i] = vyl * ampl * sin(2*PI*x);
    BaryonField[ivz  ][i] = vzl * ampl * sin(2*PI*x);
    BaryonField[ietot][i] = etotbg + etotl * ampl * sin(2*PI*x);
    if (DualEnergyFormalism) {
      fprintf(stderr, "MHD1DTestWavesInitializeGrids: Dual Energy NOT IMPLEMENTED!!\n");
      //BaryonField[ieint][i] = pl / ((Gamma-1.0)*rhol);
    }
    if (HydroMethod == MHD_RK) {
      BaryonField[iBx  ][i] = bxbg + Bxl * ampl * sin(2*PI*x);
      BaryonField[iBy  ][i] = bybg + Byl * ampl * sin(2*PI*x);
      BaryonField[iBz  ][i] = bzbg + Bzl * ampl * sin(2*PI*x);
      BaryonField[iPhi ][i] = 0.0;
    }
  }

  // the following is for a square wave test

    // BaryonField[iden ][i] = rhobg;
    // BaryonField[ivx  ][i] = 0.;
    // BaryonField[ivy  ][i] = 0.;
    // BaryonField[ivz  ][i] = 0.;
    // BaryonField[ietot][i] = etotbg;
    // if (HydroMethod == MHD_RK) {
    //   BaryonField[iBx  ][i] = bxbg;
    //   BaryonField[iBy  ][i] = bybg;
    //   BaryonField[iBz  ][i] = bzbg;
    //   BaryonField[iPhi ][i] = 0.0;
    // }
    
    // if ((x < 0.75) && (x > 0.25)) {
    //   BaryonField[iden ][i] += rhol * ampl;
    //   BaryonField[ivx  ][i] += vxl * ampl;
    //   BaryonField[ivy  ][i] += vyl * ampl;
    //   BaryonField[ivz  ][i] += vzl * ampl;
    //   BaryonField[ietot][i] += etotl * ampl;
    //   if (DualEnergyFormalism) {
    //     fprintf(stderr, "MHD1DTestWavesInitializeGrids: Dual Energy NOT IMPLEMENTED!!\n");
    //     //BaryonField[ieint][i] = pl / ((Gamma-1.0)*rhol);
    //   }
    //   if (HydroMethod == MHD_RK) {
    //     BaryonField[iBx  ][i] += Bxl * ampl;
    //     BaryonField[iBy  ][i] += Byl * ampl;
    //     BaryonField[iBz  ][i] += Bzl * ampl;
    //   }
  //}

  return SUCCESS;
}
