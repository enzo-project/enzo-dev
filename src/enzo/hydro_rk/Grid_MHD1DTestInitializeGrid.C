/***********************************************************************
/
/  GRID CLASS (INITIALIZE MHD 1D TEST)
/
/  written by: Peng Wang
/  date:       June, 2007
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

int grid::MHD1DTestInitializeGrid(float rhol, float rhor,
				  float vxl,  float vxr,
				  float vyl,  float vyr,
				  float vzl,  float vzr,
				  float pl,   float pr,
				  float Bxl,  float Bxr,
				  float Byl,  float Byr,
				  float Bzl,  float Bzr)
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
  if (UseMHD) {
    FieldType[NumberOfBaryonFields++] = Bfield1;
    FieldType[NumberOfBaryonFields++] = Bfield2;
    FieldType[NumberOfBaryonFields++] = Bfield3;
  }
  if( HydroMethod == MHD_RK ){
    FieldType[NumberOfBaryonFields++] = PhiField;
  }

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }

  int size = 1, activesize = 1, dim;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  for (dim = 0; dim < GridRank; dim++)
    activesize *= (GridDimension[dim] - 2*NumberOfGhostZones);
  
  /*
  int field;
  for (field = 0; field < NumberOfBaryonFields; field++) {
    if (BaryonField[field] == NULL) {
      BaryonField[field] = new float[size];
    }
  }

  divB = new float[activesize];
  for (int dim = 0; dim < 3; dim++) {
    gradPhi[dim] = new float[activesize];
  }

  
  for (int dim = GridRank; dim < 3; dim++) {
    for (int n = 0; n < activesize; n++) {
      gradPhi[dim][n] = 0.0;
    }
  }
  */
  this->AllocateGrids();

  
  /* transform pressure to total energy */
  float etotl, etotr, v2, B2=0;
  v2 = vxl * vxl + vyl * vyl + vzl*vzl;
  if (HydroMethod == MHD_RK || UseMHDCT) B2 = Bxl*Bxl + Byl*Byl + Bzl*Bzl;
  etotl = pl / ((Gamma-1.0)*rhol) + 0.5*v2 + 0.5*B2/rhol; 

  v2 = vxr * vxr + vyr * vyr + vzr*vzr;
  if (HydroMethod == MHD_RK || UseMHDCT) B2 = Bxr*Bxr + Byr*Byr + Bzr*Bzr; else B2 = 0.;
  etotr = pr / ((Gamma-1.0)*rhor) + 0.5*v2 + 0.5*B2/rhor;

  FLOAT x;
  int i;

#ifdef NOUSE
  /* Initial condition for Gaussian profile for the resistivity test */

  for (i = 0; i < GridDimension[0]; i++) {

    /* Compute position */

    x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
    if (x <= 0.5) {
      BaryonField[iden ][i] = 1.0;
      BaryonField[ivx  ][i] = 0.0;
      BaryonField[ivy  ][i] = 0.0;
      BaryonField[ivz  ][i] = 0.0;
      BaryonField[ietot][i] = 1.0;
      if (UseMHD) {
        BaryonField[iBx  ][i] = 0.0;
        BaryonField[iBy  ][i] = 0.0;
        BaryonField[iBz  ][i] = 1.0/sqrt(4.0*Pi*0.01)*exp(-(x-0.5)*(x-0.5)/(4.0*0.01));
        if( HydroMethod == MHD_RK ){
          BaryonField[iPhi ][i] = 0.0;
        }
      }
    } else {
      BaryonField[iden ][i] = 1.0;
      BaryonField[ivx  ][i] = 0.0;
      BaryonField[ivy  ][i] = 0.0;
      BaryonField[ivz  ][i] = 0.0;
      BaryonField[ietot][i] = 1.0;
      if (UseMHD) {
        BaryonField[iBx  ][i] = 0.0;
        BaryonField[iBy  ][i] = 0.0;
        BaryonField[iBz  ][i] = 1.0/sqrt(4.0*Pi*0.01)*exp(-(x-0.5)*(x-0.5)/(4.0*0.01));
        if( HydroMethod == MHD_RK ){
            BaryonField[iPhi ][i] = 0.0;
        }
      }
    }
  }
#endif


  for (i = 0; i < GridDimension[0]; i++) {

    /* Compute position */

    x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];  // put inerface in the middle of the Domain
    if (x <= DomainLeftEdge[0]+0.5*(DomainRightEdge[0]-DomainLeftEdge[0])) { 
      BaryonField[iden ][i] = rhol;
      BaryonField[ivx  ][i] = vxl;
      BaryonField[ivy  ][i] = vyl;
      BaryonField[ivz  ][i] = vzl;
      BaryonField[ietot][i] = etotl;
      if (DualEnergyFormalism) {
	BaryonField[ieint][i] = pl / ((Gamma-1.0)*rhol);
      }
      if (UseMHD) {
        BaryonField[iBx  ][i] = Bxl;
        BaryonField[iBy  ][i] = Byl;
        BaryonField[iBz  ][i] = Bzl;
      }
      if( HydroMethod == MHD_RK) {
            BaryonField[iPhi ][i] = 0.0;
      }
      if( UseMHDCT ){
          MagneticField[0][i] = Bxl;
          MagneticField[1][i] = Byl;
          MagneticField[2][i] = Bzl;
          for ( int k=0; k<2; k++)
          for ( int j=0; j<2; j++){
              MagneticField[1][i + MagneticDims[1][0]*(j*MagneticDims[1][1]*k)] = Byl;
              MagneticField[2][i + MagneticDims[2][0]*(j*MagneticDims[2][1]*k)] = Bzl;
          }

      }
    } else {
      BaryonField[iden ][i] = rhor;
      BaryonField[ivx  ][i] = vxr;
      BaryonField[ivy  ][i] = vyr;
      BaryonField[ivz  ][i] = vzr;
      BaryonField[ietot][i] = etotr;
      if (DualEnergyFormalism) {
	BaryonField[ieint][i] = pr / ((Gamma-1.0)*rhor);
      }
      if (UseMHD) {
        BaryonField[iBx  ][i] = Bxr;
        BaryonField[iBy  ][i] = Byr;
        BaryonField[iBz  ][i] = Bzr;
      }
      if( HydroMethod == MHD_RK ){
        BaryonField[iPhi ][i] = 0.0;
      }
      if( UseMHDCT ){
          MagneticField[0][i] = Bxr;
          MagneticField[1][i] = Byr;
          MagneticField[2][i] = Bzr;
          for ( int k=0; k<2; k++)
          for ( int j=0; j<2; j++){
              MagneticField[1][i + MagneticDims[1][0]*(j*MagneticDims[1][1]*k)] = Byr;
              MagneticField[2][i + MagneticDims[2][0]*(j*MagneticDims[2][1]*k)] = Bzr;
          }

      }
    }
  }


  return SUCCESS;
}
