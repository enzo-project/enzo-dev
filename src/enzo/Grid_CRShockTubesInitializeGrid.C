#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "phys_constants.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
		      float *TemperatureUnits, float *TimeUnits,
		      float *VelocityUnits, FLOAT Time);

int grid::CRShockTubesInitializeGrid(   float x0, 
					float rhol, float rhor,
					float vxl,  float vxr,
					float vyl,  float vyr,
					float vzl,  float vzr,
					float pl,   float pr,
					float crl,  float crr,
					float bxl,  float bxr, 
					float byl,  float byr, 
					float bzl,  float bzr
					)
{  

  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  FieldType[NumberOfBaryonFields++] = Velocity2;
  FieldType[NumberOfBaryonFields++] = Velocity3;
  FieldType[NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism) 
    FieldType[NumberOfBaryonFields++] = InternalEnergy;

  if (HydroMethod == MHD_RK){
    FieldType[NumberOfBaryonFields++] = Bfield1;
    FieldType[NumberOfBaryonFields++] = Bfield1;
    FieldType[NumberOfBaryonFields++] = Bfield3;
    FieldType[NumberOfBaryonFields++] = PhiField;
  }
  FieldType[NumberOfBaryonFields++] = CRDensity;

  
  int iCRD = FindField(CRDensity, FieldType, NumberOfBaryonFields);

  
  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }

  int size = 1, activesize = 1, dim;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  for (dim = 0; dim < GridRank; dim++)
    activesize *= (GridDimension[dim] - 2*NumberOfGhostZones);
  
  int field;
  for (field = 0; field < NumberOfBaryonFields; field++)
    if (BaryonField[field] == NULL)
      BaryonField[field] = new float[size];
  
  /* transform pressure to total energy */

  float etotl, etotr, v2, B2;
  v2 = vxl * vxl + vyl * vyl + vzl * vzl;
  etotl = pl / ((Gamma-1.0)*rhol) + 0.5*v2;

  v2 = vxr * vxr + vyr * vyr + vzr * vzr;
  etotr = pr / ((Gamma-1.0)*rhor) + 0.5*v2;

  FLOAT x;
  int i;
  for (i = 0; i < GridDimension[0]; i++) {

    x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];

    if (x <= x0) {
      BaryonField[iden ][i] = rhol;
      BaryonField[ivx  ][i] = vxl;
      BaryonField[ivy  ][i] = vyl;
      BaryonField[ivz  ][i] = vzl;
      BaryonField[ietot][i] = etotl;
      BaryonField[iCRD ][i] = crl;
      if (DualEnergyFormalism) {
	BaryonField[ieint][i] = pl / ((Gamma-1.0)*rhol);
      }
      if (HydroMethod == MHD_RK) {
	B2 = bxl*bxl + byl*byl + bzl*bzl; 
	BaryonField[ietot][i] += 0.5*B2/rhol;
	BaryonField[iBx][i] = bxl;
	BaryonField[iBy][i] = byl;
	BaryonField[iBz][i] = bzl;
	BaryonField[iPhi][i] = 0.0; 
      }
    } else {
      BaryonField[iden ][i] = rhor;
      BaryonField[ivx  ][i] = vxr;
      BaryonField[ivy  ][i] = vyr;
      BaryonField[ivz  ][i] = vzr;
      BaryonField[ietot][i] = etotr;
      BaryonField[iCRD ][i] = crr;
      if (DualEnergyFormalism) {
	BaryonField[ieint][i] = pr / ((Gamma-1.0)*rhor);
      }
      if (HydroMethod == MHD_RK) {
        B2 = bxr*bxr + byr*byr + bzr*bzr;
	BaryonField[ietot][i] += 0.5*B2/rhor;
	BaryonField[iBx][i] = bxr;
	BaryonField[iBy][i] = byr;
	BaryonField[iBz][i] = bzr;
	BaryonField[iPhi][i] = 0.0; 
      }
    }
  }

  return SUCCESS;
}

/* Version to specify three regions */

int grid::CRShockTubesInitializeGrid(   float x0,   float x1,
					float rhol, float rhor, float rhoc,
					float vxl,  float vxr,  float vxc,
					float vyl,  float vyr,  float vyc,
					float vzl,  float vzr,  float vzc,
					float pl,   float pr,   float pc,
					float crl,  float crr,  float crc, 
					float bxl,  float bxr,  float bxc, 
					float byl,  float byr,  float byc, 
					float bzl,  float bzr,  float bzc
					)
{  

  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  FieldType[NumberOfBaryonFields++] = Velocity2;
  FieldType[NumberOfBaryonFields++] = Velocity3;
  FieldType[NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism) {
    FieldType[NumberOfBaryonFields++] = InternalEnergy;
  }
  if (HydroMethod == MHD_RK){
    FieldType[NumberOfBaryonFields++] = Bfield1;
    FieldType[NumberOfBaryonFields++] = Bfield2;
    FieldType[NumberOfBaryonFields++] = Bfield3;
    FieldType[NumberOfBaryonFields++] = PhiField;
  }
  FieldType[NumberOfBaryonFields++] = CRDensity;
  

  int iCRD = FindField( CRDensity , FieldType, NumberOfBaryonFields);

  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }

  int size = 1, activesize = 1, dim;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  for (dim = 0; dim < GridRank; dim++)
    activesize *= (GridDimension[dim] - 2*NumberOfGhostZones);
  
  int field;
  for (field = 0; field < NumberOfBaryonFields; field++)
    if (BaryonField[field] == NULL)
      BaryonField[field] = new float[size];

  /* transform pressure to total energy */

  float etotl, etotr, etotc, v2, B2;
  v2 = vxl * vxl + vyl * vyl + vzl * vzl;
  etotl = pl / ((Gamma-1.0)*rhol) + 0.5*v2;

  v2 = vxr * vxr + vyr * vyr + vzr * vzr;
  etotr = pr / ((Gamma-1.0)*rhor) + 0.5*v2;

  v2 = vxc * vxc + vyc * vyc + vzc * vzc;
  etotc = pc / ((Gamma-1.0)*rhoc) + 0.5*v2;

  

  FLOAT x;
  int i;
  for (i = 0; i < GridDimension[0]; i++) {

    x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];

    if (x <= x0) {
      BaryonField[iden ][i] = rhol;
      BaryonField[ivx  ][i] = vxl;
      BaryonField[ivy  ][i] = vyl;
      BaryonField[ivz  ][i] = vzl;
      BaryonField[ietot][i] = etotl;
      BaryonField[iCRD ][i] = crl;
      if (DualEnergyFormalism) {
	BaryonField[ieint][i] = pl / ((Gamma-1.0)*rhol);
      }
      if (HydroMethod == MHD_RK) {
	B2 = bxl*bxl + byl*byl + bzl*bzl;
        BaryonField[ietot][i] += 0.5*B2/rhol;
        BaryonField[iBx][i] = bxl;
        BaryonField[iBy][i] = byl;
        BaryonField[iBz][i] = bzl;
        BaryonField[iPhi][i] = 0.0;
      }
    } else if (x <= x1) {
      BaryonField[iden ][i] = rhoc;
      BaryonField[ivx  ][i] = vxc;
      BaryonField[ivy  ][i] = vyc;
      BaryonField[ivz  ][i] = vzc;
      BaryonField[ietot][i] = etotc;
      BaryonField[iCRD ][i] = crc;
      if (DualEnergyFormalism) {
	BaryonField[ieint][i] = pc / ((Gamma-1.0)*rhoc);
      }
      if (HydroMethod == MHD_RK) {
	B2 = bxc*bxc + byc*byc + bzc*bzc;
        BaryonField[ietot][i] += 0.5*B2/rhoc;
        BaryonField[iBx][i] = bxc;
        BaryonField[iBy][i] = byc;
        BaryonField[iBz][i] = bzc;
        BaryonField[iPhi][i] = 0.0;
      }
    }
    else {
      BaryonField[iden ][i] = rhor;
      BaryonField[ivx  ][i] = vxr;
      BaryonField[ivy  ][i] = vyr;
      BaryonField[ivz  ][i] = vzr;
      BaryonField[ietot][i] = etotr;
      BaryonField[iCRD ][i] = crr;
      if (DualEnergyFormalism) {
	BaryonField[ieint][i] = pr / ((Gamma-1.0)*rhor);
      }
      if (HydroMethod == MHD_RK) {
	B2 = bxr*bxr + byr*byr + bzr*bzr;
        BaryonField[ietot][i] += 0.5*B2/rhor;
        BaryonField[iBx][i] = bxr;
        BaryonField[iBy][i] = byr;
        BaryonField[iBz][i] = bzr;
        BaryonField[iPhi][i] = 0.0;
      }
    }
    if (debug)
      printf("%"FSYM"\t%"FSYM"\t%"FSYM"\n", x , BaryonField[ieint][i], BaryonField[iCRD][i] );
  }

  // --------------- FOR DIFFUSION PROBLEM --
  // This sets of a linear ramp

  if( crc == 123.4 ){
    int i0,i1,i2;
    i0 = GridDimension[0] / 4;
    i1 = i0*2; i2 = i0*3;
	
    double mCR = ( crc - crl ) / ( i2 - i1 );
	
    for (i = 0; i < GridDimension[0]; i++) {
      if( i < i0 )
	BaryonField[iCRD][i] = crl;
      else if( i < i1 )
	BaryonField[iCRD][i] = crl + mCR*(i-i0);
      else if( i < i2 )
	BaryonField[iCRD][i] = crc - mCR*(i-i1);
      else
	BaryonField[iCRD][i] = crl;
    } // end i for
  } // end if

  /*--- FOR NEW CR DIFFUSION PROBLEM --- */
  // This sets up a Gaussian

  if ( crc == 567.8 ){
    double t0 = 1.0;	// STARTING TIME
    for (i = 0; i < GridDimension[0]; i++) {
      x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];		
      x = x - x0;	// translate gaussian to center
      BaryonField[iCRD][i] = 1.0/sqrt(4.0*pi*CRkappa*t0)
	* PEXP( -x*x/(4.0*CRkappa*t0));
    } // end i for
  } // end if

  return SUCCESS;
}
