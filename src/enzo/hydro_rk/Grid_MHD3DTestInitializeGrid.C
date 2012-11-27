/***********************************************************************
/
/  GRID CLASS (INITIALIZE MHD 3D TEST)
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
#include "EOS.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int grid::MHD3DTestInitializeGrid(int MHD3DProblemType,
				  float rhol, float rhou,
				  float vxl,  float vxu,
				  float vyl,  float vyu,
				  float pl,   float pu,
				  float Bxl,  float Bxu,
				  float Byl,  float Byu)
{  

  // Neutral
  const float HIIFraction = 1.2e-5;
  const float HeIIFraction = 1.0e-14;
  const float HeIIIFraction = 1.0e-17;
  const float HMFraction = 2.0e-9;
  const float H2IFraction = 2.0e-20;
  const float H2IIFraction = 3.0e-14;

  // Ionized
  const float HIIFractionIon = 0.99 * CoolData.HydrogenFractionByMass;
  const float HeIIFractionIon = 1.0e-14;
  const float HeIIIFractionIon = 0.99 * (1.0-CoolData.HydrogenFractionByMass);

  float HIFraction, HeIFraction, eFraction;
  float HIFractionIon, HeIFractionIon, eFractionIon;
  HIFraction = CoolData.HydrogenFractionByMass - HIIFraction;
  HIFractionIon = CoolData.HydrogenFractionByMass - HIIFractionIon;
  if (MultiSpecies > 1)
    HIFraction -= HMFraction + H2IFraction + H2IIFraction;
  HeIFraction = 1.0 - CoolData.HydrogenFractionByMass - 
    HeIIFraction - HeIIIFraction;
  HeIFractionIon = 1.0 - CoolData.HydrogenFractionByMass - 
    HeIIFractionIon - HeIIIFractionIon;
  eFraction = HIIFraction + 0.25*HeIIFraction + 0.5*HeIIIFraction;
  eFractionIon = HIIFractionIon + 0.25*HeIIFractionIon + 0.5*HeIIIFractionIon;
  if (MultiSpecies > 1)
    eFraction += 0.5*H2IIFraction - HMFraction;

  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
    VelocityUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, 
	   &VelocityUnits, Time);
  const float IonizedThreshold = 1e4 / TemperatureUnits / (Gamma-1.0);

  /* create fields */

  int GENum, TENum;

  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  FieldType[NumberOfBaryonFields++] = Velocity2;
  FieldType[NumberOfBaryonFields++] = Velocity3;
  FieldType[TENum = NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism) {
    FieldType[GENum = NumberOfBaryonFields++] = InternalEnergy;
  }

  if (HydroMethod == MHD_RK) {
    FieldType[NumberOfBaryonFields++] = Bfield1;
    FieldType[NumberOfBaryonFields++] = Bfield2;
    FieldType[NumberOfBaryonFields++] = Bfield3;
    FieldType[NumberOfBaryonFields++] = PhiField;
  }

  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, 
    H2INum, H2IINum, DINum, DIINum, HDINum;

  if (MultiSpecies) {
    FieldType[DeNum    = NumberOfBaryonFields++] = ElectronDensity;
    FieldType[HINum    = NumberOfBaryonFields++] = HIDensity;
    FieldType[HIINum   = NumberOfBaryonFields++] = HIIDensity;
    FieldType[HeINum   = NumberOfBaryonFields++] = HeIDensity;
    FieldType[HeIINum  = NumberOfBaryonFields++] = HeIIDensity;
    FieldType[HeIIINum = NumberOfBaryonFields++] = HeIIIDensity;
    if (MultiSpecies > 1) {
      FieldType[HMNum    = NumberOfBaryonFields++] = HMDensity;
      FieldType[H2INum   = NumberOfBaryonFields++] = H2IDensity;
      FieldType[H2IINum  = NumberOfBaryonFields++] = H2IIDensity;
    }
    if (MultiSpecies > 2) {
      FieldType[DINum   = NumberOfBaryonFields++] = DIDensity;
      FieldType[DIINum  = NumberOfBaryonFields++] = DIIDensity;
      FieldType[HDINum  = NumberOfBaryonFields++] = HDIDensity;
    }
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
  
  int field;
  for (field = 0; field < NumberOfBaryonFields; field++) {
    if (BaryonField[field] == NULL) {
      BaryonField[field] = new float[size];
    }
  }

  /* transform pressure to total energy */
  float etotl, etotu, v2, B2;
  v2 = vxl * vxl + vyl * vyl;
  B2 = Bxl * Bxl + Byl * Byl;
  etotl = pl / ((Gamma-1.0)*rhol) + 0.5*v2 + 0.5*B2/rhol;

  v2 = vxu * vxu + vyu * vyu;
  B2 = Bxu * Bxu + Byu * Byu;
  etotu = pu / ((Gamma-1.0)*rhou) + 0.5*v2 + 0.5*B2/rhou;

  
  int igrid;
  FLOAT x, y, z;
  if (MHD3DProblemType == 0) { // Planar shock
  float pres, eintl, eintu, h, cs, dpdrho, dpde;
  for (int k = 0; k < GridDimension[2]; k++) {
    for (int j = 0; j < GridDimension[1]; j++) {
      for (int i = 0; i < GridDimension[0]; i++) {
	/* Compute position */
	igrid = i + j*GridDimension[0] + k*GridDimension[0]*GridDimension[1];
      
	x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	if (x <= 0.1) {
	  if (k==2) {
	    printf("rhol=%"GSYM", vxl=%"GSYM", pl=%"GSYM"\n", rhol, vxl, pl);
	  }

	  BaryonField[iden ][igrid] = rhol;
	  BaryonField[ivx  ][igrid] = vxl;
	  BaryonField[ivy  ][igrid] = vyl;
	  BaryonField[ivz  ][igrid] = 0.0;
	  BaryonField[ietot][igrid] = etotl;
	  if (DualEnergyFormalism) {
	    BaryonField[ieint][igrid] = pl / ((Gamma-1.0)*rhol);
	  }
	  BaryonField[iBx  ][igrid] = Bxl;
	  BaryonField[iBy  ][igrid] = Byl;
	  BaryonField[iBz  ][igrid] = 0.0;
	  BaryonField[iPhi ][igrid] = 0.0;
	} else {
	  if (k==2) {
	    printf("rhor=%"GSYM", vxr=%"GSYM", pr=%"GSYM"\n", rhou, vxu, pu);
	  }
	  BaryonField[iden ][igrid] = rhou;
	  BaryonField[ivx  ][igrid] = vxu;
	  BaryonField[ivy  ][igrid] = vyu;
	  BaryonField[ivz  ][igrid] = 0.0;
	  BaryonField[ietot][igrid] = etotu;
	  if (DualEnergyFormalism) {
	    BaryonField[ieint][igrid] = pu / ((Gamma-1.0)*rhou);
	  }
	  BaryonField[iBx  ][igrid] = Bxu;
	  BaryonField[iBy  ][igrid] = Byu;
	  BaryonField[iBz  ][igrid] = 0.0;
	  BaryonField[iPhi ][igrid] = 0.0;
	}
      }
    }  
  }
  }

  if (MHD3DProblemType == 1) { // Uniform Density with a Shear
    float pres, eintl, eintu, h, cs, dpdrho, dpde;
    for (int k = 0; k < GridDimension[2]; k++) {
      for (int j = 0; j < GridDimension[1]; j++) {
	for (int i = 0; i < GridDimension[0]; i++) {
	  /* Compute position */
	  igrid = i + j*GridDimension[0] + k*GridDimension[0]*GridDimension[1];
	  
	  x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	  y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	  z = CellLeftEdge[2][j] + 0.5*CellWidth[2][j];
	  
	  float rho;
	  rho = 1.;
	  pres = rho/Gamma; // sound speed = 1
	  EOS(pres, rho, eintl, h, cs, dpdrho, dpde, 0, 1);
	  // impose mode perturbation
	  vxl = 0.5*cos(2.0*M_PI*y);
	  vyl = 0.;
	  etotl = eintl + 0.5*(vxl*vxl + vyl*vyl) + 0.5*(Bxl*Bxl+Byl*Byl)/rho;
	  BaryonField[iden ][igrid] = rho;
	  BaryonField[ivx  ][igrid] = vxl;
	  BaryonField[ivy  ][igrid] = vyl;
	  BaryonField[ivz  ][igrid] = 0.0;
	  
	  BaryonField[ietot][igrid] = etotl;
	  if (DualEnergyFormalism) {
	    BaryonField[ieint][igrid] = pl / ((Gamma-1.0)*rho);
	  }
	  if (HydroMethod == MHD_RK) {
	    BaryonField[iBx  ][igrid] = Bxl;
	    BaryonField[iBy  ][igrid] = Byl;
	    BaryonField[iBz  ][igrid] = 0.0;
	    BaryonField[iPhi ][igrid] = 0.0;
	  }
	  
	}
      }  
    }
  }
  
  /* Rayleigh-Taylor problem with a single-mode (2) or multiple modes
     (3) */

  int i, j, k, index, seed;
  float pres, rho, ramp, dpdrho, dpde, h, cs, vz, eintl, eintu;
  FLOAT DomainWidth[MAX_DIMENSION];
  //const float delz = 5e-3;  // range in z to apply ramp
  const float delz = 5e-10;  // range in z to apply ramp
  const float amplitude = 0.01; // perturbation amplitude

  if (MHD3DProblemType == 2 || MHD3DProblemType == 3) {
    if (HydroMethod == MHD_RK)
      ENZO_FAIL("Rayleigh-Taylor problem in 3D not setup for MHD yet.");

    seed = 123456789;
    srand(seed);

    for (i = 0; i < MAX_DIMENSION; i++)
      DomainWidth[i] = DomainRightEdge[i] - DomainLeftEdge[i];

    for (k = 0; k < GridDimension[2]; k++) {
      z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];

      // Calculate pressure from hydrostatic equilibrium
      ramp = 1.0 / (1.0 + exp(-2.0*z / delz));
      rho = rhol + ramp * (rhou-rhol);
      pres = pl + ConstantAcceleration[2] * rho * z;
      if (z <= 0)
	EOS(pres, rho, eintl, h, cs, dpdrho, dpde, 0, 1);
      else
	EOS(pres, rho, eintu, h, cs, dpdrho, dpde, 0, 1);

      for (j = 0; j < GridDimension[1]; j++) {
	y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	index = GRIDINDEX_NOGHOST(0,j,k);
	for (i = 0; i < GridDimension[0]; i++, index++) {
	  x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];

	  // Generate perturbation (2==single, 3==multiple)
	  if (MHD3DProblemType == 2)
	    vz = 0.125*amplitude * 
	      (1.0 + cos(2.0*M_PI*x / DomainWidth[0])) *
	      (1.0 + cos(2.0*M_PI*y / DomainWidth[1])) *
	      (1.0 + cos(2.0*M_PI*z / DomainWidth[2]));
	  else
	    vz = amplitude * 
	      ((float) (rand()) / (float) (RAND_MAX) - 0.5) *
	      (1.0 + cos(2.0*M_PI*z / DomainWidth[2]));

	  // Lower domain
	  if (z <= 0.0) {

	    etotl = eintl + 0.5*(vxl*vxl + vyl*vyl + vz*vz);
	    BaryonField[iden][index] = rhol;
	    BaryonField[ivx][index] = vxl;
	    BaryonField[ivy][index] = vyl;
	    BaryonField[ivz][index] = vz;
	    BaryonField[ietot][index] = etotl;

	    if (DualEnergyFormalism)
	      BaryonField[ieint][index] = pl / ((Gamma-1.0)*rho);

	  } // ENDIF (lower)

	  // Upper domain
	  else {

	    etotu = eintu + 0.5*(vxu*vxu + vyu*vyu + vz*vz);
	    BaryonField[iden][index] = rhou;
	    BaryonField[ivx][index] = vxu;
	    BaryonField[ivy][index] = vyu;
	    BaryonField[ivz][index] = vz;
	    BaryonField[ietot][index] = etotu;

	    if (DualEnergyFormalism)
	      BaryonField[ieint][index] = pu / ((Gamma-1.0)*rho);

	  } // ENDELSE (upper)
	  
	} // ENDFOR i
      } // ENDFOR j
    } // ENDFOR k
  } // ENDIF type == 2||3


  /* Magnetic explosion */
  /* Tom Abel May 2011 */


  if (MHD3DProblemType == 4) { 

    float pres, eintl, eintu, h, cs, dpdrho, dpde,rhot, bx ,by, bz;


    float *ax = new float[size];
    float *ay = new float[size];
    float *az = new float[size];
    for (int k = 0; k < GridDimension[2]; k++) {
      for (int j = 0; j < GridDimension[1]; j++) {
	for (int i = 0; i < GridDimension[0]; i++) {
	/* Compute position */
	  igrid =  GRIDINDEX_NOGHOST(i,j,k);
	  x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i] - 0.5;
	  y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j] - 0.5;
	  z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k] - 0.5;
	  float r2 = x*x+y*y+z*z;
	  float R2 = 1./16. * 1./16.;
	  float  f;
	  f = exp(-r2/2/R2);
	  ax[igrid] = - y*f * CellWidth[0][i];;  //Vector potential
	  ay[igrid] =   x*f * CellWidth[1][j];;
	  az[igrid] =   0.  * CellWidth[2][k];;
	}
      }
    }


    for (int k = 1; k < GridDimension[2]-1; k++) {
      for (int j = 1; j < GridDimension[1]-1; j++) {
	for (int i = 1; i < GridDimension[0]-1; i++) {
	  float rho;
	  rho = rhol;
	  float vx=0.,vy=0.,vz=0.;
	  int  igridyp1, igridym1, igridzp1, igridzm1;

	  igrid    = GRIDINDEX_NOGHOST(i,j,k);
	  igridyp1 = GRIDINDEX_NOGHOST(i,j+1,k);
	  igridym1 = GRIDINDEX_NOGHOST(i,j-1,k);
	  igridzp1 = GRIDINDEX_NOGHOST(i,j,k+1);
	  igridzm1 = GRIDINDEX_NOGHOST(i,j,k-1);
	  // B = Curl A              gives divergence free B-field from vector potential
	  bx = Bxl * ((az[igridyp1]-az[igridym1])/2/CellWidth[1][j] -
		      (ay[igridzp1]-ay[igridzm1])/2/CellWidth[2][k]);
	  by = Bxl * ((ax[igridzp1]-ax[igridzm1])/2/CellWidth[2][k] -
		      (az[igrid +1]-az[igrid -1])/2/CellWidth[0][i]);
	  bz = Bxl * ((ay[igrid +1]-ay[igrid -1])/2/CellWidth[0][i] -
		      (ax[igridyp1]-ax[igridym1])/2/CellWidth[1][j]);

	  pres = (EOSType > 0) ? EOSSoundSpeed*EOSSoundSpeed*rho : // isothermal sound speed 
	    EOSSoundSpeed*EOSSoundSpeed*rho ; 
	  EOS(pres, rho, eintl, h, cs, dpdrho, dpde, 0, 1);

	  etotl = eintl + 0.5*(vx*vx + vy*vy + vz*vz) + 0.5*(bx*bx+by*by+bz*bz)/rho;
	  BaryonField[iden ][igrid] = rho;
	  BaryonField[ivx  ][igrid] = 0.0 ;
	  BaryonField[ivy  ][igrid] = 0.0;
	  BaryonField[ivz  ][igrid] = 0.0;
	  
	  BaryonField[ietot][igrid] = etotl;
	  if (DualEnergyFormalism) {
	    BaryonField[ieint][igrid] = pres / ((Gamma-1.0)*rho);
	  }
	  if (HydroMethod == MHD_RK) {
	    BaryonField[iBx  ][igrid] = bx;
	    BaryonField[iBy  ][igrid] = by;
	    BaryonField[iBz  ][igrid] = bz;
	    BaryonField[iPhi ][igrid] = 0.0;
	  }
	} // endfor i
      } // endfor j 
    } // endfor k 

    delete [] ax;
    delete [] ay;
    delete [] az;

  } // if MHD3DProblemType == 4
  
  
  /* Set uniform species fractions */

  if (MultiSpecies > 0) {
    for (k = 0, index = 0; k < GridDimension[2]; k++)
      for (j = 0; j < GridDimension[1]; j++)
	for (i = 0; i < GridDimension[0]; i++, index++) {

	  // Assume the velocity perturbations are small and that
	  // thermal energy dominates
	  if (BaryonField[TENum][index] > IonizedThreshold) {
	    BaryonField[DeNum][index] = eFractionIon * BaryonField[0][index];
	    BaryonField[HINum][index] = HIFractionIon* BaryonField[0][index];
	    BaryonField[HIINum][index] = HIIFractionIon * BaryonField[0][index];
	    BaryonField[HeINum][index] = HeIFractionIon * BaryonField[0][index];
	    BaryonField[HeIINum][index] = HeIIFractionIon * BaryonField[0][index];
	    BaryonField[HeIIINum][index] = HeIIIFractionIon * BaryonField[0][index];
	  } else {
	    BaryonField[DeNum][index] = eFraction * BaryonField[0][index];
	    BaryonField[HINum][index] = HIFraction * BaryonField[0][index];
	    BaryonField[HIINum][index] = HIIFraction * BaryonField[0][index];
	    BaryonField[HeINum][index] = HeIFraction * BaryonField[0][index];
	    BaryonField[HeIINum][index] = HeIIFraction * BaryonField[0][index];
	    BaryonField[HeIIINum][index] = HeIIIFraction * BaryonField[0][index];
	  }
	}
  }

  if (MultiSpecies > 1) {
    for (k = 0, index = 0; k < GridDimension[2]; k++)
      for (j = 0; j < GridDimension[1]; j++)
	for (i = 0; i < GridDimension[0]; i++, index++) {
	  BaryonField[HMNum][index] = HMFraction * BaryonField[0][index];
	  BaryonField[H2INum][index] = H2IFraction * BaryonField[0][index];
	  BaryonField[H2IINum][index] = H2IIFraction * BaryonField[0][index];
	}
  }

  if (MultiSpecies > 2) {
    for (k = 0, index = 0; k < GridDimension[2]; k++)
      for (j = 0; j < GridDimension[1]; j++)
	for (i = 0; i < GridDimension[0]; i++, index++) {
	  BaryonField[DINum][index] = CoolData.DeuteriumToHydrogenRatio * 
	    BaryonField[HINum][index];
	  BaryonField[DIINum][index] = CoolData.DeuteriumToHydrogenRatio * 
	    BaryonField[HIINum][index];
	  BaryonField[HDINum][index] = CoolData.DeuteriumToHydrogenRatio * 
	    BaryonField[H2INum][index];
	}
  }

  return SUCCESS;
}
