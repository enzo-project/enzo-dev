/***********************************************************************
/
/  GRID CLASS (INITIALIZE MHD 2D TEST)
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

int grid::MHD2DTestInitializeGrid(int MHD2DProblemType,
				  int UseColour,
				  float RampWidth,
				  float rhol, float rhou,
				  float vxl,  float vxu,
				  float vyl,  float vyu,
				  float pl,   float pu,
				  float Bxl,  float Bxu,
				  float Byl,  float Byu)
{  

  int iZ;

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
  if (UseDivergenceCleaning) {
    FieldType[NumberOfBaryonFields++] = Phi_pField;
    //FieldType[NumberOfBaryonFields++] = DebugField;  
  }
  if (UseColour) {
    FieldType[iZ = NumberOfBaryonFields++] = Metallicity;
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
      for (int i = 0; i < size; i++) 
	BaryonField[field][i] = 0.0;
    }
  }
  */
  AllocateGrids();

  /* transform pressure to total energy */
  float etotl, etotu, v2, B2;
  v2 = vxl * vxl + vyl * vyl;
  B2 = Bxl * Bxl + Byl * Byl;
  etotl = pl / ((Gamma-1.0)*rhol) + 0.5*v2 + 0.5*B2/rhol;

  v2 = vxu * vxu + vyu * vyu;
  B2 = Bxu * Bxu + Byu * Byu;
  etotu = pu / ((Gamma-1.0)*rhou) + 0.5*v2 + 0.5*B2/rhou;

  if (MHD2DProblemType == 0 && !UseConstantAcceleration) {
    printf("Rayleigh-Taylor problem must have UseConstantAcceleration = 1\n");
    return FAIL;
  }


  int igrid, igrid2, field;
  FLOAT x, y;

  /* MHD2DProblemType == 0: Rayleigh-Taylor problem */

  if (MHD2DProblemType == 0) { 

    float pres, eintl, eintu, h, cs, dpdrho, dpde;
    float delx=0.00005; // range in y over which to apply the ramp
    for (int j = 0; j < GridDimension[1]; j++) {
      for (int i = 0; i < GridDimension[0]; i++) {
	/* Compute position */
	igrid = i + j*GridDimension[0];
	
	x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	float ramp,rho;
	ramp =  1./((1.+exp(-2/delx*(y-0.))));
	rho = rhol + ramp*(rhou-rhol);
	if (y <= 0.) {
	  // Rayleigh-Taylor problem, calculate pressure from hydro equilibrium
	  float g = ConstantAcceleration[1];
	  pres = pl+g*rhol*(y);
	  EOS(pres, rho, eintl, h, cs, dpdrho, dpde, 0, 1);
	  // impose mode perturbation
	  vyl = .01 * (1.0+cos(4.0*M_PI*(x))) * (1.0+cos(2.0/1.5*M_PI*(y))) * 0.25;
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
    if ( UseMHDCT ){
      field=0;
	    igrid2 = i+MagneticDims[field][0]*j;
      MagneticField[field][igrid2] = Bxl;
      field=1;
	    igrid2 = i+MagneticDims[field][0]*j;
      MagneticField[field][igrid2] = Byl;
      field=2;
	    igrid2 = i+MagneticDims[field][0]*j;
      MagneticField[field][igrid2] = 0;
    }
	  if (UseColour)
	    BaryonField[iZ][igrid] = 1.0;
	} else {  // endif y<0
	  /* calculate pressure from hydro equilibrium */
	  float g = ConstantAcceleration[1];
	  pres = pl+g*rho*(y);
	  EOS(pres, rho, eintu, h, cs, dpdrho, dpde, 0, 1);
	  /* impose mode perturbation */
	  vyu = 0.01 * (1.0+cos(4.0*M_PI*(x))) * (1.0+cos(2.0/1.5*M_PI*(y))) * 0.25;
	  etotu = eintu + 0.5*(vxu*vxu + vyu*vyu) + 0.5*(Bxu*Bxu+Byu*Byu)/rho;
	  BaryonField[iden ][igrid] = rho;
	  BaryonField[ivx  ][igrid] = vxu;
	  BaryonField[ivy  ][igrid] = vyu;
	  BaryonField[ivz  ][igrid] = 0.0;
	  BaryonField[ietot][igrid] = etotu;
	  if (DualEnergyFormalism) {
	    BaryonField[ieint][igrid] = pu / ((Gamma-1.0)*rho);
	  }
	  if (HydroMethod == MHD_RK) {
	    BaryonField[iBx  ][igrid] = Bxu;
	    BaryonField[iBy  ][igrid] = Byu;
	    BaryonField[iBz  ][igrid] = 0.0;
	    BaryonField[iPhi ][igrid] = 0.0;
	  }
    if ( UseMHDCT ){
      field=0;
	    igrid2 = i+MagneticDims[field][0]*j;
      MagneticField[field][igrid2] = Bxu;
      field=1;
	    igrid2 = i+MagneticDims[field][0]*j;
      MagneticField[field][igrid2] = Byu;
      field=2;
	    igrid2 = i+MagneticDims[field][0]*j;
      MagneticField[field][igrid2] = 0;
    }
	  if (UseColour)
	    BaryonField[iZ][igrid] = tiny_number;
	}
      } // endfor i
    } // endfor j 
  } // if MHD2DProblemType == 0 
  
  /* MHD2DProblemType = 1: MHD rotor 
   Reference: G. Toth, J. Comput. Phys. 161 (2000) 605. */

  if (MHD2DProblemType == 1) { 

    FLOAT r;
    FLOAT r0 = 0.1, r1 = 0.115;
    float eint, h, cs, dpdrho, dpde, etot;
//  float rho0 = 10.0, rho1 = 1.0, pres = 1.0, v0 = 2.0, Bx0 = 5.;
    float rho0 = 10.0, rho1 = 1.0, pres = 1.0, v0 = 2.0, Bx0 = 5./sqrt(4.*M_PI);

    for (int j = 0; j < GridDimension[1]; j++) {
      for (int i = 0; i < GridDimension[0]; i++) {

	igrid = i + j*GridDimension[0];

	x = CellLeftEdge[0][i] + 0.5 * CellWidth[0][i];
	y = CellLeftEdge[1][j] + 0.5 * CellWidth[1][j];	
	r = sqrt(pow(x - 0.5, 2) + pow(y - 0.5, 2));

	if (r < r0) {

	  BaryonField[iden][igrid] = rho0;
	  BaryonField[ivx ][igrid] = -v0 * (y - 0.5) / r0;
	  BaryonField[ivy ][igrid] =  v0 * (x - 0.5) / r0;
	  BaryonField[ivz ][igrid] = 0.0;
	  EOS(pres, rho0, eint, h, cs, dpdrho, dpde, 0, 1);
	  etot = eint + 
	    0.5 * (pow(BaryonField[ivx][igrid], 2) + pow(BaryonField[ivy][igrid], 2)) +
	    0.5 * Bx0 * Bx0 / rho0;
	  BaryonField[ietot][igrid] = etot;
	  if (DualEnergyFormalism) {
	    BaryonField[ieint][igrid] = eint;
	  }
	  BaryonField[iBx ][igrid] = Bx0;
	  BaryonField[iBy ][igrid] = 0.0;
	  BaryonField[iBz ][igrid] = 0.0;
	  BaryonField[iPhi][igrid] = 0.0;
	  if (UseColour)
	    BaryonField[iZ][igrid] = 1.0;

	}
	else if (r < r1) {

	  /* Toth's "Taper" function */
	  FLOAT f = (r1-r)/(r1-r0);

	  BaryonField[iden][igrid] = 1.0 + 9.0 * f;
	  BaryonField[ivx ][igrid] = -f * v0 * (y - 0.5) / r;
          BaryonField[ivy ][igrid] =  f * v0 * (x - 0.5) / r;
          BaryonField[ivz ][igrid] = 0.0;
	  EOS(pres, rho0, eint, h, cs, dpdrho, dpde, 0, 1);
          etot = eint + 
	    0.5 * (pow(BaryonField[ivx][igrid], 2) + pow(BaryonField[ivy][igrid], 2)) +
	    0.5 * Bx0 * Bx0 / (1.0 + 9.0 * f);
          BaryonField[ietot][igrid] = etot;
          if (DualEnergyFormalism) {
            BaryonField[ieint][igrid] = eint;
          }
	  BaryonField[iBx ][igrid] = Bx0;
          BaryonField[iBy ][igrid] = 0.0;
          BaryonField[iBz ][igrid] = 0.0;
	  BaryonField[iPhi][igrid] = 0.0;
	  if (UseColour)
	    BaryonField[iZ][igrid] = tiny_number;

        }
	else {

	  BaryonField[iden][igrid] = rho1;
          BaryonField[ivx ][igrid] = 0.0;
          BaryonField[ivy ][igrid] = 0.0;
          BaryonField[ivz ][igrid] = 0.0;
          EOS(pres, rho0, eint, h, cs, dpdrho, dpde, 0, 1);
          etot = eint + 0.5*Bx0*Bx0/rho1;
          BaryonField[ietot][igrid] = etot;
          if (DualEnergyFormalism) {
            BaryonField[ieint][igrid] = eint;
          }
          BaryonField[iBx ][igrid] = Bx0;
          BaryonField[iBy ][igrid] = 0.0;
          BaryonField[iBz ][igrid] = 0.0;
	  BaryonField[iPhi][igrid] = 0.0;
	  if (UseColour)
	    BaryonField[iZ][igrid] = tiny_number;

	}
      }
    }
  }
	
  /* MHD2DProblemType = 2: MHD blast wave
   Reference: T. A. Gardiner & J. M. Stone, J. Comput. Phys. 205 (2005) 509. */

  if (MHD2DProblemType == 2) { 
    FLOAT r;
    float eint, h, cs, dpdrho, dpde, etot;
    /*
    FLOAT r0 = 0.1;
    float rho0 = 1.0, pres0 = 10.0, pres1 = 1.0, Bx0 = 1.0/sqrt(2), By0 = 1.0/sqrt(2);
    */
    FLOAT r0 = 0.125;
    float rho0 = 1.0, pres0 = 100.0, pres1 = 1.0, Bx0 = 10.0/sqrt(2), By0 = 10.0/sqrt(2);

    if (HydroMethod != MHD_RK) {
      Bx0 = By0 = 0.;
    }
    
    for (int j = 0; j < GridDimension[1]; j++) {
      for (int i = 0; i < GridDimension[0]; i++) {
	
	igrid = i + j*GridDimension[0];
	
	x = CellLeftEdge[0][i] + 0.5 * CellWidth[0][i] ;
	y = CellLeftEdge[1][j] + 0.5 * CellWidth[1][j] ;	
	r = sqrt(pow(x , 2) + pow(y , 2));

	BaryonField[iden][igrid] = rho0;
	BaryonField[ivx ][igrid] = 0.0;
	BaryonField[ivy ][igrid] = 0.0;
	BaryonField[ivz ][igrid] = 0.0;
	if (HydroMethod == MHD_RK) {
	  BaryonField[iBx ][igrid] = Bx0;
	  BaryonField[iBy ][igrid] = By0;
	  BaryonField[iBz ][igrid] = 0.0;
	  BaryonField[iPhi][igrid] = 0.0;
	}
	if (r < r0) {
          EOS(pres0, rho0, eint, h, cs, dpdrho, dpde, 0, 1);
	  etot = eint + 0.5 * (Bx0 * Bx0 + By0 * By0) / rho0;
	  BaryonField[ietot][igrid] = etot;
	  if (DualEnergyFormalism) 
	    BaryonField[ieint][igrid] = eint;
	  if (UseColour)
	    BaryonField[iZ][igrid] = 1.0;

	} else {
	  EOS(pres1, rho0, eint, h, cs, dpdrho, dpde, 0, 1);
	  etot = eint + 0.5 * (Bx0 * Bx0 + By0 * By0) / rho0;
	  BaryonField[ietot][igrid] = etot;
	  if (DualEnergyFormalism) 
	    BaryonField[ieint][igrid] = eint;
	  if (UseColour)
	    BaryonField[iZ][igrid] = tiny_number;
	
	}

      }
    }
  }

  /* MHD2DProblemType = 3: MHD Kelvin-Helmholtz instability.
   Reference: T. A. Gardiner & J. M. Stone, J. Comput. Phys. 205 (2005) 509. */

  if (MHD2DProblemType == 3) { 

    float eint, h, cs, dpdrho, dpde, etot;
    float rho0 = 2.0, pres = 2.5, rho1 = 1.0, vx0 = 0.5, vx1 = -0.5, Bx = 0.5;

    srand(1564);
    
    float eps1, eps2;
    for (int j = 0; j < GridDimension[1]; j++) {
      for (int i = 0; i < GridDimension[0]; i++) {
	
	igrid = i + j*GridDimension[0];
	
	y = CellLeftEdge[1][j] + 0.5 * CellWidth[1][j];	

	eps1 = rand();
	eps1 /= RAND_MAX;
	eps1 *= 0.01;
	eps2 = rand();
	eps2 /= RAND_MAX; 
	eps2 *= 0.01;

	if (y <= 0.75 && y >= 0.25) {
	  BaryonField[iden][igrid] = rho0;
	  BaryonField[ivx ][igrid] = vx0 + eps1;
	  BaryonField[ivy ][igrid] = eps2;
	  BaryonField[ivz ][igrid] = 0.0;
	  BaryonField[iBx ][igrid] = Bx;
	  BaryonField[iBy ][igrid] = 0.0;
	  BaryonField[iBz ][igrid] = 0.0;
	  BaryonField[iPhi][igrid] = 0.0;
          EOS(pres, rho0, eint, h, cs, dpdrho, dpde, 0, 1);
	  etot = eint + 0.5 * vx0 * vx0 + 0.5 * (Bx * Bx) / rho0;
	  BaryonField[ietot][igrid] = etot;
	  if (DualEnergyFormalism) 
	    BaryonField[ieint][igrid] = eint;
	  if (UseColour)
	    BaryonField[iZ][igrid] = 1.0;
	} else {
	  BaryonField[iden][igrid] = rho1;
	  BaryonField[ivx ][igrid] = vx1 + eps1;
	  BaryonField[ivy ][igrid] = eps2;
	  BaryonField[ivz ][igrid] = 0.0;
	  BaryonField[iBx ][igrid] = Bx;
	  BaryonField[iBy ][igrid] = 0.0;
	  BaryonField[iBz ][igrid] = 0.0;
	  BaryonField[iPhi][igrid] = 0.0;
          EOS(pres, rho1, eint, h, cs, dpdrho, dpde, 0, 1);
	  etot = eint + 0.5 * vx1 * vx1 + 0.5 * (Bx * Bx) / rho1;
	  BaryonField[ietot][igrid] = etot;
	  if (DualEnergyFormalism) 
	    BaryonField[ieint][igrid] = eint;
	  if (UseColour)
	    BaryonField[iZ][igrid] = tiny_number;
	}

      }
    }
  }

  /* MHD2DProblemType = 4: MHD Kelvin-Helmholtz instability for the 
     cold/warm ISM interface. */

  if (MHD2DProblemType == 4) { 

    float DensityUnits = 1.67e-24;
    float LengthUnits = 3.0856e19;
    float TimeUnits = 1.0/sqrt(6.673e-8*DensityUnits);
    float VelocityUnits = LengthUnits/TimeUnits;

    float eint, h, cs, dpdrho, dpde, etot;
    
    float rho0 = 1.67e-22/DensityUnits,
      rho1 = 1.67e-24/DensityUnits,
      cs0 = 1e5/VelocityUnits,
      cs1 = 1e6/VelocityUnits,
      vx0 = 0,
      vx1 = 1e6/VelocityUnits,
      Bx = 0.0;
    float pres0 = rho0*cs0*cs0,
      pres1 = rho1*cs1*cs1;

    srand(1564);
    
    float eps1, eps2;
    for (int j = 0; j < GridDimension[1]; j++) {
      for (int i = 0; i < GridDimension[0]; i++) {
	
	igrid = i + j*GridDimension[0];
	
	y = CellLeftEdge[1][j] + 0.5 * CellWidth[1][j];	

	eps1 = rand();
	eps1 /= RAND_MAX;
	eps1 *= 0.01;
	eps2 = rand();
	eps2 /= RAND_MAX; 
	eps2 *= 0.01;

	if (y >= 0.5) {
	  BaryonField[iden][igrid] = rho0;
	  BaryonField[ivx ][igrid] = vx0 + eps1;
	  BaryonField[ivy ][igrid] = eps2;
	  BaryonField[ivz ][igrid] = 0.0;
	  BaryonField[iBx ][igrid] = Bx;
	  BaryonField[iBy ][igrid] = 0.0;
	  BaryonField[iBz ][igrid] = 0.0;
	  BaryonField[iPhi][igrid] = 0.0;
          EOS(pres0, rho0, eint, h, cs, dpdrho, dpde, 0, 1);
	  etot = eint + 0.5 * vx0 * vx0 + 0.5 * (Bx * Bx) / rho0;
	  BaryonField[ietot][igrid] = etot;
	  if (DualEnergyFormalism) 
	    BaryonField[ieint][igrid] = eint;
	} else {
	  BaryonField[iden][igrid] = rho1;
	  BaryonField[ivx ][igrid] = vx1 + eps1;
	  BaryonField[ivy ][igrid] = eps2;
	  BaryonField[ivz ][igrid] = 0.0;
	  BaryonField[iBx ][igrid] = Bx;
	  BaryonField[iBy ][igrid] = 0.0;
	  BaryonField[iBz ][igrid] = 0.0;
	  BaryonField[iPhi][igrid] = 0.0;
          EOS(pres1, rho1, eint, h, cs, dpdrho, dpde, 0, 1);
	  etot = eint + 0.5 * vx1 * vx1 + 0.5 * (Bx * Bx) / rho1;
	  BaryonField[ietot][igrid] = etot;
	  if (DualEnergyFormalism) 
	    BaryonField[ieint][igrid] = eint;
	}

      }
    }
  }

  /* MHD2DProblemType = 5: 
   *   Shock-vortex interaction. 
   *   Reference: Rault, Chiavassa & Donat, 2003, J. Scientific Computing, 19, 1.
   */

  if (MHD2DProblemType == 5) { 

    float rho0 = 1.0, cs = 1.0, eint, etot;
    float p0 = rho0*cs*cs;
    float eint0 = p0/((Gamma-1.0)*rho0);
    float Bx = 0.0;

    FLOAT b = 0.175, a = 0.075;
    //FLOAT b = 0, a = 0;
    
    /* parameters for the composite vortex tube */
    float Mv = 1.7;
    float a2 = a*a, b2 = b*b;
    float c2 = rho0*pow(b, 2.0*Mv*Mv*a2*b2/pow(a2-b2,2));
    float c1 = c2*exp(-0.5*Mv*Mv + Mv*Mv*a2/pow(a2-b2,2)*(0.5*a2 - 2.0*b2*log(a) - 0.5*b2*b2/a2));
    //float c1 = rho0*exp(-0.5*Mv*Mv);

    /* parameters for the shock */
    float Ms = 2;
    float rho1 = rho0/((Gamma-1.0)/(Gamma+1.0) + 2.0/((Gamma+1.0)*Ms*Ms));
    float v1 = cs*sqrt((2.0+(Gamma-1.0)*Ms*Ms)/(2.0*Gamma*Ms*Ms-Gamma+1));
    float v0 = Ms*cs;

    FLOAT x, y, r, xc = 0.2, yc = 0.5, xs = 0.4;
    int igrid;
    for (int j = 0; j < GridDimension[1]; j++) {
      for (int i = 0; i < GridDimension[0]; i++) {
	
	igrid = i + j*GridDimension[0];
	
	x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];	

	r = sqrt(pow(x-xc,2) + pow(y-yc,2));

	float sintheta = (y-yc)/r;
	float costheta = (x-xc)/r;

	if (r > b && x < xs) {
	  BaryonField[iden][igrid] = rho0;
	  BaryonField[ivx ][igrid] = v0;
	  BaryonField[ivy ][igrid] = 0.0;
	  BaryonField[ivz ][igrid] = 0.0;
	  etot = eint0 + 0.5*v0*v0;
	  BaryonField[ietot][igrid] = etot;
	  if (DualEnergyFormalism) 
	    BaryonField[ieint][igrid] = eint;
	  if (HydroMethod == MHD_RK) 
	    BaryonField[ietot][igrid] += 0.5 * (Bx * Bx) / rho0;
	} 

	if (r > b && x > xs) {
	  BaryonField[iden][igrid] = rho1;
          BaryonField[ivx ][igrid] = v1;
          BaryonField[ivy ][igrid] = 0.0;
          BaryonField[ivz ][igrid] = 0.0;
          etot = eint0 + 0.5*v1*v1;
          BaryonField[ietot][igrid] = etot;
          if (DualEnergyFormalism)
            BaryonField[ieint][igrid] = eint;
          if (HydroMethod == MHD_RK)
            BaryonField[ietot][igrid] += 0.5 * (Bx * Bx) / rho0;
        }


	if (r <= b && r >= a) {
	  BaryonField[iden][igrid] = 
	    c2*exp(Mv*Mv*a2/pow(a2-b2,2)*(0.5*r*r-2.0*b2*log(r)-0.5*b2*b2/(r*r)));
	  float vth = Mv*cs*a/(a2-b2)*(r-b2/r);
	  BaryonField[ivx][igrid] = -vth*sintheta + v0;
	  BaryonField[ivy][igrid] = vth*costheta;
	  BaryonField[ivz][igrid] = 0.0;
	  etot = eint0 + 0.5*pow(BaryonField[ivx][igrid],2)+0.5*pow(BaryonField[ivy][igrid],2);
	  BaryonField[ietot][igrid] = etot;
          if (DualEnergyFormalism)
            BaryonField[ieint][igrid] = eint;
          if (HydroMethod == MHD_RK)
            BaryonField[ietot][igrid] += 0.5*(Bx*Bx)/BaryonField[iden][igrid] ;	  
	} 
	if (r < a) {
	  BaryonField[iden][igrid] = c1*exp(0.5*Mv*Mv*r*r/a2);
	  float vth = Mv*cs*r/a;
	  BaryonField[ivx][igrid] = -vth*sintheta + v0;
          BaryonField[ivy][igrid] = vth*costheta;
          BaryonField[ivz][igrid] = 0.0;
          etot = eint0 + 0.5*pow(BaryonField[ivx][igrid],2)+0.5*pow(BaryonField[ivy][igrid],2);
          BaryonField[ietot][igrid] = etot;
          if (DualEnergyFormalism)
            BaryonField[ieint][igrid] = eint;
          if (HydroMethod == MHD_RK)
            BaryonField[ietot][igrid] += 0.5*(Bx*Bx)/BaryonField[iden][igrid] ;
	}

	/*if (x < 0.05) {
	  BaryonField[iden][igrid] = rho1;
	  BaryonField[ivx ][igrid] = v1;
          BaryonField[ivy ][igrid] = 0.0;
          BaryonField[ivz ][igrid] = 0.0;
	  BaryonField[ietot][igrid] += 0.5*v1*v1;
	  if (HydroMethod == MHD_RK) {
	    BaryonField[ietot][igrid] -= 0.5*(Bx*Bx)/rho0;
	    BaryonField[ietot][igrid] += 0.5*(Bx*Bx)/rho1;
	  }
	  }*/

	if (HydroMethod == MHD_RK) {
	  BaryonField[iBx ][igrid] = Bx;
	  BaryonField[iBy ][igrid] = 0.0;
	  BaryonField[iBz ][igrid] = 0.0;
	  BaryonField[iPhi][igrid] = 0.0;
	}

      }
    }
  }

  /* MHD2DProblemType 6: Sedov-Taylor Blast Wave 
   * Reference: Fryxell et al, 2000, ApJS, 131, 273
   * (I don't think this could be used in 2D test - comments by Ji-hoon Kim, Apr.2010)
   */

  if (MHD2DProblemType == 6) { 

    float E0 = 1.0;
    const FLOAT r0 = 0.06;
    float P0 = 3.0*(Gamma-1.0)*E0/(4.0*M_PI*pow(r0,3));
    float rho0 = 1.0;
    float P1 = 1e-5;
    float eint, h, cs, dpdrho, dpde;

    for (int j = 0; j < GridDimension[1]; j++) {
      FLOAT y = CellLeftEdge[1][j] + 0.5 * CellWidth[1][j] -0.5;	
      for (int i = 0; i < GridDimension[0]; i++) {
	
	int igrid = i + j*GridDimension[0];
	FLOAT x = CellLeftEdge[0][i] + 0.5 * CellWidth[0][i] -0.5;
	FLOAT r = sqrt(x*x + y*y);

	if (r <= r0) {
	  BaryonField[iden][igrid] = rho0;
	  BaryonField[ivx ][igrid] = 0.0;
	  BaryonField[ivy ][igrid] = 0.0;
	  BaryonField[ivz ][igrid] = 0.0;
          EOS(P0, rho0, eint, h, cs, dpdrho, dpde, 0, 1);
	  BaryonField[ietot][igrid] = eint;
	  if (DualEnergyFormalism) 
	    BaryonField[ieint][igrid] = eint;
	} else {
	  BaryonField[iden][igrid] = rho0;
	  BaryonField[ivx ][igrid] = 0.0;
	  BaryonField[ivy ][igrid] = 0.0;
	  BaryonField[ivz ][igrid] = 0.0;
          EOS(P1, rho0, eint, h, cs, dpdrho, dpde, 0, 1);
	  BaryonField[ietot][igrid] = eint;
	  if (DualEnergyFormalism) 
	    BaryonField[ieint][igrid] = eint;
	}

	if (HydroMethod == MHD_RK) {
	  BaryonField[iBx ][igrid] = Bxl;
	  BaryonField[iBy ][igrid] = Byl;
	  BaryonField[iBz ][igrid] = 0.0;
	  BaryonField[iPhi][igrid] = 0.0;
	  BaryonField[ietot][igrid] += 0.5*(Bxl*Bxl+Byl*Byl)/BaryonField[iden][igrid];
	}

      }
    }
  }

  /* MHD2DProblemType 7: Cylindrical Sedov-Taylor Blast Wave 
   * Reference: Fryxell et al, 2000, ApJS, 131, 273
   * (I don't think this could be used in 2D test - comments by
   *  Ji-hoon Kim, Apr.2010)
   * (Some parameters fixed - comments by Ji-hoon Kim, Jul.2010)
   */

  if (MHD2DProblemType == 7) { 

    float E0 = 1.0;
    const FLOAT r0 = 0.02;
    //    const FLOAT r0 = 0.1;
    float P0 = 3.0*(Gamma-1.0)*E0/(3.0*M_PI*pow(r0,2));
    float rho0 = 1.0;
    float P1 = 1e-5;
    float eint, h, cs, dpdrho, dpde;

    for (int j = 0; j < GridDimension[1]; j++) {
      FLOAT y = CellLeftEdge[1][j] + 0.5 * CellWidth[1][j] -0.5;	
      for (int i = 0; i < GridDimension[0]; i++) {
	
	int igrid = i + j*GridDimension[0];
	FLOAT x = CellLeftEdge[0][i] + 0.5 * CellWidth[0][i] -0.5;
	FLOAT r = sqrt(x*x + y*y);

	if (r <= r0) {
	  BaryonField[iden][igrid] = rho0;
	  BaryonField[ivx ][igrid] = 0.0;
	  BaryonField[ivy ][igrid] = 0.0;
	  BaryonField[ivz ][igrid] = 0.0;
          EOS(P0, rho0, eint, h, cs, dpdrho, dpde, 0, 1);
	  BaryonField[ietot][igrid] = eint;
	  if (DualEnergyFormalism) 
	    BaryonField[ieint][igrid] = eint;
	} else {
	  BaryonField[iden][igrid] = rho0;
	  BaryonField[ivx ][igrid] = 0.0;
	  BaryonField[ivy ][igrid] = 0.0;
	  BaryonField[ivz ][igrid] = 0.0;
          EOS(P1, rho0, eint, h, cs, dpdrho, dpde, 0, 1);
	  BaryonField[ietot][igrid] = eint;
	  if (DualEnergyFormalism) 
	    BaryonField[ieint][igrid] = eint;
	}

	if (HydroMethod == MHD_RK) {
	  BaryonField[iBx ][igrid] = Bxl;
	  BaryonField[iBy ][igrid] = Byl;
	  BaryonField[iBz ][igrid] = 0.0;
	  BaryonField[iPhi][igrid] = 0.0;
	}

      }
    }
  }


  /* Standing shock like in MHD2SProblemtype 5 but with small density perturbation
     downstream from the shock to look at odd-even coupling. */ 

  if (MHD2DProblemType == 8) { 

    float rho0 = 1.0, cs = 1.0, eint0,eint1, etot;
    float p0 = rho0*cs*cs;
    float Bx = Bxl;

    /* parameters for the shock */
    float Ms = 2;
    float rho1 = rho0/(1-2/(Gamma+1)*(1-1/Ms/Ms));
    float v1 = cs*sqrt((2.0+(Gamma-1.0)*Ms*Ms)/(2.0*Gamma*Ms*Ms-Gamma+1));
    float v0 = Ms*cs;
    eint0 = p0/((Gamma-1.0)*rho0);
    eint1 = eint1; // isothermal so same specific energy

    FLOAT x, y, r, xc = 0.2,  xs = 0.6, delx=0.001;
    int igrid;

    for (int j = 0; j < GridDimension[1]; j++) {
      for (int i = 0; i < GridDimension[0]; i++) {
	
	igrid = i + j*GridDimension[0];
	
	x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];	


	if (HydroMethod == MHD_RK) {
	  BaryonField[iBx ][igrid] = Bx;
	  BaryonField[iBy ][igrid] = 0.0;
	  BaryonField[iBz ][igrid] = 0.0;
	  BaryonField[iPhi][igrid] = 0.0;
	}

 	float ramp =  1./((1.+exp(-2/delx*(x-xs))));
	BaryonField[iden][igrid] = rho0 + ramp*(rho1-rho0);
	BaryonField[ivx ][igrid] = v0 + ramp*(v1-v0);
	BaryonField[ivy ][igrid] = 0.0;
	BaryonField[ivz ][igrid] = 0.0;
	etot = eint0 +  0.5*pow(v0 + ramp*(v1-v0),2);
	BaryonField[ietot][igrid] = etot;
	if (DualEnergyFormalism) 
	  BaryonField[ieint][igrid] = eint0+ramp*(eint1-eint0);
	if (HydroMethod == MHD_RK) 
	  BaryonField[ietot][igrid] += 0.5 * (Bx * Bx) / BaryonField[iden][igrid];


	// perturb ahead of shock 
	if ((y > .45 ) && (y < .55) && (x < 0.9*xs) && (x>0.8*xs)) 
	  BaryonField[iden][igrid] += 1e-4*rho0;

      }
    }
  }

  /* Kelvin Helmholtz Problem as in Robertson, Kravstov, Gnedin & Abel */ 
  if (MHD2DProblemType == 9) {
    FLOAT x,y;
    int n=0;
    float omega=vyu; // reusing UpperVelocityY as perturbation amplitude

    float ramp,eint;
    float delx=0.05; // range in y over which to apply the ramp

    for (int j = 0; j < GridDimension[1]; j++) {
      for (int i = 0; i < GridDimension[0]; i++, n++) {
	
	igrid = i + j*GridDimension[0];
	
	/* Compute position */
	
	x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	
	ramp =  1./((1.+exp(-2/delx*(y-0.25)))*(1.+exp(-2/delx*(0.75-y))));
	BaryonField[iden ][igrid] = rhou + ramp*(rhol-rhou);
	BaryonField[ivx  ][igrid] = vxu  + ramp*(vxl-vxu);
	BaryonField[ivy  ][igrid] = 0.0;
	BaryonField[ivz  ][igrid] = 0.0;
	
	// y-velocity perturbation: 
	BaryonField[ivy][igrid]  += omega*sin(2.*M_PI*x);

	eint = (pu + ramp*(pl-pu)) / ((Gamma-1.0)*BaryonField[iden][igrid]);
	BaryonField[ietot][igrid] = eint + 
	  0.5*(pow(BaryonField[ivx  ][igrid],2)+
	       pow(BaryonField[ivy  ][igrid],2)+
	       pow(BaryonField[ivz  ][igrid],2));
	if (DualEnergyFormalism) {
	  BaryonField[ieint][igrid] = eint ;
	}
	if (HydroMethod == MHD_RK) {
	  BaryonField[iBx  ][igrid] = Bxu + ramp*(Bxl-Bxu);
	  BaryonField[iBy  ][igrid] = Byu + ramp*(Byl-Byu);
	  BaryonField[iBz  ][igrid] = 0.0;
	  BaryonField[iPhi ][igrid] = 0.0;
	}
    if ( UseMHDCT ){
      field=0;
	    igrid2 = i+MagneticDims[field][0]*j;
      MagneticField[field][igrid2] = Bxu + ramp*(Bxl-Bxu);
      field=1;
	    igrid2 = i+MagneticDims[field][0]*j;
      MagneticField[field][igrid2] = Byu + ramp*(Byl-Byu);
      field=2;
	    igrid2 = i+MagneticDims[field][0]*j;
      MagneticField[field][igrid2] = 0;
    }
      }
    }
  } // endif MHD2DProblemType == 9 

  /* MHD2DProblemType == 10: Modified Rayleigh-Taylor problem */
  /* Domain goes from 0 0.1 to 0.5 0.9 */
  if (MHD2DProblemType == 10) { 

    float pres, eintl, eintu, h, cs, dpdrho, dpde;
    float delx=0.05; // range in y over which to apply the ramp
    for (int j = 0; j < GridDimension[1]; j++) {
      for (int i = 0; i < GridDimension[0]; i++) {
	/* Compute position */
	igrid = i + j*GridDimension[0];
	
	x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	float ramp,rho;
	ramp =  1./((1.+exp(-2/delx*(y-0.5))));
	rho = rhol + ramp*(rhou-rhol);
	if (y <= 0.5) {
	  // Rayleigh-Taylor problem, calculate pressure from hydro equilibrium
	  float g = ConstantAcceleration[1];
	  pres = pl+g*rho*(y-0.5);
	  EOS(pres, rho, eintl, h, cs, dpdrho, dpde, 0, 1);
	  // impose mode perturbation
	  vyl = 0.7 * (1.0+cos(4.0*M_PI*(x+.25))) * (1.0+cos(2.0/.8*M_PI*(y-0.5))) * 0.25;
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
	} else {  // endif y<0
	  /* calculate pressure from hydro equilibrium */
	  float g = ConstantAcceleration[1];
	  pres = pl+g*rho*(y-0.5);
	  EOS(pres, rho, eintu, h, cs, dpdrho, dpde, 0, 1);
	  /* impose mode perturbation */
	  vyu = 0.7 * (1.0+cos(4.0*M_PI*(x+.25))) * (1.0+cos(2.0/.8*M_PI*(y-0.5))) * 0.25;
	  etotu = eintu + 0.5*(vxu*vxu + vyu*vyu) + 0.5*(Bxu*Bxu+Byu*Byu)/rho;
	  BaryonField[iden ][igrid] = rho;
	  BaryonField[ivx  ][igrid] = vxu;
	  BaryonField[ivy  ][igrid] = vyu;
	  BaryonField[ivz  ][igrid] = 0.0;
	  BaryonField[ietot][igrid] = etotu;
	  if (DualEnergyFormalism) {
	    BaryonField[ieint][igrid] = pu / ((Gamma-1.0)*rho);
	  }
	  if (HydroMethod == MHD_RK) {
	    BaryonField[iBx  ][igrid] = Bxu;
	    BaryonField[iBy  ][igrid] = Byu;
	    BaryonField[iBz  ][igrid] = 0.0;
	    BaryonField[iPhi ][igrid] = 0.0;
	  }
	}
      } // endfor i
    } // endfor j 
  } // if MHD2DProblemType == 10

  /* MHD2DProblemType == 11: Uniform density with a sinusoidal shear velocity */
  /* Domain goes from 0 0 to 1 1 */
  // Absolutely nothing happens inthis test case. Compare to SPH ... 
  if (MHD2DProblemType == 11) { 

    float pres, eintl, eintu, h, cs, dpdrho, dpde;

    for (int j = 0; j < GridDimension[1]; j++) {
      for (int i = 0; i < GridDimension[0]; i++) {
	/* Compute position */
	igrid = i + j*GridDimension[0];
	
	x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
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
      } // endfor i
    } // endfor j 
  } // if MHD2DProblemType == 11

  /* MHD2DProblemType == 12: Uniform density with two spherical central
   overpressured regions launching sound waves. Periodic boundaries. */
  /* Domain goes from 0 0 to 1 1 */
  if (MHD2DProblemType == 12) { 
    float pres, eintl, eintu, h, cs, dpdrho, dpde;
    for (int j = 0; j < GridDimension[1]; j++) {
      for (int i = 0; i < GridDimension[0]; i++) {
	/* Compute position */
	igrid = i + j*GridDimension[0];
	
	x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];

	float rho,ramp,rad,r0=0.05, P0;
	rad = sqrt(pow(y-0.5,2)+pow(x-0.3,2));
	ramp = 1./((1. + exp((-2*( rad - r0)/RampWidth))));
	rad = sqrt(pow(y-0.5,2)+pow(x-0.7,2));
	ramp *= 1./((1. + exp((-2*( rad - r0)/RampWidth))));
	rho = 1.;
	// reused pu below as a fractional perturbation keeping sound speed =1 
	// everywhere but the central perturbations
	P0 = rho/Gamma ; // sound speed = 1
	pres = pu*P0 + ramp*(P0 - pu*P0);
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
      } // endfor i
    } // endfor j 
  } // if MHD2DProblemType == 12


  /* MHD2DProblemType == 13: A Blob Test (not finished yet) */
  /* Spherical Cloud gets overrun by Mach 10 shock */
  /* Domain goes from 0 0 to 2 1 */
  if (MHD2DProblemType == 13) { 

    float pres, eintl, eintu, h, cs, dpdrho, dpde,ramp,rhot;

    for (int j = 0; j < GridDimension[1]; j++) {
      for (int i = 0; i < GridDimension[0]; i++) {
	/* Compute position */
	igrid = i + j*GridDimension[0];
	
	x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	float rho,rad,r0=0.075, vx, vy;
	rad = sqrt(pow(y-0.5,2)+pow(x-0.5,2));
	ramp = 1./((1. + exp((-2*( rad - r0)/RampWidth))));
	rho = rhol+ramp*(rhou-rhol);
	pres = rhol/Gamma; // sound speed = 1 in medium
	EOS(pres, rho, eintl, h, cs, dpdrho, dpde, 0, 1);
	// impose mode perturbation
	vx= vxl + ramp*(vxu-vxl);
	vy = vyl*(x-0.5);
	etotl = eintl + 0.5*(vx*vx + vy*vy) + 0.5*(Bxl*Bxl+Byl*Byl)/rho;
	BaryonField[iden ][igrid] = rho;
	BaryonField[ivx  ][igrid] = vx;
	BaryonField[ivy  ][igrid] = vy;
	BaryonField[ivz  ][igrid] = 0.0;
	
	BaryonField[ietot][igrid] = etotl;
	if (DualEnergyFormalism) {
	  BaryonField[ieint][igrid] = pres / ((Gamma-1.0)*rho);
	}
	if (HydroMethod == MHD_RK) {
	  BaryonField[iBx  ][igrid] = Bxl;
	  BaryonField[iBy  ][igrid] = Byl;
	  BaryonField[iBz  ][igrid] = 0.0;
	  BaryonField[iPhi ][igrid] = 0.0;
	}
      } // endfor i
    } // endfor j 
  } // if MHD2DProblemType == 13

  /* Wengen 2 test */
  /* Tom Abel Sept 2010 */

  if (MHD2DProblemType == 14) { 

    float pres, eintl, eintu, h, cs, dpdrho, dpde,ramp,rhot, bx ,by;

    for (int j = 0; j < GridDimension[1]; j++) {
      for (int i = 0; i < GridDimension[0]; i++) {
	/* Compute position */
	igrid = i + j*GridDimension[0];
	x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	ramp =  1./(1.+exp(-2/RampWidth*(y-0.5)));
	float rho, vx, vy, f;
	rho = rhol + ramp*(rhou-rhol);

	pres = (EOSType > 0) ? EOSSoundSpeed*EOSSoundSpeed*rho : // isothermal sound speed = 0.112611
	  EOSSoundSpeed*EOSSoundSpeed*rho ; 
	EOS(pres, rho, eintl, h, cs, dpdrho, dpde, 0, 1);
	// impose mode perturbation
	//	f = cos(2.*M_PI*x*10.)*exp(-fabs(y-0.5)*10.);
	f = cos(2.*M_PI*x*10.)*exp(-fabs(y-0.5)*10.)*cos(2.*M_PI*x*3);
	vx = f * (vxl+ ramp*(vxu-vxl))  ;
	vy = vyl + ramp*(vyu - vyl);
	bx = (Bxl+ ramp*(Bxu-Bxl))  ;
	by = (Byl+ ramp*(Byu-Byl))  ;
	etotl = eintl + 0.5*(vx*vx + vy*vy) + 0.5*(bx*bx+by*by)/rho;
	BaryonField[iden ][igrid] = rho;
	BaryonField[ivx  ][igrid] = vx ;
	BaryonField[ivy  ][igrid] = vy;
	BaryonField[ivz  ][igrid] = 0.0;
	
	BaryonField[ietot][igrid] = etotl;
	if (DualEnergyFormalism) {
	  BaryonField[ieint][igrid] = pres / ((Gamma-1.0)*rho);
	}
	if (HydroMethod == MHD_RK) {
	  BaryonField[iBx  ][igrid] = bx;
	  BaryonField[iBy  ][igrid] = by;
	  BaryonField[iBz  ][igrid] = 0.0;
	  BaryonField[iPhi ][igrid] = 0.0;
	}
      } // endfor i
    } // endfor j 
  } // if MHD2DProblemType == 14

  /* Tom Abel June 2011 */
  /* Simple Sinusoidal B field */


  if (MHD2DProblemType == 15) { 

    float pres, eintl, eintu, h, cs, dpdrho, dpde,ramp,rhot, bx ,by;

    for (int j = 0; j < GridDimension[1]; j++) {
      for (int i = 0; i < GridDimension[0]; i++) {
	/* Compute position */
	igrid = i + j*GridDimension[0];
	x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	ramp =  1. ;
	float rho, vx, vy, f;
	rho = rhol;

	pres = (EOSType > 0) ? EOSSoundSpeed*EOSSoundSpeed*rho : // isothermal sound speed = 0.112611
	  rho/Gamma ; // for ideal equation of state set sounds speed = 1 
	EOS(pres, rho, eintl, h, cs, dpdrho, dpde, 0, 1);
	// impose mode perturbation
	f = (cos(2.*M_PI*x) + 1.);

	vx = 0.;
	vy = 0.;
	bx = Bxl * f ;
	by = Byl * f ;

	etotl = eintl + 0.5*(vx*vx + vy*vy) + 0.5*(bx*bx+by*by)/rho;

	BaryonField[iden ][igrid] = rho;
	BaryonField[ivx  ][igrid] = vx ;
	BaryonField[ivy  ][igrid] = vy;
	BaryonField[ivz  ][igrid] = 0.0;
	
	BaryonField[ietot][igrid] = etotl;
	if (DualEnergyFormalism) {
	  BaryonField[ieint][igrid] = pres / ((Gamma-1.0)*rho);
	}
	if (HydroMethod == MHD_RK) {
	  BaryonField[iBx  ][igrid] = bx;
	  BaryonField[iBy  ][igrid] = by;
	  BaryonField[iBz  ][igrid] = 0.0;
	  BaryonField[iPhi ][igrid] = 0.0;
	}
      } // endfor i
    } // endfor j 
  } // if MHD2DProblemType == 15
  if( UseMHDCT ){
      CenterMagneticField();
  }

  return SUCCESS;
}
