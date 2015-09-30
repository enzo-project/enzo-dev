/***********************************************************************
/
/  GRID CLASS (UPDATE MHD VARIABLES)
/
/  written by: Peng Wang
/  date:       June, 2007
/  modified1:
/
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
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
#include "EOS.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int grid::UpdateMHDPrim(float **dU, float c1, float c2)
{

  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, 
    B1Num, B2Num, B3Num, PhiNum;

  this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				   Vel3Num, TENum, B1Num, B2Num, B3Num, PhiNum);

  int i, j, k, n, n_dU, dim, igrid, field, size, activesize;
  for (dim = 0, size = 1; dim < GridRank; dim++) {
    size *= GridDimension[dim];
  }

  for (dim = 0, activesize = 1; dim < GridRank; dim++) {
    activesize *= (GridDimension[dim] - 2*NumberOfGhostZones);
  }

  float *D, *sum;
  float SmallX = 1e-20;

  // ORIGINAL
  if ( (NSpecies+NColor) > 0) {
    D = new float[activesize];
    sum = new float[activesize];
    for (i = 0; i < activesize; i++) {
      D[i] = 0.0;
      sum[i] = 0.0;
    }
  }

  float *Prim[NEQ_MHD+NSpecies+NColor];
  float *OldPrim[NEQ_MHD+NSpecies+NColor];
  this->ReturnHydroRKPointers(Prim, false);   
  this->ReturnOldHydroRKPointers(OldPrim, false);

  //##### Want to mix species and colors for renormalization?  Normally you don't
  int MixSpeciesAndColors = 0;
  int NSpecies_renorm;
  if (MixSpeciesAndColors) 
    NSpecies_renorm = NSpecies+NColor;
  else if (NoMultiSpeciesButColors) {
    NSpecies_renorm = NSpecies;
  }
  else
    switch (MultiSpecies) { //update pure species! not colours!
    case 0:  NSpecies_renorm = 0;  break;
    case 1:  NSpecies_renorm = 5;  break;
    case 2:  NSpecies_renorm = 8;  break;
    case 3:  NSpecies_renorm = 11; break;
    default: NSpecies_renorm = 0;  break;
    }
  
  // update species
  for (field = NEQ_MHD; field < NEQ_MHD+NSpecies_renorm; field++) {  
    n = 0;
    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	igrid = (k * GridDimension[1] + j) * GridDimension[0] + GridStartIndex[0];
        for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++, igrid++) {
          Prim[field][igrid] = c1*OldPrim[field][igrid] +
            (1-c1)*Prim[field][igrid]*Prim[iden][igrid] + c2*dU[field][n];
          if (NoMultiSpeciesButColors != 1)
            D[n] += Prim[field][igrid];
        }
      }
    }
  }

  // renormalize species 
  if (NoMultiSpeciesButColors != 1) {
    for (field = NEQ_MHD; field < NEQ_MHD+NSpecies_renorm; field++) {
      n = 0;
      for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
        for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
          igrid = (k * GridDimension[1] + j) * GridDimension[0] + GridStartIndex[0];
          for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++, igrid++) {
            Prim[field][igrid] = min(1.0, max((Prim[field][igrid]/D[n]), SmallX));
            sum[n] += Prim[field][igrid];
          }
        }
      }
    }
    
    for (field = NEQ_MHD; field < NEQ_MHD+NSpecies_renorm; field++) {
      n = 0;
      for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
        for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
          igrid = (k * GridDimension[1] + j) * GridDimension[0] + GridStartIndex[0];
          for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++, igrid++)
            Prim[field][igrid] /= sum[n];
        }
      }
    }
  } //close if (NoMultiSpeciesButColors == 1)
  /* Update conserved variables */

  float rho_old, vx_old, vy_old, vz_old, e_old, etot_old, Tau_old, eint_old,
    rho, vx, vy, vz, e, etot, Tau, eint, p, v2,
    D_new, S1_new, S2_new, S3_new, Tau_new, h, cs, dpdrho, dpde, Eint_new,
    Bx_old, By_old, Bz_old, Bx, By, Bz, Bx_new, By_new, Bz_new,
    Phi_old, Phi, Phi_new, B2;

  float rhou, lenu, tempu, tu, velu;
  GetUnits(&rhou, &lenu, &tempu, &tu, &velu, Time);

  n = 0;
  FLOAT x, y, z, r;

 
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {
	// first convert to conserved variables to do the update
	igrid = (k * GridDimension[1] + j) * GridDimension[0] + i;
	r = sqrt(x*x + y*y + z*z);

	rho_old  = OldBaryonField[DensNum][igrid];
	vx_old   = OldBaryonField[Vel1Num][igrid];
	vy_old   = OldBaryonField[Vel2Num][igrid];
	vz_old   = OldBaryonField[Vel3Num][igrid];
	etot_old = OldBaryonField[TENum][igrid];
	Tau_old = rho_old*etot_old;
	if (DualEnergyFormalism) {
	  eint_old = OldBaryonField[GENum][igrid];
	}
	Bx_old   = OldBaryonField[B1Num][igrid];
	By_old   = OldBaryonField[B2Num][igrid];
	Bz_old   = OldBaryonField[B3Num][igrid];
	Phi_old  = OldBaryonField[PhiNum][igrid];

	rho  = BaryonField[DensNum][igrid];
	vx   = BaryonField[Vel1Num][igrid];
	vy   = BaryonField[Vel2Num][igrid];
	vz   = BaryonField[Vel3Num][igrid];
	etot = BaryonField[TENum][igrid];
	Tau  = rho*etot;
	if (DualEnergyFormalism) {
	  eint = BaryonField[ieint][igrid];
	}
	Bx   = BaryonField[B1Num][igrid];
	By   = BaryonField[B2Num][igrid];
	Bz   = BaryonField[B3Num][igrid];
	Phi  = BaryonField[PhiNum ][igrid];

	D_new   = c1*rho_old + (1.0-c1)*rho + c2*dU[iD][n];
	S1_new  = c1*rho_old*vx_old + (1.0-c1)*rho*vx + c2*dU[iS1][n];
	S2_new  = c1*rho_old*vy_old + (1.0-c1)*rho*vy + c2*dU[iS2][n];
	S3_new  = c1*rho_old*vz_old + (1.0-c1)*rho*vz + c2*dU[iS3][n];
	Tau_new = c1*Tau_old + (1.0-c1)*Tau + c2*dU[iEtot][n];
	Bx_new  = c1*Bx_old + (1.0-c1)*Bx + c2*dU[iBx][n];
	By_new  = c1*By_old + (1.0-c1)*By + c2*dU[iBy][n];
	Bz_new  = c1*Bz_old + (1.0-c1)*Bz + c2*dU[iBz][n];
	Phi_new = c1*Phi_old + (1.0-c1)*Phi + c2*dU[iPhi][n];
	if (DualEnergyFormalism) {
	  Eint_new = c1*rho_old*eint_old + (1.0-c1)*rho*eint + c2*dU[iEint][n];
	}

	if (D_new < 0 || isnan(D_new)) {

	  printf("UpdateMHDPrim: rho <0 at %"ISYM" %"ISYM" %"ISYM": rho_old=%"FSYM", rho=%"FSYM", rho_new=%"FSYM", dU[iD]=%"FSYM"\n", 
		 i, j, k, rho_old, rho, D_new, dU[iD][n]);
	  D_new = max(rho, SmallRho);
	  printf("UpdateMHDPrim: use rho: %"FSYM"\n", D_new);
	  //	  D_new = rho;

	  return FAIL;
	}

	// convert back to primitives
	vx = S1_new/D_new;
	vy = S2_new/D_new;
	vz = S3_new/D_new;
	etot = Tau_new/D_new;

	
	if (etot < 0 && EOSType == 0) {
	  float v2_old = vx_old*vx_old + vy_old*vy_old + vz_old*vz_old;
	  float B2_old = Bx_old*vx_old + By_old*By_old + Bz_old*Bz_old;
	  printf("UpdateMHDPrim: tau < 0. etot_old=%"GSYM", etot=%"GSYM", etot_new=%"GSYM", v2=%"GSYM", v2old=%"GSYM", dU[iTau] = %"GSYM", dtFixed = %"GSYM"\n", 
		 Tau_old/rho_old, Tau/rho, Tau_new/D_new, v2, v2_old, dU[iEtot][n]*CellWidth[0][0]/dtFixed, dtFixed);
	  printf("rho_new=%"GSYM", rho=%"GSYM", rho_old=%"GSYM", B2_old/rho_old=%"GSYM"\n", D_new, rho, rho_old, B2_old/rho_old);
	  //return FAIL;
	}


	// if using polytropic EOS, calculate etot directly from density
	if (EOSType > 0) { 
	  EOS(p, D_new, eint, h, cs, dpdrho, dpde, EOSType, 0);
	  v2 = vx*vx + vy*vy + vz*vz;
	  B2 = Bx_new*Bx_new + By_new*By_new + Bz_new*Bz_new;
	  etot = eint + 0.5*v2 + 0.5*B2/D_new;
	}
	

	BaryonField[DensNum][igrid] = D_new;

	BaryonField[Vel1Num][igrid] = vx;
	BaryonField[Vel2Num][igrid] = vy;
	BaryonField[Vel3Num][igrid] = vz;
	BaryonField[TENum][igrid] = etot;
	BaryonField[B1Num][igrid] = Bx_new;
	BaryonField[B2Num][igrid] = By_new;
	BaryonField[B3Num][igrid] = Bz_new;
	BaryonField[PhiNum][igrid] = Phi_new*exp(-c1*dtFixed*pow(C_h/C_p,2));


	if (DualEnergyFormalism) {
	  v2 = vx*vx + vy*vy + vz*vz;
	  B2 = Bx_new*Bx_new + By_new*By_new + Bz_new*Bz_new;
	  float emin = SmallT/(Mu*(Gamma-1.0));

	  eint = Eint_new/D_new;
	  float eint1 = etot - 0.5*v2 - 0.5*B2/D_new;
	  if (eint1 > 0) {
	    EOS(p, D_new, eint1, h, cs, dpdrho, dpde, EOSType, 2);
	  }
	  else {
	    cs = 0.0;
	  }
	  if (cs*cs > DualEnergyFormalismEta1*v2 &&  
	      cs*cs > DualEnergyFormalismEta1*B2/D_new &&
	      eint1 > 0.5*eint) {
	    eint = eint1;
	  }
	  eint = max(eint, emin);
	  BaryonField[GENum][igrid] = eint;
	  BaryonField[TENum][igrid] = eint + 0.5*v2 + 0.5*B2/D_new;
	  if (BaryonField[GENum][igrid] < 0.0) {
	    printf("UpdateMHDPrim: eint < 0, cs2=%"GSYM", eta*v2=%"GSYM", eint=%"GSYM", etot=%"GSYM", 0.5*v2=%"GSYM", p=%"GSYM", rho=%"GSYM",0.5*B2/rho=%"GSYM"\n", 
		   cs*cs, DualEnergyFormalismEta1*v2, eint, etot, 0.5*v2, p, D_new, 0.5*B2/rho);
	    printf("dU[%"ISYM"]=%"GSYM", dU[ieint]=%"GSYM", eint_old=%"GSYM",eint1=%"GSYM"\n", iEtot, dU[iEtot][n], dU[iEint][n], eint_old, eint1);
	    return FAIL;
	  }
	}
      }
    }
  }

  // Convert species from mass fraction to density  
  if (NoMultiSpeciesButColors != 1)
    for (field = NEQ_MHD; field < NEQ_MHD+NSpecies+NColor; field++)  
      for (n = 0; n < size; n++) 
        Prim[field][n] *= Prim[iden][n];

  this->UpdateElectronDensity();

  if ( (NSpecies+NColor) > 0) {
    delete [] D;
    delete [] sum;
  }
  
  return SUCCESS;
}
