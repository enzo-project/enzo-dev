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

#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "TopGridData.h"
#include "Grid.h"
#include "EOS.h"
#include <math.h>

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

int grid::UpdateMHDPrim(float **dU, float c1, float c2)
{

  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }

  int size = 1;
  for (int dim = 0; dim < GridRank; dim++) {
    size *= GridDimension[dim];
  }

  int activesize = 1;
  for (int dim = 0; dim < GridRank; dim++) {
    activesize *= (GridDimension[dim] - 2*DEFAULT_GHOST_ZONES);
  }

  float *D, *sum;
  float SmallX = 1e-20;

  if ( (NSpecies+NColor) > 0) {
    D = new float[activesize];
    sum = new float[activesize];
    for (int i = 0; i < activesize; i++) {
      D[i] = 0.0;
      sum[i] = 0.0;
    }
  }

  /* Update species */
  int igrid;
  for (int field = NEQ_MHD; field < NEQ_MHD+NSpecies; field++) {
    int n = 0;
    for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
        for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {
          igrid = (k * GridDimension[1] + j) * GridDimension[0] + i;
          BaryonField[field][igrid] = c1*OldBaryonField[field][igrid] +
            (1-c1)*BaryonField[field][igrid]*BaryonField[iden][igrid] + c2*dU[field][n];
          D[n] += BaryonField[field][igrid];
        }
      }
    }
  }

  /* Renormalize species */

  for (int field = NEQ_MHD; field < NEQ_MHD+NSpecies; field++) {
    int n = 0;
    for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
        for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {
          igrid = (k * GridDimension[1] + j) * GridDimension[0] + i;
          BaryonField[field][igrid] = min(1.0, max((BaryonField[field][igrid]/D[n]), SmallX));	  
	  BaryonField[field][igrid] = BaryonField[field][igrid]/D[n];
          sum[n] += BaryonField[field][igrid];
        }
      }
    }
  }

  for (int field = NEQ_MHD; field < NEQ_MHD+NSpecies; field++) {
    int n = 0;
    for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
        for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {
          igrid = (k * GridDimension[1] + j) * GridDimension[0] + i;
          BaryonField[field][igrid] /= sum[n];
        }
      }
    }
  }


  /* Update conserved variables */

  float rho_old, vx_old, vy_old, vz_old, e_old, etot_old, Tau_old, eint_old,
    rho, vx, vy, vz, e, etot, Tau, eint, p, v2,
    D_new, S1_new, S2_new, S3_new, Tau_new, h, cs, dpdrho, dpde, Eint_new,
    Bx_old, By_old, Bz_old, Bx, By, Bz, Bx_new, By_new, Bz_new,
    Phi_old, Phi, Phi_new, B2;

  float rhou, lenu, tempu, tu, velu;
  GetUnits(&rhou, &lenu, &tempu, &tu, &velu, Time);

  int n = 0;
  FLOAT x, y, z, r;
  for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {
	// first convert to conserved variables to do the update
	igrid = (k * GridDimension[1] + j) * GridDimension[0] + i;
	r = sqrt(x*x + y*y + z*z);

	rho_old  = OldBaryonField[iden ][igrid];
	vx_old   = OldBaryonField[ivx  ][igrid];
	vy_old   = OldBaryonField[ivy  ][igrid];
	vz_old   = OldBaryonField[ivz  ][igrid];
	etot_old = OldBaryonField[ietot][igrid];
	Tau_old = rho_old*etot_old;
	if (DualEnergyFormalism) {
	  eint_old = OldBaryonField[ieint][igrid];
	}
	Bx_old   = OldBaryonField[iBx  ][igrid];
	By_old   = OldBaryonField[iBy  ][igrid];
	Bz_old   = OldBaryonField[iBz  ][igrid];
	Phi_old  = OldBaryonField[iPhi ][igrid];

	rho  = BaryonField[iden ][igrid];
	vx   = BaryonField[ivx  ][igrid];
	vy   = BaryonField[ivy  ][igrid];
	vz   = BaryonField[ivz  ][igrid];
	etot = BaryonField[ietot][igrid];
	Tau  = rho*etot;
	if (DualEnergyFormalism) {
	  eint = BaryonField[ieint][igrid];
	}
	Bx   = BaryonField[iBx  ][igrid];
	By   = BaryonField[iBy  ][igrid];
	Bz   = BaryonField[iBz  ][igrid];
	Phi  = BaryonField[iPhi ][igrid];


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
	  printf("UpdatePrim: rho <0 at %d %d %d: rho_old=%g, rho=%g, rho_new=%g, dU[iD]=%g\n", 
		 i, j, k, rho_old, rho, D_new, dU[iD][n]);
	  //D_new = max(D_new, SmallRho);
	  D_new = rho;
	  //D_new = rho;
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
	  printf("UpdatePrim: tau < 0. etot_old=%g, etot=%g, etot_new=%g, v2=%g, v2old=%g, dU[iTau] = %g\n", 
		 Tau_old/rho_old, Tau/rho, Tau_new/D_new, v2, v2_old, dU[iEtot][n]*CellWidth[0][0]/dtFixed);
	  printf("rho_new=%g, rho=%g, rho_old=%g, B2_old/rho_old=%g\n", D_new, rho, rho_old, B2_old/rho_old);
	  //return FAIL;
	}


	// if using polytropic EOS, calculate etot directly from density
	if (EOSType > 0) { 
	  EOS(p, D_new, eint, h, cs, dpdrho, dpde, EOSType, 0);
	  v2 = vx*vx + vy*vy + vz*vz;
	  B2 = Bx_new*Bx_new + By_new*By_new + Bz_new*Bz_new;
	  etot = eint + 0.5*v2 + 0.5*B2/D_new;
	}
	
	BaryonField[iden ][igrid] = D_new;
	BaryonField[ivx  ][igrid] = vx;
	BaryonField[ivy  ][igrid] = vy;
	BaryonField[ivz  ][igrid] = vz;
	BaryonField[ietot][igrid] = etot;
	BaryonField[iBx  ][igrid] = Bx_new;
	BaryonField[iBy  ][igrid] = By_new;
	BaryonField[iBz  ][igrid] = Bz_new;
	BaryonField[iPhi ][igrid] = Phi_new*exp(-c1*dtFixed*pow(C_h/C_p,2));
	if (DualEnergyFormalism) {
	  v2 = vx*vx + vy*vy + vz*vz;
	  B2 = Bx_new*Bx_new + By_new*By_new + Bz_new*Bz_new;
	  //float rhou, lenu, tempu, tu, velu;
	  //GetUnits(&rhou, &lenu, &tempu, &tu, &velu, Time);
	  //float eintu = lenu*lenu/tu/tu;
	  float emin = SmallT/(Mu*(Gamma-1.0));
	  //float emin = 4.0*0.48999*D_new*pow(CellWidth[0][0],2)/(Gamma*(Gamma-1.0));

	  eint = Eint_new/D_new;
	  float eint1 = etot - 0.5*v2 - 0.5*B2/D_new;
	  if (eint1 > 0) {
	    EOS(p, D_new, eint, h, cs, dpdrho, dpde, EOSType, 2);
	  }
	  else {
	    cs = 0.0;
	  }
	  if (cs*cs > DualEnergyFormalismEta1*v2 &&  
	      cs*cs > DualEnergyFormalismEta1*B2/D_new &&
	      eint1 > 0.5*eint) {
	    eint = eint1;
	  }
	  /*if (eint1 > eint) 
	    eint = eint1;*/
	  eint = max(eint, emin);
	  BaryonField[ieint][igrid] = eint;
	  BaryonField[ietot][igrid] = eint + 0.5*v2 + 0.5*B2/D_new;
	  if (BaryonField[ieint][igrid] < 0.0) {
	    printf("UpdatePrim: eint < 0, cs2=%g, eta*v2=%g, eint=%g, etot=%g, 0.5*v2=%g, p=%g, rho=%g,0.5*B2/rho=%g\n", 
		   cs*cs, DualEnergyFormalismEta1*v2, eint, etot, 0.5*v2, p, D_new, 0.5*B2/rho);
	    printf("dU[%d]=%g, dU[ieint]=%g, eint_old=%g,eint1=%g\n", iEtot, dU[iEtot][n], dU[iEint][n], eint_old, eint1);
	    return FAIL;
	  }
	}
      }
    }
  }
  
  /* Convert species from mass fraction to density */ 

  for (int field = NEQ_MHD; field < NEQ_MHD+NSpecies; field++) {
    for (int n = 0; n < size; n++) {
      BaryonField[field][n] *= BaryonField[iden][n];
    }
  }

  if ( (NSpecies+NColor) > 0) {
    delete [] D;
    delete [] sum;
  }

  
  return SUCCESS;
}
