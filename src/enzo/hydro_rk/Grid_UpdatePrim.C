/***********************************************************************
/
/  GRID CLASS (UPDATE HYDRO VARIABLES)
/
/  written by: Peng Wang
/  date:       April, 2007
/  modified1:
/
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
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

int grid::UpdatePrim(float **dU, float c1, float c2)
{

  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }

  int i, j, k, n, dim, igrid, field, size, activesize;
  for (dim = 0, size = 1; dim < GridRank; dim++) {
    size *= GridDimension[dim];
  }

  for (dim = 0, activesize = 1; dim < GridRank; dim++) {
    activesize *= (GridDimension[dim] - 2*DEFAULT_GHOST_ZONES);
  }

  float *D, *sum;
  float SmallX = 1e-20;

  if ( (NSpecies+NColor) > 0) {
    D = new float[activesize];
    sum = new float[activesize];
    for (i = 0; i < activesize; i++) {
      D[i] = 0.0;
      sum[i] = 0.0;
    }
  }

  float *Prim[NEQ_HYDRO+NSpecies+NColor];
  float *OldPrim[NEQ_HYDRO+NSpecies+NColor];
  this->ReturnHydroRKPointers(Prim, false);
  this->ReturnOldHydroRKPointers(OldPrim, false);

  // update species and colours

  for (field = NEQ_HYDRO; field < NEQ_HYDRO+NSpecies+NColor; field++) {
    n = 0;
    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	igrid = (k * GridDimension[1] + j) * GridDimension[0] + GridStartIndex[0];
        for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++, igrid++) {
          Prim[field][igrid] = c1*OldPrim[field][igrid] +
            (1-c1)*Prim[field][igrid]*Prim[iden][igrid] + c2*dU[field][n];
          D[n] += Prim[field][igrid];
        }
      }
    }
  }

  // renormalize species and colours

  for (field = NEQ_HYDRO; field < NEQ_HYDRO+NSpecies+NColor; field++) {
    n = 0;
    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	igrid = (k * GridDimension[1] + j) * GridDimension[0] + GridStartIndex[0];
        for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++, igrid++) {
          Prim[field][igrid] = min(1.0, max((Prim[field][igrid]/D[n]), SmallX));
	  Prim[field][igrid] = Prim[field][igrid]/D[n];
          sum[n] += Prim[field][igrid];
        }
      }
    }
  }

  for (field = NEQ_HYDRO; field < NEQ_HYDRO+NSpecies+NColor; field++) {
    n = 0;
    for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	igrid = (k * GridDimension[1] + j) * GridDimension[0] + GridStartIndex[0];
        for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++, igrid++)
          Prim[field][igrid] /= sum[n];
      }
    }
  }



  // update conserved variables
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, 
    B1Num, B2Num, B3Num;

  this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				   Vel3Num, TENum, B1Num, B2Num, B3Num);

  float rho_old, vx_old, vy_old, vz_old, e_old, etot_old, Tau_old, eint_old,
    rho, vx, vy, vz, e, etot, Tau, eint, p, v2,
    D_new, S1_new, S2_new, S3_new, Tau_new, h, cs, dpdrho, dpde, Eint_new;

  n = 0;
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {
	// first convert to conserved variables to do the update
	igrid = (k * GridDimension[1] + j) * GridDimension[0] + i;
	rho_old  = OldBaryonField[DensNum][igrid];
	vx_old   = OldBaryonField[Vel1Num][igrid];
	vy_old   = OldBaryonField[Vel2Num][igrid];
	vz_old   = OldBaryonField[Vel3Num][igrid];
	etot_old = OldBaryonField[TENum][igrid];
	if (DualEnergyFormalism) {
	  eint_old = OldBaryonField[GENum][igrid];
	}
	Tau_old = rho_old*etot_old;

	rho  = BaryonField[DensNum][igrid];
	vx   = BaryonField[Vel1Num][igrid];
	vy   = BaryonField[Vel2Num][igrid];
	vz   = BaryonField[Vel3Num][igrid];
	etot = BaryonField[TENum][igrid];
	if (DualEnergyFormalism) {
	  eint = BaryonField[GENum][igrid];
	}
	Tau = rho*etot;

	D_new   = c1*rho_old + (1.0-c1)*rho + c2*dU[iD][n];
	S1_new  = c1*rho_old*vx_old + (1.0-c1)*rho*vx + c2*dU[iS1][n];
	S2_new  = c1*rho_old*vy_old + (1.0-c1)*rho*vy + c2*dU[iS2][n];
	S3_new  = c1*rho_old*vz_old + (1.0-c1)*rho*vz + c2*dU[iS3][n];
	Tau_new = c1*Tau_old + (1.0-c1)*Tau + c2*dU[iEtot][n];

	if (DualEnergyFormalism) {
	  Eint_new = c1*rho_old*eint_old + (1.0-c1)*rho*eint + c2*dU[iEint][n];
	  /*if (Eint_new < 0) {
	    printf("UpdatePrim: eint < 0 in dual energy update\n");
	    printf("eint_old=%g,eint=%g,eint_new=%g, dU=%g\n",
		   eint_old, eint, Eint_new, dU[iEint][n]);
	    return FAIL;
	    }*/
	}

	if (D_new < 0 || isnan(D_new)) {
	  printf("UpdatePrim: rho <0 at %d %d %d: rho_old=%g, rho=%g, rho_new=%g, dU[iD]=%g\n", 
		 i, j, k, rho_old, rho, D_new, dU[iD][n]);
	  //D_new = max(D_new, SmallRho);
	  D_new = rho;
	  //D_new = rho;
	  return FAIL;
	}

	//D_new = max(D_new, SmallRho);

	// convert back to primitives
	vx = S1_new/D_new;
	vy = S2_new/D_new;
	vz = S3_new/D_new;
	etot = Tau_new/D_new;

	v2 = vx*vx + vy*vy + vz*vz;
	// If using polytropic EOS, calcuate etot using density
	if (EOSType > 0) { 
	  EOS(p, D_new, eint, h, cs, dpdrho, dpde, EOSType, 0);
	  etot = eint + 0.5*v2;
	}


	if (etot < 0 && EOSType == 0) {
	  float v2_old = vx_old*vx_old + vy_old*vy_old + vz_old*vz_old;
	  printf("UpdatePrim: tau < 0. etot_old=%g, etot=%g, etot_new=%g, v2=%g, v2old=%g, dU[iTau] = %g\n", 
		 Tau_old/rho_old, Tau/rho, Tau_new/D_new, v2, v2_old, dU[iEtot][n]*CellWidth[0][0]/dtFixed);
	  printf("rho_new=%g, rho=%g, rho_old=%g\n", D_new, rho, rho_old);
	  //return FAIL;
	}
	
	BaryonField[DensNum][igrid] = D_new;
	BaryonField[Vel1Num][igrid] = vx;
	BaryonField[Vel2Num][igrid] = vy;
	BaryonField[Vel3Num][igrid] = vz;
	BaryonField[TENum][igrid] = etot;

	if (DualEnergyFormalism) {
	  v2 = vx*vx + vy*vy + vz*vz;
	  eint = Eint_new/D_new;
	  float eint_du = eint;
	  float eint1 = etot - 0.5*v2;
          float emin = SmallT/(Mu*(Gamma-1.0));
	  //float emin = 4.0*0.48999*D_new*pow(CellWidth[0][0],2)/(Gamma*(Gamma-1.0));

	  if (eint1 > 0) {
	    EOS(p, D_new, eint1, h, cs, dpdrho, dpde, EOSType, 2);
	  }
	  else {
	    cs = 0.0;
	  }
	  if (cs*cs > DualEnergyFormalismEta1*v2 && eint1 > 0.5*eint) {
	    eint = eint1;
	  }
	  eint = max(eint, emin);
	  BaryonField[GENum][igrid] = eint;
	  BaryonField[TENum][igrid] = eint + 0.5*v2;
	  
	  if (eint < 0.0) {
	    printf("UpdatePrim:eint < 0, cs2=%g, eta*v2=%g, eint=%g, etot=%g, v2=%g, p=%g, rho=%g,eint1=%g\n", 
		   cs*cs, DualEnergyFormalismEta1*v2, eint_du, etot, v2, p, D_new, eint1);
	    printf("old rho=%g, old v:%g %g %g, old etot=%g oldeint=%g \n", 
		   rho_old, vx_old, vy_old, vz_old,
		   etot_old, eint_old);
	    printf("dU=%g %g %g %g %g, Tau_old=%g\n", 
		   dU[iD][n], dU[iS1][n], dU[iS2][n], dU[iS3][n], dU[iEtot][n], Tau_old);
	    //return FAIL;
	    BaryonField[GENum][igrid] = OldBaryonField[GENum][igrid];
	    BaryonField[TENum][igrid] = OldBaryonField[TENum][igrid];
	  }
	}
      }
    }
  }

  // convert species from mass fraction to density
  
  for (field = NEQ_HYDRO; field < NEQ_HYDRO+NSpecies+NColor; field++)
    for (n = 0; n < size; n++)
      Prim[field][n] *= BaryonField[DensNum][n];

  this->UpdateElectronDensity();

  if ( (NSpecies+NColor) > 0) {
    delete [] D;
    delete [] sum;
  }

  
  return SUCCESS;
}
