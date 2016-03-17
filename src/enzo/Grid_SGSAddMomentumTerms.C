#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

/* pure (unscaled) nonlinear model for TauU (full)
 * TauU = 1/12 * Delta^2 rho u_i,k u_j,k
 * see eq TODO of TODO for details
 */
void grid::SGSAddTauNLuTerm(float **Tau) {
  if (debug)
    printf("[%"ISYM"] grid::SGSAddTauNLuTerm start\n",MyProcessorNumber);

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  int B1Num, B2Num, B3Num, PhiNum;
  this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num,
      TENum, B1Num, B2Num, B3Num, PhiNum);

  float* rho;
  if (SGSFilterWidth > 1.) {
    rho = FilteredFields[0];
  } else {
    rho = BaryonField[DensNum];
  }

  int size = 1;
  int StartIndex[MAX_DIMENSION];
  int EndIndex[MAX_DIMENSION];

  for (int dim = 0; dim < MAX_DIMENSION; dim++) {
    size *= GridDimension[dim];

    /* we need Tau in the first ghost zone as well
     * as we'll take another derivative later on */
    StartIndex[dim] = GridStartIndex[dim] - 1;
    EndIndex[dim] = GridEndIndex[dim] + 1;
  }


  float CDeltaSqr = 1./12. * SGScoeffNLu * pow(SGSFilterWidth,2.) *
    pow(CellWidth[0][0]*CellWidth[1][0]*CellWidth[2][0],2./3.);

  int igrid;

  for (int k = StartIndex[2]; k <= EndIndex[2]; k++)
    for (int j = StartIndex[1]; j <= EndIndex[1]; j++)
      for (int i = StartIndex[0]; i <= EndIndex[0]; i++) {

        igrid = i + (j+k*GridDimension[1])*GridDimension[0];

        for (int l = 0; l < MAX_DIMENSION; l++) {
          Tau[XX][igrid] += CDeltaSqr * rho[igrid] * JacVel[X][l][igrid] * JacVel[X][l][igrid];
          Tau[YY][igrid] += CDeltaSqr * rho[igrid] * JacVel[Y][l][igrid] * JacVel[Y][l][igrid];
          Tau[ZZ][igrid] += CDeltaSqr * rho[igrid] * JacVel[Z][l][igrid] * JacVel[Z][l][igrid];
          Tau[XY][igrid] += CDeltaSqr * rho[igrid] * JacVel[X][l][igrid] * JacVel[Y][l][igrid];
          Tau[YZ][igrid] += CDeltaSqr * rho[igrid] * JacVel[Y][l][igrid] * JacVel[Z][l][igrid];
          Tau[XZ][igrid] += CDeltaSqr * rho[igrid] * JacVel[X][l][igrid] * JacVel[Z][l][igrid];
        }

      }

}

/* nonlinear model for TauU (full) and scaled by realiz. energy
 * TauU =  2 C Delta^2 rho |S*|^2 (u_i,k u_j,k)/(u_l,s u_l,s)
 * see eq TODO of TODO for details
 */
void grid::SGSAddTauNLuNormedEnS2StarTerm(float **Tau) {
  if (debug)
    printf("[%"ISYM"] grid::SGSAddTauNLuNormedEnS2StarTerm start\n",MyProcessorNumber);

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  int B1Num, B2Num, B3Num, PhiNum;
  this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num,
      TENum, B1Num, B2Num, B3Num, PhiNum);

  float* rho;
  if (SGSFilterWidth > 1.) {
    rho = FilteredFields[0];
  } else {
    rho = BaryonField[DensNum];
  }

  int size = 1;
  int StartIndex[MAX_DIMENSION];
  int EndIndex[MAX_DIMENSION];

  for (int dim = 0; dim < MAX_DIMENSION; dim++) {
    size *= GridDimension[dim];

    /* we need Tau in the first ghost zone as well
     * as we'll take another derivative later on */
    StartIndex[dim] = GridStartIndex[dim] - 1;
    EndIndex[dim] = GridEndIndex[dim] + 1;
  }


  float TwoCDeltaSqr = 2. * SGScoeffNLuNormedEnS2Star * pow(SGSFilterWidth,2.) *
    pow(CellWidth[0][0]*CellWidth[1][0]*CellWidth[2][0],2./3.);
  float traceSthird;
  float SStarSqr;
  float JacNorm;
  float prefactor;

  int igrid;

  for (int k = StartIndex[2]; k <= EndIndex[2]; k++)
    for (int j = StartIndex[1]; j <= EndIndex[1]; j++)
      for (int i = StartIndex[0]; i <= EndIndex[0]; i++) {

        igrid = i + (j+k*GridDimension[1])*GridDimension[0];

        traceSthird = (JacVel[X][X][igrid] + JacVel[Y][Y][igrid] + JacVel[Z][Z][igrid])/3.;

        SStarSqr = (
            2.*(pow(JacVel[X][X][igrid]-traceSthird,2.) +
              pow(JacVel[Y][Y][igrid]-traceSthird,2.) +
              pow(JacVel[Z][Z][igrid]-traceSthird,2.)
              )
            + pow(JacVel[X][Y][igrid] + JacVel[Y][X][igrid],2.)
            + pow(JacVel[Y][Z][igrid] + JacVel[Z][Y][igrid],2.)
            + pow(JacVel[X][Z][igrid] + JacVel[Z][X][igrid],2.));

        JacNorm = 0.;
        for (int l = 0; l < MAX_DIMENSION; l++)
          for (int s = 0; s < MAX_DIMENSION; s++)
            JacNorm += JacVel[l][s][igrid] * JacVel[l][s][igrid];

        prefactor = TwoCDeltaSqr * rho[igrid] * SStarSqr / JacNorm;

        for (int l = 0; l < MAX_DIMENSION; l++) {
          Tau[XX][igrid] += prefactor * JacVel[X][l][igrid] * JacVel[X][l][igrid];
          Tau[YY][igrid] += prefactor * JacVel[Y][l][igrid] * JacVel[Y][l][igrid];
          Tau[ZZ][igrid] += prefactor * JacVel[Z][l][igrid] * JacVel[Z][l][igrid];
          Tau[XY][igrid] += prefactor * JacVel[X][l][igrid] * JacVel[Y][l][igrid];
          Tau[YZ][igrid] += prefactor * JacVel[Y][l][igrid] * JacVel[Z][l][igrid];
          Tau[XZ][igrid] += prefactor * JacVel[X][l][igrid] * JacVel[Z][l][igrid];
        }

      }

}

/* pure (unscaled) nonlinear model for TauB (full)
 * TauB = 1/12 * Delta^2 B_i,k B_j,k
 * see eq TODO of TODO for details
 */
void grid::SGSAddTauNLbTerm(float **Tau) {
  if (debug)
    printf("[%"ISYM"] grid::SGSAddTauNLbTerm start\n",MyProcessorNumber);

  int size = 1;
  int StartIndex[MAX_DIMENSION];
  int EndIndex[MAX_DIMENSION];

  for (int dim = 0; dim < MAX_DIMENSION; dim++) {
    size *= GridDimension[dim];

    /* we need Tau in the first ghost zone as well
     * as we'll take another derivative later on */
    StartIndex[dim] = GridStartIndex[dim] - 1;
    EndIndex[dim] = GridEndIndex[dim] + 1;
  }


  float CDeltaSqr = 1./12. * SGScoeffNLb * pow(SGSFilterWidth,2.) * 
    pow(CellWidth[0][0]*CellWidth[1][0]*CellWidth[2][0],2./3.);

  int igrid;
  float turbMagPres;

  for (int k = StartIndex[2]; k <= EndIndex[2]; k++)
    for (int j = StartIndex[1]; j <= EndIndex[1]; j++)
      for (int i = StartIndex[0]; i <= EndIndex[0]; i++) {

        igrid = i + (j+k*GridDimension[1])*GridDimension[0];

        turbMagPres = 0.;
        for (int l = 0; l < MAX_DIMENSION; l++) {
          turbMagPres += JacB[X][l][igrid] * JacB[X][l][igrid] 
            + JacB[Y][l][igrid] * JacB[Y][l][igrid]
            + JacB[Z][l][igrid] * JacB[Z][l][igrid];

          Tau[XX][igrid] -= CDeltaSqr * JacB[X][l][igrid] * JacB[X][l][igrid];
          Tau[YY][igrid] -= CDeltaSqr * JacB[Y][l][igrid] * JacB[Y][l][igrid];
          Tau[ZZ][igrid] -= CDeltaSqr * JacB[Z][l][igrid] * JacB[Z][l][igrid];
          Tau[XY][igrid] -= CDeltaSqr * JacB[X][l][igrid] * JacB[Y][l][igrid];
          Tau[YZ][igrid] -= CDeltaSqr * JacB[Y][l][igrid] * JacB[Z][l][igrid];
          Tau[XZ][igrid] -= CDeltaSqr * JacB[X][l][igrid] * JacB[Z][l][igrid];
        }

        Tau[XX][igrid] += CDeltaSqr * turbMagPres/2.;
        Tau[YY][igrid] += CDeltaSqr * turbMagPres/2.;
        Tau[ZZ][igrid] += CDeltaSqr * turbMagPres/2.;

      }

}


/* eddy viscosity model for full tau scaled by realiz. energies
 * Tau = -2 C_1 Delta^2 rho |S*| S* + 2/3 C_2 delta_ij Delta^2 rho |S*|^2
 * see eq TODO of TODO for details
 */
void grid::SGSAddTauEVEnS2StarTerm(float **Tau) {
  if (debug)
    printf("[%"ISYM"] grid::SGSAddTauEVEnS2StarTerm start\n",MyProcessorNumber);

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  int B1Num, B2Num, B3Num, PhiNum;
  this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num,
      TENum, B1Num, B2Num, B3Num, PhiNum);

  float* rho;
  if (SGSFilterWidth > 1.) {
    rho = FilteredFields[0];
  } else {
    rho = BaryonField[DensNum];
  }

  int size = 1;
  int StartIndex[MAX_DIMENSION];
  int EndIndex[MAX_DIMENSION];

  for (int dim = 0; dim < MAX_DIMENSION; dim++) {
    size *= GridDimension[dim];

    /* we need Tau in the first ghost zone as well
     * as we'll take another derivative later on */
    StartIndex[dim] = GridStartIndex[dim] - 1;
    EndIndex[dim] = GridEndIndex[dim] + 1;
  }


  float Minus2C1DeltaSqr = -2. * SGScoeffEVStarEnS2Star * pow(SGSFilterWidth,2.) *
    pow(CellWidth[0][0]*CellWidth[1][0]*CellWidth[2][0],2./3.);
  float TwoThirdC2DeltaSqr = 2./3. * SGScoeffEnS2StarTrace * pow(SGSFilterWidth,2.) *
    pow(CellWidth[0][0]*CellWidth[1][0]*CellWidth[2][0],2./3.);

  int igrid;

  /* magic with S |S| S*... could potentially handled by
   * external function, should reduce CPU time, but increase memory usage
   */
  float traceSthird;
  float SStarSqr;

  for (int k = StartIndex[2]; k <= EndIndex[2]; k++)
    for (int j = StartIndex[1]; j <= EndIndex[1]; j++)
      for (int i = StartIndex[0]; i <= EndIndex[0]; i++) {

        igrid = i + (j+k*GridDimension[1])*GridDimension[0];

        traceSthird = (JacVel[X][X][igrid] + JacVel[Y][Y][igrid] + JacVel[Z][Z][igrid])/3.;

        SStarSqr = (
            2.*(pow(JacVel[X][X][igrid]-traceSthird,2.) +
              pow(JacVel[Y][Y][igrid]-traceSthird,2.) +
              pow(JacVel[Z][Z][igrid]-traceSthird,2.)
              )
            + pow(JacVel[X][Y][igrid] + JacVel[Y][X][igrid],2.)
            + pow(JacVel[Y][Z][igrid] + JacVel[Z][Y][igrid],2.)
            + pow(JacVel[X][Z][igrid] + JacVel[Z][X][igrid],2.));


        Tau[XX][igrid] += Minus2C1DeltaSqr * rho[igrid] * pow(SStarSqr,1./2.) * (
            JacVel[X][X][igrid] - traceSthird) + TwoThirdC2DeltaSqr * rho[igrid] * SStarSqr;
        Tau[YY][igrid] += Minus2C1DeltaSqr * rho[igrid] * pow(SStarSqr,1./2.) * (
            JacVel[Y][Y][igrid] - traceSthird) + TwoThirdC2DeltaSqr * rho[igrid] * SStarSqr;
        Tau[ZZ][igrid] += Minus2C1DeltaSqr * rho[igrid] * pow(SStarSqr,1./2.) * (
            JacVel[Z][Z][igrid] - traceSthird) + TwoThirdC2DeltaSqr * rho[igrid] * SStarSqr;

        Tau[XY][igrid] += Minus2C1DeltaSqr * rho[igrid] * pow(SStarSqr,1./2.) * (
            JacVel[X][Y][igrid] + JacVel[Y][X][igrid])/2.;
        Tau[YZ][igrid] += Minus2C1DeltaSqr * rho[igrid] * pow(SStarSqr,1./2.) * (
            JacVel[Y][Z][igrid] + JacVel[Z][Y][igrid])/2.;
        Tau[ZX][igrid] += Minus2C1DeltaSqr * rho[igrid] * pow(SStarSqr,1./2.) * (
            JacVel[Z][X][igrid] + JacVel[X][Z][igrid])/2.;

      }
}

/* scale-similarity model for TauU 
 * TauU = flt(rho) * (flt(u_i u_j) - flt(u_i) * flt(u_j))
 * see eq TODO of TODO for details
 */
void grid::SGSAddTauSSuTerm(float **Tau) {
  if (debug)
    printf("[%"ISYM"] grid::SGSAddTauSSuTerm start\n",MyProcessorNumber);

  int size = 1;
  int StartIndex[MAX_DIMENSION];
  int EndIndex[MAX_DIMENSION];

  for (int dim = 0; dim < MAX_DIMENSION; dim++) {
    size *= GridDimension[dim];

    /* we need Tau in the first ghost zone as well
     * as we'll take another derivative later on */
    StartIndex[dim] = GridStartIndex[dim] - 1;
    EndIndex[dim] = GridEndIndex[dim] + 1;
  }


  int igrid;

  for (int k = StartIndex[2]; k <= EndIndex[2]; k++)
    for (int j = StartIndex[1]; j <= EndIndex[1]; j++)
      for (int i = StartIndex[0]; i <= EndIndex[0]; i++) {

        igrid = i + (j+k*GridDimension[1])*GridDimension[0];

        Tau[XX][igrid] += SGScoeffSSu * (FltrhoUU[XX][igrid] - 
            FilteredFields[0][igrid] * FilteredFields[1][igrid] * FilteredFields[1][igrid]);
        Tau[YY][igrid] += SGScoeffSSu * (FltrhoUU[YY][igrid] - 
            FilteredFields[0][igrid] * FilteredFields[2][igrid] * FilteredFields[2][igrid]);
        Tau[ZZ][igrid] += SGScoeffSSu * (FltrhoUU[ZZ][igrid] - 
            FilteredFields[0][igrid] * FilteredFields[3][igrid] * FilteredFields[3][igrid]);
        Tau[XY][igrid] += SGScoeffSSu * (FltrhoUU[XY][igrid] - 
            FilteredFields[0][igrid] * FilteredFields[1][igrid] * FilteredFields[2][igrid]);
        Tau[YZ][igrid] += SGScoeffSSu * (FltrhoUU[YZ][igrid] - 
            FilteredFields[0][igrid] * FilteredFields[2][igrid] * FilteredFields[3][igrid]);
        Tau[XZ][igrid] += SGScoeffSSu * (FltrhoUU[XZ][igrid] - 
            FilteredFields[0][igrid] * FilteredFields[1][igrid] * FilteredFields[3][igrid]);

      }

}

/* scale-similarity model for TauB 
 * TauU = (flt(B_i B_j) - flt(B_i) * flt(B_j))
 * see eq TODO of TODO for details
 */
void grid::SGSAddTauSSbTerm(float **Tau) {
  if (debug)
    printf("[%"ISYM"] grid::SGSAddTauSSbTerm start\n",MyProcessorNumber);

  int size = 1;
  int StartIndex[MAX_DIMENSION];
  int EndIndex[MAX_DIMENSION];

  for (int dim = 0; dim < MAX_DIMENSION; dim++) {
    size *= GridDimension[dim];

    /* we need Tau in the first ghost zone as well
     * as we'll take another derivative later on */
    StartIndex[dim] = GridStartIndex[dim] - 1;
    EndIndex[dim] = GridEndIndex[dim] + 1;
  }


  int igrid;

  for (int k = StartIndex[2]; k <= EndIndex[2]; k++)
    for (int j = StartIndex[1]; j <= EndIndex[1]; j++)
      for (int i = StartIndex[0]; i <= EndIndex[0]; i++) {

        igrid = i + (j+k*GridDimension[1])*GridDimension[0];

        Tau[XX][igrid] += SGScoeffSSb * (FltBB[XX][igrid] - 
            FilteredFields[4][igrid] * FilteredFields[4][igrid]);
        Tau[YY][igrid] += SGScoeffSSb * (FltBB[YY][igrid] - 
            FilteredFields[5][igrid] * FilteredFields[5][igrid]);
        Tau[ZZ][igrid] += SGScoeffSSb * (FltBB[ZZ][igrid] - 
            FilteredFields[6][igrid] * FilteredFields[6][igrid]);
        Tau[XY][igrid] += SGScoeffSSb * (FltBB[XY][igrid] - 
            FilteredFields[4][igrid] * FilteredFields[5][igrid]);
        Tau[YZ][igrid] += SGScoeffSSb * (FltBB[YZ][igrid] - 
            FilteredFields[5][igrid] * FilteredFields[6][igrid]);
        Tau[XZ][igrid] += SGScoeffSSb * (FltBB[XZ][igrid] - 
            FilteredFields[4][igrid] * FilteredFields[6][igrid]);

      }

}

int grid::SGSAddMomentumTerms(float **dU) {
  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }

  if (Time == 0.)
    return SUCCESS;

  if (debug)
    printf("[%"ISYM"] grid::SGSAddMomentumTerms start\n",MyProcessorNumber);

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  int B1Num, B2Num, B3Num, PhiNum;
  this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num,
      TENum, B1Num, B2Num, B3Num, PhiNum);

  float* rho;
  if (SGSFilterWidth > 1.) {
    rho = FilteredFields[0];
  } else {
    rho = BaryonField[DensNum];
  }

  int size = 1;
  float *Tau[6];

  for (int dim = 0; dim < MAX_DIMENSION; dim++) {
    size *= GridDimension[dim];
  }

  for (int dim = 0; dim < 6; dim++) {
    Tau[dim] = new float[size];
    for (int i = 0; i < size; i++)
      Tau[dim][i] = 0.;
  }


  if (SGScoeffNLu != 0.) 
    SGSAddTauNLuTerm(Tau);

  if (SGScoeffNLb != 0.) 
    SGSAddTauNLbTerm(Tau);

  if ((SGScoeffEVStarEnS2Star != 0.) || (SGScoeffEnS2StarTrace != 0.))
    SGSAddTauEVEnS2StarTerm(Tau);

  if (SGScoeffNLuNormedEnS2Star != 0.)
    SGSAddTauNLuNormedEnS2StarTerm(Tau);
  
  if (SGScoeffSSu != 0.) 
    SGSAddTauSSuTerm(Tau);
  
  if (SGScoeffSSb != 0.) 
    SGSAddTauSSbTerm(Tau);


  int n = 0;
  int igrid, ip1, im1, jp1, jm1, kp1, km1;
  float MomxIncr,MomyIncr,MomzIncr,EtotIncr;

  float facX = 1. / (2. * CellWidth[0][0]);
  float facY = 1. / (2. * CellWidth[1][0]);
  float facZ = 1. / (2. * CellWidth[2][0]);

  for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
    for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
      for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, n++) {

        igrid = i + (j+k*GridDimension[1])*GridDimension[0];
        ip1 = i+1 + (j+k*GridDimension[1])*GridDimension[0];
        im1 = i-1 + (j+k*GridDimension[1])*GridDimension[0];
        jp1 = i + (j+1+k*GridDimension[1])*GridDimension[0];
        jm1 = i + (j-1+k*GridDimension[1])*GridDimension[0];
        kp1 = i + (j+(k+1)*GridDimension[1])*GridDimension[0];
        km1 = i + (j+(k-1)*GridDimension[1])*GridDimension[0];


        MomxIncr = - dtFixed * (
            (Tau[XX][ip1] - Tau[XX][im1])*facX + 
            (Tau[XY][jp1] - Tau[XY][jm1])*facY + 
            (Tau[XZ][kp1] - Tau[XZ][km1])*facZ);
        EtotIncr = BaryonField[Vel1Num][igrid] * MomxIncr + 0.5 / rho[igrid] * MomxIncr * MomxIncr;

        MomyIncr = - dtFixed * (
            (Tau[YX][ip1] - Tau[YX][im1])*facX + 
            (Tau[YY][jp1] - Tau[YY][jm1])*facY + 
            (Tau[YZ][kp1] - Tau[YZ][km1])*facZ);
        EtotIncr += BaryonField[Vel2Num][igrid] * MomyIncr + 0.5 / rho[igrid] * MomyIncr * MomyIncr;

        MomzIncr = - dtFixed * (
            (Tau[ZX][ip1] - Tau[ZX][im1])*facX + 
            (Tau[ZY][jp1] - Tau[ZY][jm1])*facY + 
            (Tau[ZZ][kp1] - Tau[ZZ][km1])*facZ);
        EtotIncr += BaryonField[Vel3Num][igrid] * MomzIncr + 0.5 / rho[igrid] * MomzIncr * MomzIncr;

        dU[ivx][n] += MomxIncr;
        dU[ivy][n] += MomyIncr;
        dU[ivz][n] += MomzIncr;
        dU[iEtot][n] += EtotIncr;
      }

  if (debug)
    printf("[%"ISYM"] grid::SGSAddMomentumTerms end, last incr: %"FSYM" %"FSYM" %"FSYM" %"FSYM"\n",
        MyProcessorNumber,MomxIncr,MomyIncr,MomzIncr,EtotIncr);

  for (int dim = 0; dim < 6; dim++) {
    delete [] Tau[dim];
  }

  return SUCCESS;
}
