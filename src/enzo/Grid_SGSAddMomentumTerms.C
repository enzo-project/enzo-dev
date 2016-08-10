/***********************************************************************
 *
 * INFORMATION This file is part of a subgrid-scale (SGS) modeling 
 * framework in order to conduct explicit large eddy simulations (LES).
 *
 * The functions in this file concern models for the turbulent 
 * stress tensor in the momentum equation.
 * It consists of the turbulent (or SGS) Reynolds stress, the SGS
 * Maxwell stress, and the SGS magnetic pressure.
 *
 * The models have been verified "a priori", i.e. in comparison to
 * reference data, in 
 * Grete et al 2015 New J. Phys. 17 023070 doi: 10.1088/1367-2630/17/2/023070
 * Grete et al 2016 Phys. Plasmas 23 062317 doi: 10.1063/1.4954304 (Grete2016a)
 * and "a posteriori", i.e. used in simulations of decaying MHD turbulence, in
 * Grete et al ... under review ... (Grete201X)
 *
 * WRITTEN BY Philipp Grete (mail@pgrete.de)
 *
 * DATE 2016
 *
************************************************************************/

#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

/* 
 * This function adds to the SGS stress tensor (Reynolds stress component)
 * the pure (unscaled) nonlinear model
 * TauU = 1/12 * Delta^2 rho u_i,k u_j,k
 *
 * See equation (35) in Grete2016a for details (such as coefficient values)
 * or Vlaykov et al 2016 Phys. Plasmas 23 062316 doi: 10.1063/1.4954303 for
 * the derivation.
 */
void grid::SGSAddTauNLuTerm(float **Tau) {
  if (debug)
    printf("[%"ISYM"] grid::SGSAddTauNLuTerm start\n",MyProcessorNumber);

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  int B1Num, B2Num, B3Num, PhiNum;
  this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num,
      TENum, B1Num, B2Num, B3Num, PhiNum);

  float* rho;
  // if an explicit filter should be used
  // (at this point the fields are already filtered, 
  // see hydro_rk/Grid_MHDSourceTerms.C and the SGSNeedJacobians switch)
  if (SGSFilterWidth > 1.) {
    rho = FilteredFields[0];
  // if the model should be calculated based on grid-scale quantities
  // (not recommended, see Grete201X)
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


  // the combined prefactor
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

/* 
 * This function adds to the SGS stress tensor (Reynolds stress component)
 * the scaled nonlinear model (magnitude is scaled by an kinetic SGS energy
 * model based on realizability conditions): *
 * TauU =  2 C Delta^2 rho |S*|^2 (u_i,k u_j,k)/(u_l,s u_l,s)
 *
 * See equation (43) and (13) in Grete2016a for details (such as coefficient values)
 */
void grid::SGSAddTauNLuNormedEnS2StarTerm(float **Tau) {
  if (debug)
    printf("[%"ISYM"] grid::SGSAddTauNLuNormedEnS2StarTerm start\n",MyProcessorNumber);

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  int B1Num, B2Num, B3Num, PhiNum;
  this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num,
      TENum, B1Num, B2Num, B3Num, PhiNum);

  float* rho;
  // if an explicit filter should be used
  // (at this point the fields are already filtered, 
  // see hydro_rk/Grid_MHDSourceTerms.C and the SGSNeedJacobians switch)
  if (SGSFilterWidth > 1.) {
    rho = FilteredFields[0];
  // if the model should be calculated based on grid-scale quantities
  // (not recommended, see Grete201X)
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

  // the combined prefactor
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

/* 
 * This function adds to the SGS stress tensor (Maxwell stress component) 
 * the pure (unscaled) nonlinear model for TauB (full)
 * TauB = 1/12 * Delta^2 B_i,k B_j,k
 *
 * See equation (36) in Grete2016a for details (such as coefficient values)
 * or Vlaykov et al 2016 Phys. Plasmas 23 062316 doi: 10.1063/1.4954303 for
 * the derivation.
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


  // the combined prefactor
  float CDeltaSqr = 1./12. * SGScoeffNLb * pow(SGSFilterWidth,2.) * 
    pow(CellWidth[0][0]*CellWidth[1][0]*CellWidth[2][0],2./3.);

  int igrid;
  float turbMagPres;

  for (int k = StartIndex[2]; k <= EndIndex[2]; k++)
    for (int j = StartIndex[1]; j <= EndIndex[1]; j++)
      for (int i = StartIndex[0]; i <= EndIndex[0]; i++) {

        igrid = i + (j+k*GridDimension[1])*GridDimension[0];

        turbMagPres = 0.;

        // the pure Maxwell stress component
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

        // the turbulent magnetic pressure component
        Tau[XX][igrid] += CDeltaSqr * turbMagPres/2.;
        Tau[YY][igrid] += CDeltaSqr * turbMagPres/2.;
        Tau[ZZ][igrid] += CDeltaSqr * turbMagPres/2.;

      }
}


/* 
 * This function adds to the SGS stress tensor (Reynolds stress component)
 * the eddy viscosity model where the strength of eddy turbulent viscosity is 
 * scaled by the kinetic SGS energy as given by a model based on 
 * realizability conditions:
 * Tau = -2 C_1 Delta^2 rho |S*| S* + 2/3 C_2 delta_ij Delta^2 rho |S*|^2
 *
 * See equation (10), (21) and (13) in Grete2016a for details (such as coefficient 
 * values) or (in practice) equations (8), (10) and (12) in Grete201X.
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
  // if an explicit filter should be used
  // (at this point the fields are already filtered, 
  // see hydro_rk/Grid_MHDSourceTerms.C and the SGSNeedJacobians switch)
    rho = FilteredFields[0];
  // if the model should be calculated based on grid-scale quantities
  // (not recommended, see Grete201X)
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


  // the combined prefactor of the deviatoric part
  float Minus2C1DeltaSqr = -2. * SGScoeffEVStarEnS2Star * pow(SGSFilterWidth,2.) *
    pow(CellWidth[0][0]*CellWidth[1][0]*CellWidth[2][0],2./3.);
  // the combined prefactor of the isotropic part
  float TwoThirdC2DeltaSqr = 2./3. * SGScoeffEnS2StarTrace * pow(SGSFilterWidth,2.) *
    pow(CellWidth[0][0]*CellWidth[1][0]*CellWidth[2][0],2./3.);

  int igrid;
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

/* 
 * This function adds to the SGS stress tensor (Reynolds stress component)
 * a scale-similarity motivated term
 * TauU = flt(rho) * (flt(u_i u_j) - flt(u_i) * flt(u_j))
 *
 * See equation (30) in Grete2016a for details (such as coefficient values)
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

/* 
 * This function adds to the SGS stress tensor (Maxwell stress component)
 * a scale-similarity motivated term
 * TauU = (flt(B_i B_j) - flt(B_i) * flt(B_j))
 *
 * See equation (31) in Grete2016a for details (such as coefficient values)
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

/*
 * This function initializes a zero stress tensor and calls the individual
 * functions that add the different terms to it.
 * Finally, the divergence of the tensor is added to the dU vector used by
 * the MUSCL framework in hydro_rk/Grid_MHDSourceTerms.C 
 */
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
  // if an explicit filter should be used
  // (at this point the fields are already filtered, 
  // see hydro_rk/Grid_MHDSourceTerms.C and the SGSNeedJacobians switch)
  if (SGSFilterWidth > 1.) {
    rho = FilteredFields[0];
  } else {
  // if the model should be calculated based on grid-scale quantities
  // (not recommended, see Grete201X)
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


  // the individual terms are added/activated by a non-zero coefficient
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

  for (int dim = 0; dim < 6; dim++) {
    delete [] Tau[dim];
  }

  return SUCCESS;
}
