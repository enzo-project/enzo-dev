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
 * Grete et al 2017 Phys. Rev. E 05 033206 (Grete2017)
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
void grid::SGS_AddMom_nonlinear_kinetic(float **Tau) {
  if (debug1)
    printf("[%"ISYM"] grid::SGS_AddMom_nonlinear_kinetic start\n",MyProcessorNumber);

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
  // (not recommended, see Grete2017)
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
  float CDeltaSqr = 1./12. * SGScoeffNLu * POW(SGSFilterWidth,2.) *
    POW(CellWidth[0][0]*CellWidth[1][0]*CellWidth[2][0],2./3.);

  int igrid;

  for (int k = StartIndex[2]; k <= EndIndex[2]; k++)
    for (int j = StartIndex[1]; j <= EndIndex[1]; j++)
      for (int i = StartIndex[0]; i <= EndIndex[0]; i++) {

        igrid = i + (j+k*GridDimension[1])*GridDimension[0];

        for (int l = 0; l < MAX_DIMENSION; l++) {
          Tau[SGSXX][igrid] += CDeltaSqr * rho[igrid] * JacVel[SGSX][l][igrid] * JacVel[SGSX][l][igrid];
          Tau[SGSYY][igrid] += CDeltaSqr * rho[igrid] * JacVel[SGSY][l][igrid] * JacVel[SGSY][l][igrid];
          Tau[SGSZZ][igrid] += CDeltaSqr * rho[igrid] * JacVel[SGSZ][l][igrid] * JacVel[SGSZ][l][igrid];
          Tau[SGSXY][igrid] += CDeltaSqr * rho[igrid] * JacVel[SGSX][l][igrid] * JacVel[SGSY][l][igrid];
          Tau[SGSYZ][igrid] += CDeltaSqr * rho[igrid] * JacVel[SGSY][l][igrid] * JacVel[SGSZ][l][igrid];
          Tau[SGSXZ][igrid] += CDeltaSqr * rho[igrid] * JacVel[SGSX][l][igrid] * JacVel[SGSZ][l][igrid];
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
void grid::SGS_AddMom_nonlinear_kinetic_scaled(float **Tau) {
  if (debug1)
    printf("[%"ISYM"] grid::SGS_AddMom_nonlinear_kinetic_scaled start\n",MyProcessorNumber);

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
  // (not recommended, see Grete2017)
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
  float TwoCDeltaSqr = 2. * SGScoeffNLuNormedEnS2Star * POW(SGSFilterWidth,2.) *
    POW(CellWidth[0][0]*CellWidth[1][0]*CellWidth[2][0],2./3.);
  float traceSthird;
  float SStarSqr;
  float JacNorm;
  float prefactor;

  int igrid;

  for (int k = StartIndex[2]; k <= EndIndex[2]; k++)
    for (int j = StartIndex[1]; j <= EndIndex[1]; j++)
      for (int i = StartIndex[0]; i <= EndIndex[0]; i++) {

        igrid = i + (j+k*GridDimension[1])*GridDimension[0];

        traceSthird = (JacVel[SGSX][SGSX][igrid] + JacVel[SGSY][SGSY][igrid] + JacVel[SGSZ][SGSZ][igrid])/3.;

        SStarSqr = (
            2.*(POW(JacVel[SGSX][SGSX][igrid]-traceSthird,2.) +
              POW(JacVel[SGSY][SGSY][igrid]-traceSthird,2.) +
              POW(JacVel[SGSZ][SGSZ][igrid]-traceSthird,2.)
              )
            + POW(JacVel[SGSX][SGSY][igrid] + JacVel[SGSY][SGSX][igrid],2.)
            + POW(JacVel[SGSY][SGSZ][igrid] + JacVel[SGSZ][SGSY][igrid],2.)
            + POW(JacVel[SGSX][SGSZ][igrid] + JacVel[SGSZ][SGSX][igrid],2.));

        JacNorm = 0.;
        for (int l = 0; l < MAX_DIMENSION; l++)
          for (int s = 0; s < MAX_DIMENSION; s++)
            JacNorm += JacVel[l][s][igrid] * JacVel[l][s][igrid];

        prefactor = TwoCDeltaSqr * rho[igrid] * SStarSqr / JacNorm;

        for (int l = 0; l < MAX_DIMENSION; l++) {
          Tau[SGSXX][igrid] += prefactor * JacVel[SGSX][l][igrid] * JacVel[SGSX][l][igrid];
          Tau[SGSYY][igrid] += prefactor * JacVel[SGSY][l][igrid] * JacVel[SGSY][l][igrid];
          Tau[SGSZZ][igrid] += prefactor * JacVel[SGSZ][l][igrid] * JacVel[SGSZ][l][igrid];
          Tau[SGSXY][igrid] += prefactor * JacVel[SGSX][l][igrid] * JacVel[SGSY][l][igrid];
          Tau[SGSYZ][igrid] += prefactor * JacVel[SGSY][l][igrid] * JacVel[SGSZ][l][igrid];
          Tau[SGSXZ][igrid] += prefactor * JacVel[SGSX][l][igrid] * JacVel[SGSZ][l][igrid];
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
void grid::SGS_AddMom_nonliner_magnetic(float **Tau) {
  if (debug1)
    printf("[%"ISYM"] grid::SGS_AddMom_nonliner_magnetic start\n",MyProcessorNumber);

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
  float CDeltaSqr = 1./12. * SGScoeffNLb * POW(SGSFilterWidth,2.) *
    POW(CellWidth[0][0]*CellWidth[1][0]*CellWidth[2][0],2./3.);

  int igrid;
  float turbMagPres;

  for (int k = StartIndex[2]; k <= EndIndex[2]; k++)
    for (int j = StartIndex[1]; j <= EndIndex[1]; j++)
      for (int i = StartIndex[0]; i <= EndIndex[0]; i++) {

        igrid = i + (j+k*GridDimension[1])*GridDimension[0];

        turbMagPres = 0.;

        // the pure Maxwell stress component
        for (int l = 0; l < MAX_DIMENSION; l++) {
          turbMagPres += JacB[SGSX][l][igrid] * JacB[SGSX][l][igrid] 
            + JacB[SGSY][l][igrid] * JacB[SGSY][l][igrid]
            + JacB[SGSZ][l][igrid] * JacB[SGSZ][l][igrid];

          Tau[SGSXX][igrid] -= CDeltaSqr * JacB[SGSX][l][igrid] * JacB[SGSX][l][igrid];
          Tau[SGSYY][igrid] -= CDeltaSqr * JacB[SGSY][l][igrid] * JacB[SGSY][l][igrid];
          Tau[SGSZZ][igrid] -= CDeltaSqr * JacB[SGSZ][l][igrid] * JacB[SGSZ][l][igrid];
          Tau[SGSXY][igrid] -= CDeltaSqr * JacB[SGSX][l][igrid] * JacB[SGSY][l][igrid];
          Tau[SGSYZ][igrid] -= CDeltaSqr * JacB[SGSY][l][igrid] * JacB[SGSZ][l][igrid];
          Tau[SGSXZ][igrid] -= CDeltaSqr * JacB[SGSX][l][igrid] * JacB[SGSZ][l][igrid];
        }

        // the turbulent magnetic pressure component
        Tau[SGSXX][igrid] += CDeltaSqr * turbMagPres/2.;
        Tau[SGSYY][igrid] += CDeltaSqr * turbMagPres/2.;
        Tau[SGSZZ][igrid] += CDeltaSqr * turbMagPres/2.;

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
 * values) or (in practice) equations (8), (10) and (12) in Grete2017.
 */
void grid::SGS_AddMom_eddy_viscosity_scaled(float **Tau) {
  if (debug1)
    printf("[%"ISYM"] grid::SGS_AddMom_eddy_viscosity_scaled start\n",MyProcessorNumber);

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
  // (not recommended, see Grete2017)
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
  float Minus2C1DeltaSqr = -2. * SGScoeffEVStarEnS2Star * POW(SGSFilterWidth,2.) *
    POW(CellWidth[0][0]*CellWidth[1][0]*CellWidth[2][0],2./3.);
  // the combined prefactor of the isotropic part
  float TwoThirdC2DeltaSqr = 2./3. * SGScoeffEnS2StarTrace * POW(SGSFilterWidth,2.) *
    POW(CellWidth[0][0]*CellWidth[1][0]*CellWidth[2][0],2./3.);

  int igrid;
  float traceSthird;
  float SStarSqr;

  for (int k = StartIndex[2]; k <= EndIndex[2]; k++)
    for (int j = StartIndex[1]; j <= EndIndex[1]; j++)
      for (int i = StartIndex[0]; i <= EndIndex[0]; i++) {

        igrid = i + (j+k*GridDimension[1])*GridDimension[0];

        traceSthird = (JacVel[SGSX][SGSX][igrid] + JacVel[SGSY][SGSY][igrid] + JacVel[SGSZ][SGSZ][igrid])/3.;

        SStarSqr = (
            2.*(POW(JacVel[SGSX][SGSX][igrid]-traceSthird,2.) +
              POW(JacVel[SGSY][SGSY][igrid]-traceSthird,2.) +
              POW(JacVel[SGSZ][SGSZ][igrid]-traceSthird,2.)
              )
            + POW(JacVel[SGSX][SGSY][igrid] + JacVel[SGSY][SGSX][igrid],2.)
            + POW(JacVel[SGSY][SGSZ][igrid] + JacVel[SGSZ][SGSY][igrid],2.)
            + POW(JacVel[SGSX][SGSZ][igrid] + JacVel[SGSZ][SGSX][igrid],2.));


        Tau[SGSXX][igrid] += Minus2C1DeltaSqr * rho[igrid] * POW(SStarSqr,1./2.) * (
            JacVel[SGSX][SGSX][igrid] - traceSthird) + TwoThirdC2DeltaSqr * rho[igrid] * SStarSqr;
        Tau[SGSYY][igrid] += Minus2C1DeltaSqr * rho[igrid] * POW(SStarSqr,1./2.) * (
            JacVel[SGSY][SGSY][igrid] - traceSthird) + TwoThirdC2DeltaSqr * rho[igrid] * SStarSqr;
        Tau[SGSZZ][igrid] += Minus2C1DeltaSqr * rho[igrid] * POW(SStarSqr,1./2.) * (
            JacVel[SGSZ][SGSZ][igrid] - traceSthird) + TwoThirdC2DeltaSqr * rho[igrid] * SStarSqr;

        Tau[SGSXY][igrid] += Minus2C1DeltaSqr * rho[igrid] * POW(SStarSqr,1./2.) * (
            JacVel[SGSX][SGSY][igrid] + JacVel[SGSY][SGSX][igrid])/2.;
        Tau[SGSYZ][igrid] += Minus2C1DeltaSqr * rho[igrid] * POW(SStarSqr,1./2.) * (
            JacVel[SGSY][SGSZ][igrid] + JacVel[SGSZ][SGSY][igrid])/2.;
        Tau[SGSZX][igrid] += Minus2C1DeltaSqr * rho[igrid] * POW(SStarSqr,1./2.) * (
            JacVel[SGSZ][SGSX][igrid] + JacVel[SGSX][SGSZ][igrid])/2.;

      }
}

/* 
 * This function adds to the SGS stress tensor (Reynolds stress component)
 * a scale-similarity motivated term
 * TauU = flt(rho) * (flt(u_i u_j) - flt(u_i) * flt(u_j))
 *
 * See equation (30) in Grete2016a for details (such as coefficient values)
 */
void grid::SGS_AddMom_scale_similarity_kinetic(float **Tau) {
  if (debug1)
    printf("[%"ISYM"] grid::SGS_AddMom_scale_similarity_kinetic start\n",MyProcessorNumber);

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

        Tau[SGSXX][igrid] += SGScoeffSSu * (FltrhoUU[SGSXX][igrid] - 
            FilteredFields[0][igrid] * FilteredFields[1][igrid] * FilteredFields[1][igrid]);
        Tau[SGSYY][igrid] += SGScoeffSSu * (FltrhoUU[SGSYY][igrid] - 
            FilteredFields[0][igrid] * FilteredFields[2][igrid] * FilteredFields[2][igrid]);
        Tau[SGSZZ][igrid] += SGScoeffSSu * (FltrhoUU[SGSZZ][igrid] - 
            FilteredFields[0][igrid] * FilteredFields[3][igrid] * FilteredFields[3][igrid]);
        Tau[SGSXY][igrid] += SGScoeffSSu * (FltrhoUU[SGSXY][igrid] - 
            FilteredFields[0][igrid] * FilteredFields[1][igrid] * FilteredFields[2][igrid]);
        Tau[SGSYZ][igrid] += SGScoeffSSu * (FltrhoUU[SGSYZ][igrid] - 
            FilteredFields[0][igrid] * FilteredFields[2][igrid] * FilteredFields[3][igrid]);
        Tau[SGSXZ][igrid] += SGScoeffSSu * (FltrhoUU[SGSXZ][igrid] - 
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
void grid::SGS_AddMom_scale_similarity_magnetic(float **Tau) {
  if (debug1)
    printf("[%"ISYM"] grid::SGS_AddMom_scale_similarity_magnetic start\n",MyProcessorNumber);

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

        Tau[SGSXX][igrid] += SGScoeffSSb * (FltBB[SGSXX][igrid] - 
            FilteredFields[4][igrid] * FilteredFields[4][igrid]);
        Tau[SGSYY][igrid] += SGScoeffSSb * (FltBB[SGSYY][igrid] - 
            FilteredFields[5][igrid] * FilteredFields[5][igrid]);
        Tau[SGSZZ][igrid] += SGScoeffSSb * (FltBB[SGSZZ][igrid] - 
            FilteredFields[6][igrid] * FilteredFields[6][igrid]);
        Tau[SGSXY][igrid] += SGScoeffSSb * (FltBB[SGSXY][igrid] - 
            FilteredFields[4][igrid] * FilteredFields[5][igrid]);
        Tau[SGSYZ][igrid] += SGScoeffSSb * (FltBB[SGSYZ][igrid] - 
            FilteredFields[5][igrid] * FilteredFields[6][igrid]);
        Tau[SGSXZ][igrid] += SGScoeffSSb * (FltBB[SGSXZ][igrid] - 
            FilteredFields[4][igrid] * FilteredFields[6][igrid]);

      }

}

/*
 * This function initializes a zero stress tensor and calls the individual
 * functions that add the different terms to it.
 * Finally, the divergence of the tensor is added to the dU vector used by
 * the MUSCL framework in hydro_rk/Grid_MHDSourceTerms.C 
 */
int grid::SGS_AddMomentumTerms(float **dU) {
  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }

  if (Time == 0.)
    return SUCCESS;

  if (debug1)
    printf("[%"ISYM"] grid::SGS_AddMomentumTerms start\n",MyProcessorNumber);

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
  // (not recommended, see Grete2017)
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
    SGS_AddMom_nonlinear_kinetic(Tau);

  if (SGScoeffNLb != 0.) 
    SGS_AddMom_nonliner_magnetic(Tau);

  if ((SGScoeffEVStarEnS2Star != 0.) || (SGScoeffEnS2StarTrace != 0.))
    SGS_AddMom_eddy_viscosity_scaled(Tau);

  if (SGScoeffNLuNormedEnS2Star != 0.)
    SGS_AddMom_nonlinear_kinetic_scaled(Tau);
  
  if (SGScoeffSSu != 0.) 
    SGS_AddMom_scale_similarity_kinetic(Tau);
  
  if (SGScoeffSSb != 0.) 
    SGS_AddMom_scale_similarity_magnetic(Tau);


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
            (Tau[SGSXX][ip1] - Tau[SGSXX][im1])*facX + 
            (Tau[SGSXY][jp1] - Tau[SGSXY][jm1])*facY + 
            (Tau[SGSXZ][kp1] - Tau[SGSXZ][km1])*facZ);
        EtotIncr = BaryonField[Vel1Num][igrid] * MomxIncr + 0.5 / rho[igrid] * MomxIncr * MomxIncr;

        MomyIncr = - dtFixed * (
            (Tau[SGSYX][ip1] - Tau[SGSYX][im1])*facX + 
            (Tau[SGSYY][jp1] - Tau[SGSYY][jm1])*facY + 
            (Tau[SGSYZ][kp1] - Tau[SGSYZ][km1])*facZ);
        EtotIncr += BaryonField[Vel2Num][igrid] * MomyIncr + 0.5 / rho[igrid] * MomyIncr * MomyIncr;

        MomzIncr = - dtFixed * (
            (Tau[SGSZX][ip1] - Tau[SGSZX][im1])*facX + 
            (Tau[SGSZY][jp1] - Tau[SGSZY][jm1])*facY + 
            (Tau[SGSZZ][kp1] - Tau[SGSZZ][km1])*facZ);
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
