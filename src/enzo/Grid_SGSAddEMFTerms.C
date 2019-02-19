/***********************************************************************
 *
 * INFORMATION This file is part of a subgrid-scale (SGS) modeling 
 * framework in order to conduct explicit large eddy simulations (LES).
 *
 * The functions in this file concern models for the turbulent 
 * electromotive force (EMF) in the induction equation.
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
 * This function adds to the EMF 
 * the pure (unscaled) full (compressible) nonlinear model:
 * EMF = 1/12 * Delta^2 * eps_ijk * (u_j,l * B_k,l - (ln rho),l u_j,l B_k)
 *
 * See equation (37) in Grete2016a for details (such as coefficient values)
 * or Vlaykov et al 2016 Phys. Plasmas 23 062316 doi: 10.1063/1.4954303 for
 * the derivation.
 */
void grid::SGS_AddEMF_nonlinear_compressive(float **EMF) {
    if (debug1)
        printf("[%"ISYM"] grid::SGS_AddEMF_nonlinear_compressive start\n",MyProcessorNumber);

    int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
    int B1Num, B2Num, B3Num, PhiNum;
    this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num,
            TENum, B1Num, B2Num, B3Num, PhiNum);

    float *rho, *Bx, *By, *Bz;

    // if an explicit filter should be used
    // (at this point the fields are already filtered, 
    // see hydro_rk/Grid_MHDSourceTerms.C)
    if (SGSFilterWidth > 1.) {
        rho = FilteredFields[0];
        Bx  = FilteredFields[4];
        By  = FilteredFields[5];
        Bz  = FilteredFields[6];
    // if the model should be calculated based on grid-scale quantities
    // (not recommended, see Grete2017)
    } else {
        rho = BaryonField[DensNum];
        Bx  = BaryonField[B1Num];
        By  = BaryonField[B2Num];
        Bz  = BaryonField[B3Num];
    }

    int size = 1;
    int StartIndex[MAX_DIMENSION];
    int EndIndex[MAX_DIMENSION];

    for (int dim = 0; dim < MAX_DIMENSION; dim++) {
        size *= GridDimension[dim];

        /* we need the EMF in the first ghost zone as well
         * as we'll take another derivative later on */
        StartIndex[dim] = GridStartIndex[dim] - 1;
        EndIndex[dim] = GridEndIndex[dim] + 1;
    }


    // the combined prefactor
    float CDeltaSqr = 1./12. * SGScoeffNLemfCompr * POW(SGSFilterWidth,2.) *
        POW(CellWidth[0][0]*CellWidth[1][0]*CellWidth[2][0],2./3.);

    int igrid, ip1, im1, jp1, jm1, kp1, km1;
    float facX = 1. / (2. * CellWidth[0][0]);
    float facY = 1. / (2. * CellWidth[1][0]);
    float facZ = 1. / (2. * CellWidth[2][0]);
    // partial derivatives of ln(rho)
    float lnRhod[MAX_DIMENSION];

    for (int k = StartIndex[2]; k <= EndIndex[2]; k++)
        for (int j = StartIndex[1]; j <= EndIndex[1]; j++)
            for (int i = StartIndex[0]; i <= EndIndex[0]; i++) {

                igrid = i + (j+k*GridDimension[1])*GridDimension[0];
                ip1 = i+1 + (j+k*GridDimension[1])*GridDimension[0];
                im1 = i-1 + (j+k*GridDimension[1])*GridDimension[0];
                jp1 = i + (j+1+k*GridDimension[1])*GridDimension[0];
                jm1 = i + (j-1+k*GridDimension[1])*GridDimension[0];
                kp1 = i + (j+(k+1)*GridDimension[1])*GridDimension[0];
                km1 = i + (j+(k-1)*GridDimension[1])*GridDimension[0];

                lnRhod[SGSX] = (rho[ip1] - rho[im1]) * facX / rho[igrid];
                lnRhod[SGSY] = (rho[jp1] - rho[jm1]) * facY / rho[igrid];
                lnRhod[SGSZ] = (rho[kp1] - rho[km1]) * facZ / rho[igrid];

                for (int l = 0; l < MAX_DIMENSION; l++) {
                    EMF[SGSX][igrid] += CDeltaSqr * (
                            JacVel[SGSY][l][igrid] * JacB[SGSZ][l][igrid] 
                            - JacVel[SGSZ][l][igrid] * JacB[SGSY][l][igrid] 
                            - lnRhod[l] * JacVel[SGSY][l][igrid] * Bz[igrid] 
                            + lnRhod[l] * JacVel[SGSZ][l][igrid] * By[igrid]);
                    EMF[SGSY][igrid] += CDeltaSqr * (
                            JacVel[SGSZ][l][igrid] * JacB[SGSX][l][igrid] 
                            - JacVel[SGSX][l][igrid] * JacB[SGSZ][l][igrid] 
                            - lnRhod[l] * JacVel[SGSZ][l][igrid] * Bx[igrid] 
                            + lnRhod[l] * JacVel[SGSX][l][igrid] * Bz[igrid]);
                    EMF[SGSZ][igrid] += CDeltaSqr * (
                            JacVel[SGSX][l][igrid] * JacB[SGSY][l][igrid] 
                            - JacVel[SGSY][l][igrid] * JacB[SGSX][l][igrid] 
                            - lnRhod[l] * JacVel[SGSX][l][igrid] * By[igrid] 
                            + lnRhod[l] * JacVel[SGSY][l][igrid] * Bx[igrid]);
                }

            }

}

/* 
 * This function adds to the EMF 
 * an eddy resistivity model where the strength of anomalous resistivity is 
 * scaled by the total SGS energy as given by a model based on 
 * realizability conditions:
 * EMF = -C * Delta^2 * sqrt(|S*|^2 + |M|^2/rho) * J
 *
 * See equation (23) and (13) in Grete2016a for details (such as coefficient values)
 */
void grid::SGS_AddEMF_eddy_resistivity(float **EMF) {
    if (debug1)
        printf("[%"ISYM"] grid::SGS_AddEMF_eddy_resistivity start\n",MyProcessorNumber);

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

        /* we need the EMF in the first ghost zone as well
         * as we'll take another derivative later on */
        StartIndex[dim] = GridStartIndex[dim] - 1;
        EndIndex[dim] = GridEndIndex[dim] + 1;
    }


    // the combined prefactor
    float MinusCDeltaSqr = -SGScoeffERS2M2Star * POW(SGSFilterWidth,2.) *
        POW(CellWidth[0][0]*CellWidth[1][0]*CellWidth[2][0],2./3.);

    int igrid;

    float traceSthird, traceMthird;
    float sqrtS2StarplusM2overRho;

    for (int k = StartIndex[2]; k <= EndIndex[2]; k++)
        for (int j = StartIndex[1]; j <= EndIndex[1]; j++)
            for (int i = StartIndex[0]; i <= EndIndex[0]; i++) {

                igrid = i + (j+k*GridDimension[1])*GridDimension[0];

                traceSthird = (JacVel[SGSX][SGSX][igrid] + JacVel[SGSY][SGSY][igrid] + JacVel[SGSZ][SGSZ][igrid])/3.;
                traceMthird = (JacB[SGSX][SGSX][igrid] + JacB[SGSY][SGSY][igrid] + JacB[SGSZ][SGSZ][igrid])/3.;

                sqrtS2StarplusM2overRho = POW(
                        2.*(POW(JacVel[SGSX][SGSX][igrid]-traceSthird,2.) +
                            POW(JacVel[SGSY][SGSY][igrid]-traceSthird,2.) +
                            POW(JacVel[SGSZ][SGSZ][igrid]-traceSthird,2.)
                           )
                        + POW(JacVel[SGSX][SGSY][igrid] + JacVel[SGSY][SGSX][igrid],2.)
                        + POW(JacVel[SGSY][SGSZ][igrid] + JacVel[SGSZ][SGSY][igrid],2.)
                        + POW(JacVel[SGSX][SGSZ][igrid] + JacVel[SGSZ][SGSX][igrid],2.)
                        + (2.*(POW(JacB[SGSX][SGSX][igrid]-traceMthird,2.) +
                                POW(JacB[SGSY][SGSY][igrid]-traceMthird,2.) +
                                POW(JacB[SGSZ][SGSZ][igrid]-traceMthird,2.)
                              )
                            + POW(JacB[SGSX][SGSY][igrid] + JacB[SGSY][SGSX][igrid],2.)
                            + POW(JacB[SGSY][SGSZ][igrid] + JacB[SGSZ][SGSY][igrid],2.)
                            + POW(JacB[SGSX][SGSZ][igrid] + JacB[SGSZ][SGSX][igrid],2.)
                          )/rho[igrid],1./2.);

                EMF[SGSX][igrid] += MinusCDeltaSqr * sqrtS2StarplusM2overRho * 
                    (JacB[SGSZ][SGSY][igrid] - JacB[SGSY][SGSZ][igrid]);
                EMF[SGSY][igrid] += MinusCDeltaSqr * sqrtS2StarplusM2overRho * 
                    (JacB[SGSX][SGSZ][igrid] - JacB[SGSZ][SGSX][igrid]);
                EMF[SGSZ][igrid] += MinusCDeltaSqr * sqrtS2StarplusM2overRho * 
                    (JacB[SGSY][SGSX][igrid] - JacB[SGSX][SGSY][igrid]);

            }
}

/* 
 * This function adds to the EMF 
 * a scale-similarity motivated term:
 * EMF = flt(u x B) - flt(u) x flt(B)
 *
 * See equation (32) in Grete2016a for details (such as coefficient values)
 */
void grid::SGS_AddEMF_scale_similarity(float **EMF) {
    if (debug1)
        printf("[%"ISYM"] grid::SGS_AddEMF_scale_similarity start\n",MyProcessorNumber);

    int size = 1;
    int StartIndex[MAX_DIMENSION];
    int EndIndex[MAX_DIMENSION];

    for (int dim = 0; dim < MAX_DIMENSION; dim++) {
        size *= GridDimension[dim];

        /* we need the EMF in the first ghost zone as well
         * as we'll take another derivative later on */
        StartIndex[dim] = GridStartIndex[dim] - 1;
        EndIndex[dim] = GridEndIndex[dim] + 1;
    }


    int igrid;

    for (int k = StartIndex[2]; k <= EndIndex[2]; k++)
        for (int j = StartIndex[1]; j <= EndIndex[1]; j++)
            for (int i = StartIndex[0]; i <= EndIndex[0]; i++) {

                igrid = i + (j+k*GridDimension[1])*GridDimension[0];

                EMF[SGSX][igrid] += SGScoeffSSemf * FltUB[SGSX][igrid] - (
                    FilteredFields[2][igrid] * FilteredFields[6][igrid] - 
                    FilteredFields[3][igrid] * FilteredFields[5][igrid]);
                EMF[SGSY][igrid] += SGScoeffSSemf * FltUB[SGSY][igrid] - (
                    FilteredFields[3][igrid] * FilteredFields[4][igrid] - 
                    FilteredFields[1][igrid] * FilteredFields[6][igrid]);
                EMF[SGSZ][igrid] += SGScoeffSSemf * FltUB[SGSZ][igrid] - (
                    FilteredFields[1][igrid] * FilteredFields[5][igrid] - 
                    FilteredFields[2][igrid] * FilteredFields[4][igrid]);
            }
}


/*
 * This function initializes a zero EMF and calls the individual
 * functions that add the different terms to the EMF.
 * Finally, the curl of the EMF is added to the dU vector used by
 * the MUSCL framework in hydro_rk/Grid_MHDSourceTerms.C 
 */
int grid::SGS_AddEMFTerms(float **dU) {
    if (ProcessorNumber != MyProcessorNumber) {
        return SUCCESS;
    }

    if (Time == 0.)
        return SUCCESS;

    if (debug1)
        printf("[%"ISYM"] grid::SGS_AddEMFTerms start\n",MyProcessorNumber);

    int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
    int B1Num, B2Num, B3Num, PhiNum;
    this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num,
            TENum, B1Num, B2Num, B3Num, PhiNum);

    int size = 1;
    float *EMF[MAX_DIMENSION];

    for (int dim = 0; dim < MAX_DIMENSION; dim++) {
        size *= GridDimension[dim];
    }

    for (int dim = 0; dim < MAX_DIMENSION; dim++) {
        EMF[dim] = new float[size];
        for (int i = 0; i < size; i++)
            EMF[dim][i] = 0.;
    }



    // the individual terms are added/activated by a non-zero coefficient
    if (SGScoeffERS2M2Star != 0.) 
        SGS_AddEMF_eddy_resistivity(EMF);

    if (SGScoeffNLemfCompr != 0.) 
        SGS_AddEMF_nonlinear_compressive(EMF);

    if (SGScoeffSSemf != 0.) 
        SGS_AddEMF_scale_similarity(EMF);

    int n = 0;
    int igrid, ip1, im1, jp1, jm1, kp1, km1;
    float BxIncr,ByIncr,BzIncr,EtotIncr;

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


                BxIncr = dtFixed * ((EMF[2][jp1]-EMF[2][jm1])*facY - (EMF[1][kp1]-EMF[1][km1])*facZ);
                EtotIncr = BaryonField[B1Num][igrid] * BxIncr + 0.5 * BxIncr * BxIncr;

                ByIncr = dtFixed * ((EMF[0][kp1]-EMF[0][km1])*facZ - (EMF[2][ip1]-EMF[2][im1])*facX);
                EtotIncr += BaryonField[B2Num][igrid] * ByIncr + 0.5 * ByIncr * ByIncr;

                BzIncr = dtFixed * ((EMF[1][ip1]-EMF[1][im1])*facX - (EMF[0][jp1]-EMF[0][jm1])*facY);
                EtotIncr += BaryonField[B3Num][igrid] * BzIncr + 0.5 * BzIncr * BzIncr;

                dU[iBx][n] += BxIncr;
                dU[iBy][n] += ByIncr;
                dU[iBz][n] += BzIncr;
                dU[iEtot][n] += EtotIncr;
            }

    if (debug1)
        printf("[%"ISYM"] grid::SGS_AddEMFTerms end, last incr: %"FSYM" %"FSYM" %"FSYM" %"FSYM"\n",
                MyProcessorNumber,BxIncr,ByIncr,BzIncr,EtotIncr);

    for (int dim = 0; dim < MAX_DIMENSION; dim++) {
        delete [] EMF[dim];
    }

    return SUCCESS;
}
