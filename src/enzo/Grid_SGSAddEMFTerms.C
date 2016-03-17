#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
/* pure (unscaled) full (compressible) nonlinear model 
 * EMF = 1/12 * Delta^2 * eps_ijk * (u_j,l * B_k,l - (ln rho),l u_j,l B_k)
 * see eq TODO of TODO for details
 */
void grid::SGSAddEMFNLemfComprTerm(float **EMF) {
    if (debug)
        printf("[%"ISYM"] grid::SGSAddEMFNLemfComprTerm start\n",MyProcessorNumber);

    int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
    int B1Num, B2Num, B3Num, PhiNum;
    this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num,
            TENum, B1Num, B2Num, B3Num, PhiNum);

    float *rho, *Bx, *By, *Bz;
    if (SGSFilterWidth > 1.) {
        rho = FilteredFields[0];
        Bx  = FilteredFields[4];
        By  = FilteredFields[5];
        Bz  = FilteredFields[6];
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


    float CDeltaSqr = 1./12. * SGScoeffNLemfCompr * pow(SGSFilterWidth,2.) *
        pow(CellWidth[0][0]*CellWidth[1][0]*CellWidth[2][0],2./3.);

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

                lnRhod[X] = (rho[ip1] - rho[im1]) * facX / rho[igrid];
                lnRhod[Y] = (rho[jp1] - rho[jm1]) * facY / rho[igrid];
                lnRhod[Z] = (rho[kp1] - rho[km1]) * facZ / rho[igrid];

                for (int l = 0; l < MAX_DIMENSION; l++) {
                    EMF[X][igrid] += CDeltaSqr * (
                            JacVel[Y][l][igrid] * JacB[Z][l][igrid] 
                            - JacVel[Z][l][igrid] * JacB[Y][l][igrid] 
                            - lnRhod[l] * JacVel[Y][l][igrid] * Bz[igrid] 
                            + lnRhod[l] * JacVel[Z][l][igrid] * By[igrid]);
                    EMF[Y][igrid] += CDeltaSqr * (
                            JacVel[Z][l][igrid] * JacB[X][l][igrid] 
                            - JacVel[X][l][igrid] * JacB[Z][l][igrid] 
                            - lnRhod[l] * JacVel[Z][l][igrid] * Bx[igrid] 
                            + lnRhod[l] * JacVel[X][l][igrid] * Bz[igrid]);
                    EMF[Z][igrid] += CDeltaSqr * (
                            JacVel[X][l][igrid] * JacB[Y][l][igrid] 
                            - JacVel[Y][l][igrid] * JacB[X][l][igrid] 
                            - lnRhod[l] * JacVel[X][l][igrid] * By[igrid] 
                            + lnRhod[l] * JacVel[Y][l][igrid] * Bx[igrid]);
                }

            }

}

/* eddy resistivity model scaled by Smagorinsky energies
 * EMF = -C * Delta^2 * sqrt(|S|^2 + |J|^2/rho) * J
 * see eq TODO of TODO for details
 */
void grid::SGSAddEMFERS2J2Term(float **EMF) {
    if (debug)
        printf("[%"ISYM"] grid::SGSAddEMFERS2J2Term start\n",MyProcessorNumber);

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

        /* we need the EMF in the first ghost zone as well
         * as we'll take another derivative later on */
        StartIndex[dim] = GridStartIndex[dim] - 1;
        EndIndex[dim] = GridEndIndex[dim] + 1;
    }


    float MinusCDeltaSqr = -SGScoeffERS2J2 * pow(SGSFilterWidth,2.) * 
        pow(CellWidth[0][0]*CellWidth[1][0]*CellWidth[2][0],2./3.);

    int igrid;
    float sqrtS2plusJ2overRho;

    for (int k = StartIndex[2]; k <= EndIndex[2]; k++)
        for (int j = StartIndex[1]; j <= EndIndex[1]; j++)
            for (int i = StartIndex[0]; i <= EndIndex[0]; i++) {

                igrid = i + (j+k*GridDimension[1])*GridDimension[0];

                sqrtS2plusJ2overRho = pow(
                        2.*(pow(JacVel[X][X][igrid],2.) +
                            pow(JacVel[Y][Y][igrid],2.) +
                            pow(JacVel[Z][Z][igrid],2.)
                           )
                        + pow(JacVel[X][Y][igrid] + JacVel[Y][X][igrid],2.)
                        + pow(JacVel[Y][Z][igrid] + JacVel[Z][Y][igrid],2.)
                        + pow(JacVel[X][Z][igrid] + JacVel[Z][X][igrid],2.)
                        + (pow(JacB[Z][Y][igrid] - JacB[Y][Z][igrid],2.) +
                            pow(JacB[X][Z][igrid] - JacB[Z][X][igrid],2.) +
                            pow(JacB[Y][X][igrid] - JacB[X][Y][igrid],2.)
                          )/rho[igrid],1./2.);

                EMF[X][igrid] += MinusCDeltaSqr * sqrtS2plusJ2overRho * 
                    (JacB[Z][Y][igrid] - JacB[Y][Z][igrid]);
                EMF[Y][igrid] += MinusCDeltaSqr * sqrtS2plusJ2overRho * 
                    (JacB[X][Z][igrid] - JacB[Z][X][igrid]);
                EMF[Z][igrid] += MinusCDeltaSqr * sqrtS2plusJ2overRho * 
                    (JacB[Y][X][igrid] - JacB[X][Y][igrid]);

            }

}

/* eddy resistivity model scaled by realiz. energies
 * EMF = -C * Delta^2 * sqrt(|S*|^2 + |M|^2/rho) * J
 * see eq TODO of TODO for details
 */
void grid::SGSAddEMFERS2M2StarTerm(float **EMF) {
    if (debug)
        printf("[%"ISYM"] grid::SGSAddEMFERS2M2StarTerm start\n",MyProcessorNumber);

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

        /* we need the EMF in the first ghost zone as well
         * as we'll take another derivative later on */
        StartIndex[dim] = GridStartIndex[dim] - 1;
        EndIndex[dim] = GridEndIndex[dim] + 1;
    }


    float MinusCDeltaSqr = -SGScoeffERS2M2Star * pow(SGSFilterWidth,2.) *
        pow(CellWidth[0][0]*CellWidth[1][0]*CellWidth[2][0],2./3.);

    int igrid;

    /* magic with S |S| S*... could potentially handled by
     * external function, should reduce CPU time, but increase memory usage
     */
    float traceSthird, traceMthird;
    float sqrtS2StarplusM2overRho;
    /* just for fun: how accurate is Dedner
     * we count the number of cells where divB is dynamically important
     * divB * Delta / |B| > 1.
     */
    int divBerror = 0;

    for (int k = StartIndex[2]; k <= EndIndex[2]; k++)
        for (int j = StartIndex[1]; j <= EndIndex[1]; j++)
            for (int i = StartIndex[0]; i <= EndIndex[0]; i++) {

                igrid = i + (j+k*GridDimension[1])*GridDimension[0];

                traceSthird = (JacVel[X][X][igrid] + JacVel[Y][Y][igrid] + JacVel[Z][Z][igrid])/3.;
                traceMthird = (JacB[X][X][igrid] + JacB[Y][Y][igrid] + JacB[Z][Z][igrid])/3.;

                if (debug && (traceMthird*pow(CellWidth[0][0]*CellWidth[1][0]*CellWidth[2][0],2./3.)*3./pow(
                                BaryonField[B1Num][igrid]*BaryonField[B1Num][igrid] +
                                BaryonField[B2Num][igrid]*BaryonField[B2Num][igrid] +
                                BaryonField[B3Num][igrid]*BaryonField[B3Num][igrid],1./2.) > 1.)) {
                    divBerror++;
                }


                sqrtS2StarplusM2overRho = pow(
                        2.*(pow(JacVel[X][X][igrid]-traceSthird,2.) +
                            pow(JacVel[Y][Y][igrid]-traceSthird,2.) +
                            pow(JacVel[Z][Z][igrid]-traceSthird,2.)
                           )
                        + pow(JacVel[X][Y][igrid] + JacVel[Y][X][igrid],2.)
                        + pow(JacVel[Y][Z][igrid] + JacVel[Z][Y][igrid],2.)
                        + pow(JacVel[X][Z][igrid] + JacVel[Z][X][igrid],2.)
                        + (2.*(pow(JacB[X][X][igrid]-traceMthird,2.) +
                                pow(JacB[Y][Y][igrid]-traceMthird,2.) +
                                pow(JacB[Z][Z][igrid]-traceMthird,2.)
                              )
                            + pow(JacB[X][Y][igrid] + JacB[Y][X][igrid],2.)
                            + pow(JacB[Y][Z][igrid] + JacB[Z][Y][igrid],2.)
                            + pow(JacB[X][Z][igrid] + JacB[Z][X][igrid],2.)
                          )/rho[igrid],1./2.);

                EMF[X][igrid] += MinusCDeltaSqr * sqrtS2StarplusM2overRho * 
                    (JacB[Z][Y][igrid] - JacB[Y][Z][igrid]);
                EMF[Y][igrid] += MinusCDeltaSqr * sqrtS2StarplusM2overRho * 
                    (JacB[X][Z][igrid] - JacB[Z][X][igrid]);
                EMF[Z][igrid] += MinusCDeltaSqr * sqrtS2StarplusM2overRho * 
                    (JacB[Y][X][igrid] - JacB[X][Y][igrid]);

            }
    if (debug)
        printf("[%"ISYM"] grid::SGSAddEMFERS2M2StarTerm divB error total: %"ISYM" |\%: %"FSYM"\n",
                MyProcessorNumber,divBerror,
                (float) divBerror / (float)((EndIndex[0] + 1 - StartIndex[0])*(EndIndex[1] + 1 - StartIndex[1])*(EndIndex[2] + 1 - StartIndex[2])));

}

/* scale-similarity model 
 * EMF = flt(u x B) - flt(u) x flt(B)
 * see eq TODO of TODO for details
 */
void grid::SGSAddEMFSSTerm(float **EMF) {
    if (debug)
        printf("[%"ISYM"] grid::SGSAddEMFSSTerm start\n",MyProcessorNumber);

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

                EMF[X][igrid] += SGScoeffSSemf * FltUB[X][igrid] - (
                    FilteredFields[2][igrid] * FilteredFields[6][igrid] - 
                    FilteredFields[3][igrid] * FilteredFields[5][igrid]);
                EMF[Y][igrid] += SGScoeffSSemf * FltUB[Y][igrid] - (
                    FilteredFields[3][igrid] * FilteredFields[4][igrid] - 
                    FilteredFields[1][igrid] * FilteredFields[6][igrid]);
                EMF[Z][igrid] += SGScoeffSSemf * FltUB[Z][igrid] - (
                    FilteredFields[1][igrid] * FilteredFields[5][igrid] - 
                    FilteredFields[2][igrid] * FilteredFields[4][igrid]);

            }

}

int grid::SGSAddEMFTerms(float **dU) {
    if (ProcessorNumber != MyProcessorNumber) {
        return SUCCESS;
    }

    if (Time == 0.)
        return SUCCESS;

    if (debug)
        printf("[%"ISYM"] grid::SGSAddEMFTerms start\n",MyProcessorNumber);

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


    if (SGScoeffERS2J2 != 0.) 
        SGSAddEMFERS2J2Term(EMF);

    if (SGScoeffERS2M2Star != 0.) 
        SGSAddEMFERS2M2StarTerm(EMF);

    if (SGScoeffNLemfCompr != 0.) 
        SGSAddEMFNLemfComprTerm(EMF);

    if (SGScoeffSSemf != 0.) 
        SGSAddEMFSSTerm(EMF);

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

    if (debug)
        printf("[%"ISYM"] grid::SGSAddEMFTerms end, last incr: %"FSYM" %"FSYM" %"FSYM" %"FSYM"\n",
                MyProcessorNumber,BxIncr,ByIncr,BzIncr,EtotIncr);

    for (int dim = 0; dim < MAX_DIMENSION; dim++) {
        delete [] EMF[dim];
    }

    return SUCCESS;
}
