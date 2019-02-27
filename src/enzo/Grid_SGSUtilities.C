/***********************************************************************
 *
 * INFORMATION This file is part of a subgrid-scale (SGS) modeling 
 * framework in order to conduct explicit large eddy simulations (LES).
 *
 * The functions in this file are considered utility functions and
 * concern the calculation of Jacobians and explicit filtering operations.
 *
 * WRITTEN BY Philipp Grete (mail@pgrete.de)
 *
 * DATE 2016
 *
************************************************************************/

#include "preincludes.h"
#include <stdlib.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

/*
 * This function conducts an explicit filtering operation on the
 * primary quantities. 
 * The original fields remain untouched and the result is stored in the
 * FilteredFields array.
 * Mass-weighted filtering applies to the velocity field.
 * The multi-dimensional filtering operation is constructed by sequential
 * application of the one-dimensional filter.
 *
 * At this point, this function is very flexible with respect to the 
 * discrete filter definition (both in effective width and weights).
 * However, this also makes the nested for-loops less efficient, but the
 * overhead is reasonable (see e.g. Grete et al 2017 Phys. Rev. E. 95 033206)
 * 
 * For additional information on how to construct multi-dimensional
 * discrete filters, see e.g. 
 * Vasilyev, Lund & Moin (1998) Journal of Comp. Physics 146, 82 or
 * Sagaut and Grohens (1999) Journal for Num. Meth. in Fluids 31, 1195
 * 
 */
int grid::SGSUtil_FilterFields() {
    if (ProcessorNumber != MyProcessorNumber) {
        return SUCCESS;
    }

    int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, 
        B1Num, B2Num, B3Num, PhiNum;
    this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, 
            TENum, B1Num, B2Num, B3Num, PhiNum);

    int size = 1;
    int StartIndex[MAX_DIMENSION];
    int EndIndex[MAX_DIMENSION];

    for (int dim = 0; dim < MAX_DIMENSION; dim++) {
        size *= GridDimension[dim];

        /* we need the filtered fields in the second ghost zone as well
         * as we need derivatives in first ghost zone */
        StartIndex[dim] = GridStartIndex[dim] - 2;
        EndIndex[dim] = GridEndIndex[dim] + 2;
    }

    int NumFilteredFields;
    if (UseMHD)
      NumFilteredFields = 7;
    else
      NumFilteredFields = 4;

    for (int m = 0; m < NumFilteredFields; m++)
        if (FilteredFields[m] == NULL) {
            FilteredFields[m] = new float[size];
            for (int o = 0; o < size; o++)
                FilteredFields[m][o] = 0.;
        }

    int N = SGSFilterStencil/2;
    int igrid, ifilter;
    float totalWeight;

    for (int k = StartIndex[2]; k <= EndIndex[2]; k++) {
      for (int j = StartIndex[1]; j <= EndIndex[1]; j++) {
        for (int i = StartIndex[0]; i <= EndIndex[0]; i++) {
        
          igrid = i + (j+k*GridDimension[1])*GridDimension[0];

          for (int l = 0; l < NumFilteredFields; l++)
            FilteredFields[l][igrid] = 0.;

          for (int  l = -N; l <= N; l++) {
            for (int  m = -N; m <= N; m++) {
              for (int  n = -N; n <= N; n++) {


                ifilter = i + l + (j + m + (k+n)*GridDimension[1])*GridDimension[0];
                totalWeight =  SGSFilterWeights[ABS(l)] * SGSFilterWeights[ABS(m)] * SGSFilterWeights[ABS(n)];
            
                // rho
                FilteredFields[0][igrid] += totalWeight * BaryonField[DensNum][ifilter];
            
                // prepare mass weighted velocity fields
                FilteredFields[1][igrid] += totalWeight * BaryonField[DensNum][ifilter]*BaryonField[Vel1Num][ifilter];
                FilteredFields[2][igrid] += totalWeight * BaryonField[DensNum][ifilter]*BaryonField[Vel2Num][ifilter];
                FilteredFields[3][igrid] += totalWeight * BaryonField[DensNum][ifilter]*BaryonField[Vel3Num][ifilter];
            
                // magnetic fields
                if (UseMHD) {
                  FilteredFields[4][igrid] += totalWeight * BaryonField[B1Num][ifilter];
                  FilteredFields[5][igrid] += totalWeight * BaryonField[B2Num][ifilter];
                  FilteredFields[6][igrid] += totalWeight * BaryonField[B3Num][ifilter];
                }
              }
            }
          } // end of innter triple for

          // now that the density is filtered, we can finalize mass-weighted filtering
          FilteredFields[1][igrid] /= FilteredFields[0][igrid];
          FilteredFields[2][igrid] /= FilteredFields[0][igrid];
          FilteredFields[3][igrid] /= FilteredFields[0][igrid];
        }
      }
    } // end of outer triple for

    return SUCCESS;
}

/*
 * This functional calculated the Jacobian of an arbitrary 3-dimensional
 * field (components given by field1, field2, and field3).
 * The result is stored in the Jac array.
 */
int grid::SGSUtil_ComputeJacobian(float *Jac[][MAX_DIMENSION],float *field1,float* field2,float* field3) {
    if (ProcessorNumber != MyProcessorNumber) {
        return SUCCESS;
    }

    int size = 1;
    int StartIndex[MAX_DIMENSION];
    int EndIndex[MAX_DIMENSION];

    for (int dim = 0; dim < MAX_DIMENSION; dim++) {
        size *= GridDimension[dim];

        /* we need the Jacobians in the first ghost zone as well
         * as we'll take another derivative later on */
        StartIndex[dim] = GridStartIndex[dim] - 1;
        EndIndex[dim] = GridEndIndex[dim] + 1;
    }

    for (int m = 0; m < MAX_DIMENSION; m++) 
        for (int n = 0; n < MAX_DIMENSION; n++) {
            if (Jac[m][n] == NULL) {
                Jac[m][n] = new float[size];
                for (int o = 0; o < size; o++)
                    Jac[m][n][o] = 0.;
            }
        }


    int igrid, ip1, im1, jp1, jm1, kp1, km1;
    float facX = 1. / (2. * CellWidth[0][0]);
    float facY = 1. / (2. * CellWidth[1][0]);
    float facZ = 1. / (2. * CellWidth[2][0]);

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

                // xdx
                Jac[SGSX][SGSX][igrid] = (field1[ip1] - field1[im1]) * facX;
                // xdy
                Jac[SGSX][SGSY][igrid] = (field1[jp1] - field1[jm1]) * facY;
                // xdz
                Jac[SGSX][SGSZ][igrid] = (field1[kp1] - field1[km1]) * facZ;

                // ydx
                Jac[SGSY][SGSX][igrid] = (field2[ip1] - field2[im1]) * facX;
                // ydy
                Jac[SGSY][SGSY][igrid] = (field2[jp1] - field2[jm1]) * facY;
                // ydz
                Jac[SGSY][SGSZ][igrid] = (field2[kp1] - field2[km1]) * facZ;

                // zdx
                Jac[SGSZ][SGSX][igrid] = (field3[ip1] - field3[im1]) * facX;
                // zdy
                Jac[SGSZ][SGSY][igrid] = (field3[jp1] - field3[jm1]) * facY;
                // zdz
                Jac[SGSZ][SGSZ][igrid] = (field3[kp1] - field3[km1]) * facZ;

            }


    return SUCCESS;
}
    
/*
 * This function conducts an explicit filter operation on mixed quantities, e.g.
 * flt(rho u_i u_j), which are required by the scale-similarity SGS model.
 * The same general comments as for SGSUtil_FilterFields (see above) apply.
 */
int grid::SGSUtil_ComputeMixedFilteredQuantities() {

    if (ProcessorNumber != MyProcessorNumber) {
        return SUCCESS;
    }
    
    int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, 
        B1Num, B2Num, B3Num, PhiNum;
    this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, 
            TENum, B1Num, B2Num, B3Num, PhiNum);

    int size = 1;
    int StartIndex[MAX_DIMENSION];
    int EndIndex[MAX_DIMENSION];

    for (int dim = 0; dim < MAX_DIMENSION; dim++) {
        size *= GridDimension[dim];

        /* we need the mixed terms in the first ghost zone as well
         * as we'll take another derivative later on */
        StartIndex[dim] = GridStartIndex[dim] - 1;
        EndIndex[dim] = GridEndIndex[dim] + 1;
    }

    for (int m = 0; m < 6; m++) {
        if (FltrhoUU[m] == NULL) {
            FltrhoUU[m] = new float[size];
            for (int o = 0; o < size; o++)
                FltrhoUU[m][o] = 0.;
        }
    }
    
    if (UseMHD) {
      for (int m = 0; m < 6; m++) {
          if (FltBB[m] == NULL) {
              FltBB[m] = new float[size];
              for (int o = 0; o < size; o++)
                  FltBB[m][o] = 0.;
          }
      }
      for (int m = 0; m < 3; m++) {
          if (FltUB[m] == NULL) {
              FltUB[m] = new float[size];
              for (int o = 0; o < size; o++)
                  FltUB[m][o] = 0.;
          }
      }
    }

    int N = SGSFilterStencil/2;
    int igrid, ifilter;
    float totalWeight;

    /* this is !highly! inefficient, just making sure it's working */
    for (int k = StartIndex[2]; k <= EndIndex[2]; k++) {
      for (int j = StartIndex[1]; j <= EndIndex[1]; j++) {
        for (int i = StartIndex[0]; i <= EndIndex[0]; i++) {
        
          igrid = i + (j+k*GridDimension[1])*GridDimension[0];

          for (int l = 0; l < 6; l++) {
            FltrhoUU[l][igrid] = 0.;
            FltBB[l][igrid] = 0.;
          }
          for (int l = 0; l < 3; l++)
            FltUB[l][igrid] = 0.;

          for (int  l = -N; l <= N; l++) {
            for (int  m = -N; m <= N; m++) {
              for (int  n = -N; n <= N; n++) {


                ifilter = i + l + (j + m + (k+n)*GridDimension[1])*GridDimension[0];
                  totalWeight =  SGSFilterWeights[ABS(l)] * SGSFilterWeights[ABS(m)] * SGSFilterWeights[ABS(n)];
            
                FltrhoUU[SGSXX][igrid] += totalWeight * BaryonField[DensNum][ifilter] *
                  BaryonField[Vel1Num][ifilter] * BaryonField[Vel1Num][ifilter];
                FltrhoUU[SGSYY][igrid] += totalWeight * BaryonField[DensNum][ifilter] *
                  BaryonField[Vel2Num][ifilter] * BaryonField[Vel2Num][ifilter];
                FltrhoUU[SGSZZ][igrid] += totalWeight * BaryonField[DensNum][ifilter] *
                  BaryonField[Vel3Num][ifilter] * BaryonField[Vel3Num][ifilter];
                FltrhoUU[SGSXY][igrid] += totalWeight * BaryonField[DensNum][ifilter] *
                  BaryonField[Vel1Num][ifilter] * BaryonField[Vel2Num][ifilter];
                FltrhoUU[SGSYZ][igrid] += totalWeight * BaryonField[DensNum][ifilter] *
                  BaryonField[Vel2Num][ifilter] * BaryonField[Vel3Num][ifilter];
                FltrhoUU[SGSXZ][igrid] += totalWeight * BaryonField[DensNum][ifilter] *
                  BaryonField[Vel1Num][ifilter] * BaryonField[Vel3Num][ifilter];

                if (!UseMHD)
                  continue;
            
                FltBB[SGSXX][igrid] += totalWeight * BaryonField[B1Num][ifilter] * BaryonField[B1Num][ifilter];
                FltBB[SGSYY][igrid] += totalWeight * BaryonField[B2Num][ifilter] * BaryonField[B2Num][ifilter];
                FltBB[SGSZZ][igrid] += totalWeight * BaryonField[B3Num][ifilter] * BaryonField[B3Num][ifilter];
                FltBB[SGSXY][igrid] += totalWeight * BaryonField[B1Num][ifilter] * BaryonField[B2Num][ifilter];
                FltBB[SGSYZ][igrid] += totalWeight * BaryonField[B2Num][ifilter] * BaryonField[B3Num][ifilter];
                FltBB[SGSXZ][igrid] += totalWeight * BaryonField[B1Num][ifilter] * BaryonField[B3Num][ifilter];
            
                FltUB[SGSX][igrid] += totalWeight * (
                  BaryonField[Vel2Num][ifilter] * BaryonField[B3Num][ifilter] -
                  BaryonField[Vel3Num][ifilter] * BaryonField[B2Num][ifilter]);
                FltUB[SGSY][igrid] += totalWeight * (
                  BaryonField[Vel3Num][ifilter] * BaryonField[B1Num][ifilter] -
                  BaryonField[Vel1Num][ifilter] * BaryonField[B3Num][ifilter]);
                FltUB[SGSZ][igrid] += totalWeight * (
                  BaryonField[Vel1Num][ifilter] * BaryonField[B2Num][ifilter] -
                  BaryonField[Vel2Num][ifilter] * BaryonField[B1Num][ifilter]);
              }
            }
          } // end of inner triple for
        }
      }
    } // end of outer triple for

    return SUCCESS;
}
