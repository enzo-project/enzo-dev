/***********************************************************************
/
/  GRID CLASS (COMPUTE CURL AND DIVERGENCE OF A VECTOR FIELD)
/
/  written by: Fen Zhao
/  date:       Sometime between 2005 and 2009
/  modified1:  Matthew Turk, November 2009
/
/  PURPOSE:
/
/   This computes the gradient and divergence of vector fields.
/   NOTE that while this ALLOCATES memory, it is up to the calling
/   function to de-allocate.
/
************************************************************************/
 
#include <hdf5.h>
#include <string.h>
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

int grid::ComputeVectorAnalysisFields(
        field_type fx, field_type fy, field_type fz,
        float* &curl_x, float* &curl_y, float* &curl_z,
        float* &div)
{
  int F1Num, F2Num, F3Num;
  F1Num=FindField(fx, FieldType, NumberOfBaryonFields);
  F2Num=FindField(fy, FieldType, NumberOfBaryonFields);
  F3Num=FindField(fz, FieldType, NumberOfBaryonFields);

  int size = 1;
  for (int dim = 0; dim < GridRank; dim++) size *= GridDimension[dim];

  if(GridRank==3){
    curl_x = new float [size];
    curl_y = new float [size];}

    curl_z = new float [size];
    div   = new float [size];

  FLOAT dx = CellWidth[0][0],
        dy = CellWidth[1][0], dz;
  if (GridRank>2)
    dz = CellWidth[2][0];

  /* Copy active part of field into grid */
  int igrid, igridyp1, igridym1, igridzp1, igridzm1;
  for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {

        igrid = (k*GridDimension[1] + j)*GridDimension[0] + i;
        igridyp1 = (k*GridDimension[1] + j + 1)*GridDimension[0] + i;
        igridym1 = (k*GridDimension[1] + j - 1)*GridDimension[0] + i;
        igridzp1 = ((k+1)*GridDimension[1]+j)*GridDimension[0] + i;
        igridzm1 = ((k-1)*GridDimension[1]+j)*GridDimension[0] + i;

        if (GridRank==3){

          div[igrid ] =
            float(
                (0.5*(BaryonField[F1Num][igrid+1]-BaryonField[F1Num][igrid-1])/dx +
                 0.5*(BaryonField[F2Num][igridyp1]-BaryonField[F2Num][igridym1])/dy +
                 0.5*(BaryonField[F3Num][igridzp1]-BaryonField[F3Num][igridzm1])/dz)
                );
          curl_x[igrid ] =
            float(
                (0.5*(BaryonField[F3Num][igridyp1]-BaryonField[F3Num][igridym1])/dy -		      
                 0.5*(BaryonField[F2Num][igridzp1]-BaryonField[F2Num][igridzm1])/dz)
                );
          curl_y[igrid ] =
            float(
                (0.5*(BaryonField[F1Num][igridzp1]-BaryonField[F1Num][igridzm1])/dz -		      
                 0.5*(BaryonField[F3Num][igrid+1]-BaryonField[F3Num][igrid-1])/dx)
                );
          curl_z[igrid ] =
            float(
                (0.5*(BaryonField[F2Num][igrid+1]-BaryonField[F2Num][igrid-1])/dx -		      
                 0.5*(BaryonField[F1Num][igridyp1]-BaryonField[F1Num][igridym1])/dy)
                );
        }

        if (GridRank==2){

          div[igrid ] =
            float(
                (0.5*(BaryonField[F1Num][igrid+1]-BaryonField[F1Num][igrid-1])/dx +
                 0.5*(BaryonField[F2Num][igridyp1]-BaryonField[F2Num][igridym1])/dy 
                ));


          curl_z[igrid] =
            float(
                (0.5*(BaryonField[F2Num][igrid+1]-BaryonField[F2Num][igrid-1])/dx -		      
                 0.5*(BaryonField[F1Num][igridyp1]-BaryonField[F1Num][igridym1])/dy)
                );
        }
      }
    }
  }

}
