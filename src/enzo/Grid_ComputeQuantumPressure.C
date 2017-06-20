/***********************************************************************
/
/  GRID CLASS (WRAPPER FOR HydroSolver)
/
/  written by: Xinyu Li
/  date:       November, 2016
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
c
c  INPUTS:
c     d       - density field (includes boundary zones)
c     dx,y,z  - zone width arrays for each dimension
c     e       - total specific energy field
c     end     - array (of dimension 3) specifying the end of the active
c               region for reach dimension (zero based)
c     eta1    - (dual) selection parameter for gas energy (typically ~0.1)
c     eta2    - (dual) selection parameter for total energy (typically ~0.001)
c     ge      - gas energy (used when idual = 1)
c     gr_x,y,zacc - gravitational acceleration fields
c     gravity - flag indicating whether or not to use gravity field (1 = yes)
c     gridvel - bulk grid velocity (vector of length 3)
c     i,j,kn  - dimensions of field arrays
c     idiff   - diffusion flag (0 = off)
c     idual   - dual energy formalism flag (0 = off)
c     ipresfree - pressure free flag (0 = off, 1 = on, i.e. p=0)
c     nhy     - cycle number (for better operator splitting)
c     rank    - dimension of problem (not currently used)
c     start   - array (of dimension 3) specifying the start of the active
c               region fo reach dimension (zero based)
c     tmp     - temporary work space (30 * largest_slice)
c     u       - x-velocity field
c     v       - y-velocity field
c     w       - z-velocity field
c     bottom  - true (1) if this is the lowest level
c     minsupecoef - coefficient for minimum pressure support (0 - not used)
/
************************************************************************/

// Solve the hydro equations with the solver, saving the subgrid fluxes
//

#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#define IDX(a,b,c) ( ((c)*jn + (b))*in + (a) )


int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int FindField(int f, int farray[], int n);



int grid::ComputeQuantumPressure(float Time)
{
    if (ProcessorNumber != MyProcessorNumber)return SUCCESS;

  /*  Locals */

  int i, j, k;
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  int QNum;
  float *d, *qpres;

  int rank = GridRank;
  int in = GridDimension[0];
  int jn = GridDimension[1];
  int kn = GridDimension[2];

  /* Allocate temporary space */

  int size = in*jn*kn;
  float *logd = new float[size];

    /* Error Check */

  //for (i = 0; i < rank; i++)
  //  if (GridDimension[i] > MAX_ANY_SINGLE_DIRECTION) {
  //    ENZO_FAIL("Quantum Pressure : A grid dimension is too long (increase max_any_single_direction.)\n");
  //  }

  //fprintf(stderr, "laplacian coefficients %lf %f %f \n", lapcoef*1e6, LengthUnits,TimeUnits);


 /* If using comoving coordinates, multiply dx by a(n+1/2).
       In one fell swoop, this recasts the equations solved by solver
       in comoving form (except for the expansion terms which are taken
       care of elsewhere). */
    FLOAT a = 1, dadt;

  if (ComovingCoordinates){
      if (CosmologyComputeExpansionFactor(Time+0.5*dtFixed, &a, &dadt) 
    == FAIL) {
  ENZO_FAIL("Error in CosmologyComputeExpansionFactors.");
      }
     }

    /* Find fields: density, total energy, velocity1-3. */

  this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum );
  QNum = FindField(QPres, FieldType, NumberOfBaryonFields);

  
  d = BaryonField[DensNum];
  qpres = BaryonField[QNum];

  /* Compute log of density */
   for (i = 0; i < size; i++) {
          //printf("%"ISYM",%"GSYM"\n",DensNum,d[i]);

    if (d[i] < 0) {
      ENZO_FAIL("Compute Log Density: Negative Density! \n");
    } else{
      logd[i] = log(d[i]);
      //fprintf(stderr, "ComputeLogd: p,a,logd=%"GSYM",%"GSYM",%"GSYM",%"ISYM",%"ISYM",%"ISYM",\n", 
      //      qpres[i],a,d[i],in,jn,kn);
    }
  }

    /* compute Quantum Pressure Field */

    for (k = 0; k < kn; k++) {
        for (j = 0; j < jn; j++) {
            for (i = 0; i < in; i++){
          //4th order differentiation
          /*qpres[IDX(i,j,k)] = (-1./12.*(logd[IDX(max(i-2,0),j,k)]+logd[IDX(min(i+2,in-1),j,k)])
                    +4./3.*(logd[IDX(max(i-1,0),j,k)]+logd[IDX(min(i+1,in-1),j,k)])
                    -5./2.*logd[IDX(i,j,k)])/(CellWidth[0][i]*CellWidth[0][i])/2.;
          qpres[IDX(i,j,k)] = qpres[IDX(i,j,k)] 
                  + pow((1./12.*(logd[IDX(max(i-2,0),j,k)]-logd[IDX(min(i+2,in-1),j,k)])
                    -2./3.*(logd[IDX(max(i-1,0),j,k)]-logd[IDX(min(i+1,in-1),j,k)]))/CellWidth[0][i],2)/4.;*/
          //2nd order differentiation
          qpres[IDX(i,j,k)] =(logd[IDX(max(i-1,0),j,k)]+logd[IDX(min(i+1,in-1),j,k)]
                             -2.*logd[IDX(i,j,k)])/(CellWidth[0][i]*CellWidth[0][i])/2.;
          qpres[IDX(i,j,k)] = qpres[IDX(i,j,k)] 
                  + pow(-1./2.*(logd[IDX(max(i-1,0),j,k)]-logd[IDX(min(i+1,in-1),j,k)])/CellWidth[0][i],2)/4.;
          //qpres[IDX(i,j,k)] = (pow(d[IDX(max(i-1,0),j,k)],0.5)-2.0*pow(d[IDX(i,j,k)],0.5)+pow(d[IDX(min(i+1,in-1),j,k)],0.5))/pow(CellWidth[0][i],2);

          if (rank > 1) {
          // 4th order
          /*qpres[IDX(i,j,k)] = qpres[IDX(i,j,k)]
          				  +(-1./12.*(logd[IDX(i,max(j-2,0),k)]+logd[IDX(i,min(j+2,jn-1),k)])
          				  +4./3.*(logd[IDX(i,max(j-1,0),k)]+logd[IDX(i,min(j+1,jn-1),k)])
          				  -5./2.*logd[IDX(i,j,k)])/(CellWidth[1][j]*CellWidth[1][j])/2.;
          qpres[IDX(i,j,k)] = qpres[IDX(i,j,k)] 
          				+ pow((1./12.*(logd[IDX(i,max(j-2,0),k)]-logd[IDX(i,min(j+2,jn-1),k)])
          				  -2./3.*(logd[IDX(i,max(j-1,0),k)]-logd[IDX(i,min(j+1,jn-1),k)]))/CellWidth[1][j],2)/4.;*/
          // 2nd order
          qpres[IDX(i,j,k)] = qpres[IDX(i,j,k)]+
                              (logd[IDX(i,max(j-1,0),k)]+logd[IDX(i,min(j+1,jn-1),k)]
                             -2.*logd[IDX(i,j,k)])/(CellWidth[1][j]*CellWidth[1][j])/2.;
          qpres[IDX(i,j,k)] = qpres[IDX(i,j,k)] 
                  + pow(-1./2.*(logd[IDX(i,max(j-1,0),k)]-logd[IDX(i,min(j+1,jn-1),k)])/CellWidth[1][j],2)/4.;
          //qpres[IDX(i,j,k)] = qpres[IDX(i,j,k)]
          //+(pow(d[IDX(i,max(j-1,0),k)],0.5)-2.0*pow(d[IDX(i,j,k)],0.5)+pow(d[IDX(i,min(j+1,jn-1),k)],0.5))/pow(CellWidth[1][j],2);


          }//endif rank >1
          if (rank > 2) {
          //4th order
          /*qpres[IDX(i,j,k)] = qpres[IDX(i,j,k)]
          				  +(-1./12.*(logd[IDX(i,j,max(k-2,0))]+logd[IDX(i,j,min(k+2,kn-1))])
          				  +4./3.*(logd[IDX(i,j,max(k-1,0))]+logd[IDX(i,j,min(k+1,kn-1))])
          				  -5./2.*logd[IDX(i,j,k)])/(CellWidth[2][k]*CellWidth[2][k])/2.;
          qpres[IDX(i,j,k)] = qpres[IDX(i,j,k)] 
          				+ pow((1./12.*(logd[IDX(i,j,max(k-2,0))]-logd[IDX(i,j,min(k+2,kn-1))])
          				  -2./3.*(logd[IDX(i,j,max(k-1,0))]-logd[IDX(i,j,min(k+1,kn-1))]))/CellWidth[2][k],2)/4.;*/
          // 2nd order
          qpres[IDX(i,j,k)] = qpres[IDX(i,j,k)]+
                              (logd[IDX(i,j,max(k-1,0))]+logd[IDX(i,j,min(k+1,kn-1))]
                             -2.*logd[IDX(i,j,k)])/(CellWidth[2][k]*CellWidth[2][k])/2.;
          qpres[IDX(i,j,k)] = qpres[IDX(i,j,k)] 
                  + pow(-1./2.*(logd[IDX(i,j,max(k-1,0))]-logd[IDX(i,j,min(k+1,kn-1))])/CellWidth[2][k],2)/4.;
          //qpres[IDX(i,j,k)] = qpres[IDX(i,j,k)]
          //+(pow(d[IDX(i,i,max(k-1,0))],0.5)-2.0*pow(d[IDX(i,j,k)],0.5)+pow(d[IDX(i,j,min(k+1,kn-1))],0.5))/pow(CellWidth[2][k],2);

          }//endif rank>2
          qpres[IDX(i,j,k)] = qpres[IDX(i,j,k)]/(a*a);
          //qpres[IDX(i,j,k)] = qpres[IDX(i,j,k)]/pow(d[IDX(i,j,k)],0.5)/(a*a);
          //qpres[IDX(i,j,k)] = 2000.;
            }//end loop over i
          }//end loop over j
        }//end loop over k
        //fprintf(stderr, "ComputeQ: p,a,logd=%"GSYM",%"GSYM",%"GSYM",%"ISYM",%"ISYM",%"ISYM",\n", 
        //    qpres[0],a,d[0],in,jn,kn);

  delete [] logd;
  
  
  return SUCCESS;

}
