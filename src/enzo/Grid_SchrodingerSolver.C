/***********************************************************************
/
/  GRID CLASS (WRAPPER FOR RK-TVD Schrodinger Eqn)
/
/  written by: Xinyu Li
/  date:       May, 2017
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

int rk3tvd(float *repsi, float *impsi, float *repsi0, float *impsi0,
         float alpha, float beta,
	       int in, int jn, int kn, int rank,
	       int is, int ie, int js, int je, int ks, int ke, 
	       float dt, float dx[], float dy[], float dz[],
	       float hmcoef,
         int gravity, double *p);

int rk4(double *repsi, double *impsi,
         int in, int jn, int kn, int rank,
         double dt, double dx[], double dy[], double dz[],
         double hmcoef);

int SchrdingerAddPotential(double *repsi, double *impsi,
         int in, int jn, int kn, int rank,
         int gin, int gjn, int gkn,
         double dt, 
         double hmcoef,
         double *p, int start1, int start2, int start3);

int GetUnits (float *DensityUnits, float *LengthUnits,
         float *TemperatureUnits, float *TimeUnits,
         float *VelocityUnits, double *MassUnits, FLOAT Time);
int FindField(int field, int farray[], int numfields);
int QuantumGetUnits(float *DensityUnits, float *LengthUnits,
       float *TemperatureUnits, float *TimeUnits,
       float *VelocityUnits, FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);


int grid::SchrodingerSolver( int nhy )
{
  /*  Locals */
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  int i, ie, is, j, je, js, k, ks, ke, n, ixyz, ret;
  //int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  //float *d, *e, *u, *v, *w;
  int dim;

  /* compute global start index for left edge of entire grid 
       (including boundary zones) */
  Elong_int GridGlobalStart[MAX_DIMENSION];

  for (dim = 0; dim < GridRank; dim++)
      GridGlobalStart[dim] = nlongint((GridLeftEdge[dim]-DomainLeftEdge[dim])/(*(CellWidth[dim]))) -
  GridStartIndex[dim];

    /* fix grid quantities so they are defined to at least 3 dims */

    for (i = GridRank; i < 3; i++) {
      GridDimension[i]   = 1;
      GridStartIndex[i]  = 0;
      GridEndIndex[i]    = 0;
      GridVelocity[i]    = 0.0;
      GridGlobalStart[i] = 0;
    }

  FLOAT a = 1, dadt;
  if (ComovingCoordinates)
      if (CosmologyComputeExpansionFactor(Time+0.5*dtFixed, &a, &dadt) 
    == FAIL) {
  ENZO_FAIL("Error in CosmologyComputeExpansionFactors.");
      }

    /* Create a cell width array to pass (and convert to absolute coords). */

    float *CellWidthTemp[MAX_DIMENSION];
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      CellWidthTemp[dim] = new float[GridDimension[dim]];
      for (i = 0; i < GridDimension[dim]; i++)
  if (dim < GridRank)
    CellWidthTemp[dim][i] = float(a*CellWidth[dim][i]);
  else
    CellWidthTemp[dim][i] = 1.0;
    }

      /* Prepare Gravity. */

  int GravityOn = 0, FloatSize = sizeof(float);
    if (SelfGravity || UniformGravity || PointSourceGravity || DiskGravity || ExternalGravity )
      GravityOn = 1;
#ifdef TRANSFER
    if (RadiationPressure)
      GravityOn = 1;
#endif  

  /* Error Check */

  for (i = 0; i < GridRank; i++)
    if (GridDimension[i] > MAX_ANY_SINGLE_DIRECTION) {
      ENZO_FAIL("Schrodinger_MAIN: A grid dimension is too long (increase max_any_single_direction.)\n");
    }

  /* Allocate temporary space for Zeus_Hydro. */

  int size = GridDimension[0]*GridDimension[1]*GridDimension[2];

  //float *repsi0 = new float[size];
  //float *impsi0 = new float[size];
  
  /* Find fields: density, total energy, velocity1-3 and set pointers to them
     Create zero fields for velocity2-3 for low-dimension runs because solver
     assumes they exist. */

  /*this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum );
  
  d = BaryonField[DensNum];
  e = BaryonField[TENum];
  u = BaryonField[Vel1Num];
  v = BaryonField[Vel2Num];
  w = BaryonField[Vel3Num]; */

  /*  Set grid start and end indicies */

  /*is = GridStartIndex[0];
  js = GridStartIndex[1];
  ks = GridStartIndex[2];
  ie = GridEndIndex[0];
  je = GridEndIndex[1];
  ke = GridEndIndex[2];*/

  /* if use more ghost zones */  
  /*ie = ie+max(0,is-3);
  je = je+max(0,js-3);
  ke = ke+max(0,ks-3);
  is = is-max(0,is-3);
  js = js-max(0,js-3);
  ks = ks-max(0,ks-3);*/

		 
  /*   1) Add source terms */  
    // get code units
  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
    VelocityUnits = 1, TimeUnits = 1;

  if (QuantumGetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
         &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.");
  }

// set poteitial offset  
int Offset[MAX_DIMENSION] = {0,0,0};
for (int dim = 0; dim < GridRank; dim++) {
    Offset[dim] = nint((CellLeftEdge[dim][0] -
      GravitatingMassFieldLeftEdge[dim])/ CellWidth[dim][0]);
    //printf("offset %d %d\n", dim , GravitatingMassFieldDimension[dim]);
  }

  // calculate hbar/m
  float hmcoef = 5.9157166856e27*TimeUnits/pow(LengthUnits,2)/FDMMass;
  //printf("hmcoef %f \n", hmcoef);

  int RePsiNum, ImPsiNum, FDMDensNum;
  RePsiNum = FindField(RePsi, FieldType, NumberOfBaryonFields);
  ImPsiNum = FindField(ImPsi, FieldType, NumberOfBaryonFields);
  FDMDensNum = FindField(FDMDensity, FieldType, NumberOfBaryonFields);

  float *repsi, *impsi, *d;
  repsi = BaryonField[RePsiNum];
  impsi = BaryonField[ImPsiNum];
  d = BaryonField[FDMDensNum];

  // Add half gravity before the TVD-RK
  if (GravityOn && (PotentialField != NULL) ){
    if (nhy == 0 ){ // Add half step of potential at the beginning
    if(SchrdingerAddPotential(repsi, impsi,
         GridDimension[0], GridDimension[1], GridDimension[2],
         GridRank, 
         GravitatingMassFieldDimension[0],GravitatingMassFieldDimension[1],GravitatingMassFieldDimension[2],
         dtFixed/2., hmcoef,
         PotentialField, 
         Offset[0], Offset[1], Offset[2]) == FAIL){
         ENZO_FAIL("Error in Add poteitial before RK!\n"); 
         }
     } else {

    if(SchrdingerAddPotential(repsi, impsi,
         GridDimension[0], GridDimension[1], GridDimension[2],
         GridRank, 
         GravitatingMassFieldDimension[0],GravitatingMassFieldDimension[1],GravitatingMassFieldDimension[2],
         dtFixed, hmcoef,
         PotentialField, 
         Offset[0], Offset[1], Offset[2]) == FAIL){
         ENZO_FAIL("Error in Add poteitial before RK!\n");
         }
    }
  }

 /* cache wave-function before TVD-RK */
  /*for (int i=0; i<size; i++){
    repsi0[i] = repsi[i];
    impsi0[i] = impsi[i];
  }

  // 1st step of rk3tvd
  
    if (rk3tvd(repsi, impsi, repsi0, impsi0, 0.0, 1.0,
     GridDimension[0], GridDimension[1], GridDimension[2],
     GridRank, 
     is, ie, js, je, ks, ke, 
     dtFixed, dx, dy, dz, hmcoef, gravity, PotentialField) == FAIL) {
    fprintf(stderr, "P(%"ISYM"): Error in rk3tvd on step %"ISYM" (dt=%"GSYM")\n", MyProcessorNumber,
      nhy, dtFixed);
    fprintf(stderr, "  grid dims = %"ISYM" %"ISYM" %"ISYM"\n", GridDimension[0], GridDimension[1], GridDimension[2]);
    ENZO_FAIL("Error in rk3tvd 1st step!\n");
    }

  // 2nd step of rk3tvd 
    if (rk3tvd(repsi, impsi, repsi0, impsi0, 0.75, 0.25,
     GridDimension[0], GridDimension[1], GridDimension[2],
     GridRank, 
     is, ie, js, je, ks, ke, 
     dtFixed, dx, dy, dz, hmcoef, gravity, PotentialField) == FAIL) {
    fprintf(stderr, "P(%"ISYM"): Error in rk3tvd on step %"ISYM" (dt=%"GSYM")\n", MyProcessorNumber,
      nhy, dtFixed);
    fprintf(stderr, "  grid dims = %"ISYM" %"ISYM" %"ISYM"\n", GridDimension[0], GridDimension[1], GridDimension[2]);
    ENZO_FAIL("Error in rk3tvd 2nd step!\n");
    }

  // final step of rk3tvd
    if (rk3tvd(repsi, impsi, repsi0, impsi0, 1./3., 2./3.,
     GridDimension[0], GridDimension[1], GridDimension[2],
     GridRank, 
     is, ie, js, je, ks, ke, 
     dtFixed, dx, dy, dz, hmcoef, gravity, PotentialField) == FAIL) {
    fprintf(stderr, "P(%"ISYM"): Error in rk3tvd on step %"ISYM" (dt=%"GSYM")\n", MyProcessorNumber,
      nhy, dtFixed);
    fprintf(stderr, "  grid dims = %"ISYM" %"ISYM" %"ISYM"\n", GridDimension[0], GridDimension[1], GridDimension[2]);
    ENZO_FAIL("Error in rk3tvd 3rd step!\n");
    }
  // Clean up 
  delete [] repsi0;
  delete [] impsi0;
    */

  // Runge-Kutta 4th order
  if ( rk4(repsi, impsi, 
          GridDimension[0], GridDimension[1], GridDimension[2], GridRank, 
          dtFixed, CellWidthTemp[0], CellWidthTemp[1], CellWidthTemp[2], hmcoef) == FAIL){
    fprintf(stderr, "P(%"ISYM"): Error in rk4 on step %"ISYM" (dt=%"GSYM")\n", MyProcessorNumber,
      nhy, dtFixed);
    fprintf(stderr, "  grid dims = %"ISYM" %"ISYM" %"ISYM"\n", GridDimension[0], GridDimension[1], GridDimension[2]);
    ENZO_FAIL("Error in 4th runge kutta advance\n");
    }


  // Add half gravity after the TVD-RK
  /*if (GravityOn && (PotentialField != NULL) ){
    if(SchrdingerAddPotential(repsi, impsi,
         GridDimension[0], GridDimension[1], GridDimension[2],
         GridRank, 
         GravitatingMassFieldDimension[0],GravitatingMassFieldDimension[1],GravitatingMassFieldDimension[2],
         is, ie, js, je, ks, ke, 
         dtFixed/2., hmcoef,
         PotentialField, 
         Offset[0], Offset[1], Offset[2]) == FAIL){
         ENZO_FAIL("Error in Add poteitial after TVD-RK!\n");
    }
  }*/


  // Update New Density

  for (int i=0; i<size; i++){
    d[i] = repsi[i]*repsi[i]+impsi[i]*impsi[i];
  }

  
  
  return SUCCESS;

}
