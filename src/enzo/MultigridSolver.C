/***********************************************************************
/
/  MULTIGRID DRIVER
/
/  written by: Greg Bryan
/  date:       January, 1998
/  modified1:
/
/  PURPOSE:
/
/  NOTE:
/
************************************************************************/
 
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
 
/* function prototypes */
 
extern "C" void FORTRAN_NAME(mg_restrict)(
				       float *source, float *dest, int *ndim,
				       int *sdim1, int *sdim2, int *sdim3,
				       int *ddim1, int *ddim2, int *ddim3);
extern "C" void FORTRAN_NAME(mg_prolong)(float *source, float *dest, int *ndim,
					 int *sdim1, int *sdim2, int *sdim3,
					 int *ddim1, int *ddim2, int *ddim3);
extern "C" void FORTRAN_NAME(mg_prolong2)(float *source,float *dest, int *ndim,
					 int *sdim1, int *sdim2, int *sdim3,
					 int *ddim1, int *ddim2, int *ddim3);
extern "C" void FORTRAN_NAME(mg_calc_defect)(
			float *solution, float *rhs, float *defect, int *ndim,
			int *sdim1, int *sdim2, int *sdim3, float *norm);
extern "C" void FORTRAN_NAME(mg_relax)(float *solution, float *rhs, int *ndim,
				       int *sdim1, int *sdim2, int *sdim3);
 
 
#define MAX_DEPTH 100
#define PRE_SMOOTH 2
#define POST_SMOOTH 3
#define NUM_CYCLES 1
 
int MultigridSolver(float *TopRHS, float *TopSolution, int Rank, int TopDims[],
		    float &norm, float &mean, int start_depth,
		    float tolerance, int max_iter)
{
 
  /* declarations. */
 
  int i, dim, MinDim, bottom, cycle, smooth,
      Dims[MAX_DIMENSION][MAX_DEPTH], Size[MAX_DEPTH];
  float *Solution[MAX_DEPTH], *RHS[MAX_DEPTH], *defect[MAX_DEPTH];
  double lmean = 0.0;

  for (Size[0] = 1, dim = 0; dim < Rank; dim++)
    Size[0] *= (Dims[dim][0] = TopDims[dim]);
  for (dim = Rank; dim < MAX_DIMENSION; dim++)
    Dims[dim][0] = 1;
 
  Solution[0] = TopSolution;
  RHS[0]      = TopRHS;
 
  /* Compute dimensions and depth of V-cycle. */
 
  int depth = 0;
  for (depth = 0; depth < MAX_DEPTH; depth++) {
 
    /* Reduce size of dimensions. */
 
    Size[depth+1] = 1;
    for (MinDim = Dims[0][depth], dim = 0; dim < Rank; dim++) {
      Dims[dim][depth+1] = (Dims[dim][depth]+1)/2;
/*      if (Dims[dim][depth+1]*2-1 != Dims[dim][depth]+1) {
	ENZO_VFAIL("Dims[%"ISYM"]=%"ISYM" not of form 2^j+1\n", dim, Dims[dim][0])
      }
*/
      MinDim = min(Dims[dim][depth+1], MinDim);
      Size[depth+1] *= Dims[dim][depth+1];
    }
    for (dim = Rank; dim < MAX_DIMENSION; dim++)
      Dims[dim][depth+1] = 1;
 
    /* Check if this is the bottom (i.e. one dimension is too small). */
 
    if (MinDim < 3)
      break;
  }
  bottom = depth;
 
  /* Error check */
 
  if (depth == MAX_DEPTH) {
    ENZO_VFAIL("Depth(%"ISYM") > MAX_DEPTH\n", depth)
  }
 
  if (start_depth > bottom) {
    ENZO_VFAIL("Start depth(%"ISYM") > bottom(%"ISYM")!\n", start_depth, bottom)
  }
 
  /* Initial smoothing of density field, if requested. */
 
  for (depth = 0; depth < start_depth; depth++) {
    RHS[depth+1]      = new float[Size[depth+1]];
    FORTRAN_NAME(mg_prolong2)(RHS[depth], RHS[depth+1], &Rank,
    		     &Dims[0][depth  ], &Dims[1][depth  ], &Dims[2][depth  ],
    		     &Dims[0][depth+1], &Dims[1][depth+1], &Dims[2][depth+1]);
  }
  for (depth = start_depth; depth > 0; depth--) {
    FORTRAN_NAME(mg_prolong2)(RHS[depth], RHS[depth-1], &Rank,
		     &Dims[0][depth  ], &Dims[1][depth  ], &Dims[2][depth  ],
    		     &Dims[0][depth-1], &Dims[1][depth-1], &Dims[2][depth-1]);
    delete [] RHS[depth];
  }
 
  //  if (start_depth == bottom)
  //    defect[bottom] = new float[Size[bottom]];
 
  /* Iterate to convergence */
 
  int iter = 0;
  float tol_check = 2*tolerance;
 
  while (iter < max_iter && tol_check > tolerance) {
 
  /* Loop over number of V-cycles. */
 
  for (cycle = 0; cycle < NUM_CYCLES; cycle++) {
 
    /* Down first half of V-cycle. */
 
    for (depth = 0; depth < bottom; depth++) {
 
      /* Allocate memory. */
 
      if (cycle == 0 && iter == 0) {
	defect[depth]     = new float[Size[depth]];
	RHS[depth+1]      = new float[Size[depth+1]];
	Solution[depth+1] = new float[Size[depth+1]];
      }
 
      /* Pre-smoothing. */
 
      for (smooth = 0; smooth < PRE_SMOOTH; smooth++)
	FORTRAN_NAME(mg_relax)(Solution[depth], RHS[depth], &Rank,
			    &Dims[0][depth], &Dims[1][depth], &Dims[2][depth]);
 
      /* Compute the defect. */
 
      FORTRAN_NAME(mg_calc_defect)(Solution[depth], RHS[depth], defect[depth],
				&Rank, &Dims[0][depth], &Dims[1][depth],
				&Dims[2][depth], &norm);
 
      /* Restrict the defect. */
 
//    printf("restricting defect %"ISYM" -> %"ISYM"\n",Dims[0][depth],Dims[0][depth+1]);
      FORTRAN_NAME(mg_restrict)(defect[depth], RHS[depth+1], &Rank,
//      FORTRAN_NAME(mg_prolong2)(defect[depth], RHS[depth+1], &Rank,
		      &Dims[0][depth  ], &Dims[1][depth  ], &Dims[2][depth  ],
		      &Dims[0][depth+1], &Dims[1][depth+1], &Dims[2][depth+1]);
 
      /* Set guess to correction equal to zero. */
 
      for (i = 0; i < Size[depth+1]; i++)
	Solution[depth+1][i] = 0;
	
    } // end of first half of V-cycle
 
    /* Smooth on bottom level. */
 
    for (smooth = 0; smooth < 3*PRE_SMOOTH; smooth++)
      FORTRAN_NAME(mg_relax)(Solution[bottom], RHS[bottom], &Rank,
		       &Dims[0][bottom], &Dims[1][bottom], &Dims[2][bottom]);
 
    /* Back up second half of V-cycle. */
 
    for (depth = bottom-1; depth >= 0; depth--) {
 
      /* Prolong coarse correction to next level. */
 
      //    printf("prolonging Solution %"ISYM" -> %"ISYM"\n",Dims[0][depth+1],Dims[0][depth]);
      FORTRAN_NAME(mg_prolong)(Solution[depth+1], defect[depth], &Rank,
		     &Dims[0][depth+1], &Dims[1][depth+1], &Dims[2][depth+1],
		     &Dims[0][depth  ], &Dims[1][depth  ], &Dims[2][depth  ]);
 
      /* Add solution correction to the one on this level. */
 
      for (i = 0; i < Size[depth]; i++)
	Solution[depth][i] += defect[depth][i];
 
      /* Post-prolongation smoothing. */
 
      for (smooth = 0; smooth < POST_SMOOTH; smooth++)
	FORTRAN_NAME(mg_relax)(Solution[depth], RHS[depth], &Rank,
		        &Dims[0][depth], &Dims[1][depth], &Dims[2][depth]);
 
    } // end loop over second-half of V-cycle
 
  } // end loop over V-cycles
 
  /* Calculate norm. */
 
  FORTRAN_NAME(mg_calc_defect)(Solution[0], RHS[0], defect[0], &Rank,
			       &Dims[0][0], &Dims[1][0], &Dims[2][0], &norm);
 
  /* Calculate mean. */
 
  lmean = 0;
  for (i = 0; i < Size[0]; i++)
    lmean += fabs(Solution[0][i]);
  lmean /= float(Size[0]);
  mean = lmean;

  iter++;
  tol_check = norm/mean;
 
//  printf("%"ISYM" (%"ISYM" %"ISYM" %"ISYM") %"GSYM" %"GSYM" %"GSYM"\n", iter, Dims[0][0], Dims[1][0],
//	 Dims[2][0], norm, mean, tol_check);
 
  } // end: iteration loop
 
  int repeat = 0;
  while (repeat < 200 && tol_check > tolerance) {
    FORTRAN_NAME(mg_relax)(Solution[0], RHS[0], &Rank,
			   &Dims[0][0], &Dims[1][0], &Dims[2][0]);
    FORTRAN_NAME(mg_calc_defect)(Solution[0], RHS[0], defect[0], &Rank,
				 &Dims[0][0], &Dims[1][0], &Dims[2][0], &norm);
    lmean = 0;
    for (i = 0; i < Size[0]; i++)
      lmean += fabs(Solution[0][i]);
    lmean /= float(Size[0]);
    mean = lmean;
    tol_check = norm/mean;
    //    printf("%"ISYM" (%"ISYM" %"ISYM" %"ISYM") %"GSYM" %"GSYM" %"GSYM"\n", repeat, Dims[0][0], Dims[1][0],
    //	   Dims[2][0], norm, mean, tol_check);
    repeat++;
  }
 
  if (tol_check > tolerance) {
    ENZO_VFAIL("Too many iterations (%"ISYM"): tol=%"GSYM", check=%"GSYM"\n", iter,
	    tolerance, tol_check)

  }
 
  /* Free allocated memory. */
 
  for (depth = 1; depth <= bottom; depth++) {
    delete [] Solution[depth];
    delete [] RHS[depth];
    delete [] defect[depth-1];
  }
 
  return SUCCESS;
}
