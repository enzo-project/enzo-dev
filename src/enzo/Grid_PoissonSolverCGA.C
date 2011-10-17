/***********************************************************************
/
/  GRID CLASS (SOLVE POISSON EQUATION IN CARTESIAN COORDINATE USING
/              CONJUGATE GRADIENT METHOD WITH JACOBI PRECONDITIONER)
/
/  written by: Peng Wang & Fen Zhao
/  date:       November, 2008
/  modified1:
/
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

int FindField(int field, int farray[], int numfields);

template <typename T>
int grid::multA(T *input, T *output, int *MatrixStartIndex, int *MatrixEndIndex )
{

  int diff[3];
  diff[0] = 1;
  diff[1] = (GridRank > 1) ? GridDimension[0] : 0;
  diff[2] = (GridRank > 2) ? GridDimension[0]*GridDimension[1] : 0;

  FLOAT dx_inv[3];
  dx_inv[0] = 1.0 / CellWidth[0][0];
  dx_inv[1] = (GridRank > 1) ? 1.0 / CellWidth[1][0] : 0.0;
  dx_inv[2] = (GridRank > 2) ? 1.0 / CellWidth[2][0] : 0.0;

  FLOAT dx_inv2[3]={dx_inv[0]*dx_inv[0], dx_inv[1]*dx_inv[1], dx_inv[2]*dx_inv[2]};
  

  int igrid;
  for (int k = MatrixStartIndex[2]; k <= MatrixEndIndex[2]; k++) {
    for (int j = MatrixStartIndex[1]; j <= MatrixEndIndex[1]; j++) {
      for (int i = MatrixStartIndex[0]; i <= MatrixEndIndex[0]; i++) {

	igrid = i + (j + k * GridDimension[1]) * GridDimension[0];

	output[igrid] = 0; 
	output[igrid] += input[igrid + diff[0]] * dx_inv2[0]; 
	output[igrid] += input[igrid - diff[0]] * dx_inv2[0]; 
	output[igrid] += input[igrid + diff[1]] * dx_inv2[1];
	output[igrid] += input[igrid - diff[1]] * dx_inv2[1];

	output[igrid] -= 2.0 * input[igrid] * dx_inv2[0];
	output[igrid] -= 2.0 * input[igrid] * dx_inv2[1];

	if (GridRank > 2){
	  output[igrid] += input[igrid + diff[2]] * dx_inv2[2];
	  output[igrid] += input[igrid - diff[2]] * dx_inv2[2];
	  output[igrid] -= 2.0 * input[igrid] * dx_inv2[2];
	}


      }
    }
  }

  return SUCCESS;

}

template <typename T>
int grid::multA2(T* input,T* output,  int *MatrixStartIndex, int *MatrixEndIndex)
{
  
  int diff[3];
  diff[0] = 1;
  diff[1] = (GridRank > 1) ? GridDimension[0] : 0;
  diff[2] = (GridRank > 2) ? GridDimension[0]*GridDimension[1] : 0;

  FLOAT dx_inv[3];
  dx_inv[0] = 1.0 / CellWidth[0][0];
  dx_inv[1] = (GridRank > 1) ? 1.0 / CellWidth[1][0] : 0.0;
  dx_inv[2] = (GridRank > 2) ? 1.0 / CellWidth[2][0] : 0.0;
  FLOAT dx_inv2[3]={dx_inv[0]*dx_inv[0]/4.0, dx_inv[1]*dx_inv[1]/4.0, dx_inv[2]*dx_inv[2]/4.0};

  int igrid;
  for (int k = MatrixStartIndex[2]; k <= MatrixEndIndex[2]; k++) {
    for (int j = MatrixStartIndex[1]; j <= MatrixEndIndex[1]; j++) {
      for (int i = MatrixStartIndex[0]; i <= MatrixEndIndex[0]; i++) {

	igrid = i + (j + k * GridDimension[1]) * GridDimension[0];
	
	output[igrid] = 0; 
	output[igrid] += input[igrid + 2*diff[0]] * dx_inv2[0]; 
	output[igrid] += input[igrid - 2*diff[0]] * dx_inv2[0]; 
	output[igrid] += input[igrid + 2*diff[1]] * dx_inv2[1];
	output[igrid] += input[igrid - 2*diff[1]] * dx_inv2[1];
	output[igrid] -= 2.0 * input[igrid] * dx_inv2[0];
	output[igrid] -= 2.0 * input[igrid] * dx_inv2[1];

	if (GridRank > 2){
	  output[igrid] += input[igrid + 2*diff[2]] * dx_inv2[2];
	  output[igrid] += input[igrid - 2*diff[2]] * dx_inv2[2];
	  output[igrid] -= 2.0 * input[igrid] * dx_inv2[2];}



	
	
      }
    }
  }

  return SUCCESS;
}

template <typename T>
T grid::dot(T *a, T *b, int size)
{
 
  T output = 0.0;

  
  for (int igrid = 0; igrid < size; igrid++) 
    output += a[igrid] * b[igrid]; 
    
    
  return output;
}

template <typename T>
int grid::setDirichletBC(T* x, int* MatrixStartIndex, int *MatrixEndIndex) {
  int ijk[3]={0,0,0}; 



  for (ijk[2] = 0; ijk[2] < GridDimension[2]; ijk[2]=ijk[2]+1) 
    for (ijk[1] = 0; ijk[1] < GridDimension[1]; ijk[1]=ijk[1]+1) 
      for (ijk[0] = 0; ijk[0] < GridDimension[0]; ijk[0]=ijk[0]+1)
	for (int i=0; i<3; i++)
	  if (ijk[i]<MatrixStartIndex[i] || ijk[i]> MatrixEndIndex[i])
	    x[ijk[0] + (ijk[1] + ijk[2] * GridDimension[1]) * GridDimension[0]]=0.0;
       
  return SUCCESS;


}

template <typename T>
int grid::setNeumannBC(T* x,  int* MatrixStartIndex, int *MatrixEndIndex, int type)
{
  int ijk[3]={0,0,0}; int index;

  int diff[3];
  diff[0] = 1;
  diff[1] = (GridRank > 1) ? GridDimension[0] : 0;
  diff[2] = (GridRank > 2) ? GridDimension[0]*GridDimension[1] : 0;
  




 //  if (type==0){

//   for (ijk[2] = MatrixStartIndex[2]; ijk[2] <= MatrixEndIndex[2]; ijk[2]=ijk[2]+1) {
//     for (ijk[1] = MatrixStartIndex[1]; ijk[1] <= MatrixEndIndex[1]; ijk[1]=ijk[1]+1) {
//       for (ijk[0] = MatrixStartIndex[0]; ijk[0] <= MatrixEndIndex[0]; ijk[0]=ijk[0]+1) {
       
// 	index =   ijk[0] + (ijk[1] + ijk[2] * GridDimension[1]) * GridDimension[0];
       
       
// 	if (ijk[0] == MatrixStartIndex[0]) {
// 	  x[index-diff[0]] = x[index+diff[0]]; 
// 	  x[index-2*diff[0]] = x[index];
// 	} else if (ijk[0] == MatrixEndIndex[0]) {
// 	  x[index+diff[0]] = x[index-diff[0]]; 
// 	  x[index+2*diff[0]] = x[index];
// 	}
       
      
// 	if (ijk[1] == MatrixStartIndex[1]) {
// 	  x[index-diff[1]] = x[index+diff[1]]; 
// 	  x[index-2*diff[1]] = x[index];
// 	} else if (ijk[1] == MatrixEndIndex[1]) {
// 	  x[index+diff[1]] = x[index-diff[1]]; 
// 	  x[index+2*diff[1]] = x[index];
// 	}

// 	if (GridRank > 2){
// 	  if (ijk[2] == MatrixStartIndex[2]) {
// 	    x[index-diff[2]] = x[index+diff[2]]; 
// 	    x[index-2*diff[2]] = x[index];
// 	  } else if (ijk[2] == MatrixEndIndex[2]) {
// 	    x[index+diff[2]] = x[index-diff[2]]; 
// 	    x[index+2*diff[2]] = x[index];
// 	  }
// 	}
//       }}}

//   T val=0.0;
//   for (ijk[2] = MatrixStartIndex[2]; ijk[2] <= MatrixEndIndex[2]; ijk[2]=ijk[2]+1) {
//     for (ijk[1] = MatrixStartIndex[1]; ijk[1] <= MatrixEndIndex[1]; ijk[1]=ijk[1]+1) {
//       for (ijk[0] = MatrixStartIndex[0]; ijk[0] <= MatrixEndIndex[0]; ijk[0]=ijk[0]+1) {
       
// 	index =   ijk[0] + (ijk[1] + ijk[2] * GridDimension[1]) * GridDimension[0];
       
       
// 	if (ijk[0] == MatrixStartIndex[0]) {
// 	  val=x[index+diff[0]];
// 	  x[index]=val;
// 	  x[index-diff[0]] = val;
// 	  x[index-2*diff[0]] = val;
// 	} else if (ijk[0] == MatrixEndIndex[0]) {
// 	  val=x[index-diff[0]];
// 	  x[index]=val;
// 	  x[index+diff[0]] = val; 
// 	  x[index+2*diff[0]] = val;
// 	}
       
 
// 	if (ijk[1] == MatrixStartIndex[1]) {
// 	  val=x[index+diff[1]];
// 	  x[index]=val;
// 	  x[index-diff[1]] = val;
// 	  x[index-2*diff[1]] = val;
// 	} else if (ijk[1] == MatrixEndIndex[1]) {
// 	  val=x[index-diff[1]];
// 	  x[index]=val;
// 	  x[index+diff[1]] = val; 
// 	  x[index+2*diff[1]] = val;
// 	}

// 	if (GridRank > 2){
// 	if (ijk[2] == MatrixStartIndex[2]) {
// 	  val=x[index+diff[2]];
// 	  x[index]=val;
// 	  x[index-diff[2]] = val;
// 	  x[index-2*diff[2]] = val;
// 	} else if (ijk[2] == MatrixEndIndex[2]) {
// 	  val=x[index-diff[2]];
// 	  x[index]=val;
// 	  x[index+diff[2]] = val; 
// 	  x[index+2*diff[2]] = val;
// 	}

// 	}
//       }}}
//  }
//  else{


  for (ijk[2] = MatrixStartIndex[2]; ijk[2] <= MatrixEndIndex[2]; ijk[2]=ijk[2]+1) 
    for (ijk[1] = MatrixStartIndex[1]; ijk[1] <= MatrixEndIndex[1]; ijk[1]=ijk[1]+1) 
      for (ijk[0] = MatrixStartIndex[0]; ijk[0] <= MatrixEndIndex[0]; ijk[0]=ijk[0]+1) {
	index =   ijk[0] + (ijk[1] + ijk[2] * GridDimension[1]) * GridDimension[0];
	if (ijk[0] == MatrixStartIndex[0] || ijk[0]-1 == MatrixStartIndex[0] ||  
	    ijk[0] == MatrixEndIndex[0] ||  ijk[0]+1 == MatrixEndIndex[0] ||
	    ijk[1] == MatrixStartIndex[1] || ijk[1]-1 == MatrixStartIndex[1] ||  
	    ijk[1] == MatrixEndIndex[1] ||  ijk[1]+1 == MatrixEndIndex[1] )
	  x[index]=0.0;
	
	if (GridRank > 2)
	  if (ijk[2] == MatrixStartIndex[2] || ijk[2]-1 == MatrixStartIndex[2] ||  
	      ijk[2] == MatrixEndIndex[2] ||  ijk[2]+1 == MatrixEndIndex[2]) 
	    x[index]=0.0;
      }//}
  
  return SUCCESS;
  
}

int grid::PoissonSolverCGA(int difftype, double *divB_p) 
{

  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }

  int size = 1;
  for (int dim = 0; dim < GridRank; dim++) {
    size *= GridDimension[dim];
  }

  int igrid;
  double threshold = PoissonApproximationThreshold;

  FLOAT dx_inv[3];
  dx_inv[0] = 1.0 / CellWidth[0][0];
  dx_inv[1] = (GridRank > 1) ? 1.0 / CellWidth[1][0] : 0.0;
  dx_inv[2] = (GridRank > 2) ? 1.0 / CellWidth[2][0] : 0.0;
 
  double *x     = new double[size];
  double *Ax    = new double[size];
  double *r_old = new double[size];
  double *r_new = new double[size];
  double *p     = new double[size];
 

  for (int i = 0; i < size; i++) {
    x[i]     = 0.0;
    Ax[i]    = 0.0;
    r_old[i] = 0.0;
    r_new[i] = 0.0;
    p[i]     = 0.0;
  }

  //Solve on the complete grid, only change the cleaning step so different boundaries

  int MatrixStartIndex[3]={ GridStartIndex[0],  GridStartIndex[1],  GridStartIndex[2]};
  int MatrixEndIndex[3]={ GridEndIndex[0],  GridEndIndex[1],  GridEndIndex[2]};

  /* Copy only the active part of BaryonField[Phi_pNum] to x so x 
     has zero values in the ghost cells. This will simplify the 
     matrix multiplication significantly. */

  int Phi_pNum = FindField(Phi_pField, FieldType, NumberOfBaryonFields);

 
    for (int k = MatrixStartIndex[2]; k <= MatrixEndIndex[2]; k++) {
      for (int j = MatrixStartIndex[1]; j <= MatrixEndIndex[1]; j++) {
	for (int i = MatrixStartIndex[0]; i <= MatrixEndIndex[0]; i++) {
	  igrid =  i + (j + k * GridDimension[1]) * GridDimension[0];
	  x[igrid] = BaryonField[Phi_pNum][igrid];
	}}}

   
    
     if (PoissonBoundaryType==Dirichlet)
      setDirichletBC(x, MatrixStartIndex, MatrixEndIndex); 
    
    
    if (difftype == 1) 
      multA(x, Ax, MatrixStartIndex, MatrixEndIndex);
    else if (difftype == 2) 
      multA2(x, Ax, MatrixStartIndex, MatrixEndIndex);

    if (PoissonBoundaryType==Neumann){
      setNeumannBC(Ax, MatrixStartIndex, MatrixEndIndex,1);
      setNeumannBC(divB_p, MatrixStartIndex, MatrixEndIndex,1);
    }

 
  /* r = divB_p - Ax */

  for (int k = MatrixStartIndex[2]; k <= MatrixEndIndex[2]; k++) 
    for (int j = MatrixStartIndex[1]; j <= MatrixEndIndex[1]; j++) 
      for (int i = MatrixStartIndex[0]; i <= MatrixEndIndex[0]; i++) {
	igrid =  i + (j + k * GridDimension[1]) * GridDimension[0]; 
	r_old[igrid] = divB_p[igrid] - Ax[igrid];	 
      }


  for (int k = MatrixStartIndex[2]; k <= MatrixEndIndex[2]; k++) 
    for (int j = MatrixStartIndex[1]; j <= MatrixEndIndex[1]; j++) 
      for (int i = MatrixStartIndex[0]; i <= MatrixEndIndex[0]; i++) {
	igrid =  i + (j + k * GridDimension[1]) * GridDimension[0]; 
	p[igrid] = r_old[igrid];
      }

  double *Ap = Ax; //reuse that allocated space for Ax
  double a, g, dotnr;

  bool converge = false;

  double r_norm = sqrt(dot(r_old, r_old, size)); 
  
 
  if (r_norm < threshold){
//    printf("r_norm1  %g < threshold %g\n", r_norm, threshold);  
    delete [] x;
    delete [] Ax;
    delete [] r_old;
    delete [] r_new;
    delete [] p;
    return SUCCESS;
  }

  int counter = 0;
  double maxnr;
  while (!converge && counter < size){

    if (PoissonBoundaryType==Dirichlet){
      setDirichletBC(p, MatrixStartIndex, MatrixEndIndex);       
    }
   

    counter++;
    
    if (difftype == 1) {
      multA(p, Ap, MatrixStartIndex, MatrixEndIndex);
    } else if (difftype == 2) {
      multA2(p, Ap, MatrixStartIndex, MatrixEndIndex);
    }

    if (PoissonBoundaryType==Neumann){
      setNeumannBC(Ap, MatrixStartIndex, MatrixEndIndex,1);
    }
        
    a = dot(r_old, r_old,size) / dot(p, Ap,size); 
  
    
    /* x = x + a * p */

  
    for (int k = MatrixStartIndex[2]; k <= MatrixEndIndex[2]; k++) 
      for (int j = MatrixStartIndex[1]; j <= MatrixEndIndex[1]; j++) 
	for (int i = MatrixStartIndex[0]; i <= MatrixEndIndex[0]; i++) {
	  igrid =  i + (j + k * GridDimension[1]) * GridDimension[0]; 
	  x[igrid] += a * p[igrid];}
   


    /* r_new = r - a * Ap */

    
    for (int k = MatrixStartIndex[2]; k <= MatrixEndIndex[2]; k++) 
      for (int j = MatrixStartIndex[1]; j <= MatrixEndIndex[1]; j++) 
	for (int i = MatrixStartIndex[0]; i <= MatrixEndIndex[0]; i++) {
	  igrid =  i + (j + k * GridDimension[1]) * GridDimension[0]; 
	  r_new[igrid] = r_old[igrid] - a * Ap[igrid];}


    dotnr = dot(r_new, r_new, size);
    g = dotnr / dot(r_old, r_old, size); 
  
    
    /* p = r_new + g*p */

    for (int k = MatrixStartIndex[2]; k <= MatrixEndIndex[2]; k++) 
      for (int j = MatrixStartIndex[1]; j <= MatrixEndIndex[1]; j++) 
	for (int i = MatrixStartIndex[0]; i <= MatrixEndIndex[0]; i++) {
	  igrid =  i + (j + k * GridDimension[1]) * GridDimension[0]; 
	  p[igrid] = r_new[igrid] + g * p[igrid]; 
	}

    /* r = r_new */

   
    for (int k = MatrixStartIndex[2]; k <= MatrixEndIndex[2]; k++) 
      for (int j = MatrixStartIndex[1]; j <= MatrixEndIndex[1]; j++) 
	for (int i = MatrixStartIndex[0]; i <= MatrixEndIndex[0]; i++) {
	  igrid =  i + (j + k * GridDimension[1]) * GridDimension[0]; 
	  r_old[igrid] = r_new[igrid];}


    /* maxnr = 0.0;
       for (int k = MatrixStartIndex[2]; k <= MatrixEndIndex[2]; k++) {
       for (int j = MatrixStartIndex[1]; j <= MatrixEndIndex[1]; j++) {
       for (int i = MatrixStartIndex[0]; i <= MatrixEndIndex[0]; i++) {
       igrid = i + (j + k * GridDimension[1]) * GridDimension[0];
       if (maxnr < abs(r_new[igrid])) 
       maxnr = abs(r_new[igrid]);
       }
       }
       } 
       printf("maxnr = %g\n", maxnr); */
    /* if (maxnr < threshold) converge = true; */

    r_norm = sqrt(dotnr);
    if (r_norm < threshold) converge = true;
    
  }


  if (debug){
//    printf("Iteration times = %d\n", counter); 
  
    if (counter == size)
      printf("Iterations reached limit, maxnr = %g\n", dotnr);
  }

 
  if (PoissonBoundaryType==Neumann)
    setNeumannBC(x, MatrixStartIndex, MatrixEndIndex,0);
  
  if (PoissonBoundaryType==Dirichlet){
    setDirichletBC(x, MatrixStartIndex, MatrixEndIndex); 
  }
    
  
  for (int igrid = 0; igrid < size; igrid++) 
    BaryonField[Phi_pNum][igrid] = x[igrid];
  
  
  delete [] x;
  delete [] Ax;
  delete [] r_old;
  delete [] r_new;
  delete [] p;

  return SUCCESS;


}














