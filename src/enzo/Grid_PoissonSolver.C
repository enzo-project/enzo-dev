/***********************************************************************
/
/  GRID CLASS (USE PROJECTION METHOD TO CLEAN DIVB)
/
/  written by: Peng Wang & Fen Zhao
/  date:       November, 2008
/  modified1:
/
/
************************************************************************/

#define DEBUG
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include <string.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

int FastFourierTransform(float *buffer, int Rank, int DimensionReal[], 
			 int Dimension[], int direction, int type);
int FindField(int field, int farray[], int numfields);

int MultigridSolver(float *TopRHS, float *TopSolution, int Rank, int TopDims[],
		    float &norm, float &mean, int start_depth, 
		    float tolerance, int max_iter);

int grid::PoissonSolver(int level) 
 /* 
     Input: type: 
     1= SOR (Successive OverRelaxation)
     2= SOR2 (with 4dx differencing)
     3= CGA (Conjugate Gradient)
     4= CGA2 (with 4dx differencing)
     5= FFT
     6= Multigrid     

     Methods 2 is sort of a failed experiment.  
     3 and 4 work very well.  
     1, 5, and 6 are placeholders for further work.
  */
{

  if (ProcessorNumber != MyProcessorNumber || UsePoissonDivergenceCleaning==0 || UseMHD==0 ) {
   
    return SUCCESS;
  }
  
  int iBx=FindField(Bfield1, FieldType, NumberOfBaryonFields);
  int iBy=FindField(Bfield2, FieldType, NumberOfBaryonFields);
  int iBz;
  if (GridRank==3) iBz=FindField(Bfield3, FieldType, NumberOfBaryonFields);

  /* Calculating divB_p */
  
  int diff[3];
  diff[0] = 1;
  diff[1] = (GridRank > 1) ? GridDimension[0] : 0;
  diff[2] = (GridRank > 2) ? GridDimension[0] * GridDimension[1] : 0;
  FLOAT dx[3];
  dx[0] = CellWidth[0][0];
  dx[1] = (GridRank > 1) ? CellWidth[1][0] : 1.0;
  dx[2] = (GridRank > 2) ? CellWidth[2][0] : 1.0;
  int size = 1;
  for (int dim = 0; dim < GridRank; dim++) {
    size *= GridDimension[dim];
  }

  float monopoleDensity = 0.0, current, Bdiff;

  double *divB_p = new double[size];

  int igrid;

  int MatrixStartIndex[3]={GridStartIndex[0],  GridStartIndex[1],  GridStartIndex[2]};
  int MatrixEndIndex[3]={GridEndIndex[0],  GridEndIndex[1],  GridEndIndex[2]};
  

  for (int igrid = 0; igrid < size; igrid++) {
    divB_p[igrid] =0.0;
  }

  for (int k = MatrixStartIndex[2]; k <= MatrixEndIndex[2]; k++) {
    for (int j = MatrixStartIndex[1]; j <= MatrixEndIndex[1]; j++) {
      for (int i = MatrixStartIndex[0]; i <= MatrixEndIndex[0]; i++) {
	igrid = GetIndex(i, j, k);
	
	divB_p[igrid] =
	  (BaryonField[iBx][igrid + diff[0]] - BaryonField[iBx][igrid - diff[0]]) / (2.0*dx[0]) +
	  (BaryonField[iBy][igrid + diff[1]] - BaryonField[iBy][igrid - diff[1]]) / (2.0*dx[1]);

	Bdiff=  
	  (BaryonField[iBx][igrid + diff[0]] - BaryonField[iBx][igrid - diff[0]]) / 2.0 +
	  (BaryonField[iBy][igrid + diff[1]] - BaryonField[iBy][igrid - diff[1]]) / 2.0;

	if (GridRank > 2){
	  divB_p[igrid] += (BaryonField[iBz][igrid + diff[2]] - BaryonField[iBz][igrid - diff[2]]) / (2.0*dx[2]);
	  Bdiff+=(BaryonField[iBz][igrid + diff[2]] - BaryonField[iBz][igrid - diff[2]]) / 2.0;
	}


	if (GridRank == 2){
	  current=sqrt(Bdiff*Bdiff/
	    (BaryonField[iBx][igrid]*BaryonField[iBx][igrid]+
	     BaryonField[iBy][igrid]*BaryonField[iBy][igrid]));}
	else if (GridRank ==3){
	  current=sqrt(Bdiff*Bdiff/
	    (BaryonField[iBx][igrid]*BaryonField[iBx][igrid]+
	     BaryonField[iBy][igrid]*BaryonField[iBy][igrid]+
	     BaryonField[iBz][igrid]*BaryonField[iBz][igrid]));}
	if (monopoleDensity<current) monopoleDensity=current;
      }
    }
  }

  if (monopoleDensity<PoissonDivergenceCleaningThreshold){
    if (debug) printf("Monopole Density %g < %g\n", monopoleDensity, PoissonDivergenceCleaningThreshold);
    return SUCCESS;
  }
//  else
//    printf("**Poisson Cleaning**\n");

  int type=UsePoissonDivergenceCleaning;

   if(debug){
   bool badDiv=false;
    float divSum = 0;
    for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
	  igrid = i + (j + k * GridDimension[1]) * GridDimension[0];
	  divSum += fabs(divB_p[igrid]/dx[0]);
	}
      }
    }
    
//    printf("Initial divB_p: %g (%g/%d) \n", divSum/size, divSum, size); 
   }
    
 

 

  if (type == 1) PoissonSolverSOR();
  else if (type == 2) PoissonSolverSOR2();
  else if (type == 3) PoissonSolverCGA(1, divB_p);
  else if (type == 4) PoissonSolverCGA(2, divB_p);
  else if (type == 5) PoissonSolverFFT();
  else if (type == 6) PoissonSolverMultigrid();   
 

  PoissonCleanStep(level);

  delete [] divB_p;
  
  return SUCCESS;
  
}


int grid::PoissonCleanStep(int level)
{

  int iBx=FindField(Bfield1, FieldType, NumberOfBaryonFields);
  int iBy=FindField(Bfield2, FieldType, NumberOfBaryonFields);
  int iBz;
  if (GridRank==3) iBz=FindField(Bfield3, FieldType, NumberOfBaryonFields);
  

  int size = 1;
  for (int dim = 0; dim < GridRank; dim++) {
    size *= GridDimension[dim];
  }

  int igrid;
  int diff[3];
  diff[0] = 1;
  diff[1] = (GridRank > 1) ? GridDimension[0] : 0;
  diff[2] = (GridRank > 2) ? GridDimension[0]*GridDimension[1] : 0;
  FLOAT dx[3];
  dx[0] = CellWidth[0][0];
  dx[1] = (GridRank > 1) ? CellWidth[1][0] : 1.0;
  dx[2] = (GridRank > 2) ? CellWidth[2][0] : 1.0;

  float Bx, By, Bz, B2, rho;
  
  int x=PoissonDivergenceCleaningBoundaryBuffer;

  int MatrixStartIndex[3]={GridStartIndex[0]+x,  GridStartIndex[1]+x,  GridStartIndex[2]+x};
  int MatrixEndIndex[3]={GridEndIndex[0]-x,  GridEndIndex[1]-x,  GridEndIndex[2]-x};

  if (GridRank==2) {
    MatrixEndIndex[2]=0;
    MatrixStartIndex[2]=0;
  }

  //int MatrixStartIndex[3]={GridStartIndex[0],  GridStartIndex[1],  GridStartIndex[2]};
  //int MatrixEndIndex[3]={GridEndIndex[0],  GridEndIndex[1],  GridEndIndex[2]};

  /* Subtract the old B energy from total energy */

 
  for (int k = MatrixStartIndex[2]; k <= MatrixEndIndex[2]; k++) {
    for (int j = MatrixStartIndex[1]; j <= MatrixEndIndex[1]; j++) {
      for (int i = MatrixStartIndex[0]; i <= MatrixEndIndex[0]; i++) {
	
	igrid = GetIndex(i, j, k);

	Bx  = BaryonField[iBx][igrid];
	By  = BaryonField[iBy][igrid];
	Bz  = (GridRank > 2) ? BaryonField[iBz][igrid]:0;
	B2  = Bx * Bx + By * By + Bz * Bz;
	rho = BaryonField[iden][igrid];
	BaryonField[ietot][igrid] -= 0.5 * B2 / rho;
      }
    }
  }  

  /* Clean B */

  int Phi_pNum = FindField(Phi_pField, FieldType, NumberOfBaryonFields);
  float *Phi_p = BaryonField[Phi_pNum];

  for (int k = MatrixStartIndex[2]; k <= MatrixEndIndex[2]; k++) {
    for (int j = MatrixStartIndex[1]; j <= MatrixEndIndex[1]; j++) {
      for (int i = MatrixStartIndex[0]; i <= MatrixEndIndex[0]; i++) {

	igrid = GetIndex(i, j, k);

	BaryonField[iBx][igrid] -= 
	  ( Phi_p[igrid + diff[0]] - Phi_p[igrid - diff[0]] ) / (2.0 * dx[0]);
	BaryonField[iBy][igrid] -= 
	  ( Phi_p[igrid + diff[1]] - Phi_p[igrid - diff[1]] ) / (2.0 * dx[1]);

	if (GridRank > 2)
	  BaryonField[iBz][igrid] -= 
	    ( Phi_p[igrid + diff[2]] - Phi_p[igrid - diff[2]] ) / (2.0 * dx[2]);

      }
    }
  }
  
  /* Add the cleaned B to total energy */

  for (int k = MatrixStartIndex[2]; k <= MatrixEndIndex[2]; k++) {
    for (int j = MatrixStartIndex[1]; j <= MatrixEndIndex[1]; j++) {
      for (int i = MatrixStartIndex[0]; i <= MatrixEndIndex[0]; i++) {

	igrid = GetIndex(i, j, k);

	Bx = BaryonField[iBx][igrid];
	By = BaryonField[iBy][igrid];
	Bz = (GridRank > 2) ? BaryonField[iBz][igrid]:0;
	B2 = Bx * Bx + By * By + Bz * Bz;
	rho = BaryonField[iden][igrid];
	BaryonField[ietot][igrid] += 0.5 * B2 / rho;

      }
    }
  }
  

  if (debug){

  double *divB_p = new double[size];
  float divSum = 0;

  for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {

	igrid = i + (j + k * GridDimension[1]) * GridDimension[0];

	divB_p[igrid] =
	  ( BaryonField[iBx][igrid + diff[0]] - BaryonField[iBx][igrid - diff[0]] ) / 2.0 +
	  ( BaryonField[iBy][igrid + diff[1]] - BaryonField[iBy][igrid - diff[1]] ) / 2.0;
	
	if (GridRank > 2)
	  divB_p[igrid] +=( BaryonField[iBz][igrid + diff[2]] - BaryonField[iBz][igrid - diff[2]] ) / 2.0;

	divSum += fabs(divB_p[igrid]/dx[0]);

      }
    }
  }  

  
//  printf("End divB_p: %g (%g/%d) \n", divSum/size, divSum, size);   //#####
  
  delete [] divB_p;
  }

  return SUCCESS;
  
}

// int grid::PoissonSolverDirichletBC(int d, double *divB_p) 

// {

  

//   const int ng = NumberOfGhostZones;
//   int igrid, igrid_rhs;

//   int Phi_pNum = FindField(Phi_pField, FieldType, NumberOfBaryonFields);

//   FLOAT dx_inv[3];
//   dx_inv[0] = 1.0 / CellWidth[0][0];
//   dx_inv[1] = (GridRank > 1) ? 1.0 / CellWidth[1][0] : 0.0;
//   dx_inv[2] = (GridRank > 2) ? 1.0 / CellWidth[2][0] : 0.0;

  
//   for (int i=0; i<3; i++)
//     dx_inv[i]=d*d*dx_inv[i]*dx_inv[i];


    

//   /* x boundary */

//   for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) 
//     for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {

//       /* left boundary */

//       igrid = GetIndex(ng - d + 1, j, k);
//       igrid_rhs = GetIndex(ng + 1, j, k);
//       divB_p[igrid_rhs] -= BaryonField[Phi_pNum][igrid] * dx_inv[0];
      
//       igrid = GetIndex(ng - d, j, k);
//       igrid_rhs = GetIndex(ng, j, k);
//       divB_p[igrid_rhs] -= BaryonField[Phi_pNum][igrid] * dx_inv[0];
      
//       /* right bounary */

//       igrid = GetIndex(GridEndIndex[0] + d - 1, j, k);
//       igrid_rhs = GetIndex(GridEndIndex[0] - 1, j, k);
//       divB_p[igrid_rhs] -= BaryonField[Phi_pNum][igrid] * dx_inv[0];
      
//       igrid = GetIndex(GridEndIndex[0] + d, j, k);
//       igrid_rhs = GetIndex(GridEndIndex[0], j, k);
//       divB_p[igrid_rhs] -= BaryonField[Phi_pNum][igrid] * dx_inv[0];      

//     }

//   /* y boundary */
  
//   if (GridRank > 1) {

//     for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
//       for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {

// 	/* left boundary */

// 	igrid = GetIndex(i, ng - d + 1, k);
// 	igrid_rhs = GetIndex(i, ng + 1, k);
// 	divB_p[igrid_rhs] -= BaryonField[Phi_pNum][igrid] * dx_inv[1];

// 	igrid = GetIndex(i, ng - d, k);
// 	igrid_rhs = GetIndex(i, ng, k);
// 	divB_p[igrid_rhs] -= BaryonField[Phi_pNum][igrid] * dx_inv[1];

// 	/* right boundary */

// 	igrid = GetIndex(i, GridEndIndex[1] + d -1, k);
// 	igrid_rhs = GetIndex(i, GridEndIndex[1] - 1, k);
// 	divB_p[igrid_rhs] -= BaryonField[Phi_pNum][igrid] * dx_inv[1];

// 	igrid = GetIndex(i, GridEndIndex[1] + d, k);
// 	igrid_rhs = GetIndex(i, GridEndIndex[1], k);
// 	divB_p[igrid_rhs] -= BaryonField[Phi_pNum][igrid] * dx_inv[1];

//       }

//   }

//   /* z boundary */

//   if (GridRank > 2) {

//     for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++)
//       for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
	
// 	/* left boundary */

// 	igrid = GetIndex(i, j, ng - d + 1);
// 	igrid_rhs = GetIndex(i, j, ng + 1);
// 	divB_p[igrid_rhs] -= BaryonField[Phi_pNum][igrid] * dx_inv[2];

// 	igrid = GetIndex(i, j, ng - d);
// 	igrid_rhs = GetIndex(i, j, ng);
// 	divB_p[igrid_rhs] -= BaryonField[Phi_pNum][igrid] * dx_inv[2];

// 	/* right boundary */

// 	igrid = GetIndex(i, j, GridEndIndex[2] + d - 1);
// 	igrid_rhs = GetIndex(i, j, GridEndIndex[2] - 1);
// 	divB_p[igrid_rhs] -= BaryonField[Phi_pNum][igrid] * dx_inv[2];

// 	igrid = GetIndex(i, j, GridEndIndex[2] + d);
//         igrid_rhs = GetIndex(i, j, GridEndIndex[2]);
// 	divB_p[igrid_rhs] -= BaryonField[Phi_pNum][igrid] * dx_inv[2];

//       }
    
//   }
  


//   return SUCCESS;

// }





//**********
//Other Solvers

int grid::PoissonSolverSOR()
{    
    return SUCCESS;
}


int grid::PoissonSolverSOR2()

{

#ifdef NOUSE
  int ijk[3];
  int index, activeIndex;
  int diffs[3]={1, GridDimension[0], GridDimension[0]*GridDimension[1]};
  int size=GridDimension[0]*GridDimension[1]*GridDimension[2];
  float dx[3]={CellWidth[0][0],CellWidth[1][0],CellWidth[2][0]};
  float dx2[3]={dx[0]*dx[0],dx[1]*dx[1],dx[2]*dx[2]};
 

   //Calculating poisson solution, SOR method

 

    float newvalue;
    //float converge=max(GridDimension[0]+1, max(GridDimension[1]+1, GridDimension[2]+1))/(Pi*2)+1;
    bool notconverge=true; float old;
    int counter=0;

    int maxdim=max(GridDimension[0]+1, max(GridDimension[1]+1, GridDimension[2]+1));
    
    float srad=(cos(Pi/GridDimension[0])+cos(Pi/GridDimension[1])+cos(Pi/GridDimension[2]))/3;
  float w=2/(1+pow(1-pow(srad,2),0.5));
  // printf(stderr, "w= %f\n", w); 
  //  w=1;

    //printf(stderr, "converge: %f \n", converge);
    //for (int n=0; n<converge; n++){
    

 
        
  while(notconverge && counter<size){
      notconverge=false;
      float maxerror=0;
      float sumerror=0;
      float aii, delta;

      //Red points
      for (ijk[2] = GridStartIndex[2]; ijk[2] <= GridEndIndex[2]; ijk[2]=ijk[2]+1){
	for (ijk[1] = GridStartIndex[2]; ijk[1] <= GridEndIndex[1]; ijk[1]=ijk[1]+1){
	  for (ijk[0] = GridStartIndex[2]; ijk[0] <= GridEndIndex[0]; ijk[0]=ijk[0]+2){
	    

	    index = GetIndex(ijk[0], ijk[1], ijk[2]);
	    activeIndex=getActiveIndex(ijk[0], ijk[1], ijk[2]);
	    
	    old=Phi_p[index];
	    
	    //For boundaries, we want dPhi=0
	    if (ijk[0]==GridStartIndex[0] || ijk[0]==GridStartIndex[0]+1) Phi_p[index-2*diffs[0]]=Phi_p[index];
	    else if (ijk[0]==GridEndIndex[0] || ijk[0]==GridEndIndex[0]-1)  Phi_p[index+2*diffs[0]]=Phi_p[index];
	    
	    if (ijk[1]==GridStartIndex[1] || ijk[1]==GridStartIndex[1]+1) Phi_p[index-2*diffs[1]]=Phi_p[index];
	    else if (ijk[1]==GridEndIndex[1] || ijk[1]==GridEndIndex[1]-1)  Phi_p[index+2*diffs[1]]=Phi_p[index];
	    
	    if (ijk[2]==GridStartIndex[2] || ijk[2]==GridStartIndex[2]+1) Phi_p[index-2*diffs[2]]=Phi_p[index];
	    else if (ijk[2]==GridEndIndex[2] || ijk[2]==GridEndIndex[2]-1)  Phi_p[index+2*diffs[2]]=Phi_p[index];
	    
	    
	    //the coefficient of A(i,i) diagonal terms
	    aii=(-2/dx2[0]-2/dx2[1]-2/dx2[2]);
	    
	    //adding the off diagonal correction terms
	    delta=0;
	    delta+=Phi_p[index+2*diffs[0]]/(4*dx2[0]);
	    delta+=Phi_p[index-2*diffs[0]]/(4*dx2[0]);
	    delta+=Phi_p[index+2*diffs[1]]/(4*dx2[1]);
	    delta+=Phi_p[index-2*diffs[1]]/(4*dx2[1]);
	    delta+=Phi_p[index+2*diffs[2]]/(4*dx2[2]);
	    delta+=Phi_p[index-2*diffs[2]]/(4*dx2[2]);
	    
	   //putting it all together
	    //Phi_p[index]=(1-w)*Phi_p[index]+w/aii*(divB_p[activeIndex]-delta);
	  

	
	    sumerror+=POW(old-Phi_p[index],2); 

	  }}}

      counter++;

      sumerror=POW(sumerror,0.5);
      if (sumerror>0.00001) notconverge=true;
	  
      if (debug&& counter%maxdim== maxdim-1) 
	 if (debug) printf("Div Cleaning %d Iterations, Error %f\n", 
		counter, sumerror);
      //  if (counter%2== 1) printf(stderr, "Div Cleaning %d Iterations, Error %f\n", counter, sumerror);



  }
  
  if (debug) printf(stderr, "Div Cleaning Takes %d Steps\n", counter);

#endif  
    
  return SUCCESS;

}


//Never Implemented FFT or Multigrid properly; Placeholders here for future work

int grid::PoissonSolverFFT(){
  return SUCCESS;
}


int grid::PoissonSolverMultigrid()
{
#ifdef NOUSE
  
  int ijk[3];
  int index, activeIndex;
  int dims[3]={GridDimension[0], GridDimension[1],GridDimension[2]};
  int size=GridDimension[0]*GridDimension[1]*GridDimension[2];
  int diff[3]={1, GridDimension[0], GridDimension[0]*GridDimension[1]};
 
  float* divB_plarger=new float [size];

   for (ijk[2] = 0; ijk[2] < GridDimension[2]; ijk[2]=ijk[2]+1){
      for (ijk[1] = 0; ijk[1] < GridDimension[1]; ijk[1]=ijk[1]+1){
	for (ijk[0] = 0; ijk[0] < GridDimension[0]; ijk[0]=ijk[0]+1){
	  divB_plarger[getIndex(ijk[0], ijk[1], ijk[2])]=0;
      }}}

  for (ijk[2] = 1; ijk[2] < GridDimension[2]-1; ijk[2]=ijk[2]+1){
    for (ijk[1] = 1; ijk[1] < GridDimension[1]-1; ijk[1]=ijk[1]+1){
      for (ijk[0] = 1; ijk[0] < GridDimension[0]-1; ijk[0]=ijk[0]+1){
	
	index=getIndex(ijk[0], ijk[1], ijk[2]);
	
	divB_plarger[index]=(
			   (BaryonField[iBx][index+diff[0]]-BaryonField[iBx][index-diff[0]])/2+
			   (BaryonField[iBy][index+diff[1]]-BaryonField[iBy][index-diff[1]])/2+
			   (BaryonField[iBz][index+diff[2]]-BaryonField[iBz][index-diff[2]])/2)*CellWidth[0][ijk[0]];
	
	
      }}}



  float tol_dim = (5.0e-4) * POW(0.1, 3-GridRank);
  tol_dim = max(sqrt(float(size))*1e-6, tol_dim);

  float norm = huge_number, mean = norm;

  MultigridSolver(divB_plarger, Phi_p, 3, dims, 
		  norm, mean, 0,
		  tol_dim, 1000);
#endif
  
  return SUCCESS;

}

// int grid::PrintToScreenBoundaries(float *field, char *display, int direction, int slice,
// 				   int check, float diffvalue)
// {
 
//   char* c= (char *) malloc(100 * sizeof(char));
//   sprintf(c, "BO_%"GOUTSYM".TNT",Time);

//   FILE *fptr=fopen(c, "w");
//   printf(fptr, "x\ty\tvalue\twidth\theight\n");
//   int xD[3]={GridDimension[0],GridDimension[1], GridDimension[2]};
//    int diffs[3]={1, xD[0], xD[1]*xD[0]}; 
 

//   bool fail=false;
//   int index;    int ijk[3];

//   if (check==0){fail=true;}

//   else if (check ==1){
//     //  printf(stderr, "Checking Grid for Badness Level %d\n", Level);
//     //int direction=0;

    
//     for (ijk[2] = GridStartIndex[2]; ijk[2] <= GridEndIndex[2]; ijk[2]=ijk[2]+1){
//       for (ijk[1] = GridStartIndex[1]; ijk[1] <= GridEndIndex[2]; ijk[1]=ijk[1]+1){
// 	for (ijk[0] = GridStartIndex[0]; ijk[0] <= GridEndIndex[2]; ijk[0]=ijk[0]+1){
// 	  index=ijk[0]+ijk[1]*xD[0]+ijk[2]*xD[1]*xD[0];	
// 	  if (fabs((float) field[index]-field[index-diffs[0]])>diffvalue ||
// 	      fabs((float) field[index]-field[index+diffs[0]])>diffvalue) fail=true;
// 	}}}}

//   if (fail){
//     printf(stderr, "\n\n*******Displaying Data********\n");
//     printf(stderr, display); 	  printf(stderr, "\n");
    
    
//     bool intertemp;
    
//     int ind1, ind2;
//     if (direction==0){ ind1=1; ind2=2;}
//     else if (direction==1){ ind1=0; ind2=2;}
//     else if (direction==2){ ind1=1; ind2=0;}
    
//      ijk[direction]=slice;
   
//      ijk[direction]=slice;
     
//      for (ijk[ind2] = 0; ijk[ind2] < GridDimension[ind2]; ijk[ind2]=ijk[ind2]+1){
//        for (ijk[ind1] = 0; ijk[ind1] < GridDimension[ind1]; ijk[ind1]=ijk[ind1]+1){

// 	 index = GetIndex(ijk[0], ijk[1], ijk[2]);
// 	 printf(stderr, "%g\t", field[index]);
// 	 printf(fptr, "%d\t%d\t%g\t1.0\t1.0\n", ijk[ind1],ijk[ind2], field[index]);
//        }
//       printf(stderr, "\n");
//     }
    
   
//   }
//   fclose(fptr);
//  return true;
// }

// int grid::PrintToScreenBoundaries(float *field, char *display, int direction, int slice){
//   PrintToScreenBoundaries(field, display, direction, (int) floor(GridDimension[0]/slice), 0, 0.0); 
//   return true;
// }


// int grid::PrintToScreenBoundaries(float *field, char *display){
//   PrintToScreenBoundaries(field, display, 2, (int) floor(GridDimension[0]/2), 0, 0.0); return true;
// }









int grid::PrintToScreenBoundaries(float *field, char *display, int direction, int slice,
				   int check, float diffvalue){

  
  // return SUCCESS;

  //if (!debug) return SUCCESS;

  
  //if (ProcessorNumber!=4) return SUCCESS;

//     if (GridLeftEdge[0]!=0.0 || GridLeftEdge[1]!=1.0){ 
//       //printf("NotGrid %g %g %g\n", GridLeftEdge[0], GridLeftEdge[1], GridLeftEdge[2]);
//       return SUCCESS;}

  if (ProcessorNumber != MyProcessorNumber) {
    printf("PrintToScreen wrong processor %d Proc != %d MyProc \n", ProcessorNumber, MyProcessorNumber);
    return SUCCESS;
  }



 //  char* c= (char *)malloc(100 * sizeof(char));
//   sprintf(c, "BO_%"GOUTSYM".TNT",Time);

//   FILE *fptr=fopen(c, "w");
//   fprintf(fptr, "x\ty\tvalue\twidth\theight\n");

  int xD[3]={GridDimension[0],GridDimension[1], GridDimension[2]};
   int diffs[3]={1, xD[0], xD[1]*xD[0]}; 
 

  bool fail=false;
  int index;    int ijk[3];

  if (check==0){fail=true;}

  else if (check ==1){
    //  printf(stdout, "Checking Grid for Badness Level %d\n", Level);
    //int direction=0;

    
    for (ijk[2] = GridStartIndex[2]; ijk[2] <= GridEndIndex[2]; ijk[2]=ijk[2]+1){
      for (ijk[1] = GridStartIndex[1]; ijk[1] <= GridEndIndex[2]; ijk[1]=ijk[1]+1){
	for (ijk[0] = GridStartIndex[0]; ijk[0] <= GridEndIndex[2]; ijk[0]=ijk[0]+1){
	  index=ijk[0]+ijk[1]*xD[0]+ijk[2]*xD[1]*xD[0];	
	//   if (abs((float) field[index]-field[index-diffs[0]])>diffvalue ||
// 	      abs((float) field[index]-field[index+diffs[0]])>diffvalue) fail=true;
	  if (field[index]<=0.0) {slice=ijk[direction]; fail=true;}
	}}}
    fail=true;
  }

   if (fail){
    printf("\n\n\n\n");
    printf("%s\n", display); 	  
    printf("Grid Edges %g %g %g\n", GridLeftEdge[0], GridLeftEdge[1], GridLeftEdge[2]);
    printf( "\n\n*******Processor # %d ********\n", ProcessorNumber);
    printf( "\n\n*******Displaying Data (Slice in %d on cell %d) (TopGrid %d)  ********\n", direction, slice, isTopGrid() );
    
    bool intertemp;
    
    int ind1, ind2;
    if (direction==0){ ind1=1; ind2=2;}
    else if (direction==1){ ind1=0; ind2=2;}
    else if (direction==2){ ind1=1; ind2=0;}
    
     ijk[direction]=slice;
    


  


     ijk[direction]=slice;
     
     for (ijk[ind2] = 0; ijk[ind2] < GridDimension[ind2]; ijk[ind2]=ijk[ind2]+1){
       for (ijk[ind1] = 0; ijk[ind1] < GridDimension[ind1]; ijk[ind1]=ijk[ind1]+1){

	 index=ijk[0]+ijk[1]*(GridDimension[0])+ ijk[2]*(GridDimension[0])*(GridDimension[1]);
	 //printf( "%2.1E \t", field[index]);
	 printf( "%.4g \t", field[index]);
	 //printf( "%d\t%d\t%g\t1.0\t1.0\n", ijk[ind1],ijk[ind2], field[index]);
       }
      printf( "\n");
     }
    
   
   }
  //fclose(fptr);
 
 return true;
}

int grid::PrintToScreenBoundaries(float *field, char *display, int direction, int slice){
  PrintToScreenBoundaries(field, display, direction, slice, 0, 0.0); 
  return true;
}


int grid::PrintToScreenBoundaries(float *field, char *display){
 if (!debug) return SUCCESS;

  //return SUCCESS;




  PrintToScreenBoundaries(field, display, 1, (int) floor(GridDimension[1]/2.0), 1, 0.0); 
  //PrintToScreenBoundaries(field, display, 1, NumberOfGhostZones, 0, 0.0);
  // PrintToScreenBoundaries(field, display, 0, NumberOfGhostZones, 0, 0.0);

  return SUCCESS;

}

int grid::PrintToScreenBoundaries(){
  //   PrintToScreenBoundaries(field, display, 2, (int) floor(GridDimension[1]/2.0), 0, 0.0); return true;
  //PrintToScreenBoundaries(field, display, 1, GridDimension[1]-1-NumberOfGhostZones, 0, 0.0);

  //if (!debug) return SUCCESS;
 
   if (GridLeftEdge[0]!=0.0 || GridLeftEdge[1]!=1.0 ||  GridLeftEdge[2]!=0.0){ 
      //printf("NotGrid %g %g %g\n", GridLeftEdge[0], GridLeftEdge[1], GridLeftEdge[2]);
      return SUCCESS;}

 
  for (int i=0; i< NumberOfBaryonFields; i++){
    printf("\n\n\n\n\n\n-------------Displaying %d (%d)\n", FieldType[i], i);
    PrintToScreenBoundaries(OldBaryonField[i], "old", 1, NumberOfGhostZones, 0, 0.0);
    PrintToScreenBoundaries(BaryonField[i], "new", 1, NumberOfGhostZones, 0, 0.0);
    PrintToScreenBoundaries(OldBaryonField[i], "old", 0, NumberOfGhostZones, 0, 0.0);
    PrintToScreenBoundaries(BaryonField[i], "new", 0, NumberOfGhostZones, 0, 0.0);
  }
  
  return SUCCESS;

}


int grid::PrintToScreenBoundaries(int field){


 
  int i=field;

  PrintToScreenBoundaries(BaryonField[i], "Density");


 return SUCCESS;
  
}
