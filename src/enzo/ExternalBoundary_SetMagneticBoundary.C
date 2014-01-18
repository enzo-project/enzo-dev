/***********************************************************************
/
/  EZTERNAL BOUNDARY CLASS (Set magnetic field on grids boundary)
/
/  written by: David Collins
/  date:       2005
/  modified1:
/
/  PURPOSE: Similar to ExternalBoundary_SetExternalBoundary.
/           Setting boundary values for the magnetic field in a cell-by-cell basis 
/           is not possible, due to our update strategy.  
/           Therefore, to minimize computation in the loop, individual loops are used
/           for each type of BC.  This is different from the other boundary value routine
/           in enzo, ExternalBoundary_SetExternalBoundary.C
/          
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"


inline int indgen(int i,int j,int k,int nx,int ny,int nz){return i+(nx)*(j+(ny)*(k));}

int ExternalBoundary::SetMagneticBoundary(int FieldRank, int GridDims[], int GridOffset[],
                          int StartIndex[], int EndIndex[],
                          float *A, int FieldType)
{
  
  //Al
  //Ai means Active Index
  //Bi means Boundary Index
  //Add[] controls the extension of the grid size along a vector field.
  //      For instance, for calls to the By face centered field, Add = {0,1,0}
  //      which corrects indexing along the y faces, but not x or z.
  int nb = NumberOfGhostZones;
  int nx = GridDims[0] - 2*nb, ny = GridDims[1] - 2*nb, nz = GridDims[2] - 2*nb;
  int i,j,k,sign,dbgflag,Bi,Ai, field;
  int nxt = nx+2*nb, nyt = ny+2*nb, nzt=nz+2*nb;
  int is,js,ks,ie, je ,ke;
  bool lx = 1, rx = 0, lxz = 0, rxz = 0, ly=0,ry=0,lz=0,rz=0;
  bool verbose = 0;
  
  //The fields centered on the boundary may require special attention.
  //For instance, on the X face, Bx on the boundary, for reflecting conditions,  must be zero.
  //Also indexing slightly changes.
  
  int Add[3] ={0,0,0};
  if(FieldType == Bfield1) {Add[0] = 1; field = 0;}
  else if(FieldType == Bfield2) {Add[1] = 1; field = 1;}
  else if(FieldType == Bfield3) {Add[2] = 1; field = 2;}
  else{
    fprintf(stderr, "Non magnetic field passed to SetMagneticBoundary.\n");
    return FAIL;
    
  }
  
  nx += Add[0]; ny += Add[1]; nz += Add[2];
  nxt += Add[0]; nyt += Add[1]; nzt += Add[2];
  is = 0; js = 0; ks = 0; ie = nxt; je = nyt; ke = nzt;
    
  //  
  //left x
  //

  if( GridOffset[0] ==0 ) 
  switch(MagneticBoundaryType[field][0][0])
    {
    case inflow: 
      if( verbose ) fprintf(stderr,"inflow?\n");
      for(k=ks;k<ke;k++)
	for(j=js;j<je;j++)
	  for(i=0;i<nb+Add[0];i++)
	    {
	      Bi = indgen(i,j,k,nxt,nyt,nzt);
	      Ai = (j+GridOffset[1]) + (k+GridOffset[2])*MagneticBoundaryDims[field][1];
              A[Bi] = 0.0;
	    }
      break;
      
    case outflow:
      if( verbose ) fprintf(stderr,"Outflow?\n");
      for(k=ks;k<ke;k++)
	for(j=js;j<je;j++)
	  for(i=0;i<nb;i++)
	    {
	      Bi = indgen(i,j,k,nxt,nyt,nzt);
	      Ai = indgen(nb,j,k,nxt,nyt,nzt);
	      A[Bi] = A[Ai];
	    }
      break;
      
    case reflecting:
      if( verbose ) fprintf(stderr,"Reflecting?\n");
      if(Add[0] == 1 )
	{
	  sign = -1;
	  for(k=ks;k<ke;k++)
	    for(j=js;j<je;j++){
	      i = nb;
	      Bi = indgen(i,j,k,nxt,nyt,nzt);
	      A[Bi] = 0;
	    }
	}
      else
	{
	  sign = 1;
	}

      for(k=ks;k<ke;k++)
	for(j=js;j<je;j++)
	  for(i=0;i<nb;i++){
	    Ai = indgen(2*nb-i-1+Add[0],j,k,nxt,nyt,nzt);
	    Bi = indgen(i,j,k,nxt,nyt,nzt);
	    A[Bi] = sign*A[Ai];
	  }

	break;

    case periodic:
      if( verbose ) fprintf(stderr,"Periodic?\n");
      break;
      if( GridOffset[0]+GridDims[0] == BoundaryDimension[0]){
	for(k=ks;k<ke;k++)
	  for(j=js;j<je;j++)
	    for(i=0;i<nb+Add[0];i++){
	      Ai = indgen(nx+i-Add[0],j,k,nxt,nyt,nzt);
	      Bi = indgen(i,j,k,nxt,nyt,nzt);
	      A[Bi] = A[Ai];
	      
	    }
      }
      break;
      
    default:
    case BoundaryUndefined:
      break;
    }
  
  // Right X face BC's
  
  if (BoundaryDimension[0] > 1 && GridOffset[0]+GridDims[0] == BoundaryDimension[0]) 
    switch(MagneticBoundaryType[field][0][1])
      {
      case inflow:
	for(k=ks;k<ke;k++)
	  for(j=js;j<je;j++)
	    for(i=0;i<nb+Add[0];i++){
	      Bi = indgen(nb+nx+i-Add[0],j,k,nxt,nyt,nzt);
	      Ai = j+GridOffset[1] + (k+GridOffset[2])*MagneticBoundaryDims[field][1];
              A[Bi] = 0.0;
	    }      
	break;
	
      case outflow:
	if( verbose ) fprintf(stderr,"Right Outflow X\n");
	for(k=ks;k<ke;k++)
	  for(j=js;j<je;j++)
	    for(i=0;i<nb;i++){
	      Bi = indgen(nb+nx+Add[0]+i,j,k,nxt,nyt,nzt);
	      Ai = indgen(nb+nx-1+Add[0]  ,j,k,nxt,nyt,nzt);
	      A[Bi] = A[Ai];
	    }
	break;
	
      case reflecting:
	if(Add[0] == 1)
	  {
	    sign = -1;
	    for(k=ks;k<ke;k++)
	      for(j=js;j<je;j++){
		Bi = indgen(nx+nb-1,j,k,nxt,nyt,nzt);
		A[Bi] = 0;
	      }
	  }
	else 
	  {sign = 1;}
	
	for(k=ks;k<ke;k++)
	  for(j=js;j<je;j++)
	    for(i=0;i<nb;i++){
	      Ai = indgen(nx+nb-i-1-Add[0],j,k,nxt,nyt,nzt);
	      Bi = indgen(nx+nb+i,j,k,nxt,nyt,nzt);
	      A[Bi] = sign*A[Ai];
	    }
	
	break;
	
      case periodic:
	break;
	if( GridOffset[0] == 0 ){
	  for(k=ks;k<ke;k++)
	    for(j=js;j<je;j++)
	      for(i=0;i<nb;i++){
		Ai = indgen(i+nb+Add[0],j,k,nxt,nyt,nzt);
		Bi = indgen(nx+nb+i,j,k,nxt,nyt,nzt);
		A[Bi] = A[Ai];
	      }
	}
	
	break;
	

      default:
      case BoundaryUndefined:
	  if(verbose) fprintf(stderr,"Undefined Boundary\n");
	break;
	
      }
  
  // Left Y face
  //
  
  if (BoundaryDimension[1] > 1 && GridOffset[1] == 0) 
    switch(MagneticBoundaryType[field][1][0])
      {
      case inflow: 
	for(k=ks;k<ke;k++)
	  for(j=0;j<nb+Add[1];j++)
	    for(i=is;i<ie;i++){
	      Bi=indgen(i,j,k,nxt,nyt,nzt);
	      Ai = i+GridOffset[0] + (k+GridOffset[2])*MagneticBoundaryDims[field][0];
              A[Bi] = 0.0;
	    }
	break;
	
      case outflow:
	for(k=ks;k<ke;k++)
	  for(j=0;j<nb;j++)
	    for(i=is;i<ie;i++){
	      Ai = indgen(i,nb,k,nxt,nyt,nzt);
	      Bi = indgen(i,j,k,nxt,nyt,nzt);
	      A[Bi] = A[Ai];
	    }
	break;
	
      case reflecting:
	
	if(Add[1] == 1)
	  {
	    sign = -1;
	    
	    for(k=ks;k<ke;k++)
	      for(i=is;i<ie;i++){
		j = nb;
		Bi = indgen(i,j,k,nxt,nyt,nzt);
		A[Bi] = 0;
		//oi(i,j,k);
	      }
	  }  
	else 
	  {sign = 1;}
	
	for(k=ks;k<ke;k++)
	  for(j=0;j<nb;j++)
	    for(i=is;i<ie;i++){
	      Ai = indgen(i,2*nb-j-1+Add[1],k,nxt,nyt,nzt);
	      Bi = indgen(i,j,k,nxt,nyt,nzt);
	      A[Bi] = sign*A[Ai];
	      
	    }
	
	break;
	
      case periodic:
	break;
	if( GridOffset[1]+GridDims[1] == BoundaryDimension[1]) {
	  for(k=ks;k<ke;k++)
	    for(j=0;j<nb+Add[1];j++)
	      for(i=is;i<ie;i++){
		Ai = indgen(i,ny+j-Add[1],k,nxt,nyt,nzt);
		Bi = indgen(i,j,k,nxt,nyt,nzt);
		
		A[Bi] = A[Ai];
		
	      }
	}
	break;
	
      default:
      case BoundaryUndefined:
	break;
	
	
      }
  
  //
  // right y 
  
  if (BoundaryDimension[1] > 1 && GridOffset[1]+GridDims[1] == BoundaryDimension[1]) 
    switch(MagneticBoundaryType[field][1][1])
      {
      case inflow:
	for(k=ks;k<ke;k++)
	  for(j=0;j<nb+Add[1];j++)
	    for(i=is;i<ie;i++){  
	      Ai =  i+GridOffset[0] + (k+GridOffset[2])*MagneticBoundaryDims[field][0];
	      Bi = indgen(i,nb+ny+j-Add[1],k,nxt,nyt,nzt);
	      A[Bi] = 0.0;
	    }
	break;
	
      case outflow: 
	for(k=ks;k<ke;k++)
	  for(j=0;j<nb;j++)
	    for(i=is;i<ie;i++){
	      Ai = indgen(i,nb+ny-1+Add[1]  ,k,nxt,nyt,nzt);
	      Bi = indgen(i,nb+ny+Add[1]+j,k,nxt,nyt,nzt);
	      A[Bi] = A[Ai];
	    }
	break;
	
      case reflecting:
	
	if( Add[1] == 1)
	  {
	    sign = -1;
	    j = 0; 
	    for(k=ks;k<ke;k++)
	      for(i=is;i<ie;i++){	      
		Bi = indgen(i,nb+ny+j-1,k,nxt,nyt,nzt);
		A[Bi] = 0;
		//oi(i,nb+ny+j,k);
	      }
	  }
	else{sign = 1;}
	
	for(k=ks;k<ke;k++)
	  for(j=0;j<nb;j++)
	    for(i=is;i<ie;i++){
	      Ai = indgen(i,nb+ny-1-j-Add[1],k,nxt,nyt,nzt);
	      Bi = indgen(i,nb+ny+j,k,nxt,nyt,nzt);
	      A[Bi] = sign *A[Ai];
	      
	    }
	break;

      case periodic:
	break;
	if( GridOffset[1] == 0){
	  for(k=ks;k<ke;k++)
	    for(j=0;j<nb;j++)
	      for(i=is;i<ie;i++){
		Ai = indgen(i,j+nb+Add[1],k,nxt,nyt,nzt);
		Bi = indgen(i,ny+nb+j,k,nxt,nyt,nzt);
		A[Bi] = A[Ai];
	      }
	}
	break;
	
      default:
      case BoundaryUndefined:
	  if(verbose) fprintf(stderr,"Undefined Boundary\n");

	break;
	
      }
  
  //
  // left z 
  //
  
  if (BoundaryDimension[2] > 1 && GridOffset[2] == 0) 
    switch(MagneticBoundaryType[field][2][0])
      {
      case inflow:
	for(k=0;k<nb+Add[2];k++)
	  for(j=js;j<je;j++)
	    for(i=is;i<ie;i++){
	      Ai = i+GridOffset[0] + (j+GridOffset[1])*MagneticBoundaryDims[field][0];
	      Bi = indgen(i,j,k,nxt,nyt,nzt);
              A[Bi] = 0.0;
	    }	      
	break;
	
      case outflow:
	for(k=0;k<nb;k++)
	  for(j=js;j<je;j++)
	    for(i=is;i<ie;i++){
	      Ai = indgen(i,j,nb,nxt,nyt,nzt);
	      Bi = indgen(i,j,k ,nxt,nyt,nzt);
	      A[Bi] = A[Ai];
	    }	      
	
	break;
	
      case reflecting:
	if( Add[2] ==1 )
	  {
	    sign = -1;
	    k = nb;
	    for(j=js;j<je;j++)
	      for(i=is;i<ie;i++)
		{
		  Bi = indgen(i,j,k,nxt,nyt,nzt);
		  A[Bi] = 0;
		  //oi(i,j,k);
		}
	  }
      else{sign = 1;}

      for(k=0;k<nb;k++)
	for(j=js;j<je;j++)
	  for(i=is;i<ie;i++){
	    Ai = indgen(i,j,2*nb-k-1+Add[2],nxt,nyt,nzt);
	    Bi = indgen(i,j,k,nxt,nyt,nzt);
	    A[Bi] = sign*A[Ai];
	    //oi(i,j,2*nb-k-1-Add[2]);
	  }
      
      break;
    case periodic:
      break;
      if( GridOffset[2]+GridDims[2] == BoundaryDimension[2]) {
	for(k=0;k<nb+Add[2];k++)
	  for(j=js;j<je;j++)
	    for(i=is;i<ie;i++){
	      Ai = indgen(i,j,nz+k-Add[2],nxt,nyt,nzt);
	      Bi = indgen(i,j,k,nxt,nyt,nzt);
	      A[Bi] = A[Ai];
	  }
      }
      break;

    default:

    case BoundaryUndefined:
	if(verbose) fprintf(stderr,"Undefined Boundary\n");
      break;

    }
  //
  // right z 
  //
  
  if (BoundaryDimension[2] > 1 && GridOffset[2]+GridDims[2] == BoundaryDimension[2]) 
    switch(MagneticBoundaryType[field][2][1])
      {
      case inflow:
	for(k=0;k<nb+Add[2];k++)
	  for(j=js;j<je;j++)
	    for(i=is;i<ie;i++){
	      Ai = i+GridOffset[0] + (j+GridOffset[1])*MagneticBoundaryDims[field][0];
	      Bi = indgen(i,j,nb+nz+k-Add[2],nxt,nyt,nzt);
              A[Bi] = 0.0;
	    }
	break;
      case outflow:
	for(k=0;k<nb;k++)
	  for(j=js;j<je;j++)
	    for(i=is;i<ie;i++){
	      Ai = indgen(i,j,nb+nz-1+Add[2]  ,nxt,nyt,nzt);
	      Bi = indgen(i,j,nb+nz+Add[2]+k,nxt,nyt,nzt);
	      A[Bi] = A[Ai];
	      
	    }
	
	break;
	
      case reflecting:
	if(Add[2] == 1)
	  {
	    sign = -1;
	    k = 0;
	  for(j=js;j<je;j++)
	    for(i=is;i<ie;i++){
	      Bi = indgen(i,j,nb+nz-1,nxt,nyt,nzt);
	      A[Bi] = 0;
	      //oi(i,j,nb+nz+k);
	    }
	  }
	else{sign = 1;}
	
	for(k=0;k<nb;k++)
	  for(j=js;j<je;j++)
	    for(i=is;i<ie;i++){
	      Ai = indgen(i,j,nb+nz-k-1-Add[2],nxt,nyt,nzt);
	      Bi = indgen(i,j,nb+nz+k,nxt,nyt,nzt);
	      A[Bi] =sign*A[Ai];
	    }
	
	break;
      case periodic:
	break;
	if( GridOffset[2] == 0 ){
	  for(k=0;k<nb;k++)
	    for(j=js;j<je;j++)
	      for(i=is;i<ie;i++){
		Ai = indgen(i,j,k+nb+Add[2],nxt,nyt,nzt);
		Bi = indgen(i,j,nz+nb+k,nxt,nyt,nzt);
		A[Bi] = A[Ai];
		
	      }
	}
	break;
	
      default:
      case BoundaryUndefined:
	break;
	
      }
  
  return SUCCESS;
}
