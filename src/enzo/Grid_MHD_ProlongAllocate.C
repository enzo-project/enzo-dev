/***********************************************************************
/
/  GRID CLASS (ALLOCATES FIELDS FOR DERIVATIVES)
/
/  written by: David Collins
/  date:       2004-2013
/  modified1:
/
/  PURPOSE:  Allocates fields for use in the divergence free interpolation
/            needed for MHDCT.  
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
#include "fortran.def"
#include "Hierarchy.h"
#include "LevelHierarchy.h"

// This is a wrapper for the prolongation routine.  Since the prolongation algorithm must
// have the most resolved information for the grid, which may not always be contained in one
// other grid, we must first calculate the derivatives that will be used in mhd_interpolate,
// then do the actual interpolation.
// For historic reasons, this is done entirely in one routine, mhd_interpolate (called from
// grid::ProlongFineGrid)
//
// Note: MHDCT only works with RefineBy = 2, thus refinement factors are hard coded.

int grid::MHD_DCheck(int * ChildDim, char * mess){

  if( ProcessorNumber!=MyProcessorNumber)
    return SUCCESS;

  int rf[3] = {2,2,2},i;

  int dXsize = (ChildDim[0]/rf[0]+2)*(ChildDim[1]/rf[1]+1)*(ChildDim[2]/rf[2]+1);
  int dYsize = (ChildDim[0]/rf[0]+1)*(ChildDim[1]/rf[1]+2)*(ChildDim[2]/rf[2]+1);
  int dZsize = (ChildDim[0]/rf[0]+1)*(ChildDim[1]/rf[1]+1)*(ChildDim[2]/rf[2]+2);


  for(i=0; i<dXsize; i++)
    if( DyBx[i] != DyBx[i] ){
      fprintf(stderr,"DCheckFail: DyBx\n");
      return FAIL;
    }
  for(i=0; i<dXsize; i++)
    if( DzBx[i] !=DzBx[i] ){
      fprintf(stderr,"Dcheck Fail: DzBx\n");
      return FAIL;
    }

  for(i=0; i<dXsize; i++)
    if(DyzBx[i] != DyzBx[i] ){
      fprintf( stderr,"Dcheck Fail: 3\n");
      return FAIL;
    }

  
  for(i=0; i<dYsize; i++)
    if(DxBy[i] !=DxBy[i] ){
      fprintf(stderr,"Dcheck Fail: 5\n");
      return FAIL;
    }

  for(i=0; i<dYsize; i++){
    if( DzBy[i] != DzBy[i]){
      fprintf(stderr,"Dcheck Fail::6\n");
      return FAIL;
    }
    if( DxzBy[i] != DxzBy[i] ){
      fprintf(stderr,"Dcheck Fail: 7\n");
      return FAIL;
    }
  }
  for(i=0; i<dZsize; i++){
    if(DxBz[i] != DxBz[i] ){
      fprintf(stderr, "Dcheck Fail 8\n");
      return FAIL;
    }
    if(DyBz[i] != DyBz[i] ){
      fprintf(stderr,"Dcheck Fail9\n");
      return FAIL;
    }
    if(DxyBz[i] !=DxyBz[i] ){
      fprintf(stderr,"Dcheck Fail10\n");
      return FAIL;
    }
  }

  //fprintf(stderr," OK \n");
  return SUCCESS;
    
}


// Note: This routine assumes that the refinement is always by a factor of two, and thus it is hard-wired.
int grid::MHD_ProlongAllocate(int * ChildDim){

  //I recognize that this may be a bad plan, hard wiring the refinement in like this.
  //However, since this routine must necessarily have a refinement=2, I'm not worried.
  int rf[3] = {2,2,2},i;

  int dXsize = (ChildDim[0]/rf[0]+2)*(ChildDim[1]/rf[1]+1)*(ChildDim[2]/rf[2]+1);
  int dYsize = (ChildDim[0]/rf[0]+1)*(ChildDim[1]/rf[1]+2)*(ChildDim[2]/rf[2]+1);
  int dZsize = (ChildDim[0]/rf[0]+1)*(ChildDim[1]/rf[1]+1)*(ChildDim[2]/rf[2]+2);

  //delete the derivatives, if they exist.  
  MHD_ProlongFree();

  DyBx = new float[ dXsize ];
  for(i=0; i<dXsize; i++)
    DyBx[i] = 0;
  DzBx = new float[ dXsize ];
  for(i=0; i<dXsize; i++)
    DzBx[i] = 0;
  DyzBx = new float[ dXsize ];
  for(i=0; i<dXsize; i++)
    DyzBx[i] = 0;

  DBxFlag = new int[ dXsize ];
  for(i=0; i<dXsize; i++)
    DBxFlag[i] = 0;

  DxBy = new float[ dYsize ];
  for(i=0; i<dYsize; i++)
    DxBy[i] = 0.01;
  DzBy = new float[ dYsize ];
  for(i=0; i<dYsize; i++)
    DzBy[i] = 0;
  DxzBy = new float[ dYsize ];
  for(i=0; i<dYsize; i++)
    DxzBy[i] = 0;
  DByFlag=new int [dYsize];
  for(i=0; i<dYsize; i++)
    DByFlag[i] = 0;


  DxBz = new float[ dZsize ];
  for(i=0; i<dZsize; i++)
    DxBz[i] = 0;
  DyBz = new float[ dZsize ];
  for(i=0; i<dZsize; i++)
    DyBz[i] = 0;
  DxyBz = new float[ dZsize ];
  for(i=0; i<dZsize; i++)
    DxyBz[i] = 0;
  DBzFlag = new int [dZsize];
  for(i=0;i<dZsize; i++)
    DBzFlag[i]=0;


  return SUCCESS;
    
}

int grid::MHD_ProlongFree(){

  if( DyBx != NULL ){
    delete [] DyBx;
    DyBx=NULL;
  }
  if( DzBx != NULL ){
    delete [] DzBx;
    DzBx=NULL;
  }
  if( DyzBx != NULL ){
    delete [] DyzBx;
    DyzBx=NULL;
  }

  if( DBxFlag != NULL ){
    delete [] DBxFlag;
    DBxFlag = NULL;
  }

  if( DxBy != NULL ){
    delete [] DxBy;
    DxBy=NULL;
  }
  if( DzBy != NULL ){
    delete [] DzBy;
    DzBy=NULL;
  }
  if( DxzBy != NULL ){
    delete [] DxzBy;
    DxzBy=NULL;
  }
  if( DByFlag != NULL ){
    delete [] DByFlag;
    DByFlag = NULL;
  }


  if( DxBz != NULL ){
    delete [] DxBz;
    DxBz=NULL;
  }
  if( DyBz != NULL ){
    delete [] DyBz;
    DyBz=NULL;
  }
  if( DxyBz != NULL ){
    delete [] DxyBz;
    DxyBz=NULL;
  }
  if( DBzFlag != NULL ){
    delete [] DBzFlag;
    DBzFlag = NULL;
  }
  return SUCCESS;
}
