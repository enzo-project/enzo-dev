//
//
// Takes the curl of the electric field, adds it to the magnetic field.
//
// MagneticField[xyz][ijk] -= Curl( ElectricField[xyz][ijk] )
// Input parameters:
//   Start and End are the CellCentered indicies to, uh... start and end.
//       Magnetic field on the faces of these cells are updated.
//       Currently, 'extra' faces are updated in Bz even for 2d (similar, By and 1d).
//   CurlOnly:  If CurlOnly = 1, then 
//       MagneticField[xyz][ijk] = Curl( ElectricField[xyz][ijk] )
//       instead of '-=' as is usual.  This is useful for taking the curl of vector 
//       potentials, or debugging.

// Local Variables:
//   Db[field][axis]: Appropriate offset in the 1d memory array for 'field' along 'axis.'
//                    Only the offsets directly relevant for a curl
//                    are stored. 
//   OK[3][2]:  A switch to eliminate components of the curl for rank!=3.    
//              Derivatives are switched out directly.
// 
//
//   Note that the methodology here for reducing dimensionality of the
//   problem is slightly different than that employed in the
//   CenterMagneticField routine, wherein the offsets are sent to zero
//   instead of individual terms being switched out.  May or may not
//   change, depending on my mood.  I can't see one way as better than
//   another as of now.

#include "performance.h"
#include "ErrorExceptions.h"
#include <math.h>
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "DebugTools.h"

//
// Method = 0: MagneticField = Curl( ElectricField )
// Method = 1: MagneticField = MagneticField - Curl( ElectricField )
// Method = 2: MagneticField = OldMagneticField - Curl(ElectricField)

//
// Method 0 is used for, say, generating a MagneticField from a Vector
//   Potential (stored in the poorly named but properly centered "ElectricField" for convinence.)
// Method 1 currently isnt' used.  
// Method 2 is what's actually done. MagneticField = OldMagnetciField + Curl( ElectricField )
//          This is done after ElectricField is updated from it's subgrids.
//          This techinique is equivalent to flux correction, but requires less machinery.
int grid::MHD_Curl(int * Start, int * End, int Method){
  

  // loop & index variables
  int i,j,k,dimX, dimY, dimZ;
  int Bdex, E1,E2,E3,E4;
  
  //Set up physical variables.
  //Note that dtUsed is 1 because dt was moved to the ElectricField, for Flux Correction.
  //It has been left here just in case...
  float dtUsed = 1;
  float dTdX[3] = {0,0,0};
  
  for( dimX=0;dimX<GridRank; dimX++ )
    dTdX[dimX] = (dtUsed*(CellWidth[dimX][0] != 0 )) ? (1.0/CellWidth[dimX][0]) : 0;
  
  // offset for the Electric field.
  int Db[3][2] = { { GridDimension[0]+1, (GridDimension[0]+1)*GridDimension[1] },
		   { GridDimension[0]*(GridDimension[1]+1), 1}, 
		   { 1, GridDimension[0] } };

  //For rank < 3, certain terms are omitted.  OK array keeps track of which
  int OK[3][2] = {{1,1},{1,1},{1,1}};

  if( GridRank < 3 ){ OK[0][1] = 0; OK[1][0] = 0;}
  if( GridRank < 2 ){ OK[0][0] = 0; OK[2][1] = 0;}

  for(dimX=0;dimX<3; dimX++){

    //in 1d, dBx/dt = 0.
    if( GridRank == 1 && dimX == 0 ) 
      continue;
    for(k=Start[2]; k<=End[2] + MHDAdd[2][dimX]; k++)
      for(j=Start[1]; j<=End[1] + MHDAdd[1][dimX]; j++)
	for(i=Start[0]; i<=End[0] + MHDAdd[0][dimX]; i++){
	  
	  
	  dimY = (dimX == 0 ) ? 1 : (dimX == 1 ) ? 2: 0;
	  dimZ = (dimX == 0 ) ? 2 : (dimX == 1 ) ? 0: 1;
	  
	  Bdex = i + MagneticDims[dimX][0]*(j+MagneticDims[dimX][1] * k);
	  E1   = i + ElectricDims[dimZ][0]*(j+ElectricDims[dimZ][1] * k);
	  E2   = E1 + Db[dimX][0];
	  E3   = i + ElectricDims[dimY][0]*(j+ElectricDims[dimY][1] * k);
	  E4   = E3 + Db[dimX][1];


	  switch( Method ){

	  case 0:
	    MagneticField[dimX][Bdex] =  
	      ( (  (OK[dimX][0] == 1 ) ? dTdX[dimY]*(ElectricField[dimZ][E2] - ElectricField[dimZ][E1]) : 0)
		-( (OK[dimX][1] == 1) ? dTdX[dimZ]*(ElectricField[dimY][E4]- ElectricField[dimY][E3]) : 0 ));
	    break;

	  case 1:
	    //fprintf(stderr,"curl 1\n");	    
	    MagneticField[dimX][Bdex] -=  
	      ( (  (OK[dimX][0] == 1 ) ? dTdX[dimY]*(ElectricField[dimZ][E2] - ElectricField[dimZ][E1]) : 0)
		-( (OK[dimX][1] == 1) ? dTdX[dimZ]*(ElectricField[dimY][E4]- ElectricField[dimY][E3]) : 0 ));
	    
	    break;
	  case 2:
	    MagneticField[dimX][Bdex] = OldMagneticField[dimX][Bdex]  -
	      ( (  (OK[dimX][0] == 1 ) ? dTdX[dimY]*(ElectricField[dimZ][E2] - ElectricField[dimZ][E1]) : 0)
		-( (OK[dimX][1] == 1) ? dTdX[dimZ]*(ElectricField[dimY][E4]- ElectricField[dimY][E3]) : 0 ));
	    break;
	  default:
	    ENZO_VFAIL(" Method = %"ISYM" isn't a valid argument to MHD_Curl.  Fix it.\n",Method)
	    break;
	  }//switch
	}//i,j,k
  }//dim
  
  return SUCCESS;
  
}






