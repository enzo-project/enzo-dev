
/***********************************************************************
/
/  GRID CLASS (CHECK THE VALIDITY OF THE MHD FIELDS AND MONITOR DIVB) 
/  modified1:
/
/  PURPOSE:  Computes DivB for the CT scheme.
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

float PreviousEnergy;

int grid::MHD_Diagnose(char * label)
{

  if( MyProcessorNumber != ProcessorNumber )
    return SUCCESS;

  // if Crash is true, MHD_Diagnose returns FAIL, and halts the function
  // of enzo on large divergences. If false, it only prints the warning.
  int Crash = FALSE;

  int size = GridDimension[0]*GridDimension[1]*GridDimension[2];
  float TotalEnergy = 0;
  float Kinetic = 0;
  float MagneticCentered = 0;
  float MagneticFace = 0;
  float Mass = 0;
  float GasEnergy = 0;
  float AbsDivB = 0;
  float TotalDivB = 0;
  float dens, v1, v2, v3, eng, b1, b2, b3, divergence;
  float MaxDivB = -1;
  float dx = 1, dy = 1, dz = 1;

  if( 0 == 0 ){
    dx = CellWidth[0][0];
    if ( GridRank > 1 )
      dy = CellWidth[1][0];
    if ( GridRank > 2 )
      dz = CellWidth[2][0];
  }

  int i,j,k, index;

  int max[3]={-1,-1,-1}, min[3] = {GridDimension[0],GridDimension[1],GridDimension[2]};
  int Convex = 1;
  float MaxTolerance = 1e-8;

  float *DivB = NULL;
  if( DivB == NULL) {
    DivB = new float[size];
    for(i=0;i<size;i++) DivB[i] = 0.0;
    
  }


  for(i=0;i<3;i++)
    if( MagneticField[i]==NULL ){
      fprintf(stderr, "MHD_Diag called with null MagField\n");
      return SUCCESS;
    }

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
    return FAIL;
  }
  
  
  //Compute total energy, gas, kinetic, and magnetic energies
  //as well as divergence.

  for(k=GridStartIndex[2];k<=GridEndIndex[2]; k++)
    for(j=GridStartIndex[1];j<=GridEndIndex[1];j++)
      for(i=GridStartIndex[0];i<=GridEndIndex[0];i++){
	
	index = i + GridDimension[0]*(j + GridDimension[1]*k);

	eng= 1.0 + ( (EquationOfState == 0 ) ? BaryonField[TENum][index] : 0.0 );
	
	v1 =  BaryonField[Vel1Num][index];
	v2 =  BaryonField[Vel2Num][index];
	v3 =  BaryonField[Vel3Num][index];
	dens = BaryonField[DensNum][index];
	b1 = CenteredB[0][index];
	b2 = CenteredB[1][index];
	b3 = CenteredB[2][index];
	
	TotalEnergy += eng;
	Kinetic += 0.5*dens*(v1*v1+v2*v2+v3*v3);
	MagneticCentered += 0.5*(b1*b1+b2*b2+b3*b3);
	Mass += dens;
	
	GasEnergy += eng - 0.5*dens*(v1*v1+v2*v2+v3*v3) - 0.5*(b1*b1+b2*b2+b3*b3);
	  
	divergence = 
	  (MagneticField[0][indexb1(i+1,j,k)] -MagneticField[0][indexb1(i,j,k)])/ (dx)+
	  ( (GridRank < 2 ) ? 0 : 
	    (MagneticField[1][indexb2(i,j+1,k)] -MagneticField[1][indexb2(i,j,k)])/ (dy )) +
	  ( (GridRank < 3 ) ? 0 : 
	    (MagneticField[2][indexb3(i,j,k+1)] -MagneticField[2][indexb3(i,j,k)])/ (dz )); 
	
	if( fabs(divergence) > MaxDivB ) MaxDivB = fabs(divergence);

	if( divergence != divergence ) {
	  fprintf(stderr,"Nan in divergence.  That's too bad...\n");
	  if( Crash == TRUE )
	    return FAIL;
	}
	DivB[index] = divergence;
	
	AbsDivB += fabs(divergence);
	TotalDivB += divergence;

	// Keep track of the structure of the divergence, and one level of its geometry (only outer bounds).

	if( fabs(divergence) > MaxTolerance ) {
	  if( i > max[0] ) max[0] = i;
	  if( i < min[0] ) min[0] = i;
	  if( j > max[1] ) max[1] = j;
	  if( j < min[1] ) min[1] = j;
	  if( k > max[2] ) max[2] = k;
	  if( k < min[2] ) min[2] = k;
	}else{
	  if( i <= max[0] && i >= min[0] &&
	      j <= max[1] && j >= min[1] &&
	      k <= max[2] && k >= min[2])
	    Convex = 0;
	}
	
      }

  //Compute current

  if( MaxDivB >= MaxTolerance || TRUE ){
    fprintf(stderr, " ++++ L1 norm             L0 norm           L1/GridSize        Max    %s\n",
	    label);
    fprintf(stderr, " ++++ %13.12e %13.12e %13.12e %13.12e\n", 
	    AbsDivB, TotalDivB, AbsDivB/size, MaxDivB);
    float IsThisNan = TotalDivB;
    if( IsThisNan != IsThisNan ) {
      fprintf(stderr,"000000000000000000000 Nan Tastic.  %13.12e \n", IsThisNan);
      return FAIL;
    }

  if( Crash == TRUE )
    return FAIL;
  else
    return SUCCESS;

  }

  return SUCCESS;

}
