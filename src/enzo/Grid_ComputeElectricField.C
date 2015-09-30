
//
// Compute the Electric Field using the algorithm presented by GS05.
// In brief, this routine computes the Electric Field from fluxes returned from the solver.
// 
// Four options exist, and are controlled by MHD_CT_Method.  Three are computed here, one is not.
// MHD_CT_Method = 0: None (Bad Idea.)
//          = 1: Balsara CT.  (Default.)
//          = 2: Athena method 1 (\alpha.  Like Lax-Fredrichs) 
//               Also needs MHD_CT_Methodparam to be set, frequently to 0.1
//          = 3: Athena method 2 (directionally switched.) (Prefered.)
//          = 4: Balsara CT, Toth Formulation.  Original method for PPM-L.
//               (Not done in the routine.)
//
// Note that switching between these methods is done through coefficients in the 
// reconstruction. These are denoted As, Bs, Cs. (See below.)

// Note that in Enzo, the field called ElectricField is actually
// dT*E.  Putting the timestep here, instead of in the curl, allows
// for our MHD analogue of flux correction.  This is accounted for in this routine.

// I exploit the cyclic permutation symmetry of the MHD equations pretty heavily in this code:  only one 
// electric field is coded, but written in a generic enough way as to be looped over.  This
// is easier to write an maintain, at the expense of slightly more complex code.


// Variables:
// -------------------
//
// Cyclic variables: dimX, dimY, dimZ 
//     Index variables for the vector quantites (V,E,B)
//     These variables are cyclicly permuted relative to eachother.  
//
//     For the E_x computation, dimX = 0, dimY = 1, dimZ = 2.
//     For the E_y computation, dimX = 1, dimY = 2, dimZ = 0.
//     For the E_z computation, dimX = 2, dimY = 0, dimZ = 1.
//
// Field index variables:  These variables relate the position in the
//   algorithm to position in the 1d memory array.
//   There are 3 different array structures used, so 3 different
//   variables are used to index the memory array.
//    C1,2,3,4:  cell centered variables:  velocity and Centered
//               Magnetic Field
//    F1,2,3,4:  Fluxes: face centered.  These also have several
//               representations, based on the solver used.
//    B1,2,3,4:  Face centere Magnetic Field.  Not that the fluxes and
//               the face centered magnetic field use different data
//               structure shapes, so require different indicies.
//    The diagram below shows the spatial location of each index,
//    relative to the electric field being calculated (at the center)
//  
//  
//    ----------------- 
//    |       |       | 
//    |       |       | 
//    |  C2   F2  C1  |  //<dbg 3>
//    |       |       | 
//    |       |       | 
//    ---F3---E---F1---   
//    |       |       | 
//    |       |       | 
//    |  C4   F4  C3  | 
//    |       |       | 
//    |       |       | 
//    ----------------- 

// Offset Variables:
//   Db[dim][direction], Dc[dim][direction]

//    These are offsets in memory for the two 'directions' transverse
//    to 'dim'.  
//    For instance, refer to the above diagram.  The point C1, at {i,j,k} is
//    indexed by index=i+nx*(j+ny*k).  If 'E' is E_z, then C2 = {i-1,j,k}
//    and is indexed by index - Dc[2][0].  C4 = {i-1,j-1,k} and is
//    indexed by index - Dc[2][0] - Dc[2][1].
//
//    Dc is the set for the cell centered fields.  Its values are:
//        { {nx   , nx*ny}, (dy, dz)
//          {nx*ny, 1    }, (dz, dx)
//          {1    , ny   }, (dx, dy)
//
//    Db is the set for the Face Centered Magnetic field, but the offsets
//       relate to the curl operator, so the size of the data structure reflects
//       the differenced component.
//        { {nx   , nx*(ny+1)}, (dy Bz, dz By)
//          {nx*ny, 1        }, (dz Bx, dx Bz)
//          {1    , nx+1     }, (dx By, dy Bx) 

//   
//    
// Macros and Pointers:
//    Enzo uses irritatingly long names for things.  These pointers
//    make my life as a programmer easier.
//

// Reconstruction Coefficients.
//    As = First Order Terms.  For Method 0, only these terms are used. Always equal to 1.
//    Bs = Linear Correction.  Equal to 1 for method 1, switched based on fluid velocity for method 2.
//    Cs = Dissipative Lax-Fredrichs term.  Only non-zero for method 2.

// EndX, EndY, EndZ
//    Usually 1.
//    Electric fields are along the edges.  End[XYZ] denotes the
//    ammount to add (to the cell centered loop indices) along [XYZ].
//    For rank<3, the 'Uniform' dimension (Z for Rank=2, Y&Z for Rank
//    = 1) End[YZ] is set to zero.

#include "performance.h"
#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"


// function SSSe 
//    Solution Sensitive Switch.
//    Switches Reconstruction Coefficients based on fluid flow properties.
//    Bs[LR] is the [Left/Right] coefficient to use, as prescribed in
//    Gardiner & Stone 2005.

inline void SSSe(float * BsL, float * BsR, float WindDirection){       
  
  float Tolerance = 1e-10;

  if( WindDirection > Tolerance ){
    *BsR = 0; *BsL=2;
  }else if( WindDirection < -1*Tolerance){
    *BsR = 2; *BsL=0;
  }else{
    *BsR = 1; *BsL=1;
  }

}

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int grid::ComputeElectricField(float dT, float ** Fluxes){  

  //loop variables.
  int i,j,k,dimX, dimY, dimZ, flag;
  int EndX=1, EndY=1,EndZ=1;

  //Permutation variables.  
  int C1, C2, C3, C4;                   // Cell centered fields.
  int B1, B2, B3, B4;                   // Face centered Mag. Field.
  int F1, F2, F3, F4;                   // Face centered flux.
  int Edex;                             // Corner centered field being updated.
  
  //inverse timestep.
  float dTi = 1.0/dT;

  //Expansion factor at t^{n+1/2}
  //dB/dt - 1/a curl(E) = 0
  //provided B = B_{half comoving} = B_{comoving} * sqrt(a)
  // E_{half comoving} = (VxB_{comoving})*sqrt{a}
  // thus 1/sqrt{a} for both a terms
  FLOAT aNpHalf=1.0, dadtNpHalf=0.0, inv_sqrt_aNpHalf=1.0;
  if(ComovingCoordinates==1 ){
      CosmologyComputeExpansionFactor(Time+0.5*dtFixed, &aNpHalf, &dadtNpHalf);
      inv_sqrt_aNpHalf = 1./sqrt(aNpHalf); //For converting to Half-Comoving
  }

  // Cell Centered Offsets.  Used for cell centered quantities as well as fluxes.
  int Dc[3][2] =       
    {{GridDimension[0],GridDimension[0]*GridDimension[1]}, 
     {GridDimension[0]*GridDimension[1],1},
     {1,GridDimension[0]}};    
  
  // Face centered offsets.  Used for the face centered magnetic field
  int Db[3][2] =        
    {{GridDimension[0],GridDimension[0]*(GridDimension[1] + 1)}, 
     {(GridDimension[0]+1)*GridDimension[1],1},
     {1,GridDimension[0]+1}};    
  
  // size of CellCentered (BaryonField) array.
  // and PPML_B indicies.  These are cyclicly permuted, just like the
  //                       S.B[0,1] indicies, but the map to the Flux array is different,
  //                       so S.B doesn't work here. (In Athena, the flux 
  //                       array doesn't keep the Longitudinal flux
  //                       (i.e. F_x has no B_x term) while PPML does.
  //                       
  // Li_B indicies are the relevant elements in the flux array. 
  
  int size = GridDimension[0]*
    ((GridRank > 1 ) ? GridDimension[1] : 1 )*
    ((GridRank > 2 ) ? GridDimension[2] : 1 );
  
  int Li_B[3][2] = { {1,1},{1,0},{0,0}};
  //Pointers (and one macro) for brevity.

  int ElectricStart[MAX_DIMENSION], ElectricEnd[MAX_DIMENSION];

  for( int dim=GridRank; dim<MAX_DIMENSION; dim++){
    ElectricStart[dim] = GridStartIndex[dim];
    ElectricEnd[dim] = GridEndIndex[dim];
  }
  for(int dim=0;dim<GridRank;dim++){
    ElectricStart[dim] = 1;
    ElectricEnd[dim] = GridDimension[dim]-2;
  }

  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum, B1Num, B2Num, B3Num;
  IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, 
                             TENum, B1Num, B2Num, B3Num);

  float *Bc[3] = {BaryonField[B1Num], BaryonField[B2Num], BaryonField[B3Num]};
  float *Bf[3] = {MagneticField[0], MagneticField[1], MagneticField[2]};
  float *Vel[3] = {BaryonField[ Vel1Num ], BaryonField[ Vel2Num ], BaryonField[ Vel3Num ] };
#define Ec(index) (Bc[ dimY ][index]*Vel[ dimZ ][index] - Bc[ dimZ ][index]*Vel[ dimY ][index] )
  
   
   //Default reconstructon coefficents.
   float As[4] = {1.0,1.0,1.0,1.0};  
   float Bs[8] = {1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0};
   float Cs[4] = {1.0,1.0,1.0,1.0};

   //
   // Alterations for 2.5, 1.5d
   // These are really just fail-safes, any fields that actually acces
   // Offsets along the 'homogenous' axii are set to zero, to avoid
   // memory access problems.
   //

   if( GridRank < 3 ){
     Dc[0][1] = 0;
     Dc[1][0] = 0;
     Db[0][1] = 0;
     Db[1][0] = 0;
     EndZ = 0;
   }
   if( GridRank < 2){
     Dc[0][0] = 0;
     Dc[2][1] = 0;
     Db[0][0] = 0;
     Db[2][1] = 0;
     EndY = 0;
   }


   for(k=ElectricStart[2]; k<=ElectricEnd[2]+EndZ; k++)
   for(j=ElectricStart[1]; j<=ElectricEnd[1]+EndY; j++)
   for(i=ElectricStart[0]; i<=ElectricEnd[0]+EndX; i++){

     for( dimX=0; dimX< 3; dimX++){

       //Double check that only the necessary variables get looped over.

       if( (dimX == 0 && i == ElectricEnd[0] +1 )||
	   (dimX == 1 && j == ElectricEnd[1] +1 )||
	   (dimX == 2 && k == ElectricEnd[2] +1 ) ) 
	 continue;

       //the cyclic variables.
       dimY = (dimX == 0 ) ? 1 : (dimX == 1 ) ? 2: 0;
       dimZ = (dimX == 0 ) ? 2 : (dimX == 1 ) ? 0: 1;

       //The index variables.
       Edex = i + ElectricDims[dimX][0] * (j + ElectricDims[dimX][1]*k);
       C1 = i + GridDimension[0] * (j + GridDimension[1]*k);
       C2 = C1 - Dc[dimX][0];
       C3 = C1 - Dc[dimX][1];
       C4 = C1 - (Dc[dimX][0] + Dc[dimX][1]);

       B1 = i + MagneticDims[ dimZ ][0]*(j + MagneticDims[ dimZ ][1]*k);
       B2 = i + MagneticDims[ dimY ][0]*(j + MagneticDims[ dimY ][1]*k);
       B3 = B1 - Db[dimX][0];
       B4 = B2 - Db[dimX][1];

       switch( HydroMethod ){
       case MHD_Li:
	 F1 = B1 +               MagneticSize[dimZ] * Li_B[dimX][0];
	 F2 = B2 +               MagneticSize[dimY] * Li_B[dimX][1];
	 F3 = B1 - Db[dimX][0] + MagneticSize[dimZ] * Li_B[dimX][0];
	 F4 = B2 - Db[dimX][1] + MagneticSize[dimY] * Li_B[dimX][1];
	 break;
       default:
	 fprintf(stderr,"Athena_ComputeElectricField and incompatable Hydro Method %d\n", HydroMethod);
	 return FAIL;
       }

       //
       // Solution/Parameter Sensitive Switches:
       // The flags are for, in order,  -d/dy, +d/dy, -d/dz, +d/dz
       //     terms in the reconstruction.

       for(flag=0; flag<4;flag++){
	 As[flag]=1;
	 switch(MHD_CT_Method){
	 case CT_BalsaraSpicer:
	   Cs[flag]= 0;
	   Bs[2*flag]  = 0;
	   Bs[2*flag+1]= 0;
	   break;
	 case CT_Athena_LF:
	   Cs[flag]=1;
	   Bs[2*flag]  = 1;
	   Bs[2*flag+1]= 1;
	   break;
	 case CT_Athena_Switch:
	   Cs[flag]= 0;

	   //The Bs for MHD_ElectricRecon = 2 is below the loop-- it just 
	   //wasn't smooth to have it in this loop.
	   break;
	 }
       }//flag loop

       // Set the electric component to use based on wind direction.

       if( MHD_CT_Method == CT_Athena_Switch ){
	 // -d/dy
	 SSSe(Bs+1,Bs+0,Vel[dimZ][C1] + Vel[dimZ][C3]);
	 // +d/dy
	 SSSe(Bs+3,Bs+2,Vel[dimZ][C2] + Vel[dimZ][C4]);
	 // -d/dz
	 SSSe(Bs+5,Bs+4,Vel[dimY][C1] + Vel[dimY][C2]);
	 // +d/dz
	 SSSe(Bs+7,Bs+6,Vel[dimY][C3] + Vel[dimY][C4]);

       }
       
       //
       //The bulk.
       //


       
       ElectricField[dimX][Edex] = 0.25*( As[0]*Fluxes[ dimZ ][F1]
					  +As[1]*Fluxes[ dimZ ][F3]
					  -As[2]*Fluxes[ dimY ][F2]
					  -As[3]*Fluxes[ dimY ][F4] );

       ElectricField[dimX][Edex] +=0.125*(
					  
					  // dimY directed derivative
					  //Note the (-) on the
					  //fluxes[dimY]. 
					  //There's an excellent
					  //reason for this, but I
					  //don't remember what it is.

					 -Bs[0]*(Ec(C1) + Fluxes[dimY][F2])  
					 -Bs[1]*(Ec(C3) + Fluxes[dimY][F4]) 
					 +Bs[2]*(-Fluxes[dimY][F2] - Ec(C2) )
					 +Bs[3]*(-Fluxes[dimY][F4] - Ec(C4) )
					 
					 // dimZ derivative.
					 -Bs[4]*(Ec(C1) - Fluxes[dimZ][F1] )
					 -Bs[5]*(Ec(C2) - Fluxes[dimZ][F3])
					 +Bs[6]*(Fluxes[dimZ][F1]- Ec(C3) )
					 +Bs[7]*(Fluxes[dimZ][F3] - Ec(C4) )

					 );


       //The dissipative derivative correction.
       ElectricField[dimX][Edex] += 
	 0.125*CT_AthenaDissipation*CellWidth[0][0]*dTi*(
						 //Dissipation term for -dimY derivative
						 +Cs[0]*(Bc[dimY][C1] - Bf[dimY][B2]
							 -(Bc[dimY][C3] - Bf[dimY][B4]))
						 //for +dimY
						 -Cs[1]*(Bf[dimY][B2] - Bc[dimY][C2]
							 -(Bf[dimY][B4] - Bc[dimY][C4]))
						 //for -dimZ (the sign change comes out of the method.)
						 -Cs[2]*(Bc[dimZ][C1] - Bf[dimZ][B1] 
							 -(Bc[dimZ][C2]-Bf[dimZ][B3]) )
						 //for +dimZ
						 +Cs[3]*(Bf[dimZ][B1] - Bc[dimZ][C3]
							 -(Bf[dimZ][B3] -Bc[dimZ][C4]))
						 );
       
       // Finally, add the dT.  While conceptually it doesn't belong
       // attached to the electric field, its  here for AMR concerns.
       // (the flux correction)
       
       ElectricField[dimX][Edex] *= dT*inv_sqrt_aNpHalf;

     }//dim
   }//i,j,k
   

   return SUCCESS;
}



//  LocalWords:  MagneticField
