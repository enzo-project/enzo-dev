#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "pout.h"
#define ABS abs

//both declared below.

float blaststyle(int i,int j, int k, int InitStyle, float BlastCenterLocal[], float Radius);
void eigen( float d, float vx, float vy, float vz, 
	    float bx, float by, float bz, float eng, 
	    float right[][7]);
float AreaWrapper(int i,int j,int k,int dim, int rank, float dx ,float dy, float dz,float Normal[]);
float VolumeWrapper(int i,int j,int k, int rank, float dx ,float dy, float dz,float Normal[]);
float AreaBelow(float pA[], float pB[], float pC[], float pD[], float Line[]);
float Triangle(float Base, float Height);
float Quadralateral(float pA[], float pB[], float pC[], float pD[]);

int grid::MHDBlastInitializeGrid(float DensityA, float DensityB,
				 float EnergyA,  float EnergyB,
				 float VelocityA[], float VelocityB[],
				 float BA[], float BB[], 
				 float Radius, float MHDBlastCenter[], int LongDimension,
				 float PerturbAmplitude, int PerturbMethod, float PerturbWavelength[],
				 int InitStyle, float Normal[])
  
{

  //<dbg>
  //fprintf(stderr,"fuck %f %f\n",DensityA,DensityB);
  //return FAIL;
  //</dbg>
  
  //Every processor needs to know this for every grid,
  //WHETHER OR NOT IT HAS THE DATA.

  float Pi = 3.14159265;

#ifndef ATHENA
   int  EquationOfState = 0 ;
   float IsothermalSoundSpeed = 0;
#endif  

  if( EquationOfState == 0 ){
    NumberOfBaryonFields = 5;
    if(DualEnergyFormalism) NumberOfBaryonFields = 6;
    fprintf(stderr,"GridDim %d %d %d\n",GridDimension[0],GridDimension[1],GridDimension[2]);
    FieldType[0] = Density;
    FieldType[1] = TotalEnergy;
    FieldType[2] = Velocity1;
    FieldType[3] = Velocity2;
    FieldType[4] = Velocity3;
    if(DualEnergyFormalism) FieldType[5] = InternalEnergy;
  }else if( EquationOfState == 1){
    NumberOfBaryonFields = 4;
    if(DualEnergyFormalism){
      NumberOfBaryonFields = 6;
      fprintf(stderr," WARNING!!!! Dual Energy Formalism and Isothermal EOS not checked.\n");
    }
    FieldType[0] = Density;
    FieldType[1] = Velocity1;
    FieldType[2] = Velocity2; 
    FieldType[3] = Velocity3;
    if(DualEnergyFormalism) FieldType[4] = InternalEnergy;
  }
  
  
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;  

  //Parameters

  int i, j, k, field, size=1, index, index2, dim, which = 0;
  float value, Xp, Yp, Zp;  //Xp is easier to searc for than x.
  
  // Declare and initialize all fields.
  int InitStartIndex[3]={0,0,0};
  int InitEndIndex[3] = {GridDimension[0]-1,GridDimension[1]-1,GridDimension[2]-1};
  float fraction; //Linear combination fraction of FieldA and FieldB
  FLOAT BlastCenterLocal[3];
  
  this->AllocateGrids();
  
  //
  // Hack for tests of random forcing.
  //
  for (dim = 0; dim < GridRank; dim++){
    size *= GridDimension[dim];
  }
  
  if(RandomForcing == TRUE){
    for(dim=0;dim<GridRank;dim++){
      RandomForcingField[dim]= new float[size];
      for(i=0;i<size;i++)
	RandomForcingField[dim][i] = 1.0;
    }
  }
  
  
  //
  //  Center Definition.  BlastCenterLocal  is in DataStructure units
  //  MHDBlastCenter is in Physical Units.
  //
  //
  
  for(i=0;i<GridRank;i++){
    BlastCenterLocal[i] = MHDBlastCenter[i];
    BlastCenterLocal[i] -= CellLeftEdge[i][0];
    BlastCenterLocal[i] /= CellWidth[i][0];
  }
  for(i=GridRank; i<3; i++)
    BlastCenterLocal[i] = 0;
  
  Radius /= CellWidth[LongDimension][0];
  
  fprintf(stderr, "Center %f %f %f, Radius, %f\n", 
	  BlastCenterLocal[0], BlastCenterLocal[1], BlastCenterLocal[2], Radius);
  

  //Variable names.
  int Eeng, Eden, Ev[3], Egas; 
  if (this->IdentifyPhysicalQuantities(Eden, Egas, Ev[0], Ev[1], 
				       Ev[2], Eeng) == FAIL) {
    fprintf(stderr, "MHDBlastInitializeGrid: Error in IdentifyPhysicalQuantities.\n");
    return FAIL;
  }
  
  //For characteristic advection.  Right[field][wave]
  //8x = Square Wave.
  //7x = Sine Wave.
  //x:
  //0 left fast
  //1 left Alfven
  //2 left Slow
  //3 contact
  //4 right slow
  //5 right alfven
  //6 right fast
  float Right[7][7], Pos, Amp, InitialPressure;

  int B2num=1, B3num= 2, wave = 1, Map[7];
  for( i=0; i<7; i++)
    for( j=0; j<7; j++)
      Right[i][j] = 10*i+j;

  if( PerturbMethod >= 70 && PerturbMethod <= 90 ){

    switch( PerturbMethod ){
    case 80: case 70: wave = 0; break;
    case 81: case 71: wave = 1; break;
    case 82: case 72: wave = 2; break;
    case 83: case 73: wave = 3; break;
    case 84: case 74: wave = 4; break;
    case 85: case 75: wave = 5; break;
    case 86: case 76: wave = 6; break;
    }

    B2num = ( InitStyle == 1 ) ? 1 : ( InitStyle == 2 ) ? 2 : ( InitStyle == 3 ) ? 0 : -2;
    B3num = ( InitStyle == 1 ) ? 2 : ( InitStyle == 2 ) ? 0 : ( InitStyle == 3 ) ? 1 : -2;

    //Assuming cyclic permutation
    switch(InitStyle){
      //Map[ baryon field index ] = EigenVector index
      //Long axis get's the "x" component of the eigen vector.
    case 1:
      B2num = 1; B3num = 2;
      Map[ Ev[0] ] = 2; Map[ Ev[1] ] = 3; Map[ Ev[2] ] = 4;
      break;
    case 2:
      B2num = 2;B3num = 0;
      Map[ Ev[0] ] = 4; Map[ Ev[1] ] = 2; Map[ Ev[2] ] = 3;
      break;
    case 3:
      B2num = 0; B3num = 1;
      Map[ Ev[0] ] = 3; Map[ Ev[1] ] = 4; Map[ Ev[2] ] = 2;
      break;
    }
    Map[ Eden ] = 0;
    if( EquationOfState == 0 )
      Map[ Eeng ] = 1;
    if( B2num == -2 ){
      fprintf(stderr, "please choose InitStyle = 1,2,3 for perturbations along x,y,z.\n");
      return FAIL;
    }

    //Eigen value order is {rho, eng, vx,vy,vz,by,bz}
    //or the proper permutation for the current problem.
    eigen(DensityA,VelocityA[0],VelocityA[1],VelocityA[2], 
	  BA[InitStyle-1], BA[B2num], BA[B3num], EnergyA, Right);

    if( EquationOfState == 0 ){
      InitialPressure = (Gamma - 1) * (EnergyA - 0.5*DensityA*
				       (VelocityA[0]*VelocityA[0]+
					VelocityA[1]*VelocityA[1]+
					VelocityA[2]*VelocityA[2])
				       -0.5*(BA[0]*BA[0]+BA[1]*BA[1]+BA[2]*BA[2]));
    }else if( EquationOfState == 1 ){
      InitialPressure = IsothermalSoundSpeed * IsothermalSoundSpeed * DensityA;
    }

    fprintf(stderr,"Eigenvector: Pressure %f\n", InitialPressure);
    for( field=0; field<7; field++)
      fprintf(stderr, "EigenVector[%d][%d] %f \n", field, wave, Right[field][wave]);
    fprintf(stderr,"EigenVector: B2 %d B3 %d \n", B2num, B3num);
  }

  //
  //Set up BarryonField and Centered Magnetic fields.
  //Add perturbatoin if necessary.
  //

  for(k=InitStartIndex[2];k<=InitEndIndex[2]; k++)
    for(j=InitStartIndex[1];j<=InitEndIndex[1];j++)
      for(i=InitStartIndex[0];i<=InitEndIndex[0];i++){
	
	index = i+GridDimension[0]*(j+GridDimension[1]*k);
	//fprintf( stderr, " %d %d %d %f ", i,j,k, Radius);

	switch( InitStyle ) {
	  
	  //most cases dealt with in the 'default' value, by the blaststyle() routine.
	case 5:
	case 6:
	case 7:
	case 8:
	  // test for passing, 'n' shit.
	  
	  which = 2;
	  if(InitStyle == 5) value = 10*(i-DEFAULT_GHOST_ZONES) + 
			       100*(j-DEFAULT_GHOST_ZONES) + 
			       1000*(k-DEFAULT_GHOST_ZONES);
	  if(InitStyle == 6) value = i;
	  if(InitStyle == 7) value = j;
	  if(InitStyle == 8) value = k;
	  
	  BaryonField[0][index] = value +  0;
	  BaryonField[1][index] = value + 1;
	  BaryonField[2][index] = value + 2;
	  BaryonField[3][index] = value + 3;
	  if( EquationOfState == 0)
	    BaryonField[4][index] = value + 4;
	  if( MHD_Used == TRUE ){
	    CenteredB[0][index] = value + 5;
	    CenteredB[1][index] = value + 6;
	    CenteredB[2][index] = value + 7;
	  }
	  
	  break;  
	case 9:
	  which = 2;
	  BaryonField[ Eden ][index] = 0;
	  BaryonField[ Ev[0] ][index] = i - DEFAULT_GHOST_ZONES;
	  BaryonField[ Ev[1] ][index] = j - DEFAULT_GHOST_ZONES;
	  BaryonField[ Ev[2] ][index] = k - DEFAULT_GHOST_ZONES;;
	  break;

	  //Volume fraction for oblique grids
	case 10:
	  which = 0;
	  fraction = VolumeWrapper(i,j,k,GridRank,CellWidth[0][0], CellWidth[1][0], CellWidth[2][0],Normal);
	  if( fraction == -1 ){
	    fprintf(stderr," MHDBlastInitializeGrid: problem with Volume Wrapper.\n");
	    return FAIL;
	  }
	  break;
	default:
	  which = 0;
	  fraction = blaststyle(i,j,k,InitStyle, BlastCenterLocal, Radius);
	  if( fraction == -1 ){
	    fprintf(stderr,"MHDBlastInitializeGrid: problem with cell centered blast style\n");
	    return FAIL;
	  }
	  break;
	  
	  
	}//switch Init Style

	if( which == 0 ) {
	  BaryonField[ Eden ][index] = (1-fraction)*DensityA + fraction*DensityB;
	  if( EquationOfState == 0 ) BaryonField[Eeng][index] = 
				       (1-fraction)*EnergyA + fraction*EnergyB;
	  BaryonField[ Ev[0] ][index] = (1-fraction)*VelocityA[0]+fraction* VelocityB[0]; 
	  if( GridRank > 1 || MHD_Used == TRUE )
	  BaryonField[ Ev[1] ][index] = (1-fraction)*VelocityA[1]+fraction* VelocityB[1]; 
	  if( GridRank > 2 || MHD_Used == TRUE )
	  BaryonField[ Ev[2] ][index] = (1-fraction)*VelocityA[2]+fraction* VelocityB[2]; 
	  if( MHD_Used == TRUE ){
	    CenteredB[0][index] = (1-fraction) * BA[0] + fraction*BB[0];
	    CenteredB[1][index] = (1-fraction) * BA[1] + fraction*BB[1];
	    CenteredB[2][index] = (1-fraction) * BA[2] + fraction*BB[2];
	  }
	}else{
	  continue;
	}

	
	switch( PerturbMethod ){

	case 1:
	  BaryonField[ Ev[0] ][index]+= PerturbAmplitude* ((float)rand()/(float)(RAND_MAX) 
						     - .5);

	  BaryonField[ Ev[1] ][index]+= PerturbAmplitude* ((float)rand()/(float)(RAND_MAX) 
						     - .5);
	  break;

	  //Plane Symmetric
	case 2: 

	    if( k == InitStartIndex[2] ){
	    BaryonField[ Ev[0] ][index]+= PerturbAmplitude* ((float)rand()/(float)(RAND_MAX) 
						       - .5);
	    BaryonField[ Ev[1] ][index]+= PerturbAmplitude* ((float)rand()/(float)(RAND_MAX) 
						       - .5);
	  }else{
	    index2 = i+GridDimension[0]*(j+GridDimension[1]*InitStartIndex[2]);
	    BaryonField[ Ev[0] ][index] = BaryonField[ Ev[0] ][index2];
	    BaryonField[ Ev[1] ][index] = BaryonField[ Ev[1] ][index2];
	  }

	  break;

	case 3:

	  if( k == InitStartIndex[2] ){

	    int modi = i % (int) PerturbWavelength[0];
	    int modj = j % (int) PerturbWavelength[1];

	    if( modi == 0 && modj == 0 ){
	      BaryonField[ Ev[0] ][index]+= PerturbAmplitude* ((float)rand()/(float)(RAND_MAX) 
							 - .5);
	      BaryonField[ Ev[1] ][index]+= PerturbAmplitude* ((float)rand()/(float)(RAND_MAX) 
							 - .5);
	    }else{
	      index2 = i-modi+GridDimension[0]*(j-modj+GridDimension[1]*k);
	      BaryonField[ Ev[0] ][index] = BaryonField[ Ev[0] ][index2];
	      BaryonField[ Ev[1] ][index] = BaryonField[ Ev[1] ][index2];
	    }

	  }else{
	    index2 = i+GridDimension[0]*(j+GridDimension[1]*InitStartIndex[2]);
	    BaryonField[ Ev[0] ][index] = BaryonField[ Ev[0] ][index2];
	    BaryonField[ Ev[1] ][index] = BaryonField[ Ev[1] ][index2];
	  }
	  break;

	  //Plane Symmetric, energy preserving
	case 7: 
	  if( EquationOfState == 0 ){
	    if( MHD_Used == TRUE ){
	      BaryonField[ Eeng ][index] -= 0.5 * (BaryonField[ Eden ][index] *
					  (BaryonField[ Ev[0] ][index] * BaryonField[ Ev[0] ][index]+
					   BaryonField[ Ev[1] ][index] * BaryonField[ Ev[1] ][index] ) );
	    }else{
	      BaryonField[ Eeng ][index] -= 0.5 * (BaryonField[ Ev[0] ][index] * BaryonField[ Ev[0] ][index]+
					  BaryonField[ Ev[1] ][index] * BaryonField[ Ev[1] ][index] );
	    }
	  }
	    if( k == InitStartIndex[2] ){
	    BaryonField[ Ev[0] ][index]+= PerturbAmplitude* ((float)rand()/(float)(RAND_MAX) 
						       - .5);
	    BaryonField[ Ev[1] ][index]+= PerturbAmplitude* ((float)rand()/(float)(RAND_MAX) 
						       - .5);
	  }else{
	    index2 = i+GridDimension[0]*(j+GridDimension[1]*InitStartIndex[2]);
	    BaryonField[ Ev[0] ][index] = BaryonField[ Ev[0] ][index2];
	    BaryonField[ Ev[1] ][index] = BaryonField[ Ev[1] ][index2];
	  }

	  if( EquationOfState == 0 ){
	    if( MHD_Used == TRUE ){
	      BaryonField[ Eeng ][index] += 0.5 * (BaryonField[ Eden ][index] *
					  (BaryonField[ Ev[0] ][index] * BaryonField[ Ev[0] ][index]+
					   BaryonField[ Ev[1] ][index] * BaryonField[ Ev[1] ][index] ) );
	    }else{
	      BaryonField[ Eeng ][index] += 0.5 * (BaryonField[ Ev[0] ][index] * BaryonField[ Ev[0] ][index]+
						   BaryonField[ Ev[1] ][index] * BaryonField[ Ev[1] ][index] );
	    }
	  }
	  break;

	case 4:	
	case 5:
	case 6:

	  Xp = (float) (i-GridStartIndex[0])/GridDimension[0] *(GridRightEdge[0]-GridLeftEdge[0]) 
	    + GridLeftEdge[0] + 0.5*CellWidth[0][i];
	  Yp = (float) (j-GridStartIndex[1])/GridDimension[1] *(GridRightEdge[1]-GridLeftEdge[1])
	    + GridLeftEdge[1] + 0.5*CellWidth[1][j];
	  Zp = (float) (k-GridStartIndex[2])/GridDimension[2] *(GridRightEdge[2]-GridLeftEdge[2])
	    + GridLeftEdge[2] + 0.5*CellWidth[2][k];
	  
	  //BaryonField[2][index] = PerturbAmplitude*sin(2*Pi * Xp / PerturbWavelength[0]);
	  BaryonField[ Ev[1] ][index] += PerturbAmplitude*sin(2*Pi * Xp / PerturbWavelength[1]);
	  break;
	  
	  //Characteristics.  Sine wave
	case 70:
	case 71:
	case 72:
	case 73:
	case 74:
	case 75:
	case 76:

	  //characteristic advection

	  Xp = (float) (i-GridStartIndex[0])/(GridDimension[0]-2*DEFAULT_GHOST_ZONES)
	    *(GridRightEdge[0]-GridLeftEdge[0]) + GridLeftEdge[0] + 0.5*CellWidth[0][i];
	  Yp = (float) (j-GridStartIndex[1])/(GridDimension[1]-2*DEFAULT_GHOST_ZONES)
	    *(GridRightEdge[1]-GridLeftEdge[1]) + GridLeftEdge[1] + 0.5*CellWidth[1][i];
	  Zp = (float) (k-GridStartIndex[2])/(GridDimension[2]-2*DEFAULT_GHOST_ZONES)
	    *(GridRightEdge[2]-GridLeftEdge[2]) + GridLeftEdge[2] + 0.5*CellWidth[2][i];


	  Pos = ( InitStyle == 1 ) ? Xp : ( InitStyle == 2 ) ? Yp : ( InitStyle == 3)? Zp : 0;
	  Amp = sin(2*Pi * Pos);
	  //Amp = (which == 0)?-1:1;
	  
	  //Make the velocity into momentum
	  for(field=0;field<3;field++)
	    BaryonField[ Ev[field] ][index] *= BaryonField[ Eden ][index];
	  
	  
	  for(field=0;field<NumberOfBaryonFields;field++){
	    BaryonField[field][index] +=PerturbAmplitude*Right[ Map[field] ][wave]*Amp;
	    //fprintf(stderr, "field %d amp %f xp %f field %f\n", field, PerturbAmplitude, Pos,
	    //BaryonField[field][index]);
	  }
	  
	  CenteredB[B2num][index] += PerturbAmplitude*Right[5][wave]*Amp;
	  CenteredB[B3num][index] += PerturbAmplitude*Right[6][wave]*Amp;
	  
	  
	  //Make the momentum into velocity 
	  for(field=0;field<3;field++)
	    BaryonField[ Ev[field] ][index] /= BaryonField[ Eden ][index];

	  break;


	    
	  //Characteristics: Step Function
	  case 80:
	  case 81:
	  case 82:
	  case 83:
	  case 84:
	  case 85:
	  case 86:

	    Amp = 1- fraction * 2;
	    //if( which == 0 ) Amp = +1;
	    //if( which == 1 ) Amp = -1;
	    
	    //Make the velocity into momentum
	    for(field=0;field<3;field++)
	      BaryonField[ Ev[field] ][index] *= BaryonField[ Eden ][index];
	    
	    
	    for(field=0; field< NumberOfBaryonFields; field++){
	      BaryonField[field][index] +=  Amp*PerturbAmplitude*Right[ Map[field] ][wave];
	      
	    }
	    CenteredB[B2num][index] +=  Amp*PerturbAmplitude*Right[5][wave];
	    CenteredB[B3num][index] +=  Amp*PerturbAmplitude*Right[6][wave];
	    
	    //Make the momentum  into velocity
	    for(field=0;field<3;field++)
	      BaryonField[ Ev[field] ][index] /= BaryonField[ Eden ][index];
	    

	    
	    break;
	    


	default:
	  break;
	}
	
      }//end barryonfield initialize
  
  //
  // MagneticField Initialize.
  //


  if( MHD_Used == TRUE )
    for(field=0;field<3;field++){
      for(k=0; k<MagneticDims[field][2]; k++)
	for(j=0; j<MagneticDims[field][1]; j++)
	  for(i=0; i<MagneticDims[field][0];i++){

	    index = i+MagneticDims[field][0]*(j+MagneticDims[field][1]*k);
	    
	    switch( InitStyle ) {
	      
	    case 5:
	    case 6:
	    case 7:
	    case 8:
	      
	      which = 2;	      
	      
	      if(InitStyle == 5) value = 10*i + 100*j + 1000*k;
	      if(InitStyle == 6) value = i; 
	      if(InitStyle == 7) value = j;
	      if(InitStyle == 8) value = 
				   ( (field==0) ? 1.0
				     : ( (field==1) ? i*j
					 : ( (field==2) ? 0.0
					     : 3000) ) );
	      /*	    
			    ( (field==0) ? (1.0*j*k)/1000 :
			    ( (field==1) ? (1.0*i*k)/1000 :
			    ( (field==2) ? (1.0*i*j)/1000 : 3000) ) );
	      */
	      //MagneticField[field][index] = (1000000*(1+field) +10000*i+100*j+k);
	      MagneticField[field][index] = value +  0* field;

	    break;

	    case 9:
	      //nothing to be done for this.
	      break;
	    case 10:
	      which = 0;
	      fraction = AreaWrapper(i,j,k,field, GridRank, CellWidth[0][0], CellWidth[1][0], CellWidth[2][0],Normal);
	      if( fraction == -1 ){
		fprintf(stderr,"MHDBlastInitializeGrid: problem with Area Wrapper.\n");
		return FAIL;
	      }

	      break;
	      
	    default:
	      which = 0;
	      fraction = blaststyle(i,j,k,InitStyle, BlastCenterLocal, Radius);

	      if( fraction == -1 ){
		fprintf(stderr,"MHDBlastInitializeGrid: problem with Face Centered blast style.\n");
		return FAIL; 
	      }
	      break;
	      
	    }//Init Style switch
	    
	    
	    if( which == 0){
	      MagneticField[field][index] = (1-fraction)*BA[field] + fraction*BB[field];
	    }else continue;


	    switch( PerturbMethod ){
	      
	      
	      
	    case 1:
	    case 2: 
	    case 3:
	    case 4:	
	    case 5:
	    case 6:
	      //These only modify the velocity.  
	      break;
	      
	      //Characteristics.  Sine wave
	    case 70:
	    case 71:
	    case 72:
	    case 73:
	  case 74:
	    case 75:
	    case 76:
	      
	      //characteristic advection
	    
	      Xp = (float) (i-GridStartIndex[0])/(GridDimension[0]-2*DEFAULT_GHOST_ZONES)
		*(GridRightEdge[0]-GridLeftEdge[0]) + GridLeftEdge[0] + 0.5*CellWidth[0][i];
	      Yp = (float) (j-GridStartIndex[1])/(GridDimension[1]-2*DEFAULT_GHOST_ZONES)
		*(GridRightEdge[1]-GridLeftEdge[1]) + GridLeftEdge[1] + 0.5*CellWidth[1][i];
	      Zp = (float) (k-GridStartIndex[2])/(GridDimension[2]-2*DEFAULT_GHOST_ZONES)
	      *(GridRightEdge[2]-GridLeftEdge[2]) + GridLeftEdge[2] + 0.5*CellWidth[2][i];
	      
	      
	      Pos = ( InitStyle == 1 ) ? Xp : ( InitStyle == 2 ) ? Yp : ( InitStyle == 3)? Zp : 0;
	      Amp = sin(2*Pi * Pos);
	      //Amp = (which == 0)?-1:1;
	      
	      
	      if( field == B2num ) MagneticField[B2num][index] += PerturbAmplitude*Right[5][wave]*Amp;
	      if( field == B3num ) MagneticField[B3num][index] += PerturbAmplitude*Right[6][wave]*Amp;
	      
	      
	      //Characteristics: Step Function
	    case 80:
	    case 81:
	    case 82:
	    case 83:
	    case 84:
	    case 85:
	    case 86:

	      Amp = 1-2*fraction;
	      //if( which == 0 ) Amp = +1;
	      //if( which == 1 ) Amp = -1;
	      
	      if( field == B2num) MagneticField[B2num][index] += Amp*PerturbAmplitude*Right[5][wave];
	      if( field == B3num) MagneticField[B3num][index] += Amp*PerturbAmplitude*Right[6][wave];
	      
	      
	      break;
	      
	      
	      
	    default:
	      break;
	    }
	    
	    
	  }//i,j,k
    }//field
  
  
  if(0==1)
    for(field=0;field<3;field++)
      for(k=MHDeStartIndex[field][2]; k<=MHDeEndIndex[field][2]; k++)
	for(j=MHDeStartIndex[field][1]; j<=MHDeEndIndex[field][1]; j++)
	  for(i=MHDeStartIndex[field][0]; i<=MHDeEndIndex[field][0];i++){
	    

	    index = i+ ElectricDims[field][0] *(j+ ElectricDims[field][1]*k);
	    switch( InitStyle ){
	    case 3:
	      ElectricField[field][index] = field + 10*i + 100*j + 1000*k;
	      break;
	      
	    default:
	      
	      ElectricField[field][index] = 0.0;
	      break;
	    }//switch
	  }

#ifdef HAOXU
    
  if(DualEnergyFormalism )
    for(index=0;index<size;index++)
    BaryonField[ Egas ][index] =
       BaryonField[ Eeng ][index]/BaryonField[ Eden ][index]
      - 0.5*(BaryonField[ Ev[0] ][index]*BaryonField[ Ev[0] ][index]+BaryonField[ Ev[1] ][index]*BaryonField[ Ev[1] ][index]+BaryonField[ Ev[2] ][index]*BaryonField[ Ev[2] ][index])
      - 0.5*(CenteredB[0][index]*CenteredB[0][index]+CenteredB[1][index]*CenteredB[1][index]+CenteredB[2][index]*CenteredB[2][index])/BaryonField[ Eden ][index];

#endif 

  
  //  this->CenterMagneticField();

  
  //<dcc>
  if(WriteInThisF(99)==TRUE ){
      
    char basename[30];
    sprintf(basename,"data9900.grid");//, LevelCycleCount[level], level, grid);
    FILE * dccptr = fopen(basename,"a");
    if( this->WriteGrid(dccptr,basename,1) == FAIL ){
      fprintf(stderr,"Robert Anton Wilson says fuck.\n");
      return FAIL;
    }
    fclose(dccptr);

  }
  //</dcc>
  
  return SUCCESS;
  
}

float blaststyle(int i,int j, int k, int InitStyle, float BlastCenterLocal[], float Radius){
  int which;
  float iF ,jF,kF;
  iF= (float) i + 0.5;
  jF= (float) j + 0.5;
  kF= (float) k + 0.5;
  switch( InitStyle ){
  case 0:
    
    if( (iF - BlastCenterLocal[0] )*(iF - BlastCenterLocal[0] )
	+(jF - BlastCenterLocal[1] )*(jF - BlastCenterLocal[1] )
	+(kF - BlastCenterLocal[2] )*(kF - BlastCenterLocal[2] ) >=
	Radius*Radius){
      
      which = 0;
      
    }else{
      which = 1;
    }
    
    break;
    
  case 1:
    
    if( fabs(iF- BlastCenterLocal[0] ) >= Radius ) {
      which = 0;
    }else{
      which = 1;
    }     
    break;
    
    
  case 2:
    
    if( fabs(jF - BlastCenterLocal[1] ) >= Radius ) {
      which = 0;
    }else{
      which = 1;
    }     
    break;
    
  case 3:
    
    if( fabs(kF - BlastCenterLocal[2] ) >= Radius ) {
      which = 0;
    }else{
      which = 1;
    }     
    break;
    
  case 40:
    
    if( (kF - BlastCenterLocal[2] )*(kF - BlastCenterLocal[2] )
	+(jF - BlastCenterLocal[1] )*(jF - BlastCenterLocal[1] )
	>= Radius*Radius){
      which = 0;
    }else{
      which = 1;
    }     
    break;
    
  case 41:
    
    if( (kF - BlastCenterLocal[2] )*(kF - BlastCenterLocal[2] )
	+(iF - BlastCenterLocal[0] )*(iF - BlastCenterLocal[0] )
	>= Radius*Radius){
      which = 0;
    }else{
      which = 1;
    }     
    break;
    
  case 42:
    
    if( (iF - BlastCenterLocal[0] )*(iF - BlastCenterLocal[0] )
	+(jF - BlastCenterLocal[1] )*(jF - BlastCenterLocal[1] )
	>= Radius*Radius){
      which = 0;
    }else{
      which = 1;
    }     
    break;
    
  default:
    which = -1;
    fprintf(stderr,"MHDBlast: Invalid Init Style %d\n",InitStyle);
    break;

  }//switch

  return which;
}//blaststyle

void eigen( float d, float vx, float vy, float vz, 
	    float bx, float by, float bz, float eng, 
	    float right[][7])
	    //, float Speeds[])
{  
  //Eigen vectors taken from Ryu & Jones, ApJ 442:228-258, 1995
  //Normalization starting on p. 231.
  float p;
  float bp = 0.5*(bx*bx + by*by + bz*bz );
  float sqrt2 = sqrt(2.0);
  float sqrt2i= 1.0/sqrt2;
  float sqrtD = sqrt(d);
  float sqrtDi = 1.0/sqrtD;
  float sbx  = sign(bx);
  float og1  = 1.0/(Gamma - 1);

#ifndef ATHENA
   int  EquationOfState = 0 ;
   float IsothermalSoundSpeed = 0;
#endif

  if( EquationOfState == 0 ){
    p = (Gamma -1 ) * (eng - 0.5* d * (vx*vx + vy*vy + vz * vz ) - 0.5*(bx*bx + by*by + bz*bz));
  }else{
    p = IsothermalSoundSpeed*IsothermalSoundSpeed * d;
  }
  
  //compute wave speeds
  float aa;
    if( EquationOfState == 0 ){
    aa = sqrt( Gamma* p/d );
  }else{
    aa = IsothermalSoundSpeed;
  }
  float cs = sqrt( 0.5*( aa*aa + 2*bp/d - sqrt( pow( (aa*aa + 2*bp/d ),2) - 4* aa*aa*bx*bx/d ) ) );
  float ca = sqrt( bx*bx/d ); 
  float cf = sqrt( 0.5*( aa*aa + 2*bp/d + sqrt( pow( (aa*aa + 2*bp/d ),2) - 4* aa*aa*bx*bx/d ) ) );

  /*
    Speeds[0] =aa;
    Speeds[1] =cs;
    Speeds[2] =ca;
    Speeds[3] =cf;
  */
  //compute ancilary values
  //The normalization of alph_f may change. This normalization uses Ryu
  //& Jones, but Balsara may be more robust.
  float betay, betaz, alph_f, alph_s, bt;
  
  if( (by == 0.0) && (bz == 0.0) ){
    betay = sqrt2i;
    betaz = sqrt2i;
    alph_f = 1;
    alph_s = 1;
  }else{
    bt = sqrt( by*by + bz*bz );
    betay = by/bt;
    betaz = bz/bt;
    alph_f = sqrt( (cf*cf-ca*ca)/(cf*cf-cs*cs) );
    alph_s = sqrt( (cf*cf-aa*aa)/(cf*cf-cs*cs) );
  }
  
  //the vectors
  
  //fast, left
  right[0][0] = alph_f;
  right[1][0] = (EquationOfState == 1 ) ? 0 :
    alph_f*0.5*(vx*vx+vy*vy+vz*vz) + 
    alph_f*cf*cf*og1 - alph_f*cf*vx + alph_s*ca*sbx*(betay*vy + betaz*vz) + 
    (Gamma-2)*og1*alph_f*(cf*cf-aa*aa);
  right[2][0] = alph_f*(vx - cf);
  right[3][0] = alph_f*vy + alph_s*betay*ca*sbx;
  right[4][0] = alph_f*vz + alph_s*betaz*ca*sbx;
  right[5][0] = alph_s*betay*cf*sqrtDi;
  right[6][0] = alph_s*betaz*cf*sqrtDi;
  
  //alfven][left
  right[0][1] = 0;
  right[1][1] = (EquationOfState == 1 ) ? 0 :
    1*(betaz*vy - betay*vz)*sbx;
  right[2][1] = 0;
  right[3][1] =  1*betaz*sbx;
  right[4][1] = -1*betay*sbx;
  right[5][1] = betaz*sqrtDi;
  right[6][1] = -betay*sqrtDi;

  
  //slow,left
  right[0][2] = alph_s;
  right[1][2] =  (EquationOfState == 1 ) ? 0 :
    alph_s*0.5*(vx*vx+vy*vy+vz*vz) + 
    alph_s*cs*cs*og1 - alph_s*cs*vx - alph_f*aa*sbx*(betay*vy + betaz*vz) +
    (Gamma-2)*og1*alph_s*(cs*cs - aa*aa );
  right[2][2] = alph_s*(vx-cs);
  right[3][2] = alph_s*vy - alph_f*betay*aa*sbx;
  right[4][2] = alph_s*vz - alph_f*betaz*aa*sbx;
  right[5][2] = -alph_f*betay*aa*aa*sqrtDi/cf;
  right[6][2] = -alph_f*betaz*aa*aa*sqrtDi/cf;
  //right[4][2]*betaz/betay; Why was this here?
  
  //entropy (no entropy wave in isothermal MHD.)(Or hydro,for that matter)
  if(EquationOfState == 1 ){
    right[0][3] = 1;
    right[1][3] = 0.5*(vx*vx+vy*vy+vz*vz);
    right[2][3] = vx;
    right[3][3] = vy;
    right[4][3] = vz;
    right[5][3] = 0;
    right[6][3] = 0;
  }else{
    right[0][3] = 0;
    right[1][3] = 0;
    right[2][3] = 0;
    right[3][3] = 0;
    right[4][3] = 0;
    right[5][3] = 0;
    right[6][3] = 0;

  }
  
  //slow,right
  right[0][4] = alph_s;
  right[1][4] = (EquationOfState == 1 ) ? 0 :
    alph_s*0.5*(vx*vx+vy*vy+vz*vz) + 
    alph_s*cs*cs*og1 + alph_s*cs*vx + alph_f*aa*sbx*(betay*vy + betaz*vz) +
    (Gamma-2)*og1*alph_s*(cs*cs - aa*aa );
  right[2][4] = alph_s*(vx+cs);
  right[3][4] = alph_s*vy + alph_f*betay*aa*sbx;
  right[4][4] = alph_s*vz + alph_f*betaz*aa*sbx;
  right[5][4] = -alph_f*betay*aa*aa*sqrtDi/cf;
  right[6][4] = -alph_f*betaz*aa*aa*sqrtDi/cf;
  
  //alfven,right
  right[0][5] = 0;
  right[1][5] = (EquationOfState == 1 ) ? 0 :
    -1*(betaz*vy - betay*vz)*sbx;
  right[2][5] = 0;
  right[3][5] = -1*betaz*sbx;
  right[4][5] = +1*betay*sbx;
  right[5][5] = betaz*sqrtDi;
  right[6][5] = -betay*sqrtDi;


  //fast, right
  right[0][6] = alph_f;
  right[1][6] = (EquationOfState == 1 ) ? 0 :
    alph_f*0.5*(vx*vx+vy*vy+vz*vz) + 
    alph_f*cf*cf*og1 + alph_f*cf*vx - alph_s*ca*sbx*(betay*vy + betaz*vz) + 
    (Gamma-2)*og1*alph_f*(cf*cf-aa*aa);
  right[2][6] = alph_f*(vx + cf);
  right[3][6] = alph_f*vy - alph_s*betay*ca*sbx;
  right[4][6] = alph_f*vz - alph_s*betaz*ca*sbx;
  right[5][6] = alph_s*betay*cf*sqrtDi;
  right[6][6] = alph_s*betaz*cf*sqrtDi;

}

