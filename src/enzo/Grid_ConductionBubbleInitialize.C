////////////////////////////////////////////////////////////////////////////////
//
//  GRID CLASS
//
//  written by: Brian O'Shea
//  date:       March 2010
//  modified1:  
//
//  PURPOSE: 
//
//  RETURNS: FAIL or SUCCESS
//
////////////////////////////////////////////////////////////////////////////////
 
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

// Function prototypes
int GetUnits (float *DensityUnits, float *LengthUnits,
	      float *TemperatureUnits, float *TimeUnits,
	      float *VelocityUnits, double *MassUnits, FLOAT Time);

int FindField(int field, int farray[], int numfields);

double g, bunch_of_constants, dKdr, 
  dKdr_cgs, K_mid, T_mid, n_mid,r_mid, 
  r_mid_cgs, r_max, r_max_cgs;

double *rad,*nofr,*Tofr;
int ncells;

#define KEV_KELVIN 1.1604e+7
#define KPC_CGS 3.0857e+21
#define DEFAULT_MU 0.6

double dndr(double T, double n);
double dtdr(double n);

void get_dens_temp(void);


// Grid Initializer
int grid::ConductionBubbleInitialize (FLOAT BubbleRadius, int PulseType, float DeltaEntropy, 
				      float MidpointEntropy, float EntropyGradient,
				      float MidpointTemperature, FLOAT BubbleCenter[MAX_DIMENSION]) {

  if (debug) {
    printf("Entering ConductionBubbleInitialize\n");
    fflush(stdout);
  }

  if (ProcessorNumber != MyProcessorNumber) 
    return SUCCESS;

  FLOAT x,y,z, r2;

  int i,j,k;

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, MetalNum;
  float TemperatureUnits = 1.0, DensityUnits = 1.0, LengthUnits = 1.0;
  float VelocityUnits = 1.0, TimeUnits = 1.0;
  double MassUnits = 1.0;

  float delta, this_delta;
  
  delta = POW(DeltaEntropy, 0.6);

  FLOAT sig2 = BubbleRadius*BubbleRadius;

  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }

  int MetallicityField = FALSE;
  if ((MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields))
      != -1)
    MetallicityField = TRUE;
  else
    MetalNum = 0;

  // Get system of units
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, 
	       &TimeUnits, &VelocityUnits, &MassUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.");
  }

  g = fabs(UniformGravityConstant)*LengthUnits/(TimeUnits*TimeUnits);

  if(UniformGravity==0) g = 0.0;  // if gravity is off make sure it's zero

  dKdr = EntropyGradient;
  K_mid = MidpointEntropy;
  T_mid = MidpointTemperature;
  r_mid = 0.5*(DomainRightEdge[0] - DomainLeftEdge[0]);
  ncells = 1024;

  rad = new double[ncells];
  nofr = new double[ncells];
  Tofr = new double[ncells];

  dKdr_cgs = dKdr * KEV_KELVIN / KPC_CGS;
  r_mid_cgs = r_mid * LengthUnits;
  r_max_cgs = 2.0*r_mid_cgs;


  printf("g, UGC = %e %e\n",g,UniformGravityConstant);
  printf("dKdr / cgs  %e  %e\n",dKdr, dKdr_cgs);
  printf("K_mid:  %e\n", K_mid);
  printf("T_mid:  %e\n",T_mid);
  printf("r_mid_cgs, r_max_cgs:  %e %e\n", r_mid_cgs, r_max_cgs);
  printf("\n");


  get_dens_temp();

  for(i=0; i<ncells; i++){
    rad[i] /= LengthUnits;  // convert to enzo distance
    nofr[i] *= DEFAULT_MU * 1.67e-24 / DensityUnits;  // convert to enzo-unit density (from number density)
    Tofr[i] /= (TemperatureUnits*(Gamma-1.0)*DEFAULT_MU);  // convert from temp to internal energy

    //printf("*** %d  %e  %e  %e\n",i,rad[i],nofr[i],Tofr[i]);
  }
  fflush(stdout);

  int GridStart[] = {0, 0, 0}, GridEnd[] = {0, 0, 0};

  for (int dim = 0; dim<GridRank; dim++) {
    GridStart[dim] = 0;
    GridEnd[dim] = GridDimension[dim]-1;
  }

  //FLOAT sig2 = PulseWidth*PulseWidth;

  int ii, small_index;
  FLOAT smallest_d, celldist;

  // loop over grid and set pulse values
  for (k = GridStart[2]; k <= GridEnd[2]; k++) 
    for (j = GridStart[1]; j <= GridEnd[1]; j++) 
      for (i = GridStart[0]; i <= GridEnd[0]; i++) {

	/* Compute position */
	x=y=z=0.0;

	/* Find distance from center. */

	// radius squared: assume we always want to be at center of 
	// box
	x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	r2 = POW(x-BubbleCenter[0], 2.0);

	if(GridRank>1){
	  y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	  r2 += POW(y-BubbleCenter[1], 2.0);
	}

	if(GridRank>2){
	  z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
	  r2 += POW(z-BubbleCenter[2], 2.0);
	}

	celldist = POW(r2,0.5);

	smallest_d = huge_number;
	for(ii=0;ii<ncells;ii++){
	  if(fabs(x-rad[ii])<smallest_d){
	    smallest_d = fabs(x-rad[ii]);
	    small_index=ii;
	  }
	}

	/*
	printf("%d %d %d   smallest_d, small_index = %e %d  (%e)   %e/%e\n",i,j,k,
	       smallest_d,small_index,
	       nofr[small_index],Tofr[small_index],x);
	*/

	this_delta = 1.0;

	if(PulseType==1){  // top hat with radius of 

	  if(celldist <= BubbleRadius){
	    this_delta = delta;
	  }
	  
	} else if (PulseType==2){

	  if(celldist <= 5.0*BubbleRadius){
	    this_delta = 1.0 + exp(-1.0*r2/sig2/2.0)*(delta-1.0);
	  }

	} else if (PulseType==3){

	  this_delta = 1.0 + (delta-1.0)*(1.0 - tanh((10.*(celldist/BubbleRadius-1.0)))) / 2.0;

	} else {
	  ENZO_VFAIL("PulseType is not specified correctly: choose 1,2 or 3 (your val %d)\n", PulseType);

	}


	BaryonField[DensNum][ELT(i,j,k)] = nofr[small_index];

	BaryonField[DensNum][ELT(i,j,k)] /= this_delta;

	if(HydroMethod==Zeus_Hydro){  // ZEUS
	  BaryonField[TENum][ELT(i,j,k)] = Tofr[small_index];  // TE = gas energy

	  BaryonField[TENum][ELT(i,j,k)] *= this_delta;

	} else{ // PPM
	  
	  BaryonField[TENum][ELT(i,j,k)] = Tofr[small_index];  // TE = total energy energy, but velocity=0 here.

	  BaryonField[TENum][ELT(i,j,k)] *= this_delta;

	  if(DualEnergyFormalism){

	    BaryonField[GENum][ELT(i,j,k)] = Tofr[small_index];  // if DEF=1, need to separately set the gas internal energy.

	    BaryonField[GENum][ELT(i,j,k)] *= this_delta;

	  } // DEF
	}

	if(TestProblemData.UseMetallicityField>0 && MetalNum != FALSE){
	  if(celldist <= BubbleRadius){
	    BaryonField[MetalNum][ELT(i,j,k)] = 
	      BaryonField[DensNum][ELT(i,j,k)]*TestProblemData.MetallicityField_Fraction;
	  } else {
	    BaryonField[MetalNum][ELT(i,j,k)] = tiny_number;
	  }

	  if(i%32==0){
	    BaryonField[MetalNum][ELT(i,j,k)] = 
	      BaryonField[DensNum][ELT(i,j,k)]*TestProblemData.MetallicityField_Fraction;
	  }

	}
      } // for(i...)

  delete [] rad;
  delete [] nofr;
  delete [] Tofr;

  if (debug) {
    printf("Exiting ConductionBubbleInitialize\n");
    fflush(stdout);}
  return SUCCESS;
}

void get_dens_temp(void){

  /*
  printf("***********\n");
  printf("g, UGC = %e %e\n",g,UniformGravityConstant);
  printf("dKdr / cgs  %e  %e\n",dKdr, dKdr_cgs);
  printf("K_mid:  %e\n", K_mid);
  printf("T_mid:  %e\n",T_mid);
  printf("r_mid_cgs, r_max_cgs:  %e %e\n", r_mid_cgs, r_max_cgs);
  */
  

  n_mid = POW(T_mid/KEV_KELVIN/K_mid, 1.5);



  // g*mu*mp/kb;

  bunch_of_constants = g*DEFAULT_MU*(1.67e-24)/(1.38e-16);

  /*
  printf("n_mid, bunch_of_constants:  %e %e %e\n",n_mid,bunch_of_constants,g);
  printf("\n");
  printf("***********\n");
  */

  double r, dr;
  double this_n, this_t, last_n, last_t;
  int i,Nsteps=200;

  this_n = last_n = n_mid;
  this_t = last_t = T_mid;
  r=r_mid_cgs;

  dr = (r_max_cgs - r_mid_cgs) / double(ncells/2);

  double this_entropy, last_entropy=0.0,dkdr=0.0;

  rad[ncells/2] = r;
  nofr[ncells/2] = this_n;
  Tofr[ncells/2] = this_t;

  double k1n,k2n,k3n,k4n,k1t,k2t,k3t,k4t;


  for(i=0;i<ncells/2-1;i++){

    last_n = this_n;
    last_t = this_t;

    // f(x) = dy/dx

    // k1 = dx*f(x,y)
    k1n = dr*dndr(this_t, this_n);
    k1t = dr*dtdr(this_n);

    // k2 = dx*f(x+0.5dx,y+k1/2)
    k2n = dr*dndr(this_t+k1t/2.0, this_n+k1n/2.0);
    k2t = dr*dtdr(this_n+k1n/2.0);

    // k3 = dx*f(x+0.5dx,y+k2/2)
    k3n = dr*dndr(this_t+k2t/2.0, this_n+k2n/2.0);
    k3t = dr*dtdr(this_n+k2n/2.0);

    // k4 = dx*f(x+dx,y+k3)
    k4n = dr*dndr(this_t+k3t, this_n+k3n);
    k4t = dr*dtdr(this_n+k3n);

    this_n += (1.0/6.0)*(k1n + 2.0*k2n + 2.0*k3n + k4n);
    this_t += (1.0/6.0)*(k1t + 2.0*k2t + 2.0*k3t + k4t);

    r += dr;

    if(this_t < 1.0e+6){
      this_n = last_n;
      this_t = last_t;
      printf("***** WARNING: FIXED TEMPERATURE AND DENSITY *****\n");
    }

    //if(this_n <= 0.0) this_n = last_n;
    //if(this_t <= 0.0) this_t = last_t;

    rad[ncells/2+i+1] = r;
    nofr[ncells/2+i+1] = this_n;
    Tofr[ncells/2+i+1] = this_t;

    //printf("+++ %e %e %e %e (%d %d)\n",this_n, this_t, r, dr,i,ncells/2+1+i);    
  }

  this_n = last_n = n_mid;
  this_t = last_t = T_mid;
  r=r_mid_cgs;

  dr = (0.0 - r_mid_cgs) / double(ncells/2);

  //printf("+-+-+-+- %e %e %e %e\n",this_n, this_t, r, dr);



  for(i=0;i<=ncells/2-1;i++){

    last_n = this_n;
    last_t = this_t;

    // f(x) = dy/dx

    // k1 = dx*f(x,y)
    k1n = dr*dndr(this_t, this_n);
    k1t = dr*dtdr(this_n);

    // k2 = dx*f(x+0.5dx,y+k1/2)
    k2n = dr*dndr(this_t+k1t/2.0, this_n+k1n/2.0);
    k2t = dr*dtdr(this_n+k1n/2.0);

    // k3 = dx*f(x+0.5dx,y+k2/2)
    k3n = dr*dndr(this_t+k2t/2.0, this_n+k2n/2.0);
    k3t = dr*dtdr(this_n+k2n/2.0);

    // k4 = dx*f(x+dx,y+k3)
    k4n = dr*dndr(this_t+k3t, this_n+k3n);
    k4t = dr*dtdr(this_n+k3n);

    this_n += (1.0/6.0)*(k1n + 2.0*k2n + 2.0*k3n + k4n);
    this_t += (1.0/6.0)*(k1t + 2.0*k2t + 2.0*k3t + k4t);

    r += dr;

    rad[ncells/2-1-i] = r;
    nofr[ncells/2-1-i] = this_n;
    Tofr[ncells/2-1-i] = this_t;

    //printf("*** %e %e %e %e (%d %d)\n",this_n, this_t, r, dr,i,ncells/2-1-i);    
  }

  rad[0] = 0.0;


  last_entropy=0;
  for(i=0;i<ncells;i++){
    this_entropy = Tofr[i] / POW(nofr[i], 2.0/3.0);
    dkdr=(this_entropy-last_entropy)/(-dr);
    printf("%d r,n,t:  %e  %e  %e  K, dkdkr:  %e  %e\n",i,rad[i],nofr[i],Tofr[i],this_entropy,dkdr);

    last_entropy = this_entropy;
  }


  return;
}

double dndr(double T, double n){
  return -0.6*(n/T)*(bunch_of_constants + dKdr_cgs*POW(n,2.0/3.0));
}

double dtdr(double n){
  return -0.4*(bunch_of_constants -1.5*dKdr_cgs*POW(n,2.0/3.0));
}
