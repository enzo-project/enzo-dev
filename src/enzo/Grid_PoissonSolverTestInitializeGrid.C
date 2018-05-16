#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"
#include "hydro_rk/EOS.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
double Gaussian(double cs);

int grid::PoissonSolverTestInitializeGrid(int TestType, float GeometryControl)
{
  /* declarations */

  int dim, i, j, k, m, sphere,B1, B2, B3, phip_num;

  int TE, IE;
  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  FieldType[NumberOfBaryonFields++] = Velocity2;
  FieldType[NumberOfBaryonFields++] = Velocity3;
  FieldType[TE=NumberOfBaryonFields++] = TotalEnergy;
  if (DualEnergyFormalism) {
    FieldType[IE=NumberOfBaryonFields++] = InternalEnergy;
  }
  FieldType[B1=NumberOfBaryonFields++] = Bfield1;
  FieldType[B2=NumberOfBaryonFields++] = Bfield2;
  FieldType[B3=NumberOfBaryonFields++] = Bfield3;
  FieldType[NumberOfBaryonFields++] = PhiField;

  if(UseDivergenceCleaning){
    FieldType[phip_num=NumberOfBaryonFields++] = Phi_pField;
  }

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }

  float rhou, lenu, tempu, tu, velu;
  
  GetUnits(&rhou, &lenu, &tempu, &tu, &velu, Time);


  this->AllocateGrids();


    

  int p,n=0;
  double x,y,z, magnitude;

  float sized[3]={DomainRightEdge[0]-DomainLeftEdge[0],
		  DomainRightEdge[1]-DomainLeftEdge[1],
		  DomainRightEdge[2]-DomainLeftEdge[2]};
  
  for (k = 0; k < GridDimension[2]; k++) {
    for (j = 0; j < GridDimension[1]; j++) {
      for (i = 0; i < GridDimension[0]; i++, n++) { 
	
	for (p = 0; p < NumberOfBaryonFields; p++){
	  
	  x=(CellLeftEdge[0][i]+0.5*CellWidth[0][i]-sized[0]/2);
	  y=(CellLeftEdge[1][j]+0.5*CellWidth[1][j]-sized[1]/2);
	  z=(CellLeftEdge[2][k]+0.5*CellWidth[2][k]-sized[2]/2);

	  //actually square of magnitude... misnomer to catch those who aren't paying attention!
	  magnitude= x*x+y*y+z*z;
	  
	  if (TestType==0){ //monopole
	    if (p==B1) BaryonField[p][n]=GeometryControl*x/pow(magnitude, 1.5)+1;
	    else if (p==B2) BaryonField[p][n]=GeometryControl*y/pow(magnitude, 1.5)+1;
	    else if (p==B3) BaryonField[p][n]=GeometryControl*z/pow(magnitude, 1.5)+1;
	    else if (p==0) BaryonField[p][n]=1;
	    else BaryonField[p][n]=0;
	  }
	  else  if (TestType==1){ //sphere

	    if (pow(magnitude, 0.5)<GeometryControl){
	      if (p==B1) BaryonField[p][n]=-1;
	      else if (p==0) BaryonField[p][n]=1;
	      else BaryonField[p][n]=0;
	    }
	    else if (p==0) BaryonField[p][n]=1;
	    else BaryonField[p][n]=0;
	  }

	  else  if (TestType==2){ //gaussian


	    if (p==B1 || p==B2 || p==B3) BaryonField[p][n]=pow(2.71828183, -1*pow(magnitude,2)/GeometryControl);
	    else if (p==0) BaryonField[p][n]=1;
	    else BaryonField[p][n]=0;
   BaryonField[TE][n]=1000;
	  }

	  else if (TestType==3){ //constant divB
	    if (p==B1) BaryonField[p][n]=-x*GeometryControl;
	    else if (p==0) BaryonField[p][n]=1;
	    else BaryonField[p][n]=0.0;
	    BaryonField[TE][n]=1;
	    BaryonField[phip_num][n]=0.0;
	  }

	  else if (TestType==4){ //constant B
	    if (p==B1) BaryonField[p][n]=GeometryControl;
	    else if (p==0) BaryonField[p][n]=1;
	    else BaryonField[p][n]=0;
	  }
	  else if (TestType==5){ //projection greens function

	    if (pow(magnitude, 0.5)<CellWidth[0][k]){
	      if (p==B2) BaryonField[p][n]=-1*sign(x);
	      else if (p==0) BaryonField[p][n]=1;
	      else BaryonField[p][n]=0;
	   
	    }
	    else if (p==0) BaryonField[p][n]=1;
	    else BaryonField[p][n]=0;
	  }
	  else  if (TestType==6){ //saw function
	    if (p==B1 && fabs(x)<0.30) BaryonField[p][n]=fabs(x)*GeometryControl;
	    else if (p==0) BaryonField[p][n]=1;
	    else BaryonField[p][n]=0;
	  }
	   






	  //---------------------------------------------------



	  else if(TestType==7){// Magnetic Rotor
	  
	    FLOAT r;
	    FLOAT r0 = 0.1, r1 = 0.115, v0=2.0;
	    float pres = 1.0, eint, h, cs, dpdrho, dpde, etot;
	    float rho0 = 10.0, rho1 = 1.0, Bx0 = 1.41047;
	    //float Bx0= 5.0;

	    /* Compute position */
	    int igrid = GetIndex(i,j,k);

	    x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	    y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];

	    r = sqrt(pow(x-0.5,2)+pow(y-0.5,2));

	    if (r < r0) {
	      BaryonField[iden][igrid] = rho0;
	      BaryonField[ivx ][igrid] = -v0*(y-0.5)/r0;
	      BaryonField[ivy ][igrid] = v0*(x-0.5)/r0;
	      BaryonField[ivz ][igrid] = 0.0;
	      EOS(pres, rho0, eint, h, cs, dpdrho, dpde, 0, 1);
	      etot = eint + 0.5*(pow(BaryonField[ivx][igrid],2)+pow(BaryonField[ivy][igrid],2))
		+ 0.5*Bx0*Bx0/rho0;
	      BaryonField[ietot][igrid] = etot;
	      if (DualEnergyFormalism) {
		BaryonField[ieint][igrid] = eint;
	      }
	      BaryonField[iBx][igrid] = Bx0;
	      BaryonField[iBy][igrid] = 0.0;
	      BaryonField[iBz][igrid] = 0.0;
	      BaryonField[iPhi][igrid] = 0.0;
	    }
	    else if (r < r1) {
	      FLOAT f = (r1-r)/(r1-r0);
	      BaryonField[iden][igrid] = 1.0+9*f;
	      BaryonField[ivx ][igrid] = -f*v0*(y-0.5)/r0;
	      BaryonField[ivy ][igrid] = f*v0*(x-0.5)/r0;
	      BaryonField[ivz ][igrid] = 0.0;
	      EOS(pres, rho0, eint, h, cs, dpdrho, dpde, 0, 1);
	      etot = eint + 0.5*(pow(BaryonField[ivx][igrid],2)+pow(BaryonField[ivy][igrid],2))+
		0.5*Bx0*Bx0/(1.0+9*f);
	      BaryonField[ietot][igrid] = etot;
	      if (DualEnergyFormalism) {
		BaryonField[ieint][igrid] = eint;
	      }
	      BaryonField[iBx][igrid] = Bx0;
	      BaryonField[iBy][igrid] = 0.0;
	      BaryonField[iBz][igrid] = 0.0;
	      BaryonField[iPhi][igrid] = 0.0;
	    }
	    else {
	      BaryonField[iden][igrid] = rho1;
	      BaryonField[ivx ][igrid] = 0.0;
	      BaryonField[ivy ][igrid] = 0.0;
	      BaryonField[ivz ][igrid] = 0.0;
	      EOS(pres, rho0, eint, h, cs, dpdrho, dpde, 0, 1);
	      etot = eint + 0.5*Bx0*Bx0/rho1;
	      BaryonField[ietot][igrid] = etot;
	      if (DualEnergyFormalism) {
		BaryonField[ieint][igrid] = eint;
	      }
	      BaryonField[iBx][igrid] = Bx0;
	      BaryonField[iBy][igrid] = 0.0;
	      BaryonField[iBz][igrid] = 0.0;
	      BaryonField[iPhi][igrid] = 0.0;
	    }
 
	  }


// 	  else if(TestType==8){// RT problem
//  /* transform pressure to total energy */
//   float etotl, etotu, v2, B2;
//   v2 = vxl * vxl + vyl * vyl;
//   B2 = Bxl * Bxl + Byl * Byl;
//   etotl = pl / ((Gamma-1.0)*rhol) + 0.5*v2 + 0.5*B2/rhol;

//   v2 = vxu * vxu + vyu * vyu;
//   B2 = Bxu * Bxu + Byu * Byu;
//   etotu = pu / ((Gamma-1.0)*rhou) + 0.5*v2 + 0.5*B2/rhou;
	  

// 	    float pres, eintl, eintu, h, cs, dpdrho, dpde;
	  
// 	    /* Compute position */
// 	    index= i + j*GridDimension[0];
	    
// 	    x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
// 	    y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
// 	    if (y <= 0.75) {
	      
// 	      // Rayleigh-Taylor problem, calculate pressure from hydro equilibrium
// 	      float g = ConstantAcceleration[1];
// 	      pres = pl+g*rhol*(y-0.75);
// 	      EOS(pres, rhol, eintl, h, cs, dpdrho, dpde, 0, 1);
// 	      // impose mode perturbation
// 	      vyl = 0.01 * (1.0+cos(4.0*Pi*(x-0.25))) * (1.0+cos(3.0*Pi*(y-0.75))) * 0.25;
// 	      etotl = eintl + 0.5*(vxl*vxl + vyl*vyl) + 0.5*(Bxl*Bxl+Byl*Byl)/rhol;
	      
	      
// 	      BaryonField[iden ][igrid] = rhol;
// 	      BaryonField[ivx  ][igrid] = vxl;
// 	      BaryonField[ivy  ][igrid] = vyl;
// 	      BaryonField[ivz  ][igrid] = 0.0;
	      
	      
// 	      BaryonField[ietot][igrid] = etotl;
// 	      if (DualEnergyFormalism) {
// 		BaryonField[ieint][igrid] = pl / ((Gamma-1.0)*rhol);
// 	      }
// 	      BaryonField[iBx  ][igrid] = Bxl;
// 	      BaryonField[iBy  ][igrid] = Byl;
// 	      BaryonField[iBz  ][igrid] = 0.0;
// 	      BaryonField[iPhi ][igrid] = 0.0;
// 	    } else {
	      
// 	      // Rayleigh-Taylor problem, calculate pressure from hydro equilibrium
// 	      float g = ConstantAcceleration[1];
// 	      pres = pu+g*rhou*(y-0.75);
// 	      EOS(pres, rhou, eintu, h, cs, dpdrho, dpde, 0, 1);
// 	      // impose mode perturbation
// 	      vyu = 0.01 * (1.0+cos(4.0*Pi*(x-0.25))) * (1.0+cos(3.0*Pi*(y-0.75))) * 0.25;
// 	      etotu = eintu + 0.5*(vxu*vxu + vyu*vyu) + 0.5*(Bxu*Bxu+Byu*Byu)/rhou;
	      
// 	      BaryonField[iden ][igrid] = rhou;
// 	      BaryonField[ivx  ][igrid] = vxu;
// 	      BaryonField[ivy  ][igrid] = vyu;
// 	      BaryonField[ivz  ][igrid] = 0.0;
// 	      BaryonField[ietot][igrid] = etotu;
// 	      if (DualEnergyFormalism) {
// 		BaryonField[ieint][igrid] = pu / ((Gamma-1.0)*rhou);
// 	      }
// 	      BaryonField[iBx  ][igrid] = Bxu;
// 	      BaryonField[iBy  ][igrid] = Byu;
// 	      BaryonField[iBz  ][igrid] = 0.0;
// 	      BaryonField[iPhi ][igrid] = 0.0;
// 	    }

// 	  }




	}}}}
	      
	     

  return SUCCESS;
}

