/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR ONE OF FOUR FDM TESTS)
/
/  written by: Xinyu Li
/  date:       May, 2019
/  modified1:
/
/  PURPOSE: set up grid for FDM tests -- type is controlled by parameter
/             LightBosonProblemType (1-4, see below for details)
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

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
#include "phys_constants.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
		      float *TemperatureUnits, float *TimeUnits,
		      float *VelocityUnits, FLOAT Time);

int grid::LightBosonInitializeGrid(float CenterPosition, int LightBosonProblemType)
{  

  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  FieldType[NumberOfBaryonFields++] = Velocity2;
  FieldType[NumberOfBaryonFields++] = Velocity3;
  FieldType[NumberOfBaryonFields++] = TotalEnergy;
  FieldType[NumberOfBaryonFields++] = RePsi;
  FieldType[NumberOfBaryonFields++] = ImPsi;
  FieldType[NumberOfBaryonFields++] = FDMDensity;

  int size = 1, index, dim;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
  int UseParticles = 1;
  FLOAT ParticleMeanDensity = .1;
  static int CollapseTestParticleCount = 0;
  int ParticleCount = 0;

  /* Return if this doesn't concern us. */

  if (ProcessorNumber != MyProcessorNumber) {
    NumberOfParticles = (UseParticles > 0) ? 1 : 0;
    for (dim = 0; dim < GridRank; dim++)
      NumberOfParticles *= (GridEndIndex[dim] - GridStartIndex[dim] + 1);
    return SUCCESS;
  }

  int field;
  for (field = 0; field < NumberOfBaryonFields; field++)
    if (BaryonField[field] == NULL)
      BaryonField[field] = new float[size];

  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
    VelocityUnits = 1, TimeUnits = 1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
         &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.");
  }

  int RePsiNum, ImPsiNum, FDMDensNum;
  RePsiNum = FindField(RePsi, FieldType, NumberOfBaryonFields);
  ImPsiNum = FindField(ImPsi, FieldType, NumberOfBaryonFields);
  FDMDensNum = FindField(FDMDensity, FieldType, NumberOfBaryonFields);

  double coef = (5.9157166856e27*TimeUnits/POW(LengthUnits,2));

  FLOAT x=0,y=0,z=0;
  float xv;

  float alpha = 1./500;
  float initialtime=-0.1;
  float sumsquare = pow(alpha,2)+POW(coef*initialtime,2);

  float a = 0.1;  // offset of Gaussian packets for (4) collision
  float expa, expb;
  float k1 = 2*pi;
  float k2 = 2*pi;

  int i,j,k;
  for (k = 0; k < GridDimension[2]; k++) {
  for (j = 0; j < GridDimension[1]; j++) {
  for (i = 0; i < GridDimension[0]; i++) {

    x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i] - CenterPosition;
    xv = CellLeftEdge[0][i] - CenterPosition;

    index = GRIDINDEX_NOGHOST(i,j,k); 

    // set up an initial Gaussian density profile, v=0

    // 1) Gaussian Density Field

    if (LightBosonProblemType == 1) {
      float r1,i1,phase,r2,i2;
      phase = -POW(x,2)/sumsquare/2*coef*initialtime;
      r1 = cos(phase);
      i1 = -sin(phase);
      r2 = sqrt((sqrt(sumsquare)+alpha)/2.);
      i2 = -sqrt((sqrt(sumsquare)-alpha)/2.);
      BaryonField[iden ][index] = 1*exp(-alpha*pow(x,2)/sumsquare)/pow(sumsquare,0.5)/5; 
      BaryonField[RePsiNum][index] = 1*exp(-alpha*pow(x,2)/sumsquare/2.)/sqrt(sumsquare)*(r1*r2-i1*i2); //sqrt(BaryonField[iden ][index]);
      BaryonField[ImPsiNum][index] = 1*exp(-alpha*pow(x,2)/sumsquare/2.)/sqrt(sumsquare)*(r1*i2+i1*r2);
      BaryonField[FDMDensNum][index] = BaryonField[RePsiNum][index]*BaryonField[RePsiNum][index] + BaryonField[ImPsiNum][index]*BaryonField[ImPsiNum][index];

      //BaryonField[ivx ][index] = 1.e-4*exp(-alpha*pow(x,2)/sumsquare); // Gaussian Density Field
      BaryonField[ivx  ][index] = xv/sumsquare*coef*coef*initialtime;
      BaryonField[ietot][index] = 100;
    }

    // 2) Fresnel solution
      
    if (LightBosonProblemType == 2) {
      if( x>0){
	BaryonField[iden ][index] = 6.;
	BaryonField[RePsiNum][index] = sqrt(6.);
	BaryonField[ImPsiNum][index] = 0;
      }else if (x==0.){
	BaryonField[iden ][index] = 5.5;
	BaryonField[RePsiNum][index] = sqrt(5.5);
	BaryonField[ImPsiNum][index] = 0;
      }else{
	BaryonField[iden ][index] = 5.;
	BaryonField[RePsiNum][index] = sqrt(5.);
	BaryonField[ImPsiNum][index] = 0;
      }
    }

    // 3) Zeldovich Test

    if (LightBosonProblemType == 3) {
      BaryonField[iden  ][index] = 1;
      if ((x<0.5)&&(x>0.4)){
	BaryonField[iden  ][index] = 1+cos(10*3.1415926536*(x-0.45));
      }else if ((x>0.5)&&(x<0.6)){
	BaryonField[iden  ][index] = 1+cos(10*3.1415926536*(x-0.55));
      }
      BaryonField[ivx  ][index] = 0;
      if ((BaryonField[iden ][index]>(1+1e-6))&&(xv<0.5)){
	BaryonField[ivx  ][index] = 1e0;
      }else if ((BaryonField[iden ][index]>(1+1e-6))&&(xv>0.5)){
	BaryonField[ivx  ][index] = -1e0;
      }
    }

    //  4) Two colliding Gaussian packets

    if (LightBosonProblemType == 4) {
      expa = exp(-alpha*pow((x+a),2)/sumsquare/2.)*50;
      expb = exp(-alpha*pow((x-a),2)/sumsquare/2.)*50;
      float rho = expa*expa + expb*expb + expa*expb*2*cos(2*k*x);

      BaryonField[RePsiNum][index] = 25+sin(16*pi*x);//expa*cos(k1*(x+a)) + expb*cos(k2*(x-a)) ;
      BaryonField[ImPsiNum][index] = 0;//expa*sin(k1*(x+a)) - expb*sin(k2*(x-a)) ;
      BaryonField[FDMDensNum][index] = BaryonField[RePsiNum][index]*BaryonField[RePsiNum][index] + BaryonField[ImPsiNum][index]*BaryonField[ImPsiNum][index];
      BaryonField[iden][index] = BaryonField[RePsiNum][index]*BaryonField[RePsiNum][index] + BaryonField[ImPsiNum][index]*BaryonField[ImPsiNum][index];

      BaryonField[ivx][index] = coef*(k1*(expa*expa-expb*expb)-2*a*alpha/sumsquare*expa*expb*sin(2*k*x))/rho;

      BaryonField[ivy  ][index] = 0;
      BaryonField[ivz  ][index] = 0;
      BaryonField[ietot][index] = 100;
    }
  }
  }
  }

  // set particles

  int SetupLoopCount, npart = 0;
  fprintf(stderr, "initialize particles \n" );
  for (SetupLoopCount = 0; SetupLoopCount < 1+min(UseParticles, 1); SetupLoopCount++) {

    /* Set particles. */

    if (UseParticles > 0 && SetupLoopCount > 0) {

      /* If particles already exist (coarse particles), then delete. */

      if (NumberOfParticles > 0) this->DeleteParticles();

      /* Use count from previous loop to set particle number. */

      NumberOfParticles = npart;
      npart = 0;

      /* Allocate space. */

      this->AllocateNewParticles(NumberOfParticles);

      /* Particle values will be set below. */

    } // end: particle initialization
            

    /* Set particles if being used (generate a number of particle
       proportional to density). */

    for (k = 0; k < GridDimension[2]; k++)
    for (j = 0; j < GridDimension[1]; j++)
    for (i = 0; i < GridDimension[0]; i++) {

      /* Compute position */
      index = GRIDINDEX_NOGHOST(i,j,k); 

      x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
      if (GridRank > 1)
	y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
      if (GridRank > 2)
	z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
      if (UseParticles)
	if (i >= GridStartIndex[0] && i <= GridEndIndex[0] &&
	    j >= GridStartIndex[1] && j <= GridEndIndex[1] &&
	    k >= GridStartIndex[2] && k <= GridEndIndex[2]  ) {
	  ParticleCount += int(BaryonField[FDMDensNum][index]/ParticleMeanDensity*0.1);
	  while (ParticleCount >= 1) {
	    if (SetupLoopCount > 0) {
	      ParticleMass[npart] = ParticleMeanDensity;
	      ParticleNumber[npart] = CollapseTestParticleCount++;
	      ParticleType[npart] = PARTICLE_TYPE_DARK_MATTER;
      //fprintf(stderr, "%d %d %d %d %f %f \n", i,j,npart,ParticleCount , BaryonField[FDMDensNum][index],ParticleMeanDensity);

      /* Set random position within cell. */

	      ParticlePosition[0][npart] = x + CellWidth[0][0]*(FLOAT(rand())/FLOAT(RAND_MAX) - 0.5);
	      if (GridRank > 1)
		ParticlePosition[1][npart] = y + CellWidth[1][0]*(FLOAT(rand())/FLOAT(RAND_MAX) - 0.5);
	      if (GridRank > 2)
		ParticlePosition[2][npart] = z + CellWidth[2][0]*(FLOAT(rand())/FLOAT(RAND_MAX) - 0.5);

      /* Set bulk velocity. */
          
	      for (dim = 0; dim < GridRank; dim++){
		ParticleVelocity[dim][npart] = (ParticlePosition[0][npart]-0.5)/sumsquare*coef*coef*initialtime; 
	      }

	    }
	    npart++;
	    ParticleCount -= 1.0;
	  }
	} // end: if statement

    } // end loop over grid
  } // end loop SetupLoopCount

  return SUCCESS;
}

