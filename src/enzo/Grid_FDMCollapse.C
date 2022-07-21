/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A COLLAPSE TEST)
/
/  written by: Xinyu Li
/  date:       May, 2019
/  modified1:
/
/  PURPOSE: set up problem to collapse a FDM halo
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
#include "phys_constants.h"
#include "CosmologyParameters.h"

#define NTHETA 1000
#define NR 1000

/********************* PROTOTYPES *********************/
int ReadFile(char *name, int Rank, int Dims[], int StartIndex[],
       int EndIndex[], int BufferOffset[], float *buffer,
       inits_type **tempbuffer, int Part, int Npart);
 
int ReadIntFile(char *name, int Rank, int Dims[], int StartIndex[],
    int EndIndex[], int BufferOffset[], int *buffer,
    int **tempbuffer, int Part, int Npart);

void ReadAttribute(hid_t dset_id, int *Attribute, char *AttributeName, FILE *log_fptr, int io_log);
int ReadAttr(char *Fname, int *Rank, int Dims[], int *NSeg, int *LSeg, FILE *log_fptr);

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

/*******************************************************/
#define READFILE ReadFile

int grid::FDMCollapseInitializeGrid(int UseParticles, float ParticleMeanDensity)
{
  /* declarations */

  int dim, i, j, k, m, field, sphere, size, iden;
  int RePsiNum, ImPsiNum, FDMDensNum;
  float xdist,ydist,zdist;

  /* create fields */

  NumberOfBaryonFields = 0;
  FieldType[NumberOfBaryonFields++] = Density;
  FieldType[NumberOfBaryonFields++] = TotalEnergy;

  int ivel = NumberOfBaryonFields;
  FieldType[NumberOfBaryonFields++] = Velocity1;
  if (GridRank > 1) 
    FieldType[NumberOfBaryonFields++] = Velocity2;
  if (GridRank > 2)
    FieldType[NumberOfBaryonFields++] = Velocity3;

  FieldType[RePsiNum = NumberOfBaryonFields++] = RePsi;
  FieldType[ImPsiNum = NumberOfBaryonFields++] = ImPsi;
  FieldType[FDMDensNum = NumberOfBaryonFields++] = FDMDensity;
  
  //printf("%d \n", NumberOfBaryonFields);
  if( WritePotential  )
    FieldType[NumberOfBaryonFields++] = GravPotential;

  /* Set various units. */
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, 
    VelocityUnits, CriticalDensity = 1, BoxLength = 1, mu = 0.6;

  FLOAT a, dadt, ExpansionFactor = 1;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, 
	   &VelocityUnits, Time);
  if (ComovingCoordinates) {
    CosmologyComputeExpansionFactor(Time, &a, &dadt);
    ExpansionFactor = a/(1.0+InitialRedshift);
    CriticalDensity = 2.78e11*pow(HubbleConstantNow, 2); // in Msolar/Mpc^3
    BoxLength = ComovingBoxSize*ExpansionFactor/HubbleConstantNow;  // in Mpc
  } else {
    CriticalDensity = 2.78e11*pow(0.74,2); // in Msolar/Mpc^3 for h=0.74
    BoxLength = LengthUnits / 3.086e24;
    HubbleConstantNow = 1.0;
    OmegaMatterNow = 1.0;
	a = 1.0;
  }

// Determine the size of the fields
 
  size = 1;
 
  for (dim = 0; dim < GridRank; dim++)
        size *= GridDimension[dim];
  int ReadData = TRUE, Offset[] = {0,0,0};
  inits_type *tempbuffer = NULL;
  
// Allocate Fields
  for (int field = 0; field < NumberOfBaryonFields; field++)
    BaryonField[field] = new float[size];

 double afloat = double(a);
 double hmcoef = 5.9157166856e27*TimeUnits/POW(LengthUnits/afloat,2)/FDMMass;

if(FDMCollapseAbsorbingBoundary){
// Read Density, use it as the absorption coefficient on the boundary
  if (READFILE("GridDensity", GridRank, GridDimension,
         GridStartIndex, GridEndIndex, Offset, BaryonField[0],
         &tempbuffer, 0, 1) == FAIL) {
    ENZO_FAIL("Error reading density.\n");}
}

// Read wavefunction
  if (READFILE("GridRePsi", GridRank, GridDimension,
         GridStartIndex, GridEndIndex, Offset, BaryonField[RePsiNum],
         &tempbuffer, 0, 1) == FAIL) {
    ENZO_FAIL("Error reading real part of wave function.\n");}
  
  if (READFILE("GridImPsi", GridRank, GridDimension,
         GridStartIndex, GridEndIndex, Offset, BaryonField[ImPsiNum],
         &tempbuffer, 0, 1) == FAIL) {
    ENZO_FAIL("Error reading imaginary part of wave function.\n");
    }    
  for (i=0; i<size; i++){
    BaryonField[FDMDensNum][i] = BaryonField[RePsiNum][i] * BaryonField[RePsiNum][i] + BaryonField[ImPsiNum][i] * BaryonField[ImPsiNum][i];
  }
 // }

  // If use particle, initial particles according to the FDM values and turn off QuantumPressure
  int CollapseTestParticleCount = 0;
  int SetupLoopCount, npart = 0;
  int ParticleCount = 0;
  int ind, indxp, indxn, indyp, indyn, indzp, indzn;
  int ip,in,jp,jn,kp,kn;
  double x,y,z,vx,vy,vz;
  double cluster_size = kpc_cm/LengthUnits;

  if (UseParticles > 0){
    if (ProcessorNumber != MyProcessorNumber) {
      NumberOfParticles = (UseParticles > 0) ? 1 : 0;
    for (dim = 0; dim < GridRank; dim++)
      NumberOfParticles *= (GridEndIndex[dim] - GridStartIndex[dim] + 1);
    return SUCCESS;
    }


      fprintf(stderr, "initialize particles \n" );
    for (SetupLoopCount = 0; SetupLoopCount < 1+min(UseParticles, 1); SetupLoopCount++) {
     if (SetupLoopCount > 0) {
      /* If particles already exist (coarse particles), then delete. */
        if (NumberOfParticles > 0) this->DeleteParticles();
      /* Use count from previous loop to set particle number. */
        NumberOfParticles = npart;
        npart = 0;
      /* Allocate space. */
        this->AllocateNewParticles(NumberOfParticles);
      /* Particle values will be set below. */
      }
	// set some test particles
	ParticleCount = 10000;

	// set many particles
	while (ParticleCount > 0) {
        if (SetupLoopCount > 0) {
		    ParticleMass[npart] = ParticleMeanDensity;
	        ParticleNumber[npart] = CollapseTestParticleCount++;
            ParticleType[npart] = PARTICLE_TYPE_DARK_MATTER;
         // Set random position within cell.
		    double theta = 3.1415927/6./1000*npart;
		    ParticlePosition[0][npart] = 0.5 + cluster_size*(FLOAT(rand())/FLOAT(RAND_MAX) - 0.5);
		    ParticlePosition[1][npart] = 0.5 + cluster_size*(FLOAT(rand())/FLOAT(RAND_MAX) - 0.5);
		    ParticlePosition[2][npart] = 0.5 + cluster_size*(FLOAT(rand())/FLOAT(RAND_MAX) - 0.5);
			ParticleVelocity[0][npart] = 0.0;//-2.5e6/(LengthUnits/TimeUnits)*sin(theta);
			ParticleVelocity[1][npart] = 0.0;// hmcoef*8*3.1415927;//2.0e6 /(LengthUnits/TimeUnits);//*cos(theta);
			ParticleVelocity[2][npart] = 0.0;
		}
		ParticleCount--;
		npart++;
	}

    // set particles following the FDM density
    /*
    for (k = 0; k < GridDimension[2]; k++)
    for (j = 0; j < GridDimension[1]; j++)
    for (i = 0; i < GridDimension[0]; i++) {
	  // Compute position 
    	  x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	    if (GridRank > 1)
	      y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	    if (GridRank > 2)
	      z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
      
      if (i >= GridStartIndex[0] && i <= GridEndIndex[0] &&
		  j >= GridStartIndex[1] && j <= GridEndIndex[1] &&
		  k >= GridStartIndex[2] && k <= GridEndIndex[2]  ) {
	      ind = GRIDINDEX_NOGHOST(i,j,k);
          indxp = GRIDINDEX_NOGHOST(i+1,j,k);
          indxn = GRIDINDEX_NOGHOST(i-1,j,k); 
          if (GridRank>1){
            indyp = GRIDINDEX_NOGHOST(i,j+1,k);
            indyn = GRIDINDEX_NOGHOST(i,j-1,k); 
          }
          if (GridRank>2){
            indzp = GRIDINDEX_NOGHOST(i,j,k+1);
            indzn = GRIDINDEX_NOGHOST(i,j,k-1);
          }
		  //printf("%d %d %d %d \n",size, i,j,k);
          //printf("x,y,z %f %f %f \n",x,y,z);
          //printf("%d %d %d \n",ind,indxp,indxn);

		  ParticleCount += int(BaryonField[FDMDensNum][ind]/ParticleMeanDensity);
	      
		  while (ParticleCount > 1) {
		      if (SetupLoopCount > 0) {
		        ParticleMass[npart] = ParticleMeanDensity;
		        ParticleNumber[npart] = CollapseTestParticleCount++;
		        ParticleType[npart] = PARTICLE_TYPE_DARK_MATTER;

		  // Set random position within cell.
      		    ParticlePosition[0][npart] = x + CellWidth[0][0]*(FLOAT(rand())/FLOAT(RAND_MAX) - 0.5);
		        ParticlePosition[1][npart] = y + CellWidth[1][0]*(FLOAT(rand())/FLOAT(RAND_MAX) - 0.5);
		        ParticlePosition[2][npart] = z + CellWidth[2][0]*(FLOAT(rand())/FLOAT(RAND_MAX) - 0.5);

		  // Set bulk velocity.
          // vx
            vx = (BaryonField[RePsiNum][ind]*(BaryonField[ImPsiNum][indxp]-BaryonField[ImPsiNum][indxn])
                  - BaryonField[ImPsiNum][ind]*(BaryonField[RePsiNum][indxp]-BaryonField[RePsiNum][indxn]))
                  *hmcoef/BaryonField[FDMDensNum][ind]/(2*CellWidth[0][i]);
            //vx = max(vx,10);
            ParticleVelocity[0][npart] = vx;
            //printf("vx %f \n",vx);
            if (GridRank>1){
              vy = (BaryonField[RePsiNum][ind]*(BaryonField[ImPsiNum][indyp]-BaryonField[ImPsiNum][indyn])
                  - BaryonField[ImPsiNum][ind]*(BaryonField[RePsiNum][indyp]-BaryonField[RePsiNum][indyn]))
                  *hmcoef/BaryonField[FDMDensNum][ind]/(2*CellWidth[1][j]);
              //vy = max(vy,10);
              ParticleVelocity[1][npart] = vy;
              //printf("vy %f \n",vy);
            }
            if (GridRank>2){
              vz = (BaryonField[RePsiNum][ind]*(BaryonField[ImPsiNum][indzp]-BaryonField[ImPsiNum][indzn])
                  - BaryonField[ImPsiNum][ind]*(BaryonField[RePsiNum][indzp]-BaryonField[RePsiNum][indzn]))
                  *hmcoef/BaryonField[FDMDensNum][ind]/(2*CellWidth[2][k]);
              //vz = max(vz,10);
              ParticleVelocity[2][npart] = vz;
              //printf("vz %f \n",vz);
            }
           } 
          npart++;
	      ParticleCount -= 1.0;
        }// end while
      }
      }// end for loop over grid */ 
   } // end loop SetupLoopCount
   NumberOfParticles = npart;
   printf("Number of Particles %d \n", NumberOfParticles);

  // turn off quantum pressure, do a pure CDM sim
  // QuantumPressure = 0;
  }
  return SUCCESS;
}
