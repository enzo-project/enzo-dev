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
#ifdef USE_MPI
#include "mpi.h"
#ifdef USE_MPE
#include "mpe.h"
#endif /* USE_MPE */
#endif /* USE_MPI */

#include <hdf5.h>
#include <stdio.h>
#include <string.h>
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
#include "fortran.def"
#include "flowdefs.h"
#include "error.h"
#include "CosmologyParameters.h"

void my_exit(int status);
 
#ifdef PROTO /* Remove troublesome HDF PROTO declaration. */
#undef PROTO
#endif

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

int grid::ParallelFDMCollapseInitializeGrid(char *FDMCollapseRePsiName, 
                                    char *FDMCollapseImPsiName,
                                    char *FDMCollapseAbsBdName,
                                    int FDMUseParticles,
                                    float FDMParticleMeanDensity,
                                    int FDMCollapseSubgridsAreStatic,
                                    int TotalRefinement)
{
  /* declarations */
  int idim, dim, vel, ibx;
  int DeNum;
 
  int ExtraField[2];
 
  inits_type *tempbuffer = NULL;
 
  FILE *log_fptr;
 
#ifdef IO_LOG
  int         io_log = 1;
#else
  int         io_log = 0;
#endif

  char pid[MAX_TASK_TAG_SIZE];
  sprintf(pid, "%"TASK_TAG_FORMAT""ISYM, MyProcessorNumber);
 
  char *logname = new char[MAX_NAME_LENGTH];
  strcpy(logname, "TSlog.");
  strcat(logname,pid);
 
  if (io_log) {
    log_fptr = fopen(logname, "a");
    fprintf(log_fptr, "\n");
    fprintf(log_fptr, "TSIG ParallelRootGridIO = %"ISYM"\n", ParallelRootGridIO);
    fprintf(log_fptr, "Processor %"ISYM", Target processor %"ISYM"\n",
        MyProcessorNumber, ProcessorNumber);
    fprintf(log_fptr, "TotalRefinement = %"ISYM"\n", TotalRefinement);
  }

  /* Determine if the data should be loaded in or not. */
 
  int ReadData = TRUE, Offset[] = {0,0,0};
 
  if (ParallelRootGridIO == TRUE && TotalRefinement == 1)
    ReadData = FALSE;
 
  if (io_log) fprintf(log_fptr, "ReadData = %"ISYM"\n", ReadData);

    /* Calculate buffer Offset (same as Grid unless doing ParallelRootGridIO
     (TotalRefinement = -1 if used as a signal that we should really load
     in the data regardless of the value of ParallelRootGridIO). */
 
  if (ParallelRootGridIO == TRUE && TotalRefinement == -1)
    for (dim = 0; dim < GridRank; dim++)
      Offset[dim] = nint((GridLeftEdge[dim] - DomainLeftEdge[dim])/CellWidth[dim][0]);

  int i, j, k, m, field, sphere, size, iden;
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

  FieldType[NumberOfBaryonFields++] = RePsi;
  FieldType[NumberOfBaryonFields++] = ImPsi;
  FieldType[NumberOfBaryonFields++] = FDMDensity;
  
  //printf("%d \n", NumberOfBaryonFields);
  if(WritePotential)
    FieldType[NumberOfBaryonFields++] = GravPotential;

  /* Set the subgrid static flag. */
 
  SubgridsAreStatic = FDMCollapseSubgridsAreStatic;

  if (ProcessorNumber == MyProcessorNumber) {
 
  /* Skip following if NumberOfBaryonFields == 0. */
 
  if (NumberOfBaryonFields > 0) {

      int DensNum = -1, GENum = -1, Vel1Num = -1, 
                         Vel2Num=-1, Vel3Num=-1, TENum=-1,
                         B1Num=-1, B2Num=-1, B3Num=-1, PhiNum=-1;
      IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, 
                           Vel2Num, Vel3Num, TENum,
                           B1Num, B2Num, B3Num, PhiNum);
      
      int RePsiNum = -1, ImPsiNum=-1, FDMDensNum=-1;
      
      RePsiNum = FindField(RePsi, FieldType, NumberOfBaryonFields);
      ImPsiNum = FindField(ImPsi, FieldType, NumberOfBaryonFields);
      FDMDensNum = FindField(FDMDensity, FieldType, NumberOfBaryonFields);

    /* Determine the size of the fields. */
 
    int size = 1;
 
    for (dim = 0; dim < GridRank; dim++)
      size *= GridDimension[dim];
 
    /* Allocate space for the fields. */
 
    if (ReadData == TRUE) {
      this->AllocateGrids();
    }
    
    if (FDMCollapseAbsBdName !=NULL && ReadData){
      if (ReadFile(FDMCollapseAbsBdName, GridRank, GridDimension,
              GridStartIndex, GridEndIndex, Offset, BaryonField[0],
              &tempbuffer, 0, 1) == FAIL) {
        ENZO_FAIL("Error reading density field.\n");
      }
    }

    if (FDMCollapseRePsiName != NULL && ReadData){
      if (ReadFile(FDMCollapseRePsiName, GridRank, GridDimension,
              GridStartIndex, GridEndIndex, Offset, BaryonField[RePsiNum],
              &tempbuffer, 0, 1) == FAIL) {
        ENZO_FAIL("Error reading real part of wave function.\n");
      }
    }
  
    if (FDMCollapseImPsiName != NULL && ReadData){
      if (ReadFile(FDMCollapseImPsiName, GridRank, GridDimension,
              GridStartIndex, GridEndIndex, Offset, BaryonField[ImPsiNum],
              &tempbuffer, 0, 1) == FAIL) {
        ENZO_FAIL("Error reading imaginary part of wave function.\n");
      }
    }

  } // end: if (NumberOfBaryonFields > 0)
  } // end: if (ProcessorNumber == MyProcessorNumber)
  OldTime = Time;
 
  if (io_log) fclose(log_fptr);

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

  double afloat = double(a);
  double hmcoef = 5.9157166856e27*TimeUnits/POW(LengthUnits/afloat,2)/FDMMass;

  // If use particle, initial particles according to the FDM values and turn off QuantumPressure
  int CollapseTestParticleCount = 0;
  int SetupLoopCount, npart = 0;
  int ParticleCount = 0;
  int ind, indxp, indxn, indyp, indyn, indzp, indzn;
  int ip,in,jp,jn,kp,kn;
  double x,y,z,vx,vy,vz;
  double r,theta,phi;
  double cluster_radius = 200*3e18/LengthUnits;

  if (FDMUseParticles > 0){
    if (ProcessorNumber != MyProcessorNumber) {
      NumberOfParticles = (FDMUseParticles > 0) ? 1 : 0;
    for (dim = 0; dim < GridRank; dim++)
      NumberOfParticles *= (GridEndIndex[dim] - GridStartIndex[dim] + 1);
    return SUCCESS;
    }

    fprintf(stdout, "FDMCollapse: initialize particles on processor %d \n", MyProcessorNumber);
    for (SetupLoopCount = 0; SetupLoopCount < 1+min(FDMUseParticles, 1); SetupLoopCount++) {
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
	ParticleCount = 1000;

	// set many particles
	while (ParticleCount > 0) {
        if (SetupLoopCount > 0) {
		    ParticleMass[npart] = FDMParticleMeanDensity;
	        ParticleNumber[npart] = CollapseTestParticleCount++;
            ParticleType[npart] = PARTICLE_TYPE_DARK_MATTER;
         // Set random position within cell.
		    theta = acos(2*(FLOAT(rand())/FLOAT(RAND_MAX) - 0.5));
			phi = 2*M_PI*(FLOAT(rand())/FLOAT(RAND_MAX));
		    r = cluster_radius*POW(FLOAT(rand())/FLOAT(RAND_MAX),1./3.);
		    ParticlePosition[0][npart] = 0.5 + r*sin(theta)*cos(phi);
		    ParticlePosition[1][npart] = 0.5 + r*sin(theta)*sin(phi);
		    ParticlePosition[2][npart] = 0.5 + r*cos(theta);
			ParticleVelocity[0][npart] = 0.0;//-2.5e6/(LengthUnits/TimeUnits)*sin(theta);
			ParticleVelocity[1][npart] = 0.0;//2.0e6 /(LengthUnits/TimeUnits);//*cos(theta);
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
   printf("FDMCollapseInitialize: Number of Particles = %d on Processor %d\n", NumberOfParticles, MyProcessorNumber);


  // turn off quantum pressure, do a pure CDM sim
  // QuantumPressure = 0;
  }
  return SUCCESS;
}
