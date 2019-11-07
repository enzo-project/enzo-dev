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


int grid::FDMCollapseInitializeGrid()
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

  if (QuantumPressure){
    FieldType[RePsiNum = NumberOfBaryonFields++] = RePsi;
    FieldType[ImPsiNum = NumberOfBaryonFields++] = ImPsi;
    FieldType[FDMDensNum = NumberOfBaryonFields++] = FDMDensity;
  }
  //printf("%d \n", NumberOfBaryonFields);

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

if(QuantumPressure){
// Read Density
  if (READFILE("GridDensity.new", GridRank, GridDimension,
         GridStartIndex, GridEndIndex, Offset, BaryonField[0],
         &tempbuffer, 0, 1) == FAIL) {
    ENZO_FAIL("Error reading density.\n");}
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

  }


  return SUCCESS;
}
