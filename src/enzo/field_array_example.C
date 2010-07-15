#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h> 
#include <hdf5.h>

#ifdef TRANSFER
#include "preincludes.h"
#endif
#include "svn_version.def"
#include "performance.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#define DEFINE_STORAGE
#include "global_data.h"
#include "units.h"
#include "flowdefs.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "TopGridData.h"
#include "CosmologyParameters.h"
#include "StarParticleData.h"
#undef DEFINE_STORAGE

#ifndef USE_MPI
#define MyProcessorNumber 0
#define ROOT_PROCESSOR 0
#endif

void my_exit(int status);

int SetDefaultGlobalValues(TopGridData &MetaData);
int CommunicationInitialize(Eint32 *argc, char **argv[]);
int CommunicationFinalize();

grid *Linear3DGrid();

Eint32 main(Eint32 argc, char *argv[]) {

  CommunicationInitialize(&argc, &argv);

  grid *agrid = Linear3DGrid();
  
  EnzoArrayFloat *dens = agrid->CreateFieldArrayFloat(Density);
  
  Eint32 index = 7 + 8*134 + 9*134*134;
  
  printf("density rank = %"ISYM"\n", dens->Rank);
  printf("density dim[0]  = %"ISYM"\n", dens->Dimension[0]);
  printf("density start[0]  = %"ISYM"\n", dens->StartIndex[0]);
  printf("density end[0]  = %"ISYM"\n", dens->EndIndex[0], 130);
  printf("density field[7 + 8*134 + 9*134*134] = %"FSYM"\n", dens->Array[index]);
  
  delete dens;    
  delete agrid;
  
  // End the overall test suite
  CommunicationFinalize();

  return 0;
}


grid *Linear3DGrid(){
  // Create a new 3D grid
  float dens = M_PI, total_energy = 0.5, internal_energy = 0.0;
  float vel[3];
  int dims[3];
  FLOAT left[3], right[3];

  grid *lineargrid = new grid;
  int i, j, k, rank = 3;
  int index;
   
  for (i = 0; i < rank; i++) {
    dims[i] = 134;
    left[i] = 0.0;
    right[i] = 1.0;
    vel[i] = (i+1) * 0.125;
  }
    
  NumberOfParticleAttributes = 0;
  lineargrid->PrepareGrid(3, dims, 
			  left, right, 2);

  int result = lineargrid->InitializeUniformGrid(dens, total_energy, internal_energy, vel);
  assert(result != FAIL);
  //  float *dens_field = lineargrid->AccessDensity();

  EnzoArrayFloat *dens_field = lineargrid->CreateFieldArrayFloat("Density");

  for (k = 3; k <= 130; k++) {
    for (j = 3; j <= 130; j++) {
      index =  k*(134)*(134) +
	j*(134) + 3;
      for (i = 3; i <= 130; i++, index++) {
	dens_field->Array[index] = (float)(i + 1000*j + 1000000*k);
      }
    }
  }
  
  delete dens_field;
  
  return lineargrid;

 }

void my_exit(int status)
{
  CommunicationFinalize();
  exit(status);
}
