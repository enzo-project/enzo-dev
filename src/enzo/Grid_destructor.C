/***********************************************************************
/
/  GRID CLASS (DESTRUCTOR)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
//
//  Grid destructor
//
#include <stdio.h>
#include <stdlib.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
 
/* function prototypes */
 
void DeleteFluxes(fluxes *Fluxes);
void WriteListOfInts(FILE *fptr, int N, int nums[]);
void DeleteStarList(Star * &Node);
#ifdef TRANSFER
PhotonPackageEntry*  DeletePhotonPackage(PhotonPackageEntry *PP);
#endif /* TRANSFER */
 
grid::~grid()
{
 
  int i;
 
  /* Error check. */
 
#ifdef UNUSED
  if (NumberOfParticles > 0) {
    fprintf(stderr, "warning: destroying live particles (%"ISYM").\n",
	    NumberOfParticles);
  /* exit(EXIT_FAILURE); */
  }
#endif /* UNUSED */
 
  for (i = 0; i < MAX_DIMENSION; i++) {
    delete [] CellLeftEdge[i];
    delete [] CellWidth[i];
    delete [] ParticlePosition[i];
    delete [] ParticleVelocity[i];
    delete [] ParticleAcceleration[i];
    delete [] AccelerationField[i];
    delete [] RandomForcingField[i];
    if (PhaseFctMultEven[i] != NULL) delete[] PhaseFctMultEven[i];
    if (PhaseFctMultOdd[i] != NULL) delete[] PhaseFctMultOdd[i];
  }
 
  if (PhaseFctInitEven != NULL) delete[] PhaseFctInitEven;
  if (PhaseFctInitOdd != NULL) delete[] PhaseFctInitOdd;

  delete ParticleAcceleration[MAX_DIMENSION];
 
  for (i = 0; i < MAX_NUMBER_OF_BARYON_FIELDS; i++) {
    delete [] BaryonField[i];
    delete [] OldBaryonField[i];
    delete [] InterpolatedField[i];
  }

#ifdef SAB
  for (i = 0; i < MAX_DIMENSION; i++) {
    if(OldAccelerationField[i] != NULL ){
      delete [] OldAccelerationField[i];
      OldAccelerationField[i] = NULL;
    }
  }
#endif
 
  DeleteFluxes(BoundaryFluxes);
  delete BoundaryFluxes;
 
  delete [] ParticleMass;
  delete [] ParticleNumber;
  delete [] ParticleType;
  delete [] PotentialField;
  delete [] GravitatingMassField;
  delete [] GravitatingMassFieldParticles;
  delete [] FlaggingField;
  delete [] MassFlaggingField;
  delete [] ParticleMassFlaggingField;
 
  for (i = 0; i < MAX_NUMBER_OF_PARTICLE_ATTRIBUTES; i++)
    delete [] ParticleAttribute[i];

  delete divB;
  for (int i=0; i<3; i++) {
    delete gradPhi[i];
  }


  DeleteStarList(Stars);

#ifdef TRANSFER
  delete PhotonPackages;
  if (FinishedPhotonPackages != NULL)
    delete FinishedPhotonPackages;
  if (PausedPhotonPackages != NULL)
    delete PausedPhotonPackages;
  delete [] SubgridMarker;
#endif

/* 
  if (debug && GridRank > 0) {
    printf("grid->destructor: deleting grid with dims = ");
    WriteListOfInts(stdout, GridRank, GridDimension);
  }
*/

  //MHD stuff 
 
  if( UseMHDCT ){
    for(i=0;i<3;i++){

      if(MagneticField[i] != NULL ){
	delete MagneticField[i];
	MagneticField[i] = NULL;
      }
      if(OldMagneticField[i] != NULL ){
	delete OldMagneticField[i];
	OldMagneticField[i] = NULL;
      }


      if(ElectricField[i] != NULL){
	delete ElectricField[i];
	ElectricField[i] = NULL;
      }

      if(AvgElectricField[i] != NULL){
	delete AvgElectricField[i];
	AvgElectricField[i] = NULL;
      }

    }

  }

}
