/***********************************************************************
/
/  GRID CLASS (CONSTRUCTOR)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
//
//  Grid constructor (Set all data to null/default state).
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
 
grid::grid()
{
 
  /* clear scalars */
 
  GridRank                              = 0;
  Time                                  = 0.0;
  OldTime                               = 0.0;
  NumberOfBaryonFields                  = 0;
  dtFixed                               = 0.0;
  NumberOfParticles                     = 0;
  GravitatingMassFieldCellSize          = FLOAT_UNDEFINED;
  GravitatingMassFieldParticlesCellSize = FLOAT_UNDEFINED;
  SubgridsAreStatic                     = FALSE;
  ProcessorNumber                       = ROOT_PROCESSOR;

  SubgridFluxStorage = NULL;
  NumberOfSubgrids = 1;
 
  /* clear MAX_DIMENSION vectors */
 
  int i, j;
  for (i = 0; i < MAX_DIMENSION; i++) {
    GridDimension[i]                 = 1;
    GridStartIndex[i]                = 0;
    GridEndIndex[i]                  = 0;
    GridLeftEdge[i]                  = DomainLeftEdge[i];
    GridRightEdge[i]                 = DomainRightEdge[i];
    CellLeftEdge[i]                  = NULL;
    CellWidth[i]                     = NULL;
    ParticlePosition[i]              = NULL;
    ParticleVelocity[i]              = NULL;
    ParticleAcceleration[i]          = NULL;
    AccelerationField[i]             = NULL;
    GravitatingMassFieldDimension[i] = 0;
    RandomForcingField[i]            = NULL;
    PhaseFctMultEven[i]              = NULL; // WS
    PhaseFctMultOdd[i]               = NULL; // WS
  }
  PhaseFctInitEven = NULL; // WS
  PhaseFctInitOdd  = NULL; // WS
 
  ParticleAcceleration[MAX_DIMENSION]      = NULL;
 
  /* clear MAX_NUMBER_OF_BARYON_FIELDS vectors & [][MAX_DIMENSION] matricies */
 
  for (i = 0; i < MAX_NUMBER_OF_BARYON_FIELDS; i++) {
    BaryonField[i]          = NULL;
    OldBaryonField[i]       = NULL;
    InterpolatedField[i]    = NULL;
    FieldType[i]            = FieldUndefined;
  }

/*
  for (i = 0; i < MAX_NUMBER_OF_BARYON_FIELDS; i++) {
    for (j = 0; j < MAX_DIMENSION; j++ ) {
      BoundaryFluxes->LeftFluxes[i][j]  = NULL;
      BoundaryFluxes->RightFluxes[i][j] = NULL;
    }
  }
*/

#ifdef SAB
  for (i = 0; i < MAX_DIMENSION; i++) {
    OldAccelerationField[i] = NULL;
  }
#endif

  AccelerationHack = FALSE;

  /* Clear miscelaneous pointers */
 
  ParticleMass                  = NULL;
  ParticleNumber                = NULL;
  ParticleType                  = NULL;
  PotentialField                = NULL;
  GravitatingMassField          = NULL;
  GravitatingMassFieldParticles = NULL;
  GravityBoundaryType           = GravityUndefined;
  for (i = 0; i < MAX_NUMBER_OF_PARTICLE_ATTRIBUTES; i++)
    ParticleAttribute[i] = NULL;

  BoundaryFluxes                = NULL;
 
  /* Clear flagging field pointers */
 
  ParticleMassFlaggingField     = NULL;
  MassFlaggingField             = NULL;
  FlaggingField                 = NULL;

  divB = NULL;
  for (int dim=0; dim<3; dim++) {
    gradPhi[dim] = NULL;
  }

#ifdef TRANSFER
  NumberOfPhotonPackages = 0;
  PhotonPackages = new PhotonPackageEntry;
  PhotonPackages->NextPackage     = NULL;
  PhotonPackages->PreviousPackage = NULL;

  FinishedPhotonPackages = new PhotonPackageEntry;
  FinishedPhotonPackages->NextPackage = NULL;
  FinishedPhotonPackages->PreviousPackage = NULL;

  PausedPhotonPackages = new PhotonPackageEntry;
  PausedPhotonPackages->NextPackage = NULL;
  PausedPhotonPackages->PreviousPackage = NULL;
 
  PhotonPackages->Photons         = 1.;
  PhotonPackages->Type            = 0;          
  PhotonPackages->Energy          = 0.;        
  PhotonPackages->EmissionTimeInterval= 0.;      
  PhotonPackages->EmissionTime    = 0.;  
  PhotonPackages->CurrentTime     = 0.;   
  PhotonPackages->Radius          = 0.;        
  PhotonPackages->ipix            = 0;         
  PhotonPackages->level           = 0;        

  sfSeed                          = 0;
  ID                              = 0;
  HasRadiation                    = FALSE;
  SubgridMarker                   = NULL;

  MaximumkphIfront                = 0;
  IndexOfMaximumkph               = INT_UNDEFINED;

  /* Initialize top level parallelism information */
  for (i=0; i<MAX_DIMENSION; i++) {
    ProcLayout[i]       = 1;
    ProcLocation[i]     = 0;
    ProcNeighbors[i][0] = 0;
    ProcNeighbors[i][1] = 0;
  }

  /* Initialize maximum radiation time step size */
  MaxRadiationDt = huge_number;

#endif

  /* Star particles */
  
  NumberOfStars = 0;
  Stars = NULL;

  for(i=0;i<3;i++){
    MagneticField[i] = NULL;
    ElectricField[i] = NULL;
    AvgElectricField[i] = NULL;
    OldMagneticField[i] = NULL;
    OldElectricField[i] = NULL;
    MHDParentTemp[i] = NULL;
  }
  dtParent = -1;

  DyBx = NULL;
  DzBx = NULL;
  DyzBx = NULL;
  DBxFlag = NULL;

  DxBy = NULL;
  DzBy = NULL;
  DxzBy = NULL;
  DByFlag = NULL;

  DxBz = NULL;
  DyBz = NULL;
  DxyBz = NULL;
  DBzFlag = NULL;

  for(int field=0;field<3;field++){
    MagneticSize[field] = -100;
    ElectricSize[field] = -100;
    for(int dim=0;dim<3;dim++){
      MHDAdd[field][dim]=(field==dim) ? 1:0;
      MagneticDims[field][dim] = -100;
    }}


  /* For once-per-rootgrid-timestep star formation, the following flag
     determines whether SF is about to occur or not. It's currently
     (April 2012) only implemented for H2REG_STAR and completely
     ignored for all other star makers. */
  MakeStars = 0;
}
