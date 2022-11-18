/***********************************************************************
/
/  Construct a fake grid for feedback algorithms based on a
/  list of active particles
/
/  written by: Nathan Goldbaum
/  date:       June 2012
/
************************************************************************/

#ifdef USE_MPI
#endif

#include "preincludes.h"

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "phys_constants.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "ActiveParticle.h"
#include "phys_constants.h"
#include "CommunicationUtilities.h"
#include "communication.h"

#include <stdio.h>
#include <stdlib.h>

int CommunicationBufferPurge(void);
int CommunicationReceiveHandler(fluxes **SubgridFluxesEstimate[] = NULL,
				int NumberOfSubgrids[] = NULL,
				int FluxFlag = FALSE,
				TopGridData* MetaData = NULL);



grid* ConstructFeedbackZone(ActiveParticleType* ThisParticle,int FeedbackRadius,
			    FLOAT dx, HierarchyEntry** Grids, int NumberOfGrids,
			    int SendField)
{
  int i,j,dim,size;
  int FeedbackZoneRank;
  FLOAT FBRdx;
  FLOAT ParticlePosition[3] = {ThisParticle->ReturnPosition()[0],
			       ThisParticle->ReturnPosition()[1],
			       ThisParticle->ReturnPosition()[2]};

  /* Build array of AP grids and check for errors */
  grid* APGrid;

  APGrid = ThisParticle->ReturnCurrentGrid();
  FBRdx = dx * FLOAT(FeedbackRadius);

  if (APGrid == NULL)
    ENZO_FAIL("Particle CurrentGrid is invalid!\n");


  std::FILE * particle_position_file = std::fopen("ActiveParticlePosition.out", "a");
  std::fprintf(particle_position_file, "%.7e, %.7e, %.7e, %.7e, %.7e, %.7e\n", 
    ParticlePosition[0]-FBRdx, 
    ParticlePosition[1]-FBRdx, 
    ParticlePosition[2]-FBRdx,
    ParticlePosition[0]+FBRdx, 
    ParticlePosition[1]+FBRdx, 
    ParticlePosition[2]+FBRdx);
  std::fclose(particle_position_file);

  std::FILE * grid_position_file = std::fopen("GridPosition.out", "a");
  std::fprintf(grid_position_file, "%.7e, %.7e, %.7e, %.7e, %.7e, %.7e\n", 
    APGrid->GetGridLeftEdge(0), 
    APGrid->GetGridLeftEdge(1), 
    APGrid->GetGridLeftEdge(2),
    APGrid->GetGridRightEdge(0), 
    APGrid->GetGridRightEdge(1), 
    APGrid->GetGridRightEdge(2));
  std::fclose(grid_position_file);


  // This should only happen if the grid pointer is invalid
  if ((APGrid->GetGridLeftEdge(0) > ParticlePosition[0]+FBRdx) ||
      (APGrid->GetGridLeftEdge(1) > ParticlePosition[1]+FBRdx) ||
      (APGrid->GetGridLeftEdge(2) > ParticlePosition[2]+FBRdx) ||
      (APGrid->GetGridRightEdge(0) < ParticlePosition[0]-FBRdx) ||
      (APGrid->GetGridRightEdge(1) < ParticlePosition[1]-FBRdx) ||
      (APGrid->GetGridRightEdge(2) < ParticlePosition[2]-FBRdx)){

        std::FILE * err_report = std::fopen("ErrorReport.out", "a");
        std::fprintf(err_report, "Truth Values:\n");
        std::fprintf(err_report, "%s, %s, %s, %s, %s, %s\n\n", 
          APGrid->GetGridLeftEdge(0) > ParticlePosition[0]+FBRdx ? "true" : "false", 
          APGrid->GetGridLeftEdge(1) > ParticlePosition[1]+FBRdx ? "true" : "false", 
          APGrid->GetGridLeftEdge(2) > ParticlePosition[2]+FBRdx ? "true" : "false",
          APGrid->GetGridRightEdge(0) < ParticlePosition[0]-FBRdx ? "true" : "false", 
          APGrid->GetGridRightEdge(1) < ParticlePosition[1]-FBRdx ? "true" : "false", 
          APGrid->GetGridRightEdge(2) < ParticlePosition[2]-FBRdx ? "true" : "false");

        std::fprintf(err_report, "FBRdx value:%.7e\n", FBRdx);
        std::fprintf(err_report, "FeedbackRadius value:%d\n", FeedbackRadius);
        std::fprintf(err_report, "Cell size:%.7e\n", dx);
        std::fprintf(err_report, "Grid->NumberOfActiveParticles :%d\n", APGrid->ReturnNumberOfActiveParticles());

        std::fclose(err_report);

        ENZO_FAIL("Particle outside own grid!\n");
      }
    


  /* Setup Feedback Zones before copying data */

  FeedbackZoneRank = APGrid->GetGridRank();

  int FeedbackZoneDimension[MAX_DIMENSION];
  FLOAT LeftCellOffset[MAX_DIMENSION],FeedbackZoneLeftEdge[MAX_DIMENSION],
    FeedbackZoneRightEdge[MAX_DIMENSION], ncells[MAX_DIMENSION];
  FLOAT CellSize, GridGZLeftEdge;

  for (dim = 0; dim < FeedbackZoneRank; dim++) {
    FeedbackZoneDimension[dim] = (2*(FeedbackRadius+NumberOfGhostZones)+1);
    CellSize = APGrid->GetCellWidth(dim,0);
    GridGZLeftEdge = APGrid->GetCellLeftEdge(dim,0);

    LeftCellOffset[dim] = MODF((ParticlePosition[dim]-GridGZLeftEdge)/CellSize,&ncells[dim]);

    FeedbackZoneLeftEdge[dim]  = GridGZLeftEdge + CellSize*(ncells[dim]-FeedbackRadius);
    FeedbackZoneRightEdge[dim] = GridGZLeftEdge + CellSize*(ncells[dim]+FeedbackRadius+1);
  }

  grid* FeedbackZone = new grid;

  FeedbackZone->InheritProperties(APGrid);

  FeedbackZone->PrepareGrid(FeedbackZoneRank, FeedbackZoneDimension,
			    FeedbackZoneLeftEdge,FeedbackZoneRightEdge,0);

  FeedbackZone->SetProcessorNumber(APGrid->ReturnProcessorNumber());

  FeedbackZone->SetTimeStep(APGrid->ReturnTimeStep());

  // This will only allocate the BaryonField on the host processor
  if (FeedbackZone->AllocateAndZeroBaryonField() == FAIL)
    ENZO_FAIL("FeedbackZone BaryonField allocation failed\n");

  if (SendField == GRAVITATING_MASS_FIELD) {
    size = 1;
    // If we're doing gravity, we need to init that field.
    FeedbackZone->InitializeGravitatingMassField(RefineBy);
    for (dim = 0; dim < FeedbackZoneRank; dim++) {
      size *= FeedbackZone->ReturnGravitatingMassFieldDimension(dim);
    }
    FeedbackZone->InitGravitatingMassField(size);
  }

  // Copy zones from this grid (which must overlap the position of the AP).
  FLOAT ZeroVector[] = {0,0,0};

  /* Post receives */

  CommunicationReceiveIndex = 0;
  CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
  CommunicationDirection = COMMUNICATION_POST_RECEIVE;

  for (j = 0; j < NumberOfGrids; j++)
    if (FeedbackZone->CopyActiveZonesFromGrid(Grids[j]->GridData,ZeroVector,SendField) == FAIL)
      ENZO_FAIL("FeedbackZone copy failed!\n");

  /* Send data */

  CommunicationDirection = COMMUNICATION_SEND;

  for (j = 0; j < NumberOfGrids; j++)
    if (FeedbackZone->CopyActiveZonesFromGrid(Grids[j]->GridData,ZeroVector,SendField) == FAIL)
      ENZO_FAIL("FeedbackZone copy failed!\n");

  /* Receive data */

  if (CommunicationReceiveHandler() == FAIL)
    ENZO_FAIL("CommunicationReceiveHandler() failed!\n");

#ifdef USE_MPI
  CommunicationBufferPurge();
#endif

  return FeedbackZone;

}
