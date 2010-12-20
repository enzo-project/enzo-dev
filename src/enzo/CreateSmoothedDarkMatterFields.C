/***********************************************************************
/
/  CALCULATE SMOOTHED PARTICLE FIELDS FOR OUTPUT
/
/  written by: John Wise
/  date:       June, 2009
/  modified1:  
/
/  PURPOSE:
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

 
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "CosmologyParameters.h"
#include "communication.h"

#include "FOF_allvars.h"
#include "FOF_nrutil.h"
#include "FOF_proto.h"

void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
void FOF_Initialize(TopGridData *MetaData, 
		    LevelHierarchyEntry *LevelArray[], 
		    FOFData &D, bool SmoothData);
void FOF_Finalize(FOFData &D, LevelHierarchyEntry *LevelArray[], 
		  TopGridData *MetaData, int FOFOnly=FALSE);
int CopyParticlesAcrossPeriodicBoundaries(FOFData &D, int TopGridResolution);
int CommunicationReceiveHandler(fluxes **SubgridFluxesEstimate[] = NULL,
				int NumberOfSubgrids[] = NULL,
				int FluxFlag = FALSE,
				TopGridData* MetaData = NULL);

int CreateSmoothedDarkMatterFields(TopGridData &MetaData, HierarchyEntry *TopGrid)
{

  if (OutputSmoothedDarkMatter == 0)
    return SUCCESS;

  if (OutputSmoothedDarkMatter < 0)
    OutputSmoothedDarkMatter = -OutputSmoothedDarkMatter;

  /* Create a LevelHierarchyEntry array */

  int level;
  LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];
  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    LevelArray[level] = NULL;
  AddLevel(LevelArray, TopGrid, 0);

  /* Initialize halo finder module, which will copy the particles to
     it.  We need this to build the tree and the nearest neighbor
     search. */

  FOFData AllVars;

  // Not used, but will fail because set_units checks if GroupMinLen >
  // DesDensityNgb
  AllVars.GroupMinLen = 50;

  // Set to be 2 cell widths for the interpolation.  This is used in
  // the copying of particles across slabs.
  AllVars.LinkLength = 2;

  set_units(AllVars);
  FOF_Initialize(&MetaData, LevelArray, AllVars, true);
  AllVars.DesDensityNgb = SmoothedDarkMatterNeighbors;

  /* Before we create the tree, we need to (1) copy the particles near
     the slab edges "shadows" to the adjacent slabs and (2) replicate
     particles on the domain boundaries to the other side for
     periodicity. */

  marking(AllVars);

  if (NumberOfProcessors > 1)
    exchange_shadow(AllVars, MetaData.TopGridDims[0], true);

  CopyParticlesAcrossPeriodicBoundaries(AllVars, MetaData.TopGridDims[0]);

  /* Create the tree */

  ngb_treeallocate(AllVars, AllVars.Nlocal, 2*AllVars.Nlocal);
  ngb_treebuild(AllVars, AllVars.Nlocal);
  set_sph_kernel(AllVars);

  /* 
     Interpolate particle quantities to the grids.  

     Communication is required because the tree domain is decomposed
     into slabs, whereas a grids->position correlation doesn't exist.
     First we interpolate to grids that overlap with the slab on this
     processor, then sum those fields on the host processor of each
     grid.
  */

  LevelHierarchyEntry *Temp;

  // Post receive calls

  CommunicationReceiveIndex = 0;
  CommunicationReceiveCurrentDependsOn = COMMUNICATION_NO_DEPENDENCE;
  CommunicationDirection = COMMUNICATION_POST_RECEIVE;

  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel)
      Temp->GridData->InterpolateParticlesToGrid(&AllVars);

  // Interpolate data and send

  CommunicationDirection = COMMUNICATION_SEND;
  if (debug)
    fprintf(stdout, "CreateSmoothedDarkMatterFields: interpolating...\n");

  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel)
      Temp->GridData->InterpolateParticlesToGrid(&AllVars);

  // Wait for the receives and sum fields

  if (CommunicationReceiveHandler() == FAIL)
    ENZO_FAIL("CommunicationReceiveHandler() failed!\n");
  

  // Reset to default
  CommunicationDirection = COMMUNICATION_SEND_RECEIVE;

  // Cleanup
  FOF_Finalize(AllVars, LevelArray, &MetaData);
  ngb_treefree();

  delete [] AllVars.Nslab;
  delete [] AllVars.Nshadow;
  delete [] AllVars.Noffset;
  delete [] AllVars.NtoLeft;
  delete [] AllVars.NtoRight;

  return SUCCESS;

}
