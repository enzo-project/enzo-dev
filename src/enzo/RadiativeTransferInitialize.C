/***********************************************************************
/
/  INITIALIZE THE RADIATIVE TRANSFER MODULE
/
/  written by: John Wise
/  date:       March, 2009
/  modified1:
/
/ PURPOSE: This routine initializes everything that is needed by the
/          ray tracer that is not included in the normal routines.
/
************************************************************************/
#include "preincludes.h"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
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
#include "LevelHierarchy.h"
#include "StarParticleData.h"
#include "RadiativeTransferHealpixRoutines.h"
#include "ImplicitProblemABC.h"
#include "FSProb.h"
#include "NullProblem.h"

int RadiativeTransferReadParameters(FILE *fptr);
int ReadPhotonSources(FILE *fptr, FLOAT CurrentTime);
int DetermineParallelism(HierarchyEntry *TopGrid, TopGridData &MetaData);
void my_exit(int status);

int RadiativeTransferInitialize(char *ParameterFile, 
				HierarchyEntry &TopGrid, 
				TopGridData &MetaData,
				ExternalBoundary &Exterior, 
				ImplicitProblemABC* &ImplicitSolver,
				LevelHierarchyEntry *LevelArray[])
{

  const char	*kphHIName     = "HI_kph";
  const char	*gammaHIName   = "HI_gamma";
  const char	*kphHeIName    = "HeI_kph";
  const char	*gammaHeIName  = "HeI_gamma";
  const char	*kphHeIIName   = "HeII_kph";
  const char	*gammaHeIIName = "HeII_gamma";
  const char	*kdissH2IName  = "H2I_kdiss";
  const char	*RadAccel1Name = "RadAccel1";
  const char	*RadAccel2Name = "RadAccel2";
  const char	*RadAccel3Name = "RadAccel3";
  const char	*MetalName     = "Metal_Density";
  const char	*ColourName    = "SN_Colour";

  int i, j, k, level;
  FILE *fptr;
  LevelHierarchyEntry *Temp;

  /* Read and set parameter values and static radiation sources */

  if ((fptr = fopen(ParameterFile, "r")) == NULL) {
    fprintf(stderr, "Error opening ParameterFile %s\n", ParameterFile);
    ENZO_FAIL("");
  }

  RadiativeTransferReadParameters(fptr);
  rewind(fptr);
  if (ProblemType == 50)
    ReadPhotonSources(fptr, MetaData.Time);

  PhotonTime = MetaData.Time;
  MetaData.FLDTime = MetaData.Time;
  MetaData.dtFLD = 0.0;

  fclose(fptr);

  if (RadiativeTransferPhotonEscapeRadius > 0) {
    PhotonEscapeFilename = new char[80];
    sprintf(PhotonEscapeFilename, "fesc%4.4d.dat", 
	    (Eint32) MetaData.DataDumpNumber-1);
  }

  /* Create all StarParticles from normal particles */

//  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
//    for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel)
//      if (Temp->GridData->FindNewStarParticles(level) == FAIL) {
//	fprintf(stderr, "Error in grid::FindStarParticles.\n");
//	ENZO_FAIL("");
//      }

  /* If we're restarting from a non-radiative run, create radiation
     fields.  We will check if the field already exists inside the
     grid routine. */

  int OldNumberOfBaryonFields = 0, FieldsToAdd = 0;
  int TypesToAdd[MAX_NUMBER_OF_BARYON_FIELDS];
  int ExistingTypes[MAX_NUMBER_OF_BARYON_FIELDS];

  for (i = 0; i < MAX_NUMBER_OF_BARYON_FIELDS; i++)
    ExistingTypes[i] = FieldUndefined;

  if (RadiativeTransfer) {
    for (i = kphHI; i <= gammaHeII; i++)
      TypesToAdd[FieldsToAdd++] = i;
    if (MultiSpecies > 1)
      TypesToAdd[FieldsToAdd++] = kdissH2I;
    if (RadiationPressure)
      for (i = RadPressure0; i <= RadPressure2; i++)
	TypesToAdd[FieldsToAdd++] = i;
    if (PopIIISupernovaUseColour)
      TypesToAdd[FieldsToAdd++] = SNColour;
    if (StarClusterUseMetalField)
      TypesToAdd[FieldsToAdd++] = Metallicity;

    for (i = FieldsToAdd; i < MAX_NUMBER_OF_BARYON_FIELDS; i++)
      TypesToAdd[i] = FieldUndefined;

    /* Check if the fields already exist. */
    
    OldNumberOfBaryonFields = LevelArray[0]->GridData->
      ReturnNumberOfBaryonFields();
    LevelArray[0]->GridData->ReturnFieldType(ExistingTypes);
    for (i = 0; i < FieldsToAdd; i++)
      for (j = 0; j < OldNumberOfBaryonFields; j++)
	if (TypesToAdd[i] == ExistingTypes[j]) {
	  for (k = i; k < FieldsToAdd; k++)
	    TypesToAdd[k] = TypesToAdd[k+1];
	  i--;
	  break;
	} // ENDIF matching type
    FieldsToAdd = 0;
    while (TypesToAdd[FieldsToAdd] != FieldUndefined)
      FieldsToAdd++;

  } // ENDIF RadiativeTransfer

  if (FieldsToAdd > 0 && debug)
    fprintf(stdout, "RadiativeTransferInitialize: Increasing baryon fields "
	    "from %"ISYM" to %"ISYM"\n", OldNumberOfBaryonFields, 
	    OldNumberOfBaryonFields+FieldsToAdd);

  if (OldNumberOfBaryonFields+FieldsToAdd > MAX_DEPTH_OF_HIERARCHY)
    ENZO_FAIL("Exceeds MAX_DEPTH_OF_HIERARCHY.  Please increase and re-compile.");

  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel)
      Temp->GridData->AddFields(TypesToAdd, FieldsToAdd);

  /* Add external boundaries */

  for (i = 0; i < FieldsToAdd; i++) {
    Exterior.AddField(TypesToAdd[i]);
  } // ENDFOR fields

  /* Assign the radiation field DataLabels */

  for (i = 0; i < FieldsToAdd; i++) {
    switch (TypesToAdd[i]) {
    case kphHI:
      DataLabel[OldNumberOfBaryonFields+i] = (char*) kphHIName;
      break;
    case gammaHI:
      DataLabel[OldNumberOfBaryonFields+i] = (char*) gammaHIName;
      break;
    case kphHeI:
      DataLabel[OldNumberOfBaryonFields+i] = (char*) kphHeIName;
      break;
    case gammaHeI:
      DataLabel[OldNumberOfBaryonFields+i] = (char*) gammaHeIName;
      break;
    case kphHeII:
      DataLabel[OldNumberOfBaryonFields+i] = (char*) kphHeIIName;
      break;
    case gammaHeII:
      DataLabel[OldNumberOfBaryonFields+i] = (char*) gammaHeIIName;
      break;
    case kdissH2I:
      DataLabel[OldNumberOfBaryonFields+i] = (char*) kdissH2IName;
      break;
    case RadPressure0:
      DataLabel[OldNumberOfBaryonFields+i] = (char*) RadAccel1Name;
      break;
    case RadPressure1:
      DataLabel[OldNumberOfBaryonFields+i] = (char*) RadAccel2Name;
      break;
    case RadPressure2:
      DataLabel[OldNumberOfBaryonFields+i] = (char*) RadAccel3Name;
      break;
    case Metallicity:
      DataLabel[OldNumberOfBaryonFields+i] = (char*) MetalName;
      break;
    case SNColour:
      DataLabel[OldNumberOfBaryonFields+i] = (char*) ColourName;
      break;
    } // ENDSWITCH
  } // ENDFOR fields

  /* Initialize SubgridMarker (do we need to do this?  it's already
     done in RebuildHierarchy) */

  /* Initialize HEALPix arrays */

  if (RadiativeTransfer) {
    pix2x = new long[1024];
    pix2y = new long[1024];
    mkPix2xy(&pix2x[0],&pix2y[0]);

    x2pix = new int[128];
    y2pix = new int[128];
    mk_xy2pix(&x2pix[0], &y2pix[0]);
  }

  // if using an implicit RT solver, declare the appropriate object here

  if (RadiativeTransferFLD == 1)
    ImplicitSolver = new FSProb; 
  else
    ImplicitSolver = new NullProblem;

  // if using the FLD solver, initialize it here
#ifdef USE_HYPRE
  if (RadiativeTransferFLD == 1) {
    ImplicitSolver = new FSProb;
    if (DetermineParallelism(&TopGrid, MetaData) == FAIL) {
      fprintf(stderr,"Error in DetermineParallelism.\n");
      my_exit(EXIT_FAILURE);
    }
  } else {
    ImplicitSolver = new NullProblem;
  }
  ImplicitSolver->Initialize(TopGrid, MetaData);
#else
  if (RadiativeTransferFLD == 1)
    ENZO_FAIL("Error: cannot use RadiativeTransferFLD without HYPRE.");
#endif

  return SUCCESS;

}
