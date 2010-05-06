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

int RadiativeTransferReadParameters(FILE *fptr);
int ReadPhotonSources(FILE *fptr, FLOAT CurrentTime);
int InitializeRadiativeTransferSpectrumTable(FLOAT Time);


int RadiativeTransferInitialize(char *ParameterFile, TopGridData &MetaData,
				ExternalBoundary &Exterior, 
				LevelHierarchyEntry *LevelArray[])
{

  const char	*kphHIName     = "HI_kph";
  const char	*gammaName     = "PhotoGamma";
  const char	*kphHeIName    = "HeI_kph";
  const char	*kphHeIIName   = "HeII_kph";
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

  if (RadiativeTransferReadParameters(fptr) == FAIL) {
    fprintf(stderr, "Error in RadiativeTransferReadParameters.\n");;
    ENZO_FAIL("");
  }
  rewind(fptr);
  if (ProblemType == 50)
    if (ReadPhotonSources(fptr, MetaData.Time) == FAIL) {
      fprintf(stderr, "Error in ReadPhotonSources.\n");;
      ENZO_FAIL("");
    }
  PhotonTime = MetaData.Time;

  fclose(fptr);

  if (RadiativeTransferPhotonEscapeRadius > 0) {
    PhotonEscapeFilename = new char[80];
    sprintf(PhotonEscapeFilename, "fesc%4.4d.dat", MetaData.DataDumpNumber-1);
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
    TypesToAdd[FieldsToAdd++] = kphHI;
    TypesToAdd[FieldsToAdd++] = PhotoGamma;
    if (RadiativeTransferHydrogenOnly == FALSE) {
      TypesToAdd[FieldsToAdd++] = kphHeI;
      TypesToAdd[FieldsToAdd++] = kphHeII;
    }
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

  if (OldNumberOfBaryonFields+FieldsToAdd > MAX_NUMBER_OF_BARYON_FIELDS)
    ENZO_FAIL("Exceeds MAX_NUMBER_OF_BARYON_FIELDS.  "
	      "Please increase and re-compile.");

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
    case PhotoGamma:
      DataLabel[OldNumberOfBaryonFields+i] = (char*) gammaName;
      break;
    case kphHeI:
      DataLabel[OldNumberOfBaryonFields+i] = (char*) kphHeIName;
      break;
    case kphHeII:
      DataLabel[OldNumberOfBaryonFields+i] = (char*) kphHeIIName;
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

  /* Check for old gammaHeI and gammaHeII fields.  Delete if they
     exist. */

  int NumberOfObsoleteFields = 2;
  int ObsoleteFields[MAX_NUMBER_OF_BARYON_FIELDS];

  ObsoleteFields[0] = gammaHeI;
  ObsoleteFields[1] = gammaHeII;
  if (RadiativeTransferHydrogenOnly == TRUE) {
    NumberOfObsoleteFields += 2;
    ObsoleteFields[2] = kphHeI;
    ObsoleteFields[3] = kphHeII;
  }

  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel)
      Temp->GridData->DeleteObsoleteFields(ObsoleteFields, 
					   NumberOfObsoleteFields);

  Exterior.DeleteObsoleteFields(ObsoleteFields, NumberOfObsoleteFields);

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

//  fprintf(stderr, "RTI: RTTS = %d, RTTST =  %s\n", 
//	  RadiativeTransferTraceSpectrum, RadiativeTransferTraceSpectrumTable); 

  /* If set, initialize spectrum table */

  if (RadiativeTransfer == TRUE &&
      RadiativeTransferTraceSpectrum == TRUE) {
    if (InitializeRadiativeTransferSpectrumTable(MetaData.Time) == FAIL) {  
      ENZO_FAIL("Error in InitializeRadiativeTransferSpectrumTable.");
    }
  }

  return SUCCESS;

}
