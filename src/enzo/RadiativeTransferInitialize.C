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
#include "RadiativeTransferHealpixRoutines.h"
#include "ImplicitProblemABC.h"
#include "gFLDProblem.h"
#include "gFLDSplit.h"
#include "FSProb.h"
// #include "MFProb.h"
// #include "MFSplit.h"
#include "NullProblem.h"

int RadiativeTransferReadParameters(FILE *fptr);
int ReadPhotonSources(FILE *fptr, FLOAT CurrentTime);
int DetermineParallelism(HierarchyEntry *TopGrid, TopGridData &MetaData);
int InitializeRadiativeTransferSpectrumTable(FLOAT Time);

int RadiativeTransferInitialize(char *ParameterFile, 
				HierarchyEntry &TopGrid, 
				TopGridData &MetaData,
				ExternalBoundary &Exterior, 
				ImplicitProblemABC* &ImplicitSolver,
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
  const char    *MetalIaName = "MetalSNIa_Density";
  const char	*ColourName    = "SN_Colour";
  const char    *RaySegName    = "Ray_Segments";
  const char    *Rad0Name      = "Radiation0";
  const char    *Rad1Name      = "Radiation1";
  const char    *Rad2Name      = "Radiation2";
  const char    *Rad3Name      = "Radiation3";

  int i, j, k, level;
  FILE *fptr;
  LevelHierarchyEntry *Temp;
  int NumberOfObsoleteFields;
  int ObsoleteFields[MAX_NUMBER_OF_BARYON_FIELDS];

  if (RadiativeTransfer == FALSE && RadiativeTransferFLD == FALSE &&  
      StarParticleFeedback == 0) {

    /* Check for radiation fields and delete them */

    NumberOfObsoleteFields = 8;
    ObsoleteFields[0] = kphHI;
    ObsoleteFields[1] = PhotoGamma;
    ObsoleteFields[2] = kphHeI;
    ObsoleteFields[3] = kphHeII;
    ObsoleteFields[4] = gammaHeI;
    ObsoleteFields[5] = gammaHeII;
    ObsoleteFields[6] = kdissH2I;
    ObsoleteFields[7] = RaySegments;

    for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
      for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel)
	Temp->GridData->DeleteObsoleteFields(ObsoleteFields, 
					     NumberOfObsoleteFields);

    Exterior.DeleteObsoleteFields(ObsoleteFields, NumberOfObsoleteFields);

    return SUCCESS;
  }

  /* Read and set parameter values and static radiation sources */

  if ((fptr = fopen(ParameterFile, "r")) == NULL) {
    ENZO_VFAIL("Error opening ParameterFile %s\n", ParameterFile)
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


  // If FLD solver handles LW radiation, set solver type
  if (RadiativeTransferFLD == 1)   ImplicitProblem = 2;

  // if using an implicit RT solver, declare the appropriate object here
  if (RadiativeTransferFLD) {
    if (ImplicitProblem == 1)
      ImplicitSolver = new gFLDProblem; 
    else if (ImplicitProblem == 2)
      ImplicitSolver = new FSProb; 
    else if (ImplicitProblem == 3)
      ImplicitSolver = new gFLDSplit; 
    // else if (ImplicitProblem == 4)    // MFProb has been removed
    //   ImplicitSolver = new MFProb; 
    // else if (ImplicitProblem == 5)    // MFSplit has been disabled
    //   ImplicitSolver = new MFSplit; 
    else
      ImplicitSolver = new NullProblem;
  }

  // if using the FLD solver, initialize it here
#ifdef USE_HYPRE
  if (RadiativeTransferFLD) {
    // first get parallelism information for implicit system
    if (DetermineParallelism(&TopGrid, MetaData) == FAIL)
      ENZO_FAIL("Error in DetermineParallelism.");
    // initialize the implicit solver
    ImplicitSolver->Initialize(TopGrid, MetaData);
  }
#else
  if (RadiativeTransferFLD)
    ENZO_FAIL("Error: cannot use RadiativeTransferFLD without HYPRE.");
#endif


  /* Create all StarParticles from normal particles */

//  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
//    for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel)
//      if (Temp->GridData->FindNewStarParticles(level) == FAIL) {
//	ENZO_FAIL("Error in grid::FindStarParticles.\n");
//      }

  /* If we're restarting from a non-radiative run, create radiation
     fields.  We will check if the field already exists inside the
     grid routine. */

  bool AddedMetallicity = false;
  int OldNumberOfBaryonFields = 0, FieldsToAdd = 0;
  int TypesToAdd[MAX_NUMBER_OF_BARYON_FIELDS];
  int ExistingTypes[MAX_NUMBER_OF_BARYON_FIELDS];

  for (i = 0; i < MAX_NUMBER_OF_BARYON_FIELDS; i++)
    ExistingTypes[i] = FieldUndefined;

  if (RadiativeTransfer || StarParticleFeedback) {

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
      if (StarClusterUseMetalField &&
	  StarParticleFeedback > 0 &&
	  StarParticleFeedback != (1 << POP3_STAR)) {
	TypesToAdd[FieldsToAdd++] = Metallicity;
	AddedMetallicity = true;
      }
      if (RadiativeTransferLoadBalance)
	TypesToAdd[FieldsToAdd++] = RaySegments;
    }
    // Add metallicity if Pop II star feedback
    if (StarParticleFeedback > 0 && 
	StarParticleFeedback != (1 << POP3_STAR) && 
	!AddedMetallicity)
      TypesToAdd[FieldsToAdd++] = Metallicity;      //#####

    if (StarMakerTypeIaSNe && StarParticleFeedback > 0)
      TypesToAdd[FieldsToAdd++] = MetalSNIaDensity;

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

  }  // ENDIF RadiativeTransfer
  else if (RadiativeTransferFLD > 1) {  // do the same for the FLD solver
    TypesToAdd[FieldsToAdd++] = RadiationFreq0;
    if (RadiativeCooling) {
      TypesToAdd[FieldsToAdd++] = kphHI;
      TypesToAdd[FieldsToAdd++] = PhotoGamma;
      if (RadiativeTransferHydrogenOnly == FALSE) {
	TypesToAdd[FieldsToAdd++] = kphHeI;
	TypesToAdd[FieldsToAdd++] = kphHeII;
      }
      if (MultiSpecies > 1)
	TypesToAdd[FieldsToAdd++] = kdissH2I;
    }

    // don't use the rest
    for (i = FieldsToAdd; i < MAX_NUMBER_OF_BARYON_FIELDS; i++)
      TypesToAdd[i] = FieldUndefined;

    // Check if the fields already exist
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

  } // ENDIF RadiativeTransferFLD


  // Add the necessary fields
  if (FieldsToAdd > 0 && debug)
    fprintf(stdout, "RadiativeTransferInitialize: Increasing baryon fields "
	    "from %"ISYM" to %"ISYM"\n", OldNumberOfBaryonFields, 
	    OldNumberOfBaryonFields+FieldsToAdd);

  // Add an extra 1 because we will need it for flagging/marking cells.
  if (OldNumberOfBaryonFields+FieldsToAdd+1 > MAX_NUMBER_OF_BARYON_FIELDS)
    ENZO_FAIL("Exceeds MAX_NUMBER_OF_BARYON_FIELDS.  "
        "Please increase and re-compile.");

  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel)
      Temp->GridData->AddFields(TypesToAdd, FieldsToAdd);

  // Add external boundaries
  for (i = 0; i < FieldsToAdd; i++) {
    Exterior.AddField(TypesToAdd[i]);
  } // ENDFOR fields


  // Assign the radiation field DataLabels
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
    case MetalSNIaDensity:
      DataLabel[OldNumberOfBaryonFields+i] = (char*) MetalIaName;
      break;
    case SNColour:
      DataLabel[OldNumberOfBaryonFields+i] = (char*) ColourName;
      break;
    case RadiationFreq0:
      DataLabel[OldNumberOfBaryonFields+i] = (char*) Rad0Name;
      break;
    case RadiationFreq1:
      DataLabel[OldNumberOfBaryonFields+i] = (char*) Rad1Name;
      break;
    case RadiationFreq2:
      DataLabel[OldNumberOfBaryonFields+i] = (char*) Rad2Name;
      break;
    case RadiationFreq3:
      DataLabel[OldNumberOfBaryonFields+i] = (char*) Rad3Name;
      break;
    case RaySegments:
      DataLabel[OldNumberOfBaryonFields+i] = (char*) RaySegName;
      break;
    } // ENDSWITCH
    if (debug) 
      printf("Adding new field %s\n", DataLabel[OldNumberOfBaryonFields+i]);
  } // ENDFOR fields

  /* Check for old gammaHeI and gammaHeII fields.  Delete if they
     exist. */

  NumberOfObsoleteFields = 2;
  ObsoleteFields[0] = gammaHeI;
  ObsoleteFields[1] = gammaHeII;
  if (RadiativeTransferHydrogenOnly == TRUE) {
    NumberOfObsoleteFields += 2;
    ObsoleteFields[2] = kphHeI;
    ObsoleteFields[3] = kphHeII;
  }
  if (RadiativeTransferLoadBalance == FALSE)
    ObsoleteFields[NumberOfObsoleteFields++] = RaySegments;

  for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
    for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel)
      Temp->GridData->DeleteObsoleteFields(ObsoleteFields, 
					   NumberOfObsoleteFields);

  Exterior.DeleteObsoleteFields(ObsoleteFields, NumberOfObsoleteFields);

  /* Initialize SubgridMarker (do we need to do this?  it's already
     done in RebuildHierarchy) */

  // Initialize HEALPix arrays
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
