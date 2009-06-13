/***********************************************************************
/
/  INITIALIZE TURBULENT CLOUD
/
/  written by: Peng Wang
/  date:       September, 2007
/  modified1:
/
/
************************************************************************/

#include <string.h>
#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "TopGridData.h"

void WriteListOfFloats(FILE *fptr, int N, float floats[]);
void WriteListOfFloats(FILE *fptr, int N, FLOAT floats[]);
void AddLevel(LevelHierarchyEntry *Array[], HierarchyEntry *Grid, int level);
int RebuildHierarchy(TopGridData *MetaData,
		     LevelHierarchyEntry *LevelArray[], int level);
int GetUnits(float *DensityUnits, float *LengthUnits,
		      float *TemperatureUnits, float *TimeUnits,
		      float *VelocityUnits, FLOAT Time);

int TurbulenceInitialize(FILE *fptr, FILE *Outfptr, 
			 HierarchyEntry &TopGrid, TopGridData &MetaData)
{
  char *DensName = "Density";
  char *TEName   = "TotalEnergy";
  char *GEName   = "GasEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  char *BxName = "Bx";
  char *ByName = "By";
  char *BzName = "Bz";
  char *PhiName = "Phi";
  char *ColourName = "colour";
  char *ElectronName = "Electron_Density";
  char *HIName    = "HI_Density";
  char *HIIName   = "HII_Density";
  char *HeIName   = "HeI_Density";
  char *HeIIName  = "HeII_Density";
  char *HeIIIName = "HeIII_Density";
  char *HMName    = "HM_Density";
  char *H2IName   = "H2I_Density";
  char *H2IIName  = "H2II_Density";
  char *DIName    = "DI_Density";
  char *DIIName   = "DII_Density";
  char *HDIName   = "HDI_Density";
  char *kphHIName    = "HI_kph";
  char *gammaHIName  = "HI_gamma";
  char *kphHeIName   = "HeI_kph";
  char *gammaHeIName = "HeI_gamma";
  char *kphHeIIName  = "HeII_kph";
  char *gammaHeIIName= "HeII_gamma";
  char *kdissH2IName = "H2I_kdiss";
  char *RadAccel1Name = "RadAccel1";
  char *RadAccel2Name = "RadAccel2";
  char *RadAccel3Name = "RadAccel3";
  char *Drive1Name = "DrivingField1";
  char *Drive2Name = "DrivingField2";
  char *Drive3Name = "DrivingField3";
  char *GravPotenName = "PotentialField";
  char *Acce1Name = "AccelerationField1";
  char *Acce2Name = "AccelerationField2";
  char *Acce3Name = "AccelerationField3";

  /* declarations */

  char  line[MAX_LINE_LENGTH];
  int   dim, ret, level, sphere, i;

  /* set default parameters */

  int RefineAtStart   = TRUE;
  int SetTurbulence = TRUE;
  int RandomSeed = 52761;
  float CloudDensity=1.0, CloudSoundSpeed=1.0, CloudMachNumber=1.0, CloudAngularVelocity = 0.0, InitialBField = 0.0;
  FLOAT CloudRadius = 0.05;
  int CloudType = 1;

  /* read input from parameter file */

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    ret += sscanf(line, "RefineAtStart = %d", &RefineAtStart);
    ret += sscanf(line, "Density = %f", &CloudDensity);
    ret += sscanf(line, "SoundVelocity = %f", &CloudSoundSpeed);
    ret += sscanf(line, "MachNumber = %f", &CloudMachNumber);
    ret += sscanf(line, "AngularVelocity = %f", &CloudAngularVelocity);
    ret += sscanf(line, "CloudRadius = %"FSYM, &CloudRadius);
    ret += sscanf(line, "SetTurbulence = %d", &SetTurbulence);
    ret += sscanf(line, "RandomSeed = %d", &RandomSeed);
    ret += sscanf(line, "InitialBfield = %f", &InitialBField);
    ret += sscanf(line, "CloudType = %d", &CloudType);

  }

  /* Convert to code units */
  
  float DensityUnits = 1.0, LengthUnits = 1.0, TemperatureUnits = 1.0, TimeUnits = 1.0, VelocityUnits = 1.0, 
    PressureUnits = 1.0, MagneticUnits = 1.0;
  if (UsePhysicalUnit) 
    GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, MetaData.Time);
  PressureUnits = DensityUnits * pow(VelocityUnits,2);
  MagneticUnits = sqrt(PressureUnits*4.0*Pi);
      
  CloudDensity /= DensityUnits;
  CloudSoundSpeed /= VelocityUnits;
  InitialBField /= MagneticUnits;
  CloudAngularVelocity *= TimeUnits;

  printf("Plasma beta=%g\n", CloudDensity*CloudSoundSpeed*CloudSoundSpeed/(InitialBField*InitialBField/2.0));
  printf("DensityUnits=%g,VelocityUnits=%g,LengthUnits=%g,TimeUnits=%g (%g yr),PressureUnits=%g\n", 
	 DensityUnits, VelocityUnits, LengthUnits, TimeUnits, TimeUnits/3.1558e7, PressureUnits);
  printf("CloudDensity=%g, CloudSoundSpeed=%g, CloudRadius=%g, CloudAngularVelocity=%g\n", 
	 CloudDensity, CloudSoundSpeed, CloudRadius, CloudAngularVelocity);

  /* Begin grid initialization */

  if (TopGrid.GridData->TurbulenceInitializeGrid(
                CloudDensity, CloudSoundSpeed, CloudRadius, CloudMachNumber, CloudAngularVelocity, InitialBField, 
		SetTurbulence, CloudType, RandomSeed, 0) == FAIL) {
    fprintf(stderr, "Error in CollapseTestInitializeGrid.\n");
    return FAIL;
  }

  /* Convert minimum initial overdensity for refinement to mass
     (unless MinimumMass itself was actually set). */

  if (MinimumMassForRefinement[0] == FLOAT_UNDEFINED) {
    MinimumMassForRefinement[0] = MinimumOverDensityForRefinement[0];
    for (int dim = 0; dim < MetaData.TopGridRank; dim++)
      MinimumMassForRefinement[0] *=(DomainRightEdge[dim]-DomainLeftEdge[dim])/
	float(MetaData.TopGridDims[dim]);
  }

  /* If requested, refine the grid to the desired level. */

  if (RefineAtStart) {

    /* Declare, initialize and fill out the LevelArray. */

    LevelHierarchyEntry *LevelArray[MAX_DEPTH_OF_HIERARCHY];
    for (level = 0; level < MAX_DEPTH_OF_HIERARCHY; level++)
      LevelArray[level] = NULL;
    AddLevel(LevelArray, &TopGrid, 0);

    /* Add levels to the maximum depth or until no new levels are created,
       and re-initialize the level after it is created. */

    for (level = 0; level < MaximumRefinementLevel; level++) {
      if (RebuildHierarchy(&MetaData, LevelArray, level) == FAIL) {
	fprintf(stderr, "Error in RebuildHierarchy.\n");
	return FAIL;
      }
      if (LevelArray[level+1] == NULL)
	break;

      LevelHierarchyEntry *Temp = LevelArray[level+1];
      while (Temp != NULL) {
	if (Temp->GridData->TurbulenceInitializeGrid(
		  CloudDensity, CloudSoundSpeed, CloudRadius, CloudMachNumber, CloudAngularVelocity, InitialBField,
		  SetTurbulence, CloudType, RandomSeed, level+1) == FAIL) {
	  fprintf(stderr, "Error in Collapse3DInitializeGrid.\n");
	  return FAIL;
	}
	Temp = Temp->NextGridThisLevel;
      }
    } // end: loop over levels

    /* Loop back from the bottom, restoring the consistency among levels. */

    for (level = MaximumRefinementLevel; level > 0; level--) {
      LevelHierarchyEntry *Temp = LevelArray[level];
      while (Temp != NULL) {
	if (Temp->GridData->ProjectSolutionToParentGrid(
				   *LevelArray[level-1]->GridData) == FAIL) {
	  fprintf(stderr, "Error in grid->ProjectSolutionToParentGrid.\n");
	  return FAIL;
	}
	Temp = Temp->NextGridThisLevel;
      }
    }

  } // end: if (RefineAtStart)

  /* set up field names and units */

  int count = 0;
  DataLabel[count++] = DensName;
  DataLabel[count++] = Vel1Name;
  DataLabel[count++] = Vel2Name;
  DataLabel[count++] = Vel3Name;
  DataLabel[count++] = TEName;
  if (DualEnergyFormalism) {
    DataLabel[count++] = GEName;
  }
  if (HydroMethod == MHD_RK) {
    DataLabel[count++] = BxName;
    DataLabel[count++] = ByName;
    DataLabel[count++] = BzName;
    DataLabel[count++] = PhiName;
  }
  if (MultiSpecies) {
    DataLabel[count++] = ElectronName;
    DataLabel[count++] = HIName;
    DataLabel[count++] = HIIName;
    DataLabel[count++] = HeIName;
    DataLabel[count++] = HeIIName;
    DataLabel[count++] = HeIIIName;
    if (MultiSpecies > 1) {
      DataLabel[count++] = HMName;
      DataLabel[count++] = H2IName;
      DataLabel[count++] = H2IIName;
    }
    if (MultiSpecies > 2) {
      DataLabel[count++] = DIName;
      DataLabel[count++] = DIIName;
      DataLabel[count++] = HDIName;
    }
  }  // if Multispecies                                                                                                      
  //if (PhotonTestUseColour)
  //DataLabel[count++] = ColourName;
#ifdef TRANSFER
  if (RadiativeTransfer)
    if (MultiSpecies) {
      DataLabel[count++]  = kphHIName;
      DataLabel[count++]  = gammaHIName;
      DataLabel[count++]  = kphHeIName;
      DataLabel[count++]  = gammaHeIName;
      DataLabel[count++]  = kphHeIIName;
      DataLabel[count++]  = gammaHeIIName;
      if (MultiSpecies > 1)
        DataLabel[count++]= kdissH2IName;
    } // if RadiativeTransfer                                                                                                

  if (RadiationPressure) {
    DataLabel[count++]  = RadAccel1Name;
    DataLabel[count++]  = RadAccel2Name;
    DataLabel[count++]  = RadAccel3Name;
  }
#endif
  if (UseDrivingField) {
    DataLabel[count++] = Drive1Name;
    DataLabel[count++] = Drive2Name;
    DataLabel[count++] = Drive3Name;
  }
  if (WritePotential) {
    DataLabel[count++] = GravPotenName;
    DataLabel[count++] = Acce1Name;
    DataLabel[count++] = Acce2Name;
    DataLabel[count++] = Acce3Name;
  }

  for (i = 0; i < count; i++) {
    DataUnits[i] = NULL;
  }

  /* If streaming movie output, write header file. */

//   FILE *header;
//   char *headerName = (char*) "movieHeader.dat";
//   int sizeOfFLOAT = sizeof(FLOAT);
//   int sizeOfRecord = (7+MAX_MOVIE_FIELDS)*sizeof(int) + sizeof(float) +
//     6*sizeof(FLOAT);
//   char *movieVersion = (char*) "1.3";
//   int nMovieFields = 0;
//   while (MovieDataField[nMovieFields] != INT_UNDEFINED &&
//          nMovieFields < MAX_MOVIE_FIELDS) nMovieFields++;

//   if (MovieSkipTimestep != INT_UNDEFINED) {
//     if ((header = fopen(headerName, "w")) == NULL) {
//       fprintf(stderr, "Error in opening movie header.\n");
//       return FAIL;
//     }
//     fprintf(header, "MovieVersion = %s\n", movieVersion);
//     fprintf(header, "RootReso = %d\n", MetaData.TopGridDims[0]);
//     fprintf(header, "FLOATSize = %d\n", sizeOfFLOAT);
//     fprintf(header, "RecordSize = %d\n", sizeOfRecord);
//     fprintf(header, "NumFields = %d\n", nMovieFields);
//     fprintf(header, "NumCPUs = %d\n", NumberOfProcessors);
//     fprintf(header, "FileStem = %s\n", NewMovieName);
//     fclose(header);
//   } /* END: write movie header file */
//   /* Open Amira Data file, if requested */

//   if (MovieSkipTimestep != INT_UNDEFINED) {
//     char *AmiraFileName = new char[80];
//     char *TextureFileName = new char[80];
//     int *RefineByArray = new int[3];
//     bool error = FALSE;
//     //    hid_t DataType = H5T_NATIVE_FLOAT;
//     char fileID[4], pid[6];
//     float root_dx = 1.0 / MetaData.TopGridDims[0];

//     //    staggering stag = CELL_CENTERED;
//     //    fieldtype field_type = SCALAR;
//     for (dim = 0; dim < MAX_DIMENSION; dim++) RefineByArray[dim] = RefineBy;

//     sprintf(pid, "_P%3.3d", MyProcessorNumber);
//     sprintf(fileID, "%3.3d", NewMovieDumpNumber);

//     strcpy(AmiraFileName, "AmiraData");
//     strcat(AmiraFileName, fileID);
//     strcat(AmiraFileName, pid);
//     strcat(AmiraFileName, ".hdf5");

//     strcpy(TextureFileName, "TextureData");
//     strcat(TextureFileName, fileID);
//     strcat(TextureFileName, pid);
//     strcat(TextureFileName, ".hdf5");

//     int field, nBaryonFields;
//     int nFields = 0;
//     while (MovieDataField[nFields] != INT_UNDEFINED)
//       nFields++;
//     nBaryonFields = TopGrid.GridData->ReturnNumberOfBaryonFields();

//     char **FieldNames = new char*[nFields];
//     for (field = 0; field < nFields; field++) {
//       FieldNames[field] = new char[64];

//       if (MovieDataField[field] != TEMPERATURE_FIELD)
//         strcpy(FieldNames[field], DataLabel[MovieDataField[field]]);
//       else
//         strcpy(FieldNames[field], "Temperature");
//     }
//     if (Movie3DVolumes > 0)
//       MetaData.AmiraGrid.AMRHDF5Create(AmiraFileName, RefineByArray,
//                                        DataType, stag, field_type,
//                                        MetaData.CycleNumber, MetaData.Time,
//                                        0, root_dx, 1,
//                                        (MyProcessorNumber == ROOT_PROCESSOR),
//                                        nFields, (NewMovieParticleOn > 0),
//                                        NumberOfParticleAttributes, FieldNames,
//                                        error);

//     if (Movie2DTextures > 0)
//       MetaData.Textures.AMRHDF5Create(TextureFileName, RefineByArray,
//                                       DataType, stag, field_type,
//                                       MetaData.CycleNumber, MetaData.Time,
//                                       0, root_dx, 1,
//                                       (MyProcessorNumber == ROOT_PROCESSOR),
//                                       nFields, (NewMovieParticleOn > 0),
//                                       NumberOfParticleAttributes, FieldNames,
//                                       error);

//     if (error) {
//       fprintf(stderr, "Error in AMRHDF5Writer.\n");
//       return FAIL;
//     }

//     delete [] AmiraFileName;
//     delete [] TextureFileName;
//     delete [] RefineByArray;
//     for (field = 0; field < nFields; field++)
//       delete [] FieldNames[field];

//   } /* ENDIF Movie */


  return SUCCESS;

}
