#include <stdlib.h>
#include <stdio.h>
#include <string.h>
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

#define RT_ENERGY_BINS 4

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
void ReadListOfFloats(FILE *fptr, int N, float floats[]);
int CreateSourceClusteringTree(int nShine, SuperSourceData *SourceList,
			       LevelHierarchyEntry *LevelArray[]);

int ReadPhotonSources(FILE *fptr, FLOAT CurrentTime)
{

  int i, source, dim, ret;

  int   PhotonTestNumberOfSources=1;
  int   PhotonTestSourceType[MAX_SOURCES];
  int   PhotonTestSourceEnergyBins[MAX_SOURCES];
  double PhotonTestSourceLuminosity[MAX_SOURCES];
  FLOAT PhotonTestSourcePosition[MAX_SOURCES][MAX_DIMENSION];
  float PhotonTestSourceLifeTime[MAX_SOURCES];
  float PhotonTestSourceCreationTime[MAX_SOURCES];
  float PhotonTestSourceRampTime[MAX_SOURCES];
  float *PhotonTestSourceSED[MAX_SOURCES];
  float *PhotonTestSourceEnergy[MAX_SOURCES];
  float PhotonTestSourceOrientation[MAX_SOURCES][MAX_DIMENSION];

  // Set defaults

  for (source = 0; source < MAX_SOURCES; source++) {
    PhotonTestSourceType[source] = Isotropic;
    PhotonTestSourceLuminosity[source] = 0.;
    PhotonTestSourceLifeTime[source] = 0.;
    if (CurrentTime > 0)
      PhotonTestSourceCreationTime[source] = 0.999*CurrentTime;
    else
      PhotonTestSourceCreationTime[source] = -1.0;
    PhotonTestSourceRampTime[source] = 0.;
    PhotonTestSourceEnergyBins[source] = 1;
    PhotonTestSourceSED[source] = NULL;
    PhotonTestSourceEnergy[source] = NULL;
    for (dim=0; dim < MAX_DIMENSION; dim++){
      PhotonTestSourcePosition[source][dim] =
	0.5*(DomainLeftEdge[dim] + DomainRightEdge[dim]);
    }
    PhotonTestSourceOrientation[source][0] =
      PhotonTestSourceOrientation[source][1] = 0.0;
    PhotonTestSourceOrientation[source][2] = 1.0;
  }

  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, 
    VelocityUnits;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, CurrentTime) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }

  char line[MAX_LINE_LENGTH];
  char *numbers;
  char *delims = (char*) " ";
  char *value;
  int count;
  bool EnergyBinsDefined = false;

  /* read input from file */
  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
     ret = 0;
    ret += sscanf(line, "PhotonTestNumberOfSources = %"ISYM,
		  &PhotonTestNumberOfSources);
     if (sscanf(line, "PhotonTestSourceType[%"ISYM"]", &source) > 0) {
      ret += sscanf(line, "PhotonTestSourceType[%"ISYM"] = %"ISYM, &source,
		    &PhotonTestSourceType[source]);
      if (debug)
	fprintf(stdout, "ReadPhotonSources: Reading Parameters of "
		"Source %"ISYM"...\n", source);
    }
    if (sscanf(line, "PhotonTestSourcePosition[%"ISYM"]", &source) > 0)
      ret += sscanf(line, "PhotonTestSourcePosition[%"ISYM"] = %"PSYM" %"PSYM" %"PSYM, 
		    &source, &PhotonTestSourcePosition[source][0],
		    &PhotonTestSourcePosition[source][1],
		    &PhotonTestSourcePosition[source][2]);
    if (sscanf(line, "PhotonTestSourceLuminosity[%"ISYM"]", &source) > 0)
      ret += sscanf(line, "PhotonTestSourceLuminosity[%"ISYM"] = %lf", &source,
		    &PhotonTestSourceLuminosity[source]);
    if (sscanf(line, "PhotonTestSourceCreationTime[%"ISYM"]", &source) > 0)
      ret += sscanf(line, "PhotonTestSourceCreationTime[%"ISYM"] = %"FSYM, &source,
		    &PhotonTestSourceCreationTime[source]);
    if (sscanf(line, "PhotonTestSourceLifeTime[%"ISYM"]", &source) > 0)
      ret += sscanf(line, "PhotonTestSourceLifeTime[%"ISYM"] = %"FSYM, &source,
		    &PhotonTestSourceLifeTime[source]);
    if (sscanf(line, "PhotonTestSourceRampTime[%"ISYM"]", &source) > 0)
      ret += sscanf(line, "PhotonTestSourceRampTime[%"ISYM"] = %"FSYM, &source,
		    &PhotonTestSourceRampTime[source]);
    if (sscanf(line, "PhotonTestSourceOrientation[%"ISYM"]", &source) > 0)
      ret += sscanf(line, "PhotonTestSourceOrientation[%"ISYM"] = %"FSYM" %"FSYM" %"FSYM, 
		    &source, &PhotonTestSourceOrientation[source][0],
		    &PhotonTestSourceOrientation[source][1],
		    &PhotonTestSourceOrientation[source][2]);
    if (sscanf(line, "PhotonTestSourceEnergyBins[%"ISYM"]", &source) > 0) {
      ret += sscanf(line, "PhotonTestSourceEnergyBins[%"ISYM"] = %"ISYM, &source,
		    &PhotonTestSourceEnergyBins[source]);
      EnergyBinsDefined = true;
    }
    if (sscanf(line, "PhotonTestSourceSED[%"ISYM"]", &source) > 0) {
      if (!EnergyBinsDefined)
	ENZO_FAIL("Must define PhotonTestSourceEnergyBins before SED!");
      PhotonTestSourceSED[source] = new float[PhotonTestSourceEnergyBins[source]+1];
      numbers = strstr(line, "=")+2;
      value = strtok(numbers, delims);
      count = 0;
      while (value != NULL) {
	PhotonTestSourceSED[source][count++] = atof(value);
	value = strtok(NULL, delims);
	ret++;
      }
    }
    if (sscanf(line, "PhotonTestSourceEnergy[%"ISYM"]", &source) > 0) {
      if (!EnergyBinsDefined)
	ENZO_FAIL("Must define PhotonTestSourceEnergyBins before Energies!");
      PhotonTestSourceEnergy[source] = new float[PhotonTestSourceEnergyBins[source]+1];
      numbers = strstr(line, "=")+2;
      value = strtok(numbers, delims);
      count = 0;
      while (value != NULL) {
	PhotonTestSourceEnergy[source][count++] = atof(value);
	value = strtok(NULL, delims);
	ret++;
      }
    }

    /* if the line is suspicious, issue a warning */

    if (ret == 0 && strstr(line, "=") &&
	strstr(line, "PhotonTestSource") && line[0] != '#')
      if (MyProcessorNumber == ROOT_PROCESSOR)
	fprintf(stderr, "warning: the following parameter line was not"
		" interpreted:\n%s\n", line);
    
  } // ENDWHILE line

  // Set default SED and energies if not user-defined

  for (source = 0; source < PhotonTestNumberOfSources; source++) {
    if (PhotonTestSourceSED[source] == NULL) {
      PhotonTestSourceSED[source] = new float[PhotonTestSourceEnergyBins[source]];
      PhotonTestSourceSED[source][0] = 1;
      for (i = 1; i < PhotonTestSourceEnergyBins[source]; i++)
	PhotonTestSourceSED[source][i] = 0;
    }
    if (PhotonTestSourceEnergy[source] == NULL) {
      PhotonTestSourceEnergy[source] = new float[PhotonTestSourceEnergyBins[source]];
      PhotonTestSourceEnergy[source][0] = 14.6;
      for (i = 1; i < PhotonTestSourceEnergyBins[source]; i++)
	PhotonTestSourceEnergy[source][i] = 0;
    }
  }

  // normalize SED

  float totSED, sum;
  for (source = 0; source < PhotonTestNumberOfSources; source++) {
    totSED = 0.;  
    for (i=0; i<PhotonTestSourceEnergyBins[source]; i++) 
      totSED += PhotonTestSourceSED[source][i];
    for (i=0; i<PhotonTestSourceEnergyBins[source]; i++) 
      PhotonTestSourceSED[source][i] /= totSED;
  }

  // list head:
  GlobalRadiationSources = new RadiationSourceEntry;
  GlobalRadiationSources->NextSource = NULL;
  GlobalRadiationSources->PreviousSource = NULL;
  for (i=0; i<PhotonTestNumberOfSources; i++) {
    if (debug) fprintf(stdout, "ReadPhotonSources: %"ISYM" %"GSYM" %"GSYM" %"GSYM"\n", 
		       i, PhotonTestSourceLuminosity[i], TimeUnits, LengthUnits);
    PhotonTestSourceLuminosity[i] *= TimeUnits/pow(LengthUnits,3);
    if (debug) fprintf(stdout, "ReadPhotonSources: %"ISYM"  %"GSYM"\n", 
		       i, PhotonTestSourceLuminosity[i]);
    RadiationSourceEntry *RadSources;
    RadSources = new RadiationSourceEntry;
    RadSources->PreviousSource = GlobalRadiationSources;
    RadSources->NextSource     = GlobalRadiationSources->NextSource;
    RadSources->SuperSource    = NULL;
    RadSources->GridID         = INT_UNDEFINED;
    RadSources->GridLevel      = INT_UNDEFINED;
    RadSources->Type           = PhotonTestSourceType[i]; 
    RadSources->Luminosity     = PhotonTestSourceLuminosity[i]; 
    RadSources->LifeTime       = PhotonTestSourceLifeTime[i]; 
    RadSources->RampTime       = PhotonTestSourceRampTime[i]; 
    RadSources->CreationTime   = PhotonTestSourceCreationTime[i]; 
    RadSources->Position       = new FLOAT[3];
    RadSources->Position[0]    = PhotonTestSourcePosition[i][0]; 
    RadSources->Position[1]    = PhotonTestSourcePosition[i][1]; 
    RadSources->Position[2]    = PhotonTestSourcePosition[i][2]; 
    RadSources->EnergyBins     = PhotonTestSourceEnergyBins[i];
    RadSources->Energy         = new float[PhotonTestSourceEnergyBins[i]];
    RadSources->SED            = new float[PhotonTestSourceEnergyBins[i]];
    RadSources->AddedEmissivity = false;
    for (int j=0; j<PhotonTestSourceEnergyBins[i]; j++){
      RadSources->Energy[j] = PhotonTestSourceEnergy[i][j];
      RadSources->SED[j]    = PhotonTestSourceSED[i][j];
    }
    if (RadSources->Type == Beamed) {
      RadSources->Orientation = new float[3];
      sum = 0;  // for normalization
      for (dim = 0; dim < MAX_DIMENSION; dim++) {
	RadSources->Orientation[dim] = PhotonTestSourceOrientation[i][dim];
	sum += PhotonTestSourceOrientation[i][dim] * PhotonTestSourceOrientation[i][dim];
      }
      sum = sqrt(sum);
      for (dim = 0; dim < MAX_DIMENSION; dim++)
	RadSources->Orientation[dim] /= sum;
    } else {
      RadSources->Orientation = NULL;
    }

    if (RadSources->Type != Isotropic && RadSources->Type != Beamed &&
	RadSources->Type != Episodic) {
      if (MyProcessorNumber == ROOT_PROCESSOR)
	fprintf(stderr, "PhotonTestSourceType must be 1, -2, -3.\n",
		"\tChanging to 1 (isotropic)\n");
      RadSources->Type = Isotropic;
    }

    GlobalRadiationSources->NextSource = RadSources;

  }  

  /* Delete allocated memory for temporary (for I/O) photon sources */

  for (source = 0; source < PhotonTestNumberOfSources; source++) {
    delete [] PhotonTestSourceEnergy[source];
    delete [] PhotonTestSourceSED[source];
  }

  /* Create tree that clusters the sources if requested */

  /* While creating tree (type SuperSource), compute position of the
     super source in each leaf. */
  
  if (RadiativeTransferSourceClustering == TRUE)
    if (CreateSourceClusteringTree(0, NULL, NULL) == FAIL) {
      ENZO_FAIL("Error in CreateSourceClusteringTree.\n");

    }
  
  return SUCCESS;

}
