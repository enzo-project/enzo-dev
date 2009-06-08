/***********************************************************************
/
/  INITS MAIN CODE
/
/  written by: Greg Bryan
/  date:       June, 1997
/  modified1:  Robert Harkness (for 64-bit precision)
/  date:       November, 2003
/
/  PURPOSE:
/    This code generates a gaussian (random phase) realization of a
/    large, discrete field, given the power spectrum of perturbations.
/
************************************************************************/
 
#include <stdlib.h>
#include <stdio.h>
 
#include "macros_and_parameters.h"
#define DEFINE_STORAGE
#include "global_data.h"
#include "CosmologyParameters.h"
#include "PowerSpectrumParameters.h"
#undef DEFINE_STORAGE
#include "Parameters.h"
 
// function prototypes
 
int InterpretCommandLine(int argc, char *argv[], char *myname,
			 char *ParameterFile[], char *SubGridParameterFile[]);
int SetParameterDefaults(parmstruct *Parameters);
int ReadParameterFile(FILE *fptr, parmstruct *Parameters);
int InitializePowerSpectrum();
int MakePowerSpectrumLookUpTable();
int GenerateRealization(parmstruct *Parameters, parmstruct *SubGridParameters);
int CosmologyReadParameters(FILE *fptr);
int ReadPowerSpectrumParameters(FILE *fptr);
 
 
Eint32 main(int argc, char *argv[])
{
 
  // Initialize
 
  debug                        = FALSE;
  char *myname                 = argv[0];
  parmstruct Parameters, *SubGridParameters = NULL;
 
  char *ParameterFile = NULL, *SubGridParameterFile = NULL;
  FILE *fptr;
 
  // Interpret command-line arguments
 
  printf("ENZO Inits v12.0 - June 27th 2004\n\n");
 
  InterpretCommandLine(argc, argv, myname, &ParameterFile, &SubGridParameterFile);
 
  // Set Parameter defaults
 
  SetParameterDefaults(&Parameters);
 
  // Open parameter file
 
  if ((fptr = fopen(ParameterFile, "r")) == NULL) {
    fprintf(stderr, "%s: error opening ParameterFile %s\n", myname, ParameterFile);
    exit(EXIT_FAILURE);
  }
 
  // Read parameters from ParameterFile
 
  if (ReadParameterFile(fptr, &Parameters) == FAIL) {
    fprintf(stderr, "Error while reading ParameterFile %s\n", ParameterFile);
    exit(EXIT_FAILURE);
  }
 
  // Read cosmology parameters
 
  rewind(fptr);
 
  if (CosmologyReadParameters(fptr) == FAIL) {
    fprintf(stderr, "Error in ReadCosmologyParameters.\n");;
    return FAIL;
  }
 
  // Read power spectrum parameters
 
  rewind(fptr);
 
  if (ReadPowerSpectrumParameters(fptr) == FAIL) {
    fprintf(stderr, "Error in ReadCosmologyParameters.\n");;
    return FAIL;
  }
 
  // Close parameter file
 
  fclose(fptr);
 
  // Read parameters from SubGridParameter file
 
  if (SubGridParameterFile != NULL) {
    if ((fptr = fopen(SubGridParameterFile, "r")) == NULL) {
      fprintf(stderr, "%s: error opening SubGridParameterFile %s\n", myname,
	      SubGridParameterFile);
      exit(EXIT_FAILURE);
    }
    SubGridParameters = new parmstruct;
    SetParameterDefaults(SubGridParameters);
    ReadParameterFile(fptr, SubGridParameters);
    fclose(fptr);
  }
 
  // Initialize the power spectrum (and set amplitude) at z=0
 
  Redshift = 0.0;
  InitializePowerSpectrum();
 
  // Generate a look-up table at the initial redshift
 
  Redshift = InitialRedshift;
  MakePowerSpectrumLookUpTable();
 
  // Generate the fields and particles
 
  GenerateRealization(&Parameters, SubGridParameters);
 
  printf("successful completion.\n");
  exit(EXIT_SUCCESS);
 
}
