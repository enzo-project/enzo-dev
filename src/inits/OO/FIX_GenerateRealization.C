/***********************************************************************
/
/  GENERATES THE FIELD AND PARTICLE REALIZATIONS
/
/  written by: Greg Bryan
/  date:       June, 1997
/  modified1:  Robert Harkness, July 2002
/
/  PURPOSE:
/
************************************************************************/
 
#include <stdlib.h>
#include <stdio.h>
 
#include "macros_and_parameters.h"
#include "global_data.h"
#include "CosmologyParameters.h"
#include "Parameters.h"
 
// function prototypes
 
extern "C" void FORTRAN_NAME(set_common)(FLOAT *lam0_in, FLOAT *omega0_in,
					 FLOAT *zri_in, FLOAT *hub_in);
extern "C" FLOAT FORTRAN_NAME(calc_f)(FLOAT *aye);
extern "C" FLOAT FORTRAN_NAME(calc_ayed)(FLOAT *aye);
 
int GenerateField(int Rank, int Dims[3], int MaxDims[3], int WaveNumberCutoff,
		  FLOAT *Field, int FieldType, int NewCenter[3],
		  int Refinement, int StartIndex[3], int Species);
int WriteField(int Rank, int Dims[3], FLOAT *Field, char *Name, int Part, int Npart,
               int GridRank, int Starts[3], int Ends[3], int RootGridDims[3]);
int WriteIntField(int Rank, int Dims[3], int *Field, char *Name, int Part, int Npart,
                  int GridRank, int Starts[3], int Ends[3], int RootGridDims[3]);
 
void fcol(FLOAT *x, int n, int m, FILE *log_fptr);
 
void RemoveSubGridParticles(FLOAT *From, FLOAT *To, int Dims[],
			   int Start[], int End[]);
 
 
 
int GenerateRealization(parmstruct *Parameters, parmstruct *SubGridParameters)
{
 
  FLOAT *ParticleField, *GridField, LeftEdge[3], RightEdge[3], ParticleOffset;
  int *ParticleTypeField;
 
  FLOAT GrowthFunction, aye = 1.0, ayed, Temp;
  int i, j, k, dim, size, index, NumberOfParticles, MaxDimsThisLevel[3];
  long_int big;
 
  FILE *dumpfile;
 
#ifdef DUMP_OK
  int         dump_ok = 1;
#else
  int         dump_ok = 0;
#endif
 
  if (dump_ok) dumpfile = fopen("DumpGR","a");
 
  // Calculate some cosmological quantities for later use
 
  FORTRAN_NAME(set_common)(&OmegaLambdaNow, &OmegaMatterNow, &InitialRedshift,
			   &HubbleConstantNow);
  GrowthFunction = FORTRAN_NAME(calc_f)(&aye);  /* dlog(D)/dlog(a) */
  ayed           = FORTRAN_NAME(calc_ayed)(&aye);
 
  if (debug) printf("GrowthFunction: dlog(D+)/dlog(a) = %"GSYM"\n", GrowthFunction);
 
  // Compute the position of the left corner of the particle grid
 
  for (dim = 0; dim < Parameters->Rank; dim++) {
    if (Parameters->NewCenter[0] != INT_UNDEFINED && Parameters->MaxDims[dim]
	!= Parameters->ParticleDims[dim]*Parameters->ParticleRefinement)
      LeftEdge[dim] =
	FLOAT((Parameters->StartIndex[dim] + Parameters->MaxDims[dim]/2-1 -
	       Parameters->NewCenter[dim] + Parameters->MaxDims[dim]      )
	      % Parameters->MaxDims[dim]) /
	FLOAT(Parameters->MaxDims[dim]);
    else
      LeftEdge[dim] = FLOAT(Parameters->StartIndex[dim]) /
	              FLOAT(Parameters->MaxDims[dim]);
    RightEdge[dim] = LeftEdge[dim] +
      FLOAT(Parameters->ParticleDims[dim]*Parameters->ParticleRefinement)/
      FLOAT(Parameters->MaxDims[dim]);
  }
 
  if (debug) {
    printf("ParticleSubgridLeftEdge  = %"FSYM" %"FSYM" %"FSYM"\n", LeftEdge[0], LeftEdge[1],
	   LeftEdge[2]);
    printf("ParticleSubgridRightEdge = %"FSYM" %"FSYM" %"FSYM"\n", RightEdge[0],
	   RightEdge[1], RightEdge[2]);
  }
 
  /* ------------------------------------------------------------------- */
 
  // Set particles
 
  if (Parameters->InitializeParticles) {
 
    // Compute required size and allocate field
 
    size = 1, NumberOfParticles = 1;
 
    for (dim = 0; dim < Parameters->Rank; dim++) {
      size *= (Parameters->ParticleDims[dim] + ((dim == 0) ? 2 : 0));
      NumberOfParticles *= Parameters->ParticleDims[dim];
    }
 
#ifdef RH_M
    big = ((long_int) (size)) * ((long_int) (sizeof(FLOAT)));
    printf("Allocate %lld bytes for ParticleField\n", big);
#endif
 
    if ((ParticleField = new FLOAT[size]) == NULL) {
      fprintf(stderr, "GenerateRealization: malloc failure (%"ISYM").\n", size);
      exit(EXIT_FAILURE);
    }
 
    // Find out where to remove particles if SubGridParameters is set
 
    int SubStart[3], SubEnd[3], ParticlesRemoved = 1, Shift[3];

    if (SubGridParameters) {
      for (dim = 0; dim < Parameters->Rank; dim++) {
        printf("Parameters->MaxDims[dim] = %d  SubGridParameters->MaxDims[dim] = %d\n", Parameters->MaxDims[dim], SubGridParameters->MaxDims[dim]);
      }
    }

    if (SubGridParameters) {
      for (dim = 0; dim < Parameters->Rank; dim++) {
	if (Parameters->MaxDims[dim] != SubGridParameters->MaxDims[dim]) {
	  fprintf(stderr, "SubGridParameter MaxDims must match.\n");
	  exit(EXIT_FAILURE);
	}
	if (SubGridParameters->StartIndex[dim] % Parameters->ParticleRefinement != 0) {
	  fprintf(stderr, "SubGrid StartIndex must be divisible by refinement.\n");
	  exit(EXIT_FAILURE);
	}
 
	// If this is the top grid, the corner is shifted if recentering
 
	if (Parameters->MaxDims[dim] ==
	    Parameters->ParticleDims[dim]*Parameters->ParticleRefinement &&
	    Parameters->NewCenter[dim] != INT_UNDEFINED)
	  Shift[dim] = Parameters->NewCenter[dim] -
	    (Parameters->MaxDims[dim]/2 - 1);
	else
	  Shift[dim] = 0;
 
	printf("Shift[%"ISYM"] = %"ISYM"\n", dim, Shift[dim]);
 
	// Compute start and end indices subgrid region
 
	SubStart[dim] = ((SubGridParameters->StartIndex[dim]-
			         Parameters->StartIndex[dim] - Shift[dim])
			 % Parameters->MaxDims[dim])/ Parameters->ParticleRefinement;
 
	MaxDimsThisLevel[dim] = Parameters->MaxDims[dim]/ Parameters->ParticleRefinement;
 
	SubStart[dim] = (SubStart[dim] + MaxDimsThisLevel[dim]) %
	                MaxDimsThisLevel[dim];
 
	SubEnd[dim] = ((SubGridParameters->StartIndex[dim] -
		               Parameters->StartIndex[dim] +
		        SubGridParameters->ParticleDims[dim]*
		        SubGridParameters->ParticleRefinement - Shift[dim])
		       % Parameters->MaxDims[dim])/ Parameters->ParticleRefinement - 1;
 
	SubEnd[dim] = (SubEnd[dim] + MaxDimsThisLevel[dim]) %
	              MaxDimsThisLevel[dim];
 
	ParticlesRemoved *= SubEnd[dim] - SubStart[dim] + 1;
 
      } // end: loop over dims
 
      NumberOfParticles -= ParticlesRemoved;
 
      if (debug)
	printf("Removing Particle Region = %"ISYM" %"ISYM" %"ISYM" -> %"ISYM" %"ISYM" %"ISYM"\n",
	       SubStart[0], SubStart[1], SubStart[2],
	       SubEnd[0], SubEnd[1], SubEnd[2]);
    }
 
    if (debug) printf("NumberOfParticles = %"ISYM"\n", NumberOfParticles);
 
      // Loop over dimensions
 
    for (dim = 0; dim < Parameters->Rank; dim++) {
      if (debug) printf("GenerateRealization: particle dim %"ISYM".\n", dim);
 
      /* 1) velocities
	    -generate the displacement field (f delta_k -i vec(k)/k^2).
	    -multiply by adot to get velocity */
	
      GenerateField(Parameters->Rank, Parameters->ParticleDims,
		    Parameters->MaxDims, Parameters->WaveNumberCutoff,
		    ParticleField, 1+dim, Parameters->NewCenter,
		    Parameters->ParticleRefinement, Parameters->StartIndex, 1);
 
      Temp = ayed * GrowthFunction;
 
      for (i = 0; i < size; i++)
	ParticleField[i] *= Temp;
 
      /* Remove subgrid particles if necessary (using temp field),
         and write out field. */
 
#ifdef RH_M
      if (SubGridParameters)
      {
        big = ((long_int) (size) * (long_int) (sizeof(FLOAT)));
        printf("Allocate %lld bytes for TempField\n", big);
      }
#endif
 
      if (SubGridParameters) {
	FLOAT *TempField = new FLOAT[size];
	RemoveSubGridParticles(ParticleField, TempField,
			       Parameters->ParticleDims, SubStart, SubEnd);
	WriteField(1, &NumberOfParticles,
		   TempField, Parameters->ParticleVelocityName, dim, 3,
                   Parameters->Rank,
                   Parameters->TopGridStart,
                   Parameters->TopGridEnd,
                   Parameters->RootGridDims);
	delete TempField;
      }
      else
	WriteField(1, &NumberOfParticles,
		   ParticleField, Parameters->ParticleVelocityName, dim, 3,
                   Parameters->Rank,
                   Parameters->TopGridStart,
                   Parameters->TopGridEnd,
                   Parameters->RootGridDims);
 
	
      /* 2) Make the position field by converting velocity to displacement
	    and adding the initial position. */
 
      Temp = 1.0/ayed;
 
      for (i = 0; i < size; i++)
	ParticleField[i] *= Temp;
 
//      ParticleOffset = (Parameters->ParticleRefinement ==
//			  Parameters->GridRefinement       ) ? 0.5 : 0.0;
 
        ParticleOffset = 0.5;
 
#ifdef SHIFT_FOR_LARS
        ParticleOffset = 1.0;
#endif /* SHIFT_FOR_LARS */
 
      if (debug) printf("ParticleOffset = %"GSYM"\n", ParticleOffset);
 
      for (k = 0; k < Parameters->ParticleDims[2]; k++)
	for (j = 0; j < Parameters->ParticleDims[1]; j++) {
	  index = (k*Parameters->ParticleDims[1] + j)*
	          Parameters->ParticleDims[0];
	  Temp = FLOAT(Parameters->ParticleRefinement) /
	         FLOAT(Parameters->MaxDims[dim]);
	  if (dim == 0)
	    for (i = 0; i < Parameters->ParticleDims[0]; i++)
	      ParticleField[index+i] += LeftEdge[dim] +
		                 (FLOAT(i)+ParticleOffset)*Temp;
	  if (dim == 1)
	    for (i = 0; i < Parameters->ParticleDims[0]; i++)
	      ParticleField[index+i] += LeftEdge[dim] +
		                 (FLOAT(j)+ParticleOffset)*Temp;
	  if (dim == 2)
	    for (i = 0; i < Parameters->ParticleDims[0]; i++)
	      ParticleField[index+i] += LeftEdge[dim] +
		                 (FLOAT(k)+ParticleOffset)*Temp;
	  for (i = 0; i < Parameters->ParticleDims[0]; i++) {
	    if (ParticleField[index+i] <  0.0) ParticleField[index+i] += 1.0;
	    if (ParticleField[index+i] >= 1.0) ParticleField[index+i] -= 1.0;
	  }
	}
 
      /* Remove subgrid particles if necessary (using temp field),
         and write out field. */
 
#ifdef RH_M
      if (SubGridParameters)
      {
        big = ((long_int) (size) * (long_int) (sizeof(FLOAT)));
        printf("Allocate %lld bytes for TempField\n", big);
      }
#endif
 
      if (SubGridParameters) {
	FLOAT *TempField = new FLOAT[size];
	RemoveSubGridParticles(ParticleField, TempField,
			       Parameters->ParticleDims, SubStart, SubEnd);
	WriteField(1, &NumberOfParticles,
		   TempField, Parameters->ParticlePositionName, dim, 3,
                   Parameters->Rank,
                   Parameters->TopGridStart,
                   Parameters->TopGridEnd,
                   Parameters->RootGridDims);
	delete TempField;
      }
      else
	WriteField(1, &NumberOfParticles,
		   ParticleField, Parameters->ParticlePositionName, dim, 3,
                   Parameters->Rank,
                   Parameters->TopGridStart,
                   Parameters->TopGridEnd,
                   Parameters->RootGridDims);
 
    } // end: loop over dims
 
    // Generate and write out particle mass field
 
    if (Parameters->ParticleMassName != NULL) {
      FLOAT ParticleMass = (OmegaMatterNow-OmegaBaryonNow)/OmegaMatterNow;
      for (dim = 0; dim < Parameters->Rank; dim++)
	ParticleMass *= Parameters->ParticleRefinement/
	  Parameters->GridRefinement;
      for (i = 0; i < NumberOfParticles; i++)
	ParticleField[i] = ParticleMass;
      WriteField(1, &NumberOfParticles,
                 ParticleField,
		 Parameters->ParticleMassName, 0, 1,
                 Parameters->Rank,
                 Parameters->TopGridStart,
                 Parameters->TopGridEnd,
                 Parameters->RootGridDims);
    }
 
    // Free field
 
    delete ParticleField;

 
    if (Parameters->ParticleTypeName != NULL) {

      // Generate and write out default (dark matter) type field

      if ((ParticleTypeField = new int[size]) == NULL) {
        fprintf(stderr, "GenerateRealization: malloc failure (%"ISYM").\n", size);
        exit(EXIT_FAILURE);
      }

      for (i = 0; i < NumberOfParticles; i++)
        ParticleTypeField[i] = PARTICLE_TYPE_DARK_MATTER;
 
      //  BOGUS PARTICLE ASSIGNMENT FOR TESTING ONLY!!!!
      //      for (i = 0; i < NumberOfParticles; i = i + 9)
      //        ParticleTypeField[i] = PARTICLE_TYPE_TRACER;
 
      WriteIntField(1, &NumberOfParticles,
                 ParticleTypeField,
                 Parameters->ParticleTypeName, 0, 1,
                 Parameters->Rank,
                 Parameters->TopGridStart,
                 Parameters->TopGridEnd,
                 Parameters->RootGridDims);

      // Free field

      delete [] ParticleTypeField;

    }
 

  } // end: if (InitializeParticles)
 
  /* ------------------------------------------------------------------- */
 
  // Set grids
 
  if (Parameters->InitializeGrids) {
 
    // Compute required size and allocate field
 
    size = 1;
 
    for (dim = 0; dim < Parameters->Rank; dim++)
      size *= (Parameters->GridDims[dim] + ((dim == 0) ? 2 : 0));
 
#ifdef RH_M
    big = ((long_int) (size) * (long_int) (sizeof(FLOAT)));
    printf("Allocate %lld bytes for GridField\n", big);
#endif
 
    if ((GridField = new FLOAT[size]) == NULL) {
      fprintf(stderr, "GenerateRealization: malloc failure (%"ISYM").\n", size);
      exit(EXIT_FAILURE);
    }
 
    /* 1) density (add one and multiply by mean density). */
 
    if (debug) printf("Generating grid densities.\n");
    GenerateField(Parameters->Rank, Parameters->GridDims,
		  Parameters->MaxDims, Parameters->WaveNumberCutoff,
		  GridField, 0, Parameters->NewCenter,
		  Parameters->GridRefinement, Parameters->StartIndex, 2);
 
    Temp = OmegaBaryonNow/OmegaMatterNow;
    for (i = 0; i < size; i++)
      GridField[i] = max(GridField[i] + 1.0, 0.1) * Temp;
 
    WriteField(Parameters->Rank, Parameters->GridDims,
	       GridField, Parameters->GridDensityName, 0, 1,
               Parameters->Rank,
               Parameters->TopGridStart,
               Parameters->TopGridEnd,
               Parameters->RootGridDims);
 
    if (dump_ok) fprintf(dumpfile, "Density, size=%"ISYM"\n", size);
    if (dump_ok) fcol(GridField, size, 8, dumpfile);
 
    /* 2) velocities. */
 
    for (dim = 0; dim < Parameters->Rank; dim++) {
      if (debug) printf("GenerateRealization: grid velocity dim %"ISYM".\n", dim);
 
      /* 1) velocities
	    -generate the displacement field (f delta_k -i vec(k)/k^2).
	    -multiply adot to get velocity */
	
      GenerateField(Parameters->Rank, Parameters->GridDims,
		    Parameters->MaxDims, Parameters->WaveNumberCutoff,
		    GridField, 1+dim, Parameters->NewCenter,
		    Parameters->GridRefinement, Parameters->StartIndex, 2);
 
      Temp = ayed * GrowthFunction;
 
      for (i = 0; i < size; i++)
	GridField[i] *= Temp;
 
      WriteField(Parameters->Rank, Parameters->GridDims,
		 GridField, Parameters->GridVelocityName, dim,3,
                 Parameters->Rank,
                 Parameters->TopGridStart,
                 Parameters->TopGridEnd,
                 Parameters->RootGridDims);
 
      if (dump_ok) fcol(GridField, size, 8, dumpfile);
 
    }
 
  }
 
  if (dump_ok) fclose(dumpfile);
 
  return SUCCESS;
}
 
 
void RemoveSubGridParticles(FLOAT *From, FLOAT *To, int Dims[],
			   int Start[], int End[])
{
  int i, j, k, findex = 0, tindex = 0;
 
  for (k = 0; k < Dims[2]; k++)
    for (j = 0; j < Dims[1]; j++)
      for (i = 0; i < Dims[0]; i++, findex++)
	if (k < Start[2] || k > End[2] ||
	    j < Start[1] || j > End[1] ||
	    i < Start[0] || i > End[0])
	  To[tindex++] = From[findex];
 
}
			
