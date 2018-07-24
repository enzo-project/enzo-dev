/***********************************************************************
/
/  GRID CLASS (DEPOSIT MUST-REFINE PARTICLES)
/
/  written by: Greg Bryan
/  date:       February, 1998
/  modified1:  May, 2009 by John Wise
/                Works with local particle storage in RebuildHierarchy,
/                which appropriately distributes memory usage.
/  modified2:  Sep, 2009 by Ji-hoon Kim
/                Now PARTICLE_TYPE_MBH also has the MUST_REFINE feature
/
/  PURPOSE:
/
************************************************************************/
 
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
#include "CosmologyParameters.h"

 
/* function prototypes */
 
extern "C" void PFORTRAN_NAME(cic_flag)(int* irefflag, FLOAT *posx, FLOAT *posy,
			FLOAT *posz, int *ndim, int *npositions,
                        int *ffield, FLOAT *leftedge,
                        int *dim1, int *dim2, int *dim3, FLOAT *cellsize,
			int *buffersize);

int grid::DepositMustRefineParticles(int pmethod, int level, bool KeepFlaggingField)
{
  /* declarations */
  //printf("grid::DepositMustRefineParticles called \n");
  int i, dim, size = 1;
  FLOAT LeftEdge[MAX_DIMENSION], CellSize;
  int ParticleBufferSize;

  ParticleBufferSize = 1;
  if (ProblemType == 106 || ProblemType ==107)
    ParticleBufferSize = 16;


  /* error check */
 
  if (ParticleMassFlaggingField == NULL) {
    fprintf(stderr, "Particle Mass Flagging Field is undefined.\n");
    return -1;
  }

  /* If refining region before supernova (to be safe in its last 5% of
     the lifetime), temporarily set particle type of star to
     PARTICLE_TYPE_MUST_REFINE. */

  if (PopIIISupernovaMustRefine == TRUE)
    this->ChangeParticleTypeBeforeSN(PARTICLE_TYPE_MUST_REFINE, level,
				     &ParticleBufferSize);

  /* Set Left edge of grid. */

  for (dim = 0; dim < GridRank; dim++) {
    LeftEdge[dim] = CellLeftEdge[dim][0];
    size *= GridDimension[dim];
  }

  CellSize = float(CellWidth[0][0]);

  /* Temporarily set the flagging field, then we will increase the
     particle mass flagging field above the refinement criteron. */

  FlaggingField = new int[size];
  for (i = 0; i < size; i++)
    FlaggingField[i] = 0;

  float UniformParticleMass = 0.0;
  if (ProblemType == 30 &&
      (MustRefineParticlesCreateParticles == 3 ||
       MustRefineParticlesCreateParticles == 4))
    UniformParticleMass = OmegaDarkMatterNow / OmegaMatterNow;

  /* Loop over all particles, marking wich ones are must refine
     To add rules, modify number of rules here and add to loop below */
  bool *rules, rule0;
  const int NumberOfRules = 2;
  rules = new bool[NumberOfRules];

  // Rules to prevent refinement, cancelling out the above rules.
  bool *antirules;
  int *AntiFlaggingField;
  int NumberOfAntiRules = 0;
  antirules = new bool[NumberOfAntiRules];

  // Add an antirule to unflag over-refined dark matter particles.
  if (MustRefineParticlesCreateParticles == 4) {
    NumberOfAntiRules++;
  }

  if (NumberOfAntiRules > 0) {
    antirules = new bool[NumberOfAntiRules];
    AntiFlaggingField = new int[size];
    for (i = 0; i < size; i++)
      AntiFlaggingField[i] = 0;
  }

  // Flag particles as must refine particles
  int *IsParticleMustRefine, *IsParticleNotMustRefine;
  IsParticleMustRefine = new int[NumberOfParticles];
  if (NumberOfAntiRules > 0) {
    IsParticleNotMustRefine = new int[NumberOfParticles];
  }
  for (i = 0; i < NumberOfParticles; i ++){
    IsParticleMustRefine[i] = 1;

    // check particle type and uniform mass
    rule0    = (ParticleType[i] == PARTICLE_TYPE_MUST_REFINE ||
                ParticleType[i] == PARTICLE_TYPE_MBH) ||
               (ParticleMass[i] < UniformParticleMass);
    rules[0] = rule0;

    // check particle mass greater than minimum mass
    rules[1] = (ParticleMass[i] > MustRefineParticlesMinimumMass);

    // add more rules here

    //

    // set flag for this particle
    for (int j = 0; j < NumberOfRules; j++)
      IsParticleMustRefine[i] *= rules[j];

    // anti-rules
    if (NumberOfAntiRules > 0) {
      IsParticleNotMustRefine[i] = 1;
      // check for over-refined dark matter particles
      antirules[0] = !rule0;
    }

    // set antiflag for this particle
    for (int j = 0; j < NumberOfAntiRules; j++)
      IsParticleNotMustRefine[i] *= antirules[j];
  }

  PFORTRAN_NAME(cic_flag)(IsParticleMustRefine,
	   ParticlePosition[0], ParticlePosition[1], ParticlePosition[2],
	   &GridRank, &NumberOfParticles, FlaggingField,
	   LeftEdge, GridDimension, GridDimension+1, GridDimension+2,
	   &CellSize, &ParticleBufferSize);

  if (NumberOfAntiRules > 0) {
    PFORTRAN_NAME(cic_flag)(IsParticleNotMustRefine,
	   ParticlePosition[0], ParticlePosition[1], ParticlePosition[2],
	   &GridRank, &NumberOfParticles, AntiFlaggingField,
	   LeftEdge, GridDimension, GridDimension+1, GridDimension+2,
	   &CellSize, &ParticleBufferSize);

    for (i = 0;i < size;i++) {
      FlaggingField[i] *= !(AntiFlaggingField[i]);
    }
  }

  /* Increase particle mass flagging field for definite refinement */

  float MustRefineMass = 
    1.001*MinimumMassForRefinement[pmethod] * 
    POW(RefineBy, level * MinimumMassForRefinementLevelExponent[pmethod]);
  if (ProblemType == 28)
    MustRefineMass = 0;

  /* Special case on level == MustRefineParticlesRefineToLevel when we
     restrict the additional AMR to regions with must-refine
     particles, and don't use the particle mass field. */
  
  int NumberOfFlaggedCells = 0;
  if (!(ProblemType == 30 &&
        (MustRefineParticlesCreateParticles == 3 ||
         MustRefineParticlesCreateParticles == 4) &&
	level == MustRefineParticlesRefineToLevel)) {
    for (i = 0; i < size; i++)
      if (FlaggingField[i]) {
	ParticleMassFlaggingField[i] = MustRefineMass;
	NumberOfFlaggedCells++;
      }
  }

  if (debug1)
    printf("DepositMRPs[%"ISYM"]: %"ISYM" flagged cells\n", 
	   level,NumberOfFlaggedCells);

  /* If refining region before supernova, change particle type back to
     its original value. */

  if (PopIIISupernovaMustRefine == TRUE)
    this->ChangeParticleTypeBeforeSN(PARTICLE_TYPE_RESET, level);

  /* Clean up */

  if (!KeepFlaggingField) {
    delete [] FlaggingField;
    FlaggingField = NULL;
  }

  delete [] IsParticleMustRefine;
  delete [] rules;

  if (NumberOfAntiRules > 0) {
    delete [] AntiFlaggingField;
    delete [] IsParticleNotMustRefine;
    delete [] antirules;
  }

  return NumberOfFlaggedCells;
 
}
