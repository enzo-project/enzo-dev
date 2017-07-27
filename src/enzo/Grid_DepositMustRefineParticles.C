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


int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);

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

  ParticleBufferSize = MustRefineParticlesBufferSize;
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

  /* Loop over all the particles, using only particles marked as
     must-refine particles. */

  float UniformParticleMass = 0.0;
  if (ProblemType == 30 && MustRefineParticlesCreateParticles == 3)
    UniformParticleMass = OmegaDarkMatterNow / OmegaMatterNow;

  /* Loop over all particles, marking wich ones are must refine
     To add rules, modify number of rules here and add to loop below */
  bool *rules;
  const int NumberOfRules = 3;
  rules = new bool[NumberOfRules];

  // Flag particles as must refine particles
  int *IsParticleMustRefine;
  IsParticleMustRefine = new int[NumberOfParticles];
  for (i = 0; i < NumberOfParticles; i ++){
    IsParticleMustRefine[i] = 1;

    // check particle type and uniform mass
    rules[0] = (ParticleType[i] == PARTICLE_TYPE_MUST_REFINE ||
                ParticleType[i] == PARTICLE_TYPE_MBH) ||
               (ParticleMass[i] < UniformParticleMass);

    // check particle mass greater than minimum mass
    rules[1] = (ParticleMass[i] > MustRefineParticlesMinimumMass);

    // add more rules here
    // Rules for individual star particles and their feedback
    //   if M < SNII thresh: Only refine when within a few dt from end 
    //   if M > SNII thresh: Always refine (if winds are on, otherwise do above)
    //   if WD and a few dt from end of life
    if (STARMAKE_METHOD(INDIVIDUAL_STAR) && STARFEED_METHOD(INDIVIDUAL_STAR)){
      rules[0] = TRUE;
      rules[1] = TRUE;

      float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits, MassUnits;

      if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
                &TimeUnits, &VelocityUnits, Time) == FAIL) {
          ENZO_FAIL("Error in GetUnits.");
      }
      MassUnits = DensityUnits * LengthUnits * LengthUnits * LengthUnits;

      float m_star      = ParticleMass[i] * this->CellWidth[0][0]*this->CellWidth[0][0]*this->CellWidth[0][0] * MassUnits / 1.989E33;
      float end_of_life = ParticleAttribute[1][i] + ParticleAttribute[0][i];

      // within some factor of the end of its life
      bool near_end_of_life = fabs(this->Time - end_of_life) < IndividualStarLifeRefinementFactor * this->dtFixed * POW(2,level); // factor of root grid, estimate root grid

      if (ParticleType[i] == PARTICLE_TYPE_INDIVIDUAL_STAR){
        rules[2] = ( ( IndividualStarStellarWinds) && (m_star > IndividualStarSNIIMassCutoff)  ) || // massive stars always on if winds are on
                   ( (!IndividualStarStellarWinds) && (near_end_of_life)  ) ||  // SNII check if no winds are on
                   ( ( IndividualStarStellarWinds) && (near_end_of_life) && (m_star < IndividualStarSNIIMassCutoff) ); // AGB wind check

      } else if (ABS(ParticleType[i]) == PARTICLE_TYPE_INDIVIDUAL_STAR_WD ||
                 ABS(ParticleType[i]) == PARTICLE_TYPE_INDIVIDUAL_STAR_REMNANT ){
        rules[2] = near_end_of_life;
      } else {
        rules[2] = FALSE; // otherwise this is NOT a must refine particle
      }

    } else{
      rules[2] = TRUE; // must set to true since using AND
    }

    // set flag for this particle
    for (int j = 0; j < NumberOfRules; j++)
      IsParticleMustRefine[i] *= rules[j];

    printf("Checked if we should refine, the answer is %"ISYM"\n", IsParticleMustRefine[i]);
    printf("rule 0 = %"ISYM" rule 1 = %"ISYM" rule 2 = %"ISYM"\n",rules[0],rules[1],rules[2]);
  }

  PFORTRAN_NAME(cic_flag)(IsParticleMustRefine,
	   ParticlePosition[0], ParticlePosition[1], ParticlePosition[2],
	   &GridRank, &NumberOfParticles, FlaggingField,
	   LeftEdge, GridDimension, GridDimension+1, GridDimension+2,
	   &CellSize, &ParticleBufferSize);

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
  if (!(ProblemType == 30 && MustRefineParticlesCreateParticles == 3 &&
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

  return NumberOfFlaggedCells;
 
}
