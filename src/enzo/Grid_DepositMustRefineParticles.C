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

#ifdef INDIVIDUALSTAR
// do a separate method to keep things clean, even though a lot will be redundant
int grid::DepositMustRefineParticles(int pmethod, int level, bool KeepFlaggingField,
                                     TopGridData *MetaData, Star *&AllStars){

  /* declarations */
  int i, dim, size = 1;
  FLOAT LeftEdge[MAX_DIMENSION], CellSize;

  if (!(AllStars)) return 0;

  float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1,
    TimeUnits = 1, VelocityUnits = 1;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
               &TimeUnits, &VelocityUnits, Time) == FAIL) {
        ENZO_FAIL("Error in GetUnits.");
  }

  /* error check */

  if (ParticleMassFlaggingField == NULL) {
    fprintf(stderr, "Particle Mass Flagging Field is undefined.\n");
    return -1;
  }

  /* Set Left edge of grid. */

  for (dim = 0; dim < GridRank; dim++) {
    LeftEdge[dim] = CellLeftEdge[dim][0];
    size *= GridDimension[dim];
  }

  CellSize = float(CellWidth[0][0]);

  /* Temporarily set the flagging field, then we will increase the
     particle mass flagging field above the refinement criteron. */

   // Have to do this since this takes place after other must refine
   // methods. Double check we are not initializing this twice
   // and that, if it is initialized, and not saving then we have
   // to reset to zero. Otherwise, it needs to persist
   if (FlaggingField == NULL){
       FlaggingField = new int[size];
   }
   if (!KeepFlaggingField){
     for (i = 0; i < size; i ++){
       FlaggingField[i] = 0;
     }
   }

  /* Loop over all the particles, using only particles marked as
     must-refine particles. */

  /* need to do the following, either:
   *    1) Make flagging field and particle arrays at size of star particles
   *    2) Check if particle position overlaps with grid z
   *    3) then check must refine rules
   *    4) save particles and number of refinement required particles -
   *       send this # to cic so only need to loop over the "yes"'s
   * ---------------------------------------------------------------------- */

  /* AJE: This is a bit of a hack --- want to make sure injection regions
          in metal mixing experiment are refined to highest level AND in the
          same was as it would be done for star particles for consistency. Hack
          this here, rather than write a new routine...
  */
  int num_events = 0;
  if (MetalMixingExperiment){
    // count the number of events that will go off this timestep:
    for (i = 0; i < MAX_TIME_ACTIONS; i++){
      if ((MetaData->Time >= (TimeActionTime[i] - 0.5*3.154E13/TimeUnits))
            && TimeActionTime[i] > 0){
        num_events++;
      }
    }
  }

  int *IsParticleMustRefine;
  FLOAT *StarPosX, *StarPosY, *StarPosZ;
  IsParticleMustRefine = new int[MetaData->NumberOfParticles + num_events];
  StarPosX = new FLOAT[MetaData->NumberOfParticles + num_events];
  StarPosY = new FLOAT[MetaData->NumberOfParticles + num_events];
  StarPosZ = new FLOAT[MetaData->NumberOfParticles + num_events];

/*
  for (i = 0; i < MetaData->NumberOfParticles; i++){
    IsParticleMustRefine[i] = 0;
    StarPosX[i] = -1.;
    StarPosY[i] = -1.;
    StarPosZ[i] = -1.;
  }
*/

  int NumberOfMustRefineStars = 0;
  FLOAT *pos;
  Star *cstar;

  for (cstar = AllStars, i = 0; cstar; cstar = cstar->NextStar) {
    float end_of_life = cstar->ReturnBirthTime() + cstar->ReturnLifetime();
    bool near_end_of_life = fabs(this->Time - end_of_life) < (IndividualStarLifeRefinementFactor * this->dtFixed * POW(2,level)); // factor of root grid, estimate root $

    IsParticleMustRefine[i] = 0;
    StarPosX[i] = StarPosY[i] = StarPosZ[i] = -1.;

    if (cstar->ReturnType() == PARTICLE_TYPE_INDIVIDUAL_STAR){
        if (( ( IndividualStarStellarWinds) && (cstar->ReturnMass() > IndividualStarSNIIMassCutoff)  ) || // massive stars always on if winds are on
                   ( (!IndividualStarStellarWinds) && (near_end_of_life)  ) ||  // SNII check if no winds are on
                   ( ( IndividualStarStellarWinds) && (near_end_of_life) && (cstar->ReturnMass() < IndividualStarSNIIMassCutoff) )){ // AGB wind check

          IsParticleMustRefine[i] = 1;
        }
    } else if (fabs(cstar->ReturnType()) == PARTICLE_TYPE_INDIVIDUAL_STAR_POPIII){
      if (near_end_of_life   &&
          (   ((cstar->ReturnBirthMass()>=TypeIILowerMass)&&(cstar->ReturnBirthMass()<=TypeIIUpperMass)) ||
              ((cstar->ReturnBirthMass()>=PISNLowerMass)  &&(cstar->ReturnBirthMass()<=PISNUpperMass  ))   )){
        IsParticleMustRefine[i] = 1;
      }
    } else if (fabs(cstar->ReturnType()) == PARTICLE_TYPE_INDIVIDUAL_STAR_WD ||
               fabs(cstar->ReturnType()) == PARTICLE_TYPE_INDIVIDUAL_STAR_REMNANT ){
      IsParticleMustRefine[i] = ((int) near_end_of_life);
    }

    if(IsParticleMustRefine[i]){
      pos = cstar->ReturnPosition();
      StarPosX[i] = pos[0];
      StarPosY[i] = pos[1];
      StarPosZ[i] = pos[2];
      i++;
    }
//    else{
//    }
  } // end for

  if (MetalMixingExperiment){

    for (int j = 0; j < MAX_TIME_ACTIONS; j++){
      if ( (MetaData->Time >= (TimeActionTime[j] - 0.5*3.154E13/TimeUnits))
          && TimeActionTime[j] > 0){

        // NOTE: known bug - indeces will not be correct b/t
        //       TimeAction and Data struct if multiple action types used....

        if (TimeActionType[j] == 4){
          IsParticleMustRefine[i] = 1; // still need to set this
          // assuming code units!!!!
          StarPosX[i] = MixingExperimentData.xpos[j];
          StarPosY[i] = MixingExperimentData.ypos[j];
          StarPosZ[i] = MixingExperimentData.zpos[j];
          i++;
        }
      }
    }


  }

  NumberOfMustRefineStars = i; // save number of stars

  PFORTRAN_NAME(cic_flag)(IsParticleMustRefine,
	                  StarPosX, StarPosY, StarPosZ,
      	                  &GridRank, &NumberOfMustRefineStars, FlaggingField,
	                  LeftEdge, GridDimension, GridDimension+1, GridDimension+2,
                          &CellSize, &IndividualStarRefineBufferSize);

  /* Increase particle mass flagging field for definite refinement */

  float MustRefineMass = 
    1.001*MinimumMassForRefinement[pmethod] * 
    POW(RefineBy, level * MinimumMassForRefinementLevelExponent[pmethod]);

  int NumberOfFlaggedCells = 0;
  for (i = 0; i < size; i++){
    if (FlaggingField[i]) {
      ParticleMassFlaggingField[i] = MustRefineMass;
      NumberOfFlaggedCells++;
    }
  }

  if (debug1)
    printf("DepositMRPs[%"ISYM"]: %"ISYM" flagged cells\n", 
	   level,NumberOfFlaggedCells);

  /* Clean up */

  if (!KeepFlaggingField) {
    delete [] FlaggingField;
    FlaggingField = NULL;
  }

  delete [] IsParticleMustRefine;
  delete [] StarPosX;
  delete [] StarPosY;
  delete [] StarPosZ;

  return NumberOfFlaggedCells;

}
#endif

//#else

int grid::DepositMustRefineParticles(int pmethod, int level, bool KeepFlaggingField)
{
  /* declarations */
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
  bool OriginalParticle;
  IsParticleMustRefine = new int[NumberOfParticles];

  if (NumberOfAntiRules > 0) {
    IsParticleNotMustRefine = new int[NumberOfParticles];
  }
  for (i = 0; i < NumberOfParticles; i ++){
    IsParticleMustRefine[i] = 1;
    if (NumberOfParticleAttributes > 0)
      OriginalParticle = (ParticleAttribute[0][i] <= 0.0);
    else
      OriginalParticle = true;

    // check particle type and uniform mass. Also check particle
    // creation time for DM particles that are positive, indicating
    // that they are either inert stellar remnants or "leftovers" from
    // particle merging

    // check particle type and uniform mass
    rule0 = (ParticleType[i] == PARTICLE_TYPE_MUST_REFINE ||
                ParticleType[i] == PARTICLE_TYPE_MBH) ||
               ((ParticleMass[i] < UniformParticleMass) && (ParticleType[i] < PARTICLE_TYPE_INDIVIDUAL_STAR)) ||
      (ParticleMass[i] < UniformParticleMass &&
       ParticleType[i] == PARTICLE_TYPE_DARK_MATTER && OriginalParticle);

    rules[0] = rule0;

    // check particle mass greater than minimum mass
    rules[1] = (ParticleMass[i] > MustRefineParticlesMinimumMass);

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

    // printf("Checked if we should refine, the answer is %"ISYM"\n", IsParticleMustRefine[i]);
    // printf("rule 0 = %"ISYM" rule 1 = %"ISYM" rule 2 = %"ISYM"\n",rules[0],rules[1],rules[2]);

  PFORTRAN_NAME(cic_flag)(IsParticleMustRefine,
	   ParticlePosition[0], ParticlePosition[1], ParticlePosition[2],
	   &GridRank, &NumberOfParticles, FlaggingField,
	   LeftEdge, GridDimension, GridDimension+1, GridDimension+2,
	   &CellSize, &ParticleBufferSize);

  if (NumberOfAntiRules > 0){
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
//#endif // ifdef INDIVIDUALSTAR
