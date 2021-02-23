/***********************************************************************
/
/  GRID CLASS (INDIVIDUAL STAR DEPOSIT MUST-REFINE PARTICLES)
/
/  written by: Andrew Emerick
/  date:       2018
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
#include "phys_constants.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);

/* function prototypes */

extern "C" void PFORTRAN_NAME(cic_flag)(int* irefflag, FLOAT *posx, FLOAT *posy,
			FLOAT *posz, int *ndim, int *npositions,
                        int *ffield, FLOAT *leftedge,
                        int *dim1, int *dim2, int *dim3, FLOAT *cellsize,
			int *buffersize);

int IndividualStarInterpolateLifetime(float &tau, const float &M,
                                      const float &metallicity, const int &mode);

//
// Separate method to keep things clean, even though a lot will be redundant
// in Grid_DepositMustRefineParticles
//
int grid::IndividualStarDepositMustRefineParticles(int pmethod, int level, bool KeepFlaggingField,
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

  /* This is a bit of a hack --- want to make sure injection regions
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

  /* Count the number of stars possible */
  int num_stars = 0;
  Star *cstar;
  for (cstar = AllStars; cstar; cstar = cstar->NextStar){
    num_stars++;
  }
  int refine_size = num_stars + num_events;
  IsParticleMustRefine = new int[refine_size];
  StarPosX = new FLOAT[refine_size];
  StarPosY = new FLOAT[refine_size];
  StarPosZ = new FLOAT[refine_size];


  for (i = 0; i < refine_size; i++){
    IsParticleMustRefine[i] = 0;
    StarPosX[i] = -1.;
    StarPosY[i] = -1.;
    StarPosZ[i] = -1.;
  }

  int NumberOfMustRefineStars = 0;
  FLOAT *pos;

  int count = 0;
  for (cstar = AllStars, i = 0; cstar; cstar = cstar->NextStar) {
    float end_of_life = 0.0;
    float lifetime = 0.0;

    /* Remnants behave differently */
    if ( !(cstar->ReturnType() == PARTICLE_TYPE_INDIVIDUAL_STAR_REMNANT)){
        //cstar->ReturnType() == PARTICLE_TYPE_INDIVIDUAL_STAR ||
        //cstar->ReturnType() == PARTICLE_TYPE_INDIVIDUAL_STAR_POPIII ||
        //cstar->ReturnType() == PARTICLE_TYPE_INDIVIDUAL_STAR_WD){

        lifetime = cstar->ReturnLifetime();

    } else{

      if (ProblemType == 260 && ChemicalEvolutionTestStarLifetime > 0){
        lifetime = ChemicalEvolutionTestStarLifetime * Myr_s / TimeUnits;
      } else {
        /* Else, we need to re-compute the main sequence (pre-remnant)
           lifetime */

        lifetime = 0.0;
        int mode = 1;
        if (IndividualStarPopIIIFormation && (cstar->IsPopIII())) mode = 3;

        IndividualStarInterpolateLifetime(lifetime, cstar->ReturnBirthMass(),
                                          cstar->ReturnMetallicity(), mode);

        lifetime /= TimeUnits;
      }
    }

    end_of_life = cstar->ReturnBirthTime() + lifetime;

    bool near_end_of_life = fabs(this->Time - end_of_life) < IndividualStarRefineTime * Myr_s / TimeUnits;
   //(IndividualStarLifeRefinementFactor * this->dtFixed * POW(2,level)); // factor of root grid, estimate root $

    IsParticleMustRefine[i] = 0;
    StarPosX[i] = StarPosY[i] = StarPosZ[i] = -1.;

    if (cstar->ReturnType() == PARTICLE_TYPE_INDIVIDUAL_STAR){
        if ((( ( IndividualStarStellarWinds) && (cstar->ReturnMass() > IndividualStarSNIIMassCutoff)  ) || // massive stars always on if winds are on
                   ( (!IndividualStarStellarWinds) && (near_end_of_life)  ) ||  // SNII check if no winds are on
                   ( ( IndividualStarStellarWinds) && (near_end_of_life) && (cstar->ReturnMass() < IndividualStarSNIIMassCutoff) )) // AGB wind check
             || (IndividualStarRefineForRadiation && (cstar->ReturnBirthMass()>=IndividualStarIonizingRadiationMinimumMass))){ // radiation check

          IsParticleMustRefine[i] = 1;
        }
    } else if (fabs(cstar->ReturnType()) == PARTICLE_TYPE_INDIVIDUAL_STAR_POPIII){
      if  ((near_end_of_life   &&
            (   ((cstar->ReturnBirthMass()>=TypeIILowerMass)&&(cstar->ReturnBirthMass()<=TypeIIUpperMass)) ||
                ((cstar->ReturnBirthMass()>=PISNLowerMass)  &&(cstar->ReturnBirthMass()<=PISNUpperMass  ))
            )) ||
            (IndividualStarRefineForRadiation && (cstar->ReturnBirthMass()>=IndividualStarIonizingRadiationMinimumMass))){

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


  int refine_buffer_size = IndividualStarRefineBufferSize;

  if (IndividualStarRefineToPhysicalRadius > 0){
    refine_buffer_size = (int)(ceil((2 * IndividualStarRefineToPhysicalRadius * pc_cm / LengthUnits) / CellSize));

    refine_buffer_size = max(refine_buffer_size, 1); // always
  }

  PFORTRAN_NAME(cic_flag)(IsParticleMustRefine,
	                  StarPosX, StarPosY, StarPosZ,
      	                  &GridRank, &(refine_size), FlaggingField,
	                  LeftEdge, GridDimension, GridDimension+1, GridDimension+2,
                          &CellSize, &refine_buffer_size);

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
