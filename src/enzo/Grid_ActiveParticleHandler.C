/***********************************************************************
/
/  GRID CLASS (HANDLE THE CREATION AND FEEDBACK OF ACTIVE PARTICLES)
/
/  written by: Greg Bryan
/  date:       March, 1997
/  modified1:  April, 2009 by JHW to have multiple types of star 
/              particles
/  modified2:  May, 2011 by MJT to be in support of active particles
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/
#include <assert.h> 
#include "preincludes.h"
#include "ErrorExceptions.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "CosmologyParameters.h"
#include "ActiveParticle.h"
#define H_DEBUG 0
/* function prototypes */
 
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int FindField(int field, int farray[], int numfields);
 
int grid::ActiveParticleHandler(HierarchyEntry* SubgridPointer, int level,
                                float dtLevelAbove, int &NumberOfNewParticles)
{

  if (EnabledActiveParticlesCount == 0) return SUCCESS;

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
 
  if (NumberOfBaryonFields == 0)
    return SUCCESS;

  /* First, set under_subgrid field */
  HierarchyEntry *Subgrid;
  this->ZeroSolutionUnderSubgrid(NULL, ZERO_UNDER_SUBGRID_FIELD);
  for (Subgrid = SubgridPointer; Subgrid; Subgrid = Subgrid->NextGridThisLevel)
    this->ZeroSolutionUnderSubgrid(Subgrid->GridData, ZERO_UNDER_SUBGRID_FIELD);

  /* initialize */

  LCAPERF_START("grid_ActiveParticleHandler");

  /* First we identify the data dependencies */

  struct ActiveParticleFormationDataFlags flags = flags_default;

  int i;
  for (i = 0; i < EnabledActiveParticlesCount; i++)
  {
    ActiveParticleType_info *ActiveParticleTypeToEvaluate = EnabledActiveParticles[i];
    ActiveParticleTypeToEvaluate->DescribeSupplementalData(flags);
  }

  struct ActiveParticleFormationData supplemental_data = data_default;
  supplemental_data.level = level;
  supplemental_data.GridID = this->ID;

  ActiveParticleType::ConstructData(this, flags, supplemental_data);

  /******************** FORMATION ********************/

  NumberOfNewParticles = 0;
  /* Now we iterate */
  for (i = 0; i < EnabledActiveParticlesCount; i++)
  {
    ActiveParticleType_info *ActiveParticleTypeToEvaluate = EnabledActiveParticles[i];
    ActiveParticleTypeToEvaluate->EvaluateFormation(
                                this, supplemental_data);
    NumberOfNewParticles += supplemental_data.NumberOfNewParticles;
    
  }

  /* Now we copy the particles from NewParticles into a statically allocated
   * array */

  if (NumberOfNewParticles > 0) {
    this->AddActiveParticles(supplemental_data.NewParticles, 0, NumberOfNewParticles);
    if (debug2)
      printf("Creating %d new active particles\n", NumberOfNewParticles);
  }

  /******************** FEEDBACK ********************/

  /* Now we iterate */
  if (NumberOfActiveParticles > 0)
    for (i = 0; i < EnabledActiveParticlesCount; i++)
      {
	ActiveParticleType_info *ActiveParticleTypeToEvaluate = EnabledActiveParticles[i];
	ActiveParticleTypeToEvaluate->EvaluateFeedback(this, supplemental_data);
      }
  
  ActiveParticleType::DestroyData(this, supplemental_data);

  //if (debug) printf("StarParticle: end\n");

  LCAPERF_STOP("grid_ActiveParticleHandler");
  return SUCCESS;
}


int grid::ActiveParticleHandler_Convert(HierarchyEntry* SubgridPointer, int level,
                                int gridnum, int &NumberOfNewParticles)
{
  if(NumberOfParticles == 0)
    return SUCCESS;
   if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS; 
  if (NumberOfBaryonFields == 0)
    return SUCCESS;
  int totalconverted = 0;  //static
  struct ActiveParticleFormationData supplemental_data = data_default;
  supplemental_data.level = level;
  supplemental_data.GridID = this->ID;
  supplemental_data.NumberOfNewParticles = 0;
   /* First, set under_subgrid field */
  HierarchyEntry *Subgrid;
  this->ZeroSolutionUnderSubgrid(NULL, ZERO_UNDER_SUBGRID_FIELD);
  for (Subgrid = SubgridPointer; Subgrid; Subgrid = Subgrid->NextGridThisLevel)
    this->ZeroSolutionUnderSubgrid(Subgrid->GridData, ZERO_UNDER_SUBGRID_FIELD);
 
  /* initialize */
  
  LCAPERF_START("grid_ActiveParticleHandler_Convert");
 
  float LifetimeFactor = 1.;
  int particles2convert = 0;

  int *hasharray = new int[NumberOfParticles];
  for (int i = 0; i < NumberOfParticles; i++)
    hasharray[i] = 1;
  NumberOfNewParticles = 0;
  for (int i = 0; i < NumberOfParticles; i++) {
    switch(ParticleType[i]) {
    case PARTICLE_TYPE_CLUSTER: 
      //      if (this->Time >= ParticleAttribute[0][i] &&
      //	  this->Time <= ParticleAttribute[0][i] + 
      //	  LifetimeFactor * ParticleAttribute[1][i]) {
	
	//printf("P%d: OK we are going to convert particle type %d (P2) on grid %d\n",
	//		 MyProcessorNumber, ParticleType[i], gridnum);
      {
      particles2convert++;
      char *dummy = new char[MAX_LINE_LENGTH];
      strcpy(dummy, "CenOstriker");
      ActiveParticleType_info *my_type =
        get_active_particle_types()[dummy];
      int enabled_num = my_type->GetEnabledParticleID();
      hasharray[i] = 0;
      /* OK we found a CenOstriker particles - make it into an active particle */
      ActiveParticleType_info *ActiveParticleTypeToEvaluate = EnabledActiveParticles[enabled_num];
      ActiveParticleTypeToEvaluate->CreateParticle(this, supplemental_data, i);
      }
      break;
	
    case PARTICLE_TYPE_SINGLE_STAR:
      //       if (this->Time >= ParticleAttribute[0][i] &&
      //	  this->Time <= ParticleAttribute[0][i] + 
      //	  LifetimeFactor * ParticleAttribute[1][i]) {
	 //printf("P%d: OK we are going to convert particle type %d (P3) on grid %d\n",
	 //     MyProcessorNumber, ParticleType[i], gridnum);
      {
      particles2convert++;
      char *dummy = new char[MAX_LINE_LENGTH];
      strcpy(dummy, "PopIII");
      ActiveParticleType_info *my_type =
        get_active_particle_types()[dummy];
      int enabled_num = my_type->GetEnabledParticleID();
      hasharray[i] = 0;
      /* OK we found a PopIII particles - make it into an active particle */
      ActiveParticleType_info *ActiveParticleTypeToEvaluate = EnabledActiveParticles[enabled_num];
      ActiveParticleTypeToEvaluate->CreateParticle(this, supplemental_data, i);
      }
      break;
    }
  }
  
  /* Now we copy the particles from NewParticles into a statically allocated
   * array */
 
  totalconverted += supplemental_data.NumberOfNewParticles;
  //if(totalconverted)
  // printf("%s: Adding %d new active particles to grid (total = %d)\n", __FUNCTION__,
  //	   supplemental_data.NumberOfNewParticles, totalconverted); fflush(stdout);
  NumberOfNewParticles = supplemental_data.NumberOfNewParticles;
  if (NumberOfNewParticles > 0) {
    this->AddActiveParticles(supplemental_data.NewParticles, 0, NumberOfNewParticles);
    //printf("Creating %d new active particles\n", NumberOfNewParticles);
  }  
#if H_DEBUG  
if(particles2convert) {
    printf("%s: P%d: particles2convert = %d\t NumberOfParticles = %d\n", __FUNCTION__, 
	   MyProcessorNumber, particles2convert, NumberOfParticles);
    printf("%s: NumberOfStars = %d\n", __FUNCTION__, ReturnNumberOfStars());
    printf("%s: NumberOfStarParticles = %d\n", __FUNCTION__, ReturnNumberOfStarParticles());
  } 
#endif
  if(particles2convert > 0) {
    int ii = 0;
    /* Destroy particles that are now active particles */
    for(int i = 0; i < NumberOfParticles; i++) {
      if(hasharray[i] == 1) {
	for(int j = 0; j < MAX_DIMENSION; j++) {
	  ParticlePosition[j][ii] = ParticlePosition[j][i];
	  ParticleVelocity[j][ii] = ParticleVelocity[j][i];
	  //ParticleAcceleration[j][ii] = ParticleAcceleration[j][i];
	}
	//if(ParticleAcceleration[i] != NULL)
	// ParticleAcceleration[MAX_DIMENSION][ii] = ParticleAcceleration[MAX_DIMENSION][i];
     
	ParticleMass[ii] = ParticleMass[i];
	ParticleNumber[ii] = ParticleNumber[i];
	ParticleType[ii] = ParticleType[i];
	for(int j = 0; j < NumberOfParticleAttributes;j++) 
	  ParticleAttribute[j][ii] = ParticleAttribute[j][i];
	ii++;
      }
    }
    SetNumberOfParticles(ii);
    NumberOfStarParticles = NumberOfStarParticles - particles2convert;
  }
 
  LCAPERF_STOP("grid_ActiveParticleHandler_Convert");
#if H_DEBUG
  if(particles2convert) {
    printf("P%d: NumberOfParticles now on grid = %d\t Number of Star Particles = %d\n", 
	   MyProcessorNumber, NumberOfParticles, NumberOfStarParticles);
  }
#endif
  return SUCCESS;
}
