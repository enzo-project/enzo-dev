
/***********************************************************************
/
/  AN EXAMPLE ACTIVE PARTICLE TYPE
/
/  written by: Matthew Turk
/  date:       May, 2011
/
/  PURPOSE:
/
************************************************************************/

#include "ActiveParticle_PopIII.h"

float CalculatePopIIILifetime(float Mass);

/* We need to make sure that we can operate on the grid, so this dance is
 * necessary to make sure that grid is 'friend' to this particle type. */

class PopIIIGrid : private grid {
  friend class ActiveParticleType_PopIII;
};

float ActiveParticleType_PopIII::OverDensityThreshold = FLOAT_UNDEFINED;
float ActiveParticleType_PopIII::MetalCriticalFraction = FLOAT_UNDEFINED;
float ActiveParticleType_PopIII::H2CriticalFraction = FLOAT_UNDEFINED;
float ActiveParticleType_PopIII::StarMass = FLOAT_UNDEFINED;

// get some parameters from the Param object
int ActiveParticleType_PopIII::InitializeParticleType() {

  OverDensityThreshold = PopIIIOverDensityThreshold;
  MetalCriticalFraction = PopIIIMetalCriticalFraction;
  H2CriticalFraction = PopIIIH2CriticalFraction;
  StarMass = PopIIIStarMass;


  ActiveParticleType::SetupBaseParticleAttributes(
    ActiveParticleType_PopIII::AttributeHandlers);

  return SUCCESS;
}

int ActiveParticleType_PopIII::EvaluateFormation
(grid *thisgrid_orig, ActiveParticleFormationData &supp_data)
{
  PopIIIGrid *tg = static_cast<PopIIIGrid *>(thisgrid_orig);

  float bmass, div, dtot, tdyn, LifetimeInYears;
  int i, j, k, dim, index, offset_y, offset_z;
  int NumberOfNewParticles = 0;

  /* Make it pretty */

  float *density = tg->BaryonField[supp_data.DensNum];
  float *velx = tg->BaryonField[supp_data.Vel1Num];
  float *vely = tg->BaryonField[supp_data.Vel2Num];
  float *velz = tg->BaryonField[supp_data.Vel3Num];

  bool HasMetalField = (supp_data.MetalNum != -1 ||
			supp_data.ColourNum != -1);

  int GridDimension[3] = {tg->GridDimension[0],
                          tg->GridDimension[1],
                          tg->GridDimension[2]};

  // Pre-calculate serialized offsets for the 3D data field.  Used for
  // the divergence.
  offset_y = tg->GridDimension[0];
  offset_z = tg->GridDimension[0] * tg->GridDimension[1];

  for (k = tg->GridStartIndex[2]; k <= tg->GridEndIndex[2]; k++) {
    for (j = tg->GridStartIndex[1]; j <= tg->GridEndIndex[1]; j++) {
      index = GRIDINDEX_NOGHOST(tg->GridStartIndex[0], j, k);
      for (i = tg->GridStartIndex[0]; i <= tg->GridEndIndex[0]; i++, index++) {
	
	// 0. If no more room for particles, quit.
	if (supp_data.NumberOfNewParticles >=
	    supp_data.MaxNumberOfNewParticles)
          continue;
	
	// 1. Finest level of refinement
	if (tg->BaryonField[tg->NumberOfBaryonFields][index] != 0.0) 
	  continue;
	
	// 2. Density greater than threshold
	if (density[index] < OverDensityThreshold)
	  continue;

	/* 3. Negative divergence: For ZEUS, the velocities are
	   face-centered, and all of the other routines have
	   cell-centered velocities. */

	if (HydroMethod == Zeus_Hydro) {
	  div = velx[index+1] - velx[index] +
	    vely[index+offset_y] - vely[index] + 
	    velz[index+offset_z] - velz[index];
	} else {
	  div = velx[index+1] - velx[index-1] + 
	    vely[index+offset_y] - vely[index-offset_y] + 
	    velz[index+offset_z] - velz[index-offset_z];
	}

	if (div > 0.0) continue;

	// 4. t_cool < t_freefall (skip if T > 11000 K)
	dtot = ( density[index] + supp_data.DarkMatterDensity[index] ) * 
	  supp_data.DensityUnits;
	tdyn = sqrt(3.0 * M_PI / 32.0 / GravConst / dtot) / supp_data.TimeUnits;
	if (tdyn < supp_data.CoolingTime[index] && 
	    supp_data.Temperature[index] > 1.1e4)
	  continue;

	// 5. If metallicity is greater than the critical metallicity
	if (HasMetalField && 
	    supp_data.TotalMetals[index] > MetalCriticalFraction)
	  continue;

	// 6. Require some H2 fraction to form a metal-free star
	if (supp_data.H2Fraction[index] < H2CriticalFraction)
	  continue;

	/*
	 * ====================================================================
	 * PARTICLE CREATION
	 * ====================================================================
	 */

	/* Compute the lifetime of a Pop III star, given its mass.
	   Note for an IMF, we have to recompute this after the random
	   sampling. */

	ActiveParticleType_PopIII *np = new ActiveParticleType_PopIII();
    supp_data.NumberOfNewParticles++;
    supp_data.NewParticles.insert(*np);
	//fprintf(stderr, "G_APH: Creating !\n");

	LifetimeInYears = CalculatePopIIILifetime(StarMass);

	// Mass of the star will be assigned by the accretion routines.
	if (RadiativeTransfer) {
	  np->Mass = 0.0;
	  np->Lifetime = LifetimeInYears * yr_s / supp_data.TimeUnits;
	} else {
	  bmass = density[index] * supp_data.MassUnits;
	  np->Mass = min(0.5 * bmass, StarMass / supp_data.MassUnits);
	  np->Lifetime = LifetimeInYears * yr_s / supp_data.TimeUnits;
	}

	np->type = PopIII;
	np->BirthTime = tg->Time;
	
	np->pos[0] = tg->CellLeftEdge[0][i] + 0.5*tg->CellWidth[0][i];
	np->pos[1] = tg->CellLeftEdge[1][j] + 0.5*tg->CellWidth[1][j];
	np->pos[2] = tg->CellLeftEdge[2][k] + 0.5*tg->CellWidth[2][k];

	/*
	  Star velocities averaged over multiple cells to avoid
	  "runaway star particle" phenomenon imethod = 2 is zeus,
	  otherwise PPM
	*/


	float *tvel = tg->AveragedVelocityAtCell(index, supp_data.DensNum,
						 supp_data.Vel1Num);
	np->vel[0] = tvel[0];
	np->vel[1] = tvel[1];
	np->vel[2] = tvel[2];

	/* Set the metallicity */

	if (HasMetalField)
	  np->Metallicity = supp_data.TotalMetals[index];
	else
	  np->Metallicity = 0.0;

      } // ENDFOR i
    } // ENDFOR j
  } // ENDFOR k

  return NumberOfNewParticles;
}


// Pop III feedback (done through rad. transfer and feedback spheres)
int ActiveParticleType_PopIII::EvaluateFeedback
(grid *thisgrid_orig, ActiveParticleFormationData &supp_data)
{
  return SUCCESS;
}

void ActiveParticleType_PopIII::DescribeSupplementalData
(ActiveParticleFormationDataFlags &flags)
{
  flags.DarkMatterDensity = true;
  flags.H2Fraction = true;
  flags.CoolingTime = true;
  flags.Temperature = true;
  flags.MetalField = true;
}

int ActiveParticleType_PopIII::SetFlaggingField(LevelHierarchyEntry *LevelArray[], int level, 
						int TopGridDims[], int PopIIIParticleID)
{
  return SUCCESS;
}

bool ActiveParticleType_PopIII::IsARadiationSource(FLOAT Time)
{
  return true;
}

/* 
 * This function can be used to reset the particle acceleration if required.
 * For example if a massless particle needs to be fixed in space. 
 * See ActiveParticle_RadiationParticle.C for details. 
 */
int ActiveParticleType_PopIII::ResetAcceleration(float *ActiveParticleAcceleration)
{
  return SUCCESS;
}

int ActiveParticleType_PopIII::CreateParticle(grid *thisgrid_orig,
					      ActiveParticleFormationData &data,
					      int particle_index)
{
  //printf("Brute force particle creation of PopIII");
  PopIIIGrid *thisGrid =
    static_cast<PopIIIGrid *>(thisgrid_orig);

  ActiveParticleType_PopIII *np = new ActiveParticleType_PopIII();
  data.NumberOfNewParticles++;
  data.NewParticles.insert(*np);

  np->level = data.level;
  np->GridID = data.GridID;
  np->CurrentGrid = thisgrid_orig;
  
  np->Mass = thisgrid_orig->ParticleMass[particle_index];
  np->type = np->GetEnabledParticleID();
  np->BirthTime = thisgrid_orig->ParticleAttribute[0][particle_index];
  np->DynamicalTime = thisgrid_orig->ParticleAttribute[1][particle_index];
  np->Metallicity = thisgrid_orig->ParticleAttribute[2][particle_index];
  np->pos[0] = thisgrid_orig->ParticlePosition[0][particle_index];
  np->pos[1] = thisgrid_orig->ParticlePosition[1][particle_index];
  np->pos[2] = thisgrid_orig->ParticlePosition[2][particle_index];
  np->vel[0] = thisgrid_orig->ParticleVelocity[0][particle_index];
  np->vel[1] = thisgrid_orig->ParticleVelocity[1][particle_index];
  np->vel[2] = thisgrid_orig->ParticleVelocity[2][particle_index];
  if (debug && data.NumberOfNewParticles > 0)
    fprintf(stderr, "AP_PopIII: Have created %"ISYM" new particles\n",
	    data.NumberOfNewParticles);

  return SUCCESS;
}
namespace {
  ActiveParticleType_info *PopIIIInfo = 
    register_ptype <ActiveParticleType_PopIII> ("PopIII");
}

std::vector<ParticleAttributeHandler*>
  ActiveParticleType_PopIII::AttributeHandlers;

