/***********************************************************************
/
/  KRAVTSOV (2003) STAR FORMATION
/
************************************************************************/

#include "ActiveParticle_Kravtsov.h"

/* We need to make sure that we can operate on the grid, so this dance is
 * necessary to make sure that grid is 'friend' to this particle type. */

class KravtsovGrid : private grid {
    friend class ActiveParticleType_Kravtsov;
};

/* Note that we only refer to SampleParticleGrid here. 
 * Given a grid object, we static case to get this:
 *
 *    SampleParticleGrid *thisgrid =
 *      static_cast<SampleParticleGrid *>(thisgrid_orig); */

float ActiveParticleType_Kravtsov::DensityThreshold = FLOAT_UNDEFINED;
float ActiveParticleType_Kravtsov::StarFormationTimeConstant = FLOAT_UNDEFINED;
float ActiveParticleType_Kravtsov::MinimumStarMass = FLOAT_UNDEFINED;

int ActiveParticleType_Kravtsov::InitializeParticleType(void) {

  DensityThreshold = StarMakerOverDensityThreshold;
  StarFormationTimeConstant = StarMakerMinimumDynamicalTime;
  MinimumStarMass = StarMakerMinimumMass;
  
  
  /* Add on the Particle Array Handlers */
  typedef ActiveParticleType_Kravtsov ap;
  AttributeVector &ah = ap::AttributeHandlers;
  ActiveParticleType::SetupBaseParticleAttributes(ah);

  /* We don't want to change the attribute this is tied with, but we do want to
  update the name. */
  
  for(AttributeVector::iterator it = ah.begin(); it != ah.end(); ++it) {
    if((*it)->name == "metallicity") (*it)->name = "metallicity_fraction";
  }

  return SUCCESS;
}


int ActiveParticleType_Kravtsov::EvaluateFormation
(grid *thisgrid_orig, ActiveParticleFormationData &supp_data)
{
  KravtsovGrid *tg =
    static_cast<KravtsovGrid *>(thisgrid_orig);


  /* Make it pretty */

  float *density = tg->BaryonField[supp_data.DensNum];
//  float *velx = tg->BaryonField[supp_data.Vel1Num];
//  float *vely = tg->BaryonField[supp_data.Vel2Num];
//  float *velz = tg->BaryonField[supp_data.Vel3Num];

  float CellWidthTemp = float(tg->CellWidth[0][0]);

  bool HasMetalField = (supp_data.MetalNum != -1 ||
			supp_data.ColourNum != -1);

  int GridDimension[3] = {tg->GridDimension[0],
                          tg->GridDimension[1],
                          tg->GridDimension[2]};

  float gasfrac, starmass, densthresh, timeconstant;
  int i,j,k,index;
  int NumberOfNewParticles = 0;


  // calculate density threshold.  odthresh is in proper particles per
  // cc and (d1/mproton) gives the mean density of the universe in
  // particles/cm^3 (assuming hydrogen is dominant)

  densthresh = DensityThreshold / (supp_data.DensityUnits / mh);


  // calculate time constant for star formation.  This assumes that
  // the user input is in units of years
  
  timeconstant = StarFormationTimeConstant * 3.156e7 / supp_data.TimeUnits;


  // for each zone, : "star" particle is created if the density
  // exceeds some threshold and this is the highest level of
  // refinement.  That's it.

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
	if (density[index] < densthresh)
	  continue;
	
	/*
	 * ====================================================================
	 * PARTICLE CREATION
	 * ====================================================================
	 */
	
	ActiveParticleType_Kravtsov *np = new ActiveParticleType_Kravtsov();
	supp_data.NewParticles.insert(*np);
    supp_data.NumberOfNewParticles++;

	// Make sure that we never give put than 90% of the cell's mass into a star particle
	gasfrac = min( 0.9, tg->dtFixed / timeconstant );
	
	// Calculate star mass in solar masses.  If this is less than
	// the user-defined threshold mass, do NOT make a star in this
	// cell.  This is not exactly in keeping with the spirit of
	// the Kravtsov algorithm, and is somewhat degenerate with the
	// density threshold, but we really don't want millions and
	// millions of star particles.

	starmass = gasfrac * density[index] * supp_data.DensityUnits * pow(supp_data.LengthUnits * CellWidthTemp, 3) /  SolarMass;

	// Do not allow stars with mass less than MinimumStarMass
	if (starmass < MinimumStarMass )
	  continue;

	np->Mass =  starmass / supp_data.MassUnits;

	np->type = Kravtsov;
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

// Feedback
int ActiveParticleType_Kravtsov::EvaluateFeedback(grid *thisgrid_orig, ActiveParticleFormationData &data)
{
  return SUCCESS;
}

void ActiveParticleType_Kravtsov::DescribeSupplementalData
(ActiveParticleFormationDataFlags &flags)
{
  flags.MetalField = true;
}

int ActiveParticleType_Kravtsov::SetFlaggingField(LevelHierarchyEntry *LevelArray[], int level,
							int TopGridDims[], int KravtsovID)
{
  return SUCCESS;
}


/* 
 * This function can be used to reset the particle acceleration if required.
 * For example if a massless particle needs to be fixed in space. 
 * See ActiveParticle_RadiationParticle.C for details. 
 */
int ActiveParticleType_Kravtsov::ResetAcceleration(float *ActiveParticleAcceleration)
{
  return SUCCESS;
}

/*
 * For brute force creation of a particle. Useful for converting from 
 * star objects to active particles.
 */
int ActiveParticleType_Kravtsov::CreateParticle(grid *thisgrid_orig,
						ActiveParticleFormationData &supp_data,
						int particle_index)
{
  return SUCCESS;
} 

namespace {
  ActiveParticleType_info *KravtsovInfo = 
    register_ptype <ActiveParticleType_Kravtsov> ("Kravtsov");
}

std::vector<ParticleAttributeHandler*>
  ActiveParticleType_Kravtsov::AttributeHandlers;

