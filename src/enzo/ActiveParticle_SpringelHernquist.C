/***********************************************************************
/
/  SPRINGEL & HERNQUIST STAR FORMATION
/
/  written by: Stephen Skory
/  date:       October, 2011
/
/  PURPOSE:
/
************************************************************************/

#include "ActiveParticle_SpringelHernquist.h"

/* We need to make sure that we can operate on the grid, so this dance is
 * necessary to make sure that grid is 'friend' to this particle type. */

class SpringelHernquistGrid : private grid {
  friend class ActiveParticleType_SpringelHernquist;
};

/* Note that we only refer to SpringelHernquistGrid here. 
 * Given a grid object, we static case to get this:
 *
 *    SpringelHernquistGrid *thisgrid =
 *      static_cast<SpringelHernquistGrid *>(thisgrid_orig); */

float ActiveParticleType_SpringelHernquist::OverDensityThreshold = FLOAT_UNDEFINED;
float ActiveParticleType_SpringelHernquist::PhysicalDensityThreshold = FLOAT_UNDEFINED;
float ActiveParticleType_SpringelHernquist::MinimumDynamicalTime = FLOAT_UNDEFINED;
float ActiveParticleType_SpringelHernquist::MinimumMass = FLOAT_UNDEFINED;

int ActiveParticleType_SpringelHernquist::InitializeParticleType() {
  // get some parameters from the Param object

  OverDensityThreshold = StarMakerOverDensityThreshold;
  PhysicalDensityThreshold = StarMakerSHDensityThreshold;
  MinimumDynamicalTime = StarMakerMinimumDynamicalTime;
  MinimumMass = StarMakerMinimumMass;

  /* Add on the Particle Array Handlers */
  typedef ActiveParticleType_SpringelHernquist ap;
  AttributeVector &ah = ap::AttributeHandlers;
  ActiveParticleType::SetupBaseParticleAttributes(ah);

  /* We don't want to change the attribute this is tied with, but we do want to
  update the name. */
  
  for(AttributeVector::iterator it = ah.begin(); it != ah.end(); ++it) {
    if((*it)->name == "metallicity") (*it)->name = "metallicity_fraction";
  }
  
  return SUCCESS;
  
}

int ActiveParticleType_SpringelHernquist::EvaluateFormation
(grid *thisgrid_orig, ActiveParticleFormationData &supp_data)
{
  SpringelHernquistGrid *tg =
    static_cast<SpringelHernquistGrid *>(thisgrid_orig);

  float bmass, tstar;
  int i, j, k, index, offset_y, offset_z;
  int NumberOfNewParticles = 0;
  float r_float, y, x, pstar, starfraction, usn;
  float msolar=SolarMass, mproton=mh, beta=0.1,
    sqrtepssn=2.0e24, kb=kboltz;

  // Define the energy output from supernovae

  usn = (1. - beta) * sqrtepssn / beta / msolar; // Below Eq (2)
  usn = usn * sqrtepssn; // this trick done to avoid floating-point issues (BWO)

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

	// Calculate star formation timescale. Eq 21.
	tstar =  yr_s * MinimumDynamicalTime * 
	  pow(density[index] * supp_data.DensityUnits / PhysicalDensityThreshold, -0.5);

	/* note: minus sign is because 'coolrate' is actual backwards: when coolrate is 
	 * negative, gas is cooling (which is the opposite of how it's defined in Springel
	 * & Hernquist 2003, MNRAS, 339, 289, eqtn. 16 */
	y = -1.0 * tstar * supp_data.CoolingRate[index] /
	  (density[index] * supp_data.DensityUnits) /
	  (beta * usn - (1.0 - beta) * 1.5 *
	   kb * supp_data.Temperature[index] / 0.6 / mproton);

	// Calculate the fraction of mass in cold clouds.
	if (y <= 0) {
	  x = 0;
	} else {
	  x = 1. + 1./ 2. / y - sqrt(1. / y + 1. / 4. / y / y); // Eq(18)
	}

	// Calculate total baryon mass in the cell.
	bmass = density[index] * supp_data.MassUnits;

	// Calculate a parameter which is compared to a random number.
	pstar = (bmass / MinimumMass) * 
	  (1. - exp(-(1. - beta) * x * tg->dtFixed * 
		    supp_data.TimeUnits / tstar)); // Eq(39)

	// Make the random number.
	r_float = float(rand())/float(RAND_MAX);

	// Finally, if this random number is smaller than pstar, make a star.
	if (r_float > pstar) continue;
	
	/*
	 * ====================================================================
	 * PARTICLE CREATION
	 * ====================================================================
	 */

	starfraction = min(MinimumMass / bmass, 0.5);

	ActiveParticleType_SpringelHernquist *np = 
      new ActiveParticleType_SpringelHernquist();
    supp_data.NumberOfNewParticles++;
    supp_data.NewParticles.insert(*np);
	//fprintf(stderr, "G_APH: Creating !\n");

	// Give the star mass and then remove the same from the grid.
	np->Mass = starfraction * density[index];
	density[index] -= starfraction * density[index];

	np->type = PARTICLE_TYPE_STAR;
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

// SH star feedback
int ActiveParticleType_SpringelHernquist::EvaluateFeedback
(grid *thisgrid_orig, ActiveParticleFormationData &data)
{
  return SUCCESS;
}

void ActiveParticleType_SpringelHernquist::DescribeSupplementalData
(ActiveParticleFormationDataFlags &flags)
{
  flags.CoolingTime = true;
  flags.Temperature = true;
  flags.MetalField = true;
}

int ActiveParticleType_SpringelHernquist::SetFlaggingField(LevelHierarchyEntry *LevelArray[], int level,
							int TopGridDims[], int SpringelHernquistID)
{
  return SUCCESS;
}

/* 
 * This function can be used to reset the particle acceleration if required.
 * For example if a massless particle needs to be fixed in space. 
 * See ActiveParticle_RadiationParticle.C for details. 
 */
int ActiveParticleType_SpringelHernquist::ResetAcceleration(float *ActiveParticleAcceleration)
{
  return SUCCESS;
}

/*
 * For brute force creation of a particle. Useful for converting from 
 * star objects to active particles.
 */
int ActiveParticleType_SpringelHernquist::CreateParticle(grid *thisgrid_orig,
							 ActiveParticleFormationData &supp_data,
							 int particle_index)
{
  return SUCCESS;
} 

namespace {
  ActiveParticleType_info *SpringelHernquistInfo = 
    register_ptype <ActiveParticleType_SpringelHernquist> ("SpringelHernquist");
}

std::vector<ParticleAttributeHandler*>
  ActiveParticleType_SpringelHernquist::AttributeHandlers;

