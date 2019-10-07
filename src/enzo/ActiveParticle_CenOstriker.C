/***********************************************************************
/
/  Cen & Ostriker star formation
/
************************************************************************/

#include "ActiveParticle_CenOstriker.h"
#include "phys_constants.h"


/* We need to make sure that we can operate on the grid, so this dance is
 * necessary to make sure that grid is 'friend' to this particle type. */

class CenOstrikerGrid : private grid {
  friend class ActiveParticleType_CenOstriker;
};

/* Defaults for parameters */

float ActiveParticleType_CenOstriker::OverdensityThreshold = FLOAT_UNDEFINED;
float ActiveParticleType_CenOstriker::MassEfficiency = FLOAT_UNDEFINED;
float ActiveParticleType_CenOstriker::MinimumDynamicalTime = FLOAT_UNDEFINED;
float ActiveParticleType_CenOstriker::MinimumStarMass = FLOAT_UNDEFINED;
float ActiveParticleType_CenOstriker::MassEjectionFraction = FLOAT_UNDEFINED;
float ActiveParticleType_CenOstriker::EnergyToThermalFeedback = FLOAT_UNDEFINED;
float ActiveParticleType_CenOstriker::MetalYield = FLOAT_UNDEFINED;

int ActiveParticleType_CenOstriker::FeedbackDistTotalCells = INT_UNDEFINED;
int ActiveParticleType_CenOstriker::FeedbackDistRadius = INT_UNDEFINED;
int ActiveParticleType_CenOstriker::FeedbackDistCellStep = INT_UNDEFINED;

bool ActiveParticleType_CenOstriker::JeansMassCriterion = true;
bool ActiveParticleType_CenOstriker::StochasticStarFormation = false;
bool ActiveParticleType_CenOstriker::UnigridVelocities = false;
bool ActiveParticleType_CenOstriker::PhysicalOverdensity = false;
bool ActiveParticleType_CenOstriker::dtDependence = true;

int ActiveParticleType_CenOstriker::InitializeParticleType() {
  
  OverdensityThreshold = StarMakerOverDensityThreshold;
  MassEfficiency = StarMakerMassEfficiency;
  MinimumDynamicalTime = StarMakerMinimumDynamicalTime;
  MinimumStarMass = StarMakerMinimumMass;
  MassEjectionFraction = StarMassEjectionFraction;
  EnergyToThermalFeedback = StarEnergyToThermalFeedback;
  MetalYield = StarMetalYield;
  FeedbackDistRadius = StarFeedbackDistRadius;
  FeedbackDistCellStep = StarFeedbackDistCellStep;
  JeansMassCriterion = true;
  StochasticStarFormation = false;
  UnigridVelocities = false;
  PhysicalOverdensity = true;


  ActiveParticleType::SetupBaseParticleAttributes(
    ActiveParticleType_CenOstriker::AttributeHandlers);

  return SUCCESS;
}

int ActiveParticleType_CenOstriker::CreateParticle(grid *thisgrid_orig,
						   ActiveParticleFormationData &data,
						   int particle_index)
{
  //printf("Brute force particle creation of CenOstriker\n");
 
  CenOstrikerGrid *thisGrid =
    static_cast<CenOstrikerGrid *>(thisgrid_orig);

  ActiveParticleType_CenOstriker *np = new ActiveParticleType_CenOstriker();
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
  
  //  if (debug && data.NumberOfNewParticles > 0)
  // fprintf(stderr, "AP_CenOstriker: Have created %"ISYM" new particles\n",
  //	    data.NumberOfNewParticles);

  return SUCCESS;
}
int ActiveParticleType_CenOstriker::EvaluateFormation(grid *thisgrid_orig, ActiveParticleFormationData &data)
{
  CenOstrikerGrid *thisGrid =
    static_cast<CenOstrikerGrid *>(thisgrid_orig);
  
  float BaryonMass,VelocityDivergence,TotalDensity,DynamicalTime,
    IsothermalSoundSpeedSquared,JeansMass,StarFraction, RandomNumber;
  float SoundSpeedConstant = 1.3095e8;
  int i, j, k, dim, index, offset_y, offset_z;
  int NumberOfNewParticles = 0;


  float *density = thisGrid->BaryonField[data.DensNum];
  float *velx = thisGrid->BaryonField[data.Vel1Num];
  float *vely = thisGrid->BaryonField[data.Vel2Num];
  float *velz = thisGrid->BaryonField[data.Vel3Num];

  FLOAT dx = data.LengthUnits * thisGrid->CellWidth[0][0];

  bool HasMetalField = (data.MetalNum != -1 || data.ColourNum != -1);

  int GridDimension[3] = {thisGrid->GridDimension[0],
                          thisGrid->GridDimension[1],
                          thisGrid->GridDimension[2]};

  // Pre-calculate serialized offsets for the 3D data field.  Used for
  // the divergence.
  offset_y = thisGrid->GridDimension[0];
  offset_z = thisGrid->GridDimension[0] * thisGrid->GridDimension[1];

  for (k = thisGrid->GridStartIndex[2]; k <= thisGrid->GridEndIndex[2]; k++) {
    for (j = thisGrid->GridStartIndex[1]; j <= thisGrid->GridEndIndex[1]; j++) {
      index = GRIDINDEX_NOGHOST(thisGrid->GridStartIndex[0], j, k);
      for (i = thisGrid->GridStartIndex[0]; i <= thisGrid->GridEndIndex[0]; i++, index++) {

	// 0. If no more room for particles, quit.
	if (data.NumberOfNewParticles >=
	    data.MaxNumberOfNewParticles)
          continue;
	
	// 1. Finest level of refinement
	if (thisGrid->BaryonField[thisGrid->NumberOfBaryonFields][index] != 0.0) 
	  continue;
	
	// 2. Density greater than threshold
	if (density[index] < OverdensityThreshold)
	  continue;

	/* 3. Negative divergence: For ZEUS, the velocities are
	   face-centered, and all of the other routines have
	   cell-centered velocities. */
	
	if (HydroMethod == Zeus_Hydro) {
	  VelocityDivergence = velx[index+1] - velx[index] +
	    vely[index+offset_y] - vely[index] + 
	    velz[index+offset_z] - velz[index];
	} else {
	  VelocityDivergence = velx[index+1] - velx[index-1] + 
	    vely[index+offset_y] - vely[index-offset_y] + 
	    velz[index+offset_z] - velz[index-offset_z];
	}

	if (VelocityDivergence > 0.0) continue;

	// 4. t_cool < t_freefall (skip if T < 11000 K)
	TotalDensity = ( density[index] + data.DarkMatterDensity[index] ) * 
	  data.DensityUnits;
	DynamicalTime = sqrt(3.0 * M_PI / 32.0 / GravConst / TotalDensity) / data.TimeUnits;

	if (DynamicalTime < data.CoolingTime[index] && 
	    data.Temperature[index] > 1.1e4)
	  continue;

	// 5. Cell mass is greater than the Jeans Mass
        BaryonMass = density[index] * data.DensityUnits * 
          POW(dx, 3) / SolarMass;
	if (JeansMassCriterion) {
	  IsothermalSoundSpeedSquared = SoundSpeedConstant * data.Temperature[index];
	  JeansMass = M_PI / (6.0 * sqrt(density[index] * data.DensityUnits)) *
	    POW(M_PI * IsothermalSoundSpeedSquared / GravConst, 1.5) / SolarMass;

	  if (BaryonMass < JeansMass)
	    continue;
	}

	// 6) Check to see if star is above threshold (given in units of M_solar)
	StarFraction = min(MassEfficiency *
                           thisGrid->ReturnTimeStep() / DynamicalTime, 0.9);
	DynamicalTime = max(DynamicalTime, MinimumDynamicalTime * yr_s /
                            data.TimeUnits);
	
	// 7) If we allow stochastic star formation, make new particles 
        //    every time the unfulfilled star formation buffer
	//    exceeds the mininimum particle mass
	if (StochasticStarFormation) {
	  if (StarFraction*BaryonMass < MinimumStarMass) {
	    UnfulfilledStarFormationMass += StarFraction*BaryonMass;
	    if (UnfulfilledStarFormationMass < MinimumStarMass) 
	      continue;
	    StarFraction = min(MinimumStarMass/BaryonMass, 0.5);
	    UnfulfilledStarFormationMass -= StarFraction*BaryonMass;
	  } 
	}

        // If not doing stochastic star formation, check if star 
        // particle mass exceeds the minimum
        else if (StarFraction * BaryonMass < MinimumStarMass) {
          continue;
        }

	/*
	 * ====================================================================
	 * PARTICLE CREATION
	 * ====================================================================
	 */

	ActiveParticleType_CenOstriker *np = new ActiveParticleType_CenOstriker();
    data.NumberOfNewParticles++;
    data.NewParticles.insert(*np);

	np->level = data.level;
	np->GridID = data.GridID;
	np->CurrentGrid = thisgrid_orig;

	np->Mass = StarFraction*density[index];
	np->type = np->GetEnabledParticleID();
	np->BirthTime = thisGrid->ReturnTime();
	np->DynamicalTime = DynamicalTime;

	np->pos[0] = thisGrid->CellLeftEdge[0][i] + 0.5*thisGrid->CellWidth[0][i];
	np->pos[1] = thisGrid->CellLeftEdge[1][j] + 0.5*thisGrid->CellWidth[1][j];
	np->pos[2] = thisGrid->CellLeftEdge[2][k] + 0.5*thisGrid->CellWidth[2][k];
	
	if (UnigridVelocities == false) {
	  if (HydroMethod != Zeus_Hydro) {
	    np->vel[0] = velx[index];
	    np->vel[1] = vely[index];
	    np->vel[2] = velz[index];
	  }
	  else {
	    np->vel[0] = 0.5*(velx[index]+velx[index+1]);
	    np->vel[1] = 0.5*(vely[index]+vely[index+offset_y]);
	    np->vel[2] = 0.5*(velz[index]+velz[index+offset_z]);
	  }
	} 
	else {
	  np->vel[0] = tiny_number;
	  np->vel[1] = tiny_number;
	  np->vel[2] = tiny_number;	  
	}

	if (HasMetalField)
	  np->Metallicity = data.TotalMetals[index];
	else
	  np->Metallicity = 0.0;

	// Remove mass from grid

	density[index] = (1.0 - StarFraction)*density[index];

      }
    }
  }

  if (debug && data.NumberOfNewParticles > 0)
    fprintf(stderr, "AP_CenOstriker: Have created %"ISYM" new particles\n",
	    data.NumberOfNewParticles);

  return SUCCESS;
}

int ActiveParticleType_CenOstriker::EvaluateFeedback
(grid *thisGrid_orig, ActiveParticleFormationData &data)
{
  CenOstrikerGrid *thisGrid =
    static_cast<CenOstrikerGrid *>(thisGrid_orig);
  
  float *density = thisGrid->BaryonField[data.DensNum];
  float *velx = thisGrid->BaryonField[data.Vel1Num];
  float *vely = thisGrid->BaryonField[data.Vel2Num];
  float *velz = thisGrid->BaryonField[data.Vel3Num];
  float *totalenergy = thisGrid->BaryonField[data.TENum];
  float *gasenergy = thisGrid->BaryonField[data.GENum];
  float *metals = thisGrid->BaryonField[data.MetalNum];
  float dt = thisGrid->dtFixed;
  float dx = float(thisGrid->CellWidth[0][0]);

  float xv1, xv2, ParticleBirthTime, ParticleDynamicalTimeAtBirth,
    ParticleMass, ParticleInitialMass, ParticleMetalFraction, 
    StarFormationDensityThisTimestep, SupernovaEnergyThisTimestep, 
    DensityToAddToEachCell, DensityRatio;

  float StellarMassFormedThisTimestepOnThisGrid = 0;

  FLOAT xpos, ypos, zpos;
  float xvel, yvel, zvel;

  FLOAT CurrentTime = thisGrid->Time;
  FLOAT xstart = thisGrid->CellLeftEdge[0][0];
  FLOAT ystart = thisGrid->CellLeftEdge[1][0];
  FLOAT zstart = thisGrid->CellLeftEdge[2][0];

  int npart = thisGrid->NumberOfActiveParticles;
  int GridXSize = thisGrid->GridDimension[0];
  int GridYSize = thisGrid->GridDimension[1];
  int GridZSize = thisGrid->GridDimension[2];
  int NumberOfGhostZones = thisGrid->GridStartIndex[0];
  int GridDimension[3] = {thisGrid->GridDimension[0],
			  thisGrid->GridDimension[1],
			  thisGrid->GridDimension[2]};
  
  int n,i,j,k,ic,kc,jc,stepk,stepj,cellstep,DistIndex,index;

  for (n=0; n < npart; n++) {
    ActiveParticleType_CenOstriker *particle = 
      static_cast<ActiveParticleType_CenOstriker*>(thisGrid->ActiveParticles[n]);

    xpos = particle->pos[0];
    ypos = particle->pos[1];
    zpos = particle->pos[2];
  
    xvel = particle->vel[0];
    yvel = particle->vel[1];
    zvel = particle->vel[2];

    ParticleBirthTime = particle->BirthTime;
    ParticleDynamicalTimeAtBirth = particle->DynamicalTime;
    ParticleMass = particle->Mass;
    ParticleMetalFraction = particle->Metallicity;
    
    // Determine how much of a given star particle would have been
    // turned into stars during this timestep.  Then, calculate the
    // mass which should have formed during this timestep dt using the
    // integral form of the Cen & Ostriker formula.
    
    xv1 = (CurrentTime - ParticleBirthTime) / ParticleDynamicalTimeAtBirth;
    if (xv1 > 12.0) continue; // current time - creation time >> dynamical time at formation, so ignore
    
    xv2 = (CurrentTime + dt - ParticleBirthTime) / ParticleDynamicalTimeAtBirth;

    // First, calculate the initial mass of the star particle in question
    
    ParticleInitialMass = ParticleMass / 
      (1.0 - MassEjectionFraction * (1.0 - (1.0 + xv1)*exp(-xv1)));
    
    // Then, calculate the amount of mass that would have formed in this timestep.

    StarFormationDensityThisTimestep = ParticleInitialMass * ((1.0 + xv1)*exp(-xv1) - 
                                                              (1.0 + xv2)*exp(-xv2));
    
    StarFormationDensityThisTimestep = max(min(StarFormationDensityThisTimestep,ParticleMass),0.0);
      
    // Calculate 3D grid indices

    i = int((xpos - xstart)/thisGrid->CellWidth[0][0]);
    j = int((ypos - ystart)/thisGrid->CellWidth[1][0]);
    k = int((zpos - zstart)/thisGrid->CellWidth[2][0]);

    // Check bounds - if star particle is outside of this grid then give a warning and continue
    
    if (i < 0 || i > GridXSize-1 || j < 0 || j > GridYSize-1 || k < 0 || k > GridZSize-1){
      fprintf(stdout, "Particle out of grid; xind, yind, zind = %"ISYM", %"ISYM", %"ISYM"\n",i,j,k);
      continue;
    }
      
    // Calculate serial index

    index = GRIDINDEX_NOGHOST(i,j,k);

    // skip if very little mass is formed

    if (StarFormationDensityThisTimestep/density[index] < 1e-10 )
      continue;

    // calculate mass added to each cell

    DensityToAddToEachCell = (StarFormationDensityThisTimestep * 
                              MassEjectionFraction) / StarFeedbackDistTotalCells;

    // If using distributed feedback, check if particle is too close 
    // to the boundary and adjust indices accordingly

    if (FeedbackDistRadius > 0)
      {
	i = max(NumberOfGhostZones + FeedbackDistRadius,
		min(GridXSize - NumberOfGhostZones - FeedbackDistRadius - 1, i));
	j = max(NumberOfGhostZones + FeedbackDistRadius,
		min(GridYSize - NumberOfGhostZones - FeedbackDistRadius - 1, j));
	k = max(NumberOfGhostZones + FeedbackDistRadius,
		min(GridZSize - NumberOfGhostZones - FeedbackDistRadius - 1, k));	
      }

    // Subtract ejected mass from particle
    
    ParticleMass -= StarFormationDensityThisTimestep *
      MassEjectionFraction;

    // Save particle mass

    particle->Mass = ParticleMass;

    // Record amount of star formation in this grid

    StellarMassFormedThisTimestepOnThisGrid += 
      StarFormationDensityThisTimestep * dt * POW(dx,3);

    // Calculate supernova energy for this timestep

    SupernovaEnergyThisTimestep = EnergyToThermalFeedback * 
      StarFormationDensityThisTimestep * POW(clight/data.VelocityUnits,2) /
      StarFeedbackDistTotalCells;

#define NO_SHARE_ENERGY
#ifdef SHARE_ENERGY
    SupernovaEnergyThisTimestep *= StarFormationDensityThisTimestep *
      MassEjectionFraction / 
      (StarFormationDensityThisTimestep * MassEjectionFraction + 
       ParticleInitialMass*exp(-xv2)*(1.0+xv2));
#endif /* SHARE_ENERGY */

    // Add energy to the energy field
    for (kc = k - FeedbackDistRadius; kc <= k + FeedbackDistRadius; kc++){
      stepk = ABS(kc - k);
      for (jc = j - FeedbackDistRadius; jc <= j + FeedbackDistRadius; jc++){
	stepj = stepk + ABS(jc - j);
	for (ic = i - FeedbackDistRadius; ic <= i + FeedbackDistRadius; ic++){
	  cellstep = stepj + ABS(ic - i);
	  DistIndex = GRIDINDEX_NOGHOST(ic,jc,kc);
	  if (cellstep <= FeedbackDistCellStep) {
	    DensityRatio = 1.0/(density[DistIndex] + DensityToAddToEachCell);
	    totalenergy[DistIndex] = ((totalenergy[DistIndex]*density[DistIndex]) + 
				      SupernovaEnergyThisTimestep)*DensityRatio;
	    if (DualEnergyFormalism == 1)
	      gasenergy[DistIndex] = ((gasenergy[DistIndex]*density[DistIndex]) + 
				      SupernovaEnergyThisTimestep)*DensityRatio;

	    // Metal feedback (note that in this function gas metal is a fraction
	    // (rho_metal/rho_gas) rather than a density.  The conversion has 
	    // been done in the handling routine)

	    // The "Cen Method".  This takes into account gas recycling:

	    if (data.MetalNum != -1) {
	      metals[DistIndex] = (metals[DistIndex] * density[DistIndex] +
		(StarFormationDensityThisTimestep / StarFeedbackDistTotalCells) *
		(MetalYield * (1.0 - ParticleMetalFraction) +
		 MassEjectionFraction * ParticleMetalFraction)) * DensityRatio;
            }

	    // Mass and momentum feedback

	    velx[DistIndex] = velx[DistIndex] * density[DistIndex] +
              DensityToAddToEachCell * xvel;
	    vely[DistIndex] = vely[DistIndex] * density[DistIndex] +
              DensityToAddToEachCell * yvel;
	    velz[DistIndex] = velz[DistIndex] * density[DistIndex] +
              DensityToAddToEachCell * zvel;
	    density[DistIndex] += DensityToAddToEachCell;
	    velx[DistIndex] /= density[DistIndex];
	    vely[DistIndex] /= density[DistIndex];
	    velz[DistIndex] /= density[DistIndex];

            // If using DualEnergyFormalism, use gas energy to 
            // set total energy
            if ((HydroMethod != Zeus_Hydro) && 
                (DualEnergyFormalism == 1)) {
              totalenergy[DistIndex] = gasenergy[DistIndex] +
                0.5 * (velx[DistIndex] * velx[DistIndex] +
                       vely[DistIndex] * vely[DistIndex] +
                       velz[DistIndex] * velz[DistIndex]);
            }

	  }
	}
      }
    }
    
  } // end loop over particles
  
  return SUCCESS;
}

void ActiveParticleType_CenOstriker::DescribeSupplementalData(ActiveParticleFormationDataFlags &flags)
{
  flags.DarkMatterDensity = true;
  flags.CoolingTime = true;
  flags.Temperature = true;
  flags.MetalField = true;
}


/* 
 * This function can be used to reset the particle acceleration if required.
 * For example if a massless particle needs to be fixed in space. 
 * See ActiveParticle_RadiationParticle.C for details. 
 */
int ActiveParticleType_CenOstriker::ResetAcceleration(float *ActiveParticleAcceleration)
{
  return SUCCESS;
}

int ActiveParticleType_CenOstriker::SetFlaggingField(LevelHierarchyEntry *LevelArray[],int level, int TopGridDims[], int ActiveParticleID)
{

  return SUCCESS;
}

namespace {
  ActiveParticleType_info *CenOstrikerInfo = 
    register_ptype <ActiveParticleType_CenOstriker> ("CenOstriker");
}

std::vector<ParticleAttributeHandler*>
  ActiveParticleType_CenOstriker::AttributeHandlers;

