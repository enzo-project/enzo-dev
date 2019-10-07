/***********************************************************************
/
/ Accreting Particle
/
************************************************************************/

#include "ActiveParticle_AccretingParticle.h"


/* We need to make sure that we can operate on the grid, so this dance is
 * necessary to make sure that grid is 'friend' to this particle type. */

class AccretingParticleGrid : private grid {
  friend class ActiveParticleType_AccretingParticle;
};

/* Note that we only refer to AccretingParticleGrid here.
 * Given a grid object, we static cast to get this:
 *
 *    AccretingParticleGrid *thisgrid =
 *      static_cast<AccretingParticleGrid *>(thisgrid_orig); */

float ActiveParticleType_AccretingParticle::OverflowFactor = FLOAT_UNDEFINED;
int ActiveParticleType_AccretingParticle::AccretionRadius = INT_UNDEFINED;
int ActiveParticleType_AccretingParticle::LinkingLength = INT_UNDEFINED;
int ActiveParticleType_AccretingParticle::RadiationParticle = INT_UNDEFINED;
double ActiveParticleType_AccretingParticle::LuminosityPerSolarMass = FLOAT_UNDEFINED;
int ActiveParticleType_AccretingParticle::RadiationSEDNumberOfBins = INT_UNDEFINED;
float* ActiveParticleType_AccretingParticle::RadiationEnergyBins = NULL;
float* ActiveParticleType_AccretingParticle::RadiationSED = NULL;
float ActiveParticleType_AccretingParticle::RadiationLifetime = FLOAT_UNDEFINED;


int ActiveParticleType_AccretingParticle::InitializeParticleType()
{
  // Leaving these defaults hardcoded for testing. NJG

  OverflowFactor = 1.01;
  LinkingLength = 4;
  AccretionRadius = 4;
  RadiationParticle = AccretingParticleRadiation;
  LuminosityPerSolarMass = AccretingParticleLuminosity;
  RadiationSEDNumberOfBins = 1;
  RadiationEnergyBins = new float[RadiationSEDNumberOfBins];
  RadiationSED = new float[RadiationSEDNumberOfBins];
  RadiationEnergyBins[0] = 21.0;
  RadiationSED[0] = 1.0;
  RadiationLifetime = 1000;


  /* Error check parameters */

  if (AccretingParticleRadiation) {
    if (RadiativeTransfer == FALSE)
    ENZO_FAIL("RadiativeTransfer must be turned with AccretingParticleRadiation.");
    if (MultiSpecies == FALSE)
      ENZO_FAIL("Multispecies must be turned with AccretingParticleRadiation.");
  }

  // Need to turn on particle mass flagging if it isn't already turned on.

  bool TurnOnParticleMassRefinement = true;
  int method;
  for (method = 0; method < MAX_FLAGGING_METHODS; method++)
    if (CellFlaggingMethod[method] == 8 || CellFlaggingMethod[method] == 4) {
      TurnOnParticleMassRefinement = false;
      break;
    }

  if (TurnOnParticleMassRefinement) {
    method = 0;
    while(CellFlaggingMethod[method] != INT_UNDEFINED)
      method++;
    CellFlaggingMethod[method] = 4;
  }

  /* Add on the Particle Array Handlers */
  typedef ActiveParticleType_AccretingParticle ap;
  AttributeVector &ah = ap::AttributeHandlers;
  ActiveParticleType::SetupBaseParticleAttributes(ah);

  ah.push_back(new Handler<ap, float, &ap::AccretionRate>("AccretionRate"));

  return SUCCESS;
}

int ActiveParticleType_AccretingParticle::EvaluateFormation
(grid *thisgrid_orig, ActiveParticleFormationData &data)
{
  // No need to do the rest if we're not on the maximum refinement level.
  if (data.level != MaximumRefinementLevel)
    return SUCCESS;

  AccretingParticleGrid *thisGrid =
    static_cast<AccretingParticleGrid *>(thisgrid_orig);

  int i,j,k,index,method,MassRefinementMethod;

  float *density = thisGrid->BaryonField[data.DensNum];

  float* velx = thisGrid->BaryonField[data.Vel1Num];
  float* vely = thisGrid->BaryonField[data.Vel2Num];
  float* velz = thisGrid->BaryonField[data.Vel3Num];

  float JeansDensityUnitConversion = (Gamma*pi*kboltz) / (Mu*mh*GravConst);
  float CellTemperature = 0;
  float JeansDensity = 0;
  float MassRefinementDensity = 0;
  float DensityThreshold = huge_number;
  float ExtraDensity = 0;

  int GridDimension[3] = {thisGrid->GridDimension[0],
                          thisGrid->GridDimension[1],
                          thisGrid->GridDimension[2]};

  FLOAT dx = thisGrid->CellWidth[0][0];

  bool HasMetalField = (data.MetalNum != -1 || data.ColourNum != -1);
  bool JeansRefinement = false;
  bool MassRefinement = false;

  const int offset[] = {1, GridDimension[0], GridDimension[0]*GridDimension[1]};

  // determine refinement criteria
  for (method = 0; method < MAX_FLAGGING_METHODS; method++) {
    if (CellFlaggingMethod[method] == 2) {
      MassRefinement = true;
      MassRefinementMethod = method;
    }
    if (CellFlaggingMethod[method] == 6)
      JeansRefinement = true;
  }

  for (k = thisGrid->GridStartIndex[2]; k <= thisGrid->GridEndIndex[2]; k++) {
    for (j = thisGrid->GridStartIndex[1]; j <= thisGrid->GridEndIndex[1]; j++) {
      for (i = thisGrid->GridStartIndex[0]; i <= thisGrid->GridEndIndex[0]; i++) {
	index = GRIDINDEX_NOGHOST(i, j, k);

	// If no more room for particles, throw an ENZO_FAIL
	if (data.NumberOfNewParticles >= data.MaxNumberOfNewParticles)
	  return FAIL;

	// Does this cell violate the Jeans condition?
	DensityThreshold = huge_number;
	if (JeansRefinement) {
	  CellTemperature = (JeansRefinementColdTemperature > 0) ? JeansRefinementColdTemperature : data.Temperature[index];
	  JeansDensity = JeansDensityUnitConversion * OverflowFactor * CellTemperature /
	    POW(data.LengthUnits*dx*4.0,2);
	  JeansDensity /= data.DensityUnits;
	  DensityThreshold = min(DensityThreshold,JeansDensity);
	}
	if (DensityThreshold == huge_number)
	  ENZO_VFAIL("Error in Accreting Particles: DensityThreshold = huge_number! \n"
		     "JeansDensity = %"GOUTSYM" \n"
		     "data.DensityUnits = %"GOUTSYM" \n"
		     "CellTemperature = %"GOUTSYM" \n",
		     JeansDensity, data.DensityUnits, CellTemperature);

	if (density[index] <= DensityThreshold)
	  continue;

	// Passed creation tests, create sink particle

	ActiveParticleType_AccretingParticle *np = new ActiveParticleType_AccretingParticle();
    data.NumberOfNewParticles++;
    data.NewParticles.insert(*np);

	ExtraDensity = density[index] - DensityThreshold;
	np->Mass = ExtraDensity;   // Particle 'masses' are actually densities
	np->type = np->GetEnabledParticleID();
	np->BirthTime = thisGrid->ReturnTime();

	np->level = data.level;
	np->GridID = data.GridID;
	np->CurrentGrid = thisGrid;

	np->pos[0] = thisGrid->CellLeftEdge[0][i] + 0.5*thisGrid->CellWidth[0][i];
	np->pos[1] = thisGrid->CellLeftEdge[1][j] + 0.5*thisGrid->CellWidth[1][j];
	np->pos[2] = thisGrid->CellLeftEdge[2][k] + 0.5*thisGrid->CellWidth[2][k];

	if (HydroMethod == PPM_DirectEuler) {
	  np->vel[0] = velx[index];
	  np->vel[1] = vely[index];
	  np->vel[2] = velz[index];
	}
	else if (HydroMethod == Zeus_Hydro) {
	  np->vel[0] = 0.5 * (velx[index] + velx[index+offset[0]]);
	  np->vel[1] = 0.5 * (vely[index] + vely[index+offset[1]]);
	  np->vel[2] = 0.5 * (velz[index] + velz[index+offset[2]]);
	} else {
	  ENZO_FAIL("AccretingParticle does not support RK Hydro or RK MHD");
	}

	if (HasMetalField)
	  np->Metallicity = data.TotalMetals[index];
	else
	  np->Metallicity = 0.0;

	np->AccretionRate = 0.0;

	// Remove mass from grid

	density[index] = DensityThreshold;

      } // i
    } // j
  } // k

  return SUCCESS;
}

int ActiveParticleType_AccretingParticle::EvaluateFeedback(grid *thisgrid_orig, ActiveParticleFormationData &data)
{
  AccretingParticleGrid *thisGrid =
    static_cast<AccretingParticleGrid *>(thisgrid_orig);

  int npart = thisGrid->NumberOfActiveParticles;

  for (int n = 0; n < npart; n++) {
    ActiveParticleType_AccretingParticle *ThisParticle =
      static_cast<ActiveParticleType_AccretingParticle*>(
              thisGrid->ActiveParticles[n]);

    ThisParticle->level = data.level;
  }

  return SUCCESS;
}

void ActiveParticleType_AccretingParticle::DescribeSupplementalData(ActiveParticleFormationDataFlags &flags)
{
  flags.DarkMatterDensity = true;
  flags.Temperature = true;
  flags.MetalField = true;
}


template <class active_particle_class>
int ActiveParticleType_AccretingParticle::BeforeEvolveLevel
(HierarchyEntry *Grids[], TopGridData *MetaData,
 int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
 int ThisLevel, bool CallEvolvePhotons, 
 int TotalStarParticleCountPrevious[],
 int AccretingParticleID)
{

  FLOAT Time = LevelArray[ThisLevel]->GridData->ReturnTime();
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
    VelocityUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, Time);

  int j, dim, ipart, nParticles;
  ActiveParticleList<ActiveParticleType> AccretingParticleList;
  if (CallEvolvePhotons)
    ActiveParticleFindAll(LevelArray, &nParticles, AccretingParticleID, 
        AccretingParticleList);

  /* Create radiation sources from particles */
  
  // Calculate conversion factor to solar masses

#ifdef TRANSFER
  if (CallEvolvePhotons) {
    RadiationSourceEntry* source;
    double dx;
    float MassConversion;
    const double LConv = (double) TimeUnits / pow(LengthUnits,3);
    const double mfactor = double(DensityUnits) / SolarMass;

    ActiveParticleType_AccretingParticle *ThisParticle = NULL;
    for (ipart = 0; ipart < nParticles; ipart++) {
      if (AccretingParticleList[ipart]->IsARadiationSource(Time)) {
	ThisParticle =
	  static_cast<ActiveParticleType_AccretingParticle*>(
          AccretingParticleList[ipart]);
	source = ThisParticle->RadiationSourceInitialize();
	dx = LengthUnits * 
      LevelArray[ThisParticle->level]->GridData->GetCellWidth(0,0);
	MassConversion = (float) (dx*dx*dx * mfactor);
	source->LifeTime       = RadiationLifetime;
	source->Luminosity = (LuminosityPerSolarMass * LConv) *
	  (ThisParticle->Mass * MassConversion); 
	source->EnergyBins = RadiationSEDNumberOfBins;
	source->Energy = new float[RadiationSEDNumberOfBins];
	source->SED = new float[RadiationSEDNumberOfBins];
	for (j = 0; j < RadiationSEDNumberOfBins; j++) {
	  source->Energy[j] = RadiationEnergyBins[j];
	  source->SED[j] = RadiationSED[j];
	}
//	printf("SRC Luminosity: L=%lg Lcode=%g M=%g Mcode=%g\n",
//	       LuminosityPerSolarMass * ThisParticle->Mass * MassConversion, 
//	       source->Luminosity, ThisParticle->Mass * MassConversion, 
//	       ThisParticle->Mass);
	if (GlobalRadiationSources->NextSource != NULL)
	  GlobalRadiationSources->NextSource->PreviousSource = source;
	GlobalRadiationSources->NextSource = source;
      }
    }
  } // ENDIF CallEvolvePhotons
#endif /* TRANSFER */
  
  return SUCCESS;
}


grid* ConstructFeedbackZone(ActiveParticleType* ThisParticle, int FeedbackRadius,
			     FLOAT dx, HierarchyEntry** Grids, int NumberOfGrids,
			     int SendField);

int DistributeFeedbackZone(grid* FeedbackZone, HierarchyEntry** Grids,
			   int NumberOfGrids, int SendField);

int ActiveParticleType_AccretingParticle::Accrete(int nParticles, 
    ActiveParticleList<ActiveParticleType>& ParticleList,
    int AccretionRadius, FLOAT dx,
    LevelHierarchyEntry *LevelArray[], int ThisLevel)
{

  /* Skip accretion if we're not on the maximum refinement level.
     This should only ever happen right after creation and then
     only in pathological cases where sink creation is happening at
     the edges of two regions at the maximum refinement level */

  if (ThisLevel < MaximumRefinementLevel)
    return SUCCESS;

  /* For each particle, loop over all of the grids and do accretion
     if the grid overlaps with the accretion zone                   */

  int i, NumberOfGrids;
  int *FeedbackRadius = NULL;
  HierarchyEntry **Grids = NULL;
  grid *sinkGrid = NULL;

  bool SinkIsOnThisProc, SinkIsOnThisGrid;

  float SubtractedMass, SubtractedMomentum[3] = {};

  NumberOfGrids = GenerateGridArray(LevelArray, ThisLevel, &Grids);

  for (i = 0; i < nParticles; i++) {
    grid* FeedbackZone = ConstructFeedbackZone(ParticleList[i], AccretionRadius,
					       dx, Grids, NumberOfGrids, ALL_FIELDS);

    if (MyProcessorNumber == FeedbackZone->ReturnProcessorNumber()) {

      float AccretionRate = 0;

      if (FeedbackZone->AccreteOntoAccretingParticle(ParticleList[i],
			      AccretionRadius*dx, &AccretionRate) == FAIL)
	return FAIL;

      // No need to communicate the accretion rate to the other CPUs since this particle is already local.
      static_cast<ActiveParticleType_AccretingParticle*>(ParticleList[i])->AccretionRate = AccretionRate;
    }

    DistributeFeedbackZone(FeedbackZone, Grids, NumberOfGrids, ALL_FIELDS);

    delete FeedbackZone;
  }

  if (AssignActiveParticlesToGrids(ParticleList, nParticles, LevelArray) == FAIL)
    return FAIL;

  delete [] Grids;
  return SUCCESS;
}

int ActiveParticleType_AccretingParticle::SetFlaggingField(
    LevelHierarchyEntry *LevelArray[], int level,
    int TopGridDims[], int AccretingParticleID)
{
  /* Generate a list of all sink particles in the simulation box */
  int i, nParticles;
  FLOAT *pos = NULL, dx=0;
  ActiveParticleList<ActiveParticleType> AccretingParticleList;
  LevelHierarchyEntry *Temp = NULL;

  ActiveParticleFindAll(LevelArray, &nParticles, AccretingParticleID, 
      AccretingParticleList);

  /* Calculate CellWidth on maximum refinement level */

  dx = (DomainRightEdge[0] - DomainLeftEdge[0]) /
    (TopGridDims[0]*POW(FLOAT(RefineBy),FLOAT(MaximumRefinementLevel)));

  for (i=0 ; i<nParticles; i++){
    pos = AccretingParticleList[i]->ReturnPosition();
    for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel)
      if (Temp->GridData->DepositRefinementZone(level,pos,AccretionRadius*dx) == FAIL) {
	ENZO_FAIL("Error in grid->DepositRefinementZone.\n")
	  }
  }

  return SUCCESS;
}

/* 
 * This function can be used to reset the particle acceleration if required.
 * For example if a massless particle needs to be fixed in space. 
 * See ActiveParticle_RadiationParticle.C for details. 
 */
int ActiveParticleType_AccretingParticle::ResetAcceleration(float *ActiveParticleAcceleration)
{
  return SUCCESS;
}

/*
 * For brute force creation of a particle. Useful for converting from 
 * star objects to active particles.
 */
int ActiveParticleType_AccretingParticle::CreateParticle(grid *thisgrid_orig,
							 ActiveParticleFormationData &supp_data,
							 int particle_index)
{
  return SUCCESS;
} 

bool ActiveParticleType_AccretingParticle::IsARadiationSource(FLOAT Time)
{
  return (RadiationParticle == TRUE) ? true : false;
}

namespace {
  ActiveParticleType_info *AccretingParticleInfo =
    register_ptype <ActiveParticleType_AccretingParticle>
    ("AccretingParticle");
}

std::vector<ParticleAttributeHandler*>
  ActiveParticleType_AccretingParticle::AttributeHandlers;
