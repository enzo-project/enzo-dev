/***********************************************************************
/
/ Smart Star Particle
/ Author: John Regan
/         Based on previous active particle code from many authors
/ 								Edited 2022 by Simone Gordon
/ Date: Early 2017
/
************************************************************************/

#include "ActiveParticle_SmartStar.h"
#include "phys_constants.h"
#define SSDEBUG 0
#define SSDEBUG_TOTALMASS 1

#define DYNAMIC_ACCRETION_RADIUS 0
#define BONDIHOYLERADIUS 1 // SG. Turning on accretion out to BH radius in Accrete.
#define MINIMUMPOTENTIAL 1
#define CALCDIRECTPOTENTIAL 0
#define JEANSREFINEMENT  0 // SG. turning off to check potential fix.
#define MASSTHRESHOLDCHECK 1  //SG. Turned on for testing. Turning off again.
#define JEANSLENGTHCALC    1
#define MASSTHRESHOLD      0.1 //Msolar in grid. SG. changed to 20 to prevent runaway SF in EvaluateFormation.
#define COOLING_TIME       0 // SG. Turn on to prevent spurious SF.Turning off again.
#define NUMSSPARTICLETYPES 4
#define JEANS_FACTOR       2
#define STELLAR_ACCRETION_OFF 1 // SG. Turns off accretion for POPIII if =1.
#define HW_BH_MASS 1   // SG. BH forms with mass according to Heger-Woosley 2002 relation.
#define SNEFEEDBACK 1
#define ONE_PARTICLE_ONLY 0 // SG. Disabling further particle formation. Must be set to 1 after forming 1 particle.
int DetermineSEDParameters(ActiveParticleType_SmartStar *SS,FLOAT Time, FLOAT dx);
float CalculatePopIIILifetime(float Mass);

/* We need to make sure that we can operate on the grid, so this dance is
 * necessary to make sure that grid is 'friend' to this particle type. */

class SmartStarGrid : private grid {
  friend class ActiveParticleType_SmartStar;
};

/* Note that we only refer to SmartStarGrid here.
 * Given a grid object, we static cast to get this:
 *
 *    SmartStarGrid *thisgrid =
 *      static_cast<SmartStarGrid *>(thisgrid_orig); */
float ActiveParticleType_SmartStar::EjectedMassThreshold = FLOAT_UNDEFINED;
double ActiveParticleType_SmartStar::LuminosityPerSolarMass = FLOAT_UNDEFINED;
int ActiveParticleType_SmartStar::RadiationSEDNumberOfBins = INT_UNDEFINED;
float* ActiveParticleType_SmartStar::RadiationEnergyBins = NULL;
float* ActiveParticleType_SmartStar::RadiationSED = NULL;
int ActiveParticleType_SmartStar::FeedbackDistRadius = INT_UNDEFINED;
int ActiveParticleType_SmartStar::FeedbackDistTotalCells = INT_UNDEFINED;
int ActiveParticleType_SmartStar::FeedbackDistCellStep = INT_UNDEFINED;
static double JeansLength(float T, float dens, float density_units);
static void UpdateAccretionRadius(ActiveParticleType*  ThisParticle, float newmass,
				  FLOAT AccretionRadius, float avgtemp,
				  float mass_units, float length_units); // SG. Took out dx argument.
static float GetStellarRadius(float cmass, float accrate);
int SmartStarPopIII_IMFInitialize(void);
int ActiveParticleType_SmartStar::InitializeParticleType()
{

  EjectedMassThreshold = 1.0;
  RadiationSEDNumberOfBins = NUMRADIATIONBINS;
  RadiationEnergyBins = new float[RadiationSEDNumberOfBins];
  RadiationSED = new float[RadiationSEDNumberOfBins];
  FeedbackDistRadius = StarFeedbackDistRadius;
  FeedbackDistCellStep = StarFeedbackDistCellStep;

#if SSDEBUG
  printf("%s: Initialising SmartStar\n", __FUNCTION__);
#endif

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
  typedef ActiveParticleType_SmartStar ap;
  AttributeVector &ah = ap::AttributeHandlers;
  ActiveParticleType::SetupBaseParticleAttributes(ah);

  ah.push_back(new Handler<ap, FLOAT, &ap::AccretionRadius>("AccretionRadius"));
  ah.push_back(new Handler<ap, int, &ap::ParticleClass>("ParticleClass"));
  ah.push_back(new Handler<ap, double, &ap::RadiationLifetime>("RadiationLifetime"));
  ah.push_back(new Handler<ap, double, &ap::StellarAge>("StellarAge"));
  ah.push_back(new Handler<ap, float, &ap::NotEjectedMass>("NotEjectedMass"));
  ah.push_back(new Handler<ap, float, &ap::MassToBeEjected>("MassToBeEjected"));
  ah.push_back(new Handler<ap, float, &ap::eta_disk>("eta_disk"));
  ah.push_back(new Handler<ap, float, &ap::beta_jet>("beta_jet"));
  ah.push_back(new Handler<ap, float, &ap::InfluenceRadius>("InfluenceRadius"));
  ah.push_back(new Handler<ap, float, &ap::epsilon_deltat>("epsilon_deltat"));
  ah.push_back(new Handler<ap, float, &ap::mass_in_accretion_sphere>("mass_in_accretion_sphere"));
  ah.push_back(new ArrayHandler<ap, float, MAX_DIMENSION, &ap::Accreted_angmom>("Accreted_angmom", 0));
  ah.push_back(new ArrayHandler<ap, float, NTIMES, &ap::AccretionRate>("AccretionRate", 0));
  ah.push_back(new ArrayHandler<ap, float, NTIMES, &ap::AccretionRateTime>("AccretionRateTime", 0));
  ah.push_back(new Handler<ap, int, &ap::TimeIndex>("TimeIndex"));
  
  /* Initialize the IMF lookup table if requested and not defined */
  /* Generate a list of all sink particles in the simulation box */

  if (PopIIIInitialMassFunction)
    SmartStarPopIII_IMFInitialize();
  
  printf("SmartStar Initialisation complete\n"); fflush(stdout);
  return SUCCESS;
}

int ActiveParticleType_SmartStar::EvaluateFormation
(grid *thisgrid_orig, ActiveParticleFormationData &data)
{
    // SG. If we just want one SmartStar in simulation.
    // SG. Temporarily disable star formation other than the current particle.
    if (ONE_PARTICLE_ONLY) {
        fprintf(stderr, "Total number of SmartStar particles is 1. Disabling further formation %s\n",
                __FUNCTION__);
        return SUCCESS;
    }
    // No need to do the rest if we're not on the maximum refinement level/ local max level of refinement
    SmartStarGrid *thisGrid =
    static_cast<SmartStarGrid *>(thisgrid_orig);

    //static int shu_collapse = 0;
    int i,j,k,index,method,MassRefinementMethod;

    float *density = thisGrid->BaryonField[data.DensNum];

    float* velx = thisGrid->BaryonField[data.Vel1Num];
    float* vely = thisGrid->BaryonField[data.Vel2Num];
    float* velz = thisGrid->BaryonField[data.Vel3Num];
    float div = 0.0, divx = 0.0, divy = 0.0, divz = 0.0, dtot = 0.0, tdyn = 0.0;
    float mass = 0.0;
    float JeansDensityUnitConversion = (Gamma*pi*kboltz) / (Mu*mh*GravConst);
    float CellTemperature = 0;
    float JeansDensity = 0, JeansMass = 0;
    float MassRefinementDensity = 0;
    float DensityThreshold = ActiveParticleDensityThreshold*mh/data.DensityUnits;   //in code density
    float ExtraDensity = 0;
    float GravitationalMinimum = 0.0;
    float ThermalEnergy = 0.0, GravitationalEnergy = 0.0, KineticEnergy = 0.0, CalcTotalEnergy = 0.0,
    TotalMass = 0.0;
    int GridDimension[3] = {thisGrid->GridDimension[0],
                          thisGrid->GridDimension[1],
                          thisGrid->GridDimension[2]};
    int size = GridDimension[0]*GridDimension[1]*GridDimension[2];
    float ConverttoSolar = data.DensityUnits*POW(data.LengthUnits, 3.0)/SolarMass;
    FLOAT dx = thisGrid->CellWidth[0][0], centralpos[3];

    bool HasMetalField = (data.MetalNum != -1 || data.ColourNum != -1);
    bool JeansRefinement = false;
    bool MassRefinement = false;
    float *PotentialField  = NULL;
    if(data.GravPotentialNum >= 0)
    PotentialField = thisGrid->BaryonField[data.GravPotentialNum];
    FLOAT dx_pc = dx*data.LengthUnits/pc_cm;   //in pc

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

#if MASSTHRESHOLDCHECK
    JeansMass = thisGrid->CalculateJeansMass(data.DensNum, data.Temperature, data.DensityUnits);  //In Msolar
#endif

// SG. Check we're on the maximum LOCAL refinement level from the get-go. 
    for (k = thisGrid->GridStartIndex[2]; k <= thisGrid->GridEndIndex[2]; k++) {
        for (j = thisGrid->GridStartIndex[1]; j <= thisGrid->GridEndIndex[1]; j++) {
            index = GRIDINDEX_NOGHOST(thisGrid->GridStartIndex[0], j, k);
            for (i = thisGrid->GridStartIndex[0]; i <= thisGrid->GridEndIndex[0]; i++, index++) {
                if (thisGrid->BaryonField[thisGrid->NumberOfBaryonFields][index] != 0.0){
                    continue;
                } else{

	DensityThreshold = ActiveParticleDensityThreshold*mh/data.DensityUnits;
	// If no more room for particles, throw an ENZO_FAIL
	if (data.NumberOfNewParticles >= data.MaxNumberOfNewParticles)
	  return FAIL;
	// 1. Check we're on the maximum refinement level - already done at start of function
	// SG. Check we're on the maximum LOCAL refinement level from the get-go.
	if (thisGrid->BaryonField[thisGrid->NumberOfBaryonFields][index] != 0.0)
	 continue;

	// 2. Does this cell violate the Jeans condition or overdensity threshold
	// Fedderath condition #1
#if JEANSREFINEMENT
	if (JeansRefinement) {
	  CellTemperature = (JeansRefinementColdTemperature > 0) ? JeansRefinementColdTemperature : data.Temperature[index];
	  /*
	   * The density threshold is exceeded here once the Jeans
	   * Density in one cell is exceeded - see Krumholz et al. (2004)
	   * Krumholz et al. use a Jeans factor of 4, I found that
	   * pushing this to 2 (i.e. making the threshold higher) is
	   * more accurate.
	   */
	  int JeansFactor = JEANS_FACTOR;
	  JeansDensity = JeansDensityUnitConversion * 1.01 * CellTemperature /
	    POW(data.LengthUnits*dx*JeansFactor,2);

	  JeansDensity /= data.DensityUnits;
	  DensityThreshold = min(DensityThreshold,JeansDensity);
	}
	if (DensityThreshold == huge_number)
	  ENZO_VFAIL("Error in Accreting Particles: DensityThreshold = huge_number! \n"
		     "JeansDensity = %"GOUTSYM" \n"
		     "data.DensityUnits = %"GOUTSYM" \n"
		     "CellTemperature = %"GOUTSYM" \n",
		     JeansDensity, data.DensityUnits, CellTemperature);
#endif

	if (density[index] <= DensityThreshold){
		continue;
	}
	mass = density[index]*dx*dx*dx;
#if SSDEBUG
	fprintf(stdout, "%s: Excellent! Density threshold exceeeded - density = %g cm^-3\n",
			__FUNCTION__, density[index]*data.DensityUnits/mh);
#endif
	/* 3. Negative divergence: For ZEUS, the velocities are
	   face-centered, and all of the other routines have
	   cell-centered velocities. */

	if (HydroMethod == Zeus_Hydro) {
	  divx = velx[index + offset[0]] - velx[index];
	  divy = vely[index + offset[1]] - vely[index];
	  divz = velz[index + offset[2]] - velz[index];
	  div = divx + divy + divz;
	} else {
	  divx = velx[index + offset[0]] - velx[index - offset[0]];
	  divy = vely[index + offset[1]] - vely[index - offset[1]];
	  divz = velz[index + offset[2]] - velz[index - offset[2]];
	  div = divx + divy + divz;
	}
	/* All three components must be negative to pass the test */
	if (divx > 0.0 || divy > 0.0 || divz > 0.0) continue;
#if SSDEBUG
        fprintf(stderr, "%s: Negative Divergence passed\n", __FUNCTION__);
#endif
	
	/* We now need to define a control volume - this is the region within
	   an accretion radius of the cell identified */
	centralpos[0] = thisGrid->CellLeftEdge[0][i] + 0.5*thisGrid->CellWidth[0][i];
	centralpos[1] = thisGrid->CellLeftEdge[1][j] + 0.5*thisGrid->CellWidth[1][j];
	centralpos[2] = thisGrid->CellLeftEdge[2][k] + 0.5*thisGrid->CellWidth[2][k];

#if COOLING_TIME
#if SSDEBUG
        fprintf(stdout, "%s: Calculate cooling time\n", __FUNCTION__);
#endif

	// 4. t_cool < t_freefall (skip if T > 11000 K)
	dtot = ( density[index] + data.DarkMatterDensity[index] ) *
	  data.DensityUnits;
	tdyn = sqrt(3.0 * M_PI / 32.0 / GravConst / dtot) / data.TimeUnits;
	fprintf(stderr, "%s: Calculate cooling time: %e, and temp = %e.\n", __FUNCTION__, data.CoolingTime[index], data.Temperature[index]);
	if (tdyn < data.CoolingTime[index] &&
	    data.Temperature[index] > 1.1e4)
	  continue;
#endif

#if MASSTHRESHOLDCHECK
#if SSDEBUG
        fprintf(stdout, "%s: Calculate Mass Threshold Check\n", __FUNCTION__);
#endif

	TotalMass = thisGrid->FindMassinGrid(data.DensNum);
	/* Mass Threshold check */
	/* The control region should contain a mass greater than the mass threshold */
	if(TotalMass*ConverttoSolar < (double)MASSTHRESHOLD) {
        continue;
	}
#if SSDEBUG
	printf("%s: Total Mass in Accretion Region = %g Msolar (Threshold = %g)\n", __FUNCTION__,
	       TotalMass*ConverttoSolar, (double)MASSTHRESHOLD);
#endif
#endif
#if MINIMUMPOTENTIAL
#if CALCDIRECTPOTENTIAL
	if(PotentialField == NULL) {
	  PotentialField = new float[size];
	  thisGrid->CalculatePotentialField(PotentialField, data.DensNum, data.DensityUnits, data.TimeUnits,data.LengthUnits);
	}
#endif
	if(PotentialField) {

#if SSDEBUG
	  fprintf(stdout, "%s: Calculate Gravitational Potential\n", __FUNCTION__);
#endif
	  /* 4. Gravitational Minimum Check */
	  /*
	   * Find the minimum of the potential over a Jeans Length.
	   */
	  double JLength = JeansLength(CellTemperature, density[index],
				       data.DensityUnits)/data.LengthUnits;
	  FLOAT search_radius = JLength;
	  GravitationalMinimum  = thisGrid->FindMinimumPotential(centralpos,
				  search_radius,
				  PotentialField);
	  if(PotentialField[index] > GravitationalMinimum) {
#if SSDEBUG
	    printf("FAILURE: GravitationalMinimum = %g\t "		\
		   "PotentialField[index] = %g\n\n", GravitationalMinimum, PotentialField[index]);
	    printf("%s: Cellwidth = %e\t JLength = %e\n", __FUNCTION__, dx, JLength);
	    printf("%s: search_radius (in cellwidths) = %f\n", __FUNCTION__, search_radius/dx);
#endif
	    continue;
	  }
#if SSDEBUG
	  printf("%s: Cellwidth = %e\t JLength = %e\n", __FUNCTION__, dx, JLength);
	  printf("%s: search_radius (in cellwidths) = %f\n", __FUNCTION__, search_radius/dx);
	  fprintf(stdout, "%s: Gravitational Potential Passed!\n", __FUNCTION__);
#endif
	}

#endif
#ifdef GRAVENERGY
	static int mincount = 0;
	/* 5. Jeans Instability Check */

	/* This is the internal energy of the gas */
	ThermalEnergy = thisGrid->FindTotalThermalEnergy(centralpos, LinkingLength*dx,
							 data.GENum);
	/* This is the internal energy plus the kinetic energy */
	CalcTotalEnergy = thisGrid->FindTotalEnergy(centralpos, LinkingLength*dx,
						data.TENum);

	GravitationalEnergy += gpot[index];

	if(fabs(GravitationalEnergy) <= 2*ThermalEnergy) continue;

	/* 6. Bound state check */
	KineticEnergy = thisGrid->FindTotalKineticEnergy(centralpos, LinkingLength*dx,
							 data.DensNum, data.Vel1Num, data.Vel2Num,
							 data.Vel3Num);
	if(KineticEnergy + ThermalEnergy + GravitationalEnergy >= 0.0) continue;

#endif

	/*
	 * Now we need to check the H2 fraction and the accretion rate now.
	 * If the H2fraction > PopIIIH2CriticalFraction then we are ok to form a PopIII star
	 * If H2fraction < PopIIIH2CriticalFraction then we should form a SMS only if the
	 * accretion rate onto the protostar (cell) is above a critical value
	 *
	 */
	float *cellvel = new float[MAX_DIMENSION];
	cellvel[0] = velx[index]; cellvel[1] = vely[index]; cellvel[2] = velz[index];
	float accrate	= thisGrid->ConvergentMassFlow(data.DensNum, data.Vel1Num,
					     dx*ACCRETIONRADIUS, centralpos,
					     cellvel, -1, -1, -1);

	ExtraDensity = density[index] - DensityThreshold;
	 /* We have the main conditions for SF now but we need to decide on the type of star 
	  *  that forms now. There are a few cases to consider:
	  * 1. If the metallicity is high then a cluster of PopII stars can form
	  * 2. If the metallicity is low but the H2 fraction is high then a PopIII star forms
	  * 3. If the metallicty is low and the H2 fraction is low but the accretion rate is high
	  *    then a SMS can form
	  * 4. Otherwise star formation is suppressed
	  * We want to avoid spurious SF in a minihalo that is being heated (e.g. by LW or dynamical 
	  * heating). Therefore we insist on the mass accretion criteria. If the mass accretion rate 
	  * is low we don't allow the SMS pathway to trigger spurious SF. 
	  */
	int stellar_type = -99999;
	// SG. Want to print out level POPIII particle forms on.
	int ThisLevel = thisGrid->GridLevel;

	//if(shu_collapse == 1)
	//  continue;
	if(ProblemType == 27) { //collapse test 
	  stellar_type = POPIII;
	}
	else if(HasMetalField &&    data.TotalMetals[index] > PopIIIMetalCriticalFraction) {
	  stellar_type = POPII;
	}
	else if(data.H2Fraction[index] >  PopIIIH2CriticalFraction) {
	  stellar_type = POPIII;
	  fprintf(stderr, "POPIII particles(%d) created and done in %s on level %"ISYM".\n", data.NumberOfNewParticles + 1,
              __FUNCTION__, ThisLevel);

	}
	else if((accrate*3.154e7*ConverttoSolar/data.TimeUnits > CRITICAL_ACCRETION_RATE*10000000.0)
		&& (dx_pc < SMS_RESOLUTION)) { // SG. Increasing x10 to x10000000 to suppress SMS formation.
	  /* 
	   * The threshold for initially forming the SMS is set at 10 times the critical rates. This 
	   * ensures we get regions of truly high accretion
	   */
	  stellar_type = SMS;
	  fprintf(stderr, "!!!!!!!!SMS Formed\t accrate = %e Msolar/yr",
		 accrate*3.154e7*ConverttoSolar/data.TimeUnits);
	}
	
	if(stellar_type < 0)
	 continue;
	
	ActiveParticleType_SmartStar *np = new ActiveParticleType_SmartStar();
	data.NumberOfNewParticles++;
	data.NewParticles.insert(*np);

	
	np->type = np->GetEnabledParticleID();
	/* 
	 * Set Birthtime to -1.0 initially. If everything is ok and the 
	 * star particle can form then we can update this to its correct value. 
	 */
	np->BirthTime = -1.0;
	np->ParticleClass = stellar_type;
	np->level = data.level;
	np->GridID = data.GridID;
	np->CurrentGrid = thisGrid;
	np->pos[0] = thisGrid->CellLeftEdge[0][i] + 0.5*thisGrid->CellWidth[0][i];
	np->pos[1] = thisGrid->CellLeftEdge[1][j] + 0.5*thisGrid->CellWidth[1][j];
	np->pos[2] = thisGrid->CellLeftEdge[2][k] + 0.5*thisGrid->CellWidth[2][k];

	float *apvel = new float[MAX_DIMENSION];
	/* Calculate AP velocity by taking average of 
	 * surrounding 125 cells 
	 */
	apvel = thisGrid->AveragedVelocityAtCell3D(index, data.DensNum, data.Vel1Num);
	np->vel[0] = apvel[0];
	np->vel[1] = apvel[1];
	np->vel[2] = apvel[2];
	if(ProblemType == 27) { //special case of pure spherical collapse
	  np->vel[0] = cellvel[0]; /* alternatively this could be set to zero */
	  np->vel[1] = cellvel[1]; /* but the velocity of the cell is more consistent */
	  np->vel[2] = cellvel[2];
	}
	
	if (HasMetalField)
	  np->Metallicity = data.TotalMetals[index];
	else
	  np->Metallicity = 0.0;

	np->TimeIndex = 0; //Start at 0 - we'll increment at the start of the update function.
	
	if (np->ParticleClass == POPIII){
		np->AccretionRadius = dx*ACCRETIONRADIUS; // SG. Turning off and on accretion radius for testing purposes.
		// np->AccretionRadius = dx*ACCRETIONRADIUS;
	}else{
		np->AccretionRadius = dx*ACCRETIONRADIUS;
	}
	for(int acc = 1; acc < NTIMES; acc++) {
	  np->AccretionRate[acc] = -11111.0;
	  np->AccretionRateTime[acc] = -11111.0;
	}
	
	/* Calculate initial accretion rate onto cell */
	np->AccretionRate[0] = accrate;
	np->AccretionRateTime[0] = np->BirthTime;
	np->RadiationLifetime= 0.0; 
	np->StellarAge = 0.0;
	np->NotEjectedMass = 0.0;
       
	/* The mass of the particles formed depends on the resolution and is handled 
	 * in a call to `RemoveMassFromGridAfterFormation` which is called from 
	 * `AfterEvolveLevel`. np->Mass is initally set here but can get overwritten 
	 *  AfterFormation. The cell values are not updated here but 
		*
	 * instead are updated in `AfterFormation`.
	 * We set the initial mass to zero so that merging works ok i.e. nothing spurious 
	 * happens in this case. 
	 */
	
	np->Mass = 0.0;
	np->oldmass = 0.0; // SG.
	fprintf(stderr,"%s: Total Mass in Accretion Region = %g Msolar (Threshold = %g)\n", __FUNCTION__,
	     TotalMass*ConverttoSolar, (double)MASSTHRESHOLD);
    } // end else
      } // i
    } // j
  } // k
  if(data.NumberOfNewParticles > 0) {
    fprintf(stderr, "Particles (%d) Created and done in %s\n", data.NumberOfNewParticles, __FUNCTION__);
  }

#if CALCDIRECTPOTENTIAL
  delete [] PotentialField;
#endif
  return SUCCESS;
}


void ActiveParticleType_SmartStar::DescribeSupplementalData(ActiveParticleFormationDataFlags &flags)
{
  flags.DarkMatterDensity = true;
  flags.Temperature = true;
  flags.MetalField = true;
  if(MultiSpecies > 1)
    flags.H2Fraction = true;
#if COOLING_TIME
  flags.CoolingTime = true;
#else
  flags.CoolingTime = false;
#endif
}


/*
 * Do feedback before evolving the level forward 
 */
template <class active_particle_class>
int ActiveParticleType_SmartStar::BeforeEvolveLevel
(HierarchyEntry *Grids[], TopGridData *MetaData,
 int NumberOfGrids, LevelHierarchyEntry *LevelArray[], 
 int ThisLevel, bool CallEvolvePhotons, 
 int TotalStarParticleCountPrevious[],
 int SmartStarID)
{
    FLOAT Time = LevelArray[ThisLevel]->GridData->ReturnTime();
    float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,VelocityUnits;
    double MassUnits;
    GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
             &TimeUnits, &VelocityUnits, Time);
    MassUnits = DensityUnits * POW(LengthUnits,3);
    int j, dim, ipart, nParticles;
    ActiveParticleList<ActiveParticleType> SmartStarList;
    if(SmartStarFeedback == FALSE)
        return SUCCESS;
    if(SmartStarBHRadiativeFeedback == FALSE && SmartStarStellarRadiativeFeedback == FALSE)
        return SUCCESS;
    if (CallEvolvePhotons)
        //fprintf(stderr,"%s: In if CallEvolvePhotons == true, the ActiveParticleFindAll triggered.\n", __FUNCTION__);
        ActiveParticleFindAll(LevelArray, &nParticles, SmartStarID,
                              SmartStarList);

#ifdef TRANSFER
  if (CallEvolvePhotons) {
			 //fprintf(stderr,"%s: TRANSFER being called. Radiative feedback should be carried out.\n", __FUNCTION__);
    RadiationSourceEntry* source;
    double dx;
    double MassConversion;
    const double LConv = (double) TimeUnits / pow(LengthUnits,3);
    const double mfactor = double(MassUnits) / SolarMass;
    
    ActiveParticleType_SmartStar *ThisParticle = NULL;
    for (ipart = 0; ipart < nParticles; ipart++) {
      if (SmartStarList[ipart]->IsARadiationSource(Time)) {
	ThisParticle =
	  static_cast<ActiveParticleType_SmartStar*>(
          SmartStarList[ipart]);
	
	if(ThisParticle->ParticleClass == BH && SmartStarBHRadiativeFeedback == FALSE)
	  continue; //No BH radiative feedback
	if( (ThisParticle->ParticleClass == POPIII || ThisParticle->ParticleClass == SMS) 
	    && SmartStarStellarRadiativeFeedback == FALSE)
	  continue; //No stellar radiative feedback

	dx = LevelArray[ThisParticle->level]->GridData->GetCellWidth(0,0);
	//fprintf(stderr, "%s: cell width on particle grid (level = %"ISYM") is %e pc.\n", __FUNCTION__, ThisParticle->level, dx*LengthUnits/pc_cm);
	MassConversion = (double) (dx*dx*dx * mfactor); //Converts to Solar Masses
	source = ThisParticle->RadiationSourceInitialize();
	double PMass = ThisParticle->Mass*MassConversion;
	float ramptime = 0.0;
	//fprintf(stderr, "%s: PMass = %e and Particle Class = %"ISYM" and OldMass = %e.\n", __FUNCTION__, PMass, ThisParticle->ParticleClass, ThisParticle->oldmass);
	if(POPIII == ThisParticle->ParticleClass ||
	   SMS == ThisParticle->ParticleClass) {
	  ramptime = yr_s * 1e5 / TimeUnits; // SG. Change to 5e4?
	}
	if(POPII == ThisParticle->ParticleClass) {
	  ramptime = yr_s * StarClusterMinDynamicalTime / TimeUnits;
	}
	/* Call Function to return SED parameters */
	if(ThisParticle->DetermineSEDParameters(Time, dx) == FAIL)
	  return FAIL;

	source->LifeTime       = ThisParticle->RadiationLifetime; 
	source->Luminosity = (ThisParticle->LuminosityPerSolarMass * LConv) * PMass;
	source->RampTime  =  ramptime;
	source->EnergyBins = RadiationSEDNumberOfBins;
	source->Energy = new float[RadiationSEDNumberOfBins];
	source->SED = new float[RadiationSEDNumberOfBins];
	for (j = 0; j < RadiationSEDNumberOfBins; j++) {
	  source->Energy[j] = ThisParticle->RadiationEnergyBins[j];
	  source->SED[j] = ThisParticle->RadiationSED[j];
	  
	}

	// SG. For debugging purposes.
	if(ThisParticle->ParticleClass == BH) {
	  printf("%s: !!!!!!!!!!!!!!!!!!!!!!!!SRC Luminosity: L=%lg Lcode=%g M=%g Mcode=%g\n", __FUNCTION__, 
		 ThisParticle->LuminosityPerSolarMass * ThisParticle->Mass * MassConversion, 
		 source->Luminosity, ThisParticle->Mass * MassConversion, ThisParticle->Mass);
	  printf("TimeIndex = %d\n", ThisParticle->TimeIndex);
	  printf("%s: BH Mass = %e Msolar AccretionRate = %e Msolar/yr\n", __FUNCTION__,
		 ThisParticle->Mass * MassConversion, 
		 (ThisParticle->AccretionRate[ThisParticle->TimeIndex]*MassConversion/TimeUnits)*yr_s);
	  printf("%s: ParticleClass = %d\t SEDs = [%f, %f, %f, %f, %f]\n", __FUNCTION__,
		 ThisParticle->ParticleClass, ThisParticle->RadiationSED[0], ThisParticle->RadiationSED[1],
		 ThisParticle->RadiationSED[2], ThisParticle->RadiationSED[3],
		 ThisParticle->RadiationSED[4]);
	}

#if SSDEBUG_RT
	printf("%s: SS (%d) Energy Bins = %1.1f (%1.1f) %1.1f (%1.1f) %1.1f (%1.1f) " \
	       "%1.1f (%1.1f) %1.1f (%1.1f)\n", __FUNCTION__,
	       ThisParticle->ParticleClass, ThisParticle->RadiationEnergyBins[0],
	       ThisParticle->RadiationSED[0],
	       ThisParticle->RadiationEnergyBins[1],ThisParticle->RadiationSED[1],
	       ThisParticle->RadiationEnergyBins[2],ThisParticle->RadiationSED[2],
	       ThisParticle->RadiationEnergyBins[3],ThisParticle->RadiationSED[3],
	       ThisParticle->RadiationEnergyBins[4],ThisParticle->RadiationSED[4]);
	

	
	if(ThisParticle->ParticleClass == SMS) {
	  printf("%s: !!!!!!!!!!!!!!!!!!!!!!!!SRC Luminosity: L=%lg Lcode=%g M=%g Mcode=%g\n", __FUNCTION__, 
		 ThisParticle->LuminosityPerSolarMass * ThisParticle->Mass * MassConversion, 
		 source->Luminosity, ThisParticle->Mass * MassConversion, ThisParticle->Mass);
	  printf("TimeIndex = %d\n", ThisParticle->TimeIndex);
	  printf("%s: SMS Mass = %e Msolar AccretionRate = %e Msolar/yr\n", __FUNCTION__,
		 ThisParticle->Mass * MassConversion, 
		 (ThisParticle->AccretionRate[ThisParticle->TimeIndex]*MassConversion/TimeUnits)*yr_s);
	  printf("%s: ParticleClass = %d\t SEDs = [%f, %f, %f, %f, %f]\n", __FUNCTION__,
		 ThisParticle->ParticleClass, ThisParticle->RadiationSED[0], ThisParticle->RadiationSED[1],
		 ThisParticle->RadiationSED[2], ThisParticle->RadiationSED[3],
		 ThisParticle->RadiationSED[4]);
	}
#endif
	if (GlobalRadiationSources->NextSource != NULL)
	  GlobalRadiationSources->NextSource->PreviousSource = source;
	GlobalRadiationSources->NextSource = source;
      }
    }
   
  } // ENDIF CallEvolvePhotons
#endif /* TRANSFER */
    return SUCCESS;
}

int ActiveParticleType_SmartStar::PopIIIFormationFromSphere(ActiveParticleType_SmartStar* SS, 
	grid* APGrid, int ThisProcessorNum, FLOAT StarLevelCellWidth, FLOAT CellVolumeStarLevel,
	FLOAT Time, LevelHierarchyEntry *LevelArray[], LevelHierarchyEntry *Temp, int ThisLevel)
	{

		/* Set the units. */
		float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits;
		double MassUnits;
		GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
			&TimeUnits, &VelocityUnits, Time);
		MassUnits = DensityUnits * POW(LengthUnits,3);
		
		/* Initialise attributes of PopIII class
		+ intstantiate sphere variables */
		float Radius = 0.0, SphereRadius, SphereMass, Age;
		float MassEnclosed = 0, ColdGasMass = 0, ColdGasFraction = 0;
		float Metallicity2 = 0, Metallicity3 = 0, Subtraction;
		int SphereContained, SphereContainedNextLevel, CellsModified = 0;
		bool MarkedSubgrids = false;

		/* Set target sphere mass to be twice the target PopIII mass */
		float TargetSphereMass = 2*PopIIIStarMass;

		/* Find radius of sphere to be accreted from */
		SS->FindAccretionSphere(LevelArray, ThisLevel, StarLevelCellWidth, Radius,
                                TargetSphereMass, MassEnclosed, Metallicity2, Metallicity3, ColdGasMass,
                                ColdGasFraction, SphereContained, MarkedSubgrids);

		/* Determine if a sphere is enclosed within the grids on next level
		If that is the case, we perform SubtractAccretedMass not here, 
		but in the EvolveLevel of the next level. */
		SphereMass = MassEnclosed; // for mass removal
		SphereRadius = Radius; // For subtraction variable
		SphereContainedNextLevel = FALSE;
		
		if (LevelArray[ThisLevel+1] != NULL) {
			fprintf(stderr, "%s: Checking if sphere can be found on next highest level.\n", __FUNCTION__);		
			SS->FindAccretionSphere(LevelArray, ThisLevel+1, StarLevelCellWidth, Radius,
                                    TargetSphereMass,MassEnclosed, Metallicity2, Metallicity3, ColdGasMass,
                                    ColdGasFraction,SphereContainedNextLevel, MarkedSubgrids);
		}

		/* Quit this routine when 
		(1) sphere is not contained, or 
		(2) sphere is contained, but the next level can contain the sphere too. */ 

		if ((SphereContained == FALSE) || 
		(SphereContained == TRUE && SphereContainedNextLevel == TRUE)){
			fprintf(stderr, "%s: Sphere not contained or sphere contained on next level.\n", 
				__FUNCTION__, SphereContained);	
			return SUCCESS;
		}

		/* Set fraction of SphereMass that will be removed */
		Subtraction = PopIIIStarMass/SphereMass;

		/* Now set cells within the radius to their values after mass subtraction. */
		for (int l = ThisLevel; l < MAX_DEPTH_OF_HIERARCHY; l++){
			for (Temp = LevelArray[l]; Temp; Temp = Temp->NextGridThisLevel){

				Temp->GridData->RemoveMassFromSphere(SS, l, SphereRadius, DensityUnits, 
					LengthUnits, VelocityUnits, TemperatureUnits, 
					TimeUnits, Subtraction, CellsModified);

			} // END Grids
		} // END Levels

		/* Exit if no mass is enclosed */
		if (SphereMass == 0){
			fprintf(stderr,"%s: No mass enclosed. No particle created.\n", __FUNCTION__);
			return SUCCESS;
		}

		fprintf(stderr,"%s: PopIII Particle created!\n", __FUNCTION__);

		/* Assign mass, radius and lifetime to particle */
		SS->Mass = PopIIIStarMass; // msun
		SS->InfluenceRadius = SphereRadius; // code units
		//SS->RadiationLifetime = CalculatePopIIILifetime(SS->Mass); // code time
		//SS->RadiationLifetime*= yr_s/TimeUnits;
		SS->RadiationLifetime = 0; // SG. Hardcoding lifetime to zero for BH-phase only runs. Replaces above two lines.
		SS->BirthTime = APGrid->ReturnTime();
		Age = Time - SS->BirthTime;

		/* Print out particle attributes */
		fprintf(stderr,"%s: Particle Mass = %1.1f Msolar\n", __FUNCTION__, SS->Mass);
		fprintf(stderr,"%s: Particle Lifetime = %1.2f Myr\n", __FUNCTION__,
                SS->RadiationLifetime*TimeUnits/Myr_s);
		fprintf(stderr,"%s: Particle Age = %1.1f Myr\n", __FUNCTION__, Age*TimeUnits/Myr_s);
		fprintf(stderr,"%s: Particle Class = %d\n", __FUNCTION__, SS->ParticleClass);
		fprintf(stderr,"%s: Removed mass from sphere of radius %e pc\n", __FUNCTION__,
                SS->InfluenceRadius*LengthUnits/pc_cm);

		/* Assign code density mass to SS and oldmass */
		SS->Mass = (PopIIIStarMass*SolarMass/MassUnits)/CellVolumeStarLevel; //code density
		SS->oldmass = 0.0;

		return SUCCESS;
} // END POPIII


grid* ConstructFeedbackZone(ActiveParticleType* ThisParticle, FLOAT FeedbackRadius,
			     FLOAT dx, HierarchyEntry** Grids, int NumberOfGrids,
			     int SendField);

int DistributeFeedbackZone(grid* FeedbackZone, HierarchyEntry** Grids,
			   int NumberOfGrids, int SendField);

int ActiveParticleType_SmartStar::EvaluateFeedback(
        grid *thisgrid_orig, ActiveParticleFormationData &data) {
    // SG. Feedback is handled elsewhere.
    return SUCCESS;
}

// SG. Altered
int ActiveParticleType_SmartStar::RemoveMassFromGridAfterFormation(int nParticles, 
    ActiveParticleList<ActiveParticleType>& ParticleList,
    LevelHierarchyEntry *LevelArray[], int ThisLevel)
{

	/* Set the units. */
	FLOAT Time = LevelArray[ThisLevel]->GridData->ReturnTime();
    float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits;
    double MassUnits;
	int SSparticles[nParticles] = {-1};
    GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
             &TimeUnits, &VelocityUnits, Time);
    MassUnits = DensityUnits * POW(LengthUnits,3);

	/*
		* Order particles in order of SMS, PopIII, PopII
		* SMS first since they have the highest accretion rates and hence 
		* should be forming first
	*/
	int k = 0, num_new_sms_stars = 0, num_new_popiii_stars = 0, num_new_popii_stars = 0;
	for (int i = 0; i < nParticles; i++) {
        grid* APGrid = ParticleList[i]->ReturnCurrentGrid();
		if (MyProcessorNumber == APGrid->ReturnProcessorNumber()) {
			ActiveParticleType_SmartStar* SS;
			SS = static_cast<ActiveParticleType_SmartStar*>(ParticleList[i]);
			if(SS->ParticleClass == SMS && SS->TimeIndex == 0) {
				SSparticles[k++] = i;
				num_new_sms_stars++;
			} 
		} 
	}

	for (int i = 0; i < nParticles; i++) {
		grid* APGrid = ParticleList[i]->ReturnCurrentGrid();
		ActiveParticleType_SmartStar* SS;
		SS = static_cast<ActiveParticleType_SmartStar*>(ParticleList[i]);
		// SG. For low resolution particles that never get accreted, the time index is never incremented
		// Hence Mass = 0 check is included. No proc check here anymore.
		if(SS->ParticleClass == POPIII && SS->TimeIndex == 0 && SS->Mass == 0) {
			SSparticles[k++] = i;
			num_new_popiii_stars++;
		} // End IF particle class POPIII
	} // End FOR
					
	for (int i = 0; i < nParticles; i++) {
		grid* APGrid = ParticleList[i]->ReturnCurrentGrid();
		if (MyProcessorNumber == APGrid->ReturnProcessorNumber()) {
			ActiveParticleType_SmartStar* SS;
			SS = static_cast<ActiveParticleType_SmartStar*>(ParticleList[i]);
			if(SS->ParticleClass == POPII && SS->TimeIndex == 0) {
				SSparticles[k++] = i;
				num_new_popii_stars++;
			}
		}
	}

	int num_new_stars = num_new_sms_stars + num_new_popiii_stars + num_new_popii_stars;
	if(num_new_stars == 0){
		return SUCCESS;
	}


    /* Instantiate common attributes of all SS particle types */
	ActiveParticleType_SmartStar* SS;
	LevelHierarchyEntry *Temp;
	grid* APGrid;
	int pindex, ThisProcessorNum;
	int cellindex_x, cellindex_y, cellindex_z, cellindex;
	float DensityThreshold,
	ParticleDensity, newcelldensity;
	FLOAT StarLevelCellWidth, CellVolumeStarLevel, dx_pc;
	float *density;

#if JEANSREFINEMENT
	bool JeansRefinement = false;
	for (int method = 0; method < MAX_FLAGGING_METHODS; method++) 
		if (CellFlaggingMethod[method] == 6)
			JeansRefinement = true;
			if (JeansRefinement) {
				int size = APGrid->GetGridSize();
				float *Temperature = new float[size]();
				APGrid->ComputeTemperatureField(Temperature);
				float CellTemperature = (JeansRefinementColdTemperature > 0) ? JeansRefinementColdTemperature : Temperature[cellindex];
				int JeansFactor = JEANS_FACTOR; 
				float JeansDensityUnitConversion = (Gamma*pi*kboltz) / (Mu*mh*GravConst);
				float JeansDensity = JeansDensityUnitConversion * 1.01 * CellTemperature /
					POW(LengthUnits*dx*JeansFactor,2);
				JeansDensity /= DensityUnits;
				DensityThreshold = min(DensityThreshold,JeansDensity);
				fprintf(stderr,"%s: Density Threshold = %e\t Jeans Density = %e\n", __FUNCTION__, DensityThreshold, JeansDensity);
				}
#endif

    /* Main loop over new particles */
	for (int k = 0; k < num_new_stars; k++){

		/* Define APGrid and SS particle */
		pindex = SSparticles[k];
		APGrid = ParticleList[pindex]->ReturnCurrentGrid();
		ThisProcessorNum = APGrid->ReturnProcessorNumber();
		SS = static_cast<ActiveParticleType_SmartStar*>(ParticleList[pindex]);

		/* Define cell width and volume on star grid */
		StarLevelCellWidth = APGrid->CellWidth[0][0];
		CellVolumeStarLevel = StarLevelCellWidth*StarLevelCellWidth*StarLevelCellWidth;
		dx_pc = StarLevelCellWidth*LengthUnits/pc_cm; // in pc

		/* Define cellindex on APGrid */
		cellindex_x = (SS->pos[0] - APGrid->CellLeftEdge[0][0])/StarLevelCellWidth,
		cellindex_y = (SS->pos[1] - APGrid->CellLeftEdge[1][0])/StarLevelCellWidth,
		cellindex_z = (SS->pos[2] - APGrid->CellLeftEdge[2][0])/StarLevelCellWidth;
		cellindex = APGrid->GetIndex(cellindex_x, cellindex_y, cellindex_z);


		/*********************************************************************
		                            	SMS CASE
		**********************************************************************/

		if(SMS == SS->ParticleClass) {
			/* 
			If resolution is sufficent, accrete as normal - 
			just remove mass from the cell 
			*/

		    /* Find DensNum */
			int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
			if (APGrid->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				Vel3Num, TENum) == FAIL){
					ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
					}

		    /* Define ParticleDensity and newcelldensity (only known to owning processor) */
			density = APGrid->BaryonField[DensNum];
			DensityThreshold = ActiveParticleDensityThreshold*mh/DensityUnits;
			ParticleDensity = density[cellindex] - DensityThreshold;
			newcelldensity = density[cellindex] - ParticleDensity;
			
			if(dx_pc <= SMS_RESOLUTION) {

				/* Update density, set SS attributes */
				density[cellindex] = newcelldensity;
				SS->BirthTime = APGrid->ReturnTime();
				SS->Mass = ParticleDensity;
				SS->oldmass = 0.0;

				/* Pathological case */
				if(ParticleDensity < 0.0) {
					printf("%s: cellindex = %d\n", __FUNCTION__, cellindex);
					printf("density[cellindex] = %e cm^-3\n", density[cellindex]*DensityUnits/mh);
					printf("DensityThreshold = %e cm^-3\n", DensityThreshold*DensityUnits/mh);
					printf("SS->ParticleClass = %d\n", SS->ParticleClass); fflush(stdout);
					ENZO_FAIL("Particle Density is negative. Oh dear.\n");
				}

			 continue;
			} // END resolution check
			continue; // SG. No low-res SMS case implemented yet.
		} // END SMS


		/*********************************************************************
		                            	POP III CASE
		**********************************************************************/

		else if (SS->ParticleClass == POPIII){
			/* 
			SG. If stellar accretion is off, pop3 star forms with final mass taken
			from sphere. The target pop3 stellar mass is set with the parameter
			'PopIIIStarMass' which is set by the user in the parameter file.
			A sphere will step out by one cell width until it reaches a radius which
			contains twice the PopIIIStarMass value. 
			*/
		    #if STELLAR_ACCRETION_OFF
			PopIIIFormationFromSphere(SS, APGrid, ThisProcessorNum, StarLevelCellWidth,
			CellVolumeStarLevel, Time, LevelArray, Temp, ThisLevel);
			continue;
		    #endif

		    /* If stellar accretion not off and resolution is sufficent, accrete as normal -
		     * just remove mass from the cell */
            if (MyProcessorNumber != APGrid->ReturnProcessorNumber()){
                continue;
            }

			/* Find DensNum */
            int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
            if (APGrid->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,Vel3Num, TENum) == FAIL){
                    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
                    }

            /* Define ParticleDensity and newcelldensity (only known to owning processor) */
            density = APGrid->BaryonField[DensNum];
            DensityThreshold = ActiveParticleDensityThreshold*mh/DensityUnits;
            ParticleDensity = density[cellindex] - DensityThreshold;
            newcelldensity = density[cellindex] - ParticleDensity;

            /* Update density, set SS attributes */
            if(dx_pc <= POPIII_RESOLUTION) {
               density[cellindex] = newcelldensity;
               SS->BirthTime = APGrid->ReturnTime();
               SS->Mass = ParticleDensity;
               SS->oldmass = 0.0;

               /* Pathological case */
               if(ParticleDensity < 0.0) {
                   printf("%s: cellindex = %d\n", __FUNCTION__, cellindex);
                   printf("density[cellindex] = %e cm^-3\n", density[cellindex]*DensityUnits/mh);
                   printf("DensityThreshold = %e cm^-3\n", DensityThreshold*DensityUnits/mh);
                   printf("SS->ParticleClass = %d\n", SS->ParticleClass); fflush(stdout);
                   ENZO_FAIL("Particle Density is negative. Oh dear.\n");
               }
                #if SSDEBUG
                fprintf(stderr, "%s: Particle with initial mass %e (%e) Msolar created\n", __FUNCTION__,
                        SS->Mass*pow(StarLevelCellWidth,3)*MassUnits/SolarMass, SS->Mass);
                #endif
                continue;
			} // END resolution check
			else {
				PopIIIFormationFromSphere(SS, APGrid, ThisProcessorNum, StarLevelCellWidth,
                                          CellVolumeStarLevel, Time, LevelArray, Temp, ThisLevel);
				continue;
            }
		} // END PopIII


		/*********************************************************************
		                            	POP II CASE
		**********************************************************************/

		else if (POPII == SS->ParticleClass) {
			/* 
			 Accrete if the mass exceeds the star cluster minimum mass.
			*/

		    /* Find DensNum */
			int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
			if (APGrid->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				Vel3Num, TENum) == FAIL){
					ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
					}

			/* Define ParticleDensity and newcelldensity (only known to owning processor) */
			density = APGrid->BaryonField[DensNum];
			DensityThreshold = ActiveParticleDensityThreshold*mh/DensityUnits;
			ParticleDensity = density[cellindex] - DensityThreshold;
			newcelldensity = density[cellindex] - ParticleDensity;

			float PopIIMass = SS->Mass*pow(StarLevelCellWidth,3)*MassUnits/SolarMass;
			if (PopIIMass > StarClusterMinimumMass){

				/* Update densities and set SS attributes */
				density[cellindex] = (1 - StarClusterFormEfficiency)*density[cellindex];
				SS->BirthTime = APGrid->ReturnTime();
				SS->Mass = StarClusterFormEfficiency*density[cellindex];
				SS->oldmass = 0.0;

				/* Pathological case */
				if(ParticleDensity < 0.0) {
					fprintf(stderr,"%s: cellindex = %d\n", __FUNCTION__, cellindex);
					fprintf(stderr,"density[cellindex] = %e cm^-3\n", density[cellindex]*DensityUnits/mh);
					fprintf(stderr,"DensityThreshold = %e cm^-3\n", DensityThreshold*DensityUnits/mh);
					fprintf(stderr,"SS->ParticleClass = %d\n", SS->ParticleClass); fflush(stdout);
					ENZO_FAIL("Particle Density is negative. Oh dear.\n");
				}
				printf("POPII: cellindex %d updated - next.\n\n", cellindex);

				continue;
			} // END mass check
		 continue; // No low-res POPII case implemented yet.
         } // END POPII
         } // END num_new_stars loop
	return SUCCESS;
} // END RemoveMassFromGridAfterFormation


int ActiveParticleType_SmartStar::Accrete(int nParticles, 
    ActiveParticleList<ActiveParticleType>& ParticleList,
    FLOAT &AccretionRadius,
    LevelHierarchyEntry *LevelArray[], int ThisLevel, int SmartStarID)
{
    /* Skip accretion if we're not on the maximum refinement level.
    This should only ever happen right after creation and then
    only in pathological cases where sink creation is happening at
    the edges of two regions at the maximum refinement level
     */

    // if (ThisLevel < MaximumRefinementLevel)
    //   return SUCCESS;

    FLOAT Time = LevelArray[ThisLevel]->GridData->ReturnTime();
    float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits;
    double MassUnits;
    GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits,
             &VelocityUnits, Time);
    MassUnits = DensityUnits * POW(LengthUnits,3);
    // SG. Find dx from grid data.
    FLOAT dx = LevelArray[ThisLevel]->GridData->CellWidth[0][0];
    FLOAT dx_pc = dx*LengthUnits/pc_cm;   //in pc
    double MassConversion = (double) (dx*dx*dx * double(MassUnits));  //convert to g

    // For each particle, loop over all the grids and do accretion if the grid overlaps with the accretion zone
    int i, NumberOfGrids;
    int *FeedbackRadius = NULL;
    HierarchyEntry **Grids = NULL;
    grid *sinkGrid = NULL;

    bool SinkIsOnThisProc, SinkIsOnThisGrid;

    float SubtractedMass, SubtractedMomentum[3] = {};

    NumberOfGrids = GenerateGridArray(LevelArray, ThisLevel, &Grids);
    float ctime = LevelArray[ThisLevel]->GridData->ReturnTime();

    /*
   * TimeDelay if we want to delay the time between 
   * turning on accretion following a supernova. 
   * In testing I didn't observe any difference. 
   * Hence, the TimeDelay is set to a short period
   * of 100,000 years which accounts for end of snowplough 
   * period
    */
    float TimeDelay = 0*yr_s/TimeUnits; //set to 100 kyr - SG change to 0 yr for BH only.
    for (i = 0; i < nParticles; i++) {
        /*
         * Accretion is only allowed if it makes sense:
         * 1. Black Holes in general are always allowed to accrete.
         * 2. SMS can accrete if there is sufficient resolution to do so
         * 3. PopIII stars can accrete if there is sufficient resolution
         * 4. PopII stars never accrete
         */
        ActiveParticleType_SmartStar* SS;
        SS = static_cast<ActiveParticleType_SmartStar*>(ParticleList[i]);
        float MassInSolar = ParticleList[i]->ReturnMass()*MassConversion/SolarMass;
        AccretionRadius = SS->AccretionRadius;
        int pclass = SS->ParticleClass;
        float Stellar_Age = SS->StellarAge;
        float p_age = ctime - SS->BirthTime;
        int MyLevel = SS->ReturnLevel(); // SG.

        //SG.
        if (ThisLevel != MyLevel){
            continue;
            }

        if(pclass == POPIII){
            /*
             * We only accrete onto POPIII stars if our maximum
             * spatial resolution is better than 1e-3 pc
             * */
            if(dx_pc > POPIII_RESOLUTION) // we don't have sufficient resolution
              continue;
            // SG.
            if (MassInSolar >= PopIIIStarMass || MassInSolar == 0){
                fprintf(stderr,"%s: Accrete is skipped.\n", __FUNCTION__);
                continue; // SG. Change to return SUCCESS? In order to skip particle assignment
            }
        }
        else if(pclass == SMS){
            /*
             * We only accrete onto SMSs if our maximum
             * spatial resolution is better than 1e-1 pc
             * */
            if(dx_pc > SMS_RESOLUTION) // we don't have sufficient resolution
                continue;
        }
        else if(pclass == BH){
            /* We always accrete onto BHs. The only restriction is that
             * we can optionally employ a time delay following a SNe explosion to
             * avoid spurious accretion.
            */
            // time delay.
            if(p_age < Stellar_Age + TimeDelay){
                fprintf(stderr, "%s: no accretion due to TimeDelay (%e yrs). p_age = %e yrs. "
                                "Stellar_Age = %e yrs.\n",
                        __FUNCTION__, TimeDelay*TimeUnits/yr_s, p_age*TimeUnits/yr_s, Stellar_Age*TimeUnits/yr_s);
                continue;
            }
        }
        else if(pclass == POPII){
            /* We never accrete onto POPII stars */
            continue;
        }

        // SG/BS Put feedback zone inside processor num
        grid* FeedbackZone = ConstructFeedbackZone(ParticleList[i], int(AccretionRadius/dx), dx, Grids,
                                                   NumberOfGrids, ALL_FIELDS);
        grid* APGrid = ParticleList[i]->ReturnCurrentGrid();

        // SG. Set to 0 before it's calculated by owning proc and then communicated with other procs in
        // CommunicateAllSumValues().
        FLOAT positions[3] = {0,0,0}; // SG. All elements initialised to zero.
        FLOAT NewAccretionRadius = 0;
        FLOAT* pos;

        // SG/BS change to continue and !=.
        if (MyProcessorNumber == FeedbackZone->ReturnProcessorNumber()) {
            float AccretionRate = 0;

            if (FeedbackZone->AccreteOntoSmartStarParticle(ParticleList[i],
			      AccretionRadius, &AccretionRate) == FAIL)
                return FAIL;

            // SG. positions is the array of de-referenced particle positions in each dim.
            positions[0] = ParticleList[i]->ReturnPosition()[0];
            positions[1] = ParticleList[i]->ReturnPosition()[1];
            positions[2] = ParticleList[i]->ReturnPosition()[2];


#if BONDIHOYLERADIUS
            /* Check what the Bondi-Hoyle radius is - we should accrete out to that if required */
            float mparticle = ParticleList[i]->ReturnMass()*dx*dx*dx;
            float *vparticle = ParticleList[i]->ReturnVelocity();
            int size = APGrid->GetGridSize();
            float *Temperature = new float[size]();
            APGrid->ComputeTemperatureField(Temperature);

            // SG print
            fprintf(stderr, "%s: mparticle = %e Msun.\n", __FUNCTION__, mparticle*MassUnits/SolarMass);

            FLOAT BondiHoyleRadius = APGrid->CalculateBondiHoyleRadius(mparticle, vparticle, Temperature);

            fprintf(stderr, "%s: BondiHoyleRadius = %e pc.\n", __FUNCTION__,
                    BondiHoyleRadius*LengthUnits/pc_cm);
            fprintf(stderr, "%s: AccretionRadius = %e pc.\n", __FUNCTION__,
                    SS->AccretionRadius*LengthUnits/pc_cm);
            /*
            SG Comment: the accretion radius will be reassigned to the Bondi radius when the accretion radius
            falls outside of the tolerance.
            */
            FLOAT tol = 0.000001*BondiHoyleRadius;
            if(BondiHoyleRadius + tol < SS->AccretionRadius || SS->AccretionRadius < BondiHoyleRadius - tol) {
                SS->AccretionRadius = BondiHoyleRadius;
                fprintf(stderr, "%s: Updating accretion radius to Bondi radius = %e pc (%f cells)\n",
                        __FUNCTION__,SS->AccretionRadius*LengthUnits/pc_cm,SS->AccretionRadius/dx);
            }
            else{
                fprintf(stderr, "%s: Accretion Radius = %e pc (%f cells). No update needed.\n",
                        __FUNCTION__,SS->AccretionRadius*LengthUnits/pc_cm, SS->AccretionRadius/dx);
            }
            AccretionRadius = SS->AccretionRadius;
            NewAccretionRadius = AccretionRadius;
            delete [] Temperature;
            Temperature = NULL;
#endif
	/* No need to communicate the accretion rate to the other CPUs since this particle is already local.
	 * Need to decide how often I update the accretion history
	 * */
        } // SG. End processor.

        // SG. Communicate with all processors the updated accretion radius.
        CommunicationAllSumValues(&NewAccretionRadius, 1);
        CommunicationAllSumValues(positions, 3);

        SS->AccretionRadius = NewAccretionRadius;
        AccretionRadius = SS->AccretionRadius;

        SS->pos[0] = positions[0];
        SS->pos[1] = positions[1];
        SS->pos[2] = positions[2];

        DistributeFeedbackZone(FeedbackZone, Grids, NumberOfGrids, ALL_FIELDS);
        delete FeedbackZone;
    } // END particles

    if (AssignActiveParticlesToGrids(ParticleList, nParticles, LevelArray) == FAIL)
    return FAIL;

    delete [] Grids;
    return SUCCESS;
} // END Accrete


int ActiveParticleType_SmartStar::SetFlaggingField(
    LevelHierarchyEntry *LevelArray[], int level,
    int TopGridDims[], int SmartStarID)
{
    /* 	SG. Get dx of grid cell here. This function is calling all grids on level.
     * 1) accrad, 2) cell width, 3) parameter for number cells to refine the accretion radius by.
     * Only if dx > dx_bondi is DepositRefinementZone triggered.
     * */

    /* Generate a list of all sink particles in the simulation box */
    int i, nParticles;
    FLOAT *pos = NULL;
    ActiveParticleList<ActiveParticleType> SmartStarList;
    LevelHierarchyEntry *Temp = NULL;
    double dx = LevelArray[level]->GridData->CellWidth[0][0]; // SG. Grid cell width.

    // SG. Get units for conversion to pc for dx in print statement.
    FLOAT Time = LevelArray[level]->GridData->ReturnTime();
    float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits;
    double MassUnits;
    GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
             &TimeUnits, &VelocityUnits, Time);

    // Find all smartstar particles
    ActiveParticleFindAll(LevelArray, &nParticles, SmartStarID,
                          SmartStarList);

    for (i=0 ; i<nParticles; i++){
        ActiveParticleType_SmartStar* SS;
        SS = static_cast<ActiveParticleType_SmartStar*>(SmartStarList[i]);
        int pclass = SS->ParticleClass;

        /* POPIII case*/
        if (pclass == POPIII){
            continue;
		}
		 /* SMS case*/
	    if (pclass == SMS) {
            continue;
		  }
        /* BH case*/
        else{
			/* Define position and accrad of BH */
			pos = SmartStarList[i]->ReturnPosition();
			double accrad = static_cast<ActiveParticleType_SmartStar*>(SmartStarList[i])->AccretionRadius;

			/* SG. Only print out accretion radius if we're on the SS grid */
			grid* APGrid = SS->ReturnCurrentGrid();
			if (MyProcessorNumber == APGrid->ReturnProcessorNumber()){
			fprintf(stderr, "%s: Accretion radius = %e (Bondi radius) and bondi factor = %e and cell_width = %e.\n",
                    __FUNCTION__, accrad, SmartStarBondiRadiusRefinementFactor, dx);
            }

			/* SG. Check for when accrad = 0 in the first 100 kyr of BH's life. */
			if (accrad < 1e-30)
			continue;

			/* SG. Calculate user-set dx_bondi and dx_bondi in pc*/
			double dx_bondi = (double) accrad/SmartStarBondiRadiusRefinementFactor;
			double dx_pc = dx*LengthUnits/pc_cm;   //in pc
			double dx_bondi_pc = dx_bondi*LengthUnits/pc_cm; //in pc
	
			/* if dx_bondi > dx, don't deposit refinement zone */
			if (dx_bondi > dx){
				//fprintf(stderr,"%s: dx_bondi = %"GSYM" pc (%"GSYM" in code units) is greater than cell width = %e pc.
                // Don't deposit refinement zone.\n",
				//__FUNCTION__, dx_bondi_pc, dx_bondi, dx_pc);
                continue;
			}
			for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel){
                // fprintf(stderr,"%s: BondiRadius/factor = %e pc is less than cell width = %e pc. Deposit refinement zone.\n",
                // 	__FUNCTION__, dx_bondi_pc, dx_pc);
                // SG Deposit refinement zone with dx_bondi
                if (Temp->GridData->DepositRefinementZone(level,pos,dx_bondi) == FAIL) {
                    ENZO_FAIL("Error in grid->DepositRefinementZone.\n")
                } // end IF
            } // end FOR
        } // end ELSE (BH)
    } // end FOR over particles
    return SUCCESS;
} // END SetFlaggingField


int ActiveParticleType_SmartStar::SmartStarParticleFeedback(int nParticles,
    ActiveParticleList<ActiveParticleType>& ParticleList, FLOAT dx, 
	LevelHierarchyEntry *LevelArray[], int ThisLevel, int SmartStarID)
{
    /* Skip if we're not on the maximum refinement level.
     This should only ever happen right after creation and then
     only in pathological cases where creation is happening at 
     the edges of two regions at the maximum refinement level- not doing this for SG. */

    /* For each particle, loop over all the grids and do feedback
     if the grid overlaps with the feedback zone */
    FLOAT Time = LevelArray[ThisLevel]->GridData->ReturnTime();
    float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits;
    double MassUnits;
    GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
             &TimeUnits, &VelocityUnits, Time);
    MassUnits = DensityUnits * POW(LengthUnits,3);
    int NumberOfGrids;
    HierarchyEntry **Grids = NULL;
    NumberOfGrids = GenerateGridArray(LevelArray, ThisLevel, &Grids);
    FLOAT dx_sg = LevelArray[ThisLevel]->GridData->CellWidth[0][0]; // SG. Grid cell width.

    for (int i = 0; i < nParticles; i++) {
        if (SmartStarFeedback == FALSE){
            continue;
        }

        ActiveParticleType_SmartStar* SS;
        SS = static_cast<ActiveParticleType_SmartStar*>(ParticleList[i]);
	    int pclass = SS->ParticleClass;
        float OldMassCheck = SS->oldmass; // SG. Checking.
        FLOAT AccretionRadius = SS->AccretionRadius;

        // SG. BH class.
        if(pclass == BH){
            grid* FeedbackZone = ConstructFeedbackZone(ParticleList[i], FLOAT(AccretionRadius/dx_sg), dx_sg,
                                                       Grids, NumberOfGrids, ALL_FIELDS);
            if (MyProcessorNumber == FeedbackZone->ReturnProcessorNumber()) {
                if (FeedbackZone->ApplySmartStarParticleFeedback(&ParticleList[i]) == FAIL)
                    return FAIL;
            }

            /* If a PISN just went off then we can safely delete the particle. */
            if(ParticleList[i]->ShouldDelete() == true) {
                printf("%s: Delete SS %d following PISN\n", __FUNCTION__,
                       static_cast<ActiveParticleType_SmartStar*>(ParticleList[i])->ReturnID());
                fflush(stdout);
                ParticleList[i]->DisableParticle(LevelArray,FeedbackZone->ReturnProcessorNumber());
                printf("%s: SS %d deleted\n", __FUNCTION__,
                       static_cast<ActiveParticleType_SmartStar*>(ParticleList[i])->ReturnID());
                fflush(stdout);
            }
	        ActiveParticleFindAll(LevelArray, &nParticles, SmartStarID,
                                  ParticleList);
	        DistributeFeedbackZone(FeedbackZone, Grids, NumberOfGrids, ALL_FIELDS);
	        delete FeedbackZone;
        } // SG. End BH class condition.
	    else if (pclass == POPIII){ // SG. Add POPIII class condition
            grid* FeedbackZone = ConstructFeedbackZone(ParticleList[i], FLOAT(AccretionRadius/dx_sg), dx_sg,
                                                       Grids, NumberOfGrids, ALL_FIELDS);

		    // SG. Set to 0 before it's calculated by owning proc and then communicated with other procs in
            // CommunicateAllSumValues().
		    FLOAT positions[3] = {0,0,0}; // SG. All elements initialised to zero.
		    FLOAT NewAccretionRadius = 0;
		    FLOAT* pos;

		    if (MyProcessorNumber == FeedbackZone->ReturnProcessorNumber()) {
                if (FeedbackZone->ApplySmartStarParticleFeedback(&ParticleList[i]) == FAIL)
                    return FAIL;

                // SG. positions is the array of dereferenced particle positions in each dim.
                positions[0] = SS->ReturnPosition()[0];
                positions[1] = SS->ReturnPosition()[1];
                positions[2] = SS->ReturnPosition()[2];
            } // END my processor

            // SG. Communicate with all procs the updated accretion radius and particle position.
            CommunicationAllSumValues(&NewAccretionRadius, 1);
            CommunicationAllSumValues(positions, 3);

            SS->AccretionRadius = NewAccretionRadius;
            AccretionRadius = SS->AccretionRadius;

            SS->pos[0] = positions[0];
            SS->pos[1] = positions[1];
            SS->pos[2] = positions[2];

            DistributeFeedbackZone(FeedbackZone, Grids, NumberOfGrids, ALL_FIELDS);
            delete FeedbackZone;

        } // SG. End POPIII class condition.
    } // SG. End FOR loop over particles
  
    if (AssignActiveParticlesToGrids(ParticleList, nParticles, LevelArray) == FAIL)
        return FAIL;

    delete [] Grids;
    return SUCCESS;
} // SG. End SmartStarParticleFeedback function.


/* 
 * This function can be used to reset the particle acceleration if required.
 * For example if a massless particle needs to be fixed in space. 
 * See ActiveParticle_RadiationParticle.C for details. 
 */
int ActiveParticleType_SmartStar::ResetAcceleration(float *ActiveParticleAcceleration)
{
  return SUCCESS;
}


/*
 * For brute force creation of a particle. Useful for converting from 
 * star objects to active particles.
 */
int ActiveParticleType_SmartStar::CreateParticle(grid *thisgrid_orig, ActiveParticleFormationData &supp_data,
                                                 int particle_index)
{
  return SUCCESS;
} 

bool ActiveParticleType_SmartStar::IsARadiationSource(FLOAT Time)
{
    if(BH == ParticleClass)
        return (SmartStarBHRadiativeFeedback == TRUE) ? true : false;
    if(POPIII == ParticleClass || SMS == ParticleClass)
        return (SmartStarStellarRadiativeFeedback == TRUE) ? true : false;
}


/*
 * Calculate the Jeans Length of a gas cell and return in cgs 
*/

static double JeansLength(float T, float dens, float density_units)
{
  float jeans_length = 0.0;

  jeans_length = 15*kboltz*T/(4.0*M_PI*GravConst*mh*dens*density_units);
  return sqrt(jeans_length);
}


/*
 * The accretion radius is updated as mass is accreted. 
 * The accretion radius is set to the gravitational radius of the star.
 */
static void UpdateAccretionRadius(ActiveParticleType*  ThisParticle, float newmass,
				  FLOAT OldAccretionRadius, float avgtemp,
				  float mass_units, float length_units)
{
    fprintf(stderr,"%s: got here.\n", __FUNCTION__); // SG. Debug comment.
    ActiveParticleType_SmartStar* SS;
    SS = static_cast<ActiveParticleType_SmartStar*>(ThisParticle);

    // SG. Finding dx and MassConversion.
    grid* APGrid = ThisParticle->ReturnCurrentGrid();
    double dx = APGrid->GetCellWidth(0,0);
    double dx_pc = dx*length_units/pc_cm;   //in pc
    float MassConversion = (float) (dx*dx*dx * double(mass_units));  //convert to g
    MassConversion = MassConversion/SolarMass; // convert to Msun

    float Temperature = avgtemp;
    FLOAT soundspeed2 = kboltz*Temperature*Gamma/(Mu*mh);
    // SG. NewAccretionRadius = Bondi Radius. Is OldAccretionRadius = 4 cell widths?
    FLOAT NewAccretionRadius = (2 * GravConst * newmass * MassConversion / soundspeed2) / length_units; //in code_length
    fprintf(stderr,"%s: NewAccretionRadius = %e and OldAccretionRadius = %e.\n", __FUNCTION__,
            NewAccretionRadius, OldAccretionRadius); // SG. Debug comment.
    SS->AccretionRadius = NewAccretionRadius;
    return;
}


// SG.JR. Take dx out of argument of function. Calculate within function.
int ActiveParticleType_SmartStar::UpdateAccretionRateStats(int nParticles,
				ActiveParticleList<ActiveParticleType>& ParticleList,
				LevelHierarchyEntry *LevelArray[], int ThisLevel)
{
    /* Get units, define current time */
    FLOAT Time = LevelArray[ThisLevel]->GridData->ReturnTime();
    float ctime = LevelArray[ThisLevel]->GridData->ReturnTime();
    float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits;
    double MassUnits;
    GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
             &TimeUnits, &VelocityUnits, Time);
    MassUnits = DensityUnits * POW(LengthUnits,3);

    /* SG. Test cell width received from GridData vs APGrid. */
    double dx_grid = LevelArray[ThisLevel]->GridData->CellWidth[0][0];
    double dx_pc1 = dx_grid*LengthUnits/pc_cm;   //in pc
    double MassConversion = (double) (dx_grid*dx_grid*dx_grid * double(MassUnits));  //convert to g
    MassConversion = MassConversion/SolarMass; // convert to Msun

    /* SG. Moved mass conversion to within loop over particles. */
    for (int i = 0; i < nParticles; i++) {
        grid* APGrid = ParticleList[i]->ReturnCurrentGrid();
        ActiveParticleType_SmartStar* SS;
        SS = static_cast<ActiveParticleType_SmartStar*>(ParticleList[i]);
        int MyLevel = SS->ReturnLevel();

        /* If not on APGrid level or mass is 0, continue */
        if (ThisLevel != MyLevel){
        return SUCCESS;
        }

        if (SS->Mass == 0){
            fprintf(stderr,"%s: SS Mass is zero. TimeIndex not incremented. Continue.\n", __FUNCTION__);
            continue;
        }

        /* Define cell width and mass conversion based on it */
        double dx = APGrid->CellWidth[0][0];
        double dx_pc = dx*LengthUnits/pc_cm;   //in pc
        double MassConversion = (double) (dx_grid*dx_grid*dx_grid * double(MassUnits));  //SG. Changed to dx_grid. convert to g
        MassConversion = MassConversion/SolarMass; // convert to Msun

        /* Only update stats if on correct processor */
        if (MyProcessorNumber == APGrid->ReturnProcessorNumber()) {
            ActiveParticleType_SmartStar* SS;
            SS = static_cast<ActiveParticleType_SmartStar*>(ParticleList[i]);

            #if SSDEBUG
            fprintf(stderr,"%s: deltatime = %f years\t TIMEGAP = %0.2f years\n",
            __FUNCTION__, (ctime - SS->AccretionRateTime[SS->TimeIndex])*TimeUnits/yr_s,
            (float)TIMEGAP);
            #endif

            /* SG. Changing TIMEGAP to 100 years as in the print statement just above. */
            if( ((ctime - SS->AccretionRateTime[SS->TimeIndex])*TimeUnits/yr_s > (float)TIMEGAP)
                || (SS->TimeIndex == 0)) {

                fprintf(stderr,"%s: level = %"ISYM" and MyLevel = %"ISYM" and NoParticles = %"ISYM".\n", __FUNCTION__,
                        ThisLevel, MyLevel, nParticles);
                fprintf(stderr,"%s: cell width = %e pc (APGrid) on level = %"ISYM".\n", __FUNCTION__, dx_pc, ThisLevel);

                float omass = SS->ReturnOldMass();
                float cmass = SS->ReturnMass();

                if(omass < -1e-10) {
                    omass = 0.0;
                    fprintf(stderr, "%s: omass less than zero. Updating to omass = %e Msun.\n",
                            __FUNCTION__, omass*MassConversion);
                }

                /* Increment time index */
                SS->TimeIndex++;

                /* timeindex and original timeindex */
                int timeindex = (SS->TimeIndex)%NTIMES;
                int otimeindex = timeindex - 1;
                if(otimeindex == -1) //loop back
                    otimeindex = NTIMES -1;
                float otime = SS->AccretionRateTime[otimeindex];
                if(otime == -1.0) {
                    otimeindex = 0;
                    otime = SS->AccretionRateTime[otimeindex];
                }

                /* Update SS accrate and timeindex */
                float deltatime = ctime - otime;
                float accrate = (cmass - omass)/deltatime;
                float Age = ctime - SS->BirthTime;
                SS->AccretionRate[timeindex] = accrate*dx*dx*dx;
                SS->AccretionRateTime[timeindex] = ctime;
                SS->TimeIndex = timeindex;

                /* Prints */
                fprintf(stderr, "old_mass = %e Msolar\t cmass = (%e code) %e Msolar\n", omass*MassConversion,
                    cmass, cmass*MassConversion);
                fprintf(stderr, "accrate = %1.2e Msolar/yr\t deltatime = %3.3f Myrs\t TimeIndex = %d\t "
                                "Particle Mass = %1.2e Msolar\t Age = %1.3f Myr\t Lifetime = %1.2f Myr\t Class = %d\n",
                    (SS->AccretionRate[timeindex]*MassUnits/TimeUnits)*yr_s/SolarMass,
                    deltatime*TimeUnits/Myr_s,
                    SS->TimeIndex,
                    SS->ReturnMass()*MassConversion,
                    Age*TimeUnits/Myr_s,
                    SS->RadiationLifetime*TimeUnits/Myr_s,
                    SS->ParticleClass);
                    /* End Prints */

                /* Set omass to cmass */
                omass = cmass;
                SS->oldmass = cmass;

                /* Change of particle class to SMS if PopIII and accretion rate is high.
                        Disabled by SG/BS for mass removal from sphere method */

                if(SS->ParticleClass > SMS){
                    continue;
                    } else {
                    if(dx_pc < SMS_RESOLUTION) {
                        /*
                        * Using the time-averaged accretion rates determine if the
                        * SMS is accreting fast enough or
                        * if it is falling onto the main sequence.
                        * This can also allow a POPIII star to change into a SMS
                        */
                        if((SS->AccretionRate[timeindex]*MassUnits/TimeUnits)*yr_s/SolarMass > CRITICAL_ACCRETION_RATE) {
                            if(SS->ParticleClass == POPIII) {
                                // SG. Disabled Particle Class switching from POPIII TO SMS.
                                printf("%s: UPDATE SG: NO CHANGE of ParticleClass.\n", __FUNCTION__);
                                }
                                // SS->ParticleClass = SMS; SG.
                                }
                                else {
                                    float Age = Time - SS->BirthTime;
                                    if(Age*TimeUnits/yr_s > 1e4 && SMS == SS->ParticleClass) {
                                        /* Don't do this at very start */
                                        printf("%s: WARNING: ParticleClass switching from SMS to "
                                               "POPIII (deltatime = %f kyrs)\n", __FUNCTION__,
                                               deltatime*TimeUnits/(yr_s*1e3));
                                        printf("%s: WARNING: Accretion Rate = %f Msolar/yr. Critical rate = "
                                               "%f Msolar/yr\n", __FUNCTION__,
                                        (SS->AccretionRate[SS->TimeIndex]*MassUnits/TimeUnits)*yr_s/SolarMass,
                                        CRITICAL_ACCRETION_RATE);
                                        SS->ParticleClass = POPIII;
                                        } // END if age & SMS
                                } // End else
                        } // End SMS_RES
                    } // END else
                } // End IF timegap exceeded
        } // End IF processor
    } // End FOR loop over particles

  return SUCCESS;
}

/* 
 * Update stellar lifetimes at each step to account for accretion and mergers
 */
int ActiveParticleType_SmartStar::UpdateRadiationLifetimes(int nParticles,
				  ActiveParticleList<ActiveParticleType>& ParticleList,
						LevelHierarchyEntry *LevelArray[],
						int ThisLevel)
{
    /* SG. Skip stellar accretion even in high-res cases. */
#if STELLAR_ACCRETION_OFF
    return SUCCESS;
#endif

  FLOAT Time = LevelArray[ThisLevel]->GridData->ReturnTime();
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,VelocityUnits;
  double MassUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
           &TimeUnits, &VelocityUnits, Time);
  MassUnits = DensityUnits * POW(LengthUnits,3);
  for (int i = 0; i < nParticles; i++) {
      grid* APGrid = ParticleList[i]->ReturnCurrentGrid();
      // SG. Finding dx and MassConversion.
      double dx = APGrid->CellWidth[0][0];
      double dx_pc = dx*LengthUnits/pc_cm;   //in pc
      double MassConversion = (double) (dx*dx*dx * double(MassUnits));  //convert to g
      MassConversion = MassConversion/SolarMass; // convert to Msun
      if (MyProcessorNumber == APGrid->ReturnProcessorNumber()) {
          ActiveParticleType_SmartStar* SS;
          SS = static_cast<ActiveParticleType_SmartStar*>(ParticleList[i]);
          if(POPIII == SS->ParticleClass) {
              double StellarMass = SS->Mass * MassConversion; //Msolar
              float logm = log10((float) StellarMass);
              // First in years, then convert to code units
              SS->RadiationLifetime = POW(10.0, (9.785 - 3.759 * logm + 1.413 * logm * logm -
                                                 0.186 * logm * logm * logm)) / (TimeUnits / yr_s);
              SS->StellarAge = SS->RadiationLifetime; //update stellar age too
          }
      }
  }
  return SUCCESS;
}


/* Find the stellar radius of the sink particle using the SAM from Smith et al. (2011) */
static float GetStellarRadius(float cmass, float accrate)
{
    float p1 = 0, p2 = 0;
    float A1 = 1.0, A2 = 1.0;
    float stellar_radius = 0.0;
    float R_ms = 0.28*POW(cmass, 0.61);
    p1 = 5.0*POW(accrate, 0.27);   //in mass units
    p2 = 7.0*POW(accrate, 0.27);   //in mass units

    if(cmass <= p1) {
        stellar_radius = 26.0*POW(cmass, 0.27)*POW(accrate/1e-3, 0.41);
    }
    else if(cmass > p1 && cmass < p2) {
        stellar_radius = A1*POW(cmass, 3.0);
    }
    else
        stellar_radius = max(A2*POW(cmass, -2.0), R_ms);
    return stellar_radius; //in solar radii
}

namespace {
  ActiveParticleType_info *SmartStarInfo =
    register_ptype <ActiveParticleType_SmartStar>
    ("SmartStar");
}

std::vector<ParticleAttributeHandler*>
  ActiveParticleType_SmartStar::AttributeHandlers;
