/***********************************************************************
/
/ Smart Star Particle
/ Author: John Regan
/         Based on previous active particle code from many authors
/ Date: Early 2017
/
************************************************************************/

#include "ActiveParticle_SmartStar.h"
#include "phys_constants.h"
#define SSDEBUG 1
#define SSDEBUG_TOTALMASS 0

#define DYNAMIC_ACCRETION_RADIUS 0
#define BONDIHOYLERADIUS 0
#define MINIMUMPOTENTIAL 1
#define CALCDIRECTPOTENTIAL 0
#define JEANSREFINEMENT  1
#define MASSTHRESHOLDCHECK 0
#define JEANSLENGTHCALC    1
#define MASSTHRESHOLD      0.1                       //Msolar in grid
#define COOLING_TIME       0
#define NUMSSPARTICLETYPES 4
int DetermineSEDParameters(ActiveParticleType_SmartStar *SS,FLOAT Time, FLOAT dx);


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
//float ActiveParticleType_SmartStar::RadiationLifetime = FLOAT_UNDEFINED;
int ActiveParticleType_SmartStar::FeedbackDistRadius = INT_UNDEFINED;
int ActiveParticleType_SmartStar::FeedbackDistTotalCells = INT_UNDEFINED;
int ActiveParticleType_SmartStar::FeedbackDistCellStep = INT_UNDEFINED;
static double JeansLength(float T, float dens, float density_units);
static void UpdateAccretionRadius(ActiveParticleType*  ThisParticle, float newmass,
				  FLOAT AccretionRadius, FLOAT dx, float avgtemp,
				  float mass_units, float length_units);
static float GetStellarRadius(float cmass, float accrate);
int StarParticlePopIII_IMFInitialize(void);
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
  ah.push_back(new Handler<ap, float, &ap::oldmass>("oldmass"));
  //ah.push_back(new ArrayHandler<ap, float, 3, &ap::acc>("particle_acceleration_x", 0));
  //ah.push_back(new ArrayHandler<ap, float, 3, &ap::acc>("particle_acceleration_y", 1));
  //ah.push_back(new ArrayHandler<ap, float, 3, &ap::acc>("particle_acceleration_z", 2));
  
  /* Initialize the IMF lookup table if requested and not defined */

  if (PopIIIInitialMassFunction)
    StarParticlePopIII_IMFInitialize();
  
  printf("SmartStar Initialisation complete\n"); fflush(stdout);
  return SUCCESS;
}

int ActiveParticleType_SmartStar::EvaluateFormation
(grid *thisgrid_orig, ActiveParticleFormationData &data)
{
  //Particle Lifetimes
  float SmartStarLifetime[4];
  SmartStarLifetime[0] = 2e6; //in yrs
  SmartStarLifetime[1] = 2e6; //in yrs
  SmartStarLifetime[2] = 1e20; //in yrs
  SmartStarLifetime[3] = 2e7; //in yrs
  
  // No need to do the rest if we're not on the maximum refinement level.
  if (data.level != MaximumRefinementLevel)
    return SUCCESS;

  SmartStarGrid *thisGrid =
    static_cast<SmartStarGrid *>(thisgrid_orig);

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
  //printf("PotentialField = %p\n", PotentialField);
  
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
#if SSDEBUG
  fprintf(stdout, "%s: Right half thinking of creating a SmartSink particle\n", __FUNCTION__); fflush(stdout);
#endif
#if MASSTHRESHOLDCHECK
  JeansMass = thisGrid->CalculateJeansMass(data.DensNum, data.Temperature, data.DensityUnits);  //In Msolar
#endif

  for (k = thisGrid->GridStartIndex[2]; k <= thisGrid->GridEndIndex[2]; k++) {
    for (j = thisGrid->GridStartIndex[1]; j <= thisGrid->GridEndIndex[1]; j++) {
      for (i = thisGrid->GridStartIndex[0]; i <= thisGrid->GridEndIndex[0]; i++) {
	index = GRIDINDEX_NOGHOST(i, j, k);
	DensityThreshold = ActiveParticleDensityThreshold*mh/data.DensityUnits; 
	// If no more room for particles, throw an ENZO_FAIL
	if (data.NumberOfNewParticles >= data.MaxNumberOfNewParticles)
	  return FAIL;
	// 1. Check we're on the maximum refinement level - already done at start of function
	
	// 2. Does this cell violate the Jeans condition or overdensity threshold
	// Fedderath condition #1
#if JEANSREFINEMENT
	if (JeansRefinement) {
	  CellTemperature = (JeansRefinementColdTemperature > 0) ? JeansRefinementColdTemperature : data.Temperature[index];
	  int JeansFactor = 8; //RefineByJeansLengthSafetyFactor
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

	if (density[index] <= DensityThreshold)
	  continue;

#if SSDEBUG
	printf("Time = %e yr\n", thisGrid->ReturnTime()*data.TimeUnits/yr_s);
	printf("Density = %g\t DensityThreshold = %g\n", density[index]*data.DensityUnits/mh, DensityThreshold*data.DensityUnits/mh);
	printf("JeansDensity = %g\t APThreshold = %g\n", JeansDensity*data.DensityUnits/mh, ActiveParticleDensityThreshold);
#endif
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
        fprintf(stdout, "%s: Negative Divergence passed\n", __FUNCTION__);
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
	  GravitationalMinimum  = thisGrid->FindMinimumPotential(centralpos,
				  JLength*RefineByJeansLengthSafetyFactor*4.0,
				  PotentialField);
	  if(PotentialField[index] > GravitationalMinimum) {
#if SSDEBUG
	    printf("FAILURE: GravitationalMinimum = %g\t "		\
		   "PotentialField[index] = %g\n\n", GravitationalMinimum, PotentialField[index]);
#endif
	    continue;
	  }
#if SSDEBUG
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


	printf("%s: Forming SS particle out of cellindex %d\n", __FUNCTION__, index); fflush(stdout);
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
	  * We want to avoid spurious SF in a minihalo that is being heated (e.g. by LW or dynamical 
	  * heating). Therefore we insist on the mass accretion criteria. If the mass accretion rate 
	  * is low we don't allow the SMS pathway to trigger spurious SF. 
	  */
	int stellar_type = -99999;
	if(HasMetalField &&    data.TotalMetals[index] > PopIIIMetalCriticalFraction) {
	  stellar_type = POPII;
	  //	  float PopIIMass = StarClusterFormEfficiency*ExtraDensity*dx*dx*dx*ConverttoSolar;
	  // if(PopIIMass < StarClusterMinimumMass)
	  //  continue;
	}
	else if(data.H2Fraction[index] >  PopIIIH2CriticalFraction) {
	  stellar_type = POPIII;
	}
	else if(accrate*3.154e7*ConverttoSolar/data.TimeUnits > 1e-2) {
	  stellar_type = SMS;
	  printf("!!!!!!!!SMS Formed\t accrate = %e Msolar/yr",
		 accrate*3.154e7*ConverttoSolar/data.TimeUnits);
	}
	
	if(stellar_type < 0)
	  continue;
	
	ActiveParticleType_SmartStar *np = new ActiveParticleType_SmartStar();
	data.NumberOfNewParticles++;
	data.NewParticles.insert(*np);

	

	np->type = np->GetEnabledParticleID();
	/* 
	 * Set Bithtime to -1.0 initially. If everything is ok and the 
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

	printf("pos = %f %f %f\n", np->pos[0], np->pos[1], np->pos[2]);
	printf("GridLeftEdge = %f %f %f\n", thisGrid->CellLeftEdge[0][0],
	       thisGrid->CellLeftEdge[1][0], thisGrid->CellLeftEdge[2][0]);
	float *apvel = new float[MAX_DIMENSION];
	/* Calculate AP velocity by taking average of 
	 * surrounding 125 cells 
	 */
	apvel = thisGrid->AveragedVelocityAtCell3D(index, data.DensNum, data.Vel1Num);
	np->vel[0] = apvel[0];
	np->vel[1] = apvel[1];
	np->vel[2] = apvel[2];
	
	
	if (HasMetalField)
	  np->Metallicity = data.TotalMetals[index];
	else
	  np->Metallicity = 0.0;

	np->TimeIndex = 0; //Start at 0 - we'll increment at the start of the update function.
	np->AccretionRadius = dx*ACCRETIONRADIUS;
	for(int acc = 1; acc < NTIMES; acc++) {
	  np->AccretionRate[acc] = -11111.0;
	  np->AccretionRateTime[acc] = -11111.0;
	}
	
	/* Calculate initial accretion rate onto cell */
	np->AccretionRate[0] = accrate;
	np->AccretionRateTime[0] = np->BirthTime;
	
	np->RadiationLifetime= SmartStarLifetime[np->ParticleClass]*yr_s/data.TimeUnits;
	np->NotEjectedMass = 0.0;
       
	/* The mass of the particles formed depends on the resolution and is handled 
	 * in a call to `RemoveMassFromGridAfterFormation` which is called from 
	 * `AfterEvolveLevel`. np->Mass is initally set here but can get overwritten 
	 *  RemoveMassFromGridAfterFormation. The cell values are not updated here but 
	 * instead are updated in `RemoveMassFromGridAfterFormation`.
	 * We set the initial mass to zero so that merging works ok i.e. nothing spurious 
	 * happens in this case. 
	 */
	np->Mass = 0.0;
	//printf("%s: Total Mass in Accretion Region = %g Msolar (Threshold = %g)\n", __FUNCTION__,
	//     TotalMass*ConverttoSolar, (double)MASSTHRESHOLD);
      } // i
    } // j
  } // k
  if(data.NumberOfNewParticles > 0) {
    printf("Particles (%d) Created and done in %s\n", data.NumberOfNewParticles, __FUNCTION__);
  }
#if CALCDIRECTPOTENTIAL
  delete [] PotentialField;
#endif
  return SUCCESS;
}

/*
 * This routine handles the stellar feedback including the 
 * accretion luminosity and supernova feedback. 
 * The jet feedback and thermal feedback are
 * calculated in SmartStarParticleFeedback and called from AfterEvolveLevel
 * The radiative feedback is called from BeforeEvolveLevel
 */
int ActiveParticleType_SmartStar::EvaluateFeedback(grid *thisgrid_orig,
						   ActiveParticleFormationData &data)
{

  /* Feedback not handled here */
  return SUCCESS; 
  
  SmartStarGrid *thisGrid =
    static_cast<SmartStarGrid *>(thisgrid_orig);
 
  float *density = thisGrid->BaryonField[data.DensNum];
  float *totalenergy = thisGrid->BaryonField[data.TENum];
  float *gasenergy = thisGrid->BaryonField[data.GENum];
  FLOAT Time = thisGrid->ReturnTime();
  float dt = thisGrid->dtFixed;
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
    VelocityUnits;
  double MassUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, Time);
  int npart = thisGrid->NumberOfActiveParticles;
  FLOAT dx = thisGrid->CellWidth[0][0];
  const double LumConvert =   POW(LengthUnits,2)/POW(TimeUnits,3);
  MassUnits = DensityUnits * POW(LengthUnits,3);
  double MassConversion =  ((double)dx*dx*dx * MassUnits);  //convert to g
  /* 
   * This routine needs to be called here. Do not return before 
   * this point!
   */
  for (int n = 0; n < npart; n++) {
    ActiveParticleType_SmartStar *ThisParticle =
      static_cast<ActiveParticleType_SmartStar*>(
              thisGrid->ActiveParticles[n]);

    ThisParticle->level = data.level;
  }

  if(SmartStarFeedback == FALSE)
    return SUCCESS;
  for (int n = 0; n < npart; n++) {
    ActiveParticleType_SmartStar *ThisParticle =
      static_cast<ActiveParticleType_SmartStar*>(
              thisGrid->ActiveParticles[n]);
    
    if(BH == ThisParticle->ParticleClass)
      continue;
    /* The first form of feedback is accretion luminosity 
     * This is given by L = G M_acc M_* / R_*
     */
    /* StarMass in units of Msolar */
    float StarMass = ThisParticle->Mass * MassConversion / SolarMass;
    /* Accretion Rate in units of Msolar/s */
    float AccretionRate = ThisParticle->AccretionRate[ThisParticle->TimeIndex]*MassUnits/(SolarMass*TimeUnits);
    if(AccretionRate*yr_s < 1e-30) /* Takes care of any negative accretion rates */
      AccretionRate = 1e-30/yr_s;
    /* Star Radius in units of Rsolar */
    float StarRadius = GetStellarRadius(StarMass, AccretionRate);
    float LThisTimestep=0.0, sn_nrg_thistimestep = 0.0;
    if((POPIII == ThisParticle->ParticleClass) || (SMS == ThisParticle->ParticleClass)) {
      /* Calculate SpecificL [ergs s^-1 g^-1] */
      float SpecificL = GravConst * AccretionRate * SolarMass / (StarRadius * SolarRadius);
      SpecificL = SpecificL/LumConvert;  /* Convert (specific) Luminosity to code units */
#if SSDEBUG
      printf("%s: dx = %e\t MassConversion = %e\n", __FUNCTION__, dx, MassConversion);
      printf("%s: AccretionRate = %e Msolar/yr (Code = %e)\n", __FUNCTION__,
	     AccretionRate*yr_s, ThisParticle->AccretionRate[ThisParticle->TimeIndex]);
      printf("%s: Mass = %e Msolar\n", __FUNCTION__, StarMass);
      printf("%s: Radius = %e Rsolar\n", __FUNCTION__, StarRadius);
      printf("%s: SpecificL Total = %e Lsolar\n", __FUNCTION__, SpecificL*LumConvert * StarMass * SolarMass/SolarLuminosity);
#endif
    /* Now calculate the energy over the timestep */
      LThisTimestep = SpecificL*dt/float(StarFeedbackDistTotalCells);
#if SSDEBUG
      printf("%s: StarFeedbackDistTotalCells = %f\n", __FUNCTION__,
	     float(StarFeedbackDistTotalCells));
      printf("%s: dt = %e yrs\n", __FUNCTION__,
	     dt*TimeUnits/yr_s);
      printf("%s: L Total = %e Lsolar [ergs]\n", __FUNCTION__,
	     LThisTimestep*(StarMass * SolarMass));
#endif   
      
      float Age = Time - ThisParticle->BirthTime;
#if SSDEBUG      
      printf("%s: Star Age = %e yrs\n", __FUNCTION__, Age*TimeUnits/yr_s);
      printf("%s: Radiation Lifetime =  %e yrs\n", __FUNCTION__, ThisParticle->RadiationLifetime*TimeUnits/yr_s);
      printf("%s: Particle Class = %d\t POPIII = %d\n", __FUNCTION__, ThisParticle->ParticleClass, POPIII);
#endif
      /* Check for star death and transition to BH */
      if(ThisParticle->RadiationLifetime < Age) {/* Star needs to go supernovae and change type */
#ifdef SNEFEEDBACK
	/* Based on the stellar mass calculate the SN energy and ejecta mass.
       Fitted from Heger & Woosley (2002) */
	if(POPIII == ThisParticle->ParticleClass) {
	  float he_core = (13.0/24.0) * (StarMass - 20.0);
	  sn_nrg_thistimestep = (5.0 + 1.304 * (he_core - 64.0)) * 1e51; //ergs
	  printf("%s: Supernova Energy = %e [ergs]\n", __FUNCTION__, sn_nrg_thistimestep);
	  fflush(stdout);
	  //code energy (specific units)
	  //sn_nrg_thistimestep /= (TimeUnits*LumConvert*ThisParticle->Mass; //specific energy
	  sn_nrg_thistimestep /= (TimeUnits*LumConvert*StarMass*SolarMass);
	}
#endif
	printf("%s: !!!!!!!!!!!!!!!Transition from %d particle to BH (Particle Type = %d)\n", 
	       __FUNCTION__, ThisParticle->ParticleClass, BH); fflush(stdout);
	printf("%s: Star Age = %e yrs\n", __FUNCTION__, Age*TimeUnits/yr_s);
	printf("%s: Radiation Lifetime =  %e yrs\n", __FUNCTION__, ThisParticle->RadiationLifetime*TimeUnits/yr_s);
	ThisParticle->ParticleClass = BH;
      }
    }

    // Calculate 3D grid indices
    FLOAT xstart = thisGrid->CellLeftEdge[0][0];
    FLOAT ystart = thisGrid->CellLeftEdge[1][0];
    FLOAT zstart = thisGrid->CellLeftEdge[2][0];
    int GridXSize = thisGrid->GridDimension[0];
    int GridYSize = thisGrid->GridDimension[1];
    int GridZSize = thisGrid->GridDimension[2];
    int GridDimension[3] = {thisGrid->GridDimension[0],
			    thisGrid->GridDimension[1],
			    thisGrid->GridDimension[2]};
    int i = int((ThisParticle->pos[0] - xstart)/thisGrid->CellWidth[0][0]);
    int j = int((ThisParticle->pos[1] - ystart)/thisGrid->CellWidth[1][0]);
    int k = int((ThisParticle->pos[2] - zstart)/thisGrid->CellWidth[2][0]);
    // Check bounds - if star particle is outside of this grid then give a warning and continue
    //printf("%s: GridDimension = %d %d %d\n", __FUNCTION__, GridDimension[0], GridDimension[1], GridDimension[2]);
    if (i < 0 || i > GridXSize-1 || j < 0 || j > GridYSize-1 || k < 0 || k > GridZSize-1){
      fprintf(stdout, "Particle out of grid; xind, yind, zind = %"ISYM", %"ISYM", %"ISYM"\n",i,j,k);
      continue;
    }
    // Calculate serial index

    int index = GRIDINDEX_NOGHOST(i,j,k);
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
    // Add energy to the energy field
    for (int kc = k - FeedbackDistRadius; kc <= k + FeedbackDistRadius; kc++){
      int stepk = ABS(kc - k);
      for (int jc = j - FeedbackDistRadius; jc <= j + FeedbackDistRadius; jc++){
	int stepj = stepk + ABS(jc - j);
	for (int ic = i - FeedbackDistRadius; ic <= i + FeedbackDistRadius; ic++){
	  int cellstep = stepj + ABS(ic - i);
	  int DistIndex = GRIDINDEX_NOGHOST(ic,jc,kc);
	  if (cellstep <= FeedbackDistCellStep) {
	    float energybefore = totalenergy[DistIndex];
	    totalenergy[DistIndex] = ((totalenergy[DistIndex]*density[DistIndex]) + 
				      LThisTimestep + sn_nrg_thistimestep)/density[DistIndex];
	    //printf("%s: Accretion Energy added = %e [code energy]\n", __FUNCTION__, LThisTimestep*(StarMass*SolarMass));
	    //printf("%s: Supernova Energy added = %e [code energy]\n", __FUNCTION__, sn_nrg_thistimestep*(StarMass*SolarMass));
	    //printf("%s: Fractional Energy increase = %e\n", __FUNCTION__,
	    // (totalenergy[DistIndex] - energybefore)/energybefore);
	    if (DualEnergyFormalism == 1) {
	      float energybefore = gasenergy[DistIndex];
	      gasenergy[DistIndex] = ((gasenergy[DistIndex]*density[DistIndex]) + 
				      LThisTimestep)/density[DistIndex];
	      //printf("%s: Fractional Gas Energy increase = %e\n", __FUNCTION__,
	      //   (gasenergy[DistIndex] - energybefore)/energybefore);
	    }
	  }
	}
      }
    }

  }
  
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
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
    VelocityUnits;
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
    ActiveParticleFindAll(LevelArray, &nParticles, SmartStarID, 
    SmartStarList);

#ifdef TRANSFER
  if (CallEvolvePhotons) {
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
	MassConversion = (double) (dx*dx*dx * mfactor); //Converts to Solar Masses
	source = ThisParticle->RadiationSourceInitialize();
	double PMass = ThisParticle->Mass*MassConversion;
	float ramptime = 0.0;
	if(POPIII == ThisParticle->ParticleClass ||
	   SMS == ThisParticle->ParticleClass) {
	  ramptime = yr_s * 1e4 / TimeUnits;
	}
	if(POPII == ThisParticle->ParticleClass) {
	  ramptime = yr_s * StarClusterMinDynamicalTime / TimeUnits;
	}
	/* Call Function to return SED parameters */
	if(ThisParticle->DetermineSEDParameters(Time, dx) == FAIL)
	  return FAIL;

	printf("%s: Mass = %e\t RadiationLifetime = %e yrs\n", __FUNCTION__, PMass,
	       ThisParticle->RadiationLifetime*TimeUnits/yr_s); fflush(stdout);

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
	//	printf("%s: source->Lifetime = %f years\n", __FUNCTION__, source->LifeTime*TimeUnits/yr_s);
	//	printf("%s: Energy = %f %f %f %f %f\n", __FUNCTION__, source->Energy[0],
	//       source->Energy[1], source->Energy[2], source->Energy[3], source->Energy[4]);
	//printf("%s: SED = %f %f %f %f %f\n", __FUNCTION__, source->SED[0],
	//       source->SED[1], source->SED[2], source->SED[3], source->SED[4]);
#if SSDEBUG
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
       exit(-99);
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


int ActiveParticleType_SmartStar::RemoveMassFromGridAfterFormation(int nParticles, 
    ActiveParticleList<ActiveParticleType>& ParticleList,
    LevelHierarchyEntry *LevelArray[], int ThisLevel)
{
  //Particle Lifetimes
  float SmartStarLifetime[4];
  SmartStarLifetime[0] = 2e6; //in yrs
  SmartStarLifetime[1] = 2e6; //in yrs
  SmartStarLifetime[2] = 1e20; //in yrs
  SmartStarLifetime[3] = 2e7; //in yrs
  
  float StellarMasstoRemove = 0.0, CellDensityAfterFormation = 0.0;
  /* Skip accretion if we're not on the maximum refinement level.
     This should only ever happen right after creation and then
     only in pathological cases where sink creation is happening at
     the edges of two regions at the maximum refinement level */

  if (ThisLevel < MaximumRefinementLevel)
    return SUCCESS;
  FLOAT Time = LevelArray[ThisLevel]->GridData->ReturnTime();
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
    VelocityUnits;
  double MassUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, Time);
  MassUnits = DensityUnits * POW(LengthUnits,3);

  float tdyn_code = StarClusterMinDynamicalTime/(TimeUnits/yr_s);
  for (int i = 0; i < nParticles; i++) {
   grid* APGrid = ParticleList[i]->ReturnCurrentGrid();
   
   if (MyProcessorNumber == APGrid->ReturnProcessorNumber()) {
     ActiveParticleType_SmartStar* SS;
     SS = static_cast<ActiveParticleType_SmartStar*>(ParticleList[i]);

     /*
      * Only interested in newly formed particles
      */
     if(SS->TimeIndex != 0)
       continue;

      FLOAT dx = APGrid->CellWidth[0][0];
      float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1,
	TimeUnits = 1, VelocityUnits = 1,
	PressureUnits = 0, GEUnits = 0, VelUnits = 0;
      double MassUnits = 1, CellVolume = 1;
      if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
		   &TimeUnits, &VelocityUnits, Time) == FAIL) {
        ENZO_FAIL("Error in GetUnits.");
      }
      MassUnits = DensityUnits * POW(LengthUnits,3);
      int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
      if (APGrid->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL)
	{
	  ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
	}
      int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
	DINum, DIINum, HDINum;
      if (MultiSpecies) 
	if (APGrid->IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, 
					HeIIINum, HMNum, H2INum, H2IINum, DINum, 
					DIINum, HDINum) == FAIL) {
	  ENZO_FAIL("Error in grid->IdentifySpeciesFields.");
	}
      FLOAT dx_pc = dx*LengthUnits/pc_cm;   //in pc
      float *density = APGrid->BaryonField[DensNum];
      int cellindex_x = (SS->pos[0] - APGrid->CellLeftEdge[0][0])/dx,
	cellindex_y = (SS->pos[1] - APGrid->CellLeftEdge[1][0])/dx,
	cellindex_z = (SS->pos[2] - APGrid->CellLeftEdge[2][0])/dx;


      printf("SS->Pos = %f %f %f\n", SS->pos[0], SS->pos[1], SS->pos[2]);
      printf("GridLeftEdge = %f %f %f\n", APGrid->CellLeftEdge[0][0], APGrid->CellLeftEdge[1][0],
	     APGrid->CellLeftEdge[2][0]);
      int cellindex = APGrid->GetIndex(cellindex_x, cellindex_y, cellindex_z);
      float DensityThreshold = ActiveParticleDensityThreshold*mh/DensityUnits; 
     
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
	int JeansFactor = 8; //RefineByJeansLengthSafetyFactor
	float JeansDensityUnitConversion = (Gamma*pi*kboltz) / (Mu*mh*GravConst);
	float JeansDensity = JeansDensityUnitConversion * 1.01 * CellTemperature /
	  POW(LengthUnits*dx*JeansFactor,2);
	JeansDensity /= DensityUnits;
	DensityThreshold = min(DensityThreshold,JeansDensity);
      }
#endif

      float ParticleDensity = density[cellindex] - DensityThreshold;
      if(ParticleDensity < 0.0) {
	 printf("%s: cellindex = %d\n", __FUNCTION__, cellindex);
	 printf("density[cellindex] = %e cm^-3\n", density[cellindex]*DensityUnits/mh);
	 printf("DensityThreshold = %e cm^-3\n", DensityThreshold*DensityUnits/mh);
	 printf("SS->ParticleClass = %d\n", SS->ParticleClass); fflush(stdout);
	 ENZO_FAIL("Particle Density is negative. Oh dear.\n");
      }
     
      float newcelldensity = density[cellindex] - ParticleDensity;
      if(SMS == SS->ParticleClass) {

	 if(dx_pc <= SMS_RESOLUTION) { /* Accrete as normal - just remove mass from the cell */
	   density[cellindex] = newcelldensity;
	   SS->BirthTime = APGrid->ReturnTime();
	   SS->Mass = ParticleDensity;
	   SS->oldmass = SS->Mass;
	   return SUCCESS;
	 }
       }
       
       else if(POPIII == SS->ParticleClass) {
	 if(dx_pc <= POPIII_RESOLUTION) { /* Accrete as normal - just remove mass from the cell */
	   density[cellindex] = newcelldensity;
	   SS->BirthTime = APGrid->ReturnTime();
	   SS->Mass = ParticleDensity;
	   SS->oldmass = SS->Mass;
	   return SUCCESS;
	 }
       }
       else if(POPII == SS->ParticleClass) {
	 /* 
	  * For PopII stars we do this if the mass exceeds the minimum mass
	  */
	 float PopIIMass = SS->Mass*dx*dx*dx*MassUnits/SolarMass;
	 if(PopIIMass > StarClusterMinimumMass) {
	   density[cellindex] = (1 - StarClusterFormEfficiency)*density[cellindex];
	   SS->BirthTime = APGrid->ReturnTime();
	   SS->Mass = StarClusterFormEfficiency*density[cellindex];
	   SS->oldmass = SS->Mass;
	   return SUCCESS;
	 }
       }   

        /* 
	 * If the formation mass is below a resolution based threshold
	 * Remove mass from grid and replace by a uniform density sphere which accounts for the 
	 * subgrid ionisation that has taken place and accounts for the mass that should have been 
	 * accreted.
	 */
       
       /***********************************************************************

            For star formation, we need to find a sphere with enough mass to
            accrete.  We step out by a cell width when searching.

       ***********************************************************************/
      printf("%s: Low resolution run.....\n", __FUNCTION__);
      FLOAT Radius = 0.0;
      int feedback_flag = -99999;
      float MassEnclosed = 0;
      float Metallicity2 = 0;
      float Metallicity3 = 0;
      float ColdGasMass = 0;
      float AvgVelocity[MAX_DIMENSION];
      for (int dim = 0; dim < MAX_DIMENSION; dim++)
	AvgVelocity[dim] = 0.0;
      bool SphereTooSmall = true;
      float ShellMass, ShellMetallicity2, ShellMetallicity3, ShellColdGasMass, 
	ShellVelocity[MAX_DIMENSION];
      while (SphereTooSmall) { 
	Radius += APGrid->CellWidth[0][0];
	bool IsSphereContained = SS->SphereContained(LevelArray, ThisLevel, Radius);

	if (IsSphereContained == false)
	  break;
	ShellMass = 0;
	ShellMetallicity2 = 0;
	ShellMetallicity3 = 0;
	ShellColdGasMass = 0;
	for (int dim = 0; dim < MAX_DIMENSION; dim++)
	  ShellVelocity[dim] = 0.0;

	bool MarkedSubgrids = false;
	LevelHierarchyEntry *Temp = NULL;
	HierarchyEntry *Temp2 = NULL;

	for (int l = ThisLevel; l < MAX_DEPTH_OF_HIERARCHY; l++) {
	  Temp = LevelArray[l];
	  while (Temp != NULL) {
	    
	    /* Zero under subgrid field */
	    
	    if (!MarkedSubgrids) {
	      Temp->GridData->
		ZeroSolutionUnderSubgrid(NULL, ZERO_UNDER_SUBGRID_FIELD);
	      Temp2 = Temp->GridHierarchyEntry->NextGridNextLevel;
	      while (Temp2 != NULL) {
		Temp->GridData->ZeroSolutionUnderSubgrid(Temp2->GridData, 
							 ZERO_UNDER_SUBGRID_FIELD);
		Temp2 = Temp2->NextGridThisLevel;
	      }
	    } // ENDIF !MarkedSubgrids
	    /* Sum enclosed mass in this grid. Mass is in Msolar*/
	    
	    Temp->GridData->GetEnclosedMassInShell(SS->pos, Radius-APGrid->CellWidth[0][0], Radius, 
						   ShellMass, ShellMetallicity2, 
						   ShellMetallicity3,
						   ShellColdGasMass, ShellVelocity,
						   -1);
	    
	    Temp = Temp->NextGridThisLevel;
	    
	  } // END: Grids
	  
	} // END: level
	MarkedSubgrids = true;
	//float values[7];
	//values[0] = ShellMetallicity2;
	//values[1] = ShellMetallicity3;
	//values[2] = ShellMass;
	//values[3] = ShellColdGasMass;
	//for (int dim = 0; dim < MAX_DIMENSION; dim++)
	//  values[4+dim] = ShellVelocity[dim];
	
	//LCAPERF_START("star_FindFeedbackSphere_Sum");
	//CommunicationAllSumValues(values, 7);
	//LCAPERF_STOP("star_FindFeedbackSphere_Sum");
	
	//ShellMetallicity2 = values[0];
	//ShellMetallicity3 = values[1];
	//ShellMass = values[2];
	//ShellColdGasMass = values[3];
	//for (int dim = 0; dim < MAX_DIMENSION; dim++)
	// ShellVelocity[dim] = values[4+dim];
	
	MassEnclosed += ShellMass;
	ColdGasMass += ShellColdGasMass;
	// Must first make mass-weighted, then add shell mass-weighted
	// (already done in GetEnclosedMassInShell) velocity and
	// metallicity.  We divide out the mass after checking if mass is
	// non-zero.
	Metallicity2 = Metallicity2 * (MassEnclosed - ShellMass) + ShellMetallicity2;
	Metallicity3 = Metallicity3 * (MassEnclosed - ShellMass) + ShellMetallicity3;
	for (int dim = 0; dim < MAX_DIMENSION; dim++)
	  AvgVelocity[dim] = AvgVelocity[dim] * (MassEnclosed - ShellMass) +
	    ShellVelocity[dim];
	printf("MassEnclosed = %e Msolar\n", MassEnclosed);
	if (MassEnclosed == 0) {
	  IsSphereContained = false;
	  return SUCCESS;
	}
	
	Metallicity2 /= MassEnclosed;
	Metallicity3 /= MassEnclosed;
	for (int dim = 0; dim < MAX_DIMENSION; dim++)
	  AvgVelocity[dim] /= MassEnclosed;
	
	
	/* Now remove mass based on star particle type */
	
	if(POPIII == SS->ParticleClass) {
	  SS->AssignMassFromIMF();

	  SphereTooSmall = MassEnclosed < (2*SS->Mass);
	  // to make the total mass PopIIIStarMass
	  StellarMasstoRemove = SS->Mass;  // [Msolar]
	
	}
	else if(SMS == SS->ParticleClass) {
	  /* I think we need a heavy seed IMF here */
	  SS->Mass = 1000.0; //hardcoding to 1000 Msolar
	  SS->RadiationLifetime =  SmartStarLifetime[SS->ParticleClass]*yr_s/TimeUnits;
	  SphereTooSmall = MassEnclosed < (2*SS->Mass);
	  StellarMasstoRemove = SS->Mass; //[Msolar]
	}
	else if(POPII == SS->ParticleClass) {
	  float AvgDensity = (float) 
	    (double(SolarMass * MassEnclosed) / 
	     double(4*pi/3.0 * pow(Radius*LengthUnits, 3))); /* cgs density */
	  float DynamicalTime = sqrt((3.0 * pi) / (32.0 * GravConst * AvgDensity)) /
	    TimeUnits;
	  float ColdGasFraction = ColdGasMass / MassEnclosed;
	  StellarMasstoRemove = ColdGasFraction * StarClusterFormEfficiency * MassEnclosed; //[Msolar]
	  SphereTooSmall = DynamicalTime < tdyn_code;
	  
	}
	// Remove the stellar mass from the sphere and distribute the
	// gas evenly in the sphere since this is what will happen once
	// the I-front passes through it.
	
	CellDensityAfterFormation = (float) 
	  (double(SolarMass * (MassEnclosed - StellarMasstoRemove)) / 
	   double(4.0*pi/3.0 * POW(Radius*LengthUnits, 3)) /
	   DensityUnits); /* converted to code density */
	
      }  /* end while(SphereTooSmall) */
#ifdef NOT_NECESSARY
       /* Don't allow the sphere to be too large (2x leeway) */
       const float epsMass = 9.0;
       float eps_tdyn;
       if ((PopIII == SS->ParticleClass || SMS == SS->ParticleClass) && LevelArray[level+1] != NULL) {
	 if (MassEnclosed > (1.0+epsMass)*(StellarMasstoRemove)) {
	   SphereContained = FALSE;
	   return SUCCESS;
	 }
       }
       else if (PopII == SS->ParticleClass && LevelArray[level+1] != NULL) {
	 eps_tdyn = sqrt(1.0+epsMass) * tdyn_code;
	 if (DynamicalTime > eps_tdyn) {
	   SphereContained = FALSE;
	   return SUCCESS;
	 }
       }
#endif
       printf("%s: Update Grid Densities (i.e. remove mass)\n", __FUNCTION__);
       /* The Radius of influence is set by the sphere over which we had to 
	* loop to find sufficient enclosed mass. 
	*/
       float InfluenceRadius = Radius;
       /* Update cell information */
       int index = 0;
       FLOAT delx = 0.0, dely = 0.0, delz = 0.0, radius2 = 0.0, DomainWidth[MAX_DIMENSION];
       int MetallicityField = FALSE;
       int SNColourNum, MetalNum, Metal2Num, MBHColourNum, Galaxy1ColourNum, 
	 Galaxy2ColourNum, MetalIaNum, MetalIINum;
       
       if (APGrid->IdentifyColourFields(SNColourNum, Metal2Num, MetalIaNum, 
				      MetalIINum, MBHColourNum, Galaxy1ColourNum, 
				      Galaxy2ColourNum) == FAIL)
	 ENZO_FAIL("Error in grid->IdentifyColourFields.\n");

       MetalNum = max(Metal2Num, SNColourNum);
       MetallicityField = (MetalNum > 0) ? TRUE : FALSE;
       float metallicity = 0.0;
       for (int dim = 0; dim < APGrid->GridRank; dim++)
	  DomainWidth[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];
       for (int k = 0; k < APGrid->GridDimension[2]; k++) {

	 delz = APGrid->CellLeftEdge[2][k] + 0.5*APGrid->CellWidth[2][k] - SS->pos[2];
	 int sz = sign(delz);
	 delz = fabs(delz);
	 delz = min(delz, DomainWidth[2]-delz);
	 
	 for (int j = 0; j < APGrid->GridDimension[1]; j++) {

	   dely = APGrid->CellLeftEdge[1][j] + 0.5*APGrid->CellWidth[1][j] - SS->pos[1];
	   int sy = sign(dely);
	   dely = fabs(dely);
	   dely = min(dely, DomainWidth[1]-dely);
	   
	   for (int i = 0; i < APGrid->GridDimension[0]; i++, index++) {
	     float ionizedFraction = 0.999;  // Assume an initial HII region
	     delx = APGrid->CellLeftEdge[0][i] + 0.5*APGrid->CellWidth[0][i] - SS->pos[0];
	     int sx = sign(delx);
	     delx = fabs(delx);
	     delx = min(delx, DomainWidth[0]-delx);
	     
	     radius2 = delx*delx + dely*dely + delz*delz;
	     if (radius2 <= InfluenceRadius*InfluenceRadius) {

	       radius2 = max(radius2, 0.0625*APGrid->CellWidth[0][i]*APGrid->CellWidth[0][i]); // (0.25*dx)^2
	       
	       if (MetallicityField == TRUE)
		 metallicity = APGrid->BaryonField[MetalNum][index] / APGrid->BaryonField[DensNum][index];
	       else
		 metallicity = 0;
	       float fh = CoolData.HydrogenFractionByMass;
	       float fhz = fh * (1-metallicity);
	       float fhez = (1-fh) * (1-metallicity);
	       SS->BirthTime = APGrid->ReturnTime();
	       double cellvolume = APGrid->CellWidth[0][i]*APGrid->CellWidth[1][j]*
		 APGrid->CellWidth[2][k];
	       SS->Mass = (StellarMasstoRemove*SolarMass/MassUnits)/cellvolume; //code density
	       SS->oldmass = SS->Mass;
	       APGrid->BaryonField[DensNum][index] = CellDensityAfterFormation;
	       if (MultiSpecies) {
		 APGrid->BaryonField[DeNum][index] = APGrid->BaryonField[DensNum][index] * ionizedFraction;
		 APGrid->BaryonField[HINum][index] = APGrid->BaryonField[DensNum][index] * (1-ionizedFraction) * fhz;
		 APGrid->BaryonField[HIINum][index] = APGrid->BaryonField[DensNum][index] * ionizedFraction * fhz;
		 APGrid->BaryonField[HeINum][index] = APGrid->BaryonField[DensNum][index] * (1-ionizedFraction) * fhez;
		 APGrid->BaryonField[HeIINum][index] = APGrid->BaryonField[DensNum][index] * (ionizedFraction) * fhez;
		 APGrid->BaryonField[HeIIINum][index] = 1e-10 * APGrid->BaryonField[DensNum][index];
	       }
	       if (MultiSpecies > 1) {
		 APGrid->BaryonField[HMNum][index] = tiny_number;
		 APGrid->BaryonField[H2INum][index] = tiny_number;
		 APGrid->BaryonField[H2IINum][index] = tiny_number;
	       }
	       if (MultiSpecies > 2) {
		 APGrid->BaryonField[DINum][index] = APGrid->BaryonField[DensNum][index] * fh *
		   CoolData.DeuteriumToHydrogenRatio * (1-ionizedFraction);
		 APGrid->BaryonField[DIINum][index] = APGrid->BaryonField[DensNum][index] * fh *
		   CoolData.DeuteriumToHydrogenRatio * ionizedFraction;
		 APGrid->BaryonField[HDINum][index] = tiny_number * APGrid->BaryonField[DensNum][index];
	       }
	       
	     }  // END if inside radius
	     
	   }  // END i-direction
	 }  // END j-direction
       }  // END k-direction
       

       
   } /*This Processor */
 } /* End loop over APs */
 return SUCCESS;
}
int ActiveParticleType_SmartStar::Accrete(int nParticles, 
    ActiveParticleList<ActiveParticleType>& ParticleList,
    FLOAT AccretionRadius, FLOAT dx,
    LevelHierarchyEntry *LevelArray[], int ThisLevel)
{

  /* Skip accretion if we're not on the maximum refinement level.
     This should only ever happen right after creation and then
     only in pathological cases where sink creation is happening at
     the edges of two regions at the maximum refinement level */

  if (ThisLevel < MaximumRefinementLevel)
    return SUCCESS;
  FLOAT Time = LevelArray[ThisLevel]->GridData->ReturnTime();
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
    VelocityUnits;
  double MassUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, Time);
  MassUnits = DensityUnits * POW(LengthUnits,3);
  double MassConversion = (double) (dx*dx*dx * double(MassUnits));  //convert to g
  /* For each particle, loop over all of the grids and do accretion
     if the grid overlaps with the accretion zone     
  */
  
  int i, NumberOfGrids;
  int *FeedbackRadius = NULL;
  HierarchyEntry **Grids = NULL;
  grid *sinkGrid = NULL;

  bool SinkIsOnThisProc, SinkIsOnThisGrid;

  float SubtractedMass, SubtractedMomentum[3] = {};

  NumberOfGrids = GenerateGridArray(LevelArray, ThisLevel, &Grids);

  for (i = 0; i < nParticles; i++) {


    /* 
     * Accretion is only allowed if it makes sense:
     * 1. Black Holes in general are always allowed to accrete.
     * 2. SMS can accrete if there is sufficient resolution to do so 
     * 3. PopIII stars can accrete if there is sufficient resolution
     * 4. PopII stars never accrete
     */

    
    float MassInSolar = ParticleList[i]->ReturnMass()*MassConversion/SolarMass;
    AccretionRadius =  static_cast<ActiveParticleType_SmartStar*>(ParticleList[i])->AccretionRadius;
    int pclass = static_cast<ActiveParticleType_SmartStar*>(ParticleList[i])->ParticleClass;
    FLOAT dx_pc = dx*LengthUnits/pc_cm;   //in pc
    if(pclass == POPIII) {
      /* 
       * We only accrete onto POPIII stars if our maximum 
       * spatial resolution is better than 1e-3 pc
       */
      if(dx_pc > POPIII_RESOLUTION) //we don't have sufficient resolution
	continue;
    }
    else if(pclass == SMS) {
       /* 
       * We only accrete onto SMSs if our maximum 
       * spatial resolution is better than 1e-1 pc
       */
      if(dx_pc > SMS_RESOLUTION) //we don't have sufficient resolution
	continue;
      
    }
    else if(pclass == BH) {
      /* We always accrete onto BHs */
      ;
    }
    else if(pclass == POPII) {
      /* We never accrete onto POPII stars */
      continue;

    }

    printf("Accreting: Particle Class %d\n", pclass);
    grid* FeedbackZone = ConstructFeedbackZone(ParticleList[i], int(AccretionRadius/dx),
					       dx, Grids, NumberOfGrids, ALL_FIELDS);
    grid* APGrid = ParticleList[i]->ReturnCurrentGrid();
    if (MyProcessorNumber == FeedbackZone->ReturnProcessorNumber()) {

      float AccretionRate = 0;

      if (FeedbackZone->AccreteOntoSmartStarParticle(ParticleList[i],
			      AccretionRadius, &AccretionRate) == FAIL)
	return FAIL;

      FLOAT *pos = ParticleList[i]->ReturnPosition();


#if BONDIHOYLERADIUS
      /* Check what the Bondi-Hoyle radius - we should accrete out to that if required */
      float mparticle = ParticleList[i]->ReturnMass()*dx*dx*dx;
      float *vparticle = new float[3];
      vparticle = ParticleList[i]->ReturnVelocity();
      int size = APGrid->GetGridSize();
      float *Temperature = new float[size]();
 
      APGrid->ComputeTemperatureField(Temperature);
      FLOAT BondiHoyleRadius = APGrid->CalculateBondiHoyleRadius(mparticle, vparticle, Temperature);
      if(static_cast<ActiveParticleType_SmartStar*>(ParticleList[i])->AccretionRadius < BondiHoyleRadius) {
	static_cast<ActiveParticleType_SmartStar*>(ParticleList[i])->AccretionRadius = BondiHoyleRadius;
	printf("%s: Updating accretion radius to Bondi-Hoyle radius = %e pc (%f cells)\n", __FUNCTION__,
	       static_cast<ActiveParticleType_SmartStar*>(ParticleList[i])->AccretionRadius*LengthUnits/pc,
	       static_cast<ActiveParticleType_SmartStar*>(ParticleList[i])->AccretionRadius/dx);
      }
#endif
	// No need to communicate the accretion rate to the other CPUs since this particle is already local.
	/* Need to decide how often I update the accretion history */
    
    }
    DistributeFeedbackZone(FeedbackZone, Grids, NumberOfGrids, ALL_FIELDS);

    delete FeedbackZone;
  }

  if (AssignActiveParticlesToGrids(ParticleList, nParticles, LevelArray) == FAIL)
    return FAIL;

  delete [] Grids;
  return SUCCESS;
}

int ActiveParticleType_SmartStar::SetFlaggingField(
    LevelHierarchyEntry *LevelArray[], int level,
    int TopGridDims[], int SmartStarID)
{
  /* Generate a list of all sink particles in the simulation box */
  int i, nParticles;
  FLOAT *pos = NULL;
  ActiveParticleList<ActiveParticleType> SmartStarList;
  LevelHierarchyEntry *Temp = NULL;

  ActiveParticleFindAll(LevelArray, &nParticles, SmartStarID, 
      SmartStarList);

  for (i=0 ; i<nParticles; i++){
    pos = SmartStarList[i]->ReturnPosition();
    FLOAT accrad = static_cast<ActiveParticleType_SmartStar*>(SmartStarList[i])->AccretionRadius;
    for (Temp = LevelArray[level]; Temp; Temp = Temp->NextGridThisLevel)
      if (Temp->GridData->DepositRefinementZone(level,pos,accrad) == FAIL) {
	ENZO_FAIL("Error in grid->DepositRefinementZone.\n")
	  }
  }

  return SUCCESS;
}




int ActiveParticleType_SmartStar::SmartStarParticleFeedback(int nParticles,
    ActiveParticleList<ActiveParticleType>& ParticleList, FLOAT dx, 
	LevelHierarchyEntry *LevelArray[], int ThisLevel)
{
  
  /* Skip if we're not on the maximum refinement level. 
     This should only ever happen right after creation and then
     only in pathological cases where creation is happening at 
     the edges of two regions at the maximum refinement level */

  if (ThisLevel < MaximumRefinementLevel || SmartStarFeedback == FALSE)
    return SUCCESS;

  /* For each particle, loop over all of the grids and do feedback 
     if the grid overlaps with the feedback zone                   */
  
  int i, NumberOfGrids;
  HierarchyEntry **Grids = NULL;
  
  NumberOfGrids = GenerateGridArray(LevelArray, ThisLevel, &Grids);
  
 
  
  for (i = 0; i < nParticles; i++) {
    FLOAT AccretionRadius =  static_cast<ActiveParticleType_SmartStar*>(ParticleList[i])->AccretionRadius;
    grid* FeedbackZone = ConstructFeedbackZone(ParticleList[i], int(AccretionRadius/dx), dx, 
					       Grids, NumberOfGrids, ALL_FIELDS);
    if (MyProcessorNumber == FeedbackZone->ReturnProcessorNumber()) {
      if (FeedbackZone->ApplySmartStarParticleFeedback(&ParticleList[i]) == FAIL)
	return FAIL;
    }

    DistributeFeedbackZone(FeedbackZone, Grids, NumberOfGrids, ALL_FIELDS);

    delete FeedbackZone;
  }

  delete [] Grids;
  return SUCCESS;
}

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
int ActiveParticleType_SmartStar::CreateParticle(grid *thisgrid_orig,
							 ActiveParticleFormationData &supp_data,
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
 * The accretion radius is to the gravitational radius of the star.
 */
static void UpdateAccretionRadius(ActiveParticleType*  ThisParticle, float newmass,
				  FLOAT OldAccretionRadius, FLOAT dx, float avgtemp,
				  float mass_units, float length_units)
{
  ActiveParticleType_SmartStar* SS;
  SS = static_cast<ActiveParticleType_SmartStar*>(ThisParticle);

  float Temperature = avgtemp;
  FLOAT soundspeed2 = kboltz*Temperature*Gamma/(Mu*mh); 
  float MassConversion = (float) (dx*dx*dx * double(mass_units));  //convert to g
  FLOAT NewAccretionRadius = (2 * GravConst * newmass * MassConversion / soundspeed2) / length_units; //in code_length
  if(NewAccretionRadius > OldAccretionRadius)
    SS->AccretionRadius = NewAccretionRadius;
  return;
}


int ActiveParticleType_SmartStar::UpdateAccretionRateStats(int nParticles,
				ActiveParticleList<ActiveParticleType>& ParticleList,
				FLOAT dx,
				LevelHierarchyEntry *LevelArray[], int ThisLevel)
{

  FLOAT Time = LevelArray[ThisLevel]->GridData->ReturnTime();
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
    VelocityUnits;
  double MassUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, Time);
  MassUnits = DensityUnits * POW(LengthUnits,3);
  double MassConversion = (double) (dx*dx*dx * double(MassUnits));  //convert to g
  MassConversion = MassConversion/SolarMass;
  float ctime = LevelArray[ThisLevel]->GridData->ReturnTime();
  FLOAT dx_pc = dx*LengthUnits/pc_cm;   //in pc
  for (int i = 0; i < nParticles; i++) {
    
    grid* APGrid = ParticleList[i]->ReturnCurrentGrid();
    
    if (MyProcessorNumber == APGrid->ReturnProcessorNumber()) {
      ActiveParticleType_SmartStar* SS;
      SS = static_cast<ActiveParticleType_SmartStar*>(ParticleList[i]);

#if SSDEBUG
      printf("%s: deltatime = %f years\t TIMEGAP = %0.2f years\n",
	     __FUNCTION__, (ctime - SS->AccretionRateTime[SS->TimeIndex])*TimeUnits/yr_s, 
	     (float)TIMEGAP);
#endif
      //We should update when the time between stored rates exceeds TIMEGAP
      if(ctime - SS->AccretionRateTime[SS->TimeIndex] > (TIMEGAP*APGrid->ReturnTimeStep())) {
	float omass = SS->oldmass;
	float cmass = ParticleList[i]->ReturnMass();
	if(cmass - omass < 0.0) { //Can happen after a restart due to rounding
	  break;
	}
	SS->TimeIndex++;
	
	int timeindex = (SS->TimeIndex)%NTIMES;
	int otimeindex = timeindex - 1;
	if(otimeindex == -1) //loop back
	  otimeindex = NTIMES -1;
	float otime = SS->AccretionRateTime[otimeindex];
	if(otime == -1.0) {
	  otimeindex = 0;
	  otime = SS->AccretionRateTime[otimeindex]; 
	}

	float deltatime = ctime - otime;
	float accrate = (cmass - omass)/deltatime;

	SS->AccretionRate[timeindex] = accrate*dx*dx*dx;
	SS->AccretionRateTime[timeindex] = ctime;
	SS->oldmass = cmass;
	SS->TimeIndex = timeindex;
	fprintf(stdout, "old_mass = %e Msolar\t cmass = %e Msolar\n", omass*MassConversion,
		cmass*MassConversion);
	fprintf(stdout, "accrate = %e Msolar/yr\t accratetime = %e yrs \t deltatime = %f yrs\t index = %d\t Particle Mass = %e Msolar\t Class = %d\n",
		(SS->AccretionRate[timeindex]*MassUnits/TimeUnits)*yr_s/SolarMass,
		(SS->AccretionRateTime[timeindex]*TimeUnits)/yr_s,
		deltatime*TimeUnits/yr_s,
		SS->TimeIndex,
		SS->ReturnMass()*MassConversion,
		SS->ParticleClass);
	if(SS->ParticleClass < BH) {
	  /*
	   * Using the time-averaged accretion rates determine if the SMS is accreting fast enough or
	   * if it is falling onto the main sequence.
	   */
	  if((SS->AccretionRate[SS->TimeIndex]*MassUnits/TimeUnits)*yr_s/SolarMass
	     > CRITICAL_ACCRETION_RATE) {
	    SS->ParticleClass = SMS;
	  }
	  else {
	    float Age = Time - SS->BirthTime;
	    if(Age*TimeUnits/yr_s > 1e4 && SMS == SS->ParticleClass) { /* Don't do this at very start */
	      printf("%s: WARNING: ParticleClass switching from SMS to POPIII (deltatime = %f kyrs)\n", __FUNCTION__,
		     deltatime*TimeUnits/(yr_s*1e3));
	      printf("%s: WARNING: Accretion Rate = %f Msolar/yr. Critical rate = %f Msolar/yr\n", __FUNCTION__,
		     (SS->AccretionRate[SS->TimeIndex]*MassUnits/TimeUnits)*yr_s/SolarMass,
		     CRITICAL_ACCRETION_RATE);
	      SS->ParticleClass = POPIII;
	    }
	  }
	}
      }	       
    }
  }
#if SSDEBUG
  printf("Done in %s\n", __FUNCTION__);
#endif
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
