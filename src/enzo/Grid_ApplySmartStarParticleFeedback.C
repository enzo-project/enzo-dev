/***********************************************************************
/
/  Apply feedback to the temporary grid for the Smart Star Particles
/ The feedback methods dealt with here are the Thermal and Kinetic modes
/
/  written by: John Regan
/  date:       December, 2017
/
/  note: Based on methods originally implemented by Stephen Skory
/ 
/  description: the functions here apply the impact of i) supernova
/  explosions, ii) the feedback from black hole thermal energy and 
/  iii) the feedback from black hole jets. 
/
/  Radiative feedback is not dealt with here. See ActiveParticle_SmartStar.C instead. 
************************************************************************/

#include "preincludes.h"

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "units.h"
#include "Fluxes.h"
#include "GridList.h"
#include "phys_constants.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "ActiveParticle.h"
#include "phys_constants.h"
#include "ActiveParticle_SmartStar.h"

#define SSFEED_DEBUG 1
#define MAX_TEMPERATURE 1e8
#define RAMPTIME 100000.0 //yrs
#define DENSITY_WEIGHTED 1
#define SINE_WAVE       0
#define IMPOSETHRESHOLD 1
#define THRESHOLDFRACTION 1   //Solarmasses ejected per jet event
#define OPENING_ANGLE pi/360.0  //pi/3.9
#define HW_BH_MASS 1   // SG. Heger-Woosley mass used for BH remanant for POPIII particles.

int search_lower_bound(float *arr, float value, int low, int high, 
		       int total);

int grid::ApplySmartStarParticleFeedback(ActiveParticleType** ThisParticle){
	
  /* Return if this doesn't involve us */
  if (MyProcessorNumber != ProcessorNumber){
		 return SUCCESS;
		} 
    
  if(SmartStarFeedback == FALSE){
   return SUCCESS;
	 }

  ActiveParticleType_SmartStar *SS = static_cast<ActiveParticleType_SmartStar*>(* ThisParticle);
  const float PISNLowerMass = 140, PISNUpperMass = 260;
  const float TypeIILowerMass = 11, TypeIIUpperMass = 40.1;
  
  // From Nomoto et al. (2006)
  const float HypernovaMetals[] = {3.36, 3.53, 5.48, 7.03, 8.59}; // SolarMass
  const float HypernovaEnergy[] = {10, 10, 20, 25, 30}; // 1e51 erg 
  const float CoreCollapseMetals[] = {3.63, 4.41, 6.71, 8.95, 11.19}; // SolarMass
  const float CoreCollapseEnergy[] = {1, 1, 1, 1, 1}; // 1e51 erg

  const float SNExplosionMass[] = {19.99, 25, 30, 35, 40.01};  // SolarMass
  const float *SNExplosionMetals = (PopIIIUseHypernova ==TRUE) ? 
    HypernovaMetals : CoreCollapseMetals;
  const float *SNExplosionEnergy = (PopIIIUseHypernova ==TRUE) ? 
    HypernovaEnergy : CoreCollapseEnergy;
  
  /* Check whether the cube that circumscribes the accretion zone intersects with this grid */

   FLOAT *pos = SS->ReturnPosition();
   FLOAT rad = SS->AccretionRadius*0.5;
   
   if ((GridLeftEdge[0] > pos[0]+rad) || (GridRightEdge[0] < pos[0]-rad) ||
       (GridLeftEdge[1] > pos[1]+rad) || (GridRightEdge[1] < pos[1]-rad) ||
       (GridLeftEdge[2] > pos[2]+rad) || (GridRightEdge[2] < pos[2]-rad))
     return SUCCESS;

  
  float dx = float(this->CellWidth[0][0]);
  FLOAT dV = POW(dx, 3.0);
  float dt = float(this->ReturnTimeStep());
  /* Set the units. */
 
  float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1,
    TimeUnits = 1, VelocityUnits = 1,
    PressureUnits = 0, GEUnits = 0, VelUnits = 0;
  double MassUnits = 1.0;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, this->ReturnTime()) == FAIL) {
        ENZO_FAIL("Error in GetUnits.");
  }
  MassUnits = DensityUnits * POW(LengthUnits,3);
  double MassConversion = (double) (dx*dx*dx * MassUnits);  //convert to g      
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
         Vel3Num, TENum) == FAIL) {
     ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
   }

  /* Find Multi-species fields. */

  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
  if (MultiSpecies) 
    if (this->IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, 
				    HeIIINum, HMNum, H2INum, H2IINum, DINum, 
				    DIINum, HDINum) == FAIL) {
        ENZO_FAIL("Error in grid->IdentifySpeciesFields.");
    }

  /***********************************************************************
                                SUPERNOVAE
  ************************************************************************/

  // Assume that the remnant is still in the free expansion stage and
  // hasn't had any radiative losses.  In this case, the ejecta will
  // be at 3/4 the radius of the shock front (see Ostriker & McKee
  // 1988 or Tenorio-Tagle 1996).
  
  float ionizedFraction = 0.999;  // Assume supernova is ionized
  FLOAT Time = this->ReturnTime();
  float Age = Time - SS->BirthTime;



  /***********************************************************************
                                SMS
  ************************************************************************/

  /* 
   * SMS don't go supernova they just directly collapse into BHs (of the same mass)
   */
  if(SS->ParticleClass == SMS) {
    if(Age > SS->RadiationLifetime) {/* SG. Two options now. Core-collapse SN option for 40Msun SMS. */
      double StellarMass = SS->Mass*MassConversion/SolarMass; /* In Msolar */
      printf("%s: StellarMass = %lf\n", __FUNCTION__, StellarMass);
      
      if (StellarMass >= TypeIILowerMass && StellarMass <= (TypeIIUpperMass + 1)) { 
	printf("%s: Star going Supernova!\n", __FUNCTION__);
	printf("%s: Age = %1.2f Myr\t RadiationLifetime = %1.2f Myr\n", __FUNCTION__,
	     Age*TimeUnits/Myr_s, SS->RadiationLifetime*TimeUnits/Myr_s);
	double SNEnergy, HeliumCoreMass, Delta_SF, MetalMass;
	FLOAT Radius = PopIIISupernovaRadius * pc_cm / LengthUnits;
	FLOAT StarLevelCellWidth = this->CellWidth[0][0];
	Radius = max(Radius, 3.5*StarLevelCellWidth);
	FLOAT EjectaVolume = 4.0/3.0 * pi * pow(Radius*LengthUnits, 3);
	FLOAT EjectaDensity = StellarMass * SolarMass / EjectaVolume / DensityUnits;
	FLOAT EjectaMetalDensity = 0.0, EjectaThermalEnergy = 0.0;
	// Core-collapse
	int bin = search_lower_bound((float*)SNExplosionMass, StellarMass, 0, 5, 5);
	float frac = (SNExplosionMass[bin+1] - StellarMass) / 
	  (SNExplosionMass[bin+1] - SNExplosionMass[bin]);
	SNEnergy = 1e51 * (SNExplosionEnergy[bin] + 
			   frac * (SNExplosionEnergy[bin+1] - SNExplosionEnergy[bin]));
	MetalMass = (SNExplosionMetals[bin] + frac * (SNExplosionMetals[bin+1] - SNExplosionMetals[bin]));
	// Heger-Woosley (2002) relation for BHMass
	HeliumCoreMass = (13./24.) * (StellarMass - 20);
	StellarMass = HeliumCoreMass; //msun
	SS->Mass = StellarMass*SolarMass/MassConversion; //code density
	
	SS->ParticleClass = BH;
	SS->StellarAge = SS->RadiationLifetime; //Record last stellar age
	SS->RadiationLifetime = 1e20;
	printf("%s: Post-SNe (SMS): ParticleClass now %d\t Lifetime = %f Myr\n", __FUNCTION__,
	       SS->ParticleClass, SS->RadiationLifetime*TimeUnits/Myr_s);
	EjectaMetalDensity = MetalMass * SolarMass / EjectaVolume / DensityUnits;
	EjectaThermalEnergy = SNEnergy / (StellarMass * SolarMass) / VelocityUnits /
	  VelocityUnits;
	
	this->ApplySphericalFeedbackToGrid(ThisParticle, EjectaDensity, EjectaThermalEnergy,
					   EjectaMetalDensity);

      } else{ /* SMS converts directly into DCBH (ORIGINAL) */
	SS->ParticleClass = BH;
	SS->StellarAge = SS->RadiationLifetime; //Record last stellar age
	SS->RadiationLifetime = 1e20;
	printf("%s: DCBH Created from SMS: ParticleClass now %d\t Stellar Age = %f Myr\t " \
	     "Lifetime = %f Myr\n", __FUNCTION__,
	     SS->ParticleClass, SS->StellarAge*TimeUnits/Myr_s,
	     SS->RadiationLifetime*TimeUnits/Myr_s);
      }
    }
    return SUCCESS;
  } 
  /***********************************************************************
                                POPIII
  ************************************************************************/
  else if(SS->ParticleClass == POPIII) {

    fprintf(stderr, "%s: age = %e Myr.\n", __FUNCTION__, Age);

    // SG. When mass is being removed from the grid during star formation.
    if (SS->Mass == 0){
      fprintf(stderr, "%s: Mass is zero, no radiation lifetime check.\n", __FUNCTION__);
      return SUCCESS;
    }
  
    if(Age > SS->RadiationLifetime) {/* Star needs to go supernovae and change type */
      fprintf(stderr, "%s: Age exceeds stellar lifetime. Transition to BH.\n", __FUNCTION__);
      
      /* We now need to convert this particle into a Black Hole if appropriate 
       * 20 Msolar - 40.1 Msolar -> Type II supernova with BH remnant 
       * 40.1 Msolar - 140 Msolar -> DCBH with Heger-Woosley relation
       * 140 Msolar - 260 Msolar -> PISN -> No remnant (delete particle)
       * 260+ Msolar - DCBH of same mass as parent star
       */
      
      printf("%s: End of star's life\n", __FUNCTION__);
      printf("%s: Age = %1.2f Myr\t RadiationLifetime = %1.2f Myr\n", __FUNCTION__,
	     Age*TimeUnits/Myr_s, SS->RadiationLifetime*TimeUnits/Myr_s);
      double StellarMass = SS->Mass*MassConversion/SolarMass; /* In Msolar */
      printf("%s: StellarMass = %lf\n", __FUNCTION__, StellarMass);
      double SNEnergy, HeliumCoreMass, Delta_SF, MetalMass;
      FLOAT Radius = PopIIISupernovaRadius * pc_cm / LengthUnits;
      FLOAT StarLevelCellWidth = this->CellWidth[0][0];
      Radius = max(Radius, 3.5*StarLevelCellWidth);
      FLOAT EjectaVolume = 4.0/3.0 * pi * pow(Radius*LengthUnits, 3);
      FLOAT EjectaDensity = StellarMass * SolarMass / EjectaVolume / DensityUnits;
      FLOAT EjectaMetalDensity = 0.0, EjectaThermalEnergy = 0.0;

      /* pair-instability SNe */
      if (StellarMass >= PISNLowerMass && StellarMass <= PISNUpperMass) {
	HeliumCoreMass = (13./24.) * (StellarMass - 20);
	SNEnergy = (5.0 + 1.304 * (HeliumCoreMass - 64)) * 1e51;
	EjectaMetalDensity = HeliumCoreMass * SolarMass / EjectaVolume / DensityUnits;
	SS->WillDelete = true;
	printf("%s: PISN detected. Particle set for deletion.\n", __FUNCTION__);
	EjectaThermalEnergy = SNEnergy / (StellarMass * SolarMass) / VelocityUnits /
	  VelocityUnits;
	
	this->ApplySphericalFeedbackToGrid(ThisParticle, EjectaDensity, EjectaThermalEnergy,
					   EjectaMetalDensity);
	printf("%s: PISN Feedback completed. Delete particle\n", __FUNCTION__);
      } 
      
      /* Normal Type II SNe: 11 <= M <= 40.1*/
      else if (StellarMass >= TypeIILowerMass && StellarMass <= TypeIIUpperMass) { 
	// For < 20msun  
	if (StellarMass < 20.0) { 
	  SNEnergy = 1e51;
	  MetalMass = 0.1077 + 0.3383 * (StellarMass - 11.0);  // Fit to Nomoto+06
	  // For > 20msun: Hypernova or Core-Collapse depending on flag (should we add the "failed" SNe?)
	  } else { 
	  // SN explosion
	  int bin = search_lower_bound((float*)SNExplosionMass, StellarMass, 0, 5, 5);
	  float frac = (SNExplosionMass[bin+1] - StellarMass) / 
	    (SNExplosionMass[bin+1] - SNExplosionMass[bin]);
	  SNEnergy = 1e51 * (SNExplosionEnergy[bin] + 
			     frac * (SNExplosionEnergy[bin+1] - SNExplosionEnergy[bin]));
	  MetalMass = (SNExplosionMetals[bin] + 
		       frac * (SNExplosionMetals[bin+1] - SNExplosionMetals[bin]));
	  // Heger-Woosley (2002) relation for BHMass
	  HeliumCoreMass = (13./24.) * (StellarMass - 20);
	  StellarMass = HeliumCoreMass; //msun
	  SS->Mass = StellarMass*SolarMass/MassConversion; //code density
	  }
	// Set to BH particle
	SS->ParticleClass = BH;
	SS->StellarAge = SS->RadiationLifetime; //Record last stellar age
	SS->RadiationLifetime = 1e20;
	
	/* SG. Set initial accretion radius of BH to something larger than cell width. 
	   High temperature after SNe results in tiny bondi radius, but no accretion can occur.
	*/
	float mparticle = SS->Mass*dx*dx*dx;
	float *vparticle = SS->ReturnVelocity();
	grid* APGrid = SS->ReturnCurrentGrid();
	int size = APGrid->GetGridSize();
	float *Temperature = new float[size]();
	APGrid->ComputeTemperatureField(Temperature);
	FLOAT BondiHoyleRadius = APGrid->CalculateBondiHoyleRadius(mparticle, vparticle, Temperature);
	SS->AccretionRadius = 1000000*BondiHoyleRadius;
	fprintf(stderr, "%s: Initial accretion radius of BH (x1e6 Bondi rad)= %e pc.\n", __FUNCTION__, 
		SS->AccretionRadius*LengthUnits/pc_cm);
	/* SG. End set accretion radius to BHL radius */
	
	printf("%s: Post-SNe: ParticleClass now %d\t Lifetime = %f Myr\n", __FUNCTION__,
	       SS->ParticleClass, SS->RadiationLifetime*TimeUnits/Myr_s);
	EjectaMetalDensity = MetalMass * SolarMass / EjectaVolume / DensityUnits;
	EjectaThermalEnergy = SNEnergy / (StellarMass * SolarMass) / VelocityUnits /
	  VelocityUnits;

	this->ApplySphericalFeedbackToGrid(ThisParticle, EjectaDensity, EjectaThermalEnergy,
					   EjectaMetalDensity);
      } // SG. End 11 <= M <= 40.1.

      /* DCBH: 40.1 msun <  Mstar < 140 msun - BH has mass set by HW relation */
	else if (TypeIIUpperMass < StellarMass && StellarMass <= PISNLowerMass){
	  // Heger-Woosley (2002) relation for BHMass
	  HeliumCoreMass = (13./24.) * (StellarMass - 20);
	  StellarMass = HeliumCoreMass; //msun
	  SS->Mass = StellarMass*SolarMass/MassConversion; //code density
	  SS->ParticleClass = BH; // Particle class change
	  SS->StellarAge = SS->RadiationLifetime; //Record last stellar age
	  SS->RadiationLifetime = 1e20;
	  printf("%s: DCBH Created from Heger-Woosley relation: ParticleClass now %d\t Stellar Age = %f Myr\t "\
		 "Lifetime = %f Myr\n", __FUNCTION__,
		 SS->ParticleClass, SS->StellarAge*TimeUnits/Myr_s,
		 SS->RadiationLifetime*TimeUnits/Myr_s);
	}
	
	/* SG. DCBH > 260msun - BH has same mass as original star */
	else {
	  SS->ParticleClass = BH;
	  SS->StellarAge = SS->RadiationLifetime; //Record last stellar age
	  SS->RadiationLifetime = 1e20;
	  printf("%s: DCBH Created from POPIII: ParticleClass now %d\t Stellar Age = %f Myr\t " \
		 "Lifetime = %f Myr\n", __FUNCTION__,
		 SS->ParticleClass, SS->StellarAge*TimeUnits/Myr_s,
		 SS->RadiationLifetime*TimeUnits/Myr_s);
	}

      }
      return SUCCESS;
    }

  /***********************************************************************
                                POPII
  ************************************************************************/
  else if(SS->ParticleClass == POPII)
    {

      /* 
       * POPII star clusters are setup to give constant SNe feedback
       * to the grid. 
       * Clusters begin giving out feedback after StarClusterSNeStart
       * They stop producing feedback after StarClusterSNeEnd
       * After each supernova event mass is deducted from the cluster
       * Clusters therefore only decrease in mass after birth. 
       */
      const float StarClusterSNeStart = 4.0;   // Myr after cluster is born
      const float StarClusterSNeEnd = 20.0; // Myr (lifetime of a 8 SolarMass star)
      
      float AgeInMyr = Age * TimeUnits / Myr_s;
      if((AgeInMyr > StarClusterSNeStart) && (AgeInMyr < StarClusterSNeEnd)) {

	printf("%s: POPII Continuous Supernova\n", __FUNCTION__);
	FLOAT StarLevelCellWidth = this->CellWidth[0][0];
	FLOAT Radius = StarClusterSNRadius * pc_cm / LengthUnits;
	if (Radius < 2*StarLevelCellWidth) {
	  Radius = 2*StarLevelCellWidth;
	}
	FLOAT BubbleVolume = (4.0 * pi / 3.0) * Radius * Radius * Radius; /* code volume */
	float dtForThisStar = this->ReturnTimeStep();
	double StellarMass = SS->Mass*MassConversion/SolarMass; /* In Msolar */
	double Delta_SF = StarMassEjectionFraction * StellarMass * dtForThisStar * 
	  TimeUnits / (16.0*Myr_s);  /* Msolar */
	printf("%s: dtForThisStar = %e Myr\n", __FUNCTION__, dtForThisStar * TimeUnits/Myr_s);
	printf("%s: OK I'm going to eject %e Msolar as Energy\n", __FUNCTION__, Delta_SF);
	
	FLOAT EjectaVolume = 4.0/3.0 * pi * pow(Radius*LengthUnits, 3);   /* cm^3 */
	FLOAT EjectaDensity = Delta_SF * SolarMass / EjectaVolume / DensityUnits;   /* code density */
	FLOAT EjectaMetalDensity = EjectaDensity * StarMetalYield; /* code density */
	FLOAT EjectaThermalEnergy = StarClusterSNEnergy / SolarMass /   
	  (VelocityUnits * VelocityUnits); 
	this->ApplySphericalFeedbackToGrid(ThisParticle, EjectaDensity, EjectaThermalEnergy,
					   EjectaMetalDensity);
	/* Remove mass from star following supernova feedback */
	double old_mass = SS->Mass;
	SS->Mass -= Delta_SF * SolarMass / MassConversion; /*Convert to code density */  
	float frac = old_mass / SS->Mass;
	float *Vel = SS->vel;
	float NewVelocity[3] =
	  {
	   (Vel[0]*frac),
	   (Vel[1]*frac),
	   (Vel[2]*frac)
	  };
	SS->SetVelocity(NewVelocity);
	printf("%s: Mass changed from %e Msolar to %e Msolar\n", __FUNCTION__,
	       old_mass*MassConversion/SolarMass, SS->Mass*MassConversion/SolarMass);
      }
      return SUCCESS;
    }
  else if(SS->ParticleClass == BH)
    {
      /***********************************************************************
                                MBH_THERMAL
      ************************************************************************/
      if(SmartStarBHFeedback == FALSE)
	return SUCCESS;
      // Similar to Supernova, but here we assume the followings:
      // EjectaDensity = 0.0
      // EjectaMetalDensity = 0.0
      float  EjectaDensity = 0.0, EjectaMetalDensity = 0.0;
      // The unit of EjectaThermalEnergy = ergs/cm^3, not ergs/g
      if (SmartStarBHThermalFeedback == TRUE) {
	printf("%s: eta_disk = %f\n", __FUNCTION__, SS->eta_disk);
	float epsilon = SS->eta_disk/(1 - SS->eta_disk);
	
	/* find mdot */
	float mdot = SS->AccretionRate[SS->TimeIndex];  //CodeMass/CodeTime
	float accrate = mdot*MassUnits/(SolarMass*TimeUnits)*3.154e7; //in Msolar/yr
	float mdot_cgs = mdot*MassUnits/TimeUnits; //g/s
	//printf("%s: dx = %e\t MassConversion = %e\n", __FUNCTION__, dx, MassConversion);
	printf("%s: AccretionRate = %e Msolar/yr %e (code) TimeIndex = %d\n", __FUNCTION__,
	       accrate, SS->AccretionRate[SS->TimeIndex], SS->TimeIndex);
	
	
	float EjectaVolumeCGS = 4.0/3.0 * PI * pow(SS->AccretionRadius*LengthUnits, 3);
	float EjectaVolume = 4.0/3.0 * PI * pow(SS->AccretionRadius, 3);
	
	float BHMass =  SS->ReturnMass()*MassConversion/SolarMass; //In solar masses
	float eddrate = 4*M_PI*GravConst*BHMass*mh/(SS->eta_disk*clight*sigma_thompson); // Msolar/s
	eddrate = eddrate*3.154e7; //in Msolar/yr
	printf("%s: Eddrate = %e Msolar/yr AccRate = %e Msolar/yr\n", __FUNCTION__, 
	       eddrate, accrate);
	if(SmartStarSuperEddingtonAdjustment == TRUE) {
	  if(accrate > eddrate) {
	    printf("%s: We are accreting at super-Eddington rates. Modifying radiative efficiency\n", __FUNCTION__);
	    float mue = 1.22, a = 0.7;
	    float Ledd = 4*M_PI*GravConst*BHMass*SolarMass*mh*mue*clight/sigma_thompson; //cgs
	    float medddot = 16.0*Ledd/(clight*clight); //cgs
	    /* Apply Madau fit to calculate Luminosity */
	    float LSuperEdd = Ledd*MadauFit(a, accrate*SolarMass/3.154e7, medddot); //cgs
	    epsilon = LSuperEdd/(mdot_cgs*clight*clight);
	    printf("%s: Using the Madau fit raditive efficiency calculated as %e\n", __FUNCTION__, epsilon);
	  }
	}
	/* When injected energy is uniform throughout the volume;
	 * The unit of EjectaThermalEnergy is CodeMass*CodeVelocity^2
	 * EjectaThermalEnergy is added to each cell normalised by the 
	 * totalEjectaVolume. Hence the units of EjectaThermalEnergy are EnergyUnits/VolumeUnits
	 * We calculate the SmartStarDiskEnergyCoupling as (v_wind/(2*c)). To do this we
	 * must fix v_wind. For v_wind we choose 0.1 c (C.-A. Faucher-Giguere, E. Quataert Arxiv:1204.2547)
	 */
	float SmartStarDiskEnergyCoupling = 0.05;
	float EjectaThermalEnergy = SmartStarDiskEnergyCoupling * epsilon * dt * 
	  mdot*clight*clight/(VelocityUnits*VelocityUnits*EjectaVolume); 
	
	/* Ramp up over RAMPTIME yrs */
	float Age = Time - SS->BirthTime;
	float BH_Age = (Age - SS->StellarAge)*TimeUnits/yr_s;
	if(BH_Age < RAMPTIME)
	  {
	    printf("BH Age = %e yrs, ramp = %e\n", BH_Age, BH_Age/(float)RAMPTIME);
	    EjectaThermalEnergy *= BH_Age/(float)RAMPTIME;
	  }
	EjectaDensity = 0.0;
	EjectaMetalDensity = 0.0;
	this->ApplySphericalFeedbackToGrid(ThisParticle, EjectaDensity, EjectaThermalEnergy,
					   EjectaMetalDensity);
	
      }
  
      /***********************************************************************
                                 MBH_JETS
      ************************************************************************/
      if(SmartStarBHJetFeedback == TRUE) {
    
	// Inject bipolar jets along the direction of the angular momentum 
	// vector L of the MBH particle (angular momentum accreted thus far)
	// or along the z-axis  - Ji-hoon Kim, Nov.2009
	int i = 0, j = 0, k = 0;
#define MAX_SUPERCELL_NUMBER 1000
	int SUPERCELL = 1; //2 for supercell of 5 cells wide = 5^3  
	int ind_cell_inside[MAX_SUPERCELL_NUMBER], ind_cell_edge[MAX_SUPERCELL_NUMBER];
	float nx_cell_edge[MAX_SUPERCELL_NUMBER], ny_cell_edge[MAX_SUPERCELL_NUMBER], 
	  nz_cell_edge[MAX_SUPERCELL_NUMBER], anglefactor[MAX_SUPERCELL_NUMBER] = {0};
	int n_cell_inside = 0, n_cell_edge = 0, ibuff = NumberOfGhostZones;
	int ii = 0, jj = 0, kk = 0, r_s = 0, ic = 0, sign = 0;
	float m_cell_inside = 0.0, m_cell_edge = 0.0;
	float L_x, L_y, L_z, L_s, nx_L = 0.0, ny_L = 0.0, nz_L = 0.0, costheta = cos(OPENING_ANGLE);
	float SSMass = SS->ReturnMass();
	float totalenergybefore = 0.0, totalenergyafter = 0.0, totalenergyadded = 0.0;
	float sumkeadded = 0.0;
	
	if (SmartStarBHJetFeedback == FALSE || SS->MassToBeEjected*MassUnits/SolarMass < 1e-10) {
	  return SUCCESS;
	}
	
	/* i, j, k are the number of cells from the edge of the grid to the smartstar*/
	i = (int)((pos[0] - this->CellLeftEdge[0][0]) / dx);
	j = (int)((pos[1] - this->CellLeftEdge[1][0]) / dx);
	k = (int)((pos[2] - this->CellLeftEdge[2][0]) / dx);
	
	/* Note that we need to inject feedback only for the finest grid the SS belongs to */
	
	if (i < ibuff || i > this->GridDimension[0]-ibuff-1 ||
	    j < ibuff || j > this->GridDimension[1]-ibuff-1 || 
	    k < ibuff || k > this->GridDimension[2]-ibuff-1 ||
	    this == NULL ||
	    SS->level < MaximumRefinementLevel) {
	  fprintf(stdout, "grid::AddFS: MBH_JETS - MBH doesn't belong to this grid.\n"); 
	  return SUCCESS;
	}
	
	
	/* find mdot */
	float mdot = SS->AccretionRate[SS->TimeIndex];
	float BHMass =  SS->ReturnMass()*MassConversion/SolarMass; //In solar masses
	float eddrate = 4*M_PI*GravConst*BHMass*SolarMass*mh/(SS->eta_disk*clight*sigma_thompson); // g/s
	eddrate = eddrate*3.154e7/SolarMass; //in Msolar/yr
	
	float AccretionRate = mdot*MassUnits/(SolarMass*TimeUnits); //in Msolar/s
	/* 
	 * Now in the case where we are subgriding the accretion formalism
	 * re-calculate the actual accretion rate and check if we are in the correct band
	 */
	//AccretionRate *= SS->epsilon_deltat;
	
	/* Debug */
	printf("%s: Eddrate = %e Msolar/yr AccRate = %e Msolar/yr\t Ratio = %f\n", __FUNCTION__,
	       eddrate, AccretionRate*3.154e7, AccretionRate*3.154e7/eddrate);
	printf("%s: dx = %e\t MassConversion = %e\n", __FUNCTION__, dx, MassConversion);
	printf("%s: AccretionRate (*deltat) = %e Msolar/yr %e (code) TimeIndex = %d\n", __FUNCTION__,
	       AccretionRate*3.154e7, SS->AccretionRate[SS->TimeIndex], SS->TimeIndex);
	float MassEjected = SS->NotEjectedMass + SS->MassToBeEjected; //code mass    
	
	
	SS->NotEjectedMass = MassEjected;
#if IMPOSETHRESHOLD
#if SSFEED_DEBUG
	printf("SSFEED_DEBUG: %s: Mass Accumulated thus far = %e Msolar (Threshold = %e Msolar)\n",
	       __FUNCTION__, SS->NotEjectedMass*MassUnits/SolarMass, SS->EjectedMassThreshold);
#endif
	if (SS->NotEjectedMass*MassUnits/SolarMass <= SS->EjectedMassThreshold) {
	  fprintf(stdout, "grid::AddFS: MBH_JETS - accumulated mass (%f Msolar) not passed threshold (%f Msolar).\n",
		  SS->NotEjectedMass*MassUnits/SolarMass, SS->EjectedMassThreshold);
	  return SUCCESS;
	}
#endif
	
	if(AccretionRate*3.154e7/eddrate < 1e-30 || AccretionRate*3.154e7/eddrate > 1.0)
	  {
	    printf("%s: AccrateionRateRatio = %f. We are in the right band to release jets\n", __FUNCTION__,
		   AccretionRate*3.154e7/eddrate);
	  }
	else
	  {
	    printf("%s: AccrateionRateRatio = %f. No jets this time\n", __FUNCTION__,
		   AccretionRate*3.154e7/eddrate);
	    return SUCCESS;
	  }
	/*end Debug*/
	
	
#if SSFEED_DEBUG
	
	printf("SSFEED_DEBUG: %s: Mass Accreted = %e Msolar\t Mass to be Ejected  = %e Msolar\n",
	       __FUNCTION__, mdot*dt*MassUnits/SolarMass, MassEjected*MassUnits/SolarMass);
#endif
	
	if (i < ibuff+SUPERCELL || i > this->GridDimension[0]-ibuff-SUPERCELL-1 || 
	    j < ibuff+SUPERCELL || j > this->GridDimension[1]-ibuff-SUPERCELL-1 ||
	    k < ibuff+SUPERCELL || k > this->GridDimension[2]-ibuff-SUPERCELL-1) {
	  fprintf(stdout, "grid::AddFS: MBH_JETS - supercell not contained; accumulated mass (%g MSolar).\n",
		  SS->NotEjectedMass*MassUnits/SolarMass); 
	  
	  // if the supercell issue hasn't allowed the jet injection for too long,
	  // issue a warning signal and output the current hierarchy at CheckForOutput
	  if (SS->NotEjectedMass*MassUnits/SolarMass > 2.0 * SS->EjectedMassThreshold) { 
	    fprintf(stdout, "grid::AddFS: MBH_JETS - jets haven't been ejected for too long!\n");
	  } 
	  
	  // otherwise, just proceed and do it later
	  return SUCCESS;
	}
	printf("SSFEED_DEBUG: %s: Lets Eject!!!!!!!!!!!\n", __FUNCTION__);
#if IMPOSETHRESHOLD
	float MassKeptInReserve = max(MassEjected - THRESHOLDFRACTION*SolarMass/MassUnits, 0.0);
	MassEjected = MassEjected - MassKeptInReserve;
#endif
	printf("Cumulative Mass to be ejected in jet will be %f Msolar\n", MassEjected*MassUnits/SolarMass);
	/* Find the directional vector n_L of angular momentum accreted thus far */
	if(SS->CalculateAccretedAngularMomentum() == FAIL) {
	  return FAIL;
	}
	L_x = SS->Accreted_angmom[0];
	L_y = SS->Accreted_angmom[1];
	L_z = SS->Accreted_angmom[2]; 
	L_s = sqrt(pow(L_x,2) + pow(L_y,2) + pow(L_z,2));
	nx_L = L_x/L_s;  //normalized directional vector
	ny_L = L_y/L_s;
	nz_L = L_z/L_s;
	L_s = sqrt(pow(nx_L,2) + pow(ny_L,2) + pow(nz_L,2));
	printf("%s: Angular momentum = %e %e %e\t L_s = %e\n", __FUNCTION__, 
	       nx_L, ny_L, nz_L, L_s);
	//nx_L = 0.0;
	//ny_L = 0.0;
	//nz_L = 1.0;
	//L_s = sqrt(pow(nx_L,2) + pow(ny_L,2) + pow(nz_L,2));
	//printf("%s: Angular momentum = %e %e %e\t L_s = %e\n", __FUNCTION__, 
	//	   nx_L, ny_L, nz_L, L_s);
	
	/* Loop over the supercell around the MBH particle (5 * 5 * 5 = 125 cells, 
	   but only the edges), and record the cells eligible for jet injection */
	int nsupercells = 0;
	for (kk = -SUPERCELL; kk <= SUPERCELL; kk++) {
	  for (jj = -SUPERCELL; jj <= SUPERCELL; jj++) {
	    for (ii = -SUPERCELL; ii <= SUPERCELL; ii++) {
	      nsupercells++;
	      r_s = sqrt(pow(ii,2) + pow(jj,2) + pow(kk,2));	    
	      if (fabs(ii) != SUPERCELL && fabs(jj) != SUPERCELL && fabs(kk) != SUPERCELL) {  //if not on edges
		//	    printf("%s: Inside: CosTheta = %f\t cellangle = %f (%f degrees)\n", __FUNCTION__, 
		//	   costheta, fabs((ii*nx_L + jj*ny_L + kk*nz_L)/r_s), 
		//	   (360/pi)*acos(fabs((ii*nx_L + jj*ny_L + kk*nz_L)/r_s)));
		ind_cell_inside[n_cell_inside] = i+ii+(j+jj+(k+kk)*this->GridDimension[1])*this->GridDimension[0];
		m_cell_inside += this->BaryonField[DensNum][ind_cell_inside[n_cell_inside]] * 
		  pow(dx, 3);
		n_cell_inside++;
		
	      } else {  //if on edges	    
		if (fabs((ii*nx_L + jj*ny_L + kk*nz_L)/r_s) > costheta) { 
		  //printf("%s: Edge: CosTheta = %f\t cellangle = %f (%f degrees)\n", __FUNCTION__, 
		  //     costheta, fabs((ii*nx_L + jj*ny_L + kk*nz_L)/r_s), 
		  //     (360/pi)*acos(fabs((ii*nx_L + jj*ny_L + kk*nz_L)/r_s)));
		  anglefactor[n_cell_edge] = fabs((ii*nx_L + jj*ny_L + kk*nz_L)/r_s);
		  ind_cell_edge[n_cell_edge] = i+ii+(j+jj+(k+kk)*this->GetGridDimension(1))*this->GetGridDimension(0);
		  nx_cell_edge[n_cell_edge]  = ii / r_s;  //directional vector
		  ny_cell_edge[n_cell_edge]  = jj / r_s;
		  nz_cell_edge[n_cell_edge]  = kk / r_s;
		  m_cell_edge += this->BaryonField[DensNum][ind_cell_edge[n_cell_edge]] * 
		    pow(this->GetCellWidth(0, 0), 3);
		  
		  totalenergybefore +=  this->BaryonField[TENum][ind_cell_edge[n_cell_edge]];
		  n_cell_edge++;
		  
		} 
		
	      }  
	      
	    }  // END ii-direction
	  }  // END jj-direction
	}  // END kk-direction
#if SSFEED_DEBUG
	printf("%s: n_cell_inside = %d\n", __FUNCTION__,  n_cell_inside);
	printf("%s: n_cell_edge = %d\n", __FUNCTION__, n_cell_edge);
	printf("%s: nsupercells = %d\n",  __FUNCTION__, nsupercells);
#endif
	/* Calculate the jet density 
	 * This is the mass of the ejected mass + the mass at the cell edges
	 */
	
	float  rho_jet = (MassEjected) / 
	  ((float)n_cell_edge * pow(this->GetCellWidth(0,0), 3)); //In code units
#if SSFEED_DEBUG
	printf("%s: rho_jet (per cell) = %e cc", __FUNCTION__, rho_jet*DensityUnits/mh);
#endif
	
	/* Calculate MBHJetsVelocity using energy conservation below:
	   SmartStarFeedbackEnergyCoupling = The fraction of feedback energy that is
	   mechanically (for MBH_JETS) coupled to the gas.
	   SmartStarFeedbackRadiativeEfficiency - The radiative efficiency of a black hole.
	   
	   MBHFeedbackEnergyCoupling * MBHFeedbackRadiativeEfficiency * Mdot * c^2 
	   = 0.5 * MBHFeedbackMassEjectionFraction * Mdot * (MBHJetsVelocity)^2                
	   
	   Note that EjectaThermalEnergy is never used; MBHFeedbackEnergyCoupling 
	   should now be calculated considering gravitational redshift (Kim et al. 2010) 
	   
	   float MBHJetsVelocity = clight * sqrt( 2 * MBHFeedbackEnergyCoupling * MBHFeedbackRadiativeEfficiency 
	   / MBHFeedbackMassEjectionFraction ) / VelocityUnits;
	*/
	
	/* This is really a bit of a cod and it may be better to set the jet velocity as some fraction of the 
	 * speed of light. */
	float MBHJetsVelocity = clight * SmartStarJetVelocity/VelocityUnits; //code velocity
	
	/* Ramp up over RAMPTIME yrs */
	float Age = this->ReturnTime() - SS->BirthTime;
	Age = Age*TimeUnits/3.154e7;
	if(Age < RAMPTIME)
	  {
	    printf("%s: Too early for jets. Age = %f yrs\n", __FUNCTION__, Age);
	    return SUCCESS;
	    MBHJetsVelocity = MBHJetsVelocity*Age/(float)RAMPTIME;
	  }
	
	float jetenergy = 0.5*MBHJetsVelocity*MBHJetsVelocity;
#if SSFEED_DEBUG
	printf("%s: Age = %f yrs\n", __FUNCTION__, Age);
	printf("%s: Jet Energy = %e (specific = %e)\n", __FUNCTION__, jetenergy*MassEjected, jetenergy);
	printf("SSFEED_DEBUG: %s: MBHJetsVelocity = %e of clight (%f km/s)\n", __FUNCTION__, 
	       MBHJetsVelocity * VelocityUnits/clight,
	       MBHJetsVelocity*VelocityUnits/1e5 );
	printf("SSFEED_DEBUG: %s: SmartStarJetVelocity = %e\n", __FUNCTION__, SmartStarJetVelocity);
#endif
	if (MBHJetsVelocity * VelocityUnits > 0.99*clight) {
	  ENZO_VFAIL("grid::AddFS: MBHJetsVelocity is ultra-relativistic! (%g/ %g/ %g/ %g c)\n",
		     MBHFeedbackEnergyCoupling, MBHFeedbackRadiativeEfficiency, 
		     MBHFeedbackMassEjectionFraction, MBHJetsVelocity * VelocityUnits / clight);
	}
	
	/* Finally, add the jet feedback at the edges (outer part of the supercell) */
	
	for (ic = 0; ic < n_cell_edge; ic++) {
	  
	  int index = ind_cell_edge[ic];
	  float angle = anglefactor[ic];
	  float cellnumberdensity = this->BaryonField[DensNum][index]*DensityUnits/mh;
	  
	  /* Update velocities and TE; note that we now have kinetic (jet) energy added, so 
	     for DualEnergyFormalism = 0 you don't have to update any energy field */
	  
	  sign = sign(nx_cell_edge[ic]*nx_L + ny_cell_edge[ic]*ny_L + nz_cell_edge[ic]*nz_L);
	  
	  /* Calculate grid velocity: the actual veloctiy injected in supercell edges.
	     This is different from MBHJetsVelocity because it is the mass-weighted average 
	     between MBHJetsVelocity and the original veloctiy in grid */  
	  float oldvel[3] = {this->BaryonField[Vel1Num][index], 
			     this->BaryonField[Vel2Num][index],
			     this->BaryonField[Vel3Num][index]};
	  float oldcellmass = this->BaryonField[DensNum][index] * pow(this->GetCellWidth(0,0), 3);
	  float energybefore = this->BaryonField[TENum][index];
	  
	  
#if DENSITY_WEIGHTED
	  
	  int dim = 0;
	  if (GENum >= 0 && DualEnergyFormalism) 
	    for (dim = 0; dim < GridRank; dim++)
	      this->BaryonField[TENum][index] -= 
		0.5 * this->BaryonField[Vel1Num+dim][index] * 
		this->BaryonField[Vel1Num+dim][index];
	  
	  this->BaryonField[Vel1Num][index] = (this->BaryonField[DensNum][index] *  this->BaryonField[Vel1Num][index] +
					       MassEjected / ((float)n_cell_edge*pow(this->GetCellWidth(0,0), 3)) *
					       sign * nx_L * MBHJetsVelocity) / 
	    (this->BaryonField[DensNum][index] + 
	     MassEjected / ((float)n_cell_edge*pow(this->GetCellWidth(0,0), 3)));  
	  this->BaryonField[Vel2Num][index] = (this->BaryonField[DensNum][index] * 
					       this->BaryonField[Vel2Num][index] +
					       MassEjected / ((float)n_cell_edge*pow(this->GetCellWidth(0,0), 3)) *
					       sign * ny_L * MBHJetsVelocity) / 
	    (this->BaryonField[DensNum][index] + 
	     MassEjected / ((float)n_cell_edge*pow(this->GetCellWidth(0,0), 3)));
	  this->BaryonField[Vel3Num][index] = (this->BaryonField[DensNum][index] * 
					       this->BaryonField[Vel3Num][index] +
					       MassEjected / ((float)n_cell_edge*pow(this->GetCellWidth(0,0), 3)) *
					       sign * nz_L * MBHJetsVelocity) / 
	    (this->BaryonField[DensNum][index] + 
	     MassEjected / ((float)n_cell_edge*pow(this->GetCellWidth(0,0), 3)));
	  
	  float newvel[3] = {this->BaryonField[Vel1Num][index], 
			     this->BaryonField[Vel2Num][index],
			     this->BaryonField[Vel3Num][index]};
	  float newvelmag = sqrt(newvel[0]*newvel[0] + newvel[1]*newvel[1] + newvel[2]*newvel[2]);
	  float energytoadd = 0.5*newvelmag*newvelmag;
	  
	  if (GENum >= 0 && DualEnergyFormalism) 
	    for (dim = 0; dim < GridRank; dim++)
	      this->BaryonField[TENum][index] += 
		0.5 * this->BaryonField[Vel1Num+dim][index] * 
		this->BaryonField[Vel1Num+dim][index];
	  
#elif SINE_WAVE
	  this->BaryonField[Vel1Num][index] += sign * nx_L * MBHJetsVelocity * angle;
	  this->BaryonField[Vel2Num][index] += sign * ny_L * MBHJetsVelocity * angle;
	  this->BaryonField[Vel3Num][index] += sign * nz_L * MBHJetsVelocity * angle;
	  float v_ejecta[3] = { sign * nx_L * MBHJetsVelocity * angle, 
				sign * ny_L * MBHJetsVelocity * angle, 
				sign * nz_L * MBHJetsVelocity * angle};
	  float v_ejectamag = sqrt(v_ejecta[0]*v_ejecta[0] + v_ejecta[1]*v_ejecta[1] + v_ejecta[2]*v_ejecta[2]);
	  float keadded = 0.5*(MassEjected/(float)n_cell_edge)*v_ejectamag*v_ejectamag;
	  
	  /* Update the new total energy. 
	   * keadded is due to the new grid velocities. 
	   * eint does not change
	   */
	  if (GENum >= 0 && DualEnergyFormalism) 
	    this->BaryonField[TENum][index] += keadded/MassEjected;
	  
	  
#else
	  
	  this->BaryonField[Vel1Num][index] += sign * nx_L * MBHJetsVelocity;
	  this->BaryonField[Vel2Num][index] += sign * ny_L * MBHJetsVelocity;
	  this->BaryonField[Vel3Num][index] += sign * nz_L * MBHJetsVelocity;
	  float newvel[3] = {this->BaryonField[Vel1Num][index], 
			     this->BaryonField[Vel2Num][index],
			     this->BaryonField[Vel3Num][index]};
	  float oldvelmag = sqrt(oldvel[0]*oldvel[0] + oldvel[1]*oldvel[1] + oldvel[2]*oldvel[2]);
	  float newvelmag = sqrt(newvel[0]*newvel[0] + newvel[1]*newvel[1] + newvel[2]*newvel[2]);
	  float kebefore = 0.5*(oldcellmass)*oldvelmag*oldvelmag;
	  float keafter = 0.5*(oldcellmass + (MassEjected/(float)n_cell_edge))*newvelmag*newvelmag;
	  float v_ejecta[3] = { sign * nx_L * MBHJetsVelocity, 
				sign * ny_L * MBHJetsVelocity, 
				sign * nz_L * MBHJetsVelocity};
	  float v_ejectamag = sqrt(v_ejecta[0]*v_ejecta[0] + v_ejecta[1]*v_ejecta[1] + v_ejecta[2]*v_ejecta[2]);
	  float keadded = 0.5*(MassEjected/(float)n_cell_edge)*v_ejectamag*v_ejectamag;
	  
	  /* Update the new total energy. 
	   * keadded is due to the new grid velocities. 
	   * eint does not change
	   */
	  totalenergyadded += keadded;
	  //#if SSFEED_DEBUG
	  //printf("%s: SSFEED_DEBUG: Adding %e energy to TE. This increases TE by a factor of %e\n",
	  //	     __FUNCTION__,  keadded/MassEjected, (energybefore + keadded/MassEjected)/energybefore);
	  //#endif
	  if (GENum >= 0 && DualEnergyFormalism) 
	    this->BaryonField[TENum][index] += keadded/MassEjected;
	  
	  totalenergyafter +=  this->BaryonField[TENum][index];
	  //#if SSFEED_DEBUG
	  // printf("EnergyConservation: Total Energy Added = %e\n", totalenergyadded);
	  //printf("EnergyConservation: Delta Grid Energy  = %e (with mass = %e)\t Jet Energy = %e\n", 
	  //     totalenergyafter - totalenergybefore, (totalenergyafter - totalenergybefore)*MassEjected,
	  //     jetenergy);
	  //#endif
	  totalenergyafter = 0; totalenergybefore = 0;totalenergyadded= 0;
#endif
	  
	  //return SUCCESS; //works
	  //printf("DeltaGrid = %e\n",  this->BaryonField[TENum][index] - energybefore);
	  
	  /* Update density, species and colour fields */
	  float OldDensity = this->BaryonField[DensNum][index];
	  float increase = (OldDensity + rho_jet) / OldDensity;
	  this->BaryonField[DensNum][index] += rho_jet;
	  //printf("%s: Increase in Density due to Jet = %e\n", __FUNCTION__, increase);
	  if (MultiSpecies) {
	    this->BaryonField[DeNum][index] *= increase;
	    this->BaryonField[HINum][index] *= increase;
	    this->BaryonField[HIINum][index] *= increase;
	    this->BaryonField[HeINum][index] *= increase;
	    this->BaryonField[HeIINum][index] *= increase;
	    this->BaryonField[HeIIINum][index] *= increase;
	  }
	  if (MultiSpecies > 1) {
	    this->BaryonField[HMNum][index] *= increase;
	    this->BaryonField[H2INum][index] *= increase;
	    this->BaryonField[H2IINum][index] *= increase;
	  }
	  if (MultiSpecies > 2) {
	    this->BaryonField[DINum][index] *= increase;
	    this->BaryonField[DIINum][index] *= increase;
	    this->BaryonField[HIINum][index] *= increase;
	    this->BaryonField[HDINum][index] *= increase;
	  }
    
	  
	}
#if IMPOSETHRESHOLD  
	SS->NotEjectedMass = MassKeptInReserve;
	printf("Mass left over for next jet = %f Msolar\n", SS->NotEjectedMass*MassUnits/SolarMass);
	SS->EjectedMassThreshold = THRESHOLDFRACTION;
	if(SS->EjectedMassThreshold < 1e-3)
	  SS->EjectedMassThreshold = 1e-3;
	SS->NotEjectedMass = 0.0;
#if SSFEED_DEBUG
	printf("%s: New threshold Mass set to %e Msolar, total BH Mass = %e Msolar\n", __FUNCTION__, 
	       SS->EjectedMassThreshold, SS->Mass*MassConversion/SolarMass);
#endif
#else
	SS->NotEjectedMass = 0.0;
#endif
      }
    }
  return SUCCESS;
}
  
