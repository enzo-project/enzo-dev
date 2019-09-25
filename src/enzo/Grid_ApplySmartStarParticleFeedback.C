/***********************************************************************
/
/  Apply feedback to the temporary grid for the Smart Star Particles
/ The feedback methods dealt with here are the Thermal and Kinetic modes
/
/  written by: John Regan
/  date:       December, 2017
/
/  note: Based on methods originally implemented by Stephen Skory
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
#define RAMPTIME 4000.0 //1e4
#define DENSITY_WEIGHTED 1
#define SINE_WAVE       0
#define IMPOSETHRESHOLD 1
#define THRESHOLDFRACTION 1   //Solarmasses ejected per jet event
#define OPENING_ANGLE pi/360.0  //pi/3.9
int grid::ApplySmartStarParticleFeedback(ActiveParticleType** ThisParticle){

  /* Return if this doesn't involve us */
  if (MyProcessorNumber != ProcessorNumber) 
    return SUCCESS;
 
  if(SmartStarBHFeedback == FALSE)
    return SUCCESS;
  if (SmartStarBHJetFeedback == FALSE && SmartStarBHThermalFeedback == FALSE) {
    return SUCCESS;
  }

  ActiveParticleType_SmartStar *SS = static_cast<ActiveParticleType_SmartStar*>(* ThisParticle);
  if(SS->ParticleClass != BH)
    return SUCCESS;
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
    TimeUnits = 1, VelocityUnits = 1, MassUnits = 1,
    PressureUnits = 0, GEUnits = 0, VelUnits = 0;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, this->ReturnTime()) == FAIL) {
        ENZO_FAIL("Error in GetUnits.");
  }
  MassUnits = DensityUnits * POW(LengthUnits,3);
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
  if(SS->AccretionRate[SS->TimeIndex] <= 0.0) {
    float mdot = SS->AccretionRate[SS->TimeIndex];  //CodeMass/CodeTime
    float MassConversion = (float) (dx*dx*dx * double(MassUnits));  //convert to g         
    float accrate = mdot*MassUnits/(SolarMass*TimeUnits)*3.154e7; //in Msolar/s      
    printf("%s: AccretionRate = %e Msolar/yr %e (code)\n", __FUNCTION__,
           accrate, SS->AccretionRate[SS->TimeIndex]);
    printf("%s: Returning until accretion rate is updated\n", __FUNCTION__);
    return SUCCESS;
  }
  /***********************************************************************
                                MBH_THERMAL
  ************************************************************************/
  int CellsModified = 0;
  // Similar to Supernova, but here we assume the followings:
  // EjectaDensity = 0.0
  // EjectaMetalDensity = 0.0
  // The unit of EjectaThermalEnergy = ergs/cm^3, not ergs/g
  if (SmartStarBHThermalFeedback == TRUE) {
    float epsilon = SS->eta_disk;
    FLOAT outerRadius2 = POW(1.2*rad, 2.0);
    /* find mdot */
    float mdot = SS->AccretionRate[SS->TimeIndex];  //CodeMass/CodeTime
    
     /* Debug */
    float MassConversion = (float) (dx*dx*dx * double(MassUnits));  //convert to g                                                                    
    float accrate = mdot*MassUnits/(SolarMass*TimeUnits)*3.154e7; //in Msolar/yr
    float mdot_cgs = mdot*MassUnits/TimeUnits; //g/s
    //printf("%s: dx = %e\t MassConversion = %e\n", __FUNCTION__, dx, MassConversion);
    printf("%s: AccretionRate = %e Msolar/yr %e (code) TimeIndex = %d\n", __FUNCTION__,
           accrate, SS->AccretionRate[SS->TimeIndex], SS->TimeIndex);
    /*end Debug*/
    float newGE = 0.0, oldGE = 0.0;
    float maxGE = MAX_TEMPERATURE / (TemperatureUnits * (Gamma-1.0) * 0.6);
   
    FLOAT radius2 = 0.0;
    float EjectaVolumeCGS = 4.0/3.0 * PI * pow(SS->AccretionRadius*LengthUnits, 3);
    float EjectaVolume = 4.0/3.0 * PI * pow(SS->AccretionRadius, 3);
    printf("%s: OuterRadius = %e\n", __FUNCTION__, sqrt(outerRadius2)*LengthUnits/pc_cm);
    CellsModified = 0;
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
    float Age = this->ReturnTime() - SS->BirthTime;
    Age = Age*TimeUnits/3.154e7 - SmartStarSMSLifetime;
    if(Age < RAMPTIME)
      {
	printf("BH Age = %f yrs, ramp = %f\n", Age, Age/(float)RAMPTIME);
	EjectaThermalEnergy *= Age/(float)RAMPTIME;
      }
    for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
      for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
	int index = GRIDINDEX_NOGHOST(GridStartIndex[0],j,k);
	for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
	
	  radius2 = POW(CellLeftEdge[0][i] + 0.5*dx - pos[0],2.0) +
	    POW(CellLeftEdge[1][j] + 0.5*dx - pos[1],2.0) +
	    POW(CellLeftEdge[2][k] + 0.5*dx - pos[2],2.0);
	  if (radius2 < outerRadius2) {
	    float r1 = sqrt(radius2) / rad;
	    float norm = 0.98;
	    float ramp = norm*(0.5 - 0.5 * tanh(10.0*(r1-1.0)));
	    /* 1/1.2^3 factor to dilute the density since we're
	       depositing a uniform ejecta in a sphere of 1.2*radius
	       without a ramp.  The ramp is only applied to the
	       energy*density factor. */
	    float factor = 0.578704;
	    
	    float Density = this->BaryonField[DensNum][index];
	    /* Get specific energy */
	    if (GENum >= 0 && DualEnergyFormalism) {

	      /* When injected energy is uniform throughout the volume;
		 EjectaThermalEnergy in EnergyUnits/VolumeUnits */
	      oldGE =  this->BaryonField[GENum][index];
	      newGE = (Density * this->BaryonField[GENum][index] +
		       ramp * factor * EjectaThermalEnergy) / Density;

	      newGE = min(newGE, maxGE);  
	      //printf("%s: Energy Before = %e\t Energy injected = %e\t Increase = %e\n", __FUNCTION__, 
	      //	     oldGE,ramp * factor * EjectaThermalEnergy / Density, (newGE - oldGE)/oldGE);
	      fflush(stdout);
	      
	      this->BaryonField[GENum][index] = newGE;
	      this->BaryonField[TENum][index] = newGE;

	      for (int dim = 0; dim < GridRank; dim++)
		this->BaryonField[TENum][index] += 
		  0.5 * this->BaryonField[Vel1Num+dim][index] * 
		  this->BaryonField[Vel1Num+dim][index];

	      //printf("%s: Increase in GE energy is %e\n", __FUNCTION__, (newGE - oldGE)/oldGE);
	      
	    } else {

	      newGE = (Density * this->BaryonField[TENum][index] +
		       ramp * factor * EjectaThermalEnergy) / Density;

	      newGE = min(newGE, maxGE);  
	      this->BaryonField[TENum][index] = newGE;

	    } //end if(GENum >= 0 && DualEnergyFormalism)

	    /* Update species and colour fields */

	    float fh = CoolData.HydrogenFractionByMass;
	    float fhez = (1-fh);
	    float ionizedFraction = 0.999;
	    if (MultiSpecies) {
	      this->BaryonField[DeNum][index] =
		this->BaryonField[DensNum][index] * ionizedFraction;
	      this->BaryonField[HINum][index] = 
		this->BaryonField[DensNum][index] * fh * (1-ionizedFraction);
	      this->BaryonField[HIINum][index] =
		this->BaryonField[DensNum][index] * fh * ionizedFraction;
	      this->BaryonField[HeINum][index] =
		0.5*this->BaryonField[DensNum][index] * fhez * (1-ionizedFraction);
	      this->BaryonField[HeIINum][index] =
		0.5*this->BaryonField[DensNum][index] * fhez * (1-ionizedFraction);
	      this->BaryonField[HeIIINum][index] =
		this->BaryonField[DensNum][index] * fhez * ionizedFraction;
	    }
#ifdef IDONTSEEWHYTHESEELEMENTSAREEFFECTED
	    if (MultiSpecies > 1) {
	      this->BaryonField[HMNum][index] = tiny_number * this->BaryonField[DensNum][index];
	      this->BaryonField[H2INum][index] = 
	      	tiny_number * this->BaryonField[DensNum][index];
	      this->BaryonField[H2IINum][index] = 
	      	tiny_number * this->BaryonField[DensNum][index];
	    }
	    if (MultiSpecies > 2) {
	      this->BaryonField[DINum][index] = this->BaryonField[DensNum][index] * fh *
		CoolData.DeuteriumToHydrogenRatio * (1-ionizedFraction);
	      this->BaryonField[DIINum][index] = this->BaryonField[DensNum][index] * fh *
		CoolData.DeuteriumToHydrogenRatio * ionizedFraction;
	      this->BaryonField[HDINum][index] = 
		tiny_number * this->BaryonField[DensNum][index];
	    }
#endif
	    CellsModified++;
	    
	  } // END if inside radius
	}  // END i-direction
      }  // END j-direction
    }  // END k-direction
    //printf("CellsModified = %d\n", CellsModified);

  }  // END MBH_THERMAL
  
  CellsModified = 0;
   /***********************************************************************
                                 MBH_JETS
  ************************************************************************/

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
    
  float MassConversion = (float) (dx*dx*dx * double(MassUnits));  //convert to g      
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

      CellsModified++;

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

    CellsModified = 0;
    return SUCCESS;
}

