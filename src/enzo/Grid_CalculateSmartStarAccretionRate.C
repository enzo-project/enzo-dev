/***********************************************************************
/
/  GRID CLASS (COMPUTE THE BONDI-HOYLE ACCRETION RATE)
/
/  written by: John Regan
/  date:       July 2016
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/
************************************************************************/
 
#include "preincludes.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "fortran.def"
#include "Grid.h"
#include "phys_constants.h"
#include "ActiveParticle_SmartStar.h"
#include "CosmologyParameters.h"

#define USEBOUNDEDNESS      0
#define TINY_NUMBER         1e-20
#define SMALL_NUMBER         1e-6
#define ACCRETION_LIMIT     1e-1
#define C_VISC              2.1e6
float bondi_alpha(float x);
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
 
float grid::CalculateSmartStarAccretionRate(ActiveParticleType* ThisParticle, 
					    FLOAT AccretionRadius, FLOAT *KernelRadius,
					    FLOAT *SumOfWeights)
{
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
  /* Find particle properties */
  FLOAT xparticle[3] = {
    ThisParticle->ReturnPosition()[0],
    ThisParticle->ReturnPosition()[1],
    ThisParticle->ReturnPosition()[2]
  };
  float vparticle[3] = {
    ThisParticle->ReturnVelocity()[0],
    ThisParticle->ReturnVelocity()[1],
    ThisParticle->ReturnVelocity()[2]
  };
  int cindex = (GridEndIndex[0] - GridStartIndex[0])/2 + GridStartIndex[0];
  int cgindex = GRIDINDEX_NOGHOST(cindex,cindex,cindex);
  float rhocell = 0.0, mcell = 0.0, CellVolume = 1.0;
  /* Calculate cell volume */
  for (int dim = 0; dim < GridRank; dim++)
  {
    CellVolume*=CellWidth[dim][0];
  }
  float mparticle = ThisParticle->ReturnMass()*CellVolume;
  ActiveParticleType_SmartStar* SS;
  SS = static_cast<ActiveParticleType_SmartStar*>(ThisParticle);
  SS->mass_in_accretion_sphere = 0.0;
  float eta_disk = SS->eta_disk;
  float WeightedSum = 0, AverageDensity = 0, RhoInfinity = 0.0;
  float AverageT=0, TotalGasMass = 0;
  float lambda_c = 0.25*exp(1.5);
  FLOAT radius2 = 0.0;
  float SmallRhoFac = 1e10, Weight = 0.0, SmallEFac = 10., SmEint = 0,  AccretedMomentum[3],
    vgas[3], etot, eint, ke,  maccreted, etotnew, rhonew, eintnew,
    kenew, mnew = 0;
  /* Find the Bondi-Hoyle radius */
  int size = this->GetGridSize();
  float *Temperature = new float[size]();
 
  this->ComputeTemperatureField(Temperature);
  /* Get indices in BaryonField for density, internal energy, thermal energy,
   * velocity */
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL)
  {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }

 
  /* Set the units. */
 
  float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1,
    TimeUnits = 1, VelocityUnits = 1,
    PressureUnits = 0, GEUnits = 0, VelUnits = 0;
  double MassUnits = 1;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
        ENZO_FAIL("Error in GetUnits.");
  }
  PressureUnits = DensityUnits * POW(LengthUnits,2) / POW(TimeUnits,2);
  GEUnits = POW(LengthUnits,2) / POW(TimeUnits, 2);
  VelUnits = LengthUnits/(TimeUnits*1e5); //convert to km/s
  MassUnits = DensityUnits * POW(LengthUnits,3);
 
  SmEint = max(
    SmallP * PressureUnits / ((Gamma - 1)*SmallRho),
    1.5 * kboltz * SmallT / (Mu*mh)
    ) / GEUnits;
  /* Estimate the relative velocity */
  float vInfinity = sqrt(pow(vparticle[0] - BaryonField[Vel1Num][cgindex],2) +
			 pow(vparticle[1] - BaryonField[Vel2Num][cgindex],2) +
			 pow(vparticle[2] - BaryonField[Vel3Num][cgindex],2));

  //float CellTemperature = Temperature[cgindex];
  float CellTemperature = FindAverageTemperatureinRegion(Temperature, xparticle, 2.0*AccretionRadius);
  if (JeansRefinementColdTemperature > 0)
    CellTemperature = JeansRefinementColdTemperature;
  float Gcode = GravConst*DensityUnits*TimeUnits*TimeUnits;
  float cInfinity = sqrt(Gamma * kboltz * CellTemperature / (Mu * mh)) /
    LengthUnits*TimeUnits;
  FLOAT BondiHoyleRadius = CalculateBondiHoyleRadius(mparticle, vparticle, Temperature); 
 
  //printf("%s:  BondiHoyleRadius = %e pc\n", __FUNCTION__,  BondiHoyleRadius*LengthUnits/pc_cm);
  //printf("%s:  AccretionRadius = %e pc\n", __FUNCTION__,  AccretionRadius*LengthUnits/pc_cm);
  /* Impose a kernel radius that regulates the weighting cells get as a function of radius */
  if (BondiHoyleRadius < CellWidth[0][0]/4.0) {  /* For BHs whose Bondi radius is not resolved */
    //printf("%s: Setting kernel radius to CellWidth, BH not resolved\n", __FUNCTION__);
    *KernelRadius = CellWidth[0][0]*2.0;
  }
  else { /*Accrete out to the BH radius */
    //printf("%s: Setting kernel radius to BondiHoyleRadius\n", __FUNCTION__);
    *KernelRadius = max(BondiHoyleRadius, AccretionRadius);
  }

  int numcells=0;
  for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      int index = GRIDINDEX_NOGHOST(GridStartIndex[0],j,k);
      for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
	radius2 =
	  POW((CellLeftEdge[0][i] + 0.5*CellWidth[0][i]) - xparticle[0],2) +
	  POW((CellLeftEdge[1][j] + 0.5*CellWidth[1][j]) - xparticle[1],2) +
	  POW((CellLeftEdge[2][k] + 0.5*CellWidth[2][k]) - xparticle[2],2);
	
	if ((AccretionRadius*AccretionRadius) > radius2) {
	  //printf("Grid:index = %d\n", index);
	  WeightedSum += BaryonField[DensNum][index]*
	    exp(-radius2/((*KernelRadius)*(*KernelRadius)));
	  (*SumOfWeights) += exp(-radius2/((*KernelRadius)*(*KernelRadius)));
	  AverageT += Temperature[index];
	  TotalGasMass += BaryonField[DensNum][index]*CellVolume;
	  numcells++;
	}
      }
    }
  }
  delete [] Temperature;
  Weight = 1.0/numcells;
  AverageT *= Weight;
  if(AverageT <= 0.0)
    AverageT = Temperature[cgindex];
  AverageDensity = WeightedSum / (*SumOfWeights);
#ifdef DEBUG_AP
  printf("AverageDensity = %g\n", AverageDensity);
#endif
  /* Calculate the accretion rate based on prescription specified */
  float AccretionRate = 0.0;
  SS->mass_in_accretion_sphere = TotalGasMass/CellVolume; //convert to density for consistency

  /* 
   * Traditional Bondi-Hoyle Prescription using the formalism
   * presented in Krumholtz et al. (2004)
   * The particle accretes according to the Bondi-Hoyle formula
   * with a density given by the average density within the 
   * accretion radius.
   */
  if(SmartStarAccretion == SPHERICAL_BONDI_HOYLE_FORMALISM ||
     SmartStarAccretion == SPHERICAL_BONDI_HOYLE_FORMALISM_WITH_VORTICITY) {
#ifdef DEBUG_AP
    printf("Doing SPHERICAL_BONDI_HOYLE_FORMALISM, SmartStarAccretion = %d\n", SmartStarAccretion);
#endif
    RhoInfinity = AverageDensity /
      bondi_alpha(1.2*CellWidth[0][0] / BondiHoyleRadius);

    /* Bondi Hoyle */
    AccretionRate = (4*pi*RhoInfinity*POW(BondiHoyleRadius,2)*
			 sqrt(POW(lambda_c*cInfinity,2) + POW(vInfinity,2)));
    /* Include Vorticity component if specified */
    if(SPHERICAL_BONDI_HOYLE_FORMALISM_WITH_VORTICITY == SmartStarAccretion) {
#ifdef DEBUG_AP
      printf("Doing SPHERICAL_BONDI_HOYLE_FORMALISM_WIDTH_VORTICITY, SmartStarAccretion = %d\n",
	     SmartStarAccretion);
#endif
      /* Include Vorticity Component */
      FLOAT vorticity[3] = {0.0, 0.0, 0.0};
      GetVorticityComponent(xparticle, vorticity);
      FLOAT vorticity_krumholz = sqrt(vorticity[0]*vorticity[0] +
				      vorticity[1]*vorticity[1] +
				      vorticity[2]*vorticity[2]) * BondiHoyleRadius / cInfinity;
      FLOAT fw = 1.0/(1 + POW(vorticity_krumholz, 0.9));
      FLOAT AccretionRate_Vorticity = 4.0*M_PI*RhoInfinity*BondiHoyleRadius*BondiHoyleRadius*cInfinity*0.34*fw;
      //float AccretionRate_Vorticity =
      //	GetVorticityComponent(BondiHoyleRadius, xparticle, cInfinity, RhoInfinity);
      //printf("AccretionRate_Vorticity = %g\n", AccretionRate_Vorticity);
      AccretionRate = POW(POW(AccretionRate, -2.0) + POW(AccretionRate_Vorticity, -2.0), -0.5);
    }
  }
  
  /* 
   * Viscous Angular momentum prescription from DeBuhr et al. (2010) 
   * The accretion rate onto the Black Hole is calculated by taking into 
   * account the angular momentum transport of gas on scales below the grid scale. 
   */
  if(SmartStarAccretion == VISCOUS_ANGULAR_MOMENTUM_TRANSPORT) {
#ifdef DEBUG_AP
    printf("Doing VISCOUS_ANGULAR_MOMENTUM_TRANSPORT, SmartStarAccretion = %d\n", SmartStarAccretion);
#endif
  float alpha = 0.1;
  float c_s = sqrt(Gamma * kboltz * AverageT / (Mu * mh)) /
    LengthUnits*TimeUnits;
  AccretionRate = 3.0 * M_PI * alpha * (c_s * c_s) * 
    (AverageDensity * AccretionRadius) /
    sqrt(Gcode * AverageDensity * 8.0 + 
	 Gcode * (TotalGasMass) / POW(CellWidth[0][0], 3.0));
  }

  /* 
   * Cen et al. (2011) 
   * A modified accretion rate based on the alpha-disk model of Shakura & Sunyaev (1973)
   * This model is designed to give reduced, physically motivated, accretion rates
   * explaining why BH feedback can not overwhelm a galaxy. 
   */
  if(SmartStarAccretion == ALPHA_DISK_CEN_2012) {
#ifdef DEBUG_AP
     printf("Doing ALPHA_DISK_CEN_2012, SmartStarAccretion = %d\n", SmartStarAccretion);
#endif
    AccretionRate = CenAccretionRate(AverageDensity, AccretionRadius,
				     xparticle, vparticle, mparticle);
  }

  /* 
   * Calculate the impact of angular momentum on the accretion flow (Rosas-Guevara (2013) 
   * The accretion rate is reduced when high angular momentum gas is located around the 
   * bondi-hoyle radius. 
   */
  if(SmartStarAccretion == ANGULAR_MOMENTUM_LIMITED_ACCRETION) {
#ifdef DEBUG_AP
    printf("Doing ANGULAR_MOMENTUM_LIMITED_ACCRETION, SmartStarAccretion = %d\n", SmartStarAccretion);
#endif
    float c_s = sqrt(Gamma * kboltz * AverageT / (Mu * mh))*TimeUnits/LengthUnits;
    float V_phi = CalculateCirculisationSpeed(Vel1Num, AccretionRadius, xparticle, vparticle);
    /* Bondi Hoyle */
    RhoInfinity = AverageDensity /
      bondi_alpha(1.2*CellWidth[0][0] / BondiHoyleRadius);
    AccretionRate = (4*pi*RhoInfinity*POW(BondiHoyleRadius,2)*
			 sqrt(POW(lambda_c*cInfinity,2) + POW(vInfinity,2)));
    //printf("AccretionRate = %g\n", AccretionRate);
    //printf("C_VISC = %g C_VISC^{1/3} = %g\t V_phi = %g\t c_s = %g km/s\n",
    //	   C_VISC, POW(C_VISC, 1.0/3.0), V_phi*VelUnits, c_s*VelUnits);
    float Vphi_CGS = V_phi*VelUnits;
    float soundspeed = c_s*VelUnits;
    if(POW(C_VISC, 1.0/3.0)*Vphi_CGS > soundspeed) {
      AccretionRate *= (POW(soundspeed/Vphi_CGS, 3.0)/C_VISC); //Rosas-Guevara et al (2013) Eq. 10
    }
  }
  
  /*
   * Take the mass flux across the accretion radius 
   * Find the radial velocity of the gas at the accretion radius and 
   * calculate the mass inflow rate. 
   *
   */
  if(SmartStarAccretion ==  CONVERGING_MASS_FLOW) {
#ifdef DEBUG_AP
    printf("Doing CONVERGING_MASS_FLOW, SmartStarAccretion = %d\n", SmartStarAccretion);
#endif
    AccretionRate = ConvergentMassFlow(DensNum, Vel1Num, AccretionRadius, xparticle, vparticle, 
				       mparticle, Gcode, GENum);
#ifdef DEBUG_AP
    printf("%s: Calculated (mass flux) accretion rate is %e Msolar/yr\n", __FUNCTION__, 
	   AccretionRate*3.154e7*MassUnits/(SolarMass*TimeUnits));
#endif
  }
  
  return AccretionRate;
}


int grid::GetVorticityComponent(FLOAT *pos, FLOAT *vorticity)
{
  if (GridRank != 3) 
    ENZO_FAIL("Devised only for three dimension.");
  int igrid[MAX_DIMENSION], dim, index, size = 1;
  float curl_x, curl_y, curl_z;
  float curl[3] = {0,0,0};
  int index_yp1, index_ym1, index_zp1, index_zm1;
  
  FLOAT dx = CellWidth[0][0], 
    dy = CellWidth[1][0], 
    dz = CellWidth[2][0];
  for (int dim = 0; dim < MAX_DIMENSION; dim++) {
    size *= GridDimension[dim];
    igrid[dim] = (int) ((pos[dim] - GridLeftEdge[dim]) / CellWidth[0][0]);
  }
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  /* Find fields: density, total energy, velocity1-3. */

  if (IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
					      Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }

  
  index_yp1 = 
    ((igrid[2] + GridStartIndex[2]) * GridDimension[1] + 
     igrid[1] + 1 + GridStartIndex[1]) * GridDimension[0] + 
    igrid[0] + GridStartIndex[0];
  index_ym1 = 
    ((igrid[2] + GridStartIndex[2]) * GridDimension[1] + 
     igrid[1] - 1 + GridStartIndex[1]) * GridDimension[0] + 
    igrid[0] + GridStartIndex[0];
  index_zp1 = 
    ((igrid[2] + 1 + GridStartIndex[2]) * GridDimension[1] + 
     igrid[1] + GridStartIndex[1]) * GridDimension[0] + 
    igrid[0] + GridStartIndex[0];
  index_zm1 = 
    ((igrid[2] - 1 + GridStartIndex[2]) * GridDimension[1] + 
     igrid[1] + GridStartIndex[1]) * GridDimension[0] + 
    igrid[0] + GridStartIndex[0];      
  curl_x = 
    float((0.5*(BaryonField[Vel3Num][index_yp1] - BaryonField[Vel3Num][index_ym1])/dy -
	   0.5*(BaryonField[Vel2Num][index_zp1] - BaryonField[Vel2Num][index_zm1])/dz));
  curl_y = 
    float((0.5*(BaryonField[Vel1Num][index_zp1] - BaryonField[Vel1Num][index_zm1])/dz -
	   0.5*(BaryonField[Vel3Num][index+1] - BaryonField[Vel3Num][index-1])/dx));
  curl_z = 
    float((0.5*(BaryonField[Vel2Num][index+1] - BaryonField[Vel2Num][index-1])/dx -
	   0.5*(BaryonField[Vel1Num][index_yp1] - BaryonField[Vel1Num][index_ym1])/dy));
  
  vorticity[0] = curl_x;
  vorticity[1] = curl_y;
  vorticity[2] = curl_z;
  return SUCCESS;
}

/*
 * Calculate the black hole accretion rate as per Cen (2012)
 */
float grid::CenAccretionRate(float density, FLOAT AccretionRadius,
			     FLOAT *pos, float *vel, float mparticle)
{
  int index_L = 0;
  int igrid[MAX_DIMENSION], dim = 0, size = 1;
  double alpha  = 0.1, kappa  = 0.4, lambda;
  double gas_angmom[] = {0.0, 0.0, 0.0}, total_gas_mass = 0.0, gas_mass = 0.0;
  float mdot = 0.0;
  FLOAT CellVolume = 1, BoxSize = 1, DensityConversion = 1, VelocityConversion = 1;
  FLOAT a = 1, dadt;
  FLOAT delx, dely, delz, velx, vely, velz;
  const double sigma_SB = 5.67e-5;

  for (int dim = 0; dim < MAX_DIMENSION; dim++) {
    size *= GridDimension[dim];
    igrid[dim] = (int) ((pos[dim] - GridLeftEdge[dim]) / CellWidth[0][0]);
  }
  /* Set the units. */
  float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1,
    TimeUnits = 1, VelocityUnits = 1,
    PressureUnits = 0, GEUnits = 0, VelUnits = 0;
  double MassUnits = 1;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.");
  }
  VelUnits = LengthUnits/(TimeUnits*1e5); //convert to km/s
  MassUnits = DensityUnits * POW(LengthUnits,3);
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  /* Find fields: density, total energy, velocity1-3. */

  if (IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
					      Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }
  
  float Gcode = GravConst*DensityUnits*TimeUnits*TimeUnits;
  for (dim = 0; dim < MAX_DIMENSION; dim++) 
    CellVolume *= CellWidth[0][0]; // in code units
  float MassConversion = CellVolume*MassUnits;
  /* Find angular momentum in 27 cells */
  int numcells = 0;
  float relvel = 0.0, sq_relvel = 0.0;
  for (int kk = -1; kk <= 1; kk++) {
    // relative position
    delz = (CellLeftEdge[2][igrid[2] + kk + GridStartIndex[2]] + 
	    0.5*CellWidth[2][0] - pos[2]); // in code units
	
    for (int jj = -1; jj <= 1; jj++) {
      dely = (CellLeftEdge[1][igrid[1] + jj + GridStartIndex[1]] + 
	      0.5*CellWidth[1][0] - pos[1]);  // in code units
	  
      for (int ii = -1; ii <= 1; ii++) {
	delx = (CellLeftEdge[0][igrid[0] + ii + GridStartIndex[0]] + 
		0.5*CellWidth[0][0] - pos[0]);  // in code units
	
	index_L = 
	  ((igrid[2] + kk + GridStartIndex[2]) * GridDimension[1] + 
	   igrid[1] + jj + GridStartIndex[1]) * GridDimension[0] + 
	  igrid[0] + ii + GridStartIndex[0];
	
	gas_mass = BaryonField[DensNum][index_L] * CellVolume; // in code mass
	    
	// relative velocity
	velx = (BaryonField[Vel1Num][index_L] - vel[0]); // in code velocity
	vely = (BaryonField[Vel2Num][index_L] - vel[1]);
	velz = (BaryonField[Vel3Num][index_L] - vel[2]);
	
	// store gas angular momentum in: M * L * L / T 
	gas_angmom[0] += gas_mass * ( vely*delz - velz*dely); 
	gas_angmom[1] += gas_mass * (-velx*delz + velz*delx);
	gas_angmom[2] += gas_mass * ( velx*dely - vely*delx);
	total_gas_mass += gas_mass;
	relvel += sqrt(velx*velx + vely*vely + velz*velz);
	sq_relvel += velx*velx + vely*vely + velz*velz;
	numcells++;
      }
    }
  }
  float mean = relvel/numcells;
  float variance = sq_relvel/numcells - mean*mean;
  float sigmavel = sqrt(variance);
  // specific gas angular momentum in: L * L / T
  for (dim = 0; dim < MAX_DIMENSION; dim++)
    gas_angmom[dim] /= total_gas_mass;

  // the specific angular momentum in Keplerian orbit in: 
  //


  double Keplerian_angmom = Gcode * mparticle /sigmavel;
  
  // now lambda = the ratio of specific angular momentum (the real vs. the Keplerian orbit)
  lambda = fabs(gas_angmom[2]) / Keplerian_angmom;
  

  fprintf(stdout, "Star::CalculateMassAccretion: method=7: lambda = %g, gas_angmom[2] = %g, "
	  "Keplerian_angmom = %g\n", lambda, gas_angmom[2], Keplerian_angmom);

  // Calculate accretion rate in Msun/s
  // mdot = 3pi*alpha * c_s^2 * Sigma / Omega
  //      = 3^(4/3)/4^(1/3) * pi * (3*k_b/m_h)^(4/3) * (3*kappa/(16*sigma_SB*G))^(1/3) * alpha^(4/3) *
  //        M^(-1/3) * Sigma^(5/3) * R * lambda^(-4/3)
  /* This needs to be cgs now */
  mdot = pow(3.0, 4.0/3.0)/pow(4.0, 1.0/3.0) * M_PI * pow(3.0*kboltz/mh, 4.0/3.0) *
    pow(3.0*kappa/16.0/sigma_SB/GravConst, 1.0/3.0) * pow(alpha, 4.0/3.0)
    * pow(mparticle * MassUnits, -1.0/3.0) *
    pow(density * DensityUnits * CellWidth[0][0]*LengthUnits, 5.0/3.0) *
    AccretionRadius*LengthUnits * pow(lambda, -7.0/3.0);
  
  //    if (mdot > 0.0) {
  // 	fprintf(stdout, "BH Accretion[%"ISYM"]: time = %"FSYM", mdot = %"GSYM" Msun/yr, "
  // 		"M_BH = %lf Msun, rho = %"GSYM" g/cm3, c_s = %"GSYM" cm/s, T = %"GSYM" K\n",
  // 		Identifier, time, mdot*yr, Mass, density*DensityUnits, c_s, temperature[index]);
  //	this->PrintInfo();  
  //    }
  return mdot*TimeUnits/MassUnits;  //units = code units
}		   


/* 
 * If requested, adandon Bondi-Hoyle formalism and
 * just use the actual converging mass into the cell to calculate mdot 
 * We just calculate the mass flux through the accretion radius of the black hole
 * This follows prescriptions used by other authors e.g. Bleuer et al. 2015
 */
float grid::ConvergentMassFlow(int DensNum, int Vel1Num, FLOAT AccretionRadius,
			       FLOAT *pos, float *vel, float SSmass, float Gcode, int GENum)
{
  int numincells = 0, numoutcells = 0;
  float mdot = 0.0;
  float epsilon = AccretionRadius*0.1;
  float *density = BaryonField[DensNum];
  float *gasvelx = BaryonField[Vel1Num];
  float *gasvely = BaryonField[Vel1Num++];
  float *gasvelz = BaryonField[Vel1Num++];
  float div = 0.0, divx = 0.0, divy = 0.0, divz = 0.0;
  const int offset[] = {1, GridDimension[0], GridDimension[0]*GridDimension[1]};
  for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      int index = GRIDINDEX_NOGHOST(GridStartIndex[0],j,k);
      for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
	/* Want to calculate the radial velocity velocity vector
	 * so that a positive radial velocity means gas travelling to the 
	 * black hole.
	 * In this case the equation is:
	 * r_1 = BH position
	 * r_2 = cell position
	 * Relpos = r_1 - r_2
	 * and the same for the velocities.
	 * The radial velocity is the dot product i.e.
	 * (v_1 - v_2).((r_1 - r_2)/|r_1 - r_2|)
	 */

	/* r-1 - r_2 */
	FLOAT relx = pos[0] - (CellLeftEdge[0][i] + 0.5*CellWidth[0][i]);
	FLOAT rely = pos[1] - (CellLeftEdge[1][j] + 0.5*CellWidth[1][j]);
	FLOAT relz = pos[2] - (CellLeftEdge[2][k] + 0.5*CellWidth[2][k]);
	FLOAT radius2 = POW(relx,2) + POW(rely,2) + POW(relz,2);
	if ((AccretionRadius*AccretionRadius) > radius2 &&
	    ((AccretionRadius-epsilon)*(AccretionRadius-epsilon)) < radius2) {
	  FLOAT relposmag = sqrt(radius2);
	  FLOAT relpos[3] = { relx/relposmag, rely/relposmag, relz/relposmag};
	  FLOAT vrel[3] = {vel[0] - gasvelx[index],
	  		   vel[1] - gasvely[index],
	  		   vel[2] - gasvelz[index]};
	  FLOAT radialvelocity = (vrel[0]*relpos[0] +
				  vrel[1]*relpos[1] +
				  vrel[2]*relpos[2] );
	  float accrate = 0.0; 
	  if(radialvelocity < 0) {
#if USEBOUNDEDNESS
	    float ke = pow(gasvelx[index], 2.0) + pow(gasvely[index], 2.0) + pow(gasvelz[index], 2.0);
	    float te = BaryonField[GENum][index];
	    FLOAT dist = sqrt(radius2);
	    float ge = Gcode*SSmass/dist;
	    if(ke+te-ge<0) { /*Only add if we are bound */
	      continue;
	    }
#endif
	    numincells++;
	    accrate = density[index]*relposmag*relposmag*radialvelocity;
	  }
	  else
	    numoutcells++;
	  //mdot += density[index]*relposmag*relposmag*radialvelocity;
	  //float accrate = density[index] * pow(AccretionRadius, 2.0) * div;
	  //printf("%s: Accretion Rate = %g\n", __FUNCTION__, accrate*4.0*M_PI);
	  mdot += accrate;
	}
      }
    }
  }
  // mdot = -4*pi*rho*R^2*V_radial
  mdot = fabs(4*M_PI*mdot); //return the accretion rate as a positive quantity
#ifdef DEBUG_AP
  printf("%s: Num InFlow cells = %d\t Num OutflowCells = %d\t mdot = %e\n", __FUNCTION__, numincells,
	 numoutcells, mdot);
#endif
  return mdot;
}

/*
 * Calculate and return the average tangential velocity at the 
 * accretion radius.
 */
float grid::CalculateCirculisationSpeed(int Vel1Num, FLOAT AccretionRadius,
					FLOAT *pos, float *vel)
{
  int numcells = 0;
  float tang_vel = 0.0;
  FLOAT total_vorticity[3] = {0.0, 0.0, 0.0};
  float epsilon = AccretionRadius*0.02;
  float *gasvelx = BaryonField[Vel1Num];
  float *gasvely = BaryonField[Vel1Num++];
  float *gasvelz = BaryonField[Vel1Num++];
  for (int k = GridStartIndex[2]; k <= GridEndIndex[2]; k++) {
    for (int j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      int index = GRIDINDEX_NOGHOST(GridStartIndex[0],j,k);
      for (int i = GridStartIndex[0]; i <= GridEndIndex[0]; i++, index++) {
	FLOAT relx = (CellLeftEdge[0][i] + 0.5*CellWidth[0][i]) - pos[0];
	FLOAT rely = (CellLeftEdge[1][j] + 0.5*CellWidth[1][j]) - pos[1];
	FLOAT relz = (CellLeftEdge[2][k] + 0.5*CellWidth[2][k]) - pos[2];
	FLOAT radius2 = POW(relx,2) + POW(rely,2) + POW(relz,2);
	if ((AccretionRadius*AccretionRadius) > radius2 &&
	    ((AccretionRadius-epsilon)*(AccretionRadius-epsilon)) < radius2) {
	  FLOAT relposmag = sqrt(relx*relx + rely*rely + relz*relz);
	  FLOAT relpos[3] = { relx/relposmag, rely/relposmag, relz/relposmag};
	  FLOAT vrel[3] = {vel[0] - gasvelx[index],
			   vel[1] - gasvely[index],
			   vel[2] - gasvelz[index]};
	  FLOAT radialvelocity = fabs(vrel[0]*relpos[0] +
					   vrel[1]*relpos[1] +
					   vrel[2]*relpos[2] );
	  /* Tangential velocity = (v x r)/|r^2| */
	  FLOAT cellpos[3] = {pos[0] + relx, pos[1] + rely, pos[2] + relz};
	  FLOAT vorticity[3] = {0.0, 0.0, 0.0};
	  GetVorticityComponent(cellpos, vorticity);
	  //tang_vel += sqrt(vrel[0]*vrel[0] + vrel[1]*vrel[1] + vrel[2]*vrel[2] -
	  //		  radialvelocity*radialvelocity);
	  total_vorticity[0] += vorticity[0];total_vorticity[1] += vorticity[1];total_vorticity[2] += vorticity[2]; 
	  numcells++;
	  //printf("%s: numcells = %d\t relpos = (%g, %g, %g)\n", __FUNCTION__, numcells, relpos[0], relpos[1], relpos[2]);
	  //printf("%s: vorticty = (%g, %g, %g)\n", __FUNCTION__, vorticity[0], vorticity[1], vorticity[2]);
	}
      }
    }
  }
  total_vorticity[0] /= numcells;
  total_vorticity[1] /= numcells;
  total_vorticity[2] /= numcells;
  //return tang_vel/numcells;
  FLOAT totalvortmag = sqrt(total_vorticity[0]*total_vorticity[0] + total_vorticity[1]*total_vorticity[1] + total_vorticity[2]*total_vorticity[2]);
  
  return totalvortmag*AccretionRadius;
}


FLOAT grid::CalculateBondiHoyleRadius(float mparticle, float *vparticle, float *Temperature)
{

  int cindex = (GridEndIndex[0] - GridStartIndex[0])/2 + GridStartIndex[0];
  int cgindex = GRIDINDEX_NOGHOST(cindex,cindex,cindex);
  float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1,
    TimeUnits = 1, VelocityUnits = 1, 
    PressureUnits = 0, GEUnits = 0, VelUnits = 0;
  double MassUnits = 1;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.");
  }
  PressureUnits = DensityUnits * POW(LengthUnits,2) / POW(TimeUnits,2);
  GEUnits = POW(LengthUnits,2) / POW(TimeUnits, 2);
  VelUnits = LengthUnits/(TimeUnits*1e5); //convert to km/s
  MassUnits = DensityUnits * POW(LengthUnits,3);
  /* Get indices in BaryonField for density, internal energy, thermal energy,
   * velocity */
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL)
  {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }

 
  /* Estimate the relative velocity */
  float vInfinity = sqrt(pow(vparticle[0] - BaryonField[Vel1Num][cgindex],2) +
			 pow(vparticle[1] - BaryonField[Vel2Num][cgindex],2) +
			 pow(vparticle[2] - BaryonField[Vel3Num][cgindex],2));

  float CellTemperature = Temperature[cgindex];
  if (JeansRefinementColdTemperature > 0)
    CellTemperature = JeansRefinementColdTemperature;

  float cInfinity = sqrt(Gamma * kboltz * CellTemperature / (Mu * mh)) /
    LengthUnits*TimeUnits;
  float Gcode = GravConst*DensityUnits*TimeUnits*TimeUnits;
#ifdef DEBUG_AP
  printf("%s: vInfinity = %f km/s\n", __FUNCTION__,  (vInfinity*LengthUnits/TimeUnits)/1e5);
  printf("%s: cInfinity = %f km/s\n", __FUNCTION__,  (cInfinity*LengthUnits/TimeUnits)/1e5);
  printf("%s: CellTemperature = %f K\n", __FUNCTION__, CellTemperature);
  printf("%s: Celllength = %e pc\n", __FUNCTION__, CellWidth[0][0]*LengthUnits/pc_cm);
#endif
  return Gcode*mparticle/
    (pow(vInfinity,2) + pow(cInfinity,2));
}
