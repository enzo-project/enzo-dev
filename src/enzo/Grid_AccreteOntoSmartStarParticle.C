/***********************************************************************
/
/  (Calculate the accretion rate and subtract accreted mass from the
/   grid.)
/
/  written by: Nathan Goldbaum
/  date:       April 2012
/
/  note:       Equation numbers refer to Krumholz McKee & Klein (2004)
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
#include "ActiveParticle_SmartStar.h"
#define UPDATE_SS_VELOCITY 1
#define NO_DEBUG_AP
#define ACCRETE_DEBUG 0
#define NO_ACCRETION 0

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
int grid::AccreteOntoSmartStarParticle(
  ActiveParticleType* ThisParticle,FLOAT AccretionRadius,
  float* AccretionRate)
{

  /* Return if this doesn't involve us */
  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;
  ActiveParticleType_SmartStar* SS;
  SS = static_cast<ActiveParticleType_SmartStar*>(ThisParticle);
  /* Check whether the cube that circumscribes the accretion zone intersects
   * with this grid */
  FLOAT xsink = ThisParticle->ReturnPosition()[0];
  FLOAT ysink = ThisParticle->ReturnPosition()[1];
  FLOAT zsink = ThisParticle->ReturnPosition()[2];
  FLOAT dx = CellWidth[0][0];
  if ((GridLeftEdge[0] > xsink+AccretionRadius) ||
      (GridLeftEdge[1] > ysink+AccretionRadius) ||
      (GridLeftEdge[2] > zsink+AccretionRadius) ||
      (GridRightEdge[0] < xsink-AccretionRadius) ||
      (GridRightEdge[1] < ysink-AccretionRadius) ||
      (GridRightEdge[2] < zsink-AccretionRadius))
    return SUCCESS;
 
  /* Delcare and initialize local variables */

  float delta_vpart[3] = {0.0, 0.0, 0.0};
  float AccretedMass = 0;  

  /* Set the units. */
 
  float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1,
    TimeUnits = 1, VelocityUnits = 1,
    PressureUnits = 0, GEUnits = 0, VelUnits = 0;
  double MassUnits = 1;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
        ENZO_FAIL("Error in GetUnits.");
  }
  float CellVolume = 1.0;
  /* Calculate cell volume */
  for (int dim = 0; dim < GridRank; dim++)
  {
    CellVolume*=CellWidth[dim][0];
  }
  MassUnits = DensityUnits * POW(LengthUnits,3);
  float mparticle = SS->ReturnMass()*CellVolume; //code mass
  float MassConversion = (float) (dx*dx*dx * double(MassUnits));  //convert to g
  FLOAT KernelRadius = 0.0, SumOfWeights = 0.0; /*Required for weighting cells for accretion */
  *AccretionRate = CalculateSmartStarAccretionRate(ThisParticle, AccretionRadius,
						   &KernelRadius, &SumOfWeights);

#if NO_ACCRETION
  *AccretionRate = 0.0;
#endif
  /*
   * 1. Calculate the accretion rate
   * 2a. Calculate how much goes into the black hole: Mdot_BH = AccRate*(1 - eta_disk) ~ 0.9*AccRate
   * 2b. Calculate how much goes into the Jet: Mdot_Jet = Mdot_BH*beta_jet ~ 90*AccRate (for Beta_Jet = 100)
   * 3. Remove the MDot_BH from the grid and also take opportunity to calculate total mass available
   *    within the accretion sphere. 
   * 4. The Thermal energy (or radiation) is composed of eta_disk*AccRate ~ 0.1*AccRate  
   */

   float MaxAccretionRate = huge_number;
   /* 
    * If we want to cap the accretion rate at the Eddington rate the most sensible thing
    * is to restrict the amount of mass that can be extracted during an accretion event.
    * Physically I don't like this and you should consider why you want to do this. 
    * It can be very useful for testing however. 
    */
   if(BH == SS->ParticleClass && SmartStarEddingtonCap) {
     float mdot_edd = 4.0 * PI * GravConst * mparticle*MassUnits * mh / 
       (clight * SS->eta_disk * sigma_thompson);
     float accrate_cgs = (*AccretionRate*MassUnits/TimeUnits);
     // printf("%s: So the accretion rate is %e Msolar/yr and the EddRate is %e Msolar/yr?\n",
     // 	    __FUNCTION__, accrate_cgs*3.154e7/SolarMass, mdot_edd*3.154e7/SolarMass);
     if(accrate_cgs > mdot_edd)
       MaxAccretionRate = mdot_edd*TimeUnits/MassUnits;
     else
       MaxAccretionRate = *AccretionRate; 
   }
  /* 
   * Remove mass from the grid in a sensible way using a kernel weighting 
   * NB: This will almost certainly cause a difference to the 
   * calculated accretion rate above since we apply a weighting kernel.
   * The weighting kernel takes account of the Bondi-Hoyle radius and whether we
   * are resolving that radius or not. The mass actually accreted (i.e. removed from 
   * the grid) can then be much less than found from the mass flux for example but 
   * is closer to what the black hole would actually accrete. 
   */
 
  RemoveMassFromGrid(ThisParticle,AccretionRadius, *AccretionRate,
		     &AccretedMass, delta_vpart,
		     KernelRadius, SumOfWeights, MaxAccretionRate);
#if  ACCRETE_DEBUG
  printf("%s: DeltaV = %e %e %e\n", __FUNCTION__,
	 delta_vpart[0], delta_vpart[1], delta_vpart[2]);
  printf("%s: AccretedMass = %e\n", __FUNCTION__, AccretedMass);
#endif
  float *Vel = SS->vel;
  /* 
   * I don't think we should update the velocity of the particle
   * This breaks conservation of momentum but the momentum is in the subgrid. The
   * velocity of the accreted gas is placed in the accretion disk and 
   * radiated back. The black hole itself probably doesn't move in nature. It
   * should not move here either
   */
#if UPDATE_SS_VELOCITY
  float NewVelocity[3] =
    {
    (Vel[0]+delta_vpart[0]),
    (Vel[1]+delta_vpart[1]),
    (Vel[2]+delta_vpart[2])
    };

  ThisParticle->SetVelocity(NewVelocity);
#endif
  /* 
   * This value is the actual accretion rate onto the SmartStar. It was initially
   * calculated according to some prescription (e.g. Bondi-Hoyle) and then
   * a kernel weighting applied to extract gas from the grid.
   */
  *AccretionRate = AccretedMass/this->dtFixed; //in units of density/time
 
  /*
   * Using the actual accretion rate we can update how much of the gas makes it onto the
   * smart star and of the mass that is available for feedback how much, if any, ends 
   * up in Jets. 
   */
  if(BH == SS->ParticleClass && SmartStarBHJetFeedback) {
    *AccretionRate *= (1.0 - SS->eta_disk);
    /* 
     * Now this is the complicated/clever bit. 
     * The MassEjected will not be a fraction of the accreted mass but a fraction of the available mass
     * from the surrounding cells. The mass in the surrounding cells is mass_in_accretion_sphere
     * beta_jet = Jet Mass Loading
     * eta_jet = Jet efficiency (comes from GRMHD simulations)
     * epsilon_deltat = subgrids the timestep and effectively restricts jet emission
     * to some small part of the timestep which we do not resolve.
     * Mdot_Jet = beta_jet * (1 - eta_disk) * mdot * epsilon_deltat
     * MassEjected = Mdot_Jet * dt 
     */
    float BHMass =  SS->ReturnMass()*MassConversion/SolarMass; //In solar masses
    float eddrate = 4*M_PI*GravConst*BHMass*SolarMass*mh/(SS->eta_disk*clight*sigma_thompson); // g/s
    eddrate = eddrate*3.154e7/SolarMass; //in Msolar/yr
    float accrate_msolar = (*AccretionRate)*3.154e7*MassConversion/(SolarMass*TimeUnits); //Msolar/yr
    /* eta_jet comes from the fit given in Sadowski et al. (2016) */
    float eta_jet = 1.3*POW(SmartStarSpin, 2.0);
    SS->beta_jet = (1.0/SmartStarJetVelocity)*(1.0/SmartStarJetVelocity)*2.0*eta_jet;
    float mdot_total = SS->mass_in_accretion_sphere/this->dtFixed; //This is the total mass available for accretion inside accretion sphere
    float accretion_ratio = mdot_total/(*AccretionRate);
    SS->epsilon_deltat = min(1.0, accretion_ratio*(1.0/(1.0 - SS->eta_disk))*(1.0/(1.0 + SS->beta_jet)));
    SS->MassToBeEjected = 0.0;
#if  ACCRETE_DEBUG
    printf("%s: Eddrate = %e Msolar/yr AccRate = %e Msolar/yr\t Ratio = %f\n", __FUNCTION__,
	   eddrate, accrate_msolar, accrate_msolar/eddrate);
#endif
    if(accrate_msolar > eddrate) {
      
      *AccretionRate *= SS->epsilon_deltat;
      SS->MassToBeEjected = SS->beta_jet*(*AccretionRate)*this->dtFixed*dx*dx*dx; //Code Mass
      *AccretionRate *= (1 - eta_jet);
#if  ACCRETE_DEBUG
      printf("%s: eta_disk = %f\n", __FUNCTION__, SS->eta_disk);
      printf("%s: eta_jet = %f\t beta_jet = %f\t ParticleClass = %d\n", __FUNCTION__, eta_jet,  SS->beta_jet,SS->ParticleClass );
      printf("%s: Mass in surrounding sphere = %e Msolar\n", __FUNCTION__, SS->mass_in_accretion_sphere*MassConversion/SolarMass);
      printf("%s: Macc = %e Msolar\t Mjet = %e Msolar\n", __FUNCTION__,
	     (*AccretionRate)*this->dtFixed*MassConversion/SolarMass,
	     SS->MassToBeEjected*MassUnits/SolarMass);
      printf("%s: accretion_ratio = %f\n", __FUNCTION__, accretion_ratio);
      printf("%s: epsilon_deltat = %e\n", __FUNCTION__, SS->epsilon_deltat);
      printf("%s: Mass to be ejected = %f Msolar (%e code)\n",  __FUNCTION__, SS->MassToBeEjected*MassUnits/SolarMass, SS->MassToBeEjected);
      printf("%s: Updated Accretion rate = %e Msolar/yr\n", __FUNCTION__,
	     (*AccretionRate)*3.154e7*MassConversion/(SolarMass*TimeUnits));
      printf("%s: No update would be = %e Msolar/yr\n", __FUNCTION__,
	     (*AccretionRate/SS->epsilon_deltat)*3.154e7*MassConversion/(SolarMass*TimeUnits));
#endif
    }
    else {
#if  ACCRETE_DEBUG
      printf("%s: AccrateionRateRatio = %f. No jets this time\n", __FUNCTION__,
	     accrate_msolar/eddrate);
#endif
      ;
    }
  }
  
  AccretedMass = (*AccretionRate)*this->dtFixed;
  ThisParticle->AddMass(AccretedMass);
#if  ACCRETE_DEBUG
  printf("%s: AccretedMass = %e Msolar\n", __FUNCTION__, AccretedMass*MassConversion/SolarMass);
  printf("%s: PrevMass = %e Msolar\t NewMass = %e Msolar\n", __FUNCTION__,
	 (SS->ReturnMass() - AccretedMass)*MassConversion/SolarMass,
	 SS->ReturnMass()*MassConversion/SolarMass);
#endif
  return SUCCESS;
}


