/***********************************************************************
/
/ Calculate the SED energy bins for POPIII & SMS particles based on particle
/ mass and age
/ Author: John Regan
/ Date: Mid 2017
/
************************************************************************/
#include "ActiveParticle_SmartStar.h"
#include "phys_constants.h"

float PopIIIEnergyBins[NUMRADIATIONBINS] = {2.0, 12.8, 14.0, 25.0, 200.0};
float PopIIISEDFracBins[NUMRADIATIONBINS] = {0.3261, 0.1073, 0.3686, 0.1965, 0.0};
//float PopIIISEDFracBins[NUMRADIATIONBINS] = {0.3261, 0.1073, 0.0, 0.0, 0.0};
float SMSEnergyBins[NUMRADIATIONBINS] = {2.0, 12.8, 14.0, 25.0, 200.0};
float SMSSEDFracBins[NUMRADIATIONBINS] = {0.9992, 0.0008, 0.0, 0.0, 0.0};
float BHEnergyBins[NUMRADIATIONBINS] = {2.0, 12.8, 19.1, 217.3, 5190.0};
float BHSEDFracBins[NUMRADIATIONBINS] = {0.0, 0.0, 0.0, 0.0, 0.0};
#define RAMPAGE 1e5
#define STELLAR 0
/* BHArray is the photonspersecond, IR Radiation, LW Radiation, Ionising Radiation, Soft X-ray, Hard X-ray */
/* 
 * The array covers black hole masses from 1 Msolar to 1e9 Msolar 
 * and accretion rates from 1e-6 to 1e3 Msolar/yr 
 */
double BHArray[100][6] = {{8.536889e+47, 1.444555e-04, 3.879133e-05, 2.669210e-03, 1.796218e-01, 7.815097e-01},
			 {6.886130e+48, 4.986472e-05, 1.346838e-05, 9.516688e-04, 1.563754e-01, 7.995463e-01},
			 {5.960135e+49, 1.890098e-05, 5.121741e-06, 3.673078e-04, 1.541222e-01, 7.943261e-01},
			 {5.303473e+50, 7.825415e-06, 2.124380e-06, 1.536224e-04, 1.622892e-01, 7.764801e-01},
			 {4.794634e+51, 3.611429e-06, 9.814051e-07, 7.130109e-05, 1.749921e-01, 7.494527e-01},
			 {4.446445e+52, 1.970246e-06, 5.357218e-07, 3.902338e-05, 1.868796e-01, 7.186502e-01},
			 {4.251238e+53, 1.315918e-06, 3.579219e-07, 2.611038e-05, 1.947393e-01, 6.941544e-01},
			 {4.151599e+54, 1.035488e-06, 2.816977e-07, 2.056684e-05, 1.991070e-01, 6.790028e-01},
			 {4.100720e+55, 9.039782e-07, 2.459464e-07, 1.796498e-05, 2.014349e-01, 6.704891e-01},
			 {4.074003e+56, 8.377752e-07, 2.279476e-07, 1.665463e-05, 2.026844e-01, 6.658005e-01},
			 {1.984489e+48, 1.464404e-03, 3.819372e-04, 2.300907e-02, 3.401724e-01, 5.979580e-01},
			 {1.213797e+49, 4.550169e-04, 1.209228e-04, 7.934018e-03, 2.403197e-01, 7.184667e-01},
			 {8.536889e+49, 1.444555e-04, 3.879133e-05, 2.669210e-03, 1.796218e-01, 7.815097e-01},
			 {6.886130e+50, 4.986472e-05, 1.346838e-05, 9.516688e-04, 1.563754e-01, 7.995463e-01},
			 {5.960135e+51, 1.890098e-05, 5.121741e-06, 3.673078e-04, 1.541222e-01, 7.943261e-01},
			 {5.303473e+52, 7.825415e-06, 2.124380e-06, 1.536224e-04, 1.622892e-01, 7.764801e-01},
			 {4.794634e+53, 3.611429e-06, 9.814051e-07, 7.130109e-05, 1.749921e-01, 7.494527e-01},
			 {4.446445e+54, 1.970246e-06, 5.357218e-07, 3.902338e-05, 1.868796e-01, 7.186502e-01},
			 {4.251238e+55, 1.315918e-06, 3.579219e-07, 2.611038e-05, 1.947393e-01, 6.941544e-01},
			 {4.151599e+56, 1.035488e-06, 2.816977e-07, 2.056684e-05, 1.991070e-01, 6.790028e-01},
			 {5.966964e+48, 1.353020e-02, 3.202910e-03, 1.284992e-01, 5.197274e-01, 2.700725e-01},
			 {3.416695e+49, 4.562143e-03, 1.149944e-03, 5.953758e-02, 4.466020e-01, 4.399171e-01},
			 {1.984489e+50, 1.464404e-03, 3.819372e-04, 2.300907e-02, 3.401724e-01, 5.979580e-01},
			 {1.213797e+51, 4.550169e-04, 1.209228e-04, 7.934018e-03, 2.403197e-01, 7.184667e-01},
			 {8.536889e+51, 1.444555e-04, 3.879133e-05, 2.669210e-03, 1.796218e-01, 7.815097e-01},
			 {6.886130e+52, 4.986472e-05, 1.346838e-05, 9.516688e-04, 1.563754e-01, 7.995463e-01},
			 {5.960135e+53, 1.890098e-05, 5.121741e-06, 3.673078e-04, 1.541222e-01, 7.943261e-01},
			 {5.303473e+54, 7.825415e-06, 2.124380e-06, 1.536224e-04, 1.622892e-01, 7.764801e-01},
			 {4.794634e+55, 3.611429e-06, 9.814051e-07, 7.130109e-05, 1.749921e-01, 7.494527e-01},
			 {4.446445e+56, 1.970246e-06, 5.357218e-07, 3.902338e-05, 1.868796e-01, 7.186502e-01},
			 {1.855450e+49, 9.198589e-02, 1.572652e-02, 2.954364e-01, 4.654125e-01, 3.036475e-02},
			 {1.050058e+50, 3.734607e-02, 7.870513e-03, 2.185504e-01, 5.319855e-01, 1.209364e-01},
			 {5.966964e+50, 1.353020e-02, 3.202910e-03, 1.284992e-01, 5.197274e-01, 2.700725e-01},
			 {3.416695e+51, 4.562143e-03, 1.149944e-03, 5.953758e-02, 4.466020e-01, 4.399171e-01},
			 {1.984489e+52, 1.464404e-03, 3.819372e-04, 2.300907e-02, 3.401724e-01, 5.979580e-01},
			 {1.213797e+53, 4.550169e-04, 1.209228e-04, 7.934018e-03, 2.403197e-01, 7.184667e-01},
			 {8.536889e+53, 1.444555e-04, 3.879133e-05, 2.669210e-03, 1.796218e-01, 7.815097e-01},
			 {6.886130e+54, 4.986472e-05, 1.346838e-05, 9.516688e-04, 1.563754e-01, 7.995463e-01},
			 {5.960135e+55, 1.890098e-05, 5.121741e-06, 3.673078e-04, 1.541222e-01, 7.943261e-01},
			 {5.303473e+56, 7.825415e-06, 2.124380e-06, 1.536224e-04, 1.622892e-01, 7.764801e-01},
			 {5.799910e+49, 3.298557e-01, 2.808534e-02, 3.445588e-01, 1.643551e-01, 8.516466e-04},
			 {3.283307e+50, 1.915271e-01, 2.361579e-02, 3.387841e-01, 3.232404e-01, 3.448759e-03},
			 {1.855450e+51, 9.198589e-02, 1.572652e-02, 2.954364e-01, 4.654125e-01, 3.036475e-02},
			 {1.050058e+52, 3.734607e-02, 7.870513e-03, 2.185504e-01, 5.319855e-01, 1.209364e-01},
			 {5.966964e+52, 1.353020e-02, 3.202910e-03, 1.284992e-01, 5.197274e-01, 2.700725e-01},
			 {3.416695e+53, 4.562143e-03, 1.149944e-03, 5.953758e-02, 4.466020e-01, 4.399171e-01},
			 {1.984489e+54, 1.464404e-03, 3.819372e-04, 2.300907e-02, 3.401724e-01, 5.979580e-01},
			 {1.213797e+55, 4.550169e-04, 1.209228e-04, 7.934018e-03, 2.403197e-01, 7.184667e-01},
			 {8.536889e+55, 1.444555e-04, 3.879133e-05, 2.669210e-03, 1.796218e-01, 7.815097e-01},
			 {6.886130e+56, 4.986472e-05, 1.346838e-05, 9.516688e-04, 1.563754e-01, 7.995463e-01},
			 {1.746463e+50, 6.558310e-01, 2.931956e-02, 1.960125e-01, 6.487986e-03, 2.781257e-04},
			 {1.016622e+51, 4.899245e-01, 3.021142e-02, 2.987737e-01, 5.118229e-02, 4.777958e-04},
			 {5.799910e+51, 3.298557e-01, 2.808534e-02, 3.445588e-01, 1.643551e-01, 8.516466e-04},
			 {3.283307e+52, 1.915271e-01, 2.361579e-02, 3.387841e-01, 3.232404e-01, 3.448759e-03},
			 {1.855450e+53, 9.198589e-02, 1.572652e-02, 2.954364e-01, 4.654125e-01, 3.036475e-02},
			 {1.050058e+54, 3.734607e-02, 7.870513e-03, 2.185504e-01, 5.319855e-01, 1.209364e-01},
			 {5.966964e+54, 1.353020e-02, 3.202910e-03, 1.284992e-01, 5.197274e-01, 2.700725e-01},
			 {3.416695e+55, 4.562143e-03, 1.149944e-03, 5.953758e-02, 4.466020e-01, 4.399171e-01},
			 {1.984489e+56, 1.464404e-03, 3.819372e-04, 2.300907e-02, 3.401724e-01, 5.979580e-01},
			 {1.213797e+57, 4.550169e-04, 1.209228e-04, 7.934018e-03, 2.403197e-01, 7.184667e-01},
			 {4.503815e+50, 9.186506e-01, 1.054325e-02, 1.742493e-02, 1.826442e-04, 1.078499e-04},
			 {2.883229e+51, 8.072469e-01, 2.280500e-02, 8.348202e-02, 4.014262e-04, 1.684695e-04},
			 {1.746463e+52, 6.558310e-01, 2.931956e-02, 1.960125e-01, 6.487986e-03, 2.781257e-04},
			 {1.016622e+53, 4.899245e-01, 3.021142e-02, 2.987737e-01, 5.118229e-02, 4.777958e-04},
			 {5.799910e+53, 3.298557e-01, 2.808534e-02, 3.445588e-01, 1.643551e-01, 8.516466e-04},
			 {3.283307e+54, 1.915271e-01, 2.361579e-02, 3.387841e-01, 3.232404e-01, 3.448759e-03},
			 {1.855450e+55, 9.198589e-02, 1.572652e-02, 2.954364e-01, 4.654125e-01, 3.036475e-02},
			 {1.050058e+56, 3.734607e-02, 7.870513e-03, 2.185504e-01, 5.319855e-01, 1.209364e-01},
			 {5.966964e+56, 1.353020e-02, 3.202910e-03, 1.284992e-01, 5.197274e-01, 2.700725e-01},
			 {3.416695e+57, 4.562143e-03, 1.149944e-03, 5.953758e-02, 4.466020e-01, 4.399171e-01},
			 {9.418773e+50, 9.943227e-01, 2.228705e-05, 2.769577e-06, 8.729875e-05, 5.157107e-05},
			 {6.675205e+51, 9.739611e-01, 1.522367e-03, 8.618425e-04, 1.231793e-04, 7.276722e-05},
			 {4.503815e+52, 9.186506e-01, 1.054325e-02, 1.742493e-02, 1.826442e-04, 1.078499e-04},
			 {2.883229e+53, 8.072469e-01, 2.280500e-02, 8.348202e-02, 4.014262e-04, 1.684695e-04},
			 {1.746463e+54, 6.558310e-01, 2.931956e-02, 1.960125e-01, 6.487986e-03, 2.781257e-04},
			 {1.016622e+55, 4.899245e-01, 3.021142e-02, 2.987737e-01, 5.118229e-02, 4.777958e-04},
			 {5.799910e+55, 3.298557e-01, 2.808534e-02, 3.445588e-01, 1.643551e-01, 8.516466e-04},
			 {3.283307e+56, 1.915271e-01, 2.361579e-02, 3.387841e-01, 3.232404e-01, 3.448759e-03},
			 {1.855450e+57, 9.198589e-02, 1.572652e-02, 2.954364e-01, 4.654125e-01, 3.036475e-02},
			 {1.050058e+58, 3.734607e-02, 7.870513e-03, 2.185504e-01, 5.319855e-01, 1.209364e-01},
			 {1.576133e+51, 9.999111e-01, 1.828011e-16, 3.751418e-20, 5.216863e-05, 3.081822e-05},
			 {1.257048e+52, 9.996359e-01, 2.584732e-09, 2.590961e-11, 6.541093e-05, 3.864101e-05},
			 {9.418773e+52, 9.943227e-01, 2.228705e-05, 2.769577e-06, 8.729875e-05, 5.157107e-05},
			 {6.675205e+53, 9.739611e-01, 1.522367e-03, 8.618425e-04, 1.231793e-04, 7.276722e-05},
			 {4.503815e+54, 9.186506e-01, 1.054325e-02, 1.742493e-02, 1.826442e-04, 1.078499e-04},
			 {2.883229e+55, 8.072469e-01, 2.280500e-02, 8.348202e-02, 4.014262e-04, 1.684695e-04},
			 {1.746463e+56, 6.558310e-01, 2.931956e-02, 1.960125e-01, 6.487986e-03, 2.781257e-04},
			 {1.016622e+57, 4.899245e-01, 3.021142e-02, 2.987737e-01, 5.118229e-02, 4.777958e-04},
			 {5.799910e+57, 3.298557e-01, 2.808534e-02, 3.445588e-01, 1.643551e-01, 8.516466e-04},
			 {3.283307e+58, 1.915271e-01, 2.361579e-02, 3.387841e-01, 3.232404e-01, 3.448759e-03},
			 {2.059903e+51, 9.999323e-01, 4.256455e-53, 3.174986e-65, 3.991679e-05, 2.358054e-05},
			 {1.853253e+52, 9.999247e-01, 1.737968e-29, 3.372028e-36, 4.436777e-05, 2.620992e-05},
			 {1.576133e+53, 9.999111e-01, 1.828011e-16, 3.751418e-20, 5.216863e-05, 3.081822e-05},
			 {1.257048e+54, 9.996359e-01, 2.584732e-09, 2.590961e-11, 6.541093e-05, 3.864101e-05},
			 {9.418773e+54, 9.943227e-01, 2.228705e-05, 2.769577e-06, 8.729875e-05, 5.157107e-05},
			 {6.675205e+55, 9.739611e-01, 1.522367e-03, 8.618425e-04, 1.231793e-04, 7.276722e-05},
			 {4.503815e+56, 9.186506e-01, 1.054325e-02, 1.742493e-02, 1.826442e-04, 1.078499e-04},
			 {2.883229e+57, 8.072469e-01, 2.280500e-02, 8.348202e-02, 4.014262e-04, 1.684695e-04},
			 {1.746463e+58, 6.558310e-01, 2.931956e-02, 1.960125e-01, 6.487986e-03, 2.781257e-04},
			 {1.016622e+59, 4.899245e-01, 3.021142e-02, 2.987737e-01, 5.118229e-02, 4.777958e-04}};
		
static int CalculateArrayIndex(float Mass, float AccRate);	
extern "C" void FORTRAN_NAME(stellar)(float *smdot, float *dt_enzo, float *slum, float *rad);
	
int DetermineSEDParameters(ActiveParticleType_SmartStar *SS, FLOAT Time, FLOAT dx)
{

  if((SmartStarBHRadiativeFeedback == FALSE) && (SmartStarStellarRadiativeFeedback == FALSE)) {
      for(int bin = 0; bin < NUMRADIATIONBINS; bin++) {
	SS->RadiationEnergyBins[bin] = 0.0;
	SS->RadiationSED[bin] = 0.0;
      }
      return SUCCESS;
    }
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
    VelocityUnits;
  double MassUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, Time);
  MassUnits = DensityUnits * POW(LengthUnits,3);
  double MassConversion = (float) (dx*dx*dx * double(MassUnits));  //convert to g from a density
  double ParticleMass =  SS->ReturnMass()*MassConversion/SolarMass; //In solar masses

#if STELLAR
  static float ptime=0.0;
  if(SS->ParticleClass == POPIII || SS->ParticleClass == SMS) {
    SS->RadiationLifetime = 1e6*yr_s/TimeUnits; //Code Time
    SS->LuminosityPerSolarMass = 6.696798e49/40.0; //In physical units

    float smdot, dt_enzo, slum, rad;
    smdot=SS->AccretionRate[SS->TimeIndex]*MassUnits/TimeUnits;
    dt_enzo=Time-ptime;

  

    if(ptime>0.0){
      FORTRAN_NAME(stellar)(&smdot, &dt_enzo, &slum, &rad);
      SS->LuminosityPerSolarMass = slum/ParticleMass;
    }
    else{
      SS->LuminosityPerSolarMass = 1e10;
    }
  }
  ptime = Time;
  return SUCCESS;
#endif
  /* 
   * The PopIII values are taken from Schaerer et al. 2002 Table 4.
   * Luminosity is NOT in ergs/s but in photons/s 
   */
  if(SS->ParticleClass == POPIII) {
    SS->RadiationLifetime = SmartStarSMSLifetime*yr_s/TimeUnits; //Code Time
    SS->LuminosityPerSolarMass = 6.696798e49/40.0; //In physical units
    for(int bin = 0; bin < NUMRADIATIONBINS; bin++) {
      SS->RadiationEnergyBins[bin] = PopIIIEnergyBins[bin];
      SS->RadiationSED[bin] = PopIIISEDFracBins[bin];
    }
  }
  /*
   * For the SMS spectrum I assume a blackbody with an effective 
   * temperature of 6000K. 
   * Luminosity is NOT in ergs/s but in photons/s 
   */
  else if (SS->ParticleClass == SMS) {
    SS->RadiationLifetime = SmartStarSMSLifetime*yr_s/TimeUnits; //Code Times
    SS->LuminosityPerSolarMass = 1.4e51/500.0;
    /* Ideally we call SLUG here. Hardcoded for now */
    for(int bin = 0; bin < NUMRADIATIONBINS; bin++) {
      SS->RadiationEnergyBins[bin] = SMSEnergyBins[bin];
      SS->RadiationSED[bin] = SMSSEDFracBins[bin];
    }
  }
  
  /*
   * The BH luminosity depends on the Black Hole mass and the accretion rate
   * The values are calculated assuming a multi-colour disc for the 
   * for accretion disk and a corona fit by a power law.
   * The values are hardcoded into the table above and valid for black
   * hole masses from 1e0 to 1e9 and accretion rates from 1e-7 to 1e2
   */
  else if(SS->ParticleClass == BH && SmartStarBHRadiativeFeedback == TRUE) {
    SS->RadiationLifetime = 1e14*yr_s/TimeUnits; //code time
    double accrate = (SS->AccretionRate[SS->TimeIndex]*(MassUnits/SolarMass)/TimeUnits)*yr_s; //Msolar/yr
    double BHMass = ParticleMass;
    float epsilon = SS->eta_disk;
    double eddrate = 4*M_PI*GravConst*BHMass*SolarMass*mh/(epsilon*clight*sigma_thompson); // g/s
    eddrate = eddrate*yr_s/SolarMass; //in Msolar/yr
    accrate = max(accrate, 1e-6);
    BHMass = min(BHMass, 1.0);
    int arrayindex = CalculateArrayIndex(BHMass, accrate);
    SS->LuminosityPerSolarMass = BHArray[arrayindex][0]/BHMass; //erg/s/msun
    
    if(SmartStarSuperEddingtonAdjustment == TRUE) {
      if(accrate > eddrate) {
	//printf("%s: Super Eddington Accretion rate detected (%f Medd) - need to modify feedback!!!!!!!!!\n", __FUNCTION__, accrate/eddrate);
	float mue = 1.22, a = 0.5;
	double Ledd = 4*M_PI*GravConst*BHMass*SolarMass*mh*mue*clight/sigma_thompson; //cgs
	double medddot = 16.0*Ledd/(clight*clight); //cgs
	double accrate_cgs = accrate*SolarMass/yr_s; 
	/* Apply Madau fit to calculate Luminosity */
	double LSuperEdd = Ledd*MadauFit(a, accrate*SolarMass/yr_s, medddot); 
	//printf("LuminosityPerSolarMass in Superdd Case is %e erg/s/msun\n", LSuperEdd/ParticleMass);
	epsilon = LSuperEdd/(accrate_cgs*clight*clight);
	//printf("epsilon updated to %f\n", epsilon);
	SS->LuminosityPerSolarMass = LSuperEdd/BHMass;
      }
    }
 
    /* Employ some ramping to stop numerical meltdown */
    float Age = (Time - SS->ReturnBirthTime())*TimeUnits/yr_s;
    //printf("%s: BH Age = %e yrs\n", __FUNCTION__, Age);
    if(Age < RAMPAGE) {
      float ramp = (Age/RAMPAGE);
      //printf("%s: ramp = %g\n", __FUNCTION__, ramp);
      SS->LuminosityPerSolarMass = SS->LuminosityPerSolarMass*ramp;
    }
    /* Use values from the BHArray to set the SED Fractions */
    for(int bin = 0; bin < NUMRADIATIONBINS; bin++) {
      SS->RadiationEnergyBins[bin] = BHEnergyBins[bin];
      SS->RadiationSED[bin] = BHArray[arrayindex][bin+1];
    }
  }
  else {
    fprintf(stderr, "%s: Particle Class = %d but no Radiative Feedback\n", __FUNCTION__, SS->ParticleClass);
    return SUCCESS;
  }

  return SUCCESS;
}


/*
 * Using the Black Hole Mass and Accretion Rate
 * calculate the index of the SED array we need.
 * Note the array is assumed to run from 1e0 -> 1e9 in Mass
 * and 1e-6 to 1e3 in msolar/yr
 * The matrix is set to be 10x10 for simplicity as it makes the 
 * indexing very easy. 
 */
static int CalculateArrayIndex(float Mass, float AccRate)
{
  float logmass = log10(Mass), logaccrate = log10(AccRate);
  int rowlength = 10;
  int accrate_offset = 5;
  int column = floor(logaccrate) + accrate_offset + 1;
  int row = floor(logmass) + 1;
  int index = rowlength*row + column;
  //printf("%s: Row: %d Column: %d Index: %d\n", __FUNCTION__, row, column, index);
  return index;
}

float MadauFit(float a, float mdot, float medddot)
{

  float A = 0.0, B = 0.0, C = 0.0;

  A = POW(0.9663 - 0.9292*a, -0.5639);
  B = POW(4.627 - 4.445*a, -0.5524);
  C = POW(827.3 - 718.1*a, -0.7060);
  float m_ratio = medddot/mdot;
  float L = A*(0.985/ (m_ratio + B) + 0.015/(m_ratio + C));

  return L;
}
