
/****************************************************************************
/
/ STAR FORMATION AND FEEDBACK ALGORITHMS FOR INDIVIDUAL STARS
/
/ written by: Andrew Emerick
/ date:       February, 2016
/ modified1:
/
/ Controls star formation for individual stars as sampled from an IMF. Stars
/ are formed stochastically following an adaptation of Goldbaum et. al. 2015
/ as in star_maker_ssn.F
/ First use case of these particles is to tie to galaxy scale chemodynamics.
/ Particle creation tags these stars with local chemical abundances and feedback
/ is tied to yield tables to deposit elements in the ISM.
*****************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "phys_constants.h"
#include "Star.h"

#include "IndividualStarProperties.h"
#include "StellarYieldsRoutines.h"


/* Following Grid_ComputeTemperatureField.C */
#ifndef MU_METAL
# define MU_METAL 16.0
#endif

/* function prototypes */

extern "C" void FORTRAN_NAME(pop3_properties)(FLOAT *mass, FLOAT* luminosity,
                                              FLOAT *lifetime);


int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);

int FindField(int f, int farray[], int n);

float SampleIMF(void);
float SampleIMF(float * data, const float & lower_mass, const float & upper_mass);
float SamplePopIII_IMF(void);

float ComputeSnIaProbability(const float &current_time, const float &formation_time, const float &lifetime, const float &TimeUnits);
unsigned_long_int mt_random(void);

void ComputeStellarWindVelocity(Star *cstar, float *v_wind);

void ComputeStellarWindMassLossRate(const float &mproj, const float &metallicity,
                                    float *dMdt);

int SetWDLifetime(float &WD_lifetime,
                    const float &current_time, const float &formation_time,
                    const float &lifetime, const float &TimeUnits);



void IndividualStarSetCoreCollapseSupernovaProperties(Star *cstar,
                                                      float &m_eject, float &E_thermal, float *metal_mass);

void IndividualStarSetPopIIISupernovaProperties(Star *cstar, float &m_eject, float &E_thermal, float *metal_mass);

void IndividualStarSetTypeIaSupernovaProperties(float &m_eject, float &E_thermal, float *metal_mass);

float ComputeOverlap(const int &i_shape, const float &radius,
                     const FLOAT &xc, const FLOAT &yc, const FLOAT &zc,
                     const FLOAT &xl, const FLOAT &yl, const FLOAT &zl,
                     const FLOAT &xr, const FLOAT &yr, const FLOAT &zr,
                     const int &nsample);


void IndividualStarSetStellarWindProperties(Star *cstar, const float &Time,
                                            const float &dtFixed, const float &TimeUnits,
                                            float &m_eject, float &E_thermal,
                                            float *metal_mass);


void ModifyStellarWindFeedback(float cell_mass, float T, float dx,
                               float MassUnits, float EnergyUnits, float &m_eject,
                               float &E_thermal, float * metal_mass,
                               float *grid_abundances);

int IndividualStarGetSETablePosition (int &i, int &j, const float &M, const float &metallicity);
int IndividualStarGetRadTablePosition(int &i, int &j, int &k,
                                      const float &Teff, const float &g, const float &metallicity);
int StellarYieldsGetYieldTablePosition(int &i, int &j, const float &M, const float &metallicity);


float IndividualStarSurfaceGravity(const float &mp, const float &R);

void IndividualStarInterpolateProperties(float &Teff, float &R,
                                         const int &i, const int &j,
                                         const float &M, const float &metallicity);



int search_lower_bound(float *arr, float value, int low, int high, int total);

int grid::GalaxySimulationInitialStars(int *nmax, int *np){

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  this->AllocateNewParticles(*nmax);

  return this->GalaxySimulationInitialStars(nmax, np, this->ParticleMass, this->ParticleType,
                                            this->ParticlePosition, this->ParticleVelocity,
                                            this->ParticleAttribute);

}

int grid::GalaxySimulationInitialStars(int *nmax, int *np, float *ParticleMass,
                                       int *ParticleType, FLOAT *ParticlePosition[],
                                       float *ParticleVelocity[], float *ParticleAttribute[]){

  /* only do this on the root processor */
  if (!(GalaxySimulationInitialStellarDist))
    return SUCCESS;

//  if !(MyProcessorNumber == ROOT_PROCESSOR)
//   return SUCCESS;

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

//  if (this->Time > 0) return SUCCESS;

  if (this->NumberOfSubgrids > 1)  // only do on highest refined grid
    return SUCCESS;

  /* get units */
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits, MassUnits;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, this->Time) == FAIL){
      ENZO_FAIL("Error in GetUnits");
  }
  MassUnits   = DensityUnits*LengthUnits*LengthUnits*LengthUnits; // mass units

  /* obtain baryon field indexes */
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, CRNum, B1Num, B2Num, B3Num;
  this->DebugCheck("StarParticleHandler");
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }

  const float msolar = 1.989E33;

  int nstar = 0;
  const int maxstar = 200000;
  FLOAT xpos[maxstar], ypos[maxstar], zpos[maxstar];
  float mass[maxstar], z[maxstar], lifetime[maxstar];

  FILE *fptr = fopen("particle_IC.in", "r");

  if (fptr == NULL){
    ENZO_FAIL("Error opening star initial positions\n");
  }

  char line[MAX_LINE_LENGTH];
  int err;
  int i = 0;

  if (IndividualStarICLifetimeMode == 2){

    while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL){
      if(line[0] != '#'){
        err = sscanf(line, "%"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM,
                           &mass[i], &z[i], &lifetime[i], &xpos[i], &ypos[i], &zpos[i]);
        i++;
      }
    }


  } else {

    while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL){
      if(line[0] != '#'){
        err = sscanf(line, "%"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM,
                           &mass[i], &z[i], &xpos[i], &ypos[i], &zpos[i]);
        i++;
      }
    }

  }

  fclose(fptr);

  nstar = i;

  int  nx = this->GridDimension[0], ny = this->GridDimension[1], nz = this->GridDimension[2];
  int  ibuff = NumberOfGhostZones;

  int count = 0;
  FLOAT cell_volume = this->CellWidth[0][0]*this->CellWidth[0][0]*this->CellWidth[0][0];
  for( i = 0; i < nstar; i++){

    xpos[i] = xpos[i]*pc_cm/LengthUnits + DiskGravityPosition[0];
    ypos[i] = ypos[i]*pc_cm/LengthUnits + DiskGravityPosition[1];
    zpos[i] = zpos[i]*pc_cm/LengthUnits + DiskGravityPosition[2];

    // make sure particle position is on this grid / processor
    if( !( (xpos[i] > this->CellLeftEdge[0][ibuff]) && (xpos[i] < this->CellLeftEdge[0][nx - ibuff] )) ||
        !( (ypos[i] > this->CellLeftEdge[1][ibuff]) && (ypos[i] < this->CellLeftEdge[1][ny - ibuff] )) ||
        !( (zpos[i] > this->CellLeftEdge[2][ibuff]) && (zpos[i] < this->CellLeftEdge[2][nz - ibuff] )) ) {
      continue;
    }


    ParticleMass[count] = mass[i] * msolar / MassUnits / (cell_volume);
    ParticleAttribute[3][count] = mass[i]; // birth mass in solar units always
    ParticleType[count] = -PARTICLE_TYPE_INDIVIDUAL_STAR;
    ParticleAttribute[0][count] = this->Time + 2.0*this->dtFixed;


    if (IndividualStarICLifetimeMode == 0){
        IndividualStarInterpolateLifetime(ParticleAttribute[1][count],
                                          mass[i], z[i], 1);
        ParticleAttribute[1][count] /= TimeUnits;

    } else if (IndividualStarICLifetimeMode == 1){
        ParticleAttribute[1][count] = 1.5 * this->dtFixed; // end life basically now

    } else if (IndividualStarICLifetimeMode == 2){ // read from file
        float Myr = 3.1556E13;
        ParticleAttribute[1][count] = lifetime[count] * Myr / TimeUnits;
    }

    ParticleAttribute[2][count]  = z[i];

    ParticlePosition[0][count] = xpos[i];
    ParticlePosition[1][count] = ypos[i];
    ParticlePosition[2][count] = zpos[i];

/*
    ParticleVelocity[0][i] = (vx[i]*1.0E5 / VelocityUnits);
    ParticleVelocity[1][i] = (vy[i]*1.0E5 / VelocityUnits);
    ParticleVelocity[2][i] = (vz[i]*1.0E5 / VelocityUnits);
*/

    int ip, jp, kp, index;
    ip = int ( (ParticlePosition[0][count] - (this->CellLeftEdge[0][0])) / (this->CellWidth[0][0]));
    jp = int ( (ParticlePosition[1][count] - (this->CellLeftEdge[1][0])) / (this->CellWidth[0][0]));
    kp = int ( (ParticlePosition[2][count] - (this->CellLeftEdge[2][0])) / (this->CellWidth[0][0]));

    index  = ip + (jp + kp * (ny)) * (nx); // flat array index
    ParticleVelocity[0][count] = BaryonField[Vel1Num][index];
    ParticleVelocity[1][count] = BaryonField[Vel2Num][index];
    ParticleVelocity[2][count] = BaryonField[Vel3Num][index];



    /* now assign metal abundnace fractions as all tiny numbers if followed */
    if (!IndividualStarOutputChemicalTags){
      if (((TestProblemData.MultiMetals == 2) || (MultiMetals == 2))){
        for (int is = 0; is < StellarYieldsNumberOfSpecies; is++){
          ParticleAttribute[4 + is][count] = tiny_number;
        }
      }
    }


    if (IndividualStarSaveTablePositions){
      /* now go trough and assign the interpolation table positions so we don't have to again */
      int tstart = ParticleAttributeTableStartIndex;

      // stellar evolution table (attr 3 = birth mass, attr 2 = metallicity)
      int t_i = -1, t_j = -1, t_k = -1;
      IndividualStarGetSETablePosition(t_i, t_j,
                                       ParticleAttribute[3][count], ParticleAttribute[2][count]);
      ParticleAttribute[tstart    ][count] = t_i;
      ParticleAttribute[tstart + 1][count] = t_j;
      // radiation properties table (only do if particle can radiate - saves time)
      if( ParticleAttribute[3][count] >= IndividualStarRadiationMinimumMass){
        float Teff, R;
        IndividualStarInterpolateProperties(Teff, R, (int)ParticleAttribute[tstart][count],
                                            (int)ParticleAttribute[tstart+1][count],
                                            ParticleAttribute[3][count], ParticleAttribute[2][count]);
        float g = IndividualStarSurfaceGravity(ParticleAttribute[3][count], R);
        t_i = -1; t_j = -1; t_k = -1;
        IndividualStarGetRadTablePosition(t_i, t_j, t_k,
                                          Teff, g, ParticleAttribute[2][count]);
        ParticleAttribute[tstart + 2][count] = t_i;
        ParticleAttribute[tstart + 3][count] = t_j;
        ParticleAttribute[tstart + 4][count] = t_k;
      } // end radiation check

      t_i = -1 ; t_j = -1;
      StellarYieldsGetYieldTablePosition(t_i, t_j,
                                         ParticleAttribute[3][count], ParticleAttribute[2][count]);
      ParticleAttribute[tstart + 5][count] = t_i;
      ParticleAttribute[tstart + 6][count] = t_j;

    } // end assign table positions

    /* Keeping this as a particle attribute */
    ParticleAttribute[NumberOfParticleAttributes-2][count] = 0.0; // wind mass ejected
    ParticleAttribute[NumberOfParticleAttributes-1][count] = 0.0; // sn mass ejected

    count++;
  }

  *np = count;

  printf("P(%"ISYM") formed %"ISYM" stars of %"ISYM"\n", MyProcessorNumber, count, nstar);

  return SUCCESS;
}


int grid::chemical_evolution_test_star_deposit(int *nmax, int *np, float *ParticleMass,
                                               int *ParticleType, FLOAT *ParticlePosition[],
                                               float *ParticleVelocity[], float *ParticleAttribute[]){


  if (ChemicalEvolutionTestNumberOfStars == 0) return SUCCESS;

  const double msolar  = 1.989e33;
  const double myr     = 3.1536e13;

  /* for convenience, rename some grid properties - will likely get optimized out */
  int  nx = this->GridDimension[0], ny = this->GridDimension[1], nz = this->GridDimension[2];
  int  ibuff = NumberOfGhostZones;

  FLOAT xstart = CellLeftEdge[0][0], ystart = CellLeftEdge[1][0], zstart = CellLeftEdge[2][0];
  float dx = this->CellWidth[0][0];

  /* identify species fields if they exist for proper computation of Mu */
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
  if ( MultiSpecies ){
    IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
                          HMNum, H2INum, H2IINum, DINum, DIINum, HDINum);
  }

  /* get metallicity tracer field number */
  int MetalNum;
  MetalNum   = FindField(Metallicity, this->FieldType, this->NumberOfBaryonFields);

  /* obtain baryon field indexes */
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, CRNum, B1Num, B2Num, B3Num;
  this->DebugCheck("StarParticleHandler");
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }

  /* get units */
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits, MassUnits;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, this->Time) == FAIL){
      ENZO_FAIL("Error in GetUnits");
  }
  MassUnits   = DensityUnits*LengthUnits*LengthUnits*LengthUnits; // mass units

//  if(this->Grid_ChemicalEvolutionTestStarFormed){
//    return SUCCESS;

 // } else
 if (ChemicalEvolutionTestNumberOfStars > 1){
    /* read in prperties from file  */

    int nstar = ChemicalEvolutionTestNumberOfStars;

    FLOAT xpos[nstar], ypos[nstar], zpos[nstar];
    float xvel[nstar], yvel[nstar], zvel[nstar];
    float mass[nstar], z[nstar];
    int pt[nstar];

    FILE *fptr = fopen("ChemicalEvolutionTest.inits", "r");
    if (fptr == NULL){
      ENZO_FAIL("Error opening star initial positions - check that you want > 1 stars and 'ChemicalEvolutionTest.inits' exists\n");
    }

    char line[MAX_LINE_LENGTH];
    int err;
    int i = 0;
    while( fgets(line, MAX_LINE_LENGTH, fptr) !=NULL){
      if(line[0] != '#'){
        err = sscanf(line, "%"FSYM " %"FSYM " %"FSYM " %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"ISYM,
                           &xpos[i], &ypos[i], &zpos[i], &xvel[i], &yvel[i], &zvel[i], &mass[i], &z[i], &pt[i]);
        i++;
      }
    }
    fclose(fptr);

    int count = 0; // total number of stars formed on this processor
    printf("nstar = %"ISYM"\n");
    for (i = 0; i < nstar; i++){

      // make sure particle position is on this grid / processor
      if( !( (xpos[i] > this->CellLeftEdge[0][ibuff]) && (xpos[i] < this->CellLeftEdge[0][nx - ibuff] )) ||
          !( (ypos[i] > this->CellLeftEdge[1][ibuff]) && (ypos[i] < this->CellLeftEdge[1][ny - ibuff] )) ||
          !( (zpos[i] > this->CellLeftEdge[2][ibuff]) && (zpos[i] < this->CellLeftEdge[2][nz - ibuff] )) ) {
        continue;
      }

      // deposit the star by hand
      ParticleMass[count] = mass[i] * msolar / MassUnits / (dx*dx*dx);
      ParticleType[count] = -pt[i];
      ParticleAttribute[0][count] = this->Time;
      ParticleNumber[count] = i; // unique ID

      // last arg tells function to return total stellar lifetime
      if(IndividualStarInterpolateLifetime(ParticleAttribute[1][i], mass[i], z[i], 1) == FAIL){
          ENZO_FAIL("Failure in stellar lifetime interpolation");
      }

      ParticleAttribute[1][count] /= TimeUnits; // convert from s to code units
      ParticleAttribute[3][count] = mass[i]; // leave in solar
      ParticleAttribute[2][count] = z[i];

      ParticlePosition[0][count] = xpos[i];
      ParticlePosition[1][count] = ypos[i];
      ParticlePosition[2][count] = zpos[i];

      ParticleVelocity[0][count] = xvel[i]*1.0E5 / VelocityUnits;
      ParticleVelocity[1][count] = yvel[i]*1.0E5 / VelocityUnits;
      ParticleVelocity[2][count] = zvel[i]*1.0E5 / VelocityUnits;


      // find grid cell and assign chemical tags
      int ip, jp, kp, n;
      ip = int ( (ParticlePosition[0][count] - (xstart)) / (dx));
      jp = int ( (ParticlePosition[1][count] - (ystart)) / (dx));
      kp = int ( (ParticlePosition[2][count] - (zstart)) / (dx));

      n  = ip + (jp + kp * (ny)) * (nx);

      /* Metal fields are all in fractions, as set in Grid_StarParticleHandler */
      if(((TestProblemData.MultiMetals == 2) || (MultiMetals == 2))){
        for( int ii = 0; ii < StellarYieldsNumberOfSpecies; ii++){
          if(StellarYieldsAtomicNumbers[ii] > 2){
            int field_num;
            this->IdentifyChemicalTracerSpeciesFieldsByNumber(field_num, StellarYieldsAtomicNumbers[ii]);

            ParticleAttribute[4 + ii][count] = BaryonField[field_num][n];
          } else if (StellarYieldsAtomicNumbers[ii] == 1){
            /* Take H and He fractions as TOTAL amount of H and He species in the cell */
            ParticleAttribute[4 + ii][count] = BaryonField[HINum][n] + BaryonField[HIINum][n];
            if (MultiSpecies > 1){
              ParticleAttribute[4 + ii][count] += BaryonField[HMNum][n] +
                                     BaryonField[H2INum][n] + BaryonField[H2IINum][n];
            }

          } else if (StellarYieldsAtomicNumbers[ii] == 2){

            ParticleAttribute[4 + ii][count] = BaryonField[HeINum][n]  +
                                               BaryonField[HeIINum][n] + BaryonField[HeIIINum][n];

          }
        } // end loop over species
      } // end species tagging

      /* now go trough and assign the interpolation table positions so we don't have to again */
      int tstart = ParticleAttributeTableStartIndex;

      // stellar evolution table (attr 3 = birth mass, attr 2 = metallicity)
      int t_i = -1, t_j = -1, t_k = -1;
      IndividualStarGetSETablePosition(t_i, t_j,
                                       ParticleAttribute[3][count], ParticleAttribute[2][count]);
      ParticleAttribute[tstart    ][count] = t_i;
      ParticleAttribute[tstart + 1][count] = t_j;
      // radiation properties table (only do if particle can radiate - saves time)
      if( ParticleAttribute[3][count] >= IndividualStarRadiationMinimumMass){
        float Teff, R;
        IndividualStarInterpolateProperties(Teff, R, (int)ParticleAttribute[tstart][count],
                                            (int)ParticleAttribute[tstart+1][count],
                                            ParticleAttribute[3][count], ParticleAttribute[2][count]);
        float g = IndividualStarSurfaceGravity(ParticleAttribute[3][count], R);

        t_i = -1; t_j = -1; t_k = -1;
        IndividualStarGetRadTablePosition(t_i, t_j, t_k,
                                         Teff, g, ParticleAttribute[2][count]);
        ParticleAttribute[tstart + 2][count] = t_i;
        ParticleAttribute[tstart + 3][count] = t_j;
        ParticleAttribute[tstart + 4][count] = t_k;
      }
       // yields table position
      t_i = -1 ; t_j = -1;
      StellarYieldsGetYieldTablePosition(t_i, t_j,
                                         ParticleAttribute[3][count], ParticleAttribute[2][count]);
      ParticleAttribute[tstart + 5][count] = t_i;
      ParticleAttribute[tstart + 6][count] = t_j;

      ParticleAttribute[NumberOfParticleAttributes-2][count] = 0.0; // wind mass ejected
      ParticleAttribute[NumberOfParticleAttributes-1][count] = 0.0; // sn mass ejected


      count++;
    } // end loop over particles
     *np = count;
     this->Grid_ChemicalEvolutionTestStarFormed = TRUE;
     return SUCCESS;
  } else {
    FLOAT xx, yy, zz;
    xx = ChemicalEvolutionTestStarPosition[0];
    yy = ChemicalEvolutionTestStarPosition[1];
    zz = ChemicalEvolutionTestStarPosition[2];

    // make sure particle position is on this grid / processor
    if( !( (xx > this->CellLeftEdge[0][ibuff ]) && (xx < this->CellLeftEdge[0][nx - ibuff ] )) ||
        !( (yy > this->CellLeftEdge[1][ibuff ]) && (yy < this->CellLeftEdge[1][ny - ibuff ] )) ||
        !( (zz > this->CellLeftEdge[2][ibuff ]) && (zz < this->CellLeftEdge[2][nz - ibuff ] )) ) {
      this->Grid_ChemicalEvolutionTestStarFormed = TRUE; // setting this here to avoid doing MPI communication
                                              // on whatever processor the star actually gets placed
      printf("P(%"ISYM") individual_star_maker: Particle not on this grid. Leaving\n", MyProcessorNumber);
      return SUCCESS;
    }
     // deposit the star by hand
    ParticleMass[0] = ChemicalEvolutionTestStarMass * msolar / MassUnits / (dx*dx*dx);
    ParticleType[0] = - PARTICLE_TYPE_INDIVIDUAL_STAR;
    ParticleAttribute[0][0] = this->Time;
     // allow user to set lifetime artificially
    if(ChemicalEvolutionTestStarLifetime > 0){
      ParticleAttribute[1][0] = ChemicalEvolutionTestStarLifetime * myr / (TimeUnits);
    } else{
      // last arg tells function to return total stellar lifetime
      if(IndividualStarInterpolateLifetime(ParticleAttribute[1][0], ChemicalEvolutionTestStarMass,
                                                                    ChemicalEvolutionTestStarMetallicity, 1) == FAIL){
          ENZO_FAIL("Failure in stellar lifetime interpolation");
      }
       ParticleAttribute[1][0] /= TimeUnits; // convert from s to code units
    }
    ParticleAttribute[3][0] = ChemicalEvolutionTestStarMass; // in solar!!!
    ParticleAttribute[2][0] = ChemicalEvolutionTestStarMetallicity;
    ParticlePosition[0][0] = ChemicalEvolutionTestStarPosition[0];
    ParticlePosition[1][0] = ChemicalEvolutionTestStarPosition[1];
    ParticlePosition[2][0] = ChemicalEvolutionTestStarPosition[2];
    ParticleVelocity[0][0] = ChemicalEvolutionTestStarVelocity[0]*1.0E5 / VelocityUnits;
    ParticleVelocity[1][0] = ChemicalEvolutionTestStarVelocity[1]*1.0E5 / VelocityUnits;
    ParticleVelocity[2][0] = ChemicalEvolutionTestStarVelocity[2]*1.0E5 / VelocityUnits;

     // find grid cell and assign chemical tags
    int ip, jp, kp, n;
    ip = int ( (ParticlePosition[0][0] - (xstart)) / (dx));
    jp = int ( (ParticlePosition[1][0] - (ystart)) / (dx));
    kp = int ( (ParticlePosition[2][0] - (zstart)) / (dx));
    n  = ip + (jp + kp * (ny)) * (nx);

    if (! IndividualStarOutputChemicalTags){
      /* Metal fields are all in fractions, as set in Grid_StarParticleHandler */
      if(((TestProblemData.MultiMetals == 2) || (MultiMetals == 2))){
        for( int ii = 0; ii < StellarYieldsNumberOfSpecies; ii++){
          if(StellarYieldsAtomicNumbers[ii] > 2){
            int field_num;
             this->IdentifyChemicalTracerSpeciesFieldsByNumber(field_num, StellarYieldsAtomicNumbers[ii]);
             ParticleAttribute[4 + ii][0] = BaryonField[field_num][n];
           } else if (StellarYieldsAtomicNumbers[ii] == 1){
            /* Take H and He fractions as TOTAL amount of H and He species in the cell */
           ParticleAttribute[4 + ii][0] = BaryonField[HINum][n] + BaryonField[HIINum][n];
             if (MultiSpecies > 1){
              ParticleAttribute[4 + ii][0] += BaryonField[HMNum][n] +
                                           BaryonField[H2INum][n] + BaryonField[H2IINum][n];
            }
           } else if (StellarYieldsAtomicNumbers[ii] == 2){
             ParticleAttribute[4 + ii][0] = BaryonField[HeINum][n]  +
                                           BaryonField[HeIINum][n] + BaryonField[HeIIINum][n];
           }
        }
      }
    } // check if we are saving chemical tags

    if(IndividualStarSaveTablePositions){

      int tstart = ParticleAttributeTableStartIndex;

      // stellar evolution table (attr 3 = birth mass, attr 2 = metallicity)
      int t_i = -1, t_j = -1, t_k = -1;
      IndividualStarGetSETablePosition(t_i, t_j,
                                       ParticleAttribute[3][0], ParticleAttribute[2][0]);
      ParticleAttribute[tstart    ][0] = t_i;
      ParticleAttribute[tstart + 1][0] = t_j;
      // radiation properties table (only do if particle can radiate - saves time)
      if( ParticleAttribute[3][0] >= IndividualStarRadiationMinimumMass){
         float Teff, R;
         IndividualStarInterpolateProperties(Teff, R, (int)ParticleAttribute[tstart][0],
                                            (int)ParticleAttribute[tstart+1][0],
                                            ParticleAttribute[3][0], ParticleAttribute[2][0]);
         float g = IndividualStarSurfaceGravity(ParticleAttribute[3][0], R);
         t_i = -1; t_j = -1; t_k = -1;
         IndividualStarGetRadTablePosition(t_i, t_j, t_k,
                                           Teff, g, ParticleAttribute[2][0]);
                                           ParticleAttribute[tstart + 2][0] = t_i;
         ParticleAttribute[tstart + 3][0] = t_j;
         ParticleAttribute[tstart + 4][0] = t_k;
      } else {
         ParticleAttribute[tstart + 2][0] = -1;
         ParticleAttribute[tstart + 3][0] = -1;
         ParticleAttribute[tstart + 4][0] = -1;
      }
       // yields table position
      t_i = -1 ; t_j = -1;
      StellarYieldsGetYieldTablePosition(t_i, t_j,
                                         ParticleAttribute[3][0], ParticleAttribute[2][0]);
      ParticleAttribute[tstart + 5][0] = t_i;
      ParticleAttribute[tstart + 6][0] = t_j;

    } // end table position save check

    /* Keeping this as a particle attribute */
    ParticleAttribute[NumberOfParticleAttributes-2][0] = 0.0; // wind mass ejected
    ParticleAttribute[NumberOfParticleAttributes-1][0] = 0.0; // sn mass ejected



    *np = 1;
    this->Grid_ChemicalEvolutionTestStarFormed = TRUE;
    printf("individual_star_maker: Formed star ChemicalEvolutionTest. M =  %"FSYM" and Z = %"FSYM". tau = %"ESYM"\n", ParticleMass[0]*(dx*dx*dx)*MassUnits/msolar, ParticleAttribute[2][0], ParticleAttribute[1][0]);
  }

  return SUCCESS;
}

int grid::individual_star_maker(float *dm, float *temp, int *nmax, float *mu, int *np,
                                float *ParticleMass,
                                int *ParticleType, FLOAT *ParticlePosition[],
                                float *ParticleVelocity[], float *ParticleAttribute[],
                                float *StellarAbundances[]){
/*-----------------------------------------------------------------------------
  INPUTS:
    dm          - dark matter density field (computed in Grid_StarParticleHandler)
    temp        - temperature field (computed in Grid_StarParticleHandler)
    nmax        - Maximum allowed number of stars that can form on a single grid
    mu          - global Mean Molecular weight of gas
    ctype       - number for desired particle type assignment

  OUTPUTS: SUCCESS or FAIL
    Creates star particle and updates all particle arrays
    modifies baryon fields during star formation
    np - number of particles created
    ParticleMass - particle masses on grid
    ParticleType - particle types on grid
    ParticlePosition - particle positions on grid
    ParticleVelocity - particle velocities
    ParticleAttribute - particle attributes
-----------------------------------------------------------------------------*/

  if (ProblemType == 31 && GalaxySimulationInitialStellarDist && this->Time <= 0.0){

    if (this->GalaxySimulationInitialStars(nmax, np, ParticleMass, ParticleType,
                                           ParticlePosition, ParticleVelocity,
                                           ParticleAttribute) == FAIL){
      return FAIL;
    }

    return SUCCESS;
  }

  const double msolar  = 1.989e33;
  const double sndspdC = 1.3095e8;
  const double myr     = 3.1536e13;

  int i, j, k, index, ii=0, istar=0, index_presf=0;
  int xo, yo, zo, rsign=1;
  float bmass, div, star_mass=0.0, sum_mass=0.0, metal_mass=0.0, H2mass=0.0;
  float pstar, mass_to_stars, mass_available, tdyn;
  float dtot, isosndsp2, jeansmass, star_fraction;
  float umean, vmean, wmean, px, py, pz, px_excess, py_excess, pz_excess;
  float rnum;
  const double m_h = 1.673e-24;
  float inv_metal_mol = 1.0 /MU_METAL;

  const int max_random = (1<<16);

  int form_method = -1; // tracker for debugging purposes

  /* for convenience, rename some grid properties - will likely get optimized out */
  int  nx = this->GridDimension[0], ny = this->GridDimension[1], nz = this->GridDimension[2];
  int  ibuff = NumberOfGhostZones;

  FLOAT xstart = CellLeftEdge[0][0], ystart = CellLeftEdge[1][0], zstart = CellLeftEdge[2][0];
  float   dx   = CellWidth[0][0];

  if (! this->isLocal()) return SUCCESS;

  if ( this->dtFixed == 0.0){
    printf("DT EQUAL TO ZERO\N");
    return FAIL;
  }

  /* obtain baryon field indexes */
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, CRNum, B1Num, B2Num, B3Num;
  this->DebugCheck("StarParticleHandler");
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }

  /* identify species fields if they exist for proper computation of Mu */
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
  if ( MultiSpecies ){
    IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
                          HMNum, H2INum, H2IINum, DINum, DIINum, HDINum);
  }

  /* get metallicity tracer field number */
  int MetalNum;
  MetalNum   = FindField(Metallicity, this->FieldType, this->NumberOfBaryonFields);

 int PopIIIMetalNum, AGBMetalNum, SNIaMetalNum, SNIIMetalNum;

  AGBMetalNum    = FindField(ExtraType0, FieldType, NumberOfBaryonFields);
  PopIIIMetalNum = FindField(ExtraType1, FieldType, NumberOfBaryonFields);
  SNIaMetalNum   = FindField(MetalSNIaDensity, FieldType, NumberOfBaryonFields);
  SNIIMetalNum   = FindField(MetalSNIIDensity, FieldType, NumberOfBaryonFields);

  if ( IndividualStarTrackAGBMetalDensity && (AGBMetalNum <= 0)){
    ENZO_FAIL("Error in finding AGB metal density field in individual_star_maker");
  }

  if ( IndividualStarPopIIIFormation && (PopIIIMetalNum <= 0)){
    ENZO_FAIL("Error in finding Pop III metal density field in individual_star_maker");
  }

  if ( IndividualStarTrackSNMetalDensity && ( (SNIIMetalNum <= 0) || (SNIaMetalNum <=0))){
    ENZO_FAIL("Error in finding SNII and SNIa metal density field in individual_star_maker.");
  }

  /* get units */
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits, MassUnits;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, this->Time) == FAIL){
      ENZO_FAIL("Error in GetUnits");
  }
  MassUnits   = DensityUnits*LengthUnits*LengthUnits*LengthUnits; // mass units


    // 3D -> 1D index
    xo = 1;
    yo = nx;
    zo = (nx) * (ny);

    // over density threshold in code units
    // if multispecies is off, assumes a value for MU
    float odthreshold, secondary_odthreshold;
    if (MultiSpecies == FALSE){
        odthreshold           = StarMakerOverDensityThreshold * m_h * (*mu) / (DensityUnits); // code density
        secondary_odthreshold = IndividualStarSecondaryOverDensityThreshold * m_h * (*mu) / (DensityUnits);
    }

    // loop over all cells, check condition, form stars stochastically
    ii = 0; index_presf = 0;

    int number_of_sf_cells = 0;
    int integer_sep = ((int) (IndividualStarCreationStencilSize + 1) / 2.0 - 1); // stencil size must be odd number

    float *ke_before;
    if (HydroMethod != 2){
      ke_before = new float[ ((int) POW(IndividualStarCreationStencilSize,3)) ];

      for (i = 0; i < ((int) POW(IndividualStarCreationStencilSize,3)); i++) ke_before[i] = 0.0;

    } else{
      ke_before = NULL;
    }

    for (k = ibuff; k < nz - ibuff; k++){
      for (j = ibuff; j < ny - ibuff; j++){
        for (i = ibuff; i < nx - ibuff; i++){

          index = i + ( j + k * (ny)) * nx;

          if (ProblemType == 30){
            if ( BaryonField[NumberOfBaryonFields][index] != 0) continue; // not highest refined zone of this region
          }

          /* if distributed star formation */
          // check center cell's SF condition, if met, do SF
          /* loop and sum over all*/

           bmass      = (BaryonField[DensNum][index]*(dx*dx*dx)) * MassUnits / msolar; // in solar masses
           metal_mass = (BaryonField[MetalNum][index] * bmass); // metalNum is converted to fraction in Grid_StarParticleHandler

               // perform the following easy checks for SF before proceeding
               // 1) Is density greater than the density threshold?
               // 2) Is temperature < the minimum temperature?
               // 3) Do not allow star formation if minimum star particle mass is
               //    above some fraction of cell mass. This is very unlikely to occur
               //    in intended use case:
               //    (i.e. baryon mass likely always going to be > ~10 solar masses)

          sum_mass = 0.0; index_presf = ii;

          /* Need to compute Mu exactly and recompute threshold */
          if ( MultiSpecies ){
            float mu_cell, number_density;

            number_density =
              0.25*(BaryonField[HeINum][index] + BaryonField[HeIINum][index] +
                    BaryonField[HeIIINum][index])                        +
                   BaryonField[HINum][index] + BaryonField[HIINum][index]    +
                   0.5 * BaryonField[DeNum][index];
            /* if H2 is present */
            if (MultiSpecies > 1){
              number_density += BaryonField[HMNum][index] +
                          0.5*(BaryonField[H2INum][index] + BaryonField[H2IINum][index]);

              H2mass = (BaryonField[H2INum][index] + BaryonField[H2IINum][index]) * bmass;
            }

            /* Metal field must be present in this star formation scheme */
            number_density += BaryonField[MetalNum][index] * inv_metal_mol;


            number_density *= BaryonField[DensNum][index] / m_h ; // now actual n density

            mu_cell = BaryonField[DensNum][index] / (number_density * m_h);
            odthreshold = StarMakerOverDensityThreshold * (mu_cell) * m_h / (DensityUnits);

            secondary_odthreshold = IndividualStarSecondaryOverDensityThreshold * (mu_cell) * m_h / DensityUnits;
          }


          if (   BaryonField[DensNum][index]      > odthreshold
              && temp[index] <= IndividualStarTemperatureThreshold){
              //&& IndividualStarMassFraction*bmass > IndividualStarIMFLowerMassCutoff
              //&& 0.5*bmass > IndividualStarIMFUpperMassCutoff){

            // star formation may be possible
            // compute values and check jeans mass unstable

            // AJE: Feb 2017 - Apparently the dm field is only computed for dm particles, and not
            //                 computed for the static background used in the isolated galaxy sims.
            //                 This *shouldn't* be an issue if the density threshold is high, as the SF
            //                 regions should be dom by self-gravity and the local DM density should be
            //                 much less than the local baryon density... this should be fixed
            //                 if used in low resolution simulations

            dtot = ( BaryonField[DensNum][index] + dm[index] ) * (DensityUnits);         // total density
            tdyn = sqrt(3.0 * pi / 32.0 / GravConst / dtot) / (TimeUnits);            // in code units

            if (StarMakerUseJeansMass){
              isosndsp2 = sndspdC * temp[index] ;
              jeansmass = pi / (6.0 * sqrt(BaryonField[DensNum][index]*DensityUnits) *
                              POW(pi * isosndsp2 / GravConst ,1.5)) / msolar; // in solar masses
            } else{
              jeansmass = -1.0; // so always < bmass
            }

            float vel_div = -1.0;
            if (IndividualStarCheckVelocityDiv){
              if (HydroMethod == 2){
                vel_div = BaryonField[Vel1Num][index + xo] - BaryonField[Vel1Num][index] +
                          BaryonField[Vel2Num][index + yo] - BaryonField[Vel2Num][index] +
                          BaryonField[Vel3Num][index + zo] - BaryonField[Vel3Num][index];
              } else{
                vel_div = BaryonField[Vel1Num][index + xo] - BaryonField[Vel1Num][index - xo] +
                          BaryonField[Vel2Num][index + yo] - BaryonField[Vel2Num][index - yo] +
                          BaryonField[Vel3Num][index + zo] - BaryonField[Vel3Num][index - zo];
              }
            }

            if (jeansmass <= bmass && vel_div < 0.0){
              float lowest_cell_mass = bmass;
              bmass = 0.0; number_of_sf_cells = 0;

              int istart, iend, jstart, jend, kstart, kend;

              // compute indeces for adjacent cells integer_sep away from center in each dir
              // stop at grid edge if cell near boundary
              istart = iend = jstart = jend = kstart = kend = 0;
              if (integer_sep > 0){
                istart   = min( i - ibuff             , integer_sep);
                iend     = min( (nx - ibuff - 1 ) - i, integer_sep);
                jstart   = min( j - ibuff             , integer_sep);
                jend     = min( (ny - ibuff - 1 ) - j, integer_sep);
                kstart   = min( k - ibuff             , integer_sep);
                kend     = min( (nz - ibuff - 1 ) - k, integer_sep);
              }

              // loop through cells and add up total amount of mass available for SF
              int l = 0;
              for (int k_loc = -kstart; k_loc <= kend; k_loc++){
                for(int j_loc = -jstart; j_loc <= jend; j_loc++){
                  for (int i_loc = -istart; i_loc <= iend; i_loc++){
                    int loc_index = (i + i_loc) + ( (j + j_loc) + (k + k_loc)*(ny))*(nx);

                    if(BaryonField[DensNum][loc_index] > secondary_odthreshold &&
                       temp[loc_index] < IndividualStarTemperatureThreshold       ){
                      float current_cell_mass = BaryonField[DensNum][loc_index]*(dx*dx*dx)*MassUnits/msolar;

                      bmass += current_cell_mass; // solar masses
                      metal_mass += BaryonField[MetalNum][loc_index] * current_cell_mass;

                      if (MultiSpecies > 1){
                        H2mass += (BaryonField[H2INum][loc_index] + BaryonField[H2IINum][loc_index])*\
                                     current_cell_mass;
                      }

                      number_of_sf_cells++;
                      if ( current_cell_mass < lowest_cell_mass){
                        lowest_cell_mass = current_cell_mass;
                      }
                    }

                    /* if PPM, need to be careful about energies */
                    if (HydroMethod != 2 && FALSE){
                      ke_before[ l ] = 0.5 * BaryonField[DensNum][loc_index] *
                                 ( BaryonField[Vel1Num][loc_index] * BaryonField[Vel1Num][loc_index] +
                                   BaryonField[Vel2Num][loc_index] * BaryonField[Vel2Num][loc_index] +
                                   BaryonField[Vel3Num][loc_index] * BaryonField[Vel3Num][loc_index]);
                    }
                    l++;
                  }
                }
              }


              // only allow star formation if IMF can fully sampled (safely) in the given region
              //   aka, make sure there is enough mass to form the most massive star + some fudge
              if ( bmass * IndividualStarMassFraction < IndividualStarIMFUpperMassCutoff){
                break;
              }

              int add_unresolved_star = FALSE; // only used when IMF mass floor is < imf lower limit
              int form_popIII_stars   = 0;

              /* Here is where star formation actually occurs */
              if( bmass*IndividualStarMassFraction > IndividualStarSFGasMassThreshold ){
                // if true, we can try and form stars. compute probability that this mass will
                // form stars this timestep
                star_fraction  = min(StarMakerMassEfficiency*(this->dtFixed)/tdyn, 1.0);
                mass_to_stars  = star_fraction * bmass;

                pstar          = mass_to_stars / IndividualStarSFGasMassThreshold;

                rnum           = (float) (random() % max_random) / ( (float) max_random);

                if ( rnum < pstar){ // form stars until mass runs out - keep star if too much is made

                  form_popIII_stars = ((IndividualStarPopIIIFormation) *\
                                       ((metal_mass/bmass) <= PopIIIMetalCriticalFraction) ); // critial metallicity

                  if (form_popIII_stars){

                    if ( (H2mass/bmass) > PopIIIH2CriticalFraction){ // must check this separately
                      float mass_counter    = IndividualStarSFGasMassThreshold;
                      float unresolved_mass = 0.0;

                      while( mass_counter > 0.0){
                        float temp_mass = SamplePopIII_IMF();

                        ParticleMass[ii] = temp_mass;
                        ParticleType[ii] = -PARTICLE_TYPE_INDIVIDUAL_STAR_POPIII;
                        ii++;
                        sum_mass += temp_mass;
                        mass_counter -= temp_mass;
                      }
                    } // else, do not form any stars!!

                  } else {
                    float mass_counter = IndividualStarSFGasMassThreshold;
                    float unresolved_mass = 0.0;
                    while( mass_counter > 0.0){
                        float temp_mass = SampleIMF();

                        // if the mass is between the lower mass cutoff and the mass floor, sum the total mass of
                        // these particles and dump into a single particle at the very end.
                        if ( (temp_mass >= IndividualStarIMFLowerMassCutoff) && (temp_mass < IndividualStarIMFMassFloor)){
                          unresolved_mass  += temp_mass;
                        } else{
                          ParticleMass[ii]  = temp_mass;
                          ParticleType[ii]  = -PARTICLE_TYPE_INDIVIDUAL_STAR;
                          ii++;
                        }
                        sum_mass         += temp_mass;
                        mass_counter     -= temp_mass;
                    }

                    if (unresolved_mass > 0.0) { // we've formed tiny stars (should always happen).. account for this
                      add_unresolved_star = TRUE;
                      ParticleMass[ii] = unresolved_mass;
                      ParticleType[ii] = -PARTICLE_TYPE_INDIVIDUAL_STAR_UNRESOLVED;
                      ii++;
                    }
                  } // check if doing popIII or normal SF

                } // endif randum number draw check

              } // endif mass threshold check

              // prepare for assigning star properties by computing the local
              // gas velocity properties (this is for velocity assignment)
              // 2 = Zeus .. otherwise PPM
              // copied from pop3_maker.F
              if (HydroMethod == 2){
                umean = (
                       0.5 * (BaryonField[Vel1Num][index   ] + BaryonField[Vel1Num][index+xo])*BaryonField[DensNum][index] +
                       0.5 * (BaryonField[Vel1Num][index-xo] + BaryonField[Vel1Num][index   ])*BaryonField[DensNum][index-xo] +
                       0.5 * (BaryonField[Vel1Num][index+xo] + BaryonField[Vel1Num][index + xo + xo])*BaryonField[DensNum][index+xo] +
                       0.5 * (BaryonField[Vel1Num][index+yo] + BaryonField[Vel1Num][index + xo + yo])*BaryonField[DensNum][index+yo] +
                       0.5 * (BaryonField[Vel1Num][index-yo] + BaryonField[Vel1Num][index + xo - yo])*BaryonField[DensNum][index-yo] +
                       0.5 * (BaryonField[Vel1Num][index+zo] + BaryonField[Vel1Num][index + xo + zo])*BaryonField[DensNum][index+zo] +
                       0.5 * (BaryonField[Vel1Num][index-zo] + BaryonField[Vel1Num][index + xo - zo])*BaryonField[DensNum][index-zo]) /
                      ( BaryonField[DensNum][index] + BaryonField[DensNum][index-xo] + BaryonField[DensNum][index+xo] +
                        BaryonField[DensNum][index-yo] + BaryonField[DensNum][index+yo] +
                        BaryonField[DensNum][index-zo] + BaryonField[DensNum][index+zo] ); //
                vmean = (
                         0.5 * (BaryonField[Vel2Num][index   ] + BaryonField[Vel2Num][index+xo])*BaryonField[DensNum][index] +
                         0.5 * (BaryonField[Vel2Num][index-xo] + BaryonField[Vel2Num][index   ])*BaryonField[DensNum][index-xo] +
                         0.5 * (BaryonField[Vel2Num][index+xo] + BaryonField[Vel2Num][index + xo + xo])*BaryonField[DensNum][index+xo] +
                         0.5 * (BaryonField[Vel2Num][index+yo] + BaryonField[Vel2Num][index + xo + yo])*BaryonField[DensNum][index+yo] +
                         0.5 * (BaryonField[Vel2Num][index-yo] + BaryonField[Vel2Num][index + xo - yo])*BaryonField[DensNum][index-yo] +
                         0.5 * (BaryonField[Vel2Num][index+zo] + BaryonField[Vel2Num][index + xo + zo])*BaryonField[DensNum][index+zo] +
                         0.5 * (BaryonField[Vel2Num][index-zo] + BaryonField[Vel2Num][index + xo - zo])*BaryonField[DensNum][index-zo]) /
                        ( BaryonField[DensNum][index] + BaryonField[DensNum][index-xo] + BaryonField[DensNum][index+xo] +
                         BaryonField[DensNum][index-yo] + BaryonField[DensNum][index+yo] +
                          BaryonField[DensNum][index-zo] + BaryonField[DensNum][index+zo] ); //
                wmean = (
                         0.5 * (BaryonField[Vel3Num][index   ] + BaryonField[Vel3Num][index+xo])*BaryonField[DensNum][index] +
                         0.5 * (BaryonField[Vel3Num][index-xo] + BaryonField[Vel3Num][index   ])*BaryonField[DensNum][index-xo] +
                         0.5 * (BaryonField[Vel3Num][index+xo] + BaryonField[Vel3Num][index + xo + xo])*BaryonField[DensNum][index+xo] +
                         0.5 * (BaryonField[Vel3Num][index+yo] + BaryonField[Vel3Num][index + xo + yo])*BaryonField[DensNum][index+yo] +
                         0.5 * (BaryonField[Vel3Num][index-yo] + BaryonField[Vel3Num][index + xo - yo])*BaryonField[DensNum][index-yo] +
                         0.5 * (BaryonField[Vel3Num][index+zo] + BaryonField[Vel3Num][index + xo + zo])*BaryonField[DensNum][index+zo] +
                         0.5 * (BaryonField[Vel3Num][index-zo] + BaryonField[Vel3Num][index + xo - zo])*BaryonField[DensNum][index-zo]) /
                        ( BaryonField[DensNum][index] + BaryonField[DensNum][index-xo] + BaryonField[DensNum][index+xo] +
                          BaryonField[DensNum][index-yo] + BaryonField[DensNum][index+yo] +
                          BaryonField[DensNum][index-zo] + BaryonField[DensNum][index+zo] ); //
              }
              else{ // PPM case
                umean = (BaryonField[Vel1Num][index]*BaryonField[DensNum][index] +
                              BaryonField[Vel1Num][index-xo]*BaryonField[DensNum][index-xo] +
                              BaryonField[Vel1Num][index+xo]*BaryonField[DensNum][index+xo] +
                              BaryonField[Vel1Num][index-yo]*BaryonField[DensNum][index-yo] +
                              BaryonField[Vel1Num][index+yo]*BaryonField[DensNum][index+yo] +
                              BaryonField[Vel1Num][index+zo]*BaryonField[DensNum][index+zo] +
                              BaryonField[Vel1Num][index-zo]*BaryonField[DensNum][index-zo] ) /
                              (BaryonField[DensNum][index] + BaryonField[DensNum][index-xo] + BaryonField[DensNum][index+xo] +
                               BaryonField[DensNum][index-yo] + BaryonField[DensNum][index+yo] +
                               BaryonField[DensNum][index-zo] + BaryonField[DensNum][index+zo]);
                vmean = (BaryonField[Vel2Num][index]*BaryonField[DensNum][index] +
                              BaryonField[Vel2Num][index-xo]*BaryonField[DensNum][index-xo] +
                              BaryonField[Vel2Num][index+xo]*BaryonField[DensNum][index+xo] +
                              BaryonField[Vel2Num][index-yo]*BaryonField[DensNum][index-yo] +
                              BaryonField[Vel2Num][index+yo]*BaryonField[DensNum][index+yo] +
                              BaryonField[Vel2Num][index+zo]*BaryonField[DensNum][index+zo] +
                              BaryonField[Vel2Num][index-zo]*BaryonField[DensNum][index-zo] ) /
                              (BaryonField[DensNum][index] + BaryonField[DensNum][index-xo] + BaryonField[DensNum][index+xo] +
                               BaryonField[DensNum][index-yo] + BaryonField[DensNum][index+yo] +
                               BaryonField[DensNum][index-zo] + BaryonField[DensNum][index+zo]);

                wmean = (BaryonField[Vel3Num][index]*BaryonField[DensNum][index] +
                              BaryonField[Vel3Num][index-xo]*BaryonField[DensNum][index-xo] +
                              BaryonField[Vel3Num][index+xo]*BaryonField[DensNum][index+xo] +
                              BaryonField[Vel3Num][index-yo]*BaryonField[DensNum][index-yo] +
                              BaryonField[Vel3Num][index+yo]*BaryonField[DensNum][index+yo] +
                              BaryonField[Vel3Num][index+zo]*BaryonField[DensNum][index+zo] +
                              BaryonField[Vel3Num][index-zo]*BaryonField[DensNum][index-zo] ) /
                              (BaryonField[DensNum][index] + BaryonField[DensNum][index-xo] + BaryonField[DensNum][index+xo] +
                               BaryonField[DensNum][index-yo] + BaryonField[DensNum][index+yo] +
                               BaryonField[DensNum][index-zo] + BaryonField[DensNum][index+zo]);
              } // imethod velocity computation


              // now assign particle properties, loop over every star
              px = 0.0; py = 0.0; pz =0.0; // initialize momentum counters
              for (istar = index_presf; istar < ii; istar++){

                ParticleAttribute[0][istar]    = this->Time;                        // formation time
                ParticleAttribute[2][istar]    = metal_mass / bmass ; //BaryonField[MetalNum][index]; // metal fraction (conv from density in Grid_StarParti$

                if (ParticleType[istar] == -PARTICLE_TYPE_INDIVIDUAL_STAR){

                  if(IndividualStarInterpolateLifetime(ParticleAttribute[1][istar], ParticleMass[istar],
                                                                                    ParticleAttribute[2][istar], 1) == FAIL){
                    printf(" %"ESYM"  %"ESYM"  %"ESYM"\n",ParticleAttribute[1][istar], ParticleMass[istar], ParticleAttribute[2][istar]);
                  ENZO_FAIL("Error in stellar lifetime interpolation");
                  }

                } else if (ParticleType[istar] == -PARTICLE_TYPE_INDIVIDUAL_STAR_UNRESOLVED){

                  ParticleAttribute[1][istar] = huge_number * TimeUnits; // make sure its VERY large

                } else if (ParticleType[istar] == -PARTICLE_TYPE_INDIVIDUAL_STAR_POPIII){

                  float temp_mass, temp_lifetime, temp_luminosity;
                  temp_mass = 1.0 * ParticleMass[istar];

                  FORTRAN_NAME(pop3_properties)(&temp_mass, &temp_luminosity, &temp_lifetime);

                  ParticleAttribute[1][istar] = temp_lifetime * 3.1556E7; // in seconds

                } // end check for particle type for assigning lifetimes

                ParticleAttribute[1][istar] /= TimeUnits; // convert from s to code units


                ParticleAttribute[3][istar]    = ParticleMass[istar]; //progenitor mass in solar (main sequence mass)
                ParticleMass[istar]            = ParticleMass[istar] * msolar / MassUnits;   // mass in code (not yet dens)

               // give the star particle a position chosen at random
                // within the cell ( so they are not all at cell center )

                rnum =  (float) (random() % max_random) / (float) (max_random);
                ParticlePosition[0][istar] = this->CellWidth[0][i]*rnum + this->CellLeftEdge[0][i];

                rnum =  (float) (random() % max_random) / (float) (max_random);
                ParticlePosition[1][istar] = this->CellWidth[1][j]*rnum + this->CellLeftEdge[1][j];

                rnum =  (float) (random() % max_random) / (float) (max_random);
                ParticlePosition[2][istar] = this->CellWidth[2][k]*rnum + this->CellLeftEdge[2][k];

                // assume velocity dispersion is isotropic in each velocity component. Multiply disp by
                // sqrt(1/3) to get disp in each component... taking above velocities as the mean
                rnum  =  (float) (random() % max_random) / (float) (max_random);
                rsign = rnum>0.5 ? 1:-1;
                ParticleVelocity[0][istar] = umean + rsign * GaussianRandomVariable() * IndividualStarVelocityDispersion * 0.577350269*1.0E5*(TimeUnits)/(LengthUnits);

                rnum  =  (float) (random() % max_random) / (float) (max_random);
                rsign = rnum>0.5 ? 1:-1;
                ParticleVelocity[1][istar] = vmean + rsign * GaussianRandomVariable() * IndividualStarVelocityDispersion * 0.577350269*1.0e5*(TimeUnits)/(LengthUnits);

                rnum  =  (float) (random() % max_random) / (float) (max_random);
                rsign = rnum>0.5 ? 1:-1;
                ParticleVelocity[2][istar] = wmean + rsign * GaussianRandomVariable() * IndividualStarVelocityDispersion * 0.577350269*1.0E5*(TimeUnits)/(LengthUnits);

                // ENSURE MOMENTUM CONSERVATION!!!!!
                // make running total of momentum in each direction
                px += ParticleVelocity[0][istar]*ParticleMass[istar];
                py += ParticleVelocity[1][istar]*ParticleMass[istar];
                pz += ParticleVelocity[2][istar]*ParticleMass[istar];

                // We did the metallicity tagging already, but now loop through and
                // do individual chemical tagging for each species tracked in the simulation
                // these are stored as particle attributes starting with attr number 5 (index 4)
                if (((TestProblemData.MultiMetals == 2) || (MultiMetals == 2))){
//                  if (IndividualStarOutputChemicalTags){
                    // output the particle formation time, birth mass (in Msun), and metallicity
//                    printf(" %"ISYM" %"ESYM" %"ESYM" %"ESYM, ParticleType[istar],
//                                                     ParticleAttribute[0][istar], ParticleAttribute[3][istar],
//                                                     ParticleAttribute[2][istar]);
//                  }
                  for( int iyield = 0; iyield < StellarYieldsNumberOfSpecies; iyield++){
                    double temp_fraction = 0.0;

                    if(StellarYieldsAtomicNumbers[iyield] > 2){
                      int field_num;
                      this->IdentifyChemicalTracerSpeciesFieldsByNumber(field_num, StellarYieldsAtomicNumbers[iyield]);

                      if (IndividualStarOutputChemicalTags){
//                        printf(" %"ESYM, BaryonField[field_num][index]);
                        StellarAbundances[iyield][istar]     = BaryonField[field_num][index];
                      } else {
                        ParticleAttribute[4 + iyield][istar] = BaryonField[field_num][index];
                      }

                    } else if (StellarYieldsAtomicNumbers[iyield] == 1){
                      // Take H and He fractions as TOTAL amount of H in the cell
                      // this is probably not needed since it should all be HI in a star forming region anyway
                      temp_fraction = BaryonField[HINum][index] + BaryonField[HIINum][index];

                      if (MultiSpecies > 1){
                        temp_fraction += BaryonField[HMNum][index] +
                                         BaryonField[H2INum][index] + BaryonField[H2IINum][index];
                      }

                      if (IndividualStarOutputChemicalTags){
//                        printf(" %"ESYM, temp_fraction);
                        StellarAbundances[iyield][istar]     = temp_fraction;
                      } else {
                        ParticleAttribute[4 + iyield][istar] = temp_fraction;
                      }
                    } else if (StellarYieldsAtomicNumbers[iyield] == 2){
                      // Again, total amount of Helium - probably not necessary, should all be HeI anyway
                      temp_fraction = BaryonField[HeINum][index] + BaryonField[HeIINum][index] +
                                      BaryonField[HeIIINum][index];

                      if (IndividualStarOutputChemicalTags){
//                        printf(" %"ESYM, temp_fraction);
                        StellarAbundances[iyield][istar]     = temp_fraction;
                      } else{
                        ParticleAttribute[4 + iyield][istar] = temp_fraction;
                      }
                    }
                  }// loop over yields

                  if (IndividualStarOutputChemicalTags){
                    int offset = 0;
                    // add two to the stellar abundances for popIII and AGB metals (if tracked)
                    if (IndividualStarTrackAGBMetalDensity){
                      StellarAbundances[StellarYieldsNumberOfSpecies + 0][istar] = BaryonField[AGBMetalNum][index];
                      offset++;
                    }

                    if (IndividualStarPopIIIFormation){
                      StellarAbundances[StellarYieldsNumberOfSpecies + offset][istar] = BaryonField[PopIIIMetalNum][index];
                      offset++;
                    }

                    if (IndividualStarTrackSNMetalDensity) {
                      StellarAbundances[StellarYieldsNumberOfSpecies + offset][istar] = BaryonField[SNIaMetalNum][index];
                      offset++;
                      StellarAbundances[StellarYieldsNumberOfSpecies + offset][istar] = BaryonField[SNIIMetalNum][index];
                      // offset++;
                    }
                  } // if output

//                  if (IndividualStarOutputChemicalTags) printf("\n");
                } // check multimetals

                if (IndividualStarSaveTablePositions && (ParticleType[istar] == -IndividualStar)){
                  int tstart = ParticleAttributeTableStartIndex;

                  // stellar evolution table (attr 3 = birth mass, attr 2 = metallicity)
                  int t_i = -1, t_j = -1, t_k = -1;
                  IndividualStarGetSETablePosition(t_i, t_j,
                                                   ParticleAttribute[3][istar], ParticleAttribute[2][istar]);
                  ParticleAttribute[tstart    ][istar] = t_i;
                  ParticleAttribute[tstart + 1][istar] = t_j;
                  // radiation properties table (only do if particle can radiate - saves time)
                  if( ParticleAttribute[3][istar] >= IndividualStarRadiationMinimumMass){
                    float Teff, R;
                    IndividualStarInterpolateProperties(Teff, R, (int)ParticleAttribute[tstart][istar],
                                                       (int)ParticleAttribute[tstart+1][istar],
                                                       ParticleAttribute[3][istar], ParticleAttribute[2][istar]);
                    float g = IndividualStarSurfaceGravity(ParticleAttribute[3][istar], R);

                    t_i = -1; t_j = -1; t_k = -1;
                    IndividualStarGetRadTablePosition(t_i, t_j, t_k,
                                                    Teff, g, ParticleAttribute[2][istar]);
                    ParticleAttribute[tstart + 2][istar] = t_i;
                    ParticleAttribute[tstart + 3][istar] = t_j;
                    ParticleAttribute[tstart + 4][istar] = t_k;
                  } else {
                    ParticleAttribute[tstart + 2][istar] = -1;
                    ParticleAttribute[tstart + 3][istar] = -1;
                    ParticleAttribute[tstart + 4][istar] = -1;
                  }

                  // yields table position
                  t_i = -1 ; t_j = -1;
                  StellarYieldsGetYieldTablePosition(t_i, t_j,
                                                   ParticleAttribute[3][istar], ParticleAttribute[2][istar]);
                  ParticleAttribute[tstart + 5][istar] = t_i;
                  ParticleAttribute[tstart + 6][istar] = t_j;

                } // end check for saving table positions

                /* Keeping this as a particle attribute */
                ParticleAttribute[NumberOfParticleAttributes-2][istar] = 0.0; // wind mass ejected
                ParticleAttribute[NumberOfParticleAttributes-1][istar] = 0.0; // sn mass ejected


              } // end while loop for assigning particle properties
/*
              if (add_unresolved_star){
                ParticleAttribute[0][ii] = this->Time;
                ParticleAttribute[2][ii] = BaryonField[MetalNum][index];

                ParticleAttribute[3][ii] = ParticleMass[ii];
                ParticleMass[ii]         = ParticleMass[ii] * msolar / MassUnits;
                ParticlePosition[0][ii]  = this->CellWidth[0][i] + this->CellLeftEdge[0][i];
                ParticlePosition[1][ii]  = this->CellWidth[1][j] + this->CellLeftEdge[1][j];
                ParticlePosition[2][ii]  = this->CellWidth[2][k] + this->CellLeftEdge[2][k];
                ParticleVelocity[0][ii]  = umean;
                ParticleVelocity[1][ii]  = vmean;
                ParticleVelocity[2][ii]  = wmean;

                px += ParticleVelocity[0][ii]*ParticleMass[ii];
                py += ParticleVelocity[1][ii]*ParticleMass[ii];
                pz += ParticleVelocity[2][ii]*ParticleMass[ii];

                ParticleAttribute[NumberOfParticleAttributes-2][ii] = 0.0;
                ParticleAttribute[NumberOfParticleAttributes-1][ii] = 0.0;

                if (((TestProblemData.MultiMetals == 2) || (MultiMetals == 2))){
                  if (IndividualStarOutputChemicalTags){
                    // output the particle formation time, birth mass (in Msun), and metallicity
                    printf(" %"ISYM" %"ESYM" %"ESYM" %"ESYM, ParticleType[ii],
                                                     ParticleAttribute[0][ii], ParticleAttribute[3][ii],
                                                     ParticleAttribute[2][ii]);
                  }
                  for( int iyield = 0; iyield < StellarYieldsNumberOfSpecies; iyield++){
                    double temp_fraction = 0.0;

                    if(StellarYieldsAtomicNumbers[iyield] > 2){
                      int field_num;
                      this->IdentifyChemicalTracerSpeciesFieldsByNumber(field_num, StellarYieldsAtomicNumbers[iyield]);

                      if (IndividualStarOutputChemicalTags){
                        printf(" %"ESYM, BaryonField[field_num][index]);
                      } else {
                        ParticleAttribute[4 + iyield][ii] = BaryonField[field_num][index];
                      }

                    } else if (StellarYieldsAtomicNumbers[iyield] == 1){
                      // Take H and He fractions as TOTAL amount of H in the cell
                      // this is probably not needed since it should all be HI in a star forming region anyway
                      temp_fraction = BaryonField[HINum][index] + BaryonField[HIINum][index];

                      if (MultiSpecies > 1){
                        temp_fraction += BaryonField[HMNum][index] +
                                         BaryonField[H2INum][index] + BaryonField[H2IINum][index];
                      }

                      if (IndividualStarOutputChemicalTags){
                        printf(" %"ESYM, temp_fraction);
                      } else {
                        ParticleAttribute[4 + iyield][ii] = temp_fraction;
                      }
                    } else if (StellarYieldsAtomicNumbers[iyield] == 2){
                      // Again, total amount of Helium - probably not necessary, should all be HeI anyway
                      temp_fraction = BaryonField[HeINum][index] + BaryonField[HeIINum][index] +
                                      BaryonField[HeIIINum][index];

                      if (IndividualStarOutputChemicalTags){
                        printf(" %"ESYM, temp_fraction);
                      } else{
                        ParticleAttribute[4 + iyield][ii] = temp_fraction;
                      }
                    }
                  }

                  if (IndividualStarOutputChemicalTags) printf("\n");
                } // check multimetals

                /////////////////////
                ii++; // because we didn't iterate this before
              }
*/
              // ---------------------------------------------------

              // ensure zero net momentum from mean velocity
              // momentum of gas converted into stars (sum_mass * umean)
              // should be equal to the total momentum of the stars
              // compute excess momentum and modify star velocity evenly (mass weighted)
              // this is not completely physical, as pre-SF and post-SF gas vel is the same
              sum_mass = sum_mass * msolar / MassUnits; // in code units


              px_excess = umean * sum_mass + px;
              py_excess = vmean * sum_mass + py;
              pz_excess = wmean * sum_mass + pz;

              // remove or add momentum evenly from each star if needed
              if ( abs(px_excess) > tiny_number) {
                for (istar = index_presf; istar < ii; istar++){
                  ParticleVelocity[0][istar] += (-1.0 * px_excess) / (ParticleMass[istar] * (float) (ii-index_presf));
                }
              }
              if ( abs(py_excess) > tiny_number){
                for (istar = index_presf; istar < ii; istar++){
                  ParticleVelocity[1][istar] = (-1.0 * py_excess) / (ParticleMass[istar] * (float) (ii-index_presf));
                }
              }
              if ( abs(pz_excess) > tiny_number){
                for (istar = index_presf; istar < ii; istar++){
                  ParticleVelocity[2][istar] = (-1.0 * pz_excess) / (ParticleMass[istar] * (float) (ii-index_presf));
                }
              }

              // now remove mass from grid - do not need to do this for tracer fields since they are kept as fractions
              // and will be modified accordingly when converted back to densities in Grid_StarParticleHandler
              l = 0;
              for (int k_loc = -kstart; k_loc <= kend; k_loc++){
                for(int j_loc = -jstart; j_loc <= jend; j_loc++){
                  for (int i_loc = -istart; i_loc <= iend; i_loc++){
                    int loc_index = (i + i_loc) + ( (j + j_loc) + (k + k_loc)*(ny))*(nx);

                    if(BaryonField[DensNum][loc_index] > secondary_odthreshold &&
                       temp[loc_index] < IndividualStarTemperatureThreshold    ){
                      // mass is removed as weighted by the previous cell mass (more mass is
                      // taken out of higher density regions). M_new = M_old - M_sf * (M_old / M_tot)
                      // where M_tot is mass of cells that meet above SF conditions. Simplifies to below eq:
                      float old_density = BaryonField[DensNum][loc_index];

                      BaryonField[DensNum][loc_index] *= (1.0 - sum_mass / bmass);

                      // adjust total energy if we are using PPM
                      if (HydroMethod != 2 && FALSE){
                          float ke_after, delta_ke;

                          ke_after = 0.5 * BaryonField[DensNum][loc_index] *
                                   ( BaryonField[Vel1Num][loc_index] * BaryonField[Vel1Num][loc_index] +
                                     BaryonField[Vel2Num][loc_index] * BaryonField[Vel2Num][loc_index] +
                                     BaryonField[Vel3Num][loc_index] * BaryonField[Vel3Num][loc_index]);

                          delta_ke = ke_after - ke_before[ l ];

                           // or TE_new = TE_old/(1.0 -sum_mass / bmass) + delta_ke / new_density;
                          BaryonField[TENum][loc_index] = (BaryonField[TENum][loc_index]*old_density +
                                                                  delta_ke) / BaryonField[DensNum][loc_index];
                      }
                      l++;
                    }
                  }
                }
              }

            } // if jeans mass unstable
          } // resolution and density


        } // enx x loop
      } // end y loop
    } // end z loop

    // Done forming stars!!! Output and exit
    if (ii > 0){
      printf("P(%"ISYM"): individual_star_maker[add]: %"ISYM" new star particles\n", MyProcessorNumber, ii);
    }
    if (ii >= *nmax){
      fprintf(stdout, "individual_star_maker: reached max new particle count!! Available: %"ISYM". Made: %"ISYM"\n", *nmax, ii);

    }

  // star masses are recorded as densities (mass / cell volume)
  // set progenitor masses in solar
  for (int counter = 0; counter < ii; counter++){
    ParticleMass[counter]   = ParticleMass[counter] / (dx*dx*dx); // code units / cell volume
  }

  *np = ii; // number of stars formed : AJE 2/29 check if this is a bug with the -1


  delete[] ke_before;

  return SUCCESS;
}


float SamplePopIII_IMF(void){

  return SampleIMF(SecondaryIMFData,
                   PopIIILowerMassCutoff, PopIIIUpperMassCutoff);
}

float SampleIMF(void){
// to retain default behavoir

  return SampleIMF(IMFData,
                   IndividualStarIMFLowerMassCutoff,
                   IndividualStarIMFUpperMassCutoff);

}

float SampleIMF(float * data, const float & lower_mass, const float & upper_mass){
/*-----------------------------------------------------------------------------
  SampleIMF

  Samples the tabulated initial mass function a la cumulative probablity density
  as created in StarParticleIndividual_IMFInitialize.

  INPUTS -
    None

  OUTPUTS -
    Mass of randomly selected star in solar masses

-----------------------------------------------------------------------------*/

  unsigned_long_int random_int = mt_random();
  const int max_random = (1<<16);
  float x = (float) (random_int%max_random) / (float) (max_random);
  float dm = log10(upper_mass / lower_mass )/ ((float) (IMF_TABLE_ENTRIES-1));
  float m;

  int bin_number;

  if (x <= data[0] ){
    bin_number = 0;
  } else if (x >= data[IMF_TABLE_ENTRIES-1]){
    bin_number = IMF_TABLE_ENTRIES - 1;
  } else{
    bin_number = search_lower_bound(data, x, 0, IMF_TABLE_ENTRIES, IMF_TABLE_ENTRIES);
  }

  m = lower_mass * POW(10.0, bin_number * dm);

  IndividualStarIMFCalls++;

  return m;
}

int grid::IndividualStarSetWDLifetime(void){
/* --------------------------------------------------
 * IndividualStarSetWDLifetime
 * --------------------------------------------------
 * Updates WD lifetimes if not yet initialized using
 * DTD SNIa model
 * --------------------------------------------------
 */

  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits, MassUnits;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
              &TimeUnits, &VelocityUnits, this->Time) == FAIL){
      ENZO_FAIL("Error in GetUnits");
  }


  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  if (NumberOfParticles == 0)
    return SUCCESS;

  for (int i = 0; i < NumberOfParticles; i++){

    if( ParticleType[i] != -PARTICLE_TYPE_INDIVIDUAL_STAR_WD ){
      continue;
    }

    //
    // lifetime is set relative to WD formation time (now)
    //
    float new_lifetime = -1;

    int result = SetWDLifetime(new_lifetime, this->Time, ParticleAttribute[0][i],
                                            ParticleAttribute[1][i], TimeUnits);

    ParticleAttribute[1][i] = new_lifetime;
    ParticleType[i]         = ABS(ParticleType[i]);

    if (ParticleAttribute[1][i] < 0){
      return FAIL;
    }
    //
    // feedback operates computing death time = lifetime + birth time
    // renormalize so as to keep birth time the original star particle birth time
    //  - original lifetime of progenitor star to WD can be backed out via postprocessing, but not birth time
    //
    if (result > 0){ // negative result means WD never exploding  -- ensure it is not this timestep
      ParticleAttribute[1][i] = fmax(new_lifetime,1.5*this->dtFixed) + (this->Time - ParticleAttribute[0][i]);
    }
  }


  return SUCCESS;
}

void ComputeStellarWindMassLossRate(const float &mproj, const float &metallicity,
                                    float *dMdt){
 /* ------------------------------------------------------------------
  * ComputeStellarWindEjectaMass
  * -------------------------------------------------------------------
  * A. Emerick - 4/22/16
  *
  * Model for stellar wind mass loss rate taken from Leitherer et. al. 1992.
  * This is the same model used in STARBURST 99 stellar wind models.
  * -------------------------------------------------------------------- */

  float L, Teff, Z, R;

  const double solar_z = 0.02; // as assumed in Leithener et. al. 1992
  const double yr      = 3.16224E7; // number of seconds in a year

  /* get properties */
  if(IndividualStarInterpolateLuminosity(L, mproj, metallicity) == FAIL){
    ENZO_FAIL("ComputeStellarWindMassLossRate: Failed to interpolate luminosity");
  }

  if(IndividualStarInterpolateProperties(Teff, R, mproj, metallicity) == FAIL){
    ENZO_FAIL("ComputeStellarWindMassLossRate: Failed to interpolate stellar properties");
  }

  /* compute logged mass loss rate */
  *dMdt = -24.06 + 2.45 * log10(L) - 1.10 * log10(mproj) + 1.31 * log10(Teff)
                                   + 0.80 * log10(metallicity / solar_z);

  *dMdt = POW(10.0, *dMdt) / yr ; // Msun / yr -> Msun / s
  return;

}

void ComputeStellarWindVelocity(Star *cstar, float *v_wind){
 /* ------------------------------------------------------------------
  * ComputeStellarWindVelocity
  * -------------------------------------------------------------------
  * A. Emerick - 4/22/16
  *
  * Model for stellar wind velocities taken from Leitherer et. al. 1992.
  * This is the same model used in STARBURST 99 stellar wind models.
  * The mass loss rate is computed elsewhere from stellar yields tables,
  * but velocity is set below using the fit function in luminosity,
  * stellar mass, effective temperature, and metallicity
  * -------------------------------------------------------------------- */

  float L, Teff, Z, R;

  const double solar_z = 0.02; // as assumed in Leithener et. al. 1992

  int* se_table_position = cstar->ReturnSETablePosition();

  /* get properties */
  IndividualStarInterpolateLuminosity(L, se_table_position[0], se_table_position[1],
                                          cstar->ReturnBirthMass(), cstar->ReturnMetallicity());
  IndividualStarInterpolateProperties(Teff, R, se_table_position[0], se_table_position[1],
                                          cstar->ReturnBirthMass(), cstar->ReturnMetallicity());

  // wind is in units of km / s
  // L - solar units
  // T - Kelvin
  // M - solar units
  // Z - solar units
  *v_wind =   1.23 - 0.30*log10(L) + 0.55 * log10(cstar->ReturnBirthMass())
            + 0.64 * log10(Teff) + 0.13*log10(cstar->ReturnMetallicity()/solar_z);
  *v_wind = POW(10.0, *v_wind);

  return;
}


int grid::IndividualStarAddFeedbackSphere(Star *cstar, float *mp, const int mode){

/*
     General function to add feedback for a given star in a spherical region

     Thermal energy injection ONLY

     mode   :   integer, (Values)
                switches between stellar wind, core collapse SN, or type ia sn

*/
  if (this->NumberOfBaryonFields == 0 || !this->isLocal() )
    return SUCCESS;

  float dx = this->CellWidth[0][0];

  float m_eject, E_thermal;

  const float mproj = cstar->ReturnBirthMass();
  const float lifetime = cstar->ReturnLifetime();

  const double msolar = 1.989E33;

  float *metal_mass; // array of individual species masses

  FLOAT * pos;
  pos = cstar->ReturnPosition();

  /* Get Units */
  float DensityUnits, LengthUnits, TemperatureUnits,
        TimeUnits, VelocityUnits, MassUnits, EnergyUnits;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
               &TimeUnits, &VelocityUnits, this->Time) == FAIL){
      ENZO_FAIL("Error in GetUnits");
  }
  MassUnits   = DensityUnits*LengthUnits*LengthUnits*LengthUnits;
  EnergyUnits = MassUnits * VelocityUnits * VelocityUnits;

  /* If we are following yields, initialize array to hold ejecta masses */
  if(IndividualStarFollowStellarYields && ((TestProblemData.MultiMetals == 2) || (MultiMetals == 2))){

    metal_mass = new float[StellarYieldsNumberOfSpecies + 1];

    for (int i = 0; i < StellarYieldsNumberOfSpecies + 1; i ++){
      metal_mass[i] = 0.0;
    }

  } else { metal_mass = NULL;}


  float cgs_lifetime = lifetime * TimeUnits;

  int stellar_wind_mode = FALSE;

  if( mode < 0 ){  // compute properties for stellar wids

    // mproj needs to be in Msun - everything else in CGS
    IndividualStarSetStellarWindProperties(cstar, this->Time, this->dtFixed, TimeUnits,
                                           m_eject, E_thermal, metal_mass);
    stellar_wind_mode = TRUE;
  } else if (mode == 1){

    // core collapse supernova
    IndividualStarSetCoreCollapseSupernovaProperties(cstar, m_eject, E_thermal, metal_mass);

    stellar_wind_mode = FALSE;
  } else if (mode == 2){

    // Type Ia supernova properties
    IndividualStarSetTypeIaSupernovaProperties(m_eject, E_thermal, metal_mass);
    // printf("m_eject  for snia = %"FSYM"\n", m_eject);
    stellar_wind_mode = FALSE;

  } else if (mode == 3){

    IndividualStarSetPopIIISupernovaProperties(cstar, m_eject, E_thermal, metal_mass);

    stellar_wind_mode = FALSE;
  }

  /* convert computed parameters to code units */
  m_eject   = m_eject*msolar / MassUnits   / (dx*dx*dx);
  E_thermal = E_thermal      / EnergyUnits / (dx*dx*dx);

  if(IndividualStarFollowStellarYields && ((TestProblemData.MultiMetals == 2) || (MultiMetals == 2))){
    for(int i = 0; i < StellarYieldsNumberOfSpecies + 1; i++){
      // printf("metal mass species %"ISYM"   = %"ESYM"\n", i, metal_mass[i]);
      metal_mass[i] = metal_mass[i] * msolar / MassUnits / (dx*dx*dx);
    }
  }

  //
  // now that we've computed the explosion properties
  // find where we should go off
  //
  if( (m_eject > 0) || (E_thermal > 0)){ // can sometimes both be zero for stellar winds due to mass corrections
    this->IndividualStarInjectSphericalFeedback(cstar, pos[0], pos[1], pos[2], m_eject, E_thermal,
                                                metal_mass, stellar_wind_mode);
  }

  float new_mass = (*mp) - m_eject * (dx*dx*dx) * MassUnits / msolar; // update mass

  if( *mp < 0 && mode != 2){ // This can happen for Type 1a since using fixed mass model
    printf("new_mass = %"ESYM" mp = %"ESYM" m_eject =%"ESYM"\n", new_mass, *mp, m_eject*dx*dx*dx*MassUnits/msolar);
    ENZO_FAIL("IndividualStarFeedback: Ejected mass greater than current particle mass - negative particle mass!!!\n");
  } else if (mode == 2){
    *mp = 0.0;
  } else{
    *mp = new_mass;
  }

  delete [] metal_mass;

  return SUCCESS;
}

int grid::IndividualStarInjectSphericalFeedback(Star *cstar,
                                                const FLOAT xp, const FLOAT yp, const FLOAT zp,
                                                const float m_eject, const float E_thermal,
                                                const float *metal_mass, const int stellar_wind_mode){

  float dx = float(this->CellWidth[0][0]); // for convenience

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, CRNum, B1Num, B2Num, B3Num;

  this->DebugCheck("StarParticleHandler");
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }

  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;

  if ( MultiSpecies ) {
    IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
                          HMNum, H2INum, H2IINum, DINum, DIINum, HDINum);
  }

  int PopIIIMetalNum, AGBMetalNum, SNIaMetalNum, SNIIMetalNum;

  AGBMetalNum    = FindField(ExtraType0, FieldType, NumberOfBaryonFields);
  PopIIIMetalNum = FindField(ExtraType1, FieldType, NumberOfBaryonFields);
  SNIaMetalNum   = FindField(MetalSNIaDensity, FieldType, NumberOfBaryonFields);
  SNIIMetalNum   = FindField(MetalSNIIDensity, FieldType, NumberOfBaryonFields);

  if ( IndividualStarTrackAGBMetalDensity && (AGBMetalNum <= 0)){
    ENZO_FAIL("Error in finding AGB metal density field in individual_star_maker");
  }

  if ( IndividualStarPopIIIFormation && (PopIIIMetalNum <= 0)){
    ENZO_FAIL("Error in finding Pop III metal density field in individual_star_maker");
  }

  if ( IndividualStarTrackSNMetalDensity && ( (SNIIMetalNum <= 0) || (SNIaMetalNum <=0))){
    ENZO_FAIL("Error in finding SNII and SNIa metal density field in individual_star_maker.");
  }

  /* Get Units */
  float DensityUnits, LengthUnits, TemperatureUnits,
        TimeUnits, VelocityUnits, MassUnits, EnergyUnits;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
               &TimeUnits, &VelocityUnits, this->Time) == FAIL){
      ENZO_FAIL("Error in GetUnits");
  }
  MassUnits   = DensityUnits*LengthUnits*LengthUnits*LengthUnits;
  EnergyUnits = MassUnits * VelocityUnits * VelocityUnits;

  FLOAT xstart = this->CellLeftEdge[0][0],
        ystart = this->CellLeftEdge[1][0],
        zstart = this->CellLeftEdge[2][0];

  const int nx = *(this->GridDimension), ny = *(this->GridDimension+1), nz = *(this->GridDimension+2);
  const int size = nx*ny*nz;

  //
  // find position and index for grid zone
  // nearest to particle
  //
  const float xpos = (xp - xstart)/dx;
  const float ypos = (yp - ystart)/dx;
  const float zpos = (zp - zstart)/dx;

  const int ic   = ((int) floor(xpos ));
  const int jc   = ((int) floor(ypos ));
  const int kc   = ((int) floor(zpos ));

  float * temperature;
  temperature = new float[size];


  if(  this->ComputeTemperatureField(temperature) == FAIL ){
    ENZO_FAIL("Error in compute temperature called from PhotoelectricHeatingFromStar");
  }

  //
  //
  // now move outward from particle computing
  // fractional volume to inject
  //
  const float radius = IndividualStarFeedbackStencilSize * dx;    // code length
  const float volume = 4.0 * pi * radius * radius * radius / 3.0; // (code length)**3
  const int   r_int  = ceil(IndividualStarFeedbackStencilSize);     // int of farthest cell in any dir.
  const float cell_volume_fraction = dx*dx*dx / volume;           // fraction of vol for each cell

  float injected_metal_mass[StellarYieldsNumberOfSpecies+1];

  if (metal_mass == NULL && (cstar)){
    injected_metal_mass[0] = cstar->ReturnMetallicity() * m_eject;
  } else if (metal_mass == NULL){
    injected_metal_mass[0] = 0.0;
  }

  // for printing stats at the end
  float total_volume_fraction = 0.0, total_grid_mass = 0.0, total_mass_injected = 0.0;
  float total_energy_injected = 0.0, max_density_on_grid = 0.0, average_density_on_grid = 0.0;
  float total_metal_mass = 0.0;
  int   cells_this_grid = 0;

  // loop over cells, compute fractional volume, inject feedback
  for(int k = kc - r_int; k <= kc + r_int; k++){
    if ( (k<0) || (k>=nz) ) continue;
    for(int j = jc - r_int; j <= jc + r_int; j++){
      if ( (j<0) || (j>=ny) ) continue;
      for(int i = ic - r_int; i <= ic + r_int; i++){
        if ( (i<0) || (i>=nx) ) continue;

        int index = i + (j + k*ny)*nx;
        float fractional_overlap = 1.0;

        // off grid - just means particle is near a boundary
        // feedback handled appropriately on corresponding grid
        if( index < 0 || index >= nx*ny*nz) continue;

        float xc,yc,zc;
        xc = (i + 0.5) * dx + xstart;
        yc = (j + 0.5) * dx + ystart;
        zc = (k + 0.5) * dx + zstart;

        // do fractional overlap calculation
        fractional_overlap = ComputeOverlap(1, radius, xp, yp, zp,
                                            xc - dx, yc - dx, zc - dx,
                                            xc + dx, yc + dx, zc + dx, IndividualStarFeedbackOverlapSample);

        if(fractional_overlap <= 0.0) continue; // cell is enirely outside sphere

        float injection_factor = cell_volume_fraction * fractional_overlap;

        float delta_mass  = m_eject   * injection_factor;
        float delta_therm = E_thermal * injection_factor;

        if (injection_factor < 0) {ENZO_FAIL("injection factor < 0");}

        if (IndividualStarFollowStellarYields && (cstar || metal_mass)){
          for(int im = 0; im < StellarYieldsNumberOfSpecies+1; im++){

            // Hack  - remove contribution from stars to experiment fields
            //          remember, im = 0 is the total metal mass field
            if ((cstar) && (im > 0) && MetalMixingExperiment){
              for (int j = 0; j < StellarYieldsNumberOfSpecies; j++){
                if (StellarYieldsAtomicNumbers[im-1] == MixingExperimentData.anums[j] ){
                  injection_factor = 0.0; // zero out these fields to limit to just mixing experiments events ONLY
                }
              }
            }

            injected_metal_mass[im] = metal_mass[im]*injection_factor;
          }
        }

        if(stellar_wind_mode){
          // Stellar winds are challenging - to avoid both very high velocities
          // on the grid and superheated gas, modify feedback accordingly
          ModifyStellarWindFeedback(BaryonField[DensNum][index],
                                    temperature[index], dx, MassUnits, EnergyUnits,
                                    delta_mass, delta_therm, injected_metal_mass,
                                    this->AveragedAbundances);
        }

        float inv_dens = 1.0 / (BaryonField[DensNum][index] + delta_mass);

        /* record statistics accumulators */
        cells_this_grid++;
        total_volume_fraction += injection_factor;
        total_mass_injected   += delta_mass;
        total_energy_injected += delta_therm;
        total_grid_mass       += BaryonField[DensNum][index];
        max_density_on_grid     = fmax( BaryonField[DensNum][index], max_density_on_grid);

        BaryonField[TENum][index] = (BaryonField[TENum][index] * BaryonField[DensNum][index]
                                     + delta_therm) * inv_dens;
        if(DualEnergyFormalism){
          BaryonField[GENum][index] = (BaryonField[GENum][index] * BaryonField[DensNum][index]
                                       + delta_therm) * inv_dens;
        }
        float old_mass = BaryonField[DensNum][index];
        BaryonField[DensNum][index] += delta_mass;

        /* add metal species if we need to */
        if(((TestProblemData.MultiMetals == 2) || (MultiMetals == 2)) && IndividualStarFollowStellarYields){
          int field_num;
          this->IdentifyChemicalTracerSpeciesFieldsByNumber(field_num, 0); // gives metallicity field


          if (cstar){ // is feedback coming from a star particle?
            BaryonField[field_num][index] += injected_metal_mass[0]; // add to metallicity field

            // Add to separate source fields if they exist
            if (IndividualStarPopIIIFormation && cstar->ReturnMetallicity() < PopIIIMetalCriticalFraction){
              BaryonField[PopIIIMetalNum][index] += injected_metal_mass[0];

            } else if (IndividualStarTrackAGBMetalDensity &&
                      (cstar->ReturnBirthMass() < IndividualStarSNIIMassCutoff) &&
                      (stellar_wind_mode)){                                 // make sure we are an AGB wind
              BaryonField[AGBMetalNum][index] += injected_metal_mass[0];

            } else if (IndividualStarTrackSNMetalDensity && !(stellar_wind_mode)){

              if (cstar->ReturnBirthMass() > IndividualStarSNIIMassCutoff){
                BaryonField[SNIIMetalNum][index] += injected_metal_mass[0];

              } else if (cstar->ReturnBirthMass() < IndividualStarSNIaMaximumMass){
                BaryonField[SNIaMetalNum][index] += injected_metal_mass[0];

              }

            } // end if tracking yield source modes

          } else if (metal_mass){
            BaryonField[field_num][index] += injected_metal_mass[0];

          } else {
            // keep same fraction if using artificial SN generaotr
            BaryonField[field_num][index] += delta_mass *
                                             BaryonField[field_num][index] / old_mass;
          }

            total_metal_mass += BaryonField[field_num][index];

          for(int im = 0; im < StellarYieldsNumberOfSpecies; im++){
            this->IdentifyChemicalTracerSpeciesFieldsByNumber(field_num,
                                                              StellarYieldsAtomicNumbers[im]);
            if (cstar || metal_mass){
              BaryonField[field_num][index] += injected_metal_mass[1 + im];
            } else { // keep same fraction if using artificial SN generator
              BaryonField[field_num][index] += delta_mass *
                                               BaryonField[field_num][index] / old_mass;
            }
          }

        } else{
          int field_num;
          this->IdentifyChemicalTracerSpeciesFieldsByNumber(field_num, 0); // gives metallicity field

          if (cstar || metal_mass){
            BaryonField[field_num][index] += injected_metal_mass[0];
          } else{
            BaryonField[field_num][index] += delta_mass *
                                             BaryonField[field_num][index]/ old_mass;
          }
        } // end yields check

      }
    }
  }


  // print SN stats to check if resolved if desired
  if (IndividualStarPrintSNStats && (!stellar_wind_mode) && (cstar || metal_mass)){
    // Column order: Grid ID, Particle ID, M_now, M_eject, Sphere Volume

    average_density_on_grid = total_grid_mass / (1.0 * cells_this_grid); // Sum Density / # cells

    /* convert to CGS */
    total_grid_mass         *= dx*dx*dx * MassUnits;
    total_mass_injected     *= dx*dx*dx * MassUnits;
    total_energy_injected   *= dx*dx*dx * EnergyUnits;
    total_metal_mass        *= dx*dx*dx * MassUnits;
    const float volume_cgs   = volume / (LengthUnits * LengthUnits * LengthUnits);
    max_density_on_grid     *= DensityUnits;
    average_density_on_grid *= DensityUnits;
    const float m_eject_cgs  = m_eject * dx*dx*dx * MassUnits;
    const float average_metallicity      = total_metal_mass / (total_mass_injected + total_grid_mass);

    /* compute Sedov-Taylor phase radius (R_dps) */

    if (cstar){
      printf("IndividualStarSNStats: %"ISYM" %"ISYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ISYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM"\n",
            this->ID, cstar->ReturnID(), this->Time, cstar->ReturnMass(), cstar->ReturnBirthMass(), cstar->ReturnMetallicity(), m_eject_cgs,
            cells_this_grid, volume_cgs, total_volume_fraction, total_mass_injected, total_energy_injected,
            total_grid_mass, max_density_on_grid, average_density_on_grid, total_metal_mass, average_metallicity);
    } else {
      printf("IndividualStarSNStats: %"ISYM" %"ISYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ISYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM"\n",
            this->ID, -1.0, this->Time, -1.0, -1.0, -1.0, m_eject_cgs,
            cells_this_grid, volume_cgs, total_volume_fraction, total_mass_injected, total_energy_injected,
            total_grid_mass, max_density_on_grid, average_density_on_grid, total_metal_mass, average_metallicity);

    }
  }


  delete [] temperature;
  // done with spherical injection feedback

  return SUCCESS;
}

void IndividualStarSetPopIIISupernovaProperties(Star *cstar, float &m_eject, float &E_thermal, float *metal_mass){
/* -------------------------------------------------------
 * IndividualStarSetPopIIISupernovaProperties
 * -------------------------------------------------------
 * A. Emerick - Oct 2018
 *
 * Set the ejected mass, energy, and metal masses for a
 * PopIII supernova explosion.
 * -------------------------------------------------------
 */

  int *yield_table_position = cstar->ReturnYieldTablePosition();
  float birth_mass = cstar->ReturnBirthMass();
  /* compute total ejected yield */

  if ( ((birth_mass >= TypeIILowerMass) && (birth_mass <= TypeIIUpperMass)) ||
       ((birth_mass >= PISNLowerMass)   && (birth_mass <= PISNUpperMass)) ){

    if ( IndividualStarFollowStellarYields && ((TestProblemData.MultiMetals == 2) || (MultiMetals == 2))){
      m_eject   = StellarYieldsInterpolatePopIIIYield(yield_table_position[0],
                                                      birth_mass, -1);
    } else{
      m_eject   = StarMassEjectionFraction * cstar->ReturnMass();
    }

    /* metal masses for tracer species */
    if(IndividualStarFollowStellarYields && ((TestProblemData.MultiMetals == 2) || (MultiMetals == 2))){
      metal_mass[0] = StellarYieldsInterpolatePopIIIYield(yield_table_position[0],
                                                          cstar->ReturnBirthMass(),
                                                          0);
      for(int i = 0; i < StellarYieldsNumberOfSpecies; i++){
        metal_mass[1+i] = StellarYieldsInterpolatePopIIIYield(yield_table_position[0],
                                                              cstar->ReturnBirthMass(),
                                                              StellarYieldsAtomicNumbers[i]);
      }
    }

    /* Set energy for normal SN */
    if (birth_mass <= TypeIIUpperMass){
      E_thermal = IndividualStarSupernovaEnergy * 1.0E51;
    } else {

      if (PopIIIPISNEnergy < 0.0){
        // taken from pop3_maker.F (heger and woosley??)
        float he_core    = (13.0 / 24.0) * (birth_mass - 20.0);
        float sne_factor = 5.0 + 1.304 * (he_core - 64.0);
        E_thermal = sne_factor * 1.0E51;
      } else {
        E_thermal = PopIIIPISNEnergy * 1.0E51;
      }

    }

  } else {
    m_eject = 0.0; // no yields -- but we should NEVER have to specify this if
                   //   feedback routines are operating correctly
    for (int i = 0; i <= StellarYieldsNumberOfSpecies; i++){
      metal_mass[i] = 0.0;
    }
    E_thermal = 0.0;
  }

  return;
}

void IndividualStarSetTypeIaSupernovaProperties(float &m_eject, float &E_thermal, float *metal_mass){
/* -------------------------------------------------------
 * IndividualStarSetTypeIaSupernovaProperties
 * -------------------------------------------------------
 * A. Emerick - Oct 2016
 *
 * Set the ejected mass, energy, and metal masses for a
 * Type Ia supernova explosion.
 *
 * Current model treats all Type Ia uniformly
 * -------------------------------------------------------
 */


  const float c_light = 2.99792458E10;

  m_eject = StellarYields_SNIaYieldsByNumber(-1); // total ejecta in solar masses

  /* set energy given user input */
  if( IndividualStarSupernovaEnergy < 0){
    E_thermal = m_eject * StarEnergyToThermalFeedback * (c_light * c_light);
  } else {
    E_thermal = IndividualStarSupernovaEnergy * 1.0E51;
  }

  /* populate metal species array if needed */
  if (IndividualStarFollowStellarYields && ((TestProblemData.MultiMetals == 2) || (MultiMetals == 2))){
    metal_mass[0] = StellarYields_SNIaYieldsByNumber(0); // total metal mass

    for( int i = 0; i < StellarYieldsNumberOfSpecies; i++){
      metal_mass[i+1] = StellarYields_SNIaYieldsByNumber(StellarYieldsAtomicNumbers[i]);
    }
  }

  return;
}

void IndividualStarSetCoreCollapseSupernovaProperties(Star *cstar,
                                                      float &m_eject, float &E_thermal, float *metal_mass){
/* -------------------------------------------------------
 * IndividualStarCoreCollapseSupernovaProperties
 * -------------------------------------------------------
 * A. Emerick - Sep 2016
 *
 * Set the ejected mass, energy, and metal masses for a
 * core collapse supernova, given star's birth mass and
 * metallicity.
 * -------------------------------------------------------
 */

  const float c_light = 2.99792458E10;

  int *yield_table_position = cstar->ReturnYieldTablePosition();

  /* compute total ejected yield */
  if ( IndividualStarFollowStellarYields && ((TestProblemData.MultiMetals == 2) || (MultiMetals == 2))){
    // 0 in first argument signifies use CC supernova yield table
    m_eject   = StellarYieldsInterpolateYield(0, yield_table_position[0], yield_table_position[1],
                                              cstar->ReturnBirthMass(), cstar->ReturnMetallicity(), -1);
  } else{
    m_eject   = StarMassEjectionFraction * cstar->ReturnMass();
  }

  /* Fail if we are injecting a second time */
  if (cstar->ReturnSNMassEjected() > 0.0){
    ENZO_FAIL("Somehow ejected SN mass twice for this particle\n");
  }

  /* set thermal energy of explosion */
  if( IndividualStarSupernovaEnergy < 0){
    E_thermal = m_eject * StarEnergyToThermalFeedback * (c_light * c_light);
  } else{
    E_thermal = IndividualStarSupernovaEnergy * 1.0E51;
  }

  /* metal masses for tracer species */
  if(IndividualStarFollowStellarYields && ((TestProblemData.MultiMetals == 2) || (MultiMetals == 2))){
    metal_mass[0] = StellarYieldsInterpolateYield(0, yield_table_position[0], yield_table_position[1],
                                                  cstar->ReturnBirthMass(), cstar->ReturnMetallicity(), 0);

    for(int i = 0; i < StellarYieldsNumberOfSpecies; i++){
      metal_mass[1+i] = StellarYieldsInterpolateYield(0, yield_table_position[0], yield_table_position[1],
                                                      cstar->ReturnBirthMass(), cstar->ReturnMetallicity(),
                                                      StellarYieldsAtomicNumbers[i]);
    }
  }

  return;
}


void ModifyStellarWindFeedback(float cell_mass, float T, float dx,
                               float MassUnits, float EnergyUnits, float &m_eject,
                               float &E_thermal, float * metal_mass,
                               float *grid_abundances){

/* ----------------------------------------------------------------------------
 * ModifyStellarWindFeedback
 * -----------------------------------------------------------------------------
 * A. Emerick - Sep 2016
 *
 * Stellar winds are very challenging to get right, even at very high resolution.
 * We assume feedback here is used on ~1-5 pc resolution, galaxy scale, Gyr sims,
 * so feedback needs to be modified to be tractable. Inject thermal energy at
 * 100% thermalization of wind KE first (E_thermal_max) using velocities from
 * full wind model. However, if cell temperature gets very hot (~10^6 K, as set
 * by IndividualStarWindTemperature) then employ maximum wind velocity cutoff,
 * changing energy injection to E_thermal_min. If E_thermal_min will large large
 * temperatures above IndividualStarWindTemperature (as happens if soruce region
 * is devoid of gas), we use an ISM mass loading model to account for shell
 * mixing, which will be heineously unresolved at 1 pc resolution (need ~0.01 pc).
 * ----------------------------------------------------------------------------- */

  if (!IndividualStarUseWindMixingModel){
     return;
  }

  const float est_mu  = 0.5; // estimated - this is approximate anyway
  const float k_boltz = 1.380658E-16;
  const float mp      = 1.6733E-24;

  float m_ism = 0.0;

  E_thermal = E_thermal * dx *dx *dx * EnergyUnits;
  m_eject   = m_eject * dx *dx *dx * MassUnits;
  cell_mass = cell_mass *dx*dx*dx*MassUnits;

  float T_final = (E_thermal + 1.5 * cell_mass * k_boltz * T / (est_mu*mp)) *
                  (2.0 * est_mu * mp/(3.0 * k_boltz * (cell_mass + m_eject)));

    if(T_final > IndividualStarWindTemperature || T > IndividualStarWindTemperature){
      /* Compute the mass that needs to be injected */
      float E_final = (3.0/2.0) * (cell_mass/(est_mu*mp))*k_boltz * T + E_thermal;
      T_final = fmax(T, IndividualStarWindTemperature);
      m_ism   = fmax( (E_final * (2.0 * est_mu * mp)/(3.0*k_boltz * T_final)) - cell_mass - m_eject, 0.0);

      /* modify metal abundances here */
      m_ism = m_ism / (dx*dx*dx) / MassUnits;
      if(((TestProblemData.MultiMetals == 2) || (MultiMetals == 2)) && IndividualStarFollowStellarYields && m_ism > 0.0){
        for (int im = 0; im < StellarYieldsNumberOfSpecies + 1; im++){

          metal_mass[im] = metal_mass[im] + m_ism * grid_abundances[im]; // ISM abundances

          if(metal_mass[im] < 0.0){
            printf("metal_mass %"ESYM" %"ISYM" %"ESYM"\n", metal_mass[im], im, grid_abundances[im]);
            ENZO_FAIL("IndividualStarFeedback: Metal mass correction < 0 and m_ism >0");
          }
        }
      } else{
        for (int im = 0; im < StellarYieldsNumberOfSpecies + 1; im++){
          if(metal_mass[im] < 0.0){
            printf("metal_mass %"ESYM" %"ISYM"\n", metal_mass[im], im);
            ENZO_FAIL("IndividualStarFeedback: Metal mass correction < 0 and m_ism < 0");
          }
        }
      }

  }

  /* make sure things aren't whacky */
  if (E_thermal < 0.0 || m_eject < 0.0 || m_ism < 0.0){
    printf("Error in stellar wind calculation. E_thermal = %"ESYM" m_eject = %"ESYM" m_ism (code) = %"ESYM"\n",E_thermal, m_eject, m_ism);
    ENZO_FAIL("IndividualStarFeedback: Negative injection values in stellar wind feedback modification\n");
  }

  /* convert back into code units */
  E_thermal  = E_thermal / (dx*dx*dx) / EnergyUnits;
  m_eject    = (m_eject)  / (dx*dx*dx) / MassUnits + m_ism;

}

float ComputeOverlap(const int &i_shape, const float &radius,
                     const FLOAT &xc, const FLOAT &yc, const FLOAT &zc,
                     const FLOAT &xl, const FLOAT &yl, const FLOAT &zl,
                     const FLOAT &xr, const FLOAT &yr, const FLOAT &zr,
                     const int &nsample){
 /* -------------------------------------------------------------------------
  * ComputeVolumeOverlap
  * -------------------------------------------------------------------------
  * A. Emerick - 9/21/16
  *
  * Computes overlap between a given rectangular grid cell and a spherical or
  * cylindrical volume using a Monte Carlo sampling method.
  *
  * Adopted from Joshua Wall's version of David Clarke's method in ZEUS-MP
  *
  * INPUTS:
  *   i_shape    - switch for shape (1 == sphere, 2 == right cylinder)
  *   xc,yc,zc   - center coordinates of sphere
  *   xl,yl,zl   - coordinates of lower left corner of grid cell
  *   xr,yr,zr   - coordinates of upper right corner of grid cell
  *   nsample    - number of Monte Carlo sample points
  *
  * OUTPUTS:
  *   overal     - fractional volume overlap for given cell
  * ------------------------------------------------------------------------- */

  float xsq[nsample], ysq[nsample], zsq[nsample];


  float factor = 1.0;
  if (i_shape == 2) factor = 0.0; // cylinder

  //
  float dx = (xr - xl) / ((float) nsample);
  float dy = (yr - yl) / ((float) nsample);
  float dz = (zr - zl) / ((float) nsample);

  for( int i = 0; i < nsample; i++){
    xsq[i] = (xl + (0.5 + i)*dx - xc)*(xl + (0.5+i)*dx - xc);
    ysq[i] = (yl + (0.5 + i)*dy - yc)*(yl + (0.5+i)*dy - yc);
    zsq[i] = (zl + (0.5 + i)*dz - zc)*(zl + (0.5+i)*dz - zc);
  }

  int inside_count = 0;

  for (int k = 0; k < nsample; k++){
    for (int j = 0; j < nsample; j++){
      for (int i = 0; i < nsample; i++){
       float r = sqrt( factor * xsq[i] + ysq[j] + zsq[k] );

       if ( r <= radius) inside_count++;
      }
    }
  }


  return ((float) inside_count) / ((float) nsample*nsample*nsample);
}

void IndividualStarSetStellarWindProperties(Star *cstar, const float &Time,
                                            const float &dtFixed, const float &TimeUnits,
                                            float &m_eject,
                                            float &E_thermal, float *metal_mass){

  float wind_lifetime, agb_start_time, wind_dt;
  const float msun = 1.989E33, k_boltz = 1.380658E-16;

  /* New variables to make code slightly cleaner + handle units */
  float mproj        = cstar->ReturnBirthMass();
  float lifetime     = cstar->ReturnLifetime() * TimeUnits;
  float metallicity  = cstar->ReturnMetallicity();
  float particle_age = (Time - cstar->ReturnBirthTime())*TimeUnits;
  float dt           = dtFixed * TimeUnits;

  int *yield_table_position = cstar->ReturnYieldTablePosition();
  int *se_table_position    = cstar->ReturnSETablePosition();

  float m_eject_total = 0.0;

  if( IndividualStarFollowStellarYields && ((TestProblemData.MultiMetals == 2) || (MultiMetals == 2))){

    // 1 = wind, -1 = return total mass
    m_eject = StellarYieldsInterpolateYield(1, yield_table_position[0], yield_table_position[1],
                                            mproj, metallicity, -1); // total ejecta mass in Msun
    m_eject_total = m_eject;

    wind_lifetime = lifetime;   // CGS units

    if( mproj < IndividualStarAGBThreshold) {

      // 2 at end of argument implies compute start time of AGB phase
      IndividualStarInterpolateLifetime(agb_start_time, se_table_position[0], se_table_position[1],
                                           mproj, metallicity, 2);
      /* sanity check */
      float temp_lifetime;
      IndividualStarInterpolateLifetime(temp_lifetime, se_table_position[0], se_table_position[1], mproj, metallicity, 1);

      wind_lifetime = lifetime - agb_start_time; // CGS Units
      if (wind_lifetime < 0.0){
        printf("WARNING LIFETIME ISSUE --- lifetime = %"ESYM" agb_start = %"ESYM" temp_lifetime = %"ESYM"\n", lifetime, agb_start_time, temp_lifetime);
      }

        //
        // To ensure total (integrated) mass ejected is accurate, make sure we don't overinject
        // mass when winds should only be "ON" for part of a timestep, either at beginning or end
        // of AGB phase, or when AGB phase is unresolved (i.e. AGB time < dt)
        //
/*
        if (particle_age > lifetime && particle_age - dt < lifetime){

          // wind_dt = fmin( fmax(0.0, lifetime - (particle_age - dt)) , lifetime - agb_start_time);
          wind_dt = fmax(0.0, lifetime - (particle_age - dt));
          wind_dt = fmin( wind_dt, lifetime - agb_start_time);

          // printf("wind lifetime mode 1\n");
        } else if (particle_age > agb_start_time && particle_age < lifetime ) {
          wind_dt = fmin(particle_age - agb_start_time,dt); // wind only occurs for part of timestep + star dies

          // printf("wind lifetime mode 2\n");
        } else if (particle_age < agb_start_time && particle_age + dt > lifetime) {
          //
          // AGB phase is unresolved. Set wind timestep to lifetime to do all ejecta this timestep
          //
          wind_dt = wind_lifetime;
          // printf("wind lifetime mode 3\n");
        } else if (particle_age < agb_start_time && particle_age + dt > agb_start_time){
          wind_dt = particle_age + dt - agb_start_time; // wind only occurs for part of timestep
          // printf("wind lifeitme mode 4\n");
        } else{
          wind_dt = fmin( lifetime - agb_start_time, dt);
          // printf("PROBLEM IN AGB WIND PHASE\n");
        }
*/

      // printf("Wind lifetime = %"ESYM" - wind_dt = %"ESYM"  %"ESYM" %"ESYM" %"ESYM" %"ESYM"\n",wind_lifetime, wind_dt, lifetime, agb_start_time, particle_age, dt);


      // end AGB check
    } else { // massive stars (constant wind)

        //
        // Check timestep to make sure we don't overinject yields at end of life
        //
        wind_dt = fmin(fmax(lifetime - particle_age, 0.0) , dt);
        if (wind_dt < 0.0){
           wind_dt = dt - (particle_age - lifetime);

           if(abs(wind_dt) > dt){
             printf("DEBUG WARNING: Something very wrong is happending at stellar wind end of life\n");
             wind_dt = 0.001*dt;
           }
        }
    }

    /* Gaurd against cases where agb phase is zero */
    wind_lifetime = (wind_lifetime < tiny_number) ? dt : wind_lifetime;
//    wind_dt       = (wind_dt       < tiny_number) ? dt : wind_dt;
    wind_dt = dt;
    //printf("corrected Wind lifetime = %"ESYM" - wind_dt = %"ESYM"  %"ESYM" %"ESYM" %"ESYM" %"ESYM"\n",wind_lifetime, wind_dt, lifetime, agb_start_time, particle_age, dt);

    if (dt == 0){
        m_eject = 0.0;
        // printf("WARNING: ZERO TIME STEP SIZE IN WIND LAUNCHING");
    } else{
        m_eject  /= wind_lifetime ; // average mass loss rate over entire wind lifetime
    }

    // end yields methods
  } else {

    // use model to compute mass loss rate instead

    ComputeStellarWindMassLossRate(mproj, metallicity, &m_eject);

    wind_dt = fmin( fmax(lifetime - particle_age, 0.0), dt);
    if (wind_dt < 0.0){
       wind_dt = dt - (particle_age - lifetime);

       if(abs(wind_dt) > dt){
         printf("DEBUG WARNING: Something very wrong is happending at stellar wind end of life\n");
         wind_dt = 0.001*dt;
       }
    }

    wind_dt = dt;
    m_eject_total = m_eject;

  } // end  checking for yields


  float v_wind;

  if(mproj < IndividualStarAGBThreshold){
    /* no good model for AGB wind - use constant user velocity */
    v_wind = IndividualStarAGBWindVelocity;

  } else  if (IndividualStarStellarWindVelocity < 0){
    ComputeStellarWindVelocity(cstar, &v_wind); // v in km /s

  } else{
    v_wind = IndividualStarStellarWindVelocity; // user chosen, in km/s
  }

  if (v_wind > IndividualStarMaximumStellarWindVelocity &&
      IndividualStarMaximumStellarWindVelocity > 0)
            v_wind = IndividualStarMaximumStellarWindVelocity;

  v_wind *= 1.0E5; // now in cgs

  /* Now that we have wind lifetime and ejected mass, compute properties of wind*/

  float correction_factor = 1.0;
  m_eject = m_eject * wind_dt; // convert Mdot to M_ej

  /* If mass injection exceeds total, correct this fractionally -
     this should rarely happen, but is here to ensure correct chemical evolution.
     Also, if this is likely the last time step the particle is alive, dump
     everything. Again, this is not physically the best solution, but
     ensures that all of the correct yields get deposited, but sacrifices
     correct temporal injection. Given that dt is generally small, this means
     the time at which yields get injected may only be off by 10^3 - 10^4 years...
     this should be irrelevant for galaxy scale (100-1000 Myr) simulations. */
  if( (cstar->ReturnWindMassEjected() + m_eject > m_eject_total) ||
      (particle_age + 2.5*dt > lifetime) ){

    float old_ejection = m_eject;
    m_eject = fmax(m_eject_total - cstar->ReturnWindMassEjected(), 0.0);
    correction_factor = m_eject / old_ejection;
  }

  float Teff, R; // Need Teff for computing thermal energy of ejecta
  IndividualStarInterpolateProperties(Teff, R,
                                      se_table_position[0], se_table_position[1],
                                      cstar->ReturnBirthMass(), cstar->ReturnMetallicity());

  E_thermal = 1.5 * Teff * (m_eject*msun / (1.0 * 1.67E-24)) * k_boltz; // current T of wind

  if( v_wind > IndividualStarMaximumStellarWindVelocity * 1.0E5){ // so we don't waste CPU
    v_wind = IndividualStarMaximumStellarWindVelocity * 1.0E5;
  }

  E_thermal = E_thermal + 0.5 * (m_eject * msun) * v_wind * v_wind; // assume 100% KE thermalization

  /* finally, compute metal masses if needed */
  float wind_scaling = wind_dt / wind_lifetime  * correction_factor;

//  if (wind_lifetime <= tiny_number){
//    wind_scaling = 0.0;
//  }

  if(IndividualStarFollowStellarYields && TestProblemData.MultiMetals==2){

    metal_mass[0] = StellarYieldsInterpolateYield(1, yield_table_position[0], yield_table_position[1],
                                                     mproj, metallicity, 0); // total metal in Msun

    for (int i = 0; i < StellarYieldsNumberOfSpecies; i++){
      metal_mass[1 + i] = StellarYieldsInterpolateYield(1, yield_table_position[0], yield_table_position[1],
                                                        mproj, metallicity,
                                                        StellarYieldsAtomicNumbers[i]);
    }

    /* scale for wind dt and lifetime */
    for(int i = 0; i < StellarYieldsNumberOfSpecies+1; i ++){
      metal_mass[i] *= wind_scaling;
    }

  }


  for(int i = 0; i < StellarYieldsNumberOfSpecies+1; i++){
      if(metal_mass[i] < 0.0){
        printf("particle age = %"ESYM" lifetim - age = %"ESYM" dt %"ESYM" %"ESYM"\n", particle_age, lifetime-particle_age, dt, wind_scaling);
        printf("metal mass = %"ESYM" wind_dt = %"ESYM" wind_lifetime = %"ESYM" eject = %"ESYM"\n",metal_mass[i], wind_dt, wind_lifetime, m_eject);
        if(i>0){
            printf("i = %"ISYM" anum = %"ISYM"\n", i, StellarYieldsAtomicNumbers[i-1]);
        } else{
            printf("i = %"ISYM"\n",i);
        }
        cstar->PrintInfo();

        ENZO_FAIL("Negative metal mass in wind setup");
      }
  }


  // done computing stellar wind properties
  return;
}
