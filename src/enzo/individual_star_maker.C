
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
int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);

int FindField(int f, int farray[], int n);

float SampleIMF(void);

float ComputeSnIaProbability(const float &current_time, const float &formation_time, const float &lifetime, const float &TimeUnits);
unsigned_long_int mt_random(void);

void ComputeStellarWindVelocity(Star *cstar, float *v_wind);

void ComputeStellarWindMassLossRate(const float &mproj, const float &metallicity,
                                    float *dMdt);

int SetWDLifetime(float &WD_lifetime,
                    const float &current_time, const float &formation_time,
                    const float &lifetime, const float &TimeUnits);



void SetFeedbackCellCenter(const FLOAT &xp, const FLOAT &yp, const FLOAT &zp,
                             const FLOAT &xstart, const FLOAT &ystart, const FLOAT &zstart,
                             const FLOAT &dx,
                             const int &nx, const int &ny, const int &nz, const int &ibuff,
                             FLOAT *xfc, FLOAT *yfc, FLOAT *zfc);

void Momentum(float *u, float *v, float *w, float *d,
              const float &up, const float &vp, const float &wp,
              const int &nx, const int &ny, const int &nz,
              const int &ic, const int &jc, const int &kc,
              const int &iface, const int &jface, const int &kface,
              const int &stencil, int idir);

void MetalConversion(float *metal_field, float *d, const float &dx,
                     const int &nx, const int &ny, const int &nz,
                     const int &ic, const int &jc, const int &kc,
                     const int &stencil, int idir);

void SumMassEnergy(float *pu, float *pv, float *pw, float *d, float *ge, float *te,
                   const int &nx, const int &ny, const int &nz,
                   const int &iface, const int &jface, const int &kface,
                   const int &ic, const int &jc, const int &kc, const int &stencil,
                   float *mass_sum, float *energy_sum, float *kin_energy_sum);

void ComputeAbcCoefficients(float *pu, float *pv, float *pw, float *d,
                            float *ge, float *pu_l, float *pv_l, float *pw_l,
                            float *d_l, const int &nx, const int &ny, const int &nz,
                            const int &iface, const int &jface, const int &kface,
                            const int &ic, const int &jc, const int &kc,
                            const int &stencil,
                            float &A, float &B, float &C);

void AddFeedbackToGridCells(float *pu, float *pv, float *pw, 
                            float *d, float *ge, float *te,
                            const int &nx, const int &ny, const int &nz,
                            const int &ic, const int &jc, const int &kc,
                            const int &iface, const int &jface, const int &kface,
                            const float &dxf, const float &dyf, const float &dzf,
                            const float &dxc, const float &dyc, const float &dzc,
                            const float &mass_per_cell, const float &mom_per_cell,
                            const float &therm_per_cell, const int &stencil,
                            const float additional_mass_factor);

void AddMetalSpeciesToGridCells(float *m, const float &mass_per_cell,
                                const int &nx, const int &ny, const int &nz,
                                const int &ic, const int &jc, const int &kc,
                                const float &dxc, const float &dyc, const float &dzc,
                                const int &stencil, const float additional_mass_factor);

void IndividualStarSetCoreCollapseSupernovaProperties(Star *cstar,
                                                      float &m_eject, float &E_thermal, float *metal_mass);

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

    xpos[i] = xpos[i]*pc/LengthUnits + DiskGravityPosition[0];
    ypos[i] = ypos[i]*pc/LengthUnits + DiskGravityPosition[1];
    zpos[i] = zpos[i]*pc/LengthUnits + DiskGravityPosition[2];

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



    /* now assign metal abundnace fractions as all tiny numbers  */
    if (TestProblemData.MultiMetals == 2){
      for (int is = 0; is < StellarYieldsNumberOfSpecies; is++){
        ParticleAttribute[4 + is][count] = tiny_number;
      }
    }

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
      if(TestProblemData.MultiMetals == 2){
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
    /* Metal fields are all in fractions, as set in Grid_StarParticleHandler */
    if(TestProblemData.MultiMetals == 2){
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
                                float *ParticleVelocity[], float *ParticleAttribute[]){
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
  float bmass, div, star_mass=0.0, sum_mass=0.0;
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

          /* if distributed star formation */
          // check center cell's SF condition, if met, do SF
          /* loop and sum over all*/

           bmass = (BaryonField[DensNum][index]*(dx*dx*dx)) * MassUnits / msolar; // in solar masses
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
            isosndsp2 = sndspdC * temp[index] ;
            jeansmass = pi / (6.0 * sqrt(BaryonField[DensNum][index]*DensityUnits) *
                            POW(pi * isosndsp2 / GravConst ,1.5)) / msolar; // in solar masses

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


              float M_max_star;
              if (IndividualStarAllowTruncatedIMF){
                  /* allow IMF to be truncated to safely allow SF in less dense regions. This may be an
                     unphysical thing to do, as it won't necessarily reproduce the input IMF. In reality
                     SF clouds will do this (low mass clouds are less likely to make very massive stars), but
                     this is all convolved together in the input IMF (IMF may vary cloud-by-cloud)........
                     use this with caution */
                  M_max_star =  min(bmass * IndividualStarMassFraction, IndividualStarIMFUpperMassCutoff);
              } else {
                  /* only allow star formation if IMF can fully sampled (safely) in the given region */
                  if ( bmass * IndividualStarMassFraction < IndividualStarIMFUpperMassCutoff){
                    break;
                  }
                  M_max_star = IndividualStarIMFUpperMassCutoff;
              }


              if(IndividualStarSFAlgorithm == 0){ /* DO NOT USE THIS */

                // calculate mass in cell that can be converted to stars in timestep
                // generally this should be small (comparable to or less than the lower mass
                // cutoff of the IMF)

                star_fraction  = min(StarMakerMassEfficiency*(this->dtFixed)/tdyn, 1.0);
                mass_to_stars  = star_fraction * bmass;
                mass_available = IndividualStarMassFraction * bmass; // AJE - fixed 6/4
                mass_to_stars  = min(mass_to_stars, mass_available);

                // If mass_to_stars greater than available mass, convert
                // all of available mass into stars
                // Frankly this is very unlikely to occur...
                // Tests as of 2/22/16 show NO SF here for at least 10^5 stars in a LMC dwarf galaxy
                if(mass_to_stars >= mass_available){
                  mass_to_stars = mass_available;
                  while( ii < *nmax && mass_to_stars > M_max_star){
                    ParticleMass[ii] = SampleIMF();
                    sum_mass        += ParticleMass[ii]; // counter for mass formed in this cell
                    mass_to_stars   -= ParticleMass[ii]; // reduce available mass
                    ii++;
                  }
                }

                // Tests (as of 2/22/16) show NO SF here for at least the first 10^5 stars
                if (mass_to_stars > M_max_star){
                  while (ii < *nmax && mass_to_stars > M_max_star){
                    ParticleMass[ii]  = SampleIMF();
                    sum_mass         += ParticleMass[ii];
                    mass_to_stars    -= ParticleMass[ii];
                    ii++;
                  }
                }

                // If mass is above IMF lower limit, star formation will happen.
                // Just form stars randomly over IMF until mass dips below lower cutoff
                if(mass_to_stars > IndividualStarIMFLowerMassCutoff){

                  // loop until mass to stars is less than 10% of smallest star particle size
                  while( ii < *nmax && mass_to_stars > 1.1 * IndividualStarIMFLowerMassCutoff){
                    float tempmass;
                    tempmass = SampleIMF();

                    if (tempmass < M_max_star){
                        ParticleMass[ii]  = SampleIMF();
                        sum_mass         += ParticleMass[ii];
                        mass_to_stars    -= ParticleMass[ii];
                        ii++;
                    } // else redraw
  
                    if (mass_to_stars < 0.0){
                      mass_to_stars = 0.0;
                    }
                  }
                } // end mass above individual star cutoff

                // now we are in the Goldbaum et. al. 2015 regime (star_maker_ssn.F)
                // Calculate probability of star forming and form stars stochastically
                if (mass_to_stars < IndividualStarIMFLowerMassCutoff && mass_to_stars > tiny_number){
                  star_mass = SampleIMF();
                  pstar     = mass_to_stars / star_mass;
                  rnum =  (float) (random() % max_random) / ((float) max_random);
                  if (rnum < pstar){
                    ParticleMass[ii]  = star_mass;
                    sum_mass         += ParticleMass[ii];
                    ii++;
                  }
                }

              } else if (IndividualStarSFAlgorithm == 1){
                /* take chunks of mass */

                if( bmass*IndividualStarMassFraction > IndividualStarSFGasMassThreshold ){ // set to ~ 2 x M_max_star
                  // if true, we can try and form stars. compute probability that this mass will
                  // form stars this timestep
                  star_fraction  = min(StarMakerMassEfficiency*(this->dtFixed)/tdyn, 1.0);
                  mass_to_stars  = star_fraction * bmass;

                  pstar          = mass_to_stars / IndividualStarSFGasMassThreshold;

                  rnum           = (float) (random() % max_random) / ( (float) max_random);

                  if ( rnum < pstar){ // form stars until mass runs out - keep star if too much is made
                      float mass_counter = IndividualStarSFGasMassThreshold;
                      while( mass_counter > 0.0){
                          ParticleMass[ii]  = SampleIMF();
                          sum_mass         += ParticleMass[ii];
                          mass_counter     -= ParticleMass[ii];
                          ii++;
                      }

                  } // endif randum number draw check

                } // endif mass threshold check

              } // endif sf algorithm for star formation

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

                ParticleType[istar]            = - PARTICLE_TYPE_INDIVIDUAL_STAR;   // negative is a "new" star
                ParticleAttribute[0][istar]    = this->Time;                        // formation time
                ParticleAttribute[2][istar]    = BaryonField[MetalNum][index]; // metal fraction (conv from density in Grid_StarParticleHandler)


                if(IndividualStarInterpolateLifetime(ParticleAttribute[1][istar], ParticleMass[istar],
                                                                                  ParticleAttribute[2][istar], 1) == FAIL){
                  ENZO_FAIL("Error in stellar lifetime interpolation");
                }
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
                if(TestProblemData.MultiMetals == 2){
                  for( int iyield = 0; iyield < StellarYieldsNumberOfSpecies; iyield++){
                    if(StellarYieldsAtomicNumbers[iyield] > 2){
                      int field_num;

                      this->IdentifyChemicalTracerSpeciesFieldsByNumber(field_num, StellarYieldsAtomicNumbers[iyield]);

                      ParticleAttribute[4 + iyield][istar] = BaryonField[field_num][index];

                    } else if (StellarYieldsAtomicNumbers[iyield] == 1){
                      // Take H and He fractions as TOTAL amount of H in the cell
                      // this is probably not needed since it should all be HI in a star forming region anyway
                      ParticleAttribute[4 + iyield][istar] = BaryonField[HINum][index] + BaryonField[HIINum][index];

                      if (MultiSpecies > 1){
                        ParticleAttribute[4 + iyield][istar] += BaryonField[HMNum][index] +
                                                         BaryonField[H2INum][index] + BaryonField[H2IINum][index];
                      }
                    } else if (StellarYieldsAtomicNumbers[iyield] == 2){
                      // Again, total amount of Helium - probably not necessary, should all be HeI anyway
                      ParticleAttribute[4 + iyield][istar] = BaryonField[HeINum][index]  +
                                           BaryonField[HeIINum][index] + BaryonField[HeIIINum][index];

                    }
                  }
                }

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

                ParticleAttribute[NumberOfParticleAttributes-2][istar] = 0.0;
                ParticleAttribute[NumberOfParticleAttributes-1][istar] = 0.0;
/*
                printf(" Mass = %"FSYM" Z = %"FSYM" ",ParticleAttribute[3][istar], ParticleAttribute[2][istar]);
                for( int ti = tstart; ti < tstart + 7; ti++){
                   printf("%"FSYM" ", ParticleAttribute[ti][istar]);
                }
                printf("\n");
*/


              } // end while loop for assigning particle properties
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

float SampleIMF(void){
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
  float dm = log10(IndividualStarIMFUpperMassCutoff / IndividualStarIMFLowerMassCutoff)/ ((float) (IMF_TABLE_ENTRIES-1));
  float m;

  int bin_number;

  if (x <= IMFData[0] ){
    bin_number = 0;
  } else if (x >= IMFData[IMF_TABLE_ENTRIES-1]){
    bin_number = IMF_TABLE_ENTRIES - 1;
  } else{
    bin_number = search_lower_bound(IMFData, x, 0, IMF_TABLE_ENTRIES, IMF_TABLE_ENTRIES);
  }

  m = IndividualStarIMFLowerMassCutoff * POW(10.0, bin_number * dm);

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


/*

 ------------------- Begin old method that should be excised ----------------------
                          --- comment date Oct 2016 ---

int grid::individual_star_feedback(int *np,
                                   float *ParticleMass, int *ParticleType, FLOAT *ParticlePosition[],
                                   float *ParticleVelocity[], float *ParticleAttribute[]){
-----------------------------------------------------------------------------
  Handles the feedback for the indivual stars formed. This includes mechanical
  feedback from stellar winds, supernovae, and (if enabled) chemical yield 
  deposition.

  INPUTS
    nx, ny, nz   - size of grid in each dimension
    dx           - current grid size (code units)
    dt           - current timestep  (code units)
    current_time - time (code units)
    DensityUnits,LengthUnits,VelocityUnits,TimeUnits  - conversion between code units and cgs
    x/y/z start  - start position of grid in each dimension
    ibuff        - size of ghost zones in each dimension
    np           - number of particles to loop over
    ParticleMass -
    ParticlePosition -
    ParticleVelocity -
    ParticleAttribute -
-----------------------------------------------------------------------------
//   copy some grid parameters for convenience


  int  nx = this->GridDimension[0];
  int  ny = this->GridDimension[1];
  int  nz = this->GridDimension[2];
  int  ibuff = NumberOfGhostZones;

  FLOAT xstart = CellLeftEdge[0][0], ystart = CellLeftEdge[1][0], zstart = CellLeftEdge[2][0];
  float   dx   = CellWidth[0][0];



  /* Get Units 
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits, MassUnits;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, this->Time) == FAIL){
      ENZO_FAIL("Error in GetUnits");
  }
  MassUnits   = DensityUnits*LengthUnits*LengthUnits*LengthUnits; // mass units

  float mp, distmass, energy, dratio;
  const double msolar = 1.989e33;                 // solar mass in cgs
  const double speed_of_light = 2.99792458e10 ;

  FLOAT particle_age, lifetime;

  const int max_random = (1<<16);

  int ip, jp, kp, index; // particle location

  int do_stellar_winds, go_supernova;

  const float fudge_factor = 1.0E10; // see comment below - hack to get WD and SNIa to work

  int IndividualStar        = PARTICLE_TYPE_INDIVIDUAL_STAR,
      IndividualStarWD      = PARTICLE_TYPE_INDIVIDUAL_STAR_WD,
      IndividualStarRemnant = PARTICLE_TYPE_INDIVIDUAL_STAR_REMNANT; // convenience


  // loop over all star particles
  for(int i = 0; i < (*np); i++){

    // where is the particle?
    ip = int ( (ParticlePosition[0][i] - (xstart)) / dx);
    jp = int ( (ParticlePosition[1][i] - (ystart)) / dx);
    kp = int ( (ParticlePosition[2][i] - (zstart)) / dx);

    float birth_mass = ParticleAttribute[3][i];

    mp = ParticleMass[i] * dx*dx*dx; // mass in code units
    mp = mp * (MassUnits) / msolar ;        // Msun

    // warning if outside current grid
    if( ip < 0 || ip > nx || jp < 0 || jp > ny || kp < 0 || kp > nz){
      printf("Warning: star particle is outside of grid\n");
      printf(" %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" %"ISYM" \n", ip, jp, kp, nx, ny, nz);
      printf(" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM"\n", dx, xstart, ystart, zstart, ParticlePosition[0][i], ParticlePosition[1][i], ParticlePosition[2][i]);
      return FAIL;
    }


    particle_age = (this->Time) - ParticleAttribute[0][i];
    lifetime     = ParticleAttribute[1][i];

    do_stellar_winds           = FALSE;
    go_supernova               = FALSE;

    if(ParticleType[i] == IndividualStar){
      if(IndividualStarStellarWinds){

        float wind_start_age = 0.0;
        if(birth_mass < IndividualStarAGBThreshold){
          // star has AGB phase, do winds only in last fraction of lifetime
          if(IndividualStarInterpolateLifetime(wind_start_age, birth_mass,
                                               ParticleAttribute[2][i], 2) == FAIL){
            ENZO_FAIL("IndividualStarFeedback: Failure in main sequence lifetime interpolation");
          }
        }

        // dt addition here so stars w/ unresolved agb phases dump all winds at last timestep
        if (particle_age  < lifetime && particle_age + this->dtFixed > wind_start_age){
          do_stellar_winds = TRUE;
        }

      } // end winds check

      // move all attribute changes until AFTER feedback has occured!!!!!
      if ( birth_mass >= IndividualStarSNIIMassCutoff && ((particle_age + this->dtFixed) > lifetime)){
        go_supernova = TRUE;
        ParticleAttribute[1][i] = lifetime*1.0E10; // lifetime set to arbitrarily large number

      } else if ( birth_mass >= IndividualStarWDMinimumMass && birth_mass <= IndividualStarWDMaximumMass
                  && particle_age + (this->dtFixed) > lifetime && IndividualStarUseSNIa){

        float wd_mass; // compute WD mass using linear fit from Salaris et. al. 2009
        if(      birth_mass  < 4.0){ wd_mass = 0.134 * birth_mass + 0.331;}
        else if (birth_mass >= 4.0){ wd_mass = 0.047 * birth_mass + 0.679;}

        ParticleMass[i] = wd_mass * (msolar / MassUnits) / (dx*dx*dx);
        ParticleType[i] = IndividualStarWD;
        /* fudge factor makes lifetime very long so particle is not deleted,
           however, SNIa scheme needs to know the main sequence lifetime. This
           confusing, but a little bit more efficient than making a new particle attribute 
        ParticleAttribute[1][i] = lifetime*fudge_factor; // make a big number

      } else if (particle_age + (this->dtFixed) > lifetime){
        ParticleMass[i] = 0.0; // kill silently
      }


    } else if (ParticleType[i] == IndividualStarWD){
      /* White Dwarf Feedback - Make SNIa 

      /* Does the progenitor mass (main sequence mass) fit within range 
      if( (birth_mass > IndividualStarSNIaMinimumMass) &&
          (birth_mass < IndividualStarSNIaMaximumMass) ){

        float formation_time = ParticleAttribute[0][i]; // formation time of main sequence star

        float PSNIa;
        float rnum;

        /* Probability that the star will explode as SNIa in this timestep 
        PSNIa  = ComputeSnIaProbability( this->Time, formation_time, lifetime/fudge_factor, TimeUnits); // units of t^((beta)) / s
        PSNIa *= this->dtFixed;

        rnum =  (float) (random() % max_random) / ((float) (max_random));

//        printf("individual_star_feedback: SNIa - M_proj, PSNIa, rnum = %"ESYM" %"ESYM"\n", PSNIa, rnum);

        if (rnum < PSNIa){
          go_supernova = TRUE;
        }
      } // end SNIa progenitor + WD check

    } // end type check

    if(do_stellar_winds || go_supernova){
        float sum_dens = 0.0;

        // call feedback function to do stellar winds
        if(do_stellar_winds){


//          printf("ISF: Calling feedback general to do stellar winds\n");
//          printf("ISF: Current Mass = %"ESYM" Particle aatribute 3 = %"ESYM" mproj = %"ESYM"\n", mp,ParticleAttribute[3][i], birth_mass);

          /* Apply stellar wind feedback. Determined by setting last arguemtn to -1 
          this->IndividualStarAddFeedbackGeneral(ParticlePosition[0][i], ParticlePosition[1][i], ParticlePosition[2][i],
                                                 ParticleVelocity[0][i], ParticleVelocity[1][i], ParticleVelocity[2][i],
                                                 birth_mass, ParticleAttribute[1][i], particle_age, ParticleAttribute[2][i], &mp, -1);

          ParticleMass[i] = mp* (msolar/MassUnits) / (dx*dx*dx); // update particle mass and put in code units
        }

        // call feedback function to do supernova feedback (either SNIa or core collapse)
        if(go_supernova){

          if( ParticleType[i] != IndividualStarWD){
//            printf("Calling feedback to do cc supernova");
            /* do core collapse supernova feedback - set by last value == 1 
            this->IndividualStarAddFeedbackGeneral(ParticlePosition[0][i], ParticlePosition[1][i], ParticlePosition[2][i],
                                                   ParticleVelocity[0][i], ParticleVelocity[1][i], ParticleVelocity[2][i],
                                                   birth_mass, ParticleAttribute[1][i], particle_age, ParticleAttribute[2][i], &mp, 1);

            ParticleMass[i] = mp * (msolar/MassUnits) / (dx*dx*dx); // update particle mass and put in code units
            ParticleType[i] = IndividualStarRemnant; // change type
            ParticleAttribute[1][i] = 1.0E10 * ParticleAttribute[1][i];
          } else{
            printf("calling feedback to do supernova 1a\n");
            /* do SNIa supernova feedback - set by last value == 1 
            this->IndividualStarAddFeedbackGeneral(ParticlePosition[0][i], ParticlePosition[1][i], ParticlePosition[2][i],
                                             ParticleVelocity[0][i], ParticleVelocity[1][i], ParticleVelocity[2][i],
                                             birth_mass, ParticleAttribute[1][i], particle_age, ParticleAttribute[2][i], &mp, 2);

            ParticleMass[i]     = 0.0;                             // make particle mass zero - now a masless tracer
            ParticleAttribute[1][i] = 1.0E10 * ParticleAttribute[1][i];    // make lifetime infinite  -

          }






          // put attriute changes here
        }

    } // if do feedback


  } // loop over particles

  return SUCCESS;
}

 ------------------- End old method that should be excised ----------------------
                          --- comment date Oct 2016 ---

*/


/*

 ------------------- Start old method that should be excised ----------------------
                          --- comment date Oct 2016 ---


float ComputeSnIaProbability(const float &current_time, const float &formation_time,
                             const float &lifetime, const float &TimeUnits){
 /* -----------------------------------------------------------------
  * ComputeSnIaProbability
  *
  * Computes dPdt for a given white dwarf that might go supernova. The
  * probability is such that the integral over dP from the time the WD
  * was born for a hubble time afterwards is a free parameter on order of
  * a few percent. This is IndividualStarSNIaFraction, or fraction of WD's
  * over a certain progenitor mass range that will go supernova in a hubble
  * time.
  * ------------------------------------------------------------------- 


 float dPdt;
 const float hubble_time = 4.382E17; // need to allow for cosmology units AJE TO DO
                                     // currently assumes H_o = 70.4

 dPdt = IndividualStarSNIaFraction;

 /* conmpute normalized probability - normalized by integral over WD formation time to hubble time 
 if (IndividualStarDTDSlope == 1.0){
   dPdt /= log( ((hubble_time / TimeUnits) + lifetime) / lifetime );
 } else{
   dPdt *= (-IndividualStarDTDSlope + 1.0);
   dPdt /= ( POW( (hubble_time / TimeUnits) + lifetime   , -IndividualStarDTDSlope + 1) -
             POW( (lifetime)                      , -IndividualStarDTDSlope + 1));
 }

 dPdt = dPdt * POW( ((current_time) - (formation_time)), -IndividualStarDTDSlope);


 return dPdt;
}

 ------------------- End old method that should be excised ----------------------
                          --- comment date Oct 2016 ---


*/

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

/*
--------------------- Being Old method that needs to be excised -----------------------------

int grid::IndividualStarAddFeedbackGeneral(const FLOAT &xp, const FLOAT &yp, const FLOAT &zp,
                       const float &up, const float &vp, const float &wp,
                       const float &mproj, const float &lifetime, const float &particle_age,
                       const float &metallicity, float *mp, int mode     ){

// mode < 0 - wind .... mode = 1 SNII ... mode = 2 SNIA

  float m_eject, E_thermal, E_kin, f_kinetic, v_wind, p_feedback;
  const double c_light = 2.99792458E10; const double msolar = 1.989E33;

  float *metal_mass;

  /* Get Units 
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits, MassUnits;
  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, this->Time) == FAIL){
      ENZO_FAIL("Error in GetUnits");
  }
  MassUnits   = DensityUnits*LengthUnits*LengthUnits*LengthUnits; // mass units


  /* rename some grid parameters for convenience 
  int  nx = this->GridDimension[0], ny = this->GridDimension[1], nz = this->GridDimension[2];
  int  ibuff = NumberOfGhostZones;

  FLOAT xstart = CellLeftEdge[0][0], ystart = CellLeftEdge[1][0], zstart = CellLeftEdge[2][0];
  float   dx   = CellWidth[0][0];


  //printf("ISF: Assiging initital values to metal mass\n");
  if(IndividualStarFollowStellarYields && TestProblemData.MultiMetals == 2){

    metal_mass = new float[StellarYieldsNumberOfSpecies + 1];

    for (int i = 0; i < StellarYieldsNumberOfSpecies + 1; i ++){
      metal_mass[i] = 0.0;
    }

  } else { metal_mass = NULL; }


  float wind_dt = 0.0, wind_lifetime = lifetime;
  int check_mass_in_region = 0;
  float E_thermal_wind = 0.0;
  if (mode < 0){ // computing stellar winds

    /* -------------------------------------------
     *
     * Compute the mass ejecta and energetics for
     * stellar winds
     *
     * -------------------------------------------

    /* compute ejecta mass - use yield tables if yields are tracked, otherwise use model 
    if ( IndividualStarFollowStellarYields && TestProblemData.MultiMetals == 2){
      m_eject   = StellarYieldsInterpolateYield(1, mproj, metallicity, -1) *msolar/ MassUnits; /* first arg, 1 = wind ; last -1 = tot mass 
      //printf("TOTAL INTERP: m_eject msolar MassUnits  %"ESYM" %"ESYM" %"ESYM"\n", m_eject*MassUnits/msolar, msolar, MassUnits);
      wind_lifetime = lifetime;
      if (mproj < IndividualStarAGBThreshold){

        float agb_start_time ; // point in star's life where AGB phase occurs

        if(IndividualStarInterpolateLifetime(agb_start_time, mproj, metallicity, 2) == FAIL){
          ENZO_FAIL("Individual star add feedback general: failure in interpolating lifetime");
        }

        agb_start_time /= TimeUnits;
        wind_lifetime   = lifetime - agb_start_time; // wind lifetime is particle life - start time of AGB

        //
        // To ensure total (integrated) mass ejected is accurate, make sure we don't overinject
        // mass when winds should only be "ON" for part of a timestep, either at beginning or end
        // of AGB phase, or when AGB phase is unresolved (i.e. AGB time < dt)
        //
        if (particle_age > agb_start_time && particle_age + this->dtFixed > lifetime) {
          wind_dt = lifetime - particle_age; // wind only occurs for part of timestep + star dies

        } else if (particle_age < agb_start_time && particle_age + this->dtFixed > lifetime) {
          //
          // AGB phase is unresolved. Set wind timestep to lifetime to do all ejecta this timestep
          //
          wind_dt = wind_lifetime; 

        } else if (particle_age < agb_start_time && particle_age + this->dtFixed > agb_start_time){
          wind_dt = particle_age + this->dtFixed - agb_start_time; // wind only occurs for part of timestep
        }

      } else { // massive stars (constant winds over lifetime)

        //
        // Check timestep to make sure we don't overinject yields at end of life
        //
        wind_dt = fmin(lifetime - particle_age, this->dtFixed);
        if (wind_dt < 0.0){
           wind_dt = this->dtFixed - (particle_age - lifetime);

           if(abs(wind_dt) > this->dtFixed){
             printf("DEBUG WARNING: Something very wrong is happending at stellar wind end of life\n");
             wind_dt = 0.001*this->dtFixed;
           }
        }
      }

      m_eject  /= wind_lifetime ; // average mass loss rate over entire wind lifetime

    } else{
      // gives m_eject as Msun / s
      ComputeStellarWindMassLossRate(mproj, metallicity, &m_eject);
      m_eject = m_eject * msolar / MassUnits * TimeUnits;  // convert to code mass / code time

      // make sure we don't over-inject mass (i.e. partial timestep)
      wind_dt = fmin( lifetime - particle_age, this->dtFixed);
    }

    // Finally, compute total amount of mass ejected this timestep
    m_eject = m_eject * wind_dt;

    E_thermal = 1.5 * 2.0E5 * (m_eject / (1.0*1.67E-24/MassUnits)) * (1.380658E-16) / ( MassUnits*VelocityUnits*VelocityUnits);
    E_thermal_wind = E_thermal; // save this separate component

    /* set wind velocity depending on mode
    if ( IndividualStarStellarWindVelocity < 0){
      ComputeStellarWindVelocity(mproj, metallicity, &v_wind); // compute wind velocity in km / s using model
    } else {
      v_wind = IndividualStarStellarWindVelocity; // wind velocity in km / s
    }

//    printf("ISF: Stellar wind mass in Msun = %"ESYM" in code units %"ESYM"\n", m_eject *MassUnits/msolar, m_eject);
    //printf("ISF: Total Expected momentum in cgs %"ESYM" and in code units %"ESYM"\n", v_wind*1.0E5*m_eject*MassUnits/msolar, m_eject * v_wind*1.0E5/VelocityUnits);
//    printf("ISF: Stellar wind in km / s %"ESYM" in code %"ESYM" and code vel = %"ESYM"\n", v_wind , v_wind *1.0E5/ VelocityUnits, VelocityUnits);

    v_wind     = (v_wind*1.0E5) / VelocityUnits; // convert from km/s to cm/s then to code units

    //p_feedback = m_eject * v_wind;
    float IndividualStarMaximumWindVelocity = 500.0, IndividualStarWindThermEfficiency = 0.9;
    float v_max = IndividualStarMaximumWindVelocity * 1.0E5 / VelocityUnits;
    //p_feedback = m_wind * IndividualStarMaximumWindVelocity;


    //E_thermal  = E_thermal * 1.0;
  //  p_feedback = m_eject * (IndividualStarMaximumWindVelocity*1.0E5 / VelocityUnits);
//    E_kin = 0.5 * m_eject * (IndividualStarMaximumWindVelocity*IndividualStarMaximumWindVelocity*1.0E10 / VelocityUnits / VelocityUnits);
    //E_kin = 0.0;

    v_wind     = IndividualStarMaximumWindVelocity*1.0E5 / VelocityUnits;

    E_thermal  = E_thermal + 0.5 * m_eject * (v_wind*v_wind) * IndividualStarWindThermEfficiency;
    E_kin      = 0.5 * m_eject * v_wind*v_wind * (1.0 - IndividualStarWindThermEfficiency);
    check_mass_in_region = 1; // thermal to KE switch
    E_kin = 0.0;
    p_feedback = m_eject * v_wind;

  } else if (mode == 1) {

    /* -------------------------------------------
     *
     * Compute the mass ejecta and energetics for
     * core collapse supernova
     *
     * -------------------------------------------
     */

    /* use yield tables to compute supernova ejecta mass - otherwise just eject the entire star
    if ( IndividualStarFollowStellarYields && TestProblemData.MultiMetals == 2){
      m_eject   = StellarYieldsInterpolateYield(0, mproj, metallicity, -1) * msolar / MassUnits;
    } else{
      m_eject   = StarMassEjectionFraction * mproj * msolar / MassUnits;
    }

    if( IndividualStarSupernovaEnergy < 0){
      E_thermal = m_eject * StarEnergyToThermalFeedback * (c_light * c_light/(VelocityUnits*VelocityUnits));
    } else {
      E_thermal = IndividualStarSupernovaEnergy * 1.0E51 / (MassUnits*VelocityUnits*VelocityUnits);
    }

    v_wind     = 0.0;
    p_feedback = 0.0;
    E_kin      = 0.0;
    printf("AddFeedbackGeneral: M_proj %"FSYM" Z = %"FSYM", M_eject = %"ESYM" E_thermal = %"ESYM"\n", mproj, metallicity, m_eject*MassUnits/msolar, E_thermal*VelocityUnits*VelocityUnits*MassUnits);
  } else if ( mode == 2){ // Type Ia supernova

    /* -------------------------------------------
     *
     * Compute the mass ejecta and energetics for
     * Type Ia supernova
     *
     * -------------------------------------------

    m_eject = StellarYields_SNIaYieldsByNumber(-1) * msolar / MassUnits; // total ejected mass

    if( IndividualStarSupernovaEnergy < 0){
      E_thermal = m_eject * StarEnergyToThermalFeedback * (c_light * c_light/(VelocityUnits*VelocityUnits));
    } else {
      E_thermal = IndividualStarSupernovaEnergy * 1.0E51 / (MassUnits*VelocityUnits*VelocityUnits);
    }

    v_wind     = 0.0;
    p_feedback = 0.0;
    E_kin      = 0.0;
  }


  /* if we are tracking yeilds, interpolat the ejecta mass for each species
  //printf("Tabulating metal ejecta mass fields \n");
  if(TestProblemData.MultiMetals == 2 && IndividualStarFollowStellarYields){
    int interpolation_mode;     // switch for stellar winds vs cc SN interpolation

    //
    // Switch around modes for interpolating either winds or supernova
    // to wrapper functions
    //
    if (mode == 1){
      interpolation_mode = 0;
    } else if (mode < 0){
      interpolation_mode = 1;
    }

    // for each metal species, compute the total metal mass ejected depending on supernova type
    if (mode < 0){

      metal_mass[0] = StellarYieldsInterpolateYield(1, mproj, metallicity, 0) *msolar / MassUnits / (dx*dx*dx);

      for (int i = 0; i < StellarYieldsNumberOfSpecies; i++){
        metal_mass[1 + i] = StellarYieldsInterpolateYield(1, mproj, metallicity, StellarYieldsAtomicNumbers[i]) * msolar / MassUnits / (dx*dx*dx);
      }

      // metal_mass now contains total mass ejected over wind lifetime. Adjust using wind loss rate and 
      // finite timestep check performed above
      for (int i = 0; i < StellarYieldsNumberOfSpecies + 1; i++){
        metal_mass[i] *= wind_dt / wind_lifetime;
      }

    } else if (mode == 1){
      // Core collapse supernova

      // First argument teslls interpolations to look at supernova yield tables
      // last argument (atomic number) of zero means get total metal mass
      metal_mass[0] = StellarYieldsInterpolateYield(0, mproj, metallicity, 0) * msolar / MassUnits / (dx*dx*dx);

      for(int i = 0; i < StellarYieldsNumberOfSpecies; i ++){
        metal_mass[1 + i] = StellarYieldsInterpolateYield(0, mproj, metallicity, StellarYieldsAtomicNumbers[i]) * msolar / MassUnits / (dx*dx*dx);
      }
    } else if (mode == 2){
      metal_mass[0] = StellarYields_SNIaYieldsByNumber(0) * msolar / MassUnits / (dx*dx*dx);

      for(int i = 0; i < StellarYieldsNumberOfSpecies + 1; i++){
        metal_mass[1 + i] = StellarYields_SNIaYieldsByNumber( StellarYieldsAtomicNumbers[i] ) * msolar / MassUnits / (dx * dx * dx);
      }
    }


    //printf("Metal masses in array ");
//    printf("ANUM = %"ISYM" %"ESYM " %"ESYM " -- \n", 0, (wind_lifetime/wind_dt) *metal_mass[0] * dx*dx*dx *MassUnits / msolar , metal_mass[0]);
//    for(int i = 0; i < StellarYieldsNumberOfSpecies; i++){
//      printf("ANUM = %"ISYM" %"ESYM " %"ESYM " -- \n", StellarYieldsAtomicNumbers[i], (wind_lifetime/wind_dt) *metal_mass[i+1] * dx*dx*dx *MassUnits / msolar , metal_mass[i+1]);
//    }
//    printf("\n");
  }


  m_eject    = m_eject    / (dx*dx*dx);   // now in code units (code mass / code volume)
  p_feedback = p_feedback / (dx*dx*dx);
  E_thermal  = E_thermal  / (dx*dx*dx);
  E_kin      = E_kin      / (dx*dx*dx);

  /* find coordinates of feedback center. This is nominally the particle position
     but is shifted if needed if particle is too close to grid boundary.
     This is taken care of below (xfc = x feedback center)
  float xfc, yfc, zfc;
  SetFeedbackCellCenter( xp, yp, zp, xstart, ystart, zstart, dx,
                           nx, ny, nz, ibuff, &xfc, &yfc, &zfc);

  //printf("ISF: injecting feedback to grid\n");
  this->IndividualStarInjectFeedbackToGrid(xfc, yfc, zfc,
                                           up, vp, wp,
                                           m_eject, E_thermal, E_kin, p_feedback, metal_mass,
                                           check_mass_in_region, E_thermal_wind); // function call

  // subtract mass from particle (in solar units)
  float new_mass = (*mp) - m_eject * (dx*dx*dx)*MassUnits/msolar;

  if( *mp < 0 && mode != 2){ // ignore this check for SNIa yields
      printf("new_mass = %"ESYM" mp = %"ESYM" m_eject =%"ESYM"\n", new_mass, *mp, m_eject*dx*dx*dx*MassUnits/msolar);
      ENZO_FAIL("IndividualStarFeedback: Ejected mass greater than current particle mass - negative particle mass!!!\n");
  }

  *mp = new_mass;

  delete [] metal_mass;

  return SUCCESS;
}

------------------------------------------ End old method tht needs to be excsied ---------------------------------


*/


void SetFeedbackCellCenter(const FLOAT &xp, const FLOAT &yp, const FLOAT &zp,
                             const FLOAT &xstart, const FLOAT &ystart, const FLOAT &zstart,
                             const FLOAT &dx,
                             const int &nx, const int &ny, const int &nz,
                             const int &ibuff,
                             FLOAT *xfc, FLOAT *yfc, FLOAT *zfc){

  /* checks if cell center is O.K. and rescales if needed */

  int fbuff = ibuff + ((int) (IndividualStarFeedbackStencilSize+1)/2.0 - 1); // number of cells away from edge
  float rnum;
  const int max_random = (1<<16);


  *xfc = xp;
  *yfc = yp;
  *zfc = zp;

  return;

/*
  if ( xp < xstart +         fbuff*dx ||
       xp > xstart + dx*nx - fbuff*dx ||
       yp < ystart +         fbuff*dx ||
       yp > ystart + dx*ny - fbuff*dx ||
       zp < zstart +         fbuff*dx ||
       zp > zstart + dx*nz - fbuff*dx   ){

    *xfc      = fmin( fmax(xp, xstart + fbuff*dx), xstart + dx*nx - fbuff*dx);
    *yfc      = fmin( fmax(yp, ystart + fbuff*dx), ystart + dx*ny - fbuff*dx);
    *zfc      = fmin( fmax(zp, zstart + fbuff*dx), zstart + dx*nz - fbuff*dx);

    xfcshift = xp - *xfc;
    yfcshift = yp - *yfc;
    zfcshift = zp - *zfc;

    printf("Warning: Shifting feedback zone %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM"\n", *xfc, *yfc, *zfc, xfcshift, yfcshift, zfcshift);
  } else{
    *xfc = xp;
    *yfc = yp;
    *zfc = zp;
  }
*/
}


int grid::IndividualStarInjectFeedbackToGrid(const FLOAT &xfc, const FLOAT &yfc, const FLOAT &zfc,
                               float up, float wp, float vp,
                               float m_eject, float E_thermal, float E_kin, float p_feedback, float *metal_mass,
                               int check_mass_in_region, float E_thermal_wind){

  float dx = float(CellWidth[0][0]);

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


  FLOAT xstart = this->CellLeftEdge[0][0],
        ystart = this->CellLeftEdge[1][0],
        zstart = this->CellLeftEdge[2][0];

  int nx = *(GridDimension), ny = *(GridDimension+1), nz = *(GridDimension+2);
  int size = nx*ny*nz;
  int number_of_cells;
  int stencil = IndividualStarFeedbackStencilSize; //renaming in case variable size
                                                   // is added in future
  /* for now, check total gas mass in injection region and print warning if too low */

  number_of_cells = POW(stencil, 3);
  /* scale everything to be the injected values in each cell */
  m_eject    = m_eject   / ((float) number_of_cells);
  E_thermal  = E_thermal / ((float) number_of_cells);
  p_feedback = p_feedback / ((float) number_of_cells - 1); // no momentum in center cell
  E_kin      = E_kin / ((float) number_of_cells - 1); // E_kin is totoal?   //  / ((float) number_of_cells);

  /* */

  float face_shift;
  float xface, yface, zface, dxf, dyf, dzf;
  float xpos, ypos, zpos, dxc, dyc, dzc;
  int iface, jface, kface, ic, jc, kc;

  /* check hydro method and shift things around
     to face-centered */
  face_shift = 0.0;
  if (HydroMethod == 2){ face_shift = 0.5; }

  /* AJE: I Suspect -0.5 should not be there in xface/xpos - June 2016*/

  xface = (xfc - xstart)/dx  - face_shift;
  yface = (yfc - ystart)/dx  - face_shift;
  zface = (zfc - zstart)/dx  - face_shift;

  iface = ((int) floor(xface + 0.5));
  jface = ((int) floor(yface + 0.5));
  kface = ((int) floor(zface + 0.5));

  dxf = iface + 0.5 - xface;
  dyf = jface + 0.5 - yface;
  dzf = kface + 0.5 - zface;

  /* we now need the index of the cell to add mass */
  xpos = (xfc - xstart)/dx;
  ypos = (yfc - ystart)/dx;
  zpos = (zfc - zstart)/dx;

  ic   = ((int) floor(xpos + 0.5));
  jc   = ((int) floor(ypos + 0.5));
  kc   = ((int) floor(zpos + 0.5));

  dxc  = ic + 0.5 - xpos;
  dyc  = jc + 0.5 - ypos;
  dzc  = kc + 0.5 - zpos;

  /* allocate local field stencil - bigger than feedback stencil */
  float *u_local, *v_local, *w_local, *d_local, *ge_local, *te_local;
  float *ke_before;

  float additional_mass_factor = 1.0;

  int local_number_of_cells = POW(stencil + 1,3);

  u_local  = new float[local_number_of_cells];
  v_local  = new float[local_number_of_cells];
  w_local  = new float[local_number_of_cells];
  d_local  = new float[local_number_of_cells];
  ge_local = new float[local_number_of_cells];
  te_local = new float[local_number_of_cells];


  if (HydroMethod != 2){
    ke_before = new float[local_number_of_cells];
  } else {ke_before = NULL;}

  /* should go up to stencil + 1 I think (stencil = 3, 0, 1, 2, 3) */
  for(int k = 0; k < stencil + 1; k++){
    for(int j = 0; j < stencil + 1; j++){
      for(int i = 0; i < stencil + 1; i++){
        int index = i + (j + k * (stencil + 1)) * (stencil + 1);
        u_local[index] = 0.0;
        v_local[index] = 0.0;
        w_local[index] = 0.0;
        d_local[index] = 0.0;
        ge_local[index] = 0.0;
        te_local[index] = 0.0;

        if( HydroMethod != 2){
          ke_before[index] = 0.0;
        }
      }
    }
  }
  /* done allocating the zeroed arrays */

  /* Check mass in region - if too small, switch to kinetic only */
  if (check_mass_in_region) {
      printf("Checking mass in region\n");
      float * temperature;
      temperature = new float[size];

      if(  this->ComputeTemperatureField(temperature) == FAIL ){
        ENZO_FAIL("Error in compute temperature called from PhotoelectricHeatingFromStar");
      }


      float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits, MassUnits, EnergyUnits;
      if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits, &TimeUnits, &VelocityUnits, this->Time) == FAIL){
          ENZO_FAIL("Error in GetUnits");
      }
      MassUnits   = DensityUnits*LengthUnits*LengthUnits*LengthUnits; // mass units
      EnergyUnits = MassUnits * VelocityUnits * VelocityUnits;

      float cell_mass_safety_factor = 50.0 ; // MAKE PARAMETER

      int integer_sep = ((int) (stencil + 1)/2.0 - 1);
      float cell_mass = 1.0E99;
      float max_temperature = 1.0;
      int min_index = -1;
      for(int k = kc - integer_sep - 1 ; k <= kc + integer_sep; k ++){
        for(int j = jc - integer_sep - 1 ; j <= jc + integer_sep; j++){
          for(int i = ic - integer_sep - 1; i <= ic + integer_sep; i++){
            /* check this index, I think it should be nx ny not stencil+1 */
            int index = i + (j + k * (ny))*(nx);
            if ( temperature[index] > max_temperature){
              min_index = index;
              cell_mass = BaryonField[DensNum][index];
              max_temperature = temperature[index];
            }

          }
        }
      }

      /* now adjust feedback appropriately */
      float est_mu = 0.4; // estimated - this is approximate anyway
      float k_boltz = 1.380658E-16;
      float mp     = 1.6733E-24;
      float T = temperature[min_index]; // * TemperatureUnits;

      E_thermal = E_thermal  * dx *dx *dx * EnergyUnits;
      cell_mass = cell_mass* dx *dx *dx * MassUnits;
      m_eject   = m_eject * dx *dx *dx * MassUnits;
      E_kin     = E_kin * dx *dx *dx *EnergyUnits;

      printf("Minimum Index %i - T = %E - M = %E \n", min_index, T, cell_mass);
      printf("E_thermal = %E - E_kin = %E\n", E_thermal, E_kin);

      float T_final = (E_thermal + 1.5 * cell_mass * k_boltz * T / (est_mu*mp)) *
                      (2.0 * est_mu * mp/(3.0 * k_boltz * (cell_mass + m_eject)));

      float T_wind_threshold = 1.0E6;
      float maximum_velocity = 300.0 * 1.0E5;

      if (T > T_wind_threshold){ // already hot / too hot - momentum!
        additional_mass_factor = 1.85;

      } else if (T_final > T_wind_threshold){
        additional_mass_factor = 1.85;
     }



     E_kin = E_kin / (dx*dx*dx) / EnergyUnits;
     E_thermal = E_thermal / (dx*dx*dx) / EnergyUnits;
     m_eject = m_eject / (dx*dx*dx) / MassUnits;
     p_feedback = p_feedback / (dx*dx*dx) / (MassUnits * VelocityUnits);
      // maybe switch to momentum only at some point???



    delete[] temperature;
  }


  /* add up kinetic energy in the clocal region */
  if( HydroMethod != 2 ){
    int integer_sep = ((int) (stencil+1)/2.0 - 1); //floor((stencil + 1) / 2.0);

    int local_index = 0; // AJE 5 - 10 - 16
/* changing + 1 to - 1 here - - */

    for(int k = kc - integer_sep - 1 ; k <= kc + integer_sep; k ++){
      for(int j = jc - integer_sep - 1 ; j <= jc + integer_sep; j++){
        for(int i = ic - integer_sep - 1; i <= ic + integer_sep; i++){
          /* check this index, I think it should be nx ny not stencil+1 */
          int index = i + (j + k * (ny))*(nx);


          if( index >= 0 && index < size){
            ke_before[local_index] = 0.5 * BaryonField[DensNum][index] *
                                   ( (BaryonField[Vel1Num][index] * BaryonField[Vel1Num][index]) +
                                     (BaryonField[Vel2Num][index] * BaryonField[Vel2Num][index]) +
                                     (BaryonField[Vel3Num][index] * BaryonField[Vel3Num][index]));
          } else{
            ke_before[local_index] = 0.0; // off grid but give some value
          }

          local_index++;
        }
      }
    }
  } // end kinetic energy sum for Zeus


  /* function call to convert velocities to momentum in particle frame */
  Momentum(BaryonField[Vel1Num], BaryonField[Vel2Num], BaryonField[Vel3Num], BaryonField[DensNum],
           up, vp, wp, nx, ny, nz, ic, jc, kc, iface, jface, kface, stencil, 1);
  /* end convert to momenta */


/*
  AJE 9/13/16: Should not need this anymore since we moved feedback out of
               Grid_StarParticleHandler


  if(IndividualStarFollowStellarYields && TestProblemData.MultiMetals == 2){
    int field_num;
    this->IdentifyChemicalTracerSpeciesFieldsByNumber(field_num, 0);
    MetalConversion( BaryonField[field_num], BaryonField[DensNum],
                     dx, nx, ny, nz, ic, jc, kc, stencil, 1);

    for(int m = 0; m < StellarYieldsNumberOfSpecies; m++){
      this->IdentifyChemicalTracerSpeciesFieldsByNumber(field_num, StellarYieldsAtomicNumbers[m]);

      MetalConversion( BaryonField[field_num], BaryonField[DensNum],
                       dx, nx, ny, nz, ic, jc, kc, stencil, 1);
    }
  }


*/


  /* compute the total mass and energy in cells before the explosion */
  float mass_before, energy_before, kin_energy_before;

  SumMassEnergy(BaryonField[Vel1Num], BaryonField[Vel2Num], BaryonField[Vel3Num], BaryonField[DensNum],
                  BaryonField[GENum],   BaryonField[TENum],
                  nx, ny, nz, iface, jface, kface, ic, jc, kc, stencil,
                  &mass_before, &energy_before, &kin_energy_before);


  /* Now add mass and momentum terms to the dummy fields constructed earlier */
  AddFeedbackToGridCells(u_local, v_local, w_local, d_local, ge_local, te_local,
                         stencil+1, stencil+1, stencil+1, 1, 1, 1, 1, 1, 1,
                         dxf, dyf, dzf, dxc, dyc, dzc, m_eject, 1.0, 0.0, stencil, additional_mass_factor);

  /* quadratic equation in p - momentum energy injection */
  float mom_per_cell = 0.0;
  if (E_kin > 0) {
    float A, B, C;
    ComputeAbcCoefficients( BaryonField[Vel1Num], BaryonField[Vel2Num], BaryonField[Vel3Num],
                            BaryonField[DensNum], BaryonField[GENum]  ,
                            u_local, v_local, w_local, d_local,
                            nx, ny, nz, iface, jface, kface, ic, jc, kc, stencil,
                            A, B, C);

    A = A - (kin_energy_before + E_kin);
    mom_per_cell  = (-B + sqrt(B*B - 4.0 * A * C)) / (2.0 * C);

    //printf("ISF: Coeffs mom_per_cell %"ESYM"\n", mom_per_cell);
  } else { // no kinetic energy injection - add feedback will add mass directly
    mom_per_cell = p_feedback;
  }

  /* add metal feedback - mass in cells */
  //printf("ISF: Starting metal injection feedback calls\n");
  if(TestProblemData.MultiMetals == 2 && IndividualStarFollowStellarYields && (metal_mass)){
    /* For the first call, add in general metallicity field */
    int field_num;
    this->IdentifyChemicalTracerSpeciesFieldsByNumber(field_num, 0); // when atomic number is zero, gives metallicity field

    AddMetalSpeciesToGridCells(BaryonField[field_num], metal_mass[0] / ((float) number_of_cells),
                               nx, ny, nz, ic, jc, kc, dxc, dyc, dzc, stencil, additional_mass_factor);

    for(int ii = 0; ii < StellarYieldsNumberOfSpecies; ii++){
      this->IdentifyChemicalTracerSpeciesFieldsByNumber(field_num, StellarYieldsAtomicNumbers[ii]);

       //printf("ISF: Adding metal feedaback for field %"ISYM" and atomic number %"ISYM" %"ISYM"\n", field_num, StellarYieldsAtomicNumbers[ii], ii);
       AddMetalSpeciesToGridCells(BaryonField[field_num], metal_mass[1 + ii] / ((float) number_of_cells),
                                  nx, ny, nz, ic, jc, kc, dxc, dyc, dzc, stencil, additional_mass_factor);

    }

  } // 


  /* Now call add feedback again to add the feedback into the grid cells */
  AddFeedbackToGridCells(BaryonField[Vel1Num], BaryonField[Vel2Num], BaryonField[Vel3Num],
                         BaryonField[DensNum], BaryonField[GENum],   BaryonField[TENum]  ,
                         nx, ny, nz, ic, jc, kc, iface, jface, kface,
                         dxf, dyf, dzf, dxc, dyc, dzc, m_eject, mom_per_cell, E_thermal,
                         stencil, additional_mass_factor);



  /* compute total mass and energy after feedback has been added */
  float mass_after, energy_after, kin_energy_after;
  SumMassEnergy(BaryonField[Vel1Num], BaryonField[Vel2Num], BaryonField[Vel3Num],
                    BaryonField[DensNum], BaryonField[GENum],   BaryonField[TENum],
                    nx, ny, nz, iface, jface, kface, ic, jc, kc, stencil,
                    &mass_after, &energy_after, &kin_energy_after);

  /* error checking statments go here */
  //printf("ISF energy cheks, before %"ESYM" after %"ESYM" eject %"ESYM"\n", energy_before, energy_after, E_thermal);
  //printf("ISF Mass checks, before %"ESYM" after %"ESYM" eject %"ESYM"\n", mass_before, mass_after, m_eject);
  //printf("ISF Kinetic energy checks, before, after, E_kin %"ESYM" %"ESYM" %"ESYM"\n", kin_energy_before, kin_energy_after, E_kin);
  /*           -------------           */


  /* Now, reset the velocity fields to simulation frame */
  Momentum(BaryonField[Vel1Num], BaryonField[Vel2Num], BaryonField[Vel3Num],
           BaryonField[DensNum], up, vp, wp, nx, ny, nz, ic, jc, kc,
           iface, jface, kface, stencil, -1);

/*

  AJE 9/13/16: Should not need this anymore since we moved feedback out of
               Grid_StarParticleHandler

  // If needed, convert metals back to metal fractions
  if(IndividualStarFollowStellarYields && TestProblemData.MultiMetals == 2){
    int field_num;
    this->IdentifyChemicalTracerSpeciesFieldsByNumber(field_num, 0);
    MetalConversion( BaryonField[field_num], BaryonField[DensNum],
                     dx, nx, ny, nz, ic, jc, kc, stencil, -1);


    for(int m = 0; m < StellarYieldsNumberOfSpecies; m++){
      int field_num;
      this->IdentifyChemicalTracerSpeciesFieldsByNumber(field_num, StellarYieldsAtomicNumbers[m]);

      MetalConversion( BaryonField[field_num], BaryonField[DensNum],
                       dx, nx, ny, nz, ic, jc, kc, stencil, -1);
    }
  }
*/
  /* Adjust the total energy field if we are using PPM */
  if (HydroMethod != 2){
    float ke_injected = 0.0;
    float delta_ke = 0.0;
    float ke_after = 0.0;
    int integer_sep = ((int) (stencil+1)/2.0 - 1); // floor((stencil + 1) / 2.0);
    //printf("ISF kinetic feedback: integer_separation = %"ISYM"\n",integer_sep);

/* CHanging + 1 to -1 AJE 8/17 */

    int local_index = 0;
    for(int k = kc - integer_sep - 1; k <= kc + integer_sep; k ++){
      for(int j = jc - integer_sep - 1; j <= jc + integer_sep ; j++){
        for(int i = ic - integer_sep - 1; i <= ic + integer_sep; i++){


//    for(int k = -integer_sep; k < integer_sep; k++){
//      for(int j = -integer_sep; j < integer_sep; j++){
//        for(int i = -integer_sep; i < integer_sep; i++){

          int index   = i + ( j + k*ny)*nx;

          if (index >= 0 && index < size){

            ke_after = 0.5 * BaryonField[DensNum][index] *
                      ( BaryonField[Vel1Num][index] * BaryonField[Vel1Num][index] +
                        BaryonField[Vel2Num][index] * BaryonField[Vel2Num][index] +
                        BaryonField[Vel3Num][index] * BaryonField[Vel3Num][index]);

            delta_ke = ke_after - ke_before[local_index];

            BaryonField[TENum][index] += delta_ke/BaryonField[DensNum][index];

            ke_injected += delta_ke;

          } // off grid if not true... do nothing

          local_index++;
        }
      }
    }
    //printf("IndividualStarFeedback: change in kinetic energy %"ESYM"\n", ke_injected);
  } // endif


  /* free up memory */
  delete[] u_local;
  delete[] v_local;
  delete[] w_local;
  delete[] d_local;
  delete[] ge_local;
  delete[] te_local;
  delete[] ke_before;

  return SUCCESS;
}


void Momentum(float *u, float *v, float *w, float *d,
              const float &up, const float &vp, const float &wp,
              const int &nx, const int &ny, const int &nz,
              const int &ic, const int &jc, const int &kc,
              const int &iface, const int &jface, const int &kface,
              const int &stencil, int idir){

  int xo, yo, zo;
  // unit one shift in each direction
  xo = 1;
  yo = nx;
  zo = (nx * ny);

  /* making this a bigger region.... for stencil = 3, do -2 -1 0 1 2  -- 4/20/16 */
//  printf("Momentum Conversion: Direction %"ISYM"\n",idir);
//  int integer_sep = floor((stencil+1)/2.0) + 1; - AJE 5/10/16
/* CHanging +1 to -1 AJE 8/17 */
  int integer_sep = ( (int) (stencil+1)/2.0 -1);
  for(int k = -integer_sep -1 ; k <= integer_sep; k++){
    for(int j = -integer_sep -1 ; j <= integer_sep; j++){
      for(int i = -integer_sep -1 ; i <= integer_sep; i++){

        int index=0, x_index=0, y_index = 0, z_index = 0;
        if ( idir == 1.0) {
          if ( HydroMethod == 2 ){

            x_index = (iface + i) + ( (jc    + j) + (kc    + k) * ny) * nx;
            y_index = (ic    + i) + ( (jface + j) + (kc    + k) * ny) * nx;
            z_index = (ic    + i) + ( (jc    + j) + (kface + k) * ny) * nx;

  //          printf("Input u = %"ESYM" v = %"ESYM" w = %"ESYM"\n", u[x_index], v[y_index], w[z_index]);

            u[x_index] = (u[x_index] - up)*( 0.5 *(d[x_index] + d[x_index + xo]));
            v[y_index] = (v[y_index] - vp)*( 0.5 *(d[y_index] + d[y_index + yo]));
            w[z_index] = (w[z_index] - wp)*( 0.5 *(d[z_index] + d[z_index + zo]));

    //        printf("Output u = %"ESYM" v = %"ESYM" w = %"ESYM"\n", u[x_index], v[y_index], w[z_index]);

          } else {
            index = (ic + i) + ( (jc + j) + (kc + k)*ny)*nx;

            if (index < 0 || index >= nx*ny*nz){
              continue; // off grid
            }

            u[index] = (u[index] - up)*d[index];
            v[index] = (v[index] - vp)*d[index];
            w[index] = (w[index] - wp)*d[index];
          } // hydro method check


        } else { // reverse
          if( HydroMethod == 2){

            x_index = (iface + i) + ( (jc    + j) + (kc    + k) * ny) * nx;
            y_index = (ic    + i) + ( (jface + j) + (kc    + k) * ny) * nx;
            z_index = (ic    + i) + ( (jc    + j) + (kface + k) * ny) * nx;


      //      printf("Input u = %"ESYM" v = %"ESYM" w = %"ESYM"\n", u[x_index], v[y_index], w[z_index]);

            u[x_index] = u[x_index] / (0.5 * (d[x_index] + d[x_index + xo])) + up;
            v[y_index] = v[y_index] / (0.5 * (d[y_index] + d[y_index + yo])) + vp;
            w[z_index] = w[z_index] / (0.5 * (d[z_index] + d[z_index + zo])) + wp;

        //    printf("Output u = %"ESYM" v = %"ESYM" w = %"ESYM"\n", u[x_index], v[y_index], w[z_index]);
          } else{
            index = (ic + i) + ( (jc + j) + (kc + k) *ny)*nx;

            if( index < 0 || index >= nx*ny*nz){
              continue; // off grid
            }

            u[index] = u[index] / d[index] + up;
            v[index] = v[index] / d[index] + vp;
            w[index] = w[index] / d[index] + wp;
          }
        } // end idir check

        //printf("momentum_indexes: x y z %"ISYM" %"ISYM" %"ISYM"\n",x_index, y_index, z_index);
      }
    }
  } // end loop


// done
}

void MetalConversion(float *m, float *d, const float &dx,
                     const int &nx, const int &ny, const int &nz,
                     const int &ic, const int &jc, const int &kc,
                     const int &stencil, int idir){
 /* -----------------------------------------------------------------
  * MetalConversion
  * -----------------------------------------------------------------
  * A. Emerick - 04/19/16
  * ----------------------
  * Converts an arbitrary metal field from metal density to metal
  * density (a proxy for mass since we are operating on a grid with
  * uniform dx) and vice versa depending on idir. Feedback injection
  * occurs in terms of metall density (mass), field stored as fraction
  * ------------------------------------------------------------------- */

  /* any add all changes to momentum above should (probably) be reflected here as well */
  /* Metal fields give metal density in a given cell - convert to mass */
//  int integer_sep = floor((stencil+1)/2.0) + 1;
  int integer_sep = ((int) (stencil+1)/2.0 - 1);
/* changing +1 to -1 AJE 8/17 */
  for(int k = -integer_sep - 1; k <= integer_sep ; k++){
    for(int j = -integer_sep - 1; j <= integer_sep ; j++){
      for(int i = -integer_sep - 1; i <= integer_sep ; i++){

        int index = (ic + i) + ( (jc + j) + (kc + k) * ny) * nx;

        if (index < 0 || index >= nx*ny*nz){
          continue; // off grid
        }


        if (idir == 1.0){
          m[index] = m[index] / d[index];
        } else {
          m[index] = m[index] * d[index];
        }

      } // end i
    } // end j
  } // end k
} // done with metal conversion

void SumMassEnergy(float *pu, float *pv, float *pw, float *d, float *ge, float *te,
                   const int &nx, const int &ny, const int &nz, 
                   const int &iface, const int &jface, const int &kface,
                   const int &ic, const int &jc, const int &kc, const int &stencil,
                   float *mass_sum, float *energy_sum, float *kin_energy_sum){

  *mass_sum = 0.0; *energy_sum = 0.0; *kin_energy_sum = 0.0;

  int xo, yo, zo;
  // unit one shift in each direction
  xo = 1;
  yo = nx;
  zo = (nx * ny);

//  int integer_sep = floor((stencil+1)/2.0);
  int integer_sep = ((int) (stencil+1)/2.0 - 1);
/*changing +1 to -1 AJE 8/17 */
  for(int k = -integer_sep -1; k <= integer_sep; k++){
    for(int j = -integer_sep -1; j <= integer_sep; j++){
      for(int i = -integer_sep -1; i <= integer_sep; i++){
        float mass_term, mom_term, gas_energy = 0.0, kinetic_energy;

        int index   = (ic + i) + ( (jc + j) + (kc + k)*ny)*nx;
        int x_index = (iface + i) + ( (jc    + j) + (kc    + k)*ny)*nx;
        int y_index = (ic    + i) + ( (jface + j) + (kc    + k)*ny)*nx;
        int z_index = (ic    + i) + ( (jc    + j) + (kface + k)*ny)*nx;

        if ( (index   < 0 || index >= nx*ny*nz) ||
             (x_index < 0 || x_index >= nx*ny*nz) ||
             (y_index < 0 || y_index >= nx*ny*nz) ||
             (z_index < 0 || z_index >= nx*ny*nz) ) {
          continue; // off grid
        }

        mass_term = d[index];
        mom_term  = pu[x_index]*pu[x_index] +
                    pv[y_index]*pv[y_index] +
                    pw[z_index]*pw[z_index] ;

        /* total mass and energy */
        kinetic_energy   = mom_term / (2.0 * mass_term);
        *mass_sum        = *mass_sum + mass_term;
        *kin_energy_sum  = *kin_energy_sum + kinetic_energy;

        /* account for thermal energy depending on Hydro solver */
        if (HydroMethod == 2){
          gas_energy = te[index] * d[index];
        }
        if (DualEnergyFormalism) {
          gas_energy = ge[index]*d[index];
        }
        if (HydroMethod != 2 && DualEnergyFormalism == 0) {
          gas_energy = te[index]*d[index] - kinetic_energy;
        }

        *energy_sum = *energy_sum + kinetic_energy + gas_energy;

      }
    }
  }// end loop


// done with routine
}


void ComputeAbcCoefficients(float *pu, float *pv, float *pw, float *d,
                            float *ge, float *pu_l, float *pv_l, float *pw_l,
                            float *d_l, const int &nx, const int &ny, const int &nz,
                            const int &iface, const int &jface, const int &kface,
                            const int &ic, const int &jc, const int &kc,
                            const int &stencil, float &A, float &B, float &C){

  float mass_sum, energy_sum;

  A = 0.0; B = 0.0; C = 0.0; mass_sum = 0.0; energy_sum = 0.0;

  int xo, yo, zo;
  // unit one shift in each direction
  xo = 1;
  yo = nx;
  zo = (nx * ny);

  float mass_term=0.0, mom_term=0.0, b_term=0.0, c_term = 0.0;

  int integer_sep = ((int) (stencil+1)/2.0 - 1);
  int l_index     = 0; int index = 0;

/* changing + 1 to -1 AJE 8/17 */

  for(int k = -integer_sep - 1; k <= integer_sep ; k++){
    for(int j = -integer_sep - 1; j <= integer_sep ; j++){
      for(int i = -integer_sep -1; i <= integer_sep ; i++){
        int x_index, y_index, z_index;


        index   = (ic + i) + ( (jc + j) + (kc + k)*ny)*nx;
        x_index = (iface + i) + ( (jc    + j) + (kc    + k)*ny)*nx;
        y_index = (ic    + i) + ( (jface + j) + (kc    + k)*ny)*nx;
        z_index = (ic    + i) + ( (jc    + j) + (kface + k)*ny)*nx;

        if ( (index   < 0 || index >= nx*ny*nz) ||
             (x_index < 0 || x_index >= nx*ny*nz) ||
             (y_index < 0 || y_index >= nx*ny*nz) ||
             (z_index < 0 || z_index >= nx*ny*nz) ) {
          continue; // off grid
        }


        mass_term = d[index];
        mom_term  = pu[x_index]*pu[x_index] +
                    pv[y_index]*pv[y_index] +
                    pw[z_index]*pw[z_index] ;


        mass_term += d_l[l_index];

        b_term     = pu[x_index]*pu_l[l_index] +
                     pv[y_index]*pv_l[l_index] +
                     pw[z_index]*pw_l[l_index]  ;

        c_term     = pu_l[l_index] * pu_l[l_index] +
                     pv_l[l_index] * pv_l[l_index] +
                     pw_l[l_index] * pw_l[l_index] ;


        A         += mom_term/(2.0 * mass_term);
        B         += b_term/mass_term;
        C         += c_term/(2.0 * mass_term);

        l_index++;
      }
    }
  }
  //printf("ComputeAbcCoefficients: local_index = %"ISYM" integer_sep = %"ISYM" A = %"ESYM" B = %"ESYM" C = %"ESYM"\n", l_index, index, A, B, C);


} // end comput coeff

void AddMetalSpeciesToGridCells(float *m, const float &mass_per_cell,
                                const int &nx, const int &ny, const int &nz,
                                const int &ic, const int &jc, const int &kc,
                                const float &dxc, const float &dyc, const float &dzc,
                                const int &stencil, const float additional_mass_factor){
 /* -------------------------------------------------------------------------
  * AddMetalSpeciesToGridCells
  * --------------------------------------------------------------------------
  * Adds in metal species deposition for a given metal field and mass ejection.
  * This is a copy / adaptation of the below algorithm (AddFeedbackToGridCells).
  * Any modification to that function should be reflected here as well.
  * I know this is gross, but it is somewhat more efficient.
  * ------------------------------------------------------------------------- */

  int integer_sep = ((int) (stencil+1)/2.0 - 1);
  float delta_mass = 0.0, total_mass = 0.0;

  for (int k = -integer_sep; k <= integer_sep; k++){
    for (int j = -integer_sep; j <= integer_sep; j++){
      for (int i = -integer_sep; i <= integer_sep; i++){

        /* + 1 to -1 AJE 8/17 be careful here */
        for(int i_loc = i - 1; i_loc <= i; i_loc++){
          float dxc_loc = ( (i_loc == i) ? dxc : 1.0 - dxc);

          for(int j_loc = j - 1; j_loc <= j ; j_loc++){
            float dyc_loc = ( (j_loc == j) ? dyc : 1.0 - dyc);

            for( int k_loc = k - 1; k_loc <= k ; k_loc++){
              float dzc_loc = ( (k_loc == k) ? dzc : 1.0 - dzc);


              int index = (ic + i_loc) + ( (jc + j_loc) + (kc + k_loc)*ny) * nx;

              if (index < 0 || index >= nx*ny*nz){
                continue; // off grid
              }

              delta_mass = mass_per_cell * dxc_loc * dyc_loc * dzc_loc * additional_mass_factor;

              m[index] = m[index] + delta_mass;

              total_mass += delta_mass;
            } //

          } //

        } //

      }
    }
  } // end k loop 

  //printf("MetalFeedback: Deposited total metal mass (density) %"ESYM"\n", total_mass);
}

void AddFeedbackToGridCells(float *pu, float *pv, float *pw,
                            float *d, float * ge, float *te,
                            const int &nx, const int &ny, const int &nz,
                            const int &ic, const int &jc, const int &kc,
                            const int &iface, const int &jface, const int &kface,
                            const float &dxf, const float &dyf, const float &dzf,
                            const float &dxc, const float &dyc, const float &dzc,
                            const float &mass_per_cell, const float &mom_per_cell,
                            const float &therm_per_cell, const int &stencil,
                            const float additional_mass_factor){


  int integer_sep = ((int) (stencil+1)/2.0 - 1);
  float total_mass = 0.0, delta_therm =0.0;

  int on_grid, off_grid;
  on_grid = 0; off_grid = 0;

  // should go over a stencilxstencilxstencil region (if stencil is 3, -1, 0, 1)
  for(int k = -integer_sep; k <= integer_sep; k++){
    for(int j = -integer_sep; j <= integer_sep; j++){
      for(int i = -integer_sep; i <= integer_sep; i++){

/* + 1 to - 1 here - be careful */

        for(int i_loc = i - 1; i_loc <= i; i_loc++){
          float dxf_loc, dxc_loc;

          dxf_loc = dxf;  dxc_loc = dxc;
          if( i_loc == i){
            dxf_loc = 1.0 - dxf;
            dxc_loc = 1.0 - dxc;
          }

          for(int j_loc = j - 1; j_loc <= j; j_loc++){
            float dyf_loc, dyc_loc;

            dyf_loc = dyf;  dyc_loc = dyc;
            if( j_loc == j){
              dyf_loc = 1.0 - dyf;
              dyc_loc = 1.0 - dyc;
            }

            for(int k_loc = k - 1; k_loc <= k; k_loc++){
              float dzf_loc, dzc_loc;

              dzf_loc = dzf;    dzc_loc = dzc;
              if( k_loc == k){
                dzf_loc = 1.0 - dzf;
                dzc_loc = 1.0 - dzc;
              }
              /* now we do the feedback */
              int index, x_index, y_index, z_index;
              float delta_mass, delta_pu, delta_pv, delta_pw, inv_dens;

              index = (ic + i_loc) + ( (jc + j_loc) + (kc + k_loc)*ny)*nx;

              x_index = (iface + i_loc) + ( (jc    + j_loc) + (kc    + k_loc)*ny)*nx;
              y_index = (ic    + i_loc) + ( (jface + j_loc) + (kc    + k_loc)*ny)*nx;
              z_index = (ic    + i_loc) + ( (jc    + j_loc) + (kface + k_loc)*ny)*nx;

              if ( (index   < 0 || index >= nx*ny*nz) ||
                 (x_index < 0 || x_index >= nx*ny*nz) ||
                 (y_index < 0 || y_index >= nx*ny*nz) ||
                 (z_index < 0 || z_index >= nx*ny*nz) ) {
                off_grid++;
                continue; // off grid
              }
              on_grid++;

              /* do things here finally */
              delta_mass  = mass_per_cell * dxc_loc * dyc_loc * dzc_loc * additional_mass_factor;

              // add momentum, then do sign changing later
              delta_pu    = mom_per_cell * dxf_loc * dyc_loc * dzc_loc;
              delta_pv    = mom_per_cell * dxc_loc * dyf_loc * dzc_loc;
              delta_pw    = mom_per_cell * dxc_loc * dyc_loc * dzf_loc;

              // change sign to point away from central sell
              // +1 or -1 if index is pos or neg respectively
              // multiply by zero if we are the central cell
              delta_pu *=  ( (i > 0) ? 1.0 : ( i < 0 ? -1.0 : 0));
              delta_pv *=  ( (j > 0) ? 1.0 : ( j < 0 ? -1.0 : 0));
              delta_pw *=  ( (k > 0) ? 1.0 : ( k < 0 ? -1.0 : 0));

              delta_therm = therm_per_cell * dxc_loc * dyc_loc * dzc_loc;

              /* add mass momentum and thermal energy */
              inv_dens   = 1.0 / (d[index] + delta_mass);

              int mom_norm = 3;
              if ( i == 0) mom_norm--;
              if ( j == 0) mom_norm--;
              if ( k == 0) mom_norm--;
              if (mom_norm == 0) mom_norm = 1;
              mom_norm = 1;

              pu[x_index] = pu[x_index] + delta_pu / (sqrt((float) mom_norm));
              pv[y_index] = pv[y_index] + delta_pv / (sqrt((float) mom_norm));
              pw[z_index] = pw[z_index] + delta_pw / (sqrt((float) mom_norm));

              total_mass  += delta_mass;

              te[index] = (te[index]*d[index] + delta_therm) * inv_dens;

              if (DualEnergyFormalism) {
                ge[index] = (ge[index]*d[index] + delta_therm) * inv_dens;
              }
              d[index] = d[index] + delta_mass;

              //if( d[index] < 0.0 || delta_mass < 0.0){
              //  printf("individual_star: Feedback producing negative densities %"ESYM" %"ESYM" %"ESYM" %"ESYM"\n", d[index], delta_mass, te[index], delta_therm);
              //  printf("--------------------------------- %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM"\n",mass_per_cell,dxf_loc, dxc_loc, dyf_loc, dyc_loc, dzf_loc, dzc_loc);
              //}
          //    printf("feedback_indexes: i x y z %"ISYM" %"ISYM" %"ISYM" %"ISYM"\n",index, x_index, y_index, z_index);

            } // k loc
          } // j loc
        } // i loc


      }
    }
  }//end loop
  //printf("AddFeedbackToGridCells: on_grid = %"ISYM" off_grid = %"ISYM" total = %"ISYM"\n", on_grid, off_grid, on_grid+off_grid);
  //printf("AddedFeedbackToGridCells: mom_per_cell = %"ESYM" therm_per_cell = %"ESYM"\n", mom_per_cell, therm_per_cell);
  //printf("AddedFeedbackToGridCells: M_tot =  %"ESYM" therm = %"ESYM"\n", total_mass, delta_therm);

}

//
//
// TODO (to do, To Do): Oct 2016 - clean up function call and variables now that it gets star class object
//


int grid::IndividualStarAddFeedbackSphere(Star *cstar, const FLOAT &xp, const FLOAT &yp, const FLOAT &zp,
                                          const float &up, const float &vp, const float &wp, // might not need vel
                                          const float &mproj, const float &lifetime, const float &particle_age,
                                          const float &metallicity, float *mp, int mode){

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
  const double msolar = 1.989E33;

  float *metal_mass; // array of individual species masses

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
  if(IndividualStarFollowStellarYields && TestProblemData.MultiMetals == 2){

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
  }

  /* convert computed parameters to code units */
  m_eject   = m_eject*msolar / MassUnits   / (dx*dx*dx);
  E_thermal = E_thermal      / EnergyUnits / (dx*dx*dx);

  if(IndividualStarFollowStellarYields && TestProblemData.MultiMetals == 2){
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
    this->IndividualStarInjectSphericalFeedback(cstar, xp, yp, zp, m_eject, E_thermal,
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
                                                const FLOAT &xp, const FLOAT &yp, const FLOAT &zp,
                                                float m_eject, float E_thermal,
                                                float *metal_mass, int stellar_wind_mode){

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

  int nx = *(this->GridDimension), ny = *(this->GridDimension+1), nz = *(this->GridDimension+2);
  int size = nx*ny*nz;
  int number_of_cells;

  float xpos, ypos, zpos;
  int ic, jc, kc;

  //
  // find position and index for grid zone
  // nearest to particle
  //
  xpos = (xp - xstart)/dx;
  ypos = (yp - ystart)/dx;
  zpos = (zp - zstart)/dx;

  ic   = ((int) floor(xpos ));
  jc   = ((int) floor(ypos ));
  kc   = ((int) floor(zpos ));

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
  float radius = IndividualStarFeedbackStencilSize * dx;    // code length
  float volume = 4.0 * pi * radius * radius * radius / 3.0; // (code length)**3
  int   r_int  = ceil(IndividualStarFeedbackStencilSize);     // int of farthest cell in any dir.
  float cell_volume_fraction = dx*dx*dx / volume;           // fraction of vol for each cell

  float injected_metal_mass[StellarYieldsNumberOfSpecies+1];

  if (metal_mass == NULL && (cstar)){
    injected_metal_mass[0] = cstar->ReturnMetallicity() * m_eject;
  } else {
    injected_metal_mass[0] = 0.0;
  }

  // for printing stats at the end
  float total_volume_fraction = 0.0, total_grid_mass = 0.0, total_mass_injected = 0.0;
  float total_energy_injected = 0.0, max_density_on_grid = 0.0, average_density_on_grid = 0.0;
  float total_metal_mass = 0.0;
  int   cells_this_grid = 0;

  // loop over cells, compute fractional volume, inject feedback
  for(int k = kc - r_int; k <= kc + r_int; k++){
    for(int j = jc - r_int; j <= jc + r_int; j++){
      for(int i = ic - r_int; i <= ic + r_int; i++){

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

        if (IndividualStarFollowStellarYields && cstar){
          for(int im = 0; im < StellarYieldsNumberOfSpecies+1; im++){
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
        if(TestProblemData.MultiMetals == 2 && IndividualStarFollowStellarYields){
          int field_num;
          this->IdentifyChemicalTracerSpeciesFieldsByNumber(field_num, 0); // gives metallicity field


          if (cstar){
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
            if (cstar){
              BaryonField[field_num][index] += injected_metal_mass[1 + im];
            } else { // keep same fraction if using artificial SN generator
              BaryonField[field_num][index] += delta_mass *
                                               BaryonField[field_num][index] / old_mass;
            }
          }

        } else{
          int field_num;
          this->IdentifyChemicalTracerSpeciesFieldsByNumber(field_num, 0); // gives metallicity field

          if (cstar){
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
  if (IndividualStarPrintSNStats && (!stellar_wind_mode) && (cstar)){
    // Column order: Grid ID, Particle ID, M_now, M_eject, Sphere Volume

    float average_metallicity;

    average_density_on_grid = total_grid_mass / (1.0 * cells_this_grid); // Sum Density / # cells

    /* convert to CGS */
    total_grid_mass         *= dx*dx*dx * MassUnits;
    total_mass_injected     *= dx*dx*dx * MassUnits;
    total_energy_injected   *= dx*dx*dx * EnergyUnits;
    total_metal_mass        *= dx*dx*dx * MassUnits;
    volume                  /= (LengthUnits * LengthUnits * LengthUnits);
    max_density_on_grid     *= DensityUnits;
    average_density_on_grid *= DensityUnits;
    m_eject                 *= dx*dx*dx * MassUnits;
    average_metallicity      = total_metal_mass / (total_mass_injected + total_grid_mass);

    /* compute Sedov-Taylor phase radius (R_dps) */

    printf("IndividualStarSNStats: %"ISYM" %"ISYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ISYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM" %"ESYM"\n",
          this->ID, cstar->ReturnID(), this->Time, cstar->ReturnMass(), cstar->ReturnBirthMass(), cstar->ReturnMetallicity(), m_eject,
          cells_this_grid, volume, total_volume_fraction, total_mass_injected, total_energy_injected,
          total_grid_mass, max_density_on_grid, average_density_on_grid, total_metal_mass, average_metallicity);
  }


  delete [] temperature;
  // done with spherical injection feedback

  return SUCCESS;
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
  if (IndividualStarFollowStellarYields && TestProblemData.MultiMetals == 2){
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
  if ( IndividualStarFollowStellarYields && TestProblemData.MultiMetals == 2){
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
  if(IndividualStarFollowStellarYields && TestProblemData.MultiMetals == 2){
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
      if(TestProblemData.MultiMetals == 2 && IndividualStarFollowStellarYields && m_ism > 0.0){
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

  if( IndividualStarFollowStellarYields && TestProblemData.MultiMetals == 2){

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
