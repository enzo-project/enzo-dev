/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR A TEST OF CHEMICAL EVOLUTION MODELS)
/
/  written by: Andrew Emerick
/  date:       Feb, 2016
/
/  PURPOSE: Plain, uniform gas grid. Star is deposited in individual_star_maker
/           this is a boring set up.
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "EnzoTiming.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "CosmologyParameters.h"
#include "phys_constants.h"
#include "IndividualStarProperties.h"
#include "StellarYieldsRoutines.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);


int grid::ChemicalEvolutionTestInitializeGrid(float GasDensity, float GasTemperature,
                                              float GasMetallicity){

  if (ProcessorNumber != MyProcessorNumber){
    return SUCCESS;
  }


  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
                              &TimeUnits, &VelocityUnits, Time) == FAIL){
      ENZO_FAIL("Error in GetUnits.");
  }

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                       Vel3Num, TENum) == FAIL){
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }

  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum, MetalNum;

  if(MultiSpecies){
    if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
                              HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL){
      ENZO_FAIL("Error in grid->IdentifySpeciesFields.");
    }
  }

  int MetallicityField = FALSE;
  if ((MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields)) != -1){
    MetallicityField = TRUE;
  } else {
    MetalNum = 0;
    printf("ChemicalEvolutionTest: Metallicity Field not found.\n");
  }

  int size = 1, i;

  for (int dim = 0; dim < GridRank; dim++){
    size *= GridDimension[dim];
  }

  if (TestProblemData.UseMetallicityField){
    for ( i = 0; i < size; i++){
      BaryonField[MetalNum][i] = GasMetallicity * GasDensity;
    }
  }


  if (MultiSpecies){
    // set set background to primordial and 100% ionized (only HII and HeIII)
    for( i = 0; i < size; i++){
      BaryonField[HIINum][i]   = GasDensity * TestProblemData.HydrogenFractionByMass
                                            * TestProblemData.HII_Fraction ;
      BaryonField[HeIINum][i]  = GasDensity * TestProblemData.HeII_Fraction *
                                          (1.0 - TestProblemData.HydrogenFractionByMass);
      BaryonField[HeIIINum][i] = GasDensity * TestProblemData.HeIII_Fraction
                                            * (1.0 - TestProblemData.HydrogenFractionByMass);
      BaryonField[HeINum][i]   = (1.0 - TestProblemData.HydrogenFractionByMass) * GasDensity -
                                 BaryonField[HeIINum][i] - BaryonField[HeIIINum][i];
      if(MultiSpecies > 1){
        BaryonField[HMNum][i]  = TestProblemData.HM_Fraction  * BaryonField[HIINum][i];
        BaryonField[H2INum][i] = TestProblemData.H2I_Fraction * GasDensity
                                                   * TestProblemData.HydrogenFractionByMass;
        BaryonField[H2IINum][i] = TestProblemData.H2II_Fraction * 2.0 * BaryonField[HIINum][i];
      }

      BaryonField[HINum][i] = TestProblemData.HydrogenFractionByMass *GasDensity
                                                - BaryonField[HIINum][i];

      if( MultiSpecies > 1){
        BaryonField[HINum][i] -= (BaryonField[HMNum][i] + BaryonField[H2IINum][i] +
                                 BaryonField[H2INum][i]);
      }

      // Electron density: sum up ionized species
      BaryonField[DeNum][i] = BaryonField[HIINum][i] + 0.25 * BaryonField[HeIINum][i] +
                                                      0.50 * BaryonField[HeIIINum][i];
      if (MultiSpecies > 1){
        BaryonField[DeNum][i] += 0.5*BaryonField[H2IINum][i] - BaryonField[HMNum][i];
      }

      if (MultiSpecies > 2){
        BaryonField[DINum ][i] = TestProblemData.DeuteriumToHydrogenRatio * BaryonField[HINum][i];
        BaryonField[DIINum][i] = TestProblemData.DeuteriumToHydrogenRatio * BaryonField[HIINum][i];
        BaryonField[HDINum][i] = 0.75 * TestProblemData.DeuteriumToHydrogenRatio * BaryonField[H2INum][i];
      }

    }// loop over cells
  } // Multispecies



  /* Loop over all requested stellar yields species and assign initial values */
  if (TestProblemData.MultiMetals == 2 || MultiMetals == 2){
    for( int sp = 0; sp < StellarYieldsNumberOfSpecies; sp++){
      if(StellarYieldsAtomicNumbers[sp] > 2){
        int   field_num;
        float fraction;

        fraction = TestProblemData.ChemicalTracerSpecies_Fractions[sp];

        IdentifyChemicalTracerSpeciesFieldsByNumber(field_num, StellarYieldsAtomicNumbers[sp]);

        for(i = 0; i < size; i++){ // assign initial values
          BaryonField[field_num][i] = fraction * GasDensity;
        }
      }
    } // loop over yields

  } // MULTI METALS

  /* now go through and initialize the particles */
  int MaximumNumberOfNewParticles = 10000;
  int NumberOfNewParticles = 0;
  this->AllocateNewParticles(MaximumNumberOfNewParticles);

  this->chemical_evolution_test_star_deposit(&MaximumNumberOfNewParticles,
                                             &NumberOfNewParticles, this->ParticleMass,
                                             this->ParticleType, this->ParticlePosition,
                                             this->ParticleVelocity, this->ParticleAttribute);

  if (NumberOfNewParticles > 0) {
    this->NumberOfParticles = NumberOfNewParticles;
  }
//  } else if (ChemicalEvolutionTestNumberOfStars > 0) {
//    ENZO_FAIL("Was not able to deposit stars in chemical evolution test\n");
//  } // end: if (NumberOfNewParticles > 0)



  return SUCCESS;
}


int grid::chemical_evolution_test_star_deposit(int *nmax, int *np, float *ParticleMass,
                                               int *ParticleType, FLOAT *ParticlePosition[],
                                               float *ParticleVelocity[], float *ParticleAttribute[]){


  if (ChemicalEvolutionTestNumberOfStars == 0) return SUCCESS;

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
      ParticleMass[count] = mass[i] * SolarMass / MassUnits / (dx*dx*dx);
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

      ParticleVelocity[0][count] = xvel[i]*kpc_cm / VelocityUnits;
      ParticleVelocity[1][count] = yvel[i]*kpc_cm / VelocityUnits;
      ParticleVelocity[2][count] = zvel[i]*kpc_cm / VelocityUnits;


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
    ParticleMass[0] = ChemicalEvolutionTestStarMass * SolarMass / MassUnits / (dx*dx*dx);
    ParticleType[0] = - PARTICLE_TYPE_INDIVIDUAL_STAR;
    ParticleAttribute[0][0] = this->Time;
     // allow user to set lifetime artificially
    if(ChemicalEvolutionTestStarLifetime > 0){
      ParticleAttribute[1][0] = ChemicalEvolutionTestStarLifetime * Myr_s / (TimeUnits);
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
    ParticleVelocity[0][0] = ChemicalEvolutionTestStarVelocity[0]*kpc_cm / VelocityUnits;
    ParticleVelocity[1][0] = ChemicalEvolutionTestStarVelocity[1]*kpc_cm / VelocityUnits;
    ParticleVelocity[2][0] = ChemicalEvolutionTestStarVelocity[2]*kpc_cm / VelocityUnits;

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
    printf("individual_star_maker: Formed star ChemicalEvolutionTest. M =  %"FSYM" and Z = %"FSYM". tau = %"ESYM"\n", ParticleMass[0]*(dx*dx*dx)*MassUnits/SolarMass, ParticleAttribute[2][0], ParticleAttribute[1][0]);
  }

  return SUCCESS;
}
