/*
    Algorithm creates the stars for Mechanical star maker
        Formation follows Hopkins 2017 "How to model supernovae"

        07/2019: Azton Wells
 */
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include <vector>
#include <limits.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "CosmologyParameters.h"
#include "StarParticleData.h"
#include "phys_constants.h"
//Prototypes
    int checkCreationCriteria(float* Density,
                        float* Metals,
                        float* Temperature,float* DMField,
                        float* Vel1, float* Vel2, float* Vel3, float* TotE,
                        float* CoolingTime, int* GridDim,
                        float* shieldedFraction, float* freeFallTime,
                        float* dynamicalTime, int i, int j, int k,
                        float Time, float* RefinementField, float CellWidth,
                        bool* gridShouldFormStars, bool* notEnoughMetals, 
                        int continuingFormation, int* seedIndex);
    int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, float *MassUnits, float Time);
    int FindField(int field, int farray[], int numfields);

    void mt_init(unsigned_int seed);
    unsigned_long_int mt_random();


/* Creation Routine */
int grid::MechStars_Creation(grid* ParticleArray, float* Temperature,
        float *DMField, float* totalMetal, int level, float* CoolingTime,
        int MaximumNumberOfNewParticles, int* NumberOfParticlesSoFar)
{
    bool use_F2 = true; // use FIRE-2 methods of self-shielded fractionn and virial parameter
    float conversion_fraction = StarMakerMassEfficiency; // max fraction of baryon mass that gets converted to stars
                                   // on successful formation check
    if (MyProcessorNumber != ProcessorNumber)
        return 0;
    /* If limiting timesteps for SNe, return if not the right level */
    if (level < StarMakeLevel) return 0;
    /* 
        these flags are used to determine if we should set off a seed SNe 
        if (gridShouldFormStars && !notEnoughMetals) at the end, we'll seed 
        with P3 SNe if requested
    */
    bool gridShouldFormStars=false, notEnoughMetals = true;

    bool debug = true; // local debug flag; theres a lot of printing in here 
    mt_init(clock());


    //get field numbers
    int DensNum, GENum, Vel1Num, Vel2Num,Vel3Num, TENum;
    if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
        fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
        return FAIL;
    }

    int MetallicityField = FALSE, MetalNum, MetalIaNum, MetalIINum, SNColourNum;
    if ((MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields))
             != -1)
    {
        MetallicityField = TRUE;
               float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1,
                        TimeUnits = 1, VelocityUnits = 1, MassUnits = 1;
                if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
                        &TimeUnits, &VelocityUnits, &MassUnits, this->Time) == FAIL) {
                    fprintf(stderr, "Error in GetUnits.\n");
                    return FAIL;
                    }
    }
    else
        MetalNum = 0;
    if (!MetalNum){
        fprintf(stderr, "Can only use mechanical star maker routines with metallicity enabled!");
        return FAIL;
    }
    MetalIaNum = FindField(MetalSNIaDensity, FieldType, NumberOfBaryonFields);
    MetalIINum = FindField(MetalSNIIDensity, FieldType, NumberOfBaryonFields);
    SNColourNum = FindField(SNColour, FieldType, NumberOfBaryonFields);
    int size =1;
    for (int dim = 0; dim < GridRank; dim ++)
        size *= GridDimension[dim];


    float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1,
                TimeUnits = 1, VelocityUnits = 1, MassUnits = 1;
    if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
                &TimeUnits, &VelocityUnits, &MassUnits, this->Time) == FAIL) {
        fprintf(stderr, "Error in GetUnits.\n");
    return FAIL;
    }

    /*      
        ReDefine MassUnits so that code * mass = physical (Msun)
            and physical(Msun) / MassUnits = code */
    MassUnits = DensityUnits*pow(LengthUnits*CellWidth[0][0], 3)/SolarMass;

    /* Index of last cell that was capable of star formation but has no metals */
    int *seedIndex = new int [3];

    FLOAT dx = CellWidth[0][0];
    FLOAT cell_vol = dx*dx*dx;
    int GZ = int(NumberOfGhostZones);
    int nCreated = *NumberOfParticlesSoFar;
    // first, loop over particles to find the stars on this grid;
    // if a cell qualifies for star formation, we want to add to 
    // existing stars before forming new ones!
    // std::vector<int> ptcl_inds; // vector or xyz indices of particle positions
    // for (int pIndex = 0; pIndex < NumberOfParticles; pIndex++){
    //     if (ParticleType[pIndex] ==2){
    //        ptcl_inds.push_back(ParticleNumber[pIndex]) ;
        
    //     }
    // }
    // int nPriorStars = ptcl_inds.size();
    //fprintf(stdout, "Starting creation with %d prior particles\n",nCreated);
    for (int k = GZ; k< GridDimension[2]-GZ; k++){
        for(int j = GZ; j < GridDimension[1]-GZ; j++){
            for (int i = GZ; i < GridDimension[0]-GZ; i++){
                /*
                Particle creation has several criteria:
                    0. is this the finest level for this cell
                    1. Density > overdensity
                    2. is flow converging by finite differences
                    3. is cooling time < dynamical time
                    4. is gas mass > jean critical mass
                        Fire-2 motivated checks:
                    5. is gas self shielded by Krumholz & Gnedin 2007 criteria
                    6. is virial parameter < 1

                 */
                int index = i+ j*GridDimension[0]+k*GridDimension[0]*GridDimension[1];
                if (BaryonField[NumberOfBaryonFields][index] != 0.0) continue;
                float shieldedFraction = 0;
                float freeFallTime = 0;
                float dynamicalTime = 0;
                float Time = this->Time;
                int createStar = false;
        

                    createStar = checkCreationCriteria(BaryonField[DensNum],
                        totalMetal, Temperature, DMField,
                        BaryonField[Vel1Num], BaryonField[Vel2Num],
                        BaryonField[Vel3Num], BaryonField[TENum],
                        CoolingTime, GridDimension, &shieldedFraction,
                        &freeFallTime, &dynamicalTime, i,j,k,Time,
                        BaryonField[NumberOfBaryonFields], CellWidth[0][0],
                        &gridShouldFormStars, &notEnoughMetals, 0, seedIndex);


                    if (createStar!=FAIL){

                        /* Determine Mass of new particle 
                            WARNING: this removes the mass of the formed particle from the 
                            host cell.  If your simulation has very small (>15 Msun) baryon mass
                            per cell, it will break your sims! - AIW
                        */
                        float divisor = max( 1.0, freeFallTime * TimeUnits / Myr_s);
                        float MaximumStarMass = StarMakerMaximumFormationMass;
                        if (MaximumStarMass < 0)
                            MaximumStarMass = conversion_fraction * BaryonField[DensNum][index] * MassUnits;
                        float MassShouldForm = 0.0;
                        // if (use_F2)
                            MassShouldForm = min(shieldedFraction * BaryonField[DensNum][index]
                                        * MassUnits / divisor, conversion_fraction * BaryonField[DensNum][index] * MassUnits / divisor);
                        // else
                        //     MassShouldForm = (MaximumStarMass/divisor);
                        
                        // Probability has the last word
                        // FIRE-2 uses p = 1 - exp (-MassShouldForm*dt / M_gas_particle) to convert a whole particle to star particle
                        //  We convert a fixed portion of the baryon mass (or the calculated amount)
                        float p_form = 1.0 - exp(-1*MassShouldForm * this->dtFixed 
						                / (MaximumStarMass)); 
                        
                        float random = float(mt_random())/float(UINT_MAX);
                        
                        if (debug && MassShouldForm > 0)
			                printf("[t=%f, np = %d] Expected Mass = %12.5e; Cell_mass = %12.5e; f_s = %12.5e; t_ff = %12.5e;  time-factor = %12.5e; nb = %12.5e; tdyn = %3.3f; t_cool = %3.3f; pform = %12.5e, rand = %12.5e; rand_max = %ld\n", Time*TimeUnits/Myr_s, 0,
                            MassShouldForm, BaryonField[DensNum][index] * MassUnits, 
                            shieldedFraction, freeFallTime*TimeUnits / Myr_s, 1.0/(freeFallTime*TimeUnits) * Myr_s,
                            BaryonField[DensNum][index]*DensityUnits/(mh/0.6), dynamicalTime / Myr_s, 
                            CoolingTime[index]*TimeUnits/Myr_s, p_form, random, RAND_MAX);
                        
                        if (MassShouldForm < 0 && !use_F2){
                             printf("Negative formation mass: %f %f\n",shieldedFraction, freeFallTime);
                             continue;
                        }

                        /* New star is MassShouldForm up to `conversion_fraction` * baryon mass of the cell, but at least 15 msun */
                        float newMass = min(MassShouldForm/MassUnits, MaximumStarMass/MassUnits); 

                        if ((newMass*MassUnits < StarMakerMinimumMass) /* too small */
                                || (random > p_form) /* too unlikely */
                                || (newMass > BaryonField[DensNum][index])) /* too big compared to cell 
                                                                            (make sure min/max mass is within reasonable range for gas mass)*/
                                { 
                                    continue;
                                }
 

                        ParticleArray->ParticleMass[nCreated] = newMass;
                        if (StarParticleRadiativeFeedback)
                            ParticleArray->ParticleAttribute[1][nCreated] = 25.0 * Myr_s / TimeUnits; // need 25 Myr lifetime for ray tracing feedback
                        else
                            ParticleArray->ParticleAttribute[1][nCreated] = 0.0; // Tracking SNE in TDP field dynamicalTime/TimeUnits;
                        ParticleArray->ParticleAttribute[0][nCreated] = Time;

                        ParticleArray->ParticleAttribute[2][nCreated] = totalMetal[index]
                                    /BaryonField[DensNum][index];
                        if (StarMakerTypeIaSNe)
                            ParticleArray->ParticleAttribute[3][nCreated] = BaryonField[MetalIaNum][index];

                        ParticleArray->ParticleType[nCreated] = PARTICLE_TYPE_STAR;
                        BaryonField[DensNum][index] -= newMass;
                        BaryonField[MetalNum][index] -= newMass*totalMetal[index]/BaryonField[DensNum][index];
                        if (SNColourNum > 0)
                            BaryonField[SNColourNum][index] -= newMass/BaryonField[DensNum][index]*BaryonField[SNColourNum][index]/BaryonField[DensNum][index];


        // type IA metal field isnt here right now
                        // if (StarMakerTypeIASNe)
                        //     ParticleArray->ParticleAttribute[3][nCreated] = BaryonField[MetalIANum];
                        float vX = 0.0;
                        float vY = 0.0;
                        float vZ = 0.0;
                        if (HydroMethod != 0)
                            fprintf(stderr,"Mechanical star maker not tested for anything except HydroMethod = 0\n");
                        /* average particle velocity over many cells to prevent runaway */
                        float cnter = 0;
                        float msum = 0;
                        for (int kp = max(0,k-3); kp <= min(k+3, GridDimension[2]); kp++)
                            for (int jp = max(0,j-3); jp <= min(j+3, GridDimension[1]) ; jp++)
                                for (int ip = max(0,i-3); ip <= min(i + 3, GridDimension[0]) ; ip++)
                                {
                                    cnter ++;
                                    int ind = ip + jp*GridDimension[0]+kp*GridDimension[0]+GridDimension[1];
                                    vX += BaryonField[Vel1Num][ind]*BaryonField[DensNum][ind];
                                    vY += BaryonField[Vel2Num][ind]*BaryonField[DensNum][ind];
                                    vZ += BaryonField[Vel3Num][ind]*BaryonField[DensNum][ind];
                                    msum += BaryonField[DensNum][ind];
                                }
                        vX = vX / msum;
                        vY = vY / msum;
                        vZ = vZ / msum;
                        float MaxVelocity = 150.*1.0e5/VelocityUnits;
                        ParticleArray->ParticleVelocity[0][nCreated] = 
                            (abs(vX) > MaxVelocity)?(MaxVelocity*((vX > 0)?(1):(-1))):(vX);
                        ParticleArray->ParticleVelocity[1][nCreated] = 
                            (abs(vY) > MaxVelocity)?(MaxVelocity*((vY > 0)?(1):(-1))):(vY);
                        ParticleArray->ParticleVelocity[2][nCreated] = 
                            (abs(vZ) > MaxVelocity)?(MaxVelocity*((vZ > 0)?(1):(-1))):(vZ);

                        /* give it position at center of host cell */

                        ParticleArray->ParticlePosition[0][nCreated] = CellLeftEdge[0][0]
                                                +(dx*(FLOAT(i)-0.5));
                        ParticleArray->ParticlePosition[1][nCreated] = CellLeftEdge[1][0]
                                                +(dx*(FLOAT(j)-0.5));
                        ParticleArray->ParticlePosition[2][nCreated] = CellLeftEdge[2][0]
                                                +(dx*(FLOAT(k)-0.5));
                        if (nCreated >= MaximumNumberOfNewParticles) return nCreated;
                        // if (debug)
                            fprintf(stdout, "\t\tCreated star: [%f Myr] M_cell = %e T_cell = %e\n\t\t\tL = %d  N = %d Type = %d\n\t\t\tM_* = %e Tdyn = %e Attr1 = %f Attr2 = %e Attr3 = %e\n\t\t\tPosition = [%f %f %f]\n\t\t\tvelocity = [%f %f %f]\n\t\t\tgrid_index = %d cell_index = %d\n\t\t\ti = %d j = %d k = %d\n",
                                Time*TimeUnits/3.1557e13,
                                BaryonField[DensNum][index]*MassUnits, 
                                Temperature[index],
                                level, nCreated+1,
                                ParticleArray->ParticleType[nCreated],
                                ParticleArray->ParticleMass[nCreated]*MassUnits,
                                dynamicalTime / yr_s,
                                ParticleArray->ParticleAttribute[0][nCreated],
                                ParticleArray->ParticleAttribute[1][nCreated],
                                ParticleArray->ParticleAttribute[2][nCreated],
                                ParticleArray->ParticlePosition[0][nCreated],
                                ParticleArray->ParticlePosition[1][nCreated],
                                ParticleArray->ParticlePosition[2][nCreated],
                                ParticleArray->ParticleVelocity[0][nCreated]*VelocityUnits/1e5,
                                ParticleArray->ParticleVelocity[1][nCreated]*VelocityUnits/1e5,
                                ParticleArray->ParticleVelocity[2][nCreated]*VelocityUnits/1e5,
                                index, GridDimension[0]*GridDimension[2]*GridDimension[3], i,j,k);
                        nCreated ++;
                }
            }//end for k
        }//end for j
    } // end for i
    /*
        If the grid has met star formation criteria, but stars are not formed due to lack of metals,
        we set off a P3 SN event at the last grid cell that could have hosted star formation. Can and will set
        off multiple events per grid cell if the same cell meets criteria on the next iteration!
     */
    if (gridShouldFormStars && MechStarsSeedField && (nCreated == 0)){
        // set off a p3 supernova at at the last cell that could 
        // host star formation in the grid if the
        // grid can host star formation but has no appreciable metals
        fprintf(stdout, "\n\n\n[%d] %d %d %d Creating seed field!\n\n\n\n", 
                ID,seedIndex[0], seedIndex[1], seedIndex[2]) ;
        MechStars_SeedSupernova(&totalMetal[0], Temperature, seedIndex);
        
    }
    //if (nCreated > 0 && debug){
    //  fprintf(stdout, "Created %d star particles\n",nCreated);
    //    }
    delete [] seedIndex;
    // delete [] totalMetal;
    *NumberOfParticlesSoFar = nCreated;
    return nCreated;
}
