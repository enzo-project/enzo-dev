/*
    Algorithm creates the stars for Mechanical star maker
        Formation follows Hopkins 2017 "How to model supernovae"

        07/2019: Azton Wells
 */
#include <stdio.h>
#include <math.h>
#include <mpi.h>
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


/* Creation Routine */
int grid::MechStars_Creation(grid* ParticleArray, float* Temperature,
        float *DMField, float* totalMetal, int level, float* CoolingTime,
        int MaximumNumberOfNewParticles, int* NumberOfParticlesSoFar)
{
    /* If limiting timesteps for SNe, return if not the right level */
    if (level < StarMakeLevel) return 0;
    /* 
        these flags are used to determine if we should set off a seed SNe 
        if (gridShouldFormStars && !notEnoughMetals) at the end, we'll seed 
        with P3 SNe if requested
    */
    bool gridShouldFormStars=false, notEnoughMetals = true;

    bool debug = true; // local debug flag; theres a lot of printing in here 


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


                    if (createStar){

                        /* Determine Mass of new particle */

                        float MassShouldForm = (shieldedFraction * BaryonField[DensNum][index]
                                        * MassUnits / (freeFallTime * TimeUnits) * Myr_s);
                        // Probability has the last word here
                        // FIRE-2 uses p = 1 - exp (-MassShouldForm*dt / M_gas_particle) to convert a whole particle to star particle
                        float p_form = 1.0 - exp(-1*MassShouldForm * this->dtFixed / (StarMakerMaximumFormationMass)); 
                        float random = (float) rand() / (float)(RAND_MAX);

                        printf("Expected Mass = %12.5e; Cell_mass = %12.5e; f_s = %12.5e; t_ff = %12.5e;  time-factor = %12.5e; nb = %12.5e; pform = %12.5e, rand = %12.5e; rand_max = %ld\n", 
                            MassShouldForm, BaryonField[DensNum][index] * MassUnits, shieldedFraction, freeFallTime*TimeUnits, 1.0/(freeFallTime*TimeUnits) * Myr_s,
                            BaryonField[DensNum][index]*DensityUnits/(mh/0.6), p_form, random, RAND_MAX);
                        if (MassShouldForm < 0){
                            printf("Negative formation mass: %f %f",shieldedFraction, freeFallTime);
                            continue;
                        }

                        /* limit new mass to 1/2 gas in cell, or StarMakerMaximumFormationMass */
                        float newMass = min(min(MassShouldForm/MassUnits, 0.5*BaryonField[DensNum][index]), StarMakerMaximumFormationMass/MassUnits);
                        
                        // only if mass is large enough

                        if ((newMass*MassUnits < StarMakerMinimumMass) || (random > p_form)){ // too small mass or probability is not in your favor, young one.
                                continue;
                        }
                        float totalDensity = (BaryonField[DensNum][index]+DMField[index])*DensityUnits;
                        dynamicalTime = pow(3.0*pi/32.0/GravConst/totalDensity, 0.5);
                        // fprintf(stdout, "DynamicalTime = %e\n", dynamicalTime);
                        ParticleArray->ParticleMass[nCreated] = newMass;
                        if (StarParticleRadiativeFeedback)
                            ParticleArray->ParticleAttribute[1][nCreated] = 25.0 * Myr_s * TimeUnits; // need 25 Myr lifetime for ray tracing feedback
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
                        /* average particle velocity over 125 cells to prevent runaway */
                        for (int kp = k-2; kp <= k+2; kp++)
                            for (int jp = j-2; jp <= j+2 ; jp++)
                                for (int ip = i-2; ip <= i + 2 ; ip++)
                                {
                                    int ind = ip + jp*GridDimension[0]+kp*GridDimension[0]+GridDimension[1];
                                    vX += BaryonField[Vel1Num][ind];
                                    vY += BaryonField[Vel2Num][ind];
                                    vZ += BaryonField[Vel3Num][ind];
                                }
                        float MaxVelocity = 250.*1.0e5/VelocityUnits;
                        ParticleArray->ParticleVelocity[0][nCreated] = 
                            (abs(vX/125.) > MaxVelocity)?(MaxVelocity*((vX > 0)?(1):(-1))):(vX/125.);
                        ParticleArray->ParticleVelocity[1][nCreated] = 
                            (abs(vY/125.) > MaxVelocity)?(MaxVelocity*((vY > 0)?(1):(-1))):(vY/125.);
                        ParticleArray->ParticleVelocity[2][nCreated] = 
                            (abs(vZ/125.) > MaxVelocity)?(MaxVelocity*((vZ > 0)?(1):(-1))):(vZ/125.);

                        /* give it position at center of host cell */

                        ParticleArray->ParticlePosition[0][nCreated] = CellLeftEdge[0][0]
                                                +(dx*(FLOAT(i)-0.5));
                        ParticleArray->ParticlePosition[1][nCreated] = CellLeftEdge[1][0]
                                                +(dx*(FLOAT(j)-0.5));
                        ParticleArray->ParticlePosition[2][nCreated] = CellLeftEdge[2][0]
                                                +(dx*(FLOAT(k)-0.5));
                        if (nCreated >= MaximumNumberOfNewParticles) return nCreated;
                        if (debug)
                            fprintf(stdout,"Created star: [%f] %e %e ::: %d %d %d ::: %e %f %e %e::: %f %f %f ::: %f %f %f \t::: %d %d ::: %d %d %d\n",
                                Time*TimeUnits/3.1557e13,
                                BaryonField[DensNum][index],
                                BaryonField[DensNum][index]*MassUnits, 
                                level, nCreated+1,
                                ParticleArray->ParticleType[nCreated],
                                ParticleArray->ParticleMass[nCreated]*MassUnits,
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
