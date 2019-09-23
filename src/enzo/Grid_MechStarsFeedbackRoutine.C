/*
    Routine to determine feedback quantities and couple them to the
    Grid.  Coupling follows the implementation of Hopkins 2017 with
    modification for Enzo's fixed grid

    07/2019: Azton Wells
 */
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include "ErrorExceptions.h"
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


    int determineSN(float age, int* nSNII, int* nSNIA, float massMsun,
                    float TimeUnits, float dtFixed);
    int determineWinds(float age, float* eWinds, float* zWinds, float* mWinds,
                        float massMsun, float zZsun, float TimeUnits, float dtFixed);
int checkCreationCriteria(float* Density, float* Metals,
                        float* Temperature,float* DMField,
                        float* Vel1, float* Vel2, float* Vel3, 
                        float* CoolingTime, int* GridDim,
                        float* shieldedFraction, float* freeFallTime, 
                        float* dynamicalTime, int i, int j, int k, 
                        float Time, float* RefinementField, float CellWidth,
                        bool* gridShouldFormStars, bool* notEnoughMetals, 
                        int continuingFormation, int* seedIndex);
    int FindField(int field, int farray[], int numfields);
    int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, float *MassUnits, FLOAT Time);





int grid::MechStars_FeedbackRoutine(int level, float* mu_field, 
            float* Temperature, float* CoolingTime, float* DMField)
{

    //fprintf(stdout,"IN FEEDBACK ROUTINE\n  %d   %d   %d\n",
        //SingleSN, StellarWinds, UnrestrictedSN);
    float stretchFactor = 1.0;//1/sqrt(2) to cover cell diagonal
    /* Get units to use */
    bool SingleWinds = true; // flag to consolidate wind feedback into one event centered on most massive cell in grid
            // I wouldn't recommend this for unigrid runs...
    bool debug = true;
    float startFB = MPI_Wtime();
    int dim, i, j, k, index, size, field, GhostZones = NumberOfGhostZones;
    int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;

    /* Compute size (in floats) of the current grid. */

    size = 1;
    for (dim = 0; dim < GridRank; dim++)
        size *= GridDimension[dim];
    int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
                        HMNum, H2INum, H2IINum, DINum, DIINum, HDINum;
    /* Find fields: density, total energy, velocity1-3. */

    if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                        Vel3Num, TENum) == FAIL) {
        fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
        return FAIL;
    }
    /* Set the units */
    float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1,
            TimeUnits = 1, VelocityUnits = 1, MassUnits = 1;
    if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
            &TimeUnits, &VelocityUnits, &MassUnits, this->Time) == FAIL) {
        fprintf(stderr, "Error in GetUnits.\n");
        return FAIL;
    }
    double dx = CellWidth[0][0];
    MassUnits = DensityUnits*pow(LengthUnits*dx, 3)/SolarMass;
    /*
        get metallicity field and set flag; assumed true thoughout feedback
        since so many quantities are metallicity dependent
     */
    int MetallicityField = FALSE, MetalNum, SNColourNum=-1;
    if ((MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields))
        != -1)
        MetallicityField = TRUE;
    else
        MetalNum = 0;
    if (MechStarsSeedField) 
        SNColourNum = FindField(SNColour, FieldType, NumberOfBaryonFields);
    float* totalMetal = new float [size];
    for (int i = 0; i < size; i++){
        totalMetal[i] = BaryonField[MetalNum][i];
        if (MechStarsSeedField)
            totalMetal[i] += BaryonField[SNColourNum][i];
    }
    int numSN = 0; // counter of events
    int c = 0; // counter of particles
    float maxD = 0.0;
    int maxI=0, maxJ=0, maxK=0;
    int maxindex=0;
    /* Begin Iteration of all particles */
    // printf("\nIterating all particles  ");
    for (int pIndex=0; pIndex < NumberOfParticles; pIndex++){
        // if (ParticleType[pIndex] != 1 && debug)
        // fprintf(stdout,"PARTICLE: %d %d %e %f\n", ParticleType[pIndex],
            // ParticleNumber[pIndex],
            // ParticleMass[pIndex],
            // ParticleAttribute[0][pIndex]);
        /* Selection criteria */

        if (ParticleType[pIndex] == PARTICLE_TYPE_STAR
                && ParticleMass[pIndex] > 0.0
                && ParticleAttribute[0][pIndex] > 0.0){
            c++;
            // if (StarMakerAgeCutoff)
            //     if ((Time-ParticleAttribute[0][pIndex])
            //         *TimeUnits/(150*3.1557e7) > 150)
            //         continue;

            /* get index of cell hosting particle */
            float xp = ParticlePosition[0][pIndex];
            float yp = ParticlePosition[1][pIndex];
            float zp = ParticlePosition[2][pIndex];

            int ip = (xp-CellLeftEdge[0][0]-0.5*dx)/dx;
            int jp = (yp-CellLeftEdge[1][0]-0.5*dx)/dx;
            int kp = (zp-CellLeftEdge[2][0]-0.5*dx)/dx;



            /* error check particle position; Cant be on the border or outside grid
                If on border, reposition to within grid for CIC deposit */
            FLOAT age = (Time-ParticleAttribute[0][pIndex])*TimeUnits/3.1557e13;// Myr

            float gridDx = GridDimension[0]*dx;
            float gridDy = GridDimension[1]*dx;
            float gridDz = GridDimension[2]*dx;
          /* Keep particle 2 cells from edge since we cant transfer to
                neighboring grids */
            float borderDx = (stretchFactor+1)*dx;
            if (xp > CellLeftEdge[0][0]+gridDx
                || xp < CellLeftEdge[0][0]
                || yp > CellLeftEdge[1][0]+gridDy
                || yp < CellLeftEdge[1][0]
                || zp > CellLeftEdge[2][0]+gridDz
                || zp < CellLeftEdge[2][0]){
                fprintf(stderr, "Particle %d out of grid!\nage: %d, pos: %f, %f, %f GridEdge: %f %f %f", pIndex,
                    age, xp, yp, zp, CellLeftEdge[0][0], CellLeftEdge[1][0], CellLeftEdge[2][0]
                    );
                EnzoFatalException("Star Maker Mechanical: particle out of grid!\n");
                }
            int shifted = 0;

            if (xp < CellLeftEdge[0][0]+borderDx){
                xp = CellLeftEdge[0][0]+borderDx+0.5*dx;
                shifted ++;
            }
            if (xp > CellLeftEdge[0][0]+gridDx-borderDx){
                xp = CellLeftEdge[0][0]+gridDx-borderDx-0.5*dx;
                shifted = 1;
            }
            if (yp < CellLeftEdge[1][0]+borderDx){
                yp = CellLeftEdge[1][0]+borderDx+0.5*dx;
                shifted = 1;
            }
            if (yp > CellLeftEdge[1][0]+gridDx-borderDx){
                yp = CellLeftEdge[1][0]+gridDx-borderDx-0.5*dx;
                shifted = 1;
            }
            if (zp < CellLeftEdge[2][0]+borderDx){
                zp = CellLeftEdge[2][0]+borderDx+0.5*dx;
                shifted = 1;
            }
            if (zp > CellLeftEdge[2][0]+gridDx-borderDx){
                zp = CellLeftEdge[2][0]+gridDx-borderDx-0.5*dx;
                shifted = 1;
            }
            if (shifted > 0){
            if (debug)
                    fprintf(stderr, "Particle position shifted away from edge: %e: %f %f %f\n%f %f %f\n",
                            age,xp, yp, zp, CellLeftEdge[0][0]+borderDx, CellLeftEdge[1][0]+borderDx, CellLeftEdge[2][0]+borderDx);
            int ip = (xp-CellLeftEdge[0][0]-0.5*dx)/dx;
            int jp = (yp-CellLeftEdge[1][0]-0.5*dx)/dx;
            int kp = (zp-CellLeftEdge[2][0]-0.5*dx)/dx;
            }
            /* Check for continual formation.  Continually forming new mass allows the 
                star particle count to stay lower, ultimately reducing runtime by having 
                fewer particles to iterate. 
            */
            index = ip+jp*GridDimension[0]+kp*GridDimension[0]*GridDimension[1];

            float shieldedFraction = 0, dynamicalTime = 0, freeFallTime = 0;
            bool gridShouldFormStars = true, notEnoughMetals=false;
            float zFraction = BaryonField[MetalNum][index];
            if (MechStarsSeedField){
                    zFraction += BaryonField[SNColourNum][index];
            }
            if (ParticleMass[pIndex]*MassUnits < StarMakerMaximumMass){
                int createStar = checkCreationCriteria(BaryonField[DensNum],
                        &zFraction, Temperature, DMField,
                        BaryonField[Vel1Num], BaryonField[Vel2Num],
                        BaryonField[Vel3Num],
                        CoolingTime, GridDimension, &shieldedFraction,
                        &freeFallTime, &dynamicalTime, ip,jp,kp,Time,
                        BaryonField[NumberOfBaryonFields], CellWidth[0][0],
                        &gridShouldFormStars, &notEnoughMetals, 1, NULL);
                if(createStar){
                    float MassShouldForm =min((shieldedFraction * BaryonField[DensNum][index]
                                        * MassUnits / freeFallTime * this->dtFixed*TimeUnits/3.1557e13),
                                        0.5*BaryonField[DensNum][index]*MassUnits);
                    printf("Adding new mass %e\n",MassShouldForm);
                    /* Dont allow negative mass, or taking all gas in cell */
                    if (MassShouldForm < 0 )
                        MassShouldForm = 0;
                    if (MassShouldForm > 0.5*BaryonField[DensNum][index]*MassUnits)
                        MassShouldForm = 0.5*BaryonField[DensNum][index]*MassUnits;

                    // Set units and modify particle
                    MassShouldForm /= MassUnits;
                    if (MassShouldForm > 0){
                        float delta = MassShouldForm/(ParticleMass[pIndex]+MassShouldForm);

                        /* modify metallicity */
                        zFraction /= BaryonField[DensNum][index];
                        // update mass-weighted metallicity fraction of star particle
                        ParticleAttribute[2][pIndex] = (ParticleAttribute[2][pIndex]*(1.-delta)+zFraction*delta); 
                        // update mass-weighted age of star particle
                        if (age > 3.5) // only update if particle is old enough for SNe
                            ParticleAttribute[0][pIndex] = (ParticleAttribute[0][pIndex]*(1.-delta)+Time*delta);
                        /* Add new formation mass to particle */
                        ParticleMass[pIndex] += MassShouldForm;   
                        printf("[%f] added new mass %e + %e = %e newZ = %f newAge = %f\n", 
                            Time*TimeUnits/3.1557e13, (ParticleMass[pIndex]-MassShouldForm)*MassUnits, 
                            MassShouldForm*MassUnits, ParticleMass[pIndex]*MassUnits,
                            ParticleAttribute[2][pIndex],(Time- ParticleAttribute[0][pIndex])*TimeUnits/3.1557e13);
                        /* Take formed mass out of grid cell */
                        BaryonField[DensNum][index] -= MassShouldForm;
                        /* Take metals out of host cell too! */
                        BaryonField[MetalNum][index] -= BaryonField[MetalNum][index]/BaryonField[DensNum][index]*MassShouldForm;
                        if (MechStarsSeedField)
                            BaryonField[SNColourNum][index] -= BaryonField[SNColourNum][index]/BaryonField[DensNum][index]*MassShouldForm;             
                    }
                }
            }
            
            /* Start actual feedback: Supernova calculations */
            int nSNII = 0;
            int nSNIA = 0;
            float SNMassEjected = 0, SNMetalEjected = 0;
            //fprintf(stdout, "Checking particle:  %f    %e \n", age, ParticleMass[pIndex]*MassUnits);
            /* determine how many supernova events */
            if (SingleSN)
            {
                // fprintf(stdout,"Checking for SN age = %f\n", age);
                determineSN(age, &nSNII, &nSNIA, ParticleMass[pIndex]*MassUnits,
                            TimeUnits, dtFixed);
                numSN += nSNII+nSNIA;
		        if ((nSNII > 0 || nSNIA > 0) && debug)
                    fprintf(stdout,"SUPERNOVAE!!!! %d %d level = %d\n", nSNII, nSNIA, level);
                if (nSNII > 0 || nSNIA > 0){
                    /* set feedback qtys based on number and types of events */
                        /* 1e51 erg per sn */
                    ParticleAttribute[1][pIndex] += nSNII + nSNIA;
                    float energySN = (nSNII + nSNIA)*1e51;

                        /*10.5 Msun ejecta for type II and IA*/
                    SNMassEjected = (nSNII+nSNIA)*10.5;
                        /* Metal yeilds from starburst 99 */
                    float starMetal = (ParticleAttribute[2][pIndex]/0.02); //determines metal content of SNeII
                    // if (StarMakerTypeIISNeMetalField)
                    //     starMetal += ParticleAttribute[4][pIndex]/0.02;
                    /* Couple these in the deposit routine */
                    // SNMetalEjected = nSNII*(1.91+0.0479*max(zZsun, 1.65));
                    // SNMetalEjected += nSNIA*(1.4); // this metal should get coupled to SNIA field if its being used
                    MechStars_DepositFeedback(energySN, SNMassEjected, SNMetalEjected, totalMetal,
                                &ParticleVelocity[0][pIndex], &ParticleVelocity[1][pIndex], &ParticleVelocity[2][pIndex],
                                &ParticlePosition[0][pIndex], &ParticlePosition[1][pIndex], &ParticlePosition[2][pIndex],
                                ip, jp, kp, size, mu_field, 0, nSNII, nSNIA, starMetal, 0);
                    ParticleAttribute[1][pIndex] += nSNII+nSNIA;
                }
            }

            /* Do the same for winds. Cooling Radius is very small,
            So almost no energy is coupled, but some mass may be. */
            float windEnergy=0, windMass=0, windMetals=0;

            /*
                Ignore very old stars, veryvery young stars, and ones whose mass is depleted
             */
            if (StellarWinds && age > 0.001 && ParticleMass[pIndex]*MassUnits > 1)
            {
                // printf("Checking Winds\n");
                float zZsun = min(ParticleAttribute[2][pIndex]/0.02, MechStarsCriticalMetallicity);

                determineWinds(age, &windEnergy, &windMass, &windMetals,
                                ParticleMass[pIndex]*MassUnits, zZsun,
                                TimeUnits, dtFixed);
                if (windMass > 10) fprintf(stdout,"Really High Wind Mass!!\n");
                if (windEnergy > 1e5){
                    printf("Winds: M = %e E=%e\n", windMass, windEnergy);
                    MechStars_DepositFeedback(windEnergy, windMass, windMetals, totalMetal,
                                        &ParticleVelocity[0][pIndex], &ParticleVelocity[1][pIndex], &ParticleVelocity[2][pIndex],
                                        &ParticlePosition[0][pIndex], &ParticlePosition[1][pIndex], &ParticlePosition[2][pIndex],
                                        ip, jp, kp, size, mu_field, 1, 0, 0, 0.0, 0);
                }

            }
            if (windMass > 0.0 || SNMassEjected > 0){
                //if (debug) printf("Subtracting off mass %e\n",(windMass+SNMassEjected));
                ParticleMass[pIndex] -= (windMass+SNMassEjected)/MassUnits;
            }
            // printf("Post-feedback MP = %e\n", ParticleMass[pIndex]*MassUnits);
        }
    }// end iteration over particles
    if (c > 0){
        fprintf(stdout, "Ptcl Number = %d Events = %d FeedbackTime = %e Size = %d\n",
            c, numSN, MPI_Wtime()-startFB, GridDimension[0]*GridDimension[1]*GridDimension[2]);
    }
    /* to avoid iterating deposition over 10k particles, do ONE winds feedback that is summed all wind feedback in the region, 
        centered on the most massive cell in the grid*/
    delete [] totalMetal;
    return SUCCESS;
}
