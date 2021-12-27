/*
    Routine to determine feedback quantities and couple them to the
    Grid.  Coupling follows the implementation of Hopkins 2017 with
    modification for Enzo's fixed grid

    07/2019: Azton Wells
 */
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <signal.h>
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

int determineSN(float age, int *nSNII, int *nSNIA, float massMsun,
                float TimeUnits, float dtFixed);
int determineWinds(float age, float *eWinds, float *zWinds, float *mWinds,
                   float massMsun, float zZsun, float TimeUnits, float dtFixed);
int checkCreationCriteria(float *Density, float *Metals,
                          float *Temperature, float *DMField,
                          float *Vel1, float *Vel2, float *Vel3, float* TotE,
                          float *CoolingTime, int *GridDim,
                          float *shieldedFraction, float *freeFallTime,
                          float *dynamicalTime, int i, int j, int k,
                          float Time, float *RefinementField, float CellWidth,
                          bool *gridShouldFormStars, bool *notEnoughMetals,
                          int continuingFormation, int *seedIndex);
int FindField(int field, int farray[], int numfields);
int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, float *MassUnits, FLOAT Time);
int MechStars_depositEmissivityField(int index, float cellwidth,
                                     float *emissivity0, float age, float mass,
                                     float TimeUnits, float dt);

int grid::MechStars_FeedbackRoutine(int level, float *mu_field, float *totalMetal,
                                    float *Temperature, float *CoolingTime, float *DMField)
{

    //fprintf(stdout,"IN FEEDBACK ROUTINE\n  %d   %d   %d\n",
    //SingleSN, StellarWinds, UnrestrictedSN);
    if (MyProcessorNumber != ProcessorNumber)
        return 0;

    float Zsolar = 0.02;
    float stretchFactor = 2.0; // distance in cells required from grid edge for FB to function.
    bool debug = false;
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
                                         Vel3Num, TENum) == FAIL)
    {
        fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
        return FAIL;
    }
    /* Set the units */
    float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1,
          TimeUnits = 1, VelocityUnits = 1, MassUnits = 1;
    if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
                 &TimeUnits, &VelocityUnits, &MassUnits, this->Time) == FAIL)
    {
        fprintf(stderr, "Error in GetUnits.\n");
        return FAIL;
    }
    FLOAT dx = CellWidth[0][0];
    MassUnits = DensityUnits * pow(LengthUnits * dx, 3) / SolarMass;
    /*
        get metallicity field and set flag; assumed true thoughout feedback
        since so many quantities are metallicity dependent
     */
    int MetallicityField = FALSE, MetalNum, SNColourNum = -1;
    if ((MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields)) != -1)
        MetallicityField = TRUE;
    else
        MetalNum = 0;
    if (MechStarsSeedField)
        SNColourNum = FindField(SNColour, FieldType, NumberOfBaryonFields);

    /* Find and set emissivity field if being used */
    int EmisNum = -1;
    if (StarMakerEmissivityField)
    {
        EmisNum = FindField(Emissivity0, FieldType, NumberOfBaryonFields);
    }
    if (EmisNum > 0)
    {
        /* Zero the emissivity first, as we dont want to to accumulate */
        for (int i = 0; i < size; ++i)
            BaryonField[EmisNum][i] = 0.0;
    }

    int numSN = 0; // counter of events
    int c = 0;     // counter of particles
    // printf("\nIterating all particles  ");
    for (int pIndex = 0; pIndex < NumberOfParticles; pIndex++)
    {

        if (ParticleType[pIndex] == PARTICLE_TYPE_STAR && ParticleMass[pIndex] > 00.0 && ParticleAttribute[0][pIndex] > 0.0)
        {
            float age = (Time - ParticleAttribute[0][pIndex]) * TimeUnits / Myr_s; // Myr
            c++;

            /* get index of cell hosting particle */
            FLOAT xp = ParticlePosition[0][pIndex];
            FLOAT yp = ParticlePosition[1][pIndex];
            FLOAT zp = ParticlePosition[2][pIndex];


            /* error check particle position; Cant be on the border or outside grid
                If on border, reposition to within grid for CIC deposit */

            FLOAT gridDx = GridDimension[0] * dx;
            FLOAT gridDy = GridDimension[1] * dx;
            FLOAT gridDz = GridDimension[2] * dx;
            /* Keep particle 2 cells from edge since we cant transfer to
                neighboring grids */
            FLOAT borderDx = (stretchFactor)*dx;
                if (xp > GridRightEdge[0] 
                    || xp < GridLeftEdge[0] 
                    || yp > GridRightEdge[1]
                    || yp < GridLeftEdge[1] 
                    || zp > GridRightEdge[2] 
                    || zp < GridLeftEdge[2])
                    {
                    if (debug){
                        fprintf(stderr, "[%d--%e] Particle %" ISYM "  with type %" ISYM " out of grid! Mass: %" GSYM " age: %" FSYM "\npos: %" FSYM ", %" FSYM " %" FSYM ", Vel: %" FSYM " %" FSYM " %f\nLeftEdge: %" FSYM " %" FSYM " %" FSYM "\nRightEdge: %" FSYM " %" FSYM " %" FSYM "\n",
                                level, dx, pIndex, ParticleType[pIndex], ParticleMass[pIndex] * MassUnits,
                                age, xp, yp, zp,
                                ParticleVelocity[0][pIndex],
                                ParticleVelocity[1][pIndex],
                                ParticleVelocity[2][pIndex],
                                GridLeftEdge[0], GridLeftEdge[1], GridLeftEdge[2],
                                GridRightEdge[0], GridRightEdge[1], GridRightEdge[2]);
                        fprintf(stderr, "Particle %d: Ct = %f Zzsun = %f TDP = %f\n",
                                pIndex, ParticleAttribute[0][pIndex],
                                ParticleAttribute[2][pIndex],
                                ParticleAttribute[1][pIndex]);
                        fprintf(stderr, "Accelerations: %"FSYM", %"FSYM", %"FSYM"\n",
                                ParticleAcceleration[0][pIndex], ParticleAcceleration[1][pIndex],
                                ParticleAcceleration[2][pIndex]);
                    }
                    continue;
                        // EnzoFatalException("Star Maker Mechanical: particle out of grid!\n");
                        //raise(SIGABRT); // fails more quickly and puts out a core dump for analysis
                    }
            int shifted = 0;

            if (xp < CellLeftEdge[0][0] + borderDx)
            {
                printf("Shifting particle X %e -> %e\n", xp, xp + borderDx+0.5*dx);
                xp = CellLeftEdge[0][0] + borderDx + 0.5 * dx;
                shifted++;
            }
            if (xp > CellLeftEdge[0][0] + gridDx - borderDx)
            {
                printf("Shifting particle X %e -> %e\n", xp, xp - borderDx - 0.5*dx);
                xp = CellLeftEdge[0][0] + gridDx - borderDx - 0.5*dx;
                shifted = 1;
            }
            if (yp < CellLeftEdge[1][0] + borderDx)
            {
                printf("Shifting particle Y %e -> %e\n", yp, yp + borderDx+0.5*dx);
                yp = CellLeftEdge[1][0] + borderDx + 0.5 * dx;
                shifted = 1;
            }
            if (yp > CellLeftEdge[1][0] + gridDy - borderDx)
            {
                printf("Shifting particleY %e -> %e\n", yp, yp - borderDx - 0.5*dx);
                yp = CellLeftEdge[1][0] + gridDy - borderDx - 0.5 * dx;
                shifted = 1;
            }
            if (zp < CellLeftEdge[2][0] + borderDx)
            {
                printf("Shifting particle Z %e -> %e\n", zp, zp + borderDx+0.5*dx);
                zp = CellLeftEdge[2][0] + borderDx + 0.5 * dx;
                shifted = 1;
            }
            if (zp > CellLeftEdge[2][0] + gridDz - borderDx)
            {
                printf("Shifting particle Z %e -> %e\n", zp, zp - borderDx-0.5*dx);
                zp = CellLeftEdge[2][0] + gridDz - borderDx - 0.5 * dx;
                shifted = 1;
            }

            int ip = (xp - CellLeftEdge[0][0] - 0.5 * dx) / dx;
            int jp = (yp - CellLeftEdge[1][0] - 0.5 * dx) / dx;
            int kp = (zp - CellLeftEdge[2][0] - 0.5 * dx) / dx;

            index = ip + jp * GridDimension[0] + kp * GridDimension[0] * GridDimension[1];


            float shieldedFraction = 0, dynamicalTime = 0, freeFallTime = 0;
            bool gridShouldFormStars = true, notEnoughMetals = false;
            float zFraction = totalMetal[index];
            float pmassMsun = ParticleMass[pIndex] * MassUnits;
            

            /* Start actual feedback: Supernova calculations */

            int nSNII = 0;
            int nSNIA = 0;
            float SNMassEjected = 0, SNMetalEjected = 0;

            /* determine how many supernova events */

            if (SingleSN)
            {
                /* 
                    Determine SN events from rates (currently taken from
                    Hopkins, 2018)
                 */
                determineSN(age, &nSNII, &nSNIA, pmassMsun,
                            TimeUnits, dtFixed);
                numSN += nSNII + nSNIA;
                if ((nSNII > 0 || nSNIA > 0) && debug)
                    fprintf(stdout, "SUPERNOVAE!!!! %d %d level = %d age = %f\n", nSNII, nSNIA, level, age);
                if (nSNII > 0 || nSNIA > 0)
                {
                    /* set feedback qtys based on number and types of events */
                    /* 1e51 erg per sn */
                    float energySN = (nSNII + nSNIA) * 1e51;

                    /*10.5 Msun ejecta for type II and IA*/
                    SNMassEjected = nSNII*10.5 + nSNIA * 1.5;
                    float starMetal = (ParticleAttribute[2][pIndex] / Zsolar); //determines metal content of SNeII

                    MechStars_DepositFeedback(energySN, SNMassEjected, SNMetalEjected, totalMetal, Temperature,
                                              &ParticleVelocity[0][pIndex], &ParticleVelocity[1][pIndex], &ParticleVelocity[2][pIndex],
                                              &xp, &yp, &zp,
                                              ip, jp, kp, size, mu_field, 0, nSNII, nSNIA, starMetal, 0);

                    // can only track number of events in dynamical time if not using it to determine lifetime

                    if (!StarParticleRadiativeFeedback)
                        ParticleAttribute[1][pIndex] += nSNII + nSNIA;
                    // else we can log to file for tracking...
                    else 
                        fprintf(stdout, "***%d -- %f -- SN %d -- SNIa %d\n", ParticleNumber[pIndex], Time * TimeUnits/Myr_s, nSNII, nSNIA);
                    ParticleMass[pIndex] = ParticleMass[pIndex] - SNMassEjected/MassUnits;
                }
            }

            /* Do the same for winds. Cooling Radius is very small,
            So almost no energy is coupled, but some mass may be. */
            float windEnergy = 0, windMass = 0, windMetals = 0;

            /*
                Ignore very old stars, veryvery young stars, and ones whose mass is depleted
             */
            if (StellarWinds && age > 0.001 && ParticleMass[pIndex] * MassUnits > 1)
            {
                // printf("Checking Winds\n");
                float zZsun = min(ParticleAttribute[2][pIndex] / Zsolar, MechStarsCriticalMetallicity);

                determineWinds(age, &windEnergy, &windMass, &windMetals,
                               ParticleMass[pIndex] * MassUnits, zZsun,
                               TimeUnits, dtFixed);
                if (windMass > 100)
                    fprintf(stdout, "Really High Wind Mass = %e\n", windMass);
                if (windEnergy > 1e5)
                {
                    //printf("Winds: M = %e E=%e\n", windMass, windEnergy);
                    MechStars_DepositFeedback(windEnergy, windMass, windMetals, totalMetal, Temperature,
                                              &ParticleVelocity[0][pIndex], &ParticleVelocity[1][pIndex], &ParticleVelocity[2][pIndex],
                                              &xp, &yp, &zp,
                                              ip, jp, kp, size, mu_field, 1, 0, 0, 0.0, 0);
                    ParticleMass[pIndex] = ParticleMass[pIndex] - windMass/MassUnits;
                
                }
            }
            
            /* 
            if these stars are used in conjunction with FLDImplicit or FLDSplit.  This functionality has not been verified
            */
            if (StarMakerEmissivityField)
            {
                MechStars_depositEmissivityField(index, CellWidth[0][0], BaryonField[EmisNum],
                                                 age, ParticleMass[pIndex] * MassUnits, TimeUnits, dtFixed);
            }
        }
    } // end iteration over particles
    if (numSN > 0)
    {
        fprintf(stdout, "Ptcl Number = %d Events = %d FeedbackTime = %e Size = %d\n",
                c, numSN, MPI_Wtime() - startFB, GridDimension[0] * GridDimension[1] * GridDimension[2]);
    }
    return SUCCESS;
}
