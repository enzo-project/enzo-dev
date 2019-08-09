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
                        float* Vel1, float* Vel2, float* Vel3, 
                        float* CoolingTime, int* GridDim,
                        float* shieldedFraction, float* freeFallTime, 
                        float* dynamicalTime, int i, int j, int k, 
                        float Time, float* RefinementField, float CellWidth);
    int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, float *MassUnits, float Time);
    int FindField(int field, int farray[], int numfields);


/* Creation Routine */
int grid::MechStars_Creation(grid* ParticleArray, float* Temperature, 
        float *DMField, int level, float* CoolingTime, 
        int MaximumNumberOfNewParticles, int* NumberOfParticlesSoFar)
{
    if (level < StarMakeLevel) return 0;
    float stretchFactor=1.4;
    bool debug = true;
    // fprintf(stdout, "Preparing to check grids\n");
    // limit creation to level specified in parameter file

    //get field numbers
    int DensNum, GENum, Vel1Num, Vel2Num,Vel3Num, TENum;
    if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
        fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
        return FAIL;
    }
    
    int MetallicityField = FALSE, MetalNum;
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

    
    int rank = GridRank;

    
    float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1,
                TimeUnits = 1, VelocityUnits = 1, MassUnits = 1;
    if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
                &TimeUnits, &VelocityUnits, &MassUnits, this->Time) == FAIL) {
        fprintf(stderr, "Error in GetUnits.\n");
    return FAIL;    
    } 

/*      Define MassUnits so that code * mass = physical (Msun)
            and physical(Msun) / MassUnits = code */
    MassUnits = DensityUnits*pow(LengthUnits*CellWidth[0][0], 3)/SolarMass;
    
    float dx = CellWidth[0][0];
    int GZ = int(NumberOfGhostZones);
    int nCreated = *NumberOfParticlesSoFar;
    fprintf(stdout, "Starting creation with %d prior particles\n",nCreated);
    for (int i = GZ; i < GridDimension[0]-GZ; i++){
        for(int j = GZ; j < GridDimension[1]-GZ; j++){
            for (int k = GZ; k< GridDimension[2]-GZ; k++){   
                /*
                Particle creation has several criteria:

                    1. Density > overdensity
                    2. is flow converging by finite differences
                    3. is cooling time < dynamical time
                    4. is gas mass > jean critical mass
                    5. is gas self shielded by Krumholz & Gnedin 2007 criteria
 
                 */
                float shieldedFraction = 0;
                float freeFallTime = 0;
                float dynamicalTime = 0;
                float Time = this->Time;

                int stat = checkCreationCriteria(BaryonField[DensNum],
                    BaryonField[MetalNum], Temperature, DMField, 
                    BaryonField[Vel1Num], BaryonField[Vel2Num], 
                    BaryonField[Vel3Num],
                    CoolingTime, GridDimension, &shieldedFraction,
                    &freeFallTime, &dynamicalTime, i,j,k,Time, 
                    BaryonField[NumberOfBaryonFields], CellWidth[0][0]);
                              
                
                
                if (stat){
                    int index = i+ j*GridDimension[0]+k*GridDimension[0]*GridDimension[1];

                    /* Determine Mass of new particle */

                    float MassShouldForm =(shieldedFraction * BaryonField[DensNum][index] 
                                    * MassUnits / freeFallTime * this->dtFixed*TimeUnits/3.1557e13);

                    if (MassShouldForm < 0){
                        printf("Negative formation mass: %f %f",shieldedFraction, freeFallTime);
                        continue;
                    }
                    float newMass = min(MassShouldForm,
                                    StarMakerMaximumFormationMass);
                    newMass = newMass/MassUnits;
                    if (newMass > BaryonField[DensNum][index]) exit(136);
                    nCreated ++;
                    float totalDensity = (BaryonField[DensNum][index]+DMField[index])*DensityUnits;
                    dynamicalTime = pow(3.0*pi/32.0/GravConst/totalDensity, 0.5);
                    fprintf(stdout, "DynamicalTime = %e\n", dynamicalTime);
                    ParticleArray->ParticleMass[nCreated] = newMass;
                    ParticleArray->ParticleAttribute[1][nCreated] = dynamicalTime/TimeUnits;
                    ParticleArray->ParticleAttribute[0][nCreated] = Time;
                    ParticleArray->ParticleAttribute[2][nCreated] = BaryonField[MetalNum][index]
                                /BaryonField[DensNum][index]/0.02;
                    ParticleArray->ParticleType[nCreated] = PARTICLE_TYPE_STAR;
                    BaryonField[DensNum][index] -= newMass;

    // type IA metal field isnt here right now
                    // if (StarMakerTypeIASNe) 
                    //     ParticleArray->ParticleAttribute[3][nCreated] = BaryonField[MetalIANum];
                    int preI = index-1;
                    int postI = index+1;
                    int preJ = i+ (j-1)*GridDimension[0]+k*GridDimension[0]*GridDimension[1];
                    int postJ = i+ (j+1)*GridDimension[0]+k*GridDimension[0]*GridDimension[1];
                    int preK = i+ j*GridDimension[0]+(k-1)*GridDimension[0]*GridDimension[1];
                    int postK = i+ j*GridDimension[0]+(k+1)*GridDimension[0]*GridDimension[1]; 
                    if (HydroMethod != 0) 
                        fprintf(stderr,"Mechanical star maker not tested for anything except HydroMethod = 0\n");
                    ParticleArray->ParticleVelocity[0][nCreated] = 
                        0.33*(BaryonField[Vel1Num][preI]+BaryonField[Vel1Num][postI]
                            +BaryonField[Vel1Num][index]);
                    ParticleArray->ParticleVelocity[1][nCreated] = 
                        0.33*(BaryonField[Vel2Num][preI]+BaryonField[Vel2Num][postI]
                            +BaryonField[Vel2Num][index]);
                    ParticleArray->ParticleVelocity[2][nCreated] = 
                        0.33*(BaryonField[Vel3Num][preI]+BaryonField[Vel3Num][postI]
                            +BaryonField[Vel3Num][index]);

                    /* give it position at center of host cell */

                    ParticleArray->ParticlePosition[0][nCreated] = CellLeftEdge[0][0]
                                            +(dx*(float(i)-0.5));
                    ParticleArray->ParticlePosition[1][nCreated] = CellLeftEdge[1][0]
                                            +(dx*(float(j)-0.5));
                    ParticleArray->ParticlePosition[2][nCreated] = CellLeftEdge[2][0]
                                            +(dx*(float(k)-0.5));
                    if (nCreated >= MaximumNumberOfNewParticles) return nCreated;
                    fprintf(stdout,"Created star: %d %d ::: %e %f %e %e::: %f %f %f ::: %f %f %f ::: %d %d ::: %d %d %d\n",
                        nCreated, 
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
                }
            }//end for k
        }//end for j
    } // end for i
    if (nCreated > 0){
        fprintf(stdout, "Created %d star particles\n",nCreated);
    }
    *NumberOfParticlesSoFar = nCreated;
    return nCreated;
}