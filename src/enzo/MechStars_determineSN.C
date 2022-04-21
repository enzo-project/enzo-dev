/*
    Probabilistically determines supernova based on analytic
    starburst99 simulations.  fits taken from Hopkins 2017

    07/2019: Azton Wells
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include <random>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "StarParticleData.h"
#include "phys_constants.h"

void mt_init(unsigned_int seed);
unsigned_long_int mt_random();
int determineSN(float age, int* nSNII, int* nSNIA, 
                float massMsun, float TimeUnits, float dt){
    if (NEvents > 0){
        *nSNII = 1;
        NEvents -= 1;
        return SUCCESS;
    }

    unsigned_long_int max_random = 2e9; // 1'000'000;
    mt_init(clock());
    
    std::random_device rd;
    std::mt19937 e2(rd());
    *nSNII = 0;
    *nSNIA = 0;
    float RII=0, RIA=0, PII=0, PIA=0;
    float random;
    if (SingleSN == 1 && NEvents < 0)
    {   
    
        if (age < 3.401)
        {
            RII = 0.0;
            RIA = 0.0;
        }
        else if (3.401 <= age && age< 10.37)
        {
                RII = 5.408e-4;
                RIA = 0.0;
        }
        else if (10.37 <= age && age < 37.53)
        {
                RII = 2.516e-4;
                RIA = 0.0;
        }
        else if (37.53 <= age)
        {
                RII = 0.0;
                RIA = 5.2e-8+1.6e-5*exp(-1.0*pow((age-50.0)/10.0, 2)/2.0);
        }
        /* rates -> probabilities */
        if (RII > 0){
            PII = RII * massMsun * (TimeUnits*dt) / Myr_s ;
            if (PII > 1 && UnrestrictedSN)
                while (PII > 1){
                    *nSNII += 1;
                    PII -= 1;

            }
            std::binomial_distribution<> dist(1, PII);
            random = dist(e2);
            // fprintf(stdout, "PII =%f -- %f %e %f dt = %f rnd = %e; URS = %lld\n", 
            //                     PII, RII, massMsun, age, dt * TimeUnits/Myr_s, random, UnrestrictedSN);

            if (PII > 1.0 && !UnrestrictedSN){
                ENZO_FAIL("PII too large!");
            }
            // int psn = *nSNII;
            if (random){
                *nSNII += 1;
            }
        }
        if (*nSNII > 0)
            fprintf(stdout, "Positive SN predicted: RII = %e; PII = %e; nSN = %d; dt = %e (%e Myr); M* = %e; A* = %f; Rand = %e\n", \
                            RII, PII, *nSNII, dt, dt * TimeUnits / Myr_s, massMsun, age, random);
        
        if (RIA > 0){
            PIA = RIA*massMsun / Myr_s * TimeUnits * dt;
            if (PIA > 1.0 && UnrestrictedSN)
                while (PIA > 1.0){
                    *nSNIA += 1;
                    PIA -= 1;
                }            
            std::binomial_distribution<> dist(1, PIA);
            random = dist(e2); 
            if (PIA < 1.0 / float(max_random)){
                fprintf(stdout, "WARNING: SN probability < minimum random!  Edit modulo factors of random!\n");
            }            
            if (random)
                *nSNIA += 1;
        }
    }
        return SUCCESS;
    
}
