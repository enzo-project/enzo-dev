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
    // const int max_random = (1<<16);
    if (NEvents > 0){
        *nSNII = 1;
        NEvents -= 1;
        return SUCCESS;
    }
    /* else, calculate SN rate, probability and determine number of events */
    mt_init(clock());
    *nSNII = 0;
    *nSNIA = 0;
    float RII=0, RIA=0, PII=0, PIA=0;
    float random;
    if (SingleSN == 1 && NEvents < 0)
    {   
        // printf("Calculating rates\n");
        /* age-dependent rates */
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
	    //    fprintf(stdout, "Rates: For age %f Myr, RII = %f; RIA = %f ", age, RII, RIA);
        /* rates -> probabilities */
        if (RII > 0){
        // printf("Zcpl = %e", zCouple);
            PII = RII * massMsun / Myr_s *TimeUnits*dt;
            random = float(mt_random())/float(UINT_MAX);
            // fprintf(stdout, "PII =%f -- %f %e %f rnd = %e\n", PII, RII, massMsun, age, random);
            if (PII > 1.0 && UnrestrictedSN == TRUE){
                int round = (int)PII;
                *nSNII = round;
                PII -= round;
            }
            if (PII > 1.0 && !UnrestrictedSN){
                ENZO_FAIL("PII too large!");
            }
            int psn = *nSNII;
            if (random < PII){
                *nSNII = psn+1;
            }
        }
        // if (*nSNII > 0)
            // fprintf(stdout, "Positive SN predicted: RII = %e; PII = %e; dt = %e (%e Myr); M* = %e; A* = %f; Rand = %e\n", \
            //                 RII, PII, dt, dt * TimeUnits / Myr_s, massMsun, age, random);
        // printf("RANDOM = %f\n", random);            
        // printf("N SNII=%d\n",*nSNII);
        
        if (RIA > 0){
            PIA = RIA*massMsun / Myr_s * TimeUnits * dt;
            random = float(mt_random())/float(UINT_MAX);
            
            if (PIA > 1.0 && UnrestrictedSN == TRUE)
            {
                int round = int(PIA);
                *nSNIA = round;
                PIA -= round;
            }
            int psn = *nSNIA;
            
            if (random < PIA)
                *nSNIA = psn+1;
            // if (*nSNIA > 0)
            //     fprintf(stdout, "PIA = %f\n", PIA);
        }
    }
        return SUCCESS;
    
}
