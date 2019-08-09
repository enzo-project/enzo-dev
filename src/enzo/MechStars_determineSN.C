/*
    Probabilistically determines supernova based on analytic
    starburst99 simulations.  fits taken from Hopkins 2017

    07/2019: Azton Wells
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "StarParticleData.h"
#include "phys_constants.h"

int determineSN(float age, int* nSNII, int* nSNIA, 
                float massMsun, float TimeUnits, float dt){

    if (NEvents > 0){
        *nSNII = 1;
        NEvents -= 1;
        return SUCCESS;
    }
    /* else, calculate SN rate, probability and determine number of events */
    int seed = clock();
    float RII=0, RIA=0, PII=0, PIA=0, random = 0;
    if (SingleSN == 1 && NEvents < 0)
    {   
        // printf("Calculating rates\n");
        /* age-dependent rates */
        if (age < 3.401)
        {
            RII = 0.0;
            RIA = 0.0;
        }
        if (3.401 <= age < 10.37)
        {
                RII = 5.408e-4;
                RIA = 0.0;
        }
        if (10.37 <= age < 37.53)
        {
                RII = 2.516e-4;
                RIA = 0.0;
        }
        if (37.52 <= age)
        {
                RII = 0.0;
                RIA = 5.2e-8+1.6e-5*exp(-1.0*pow((age-50.0)/10.0, 2)/2.0);
        }
        // printf("Rates: %f %f %f\n", age, RII, RIA);
        /* rates -> probabilities */
        if (RII > 0){
            srand(seed);
            PII = RII * massMsun / 3.1557e13 *TimeUnits*dt;
            // printf("PII =%f\n %f %e %f\n", PII, RII, massMsun, age);
            random = float(rand())/float(RAND_MAX);
            if (PII > 1.0 && UnrestrictedSN == TRUE){
                int round = (int)PII;
                *nSNII = round;
                PII -= round;
            }
            int psn = *nSNII;
            if (random < PII){
                *nSNII = psn+1;
            }
        }
        // printf("RANDOM = %f\n", random);            
        // printf("N SNII=%d\n",*nSNII);
        
        if (RIA > 0){
            srand(seed);
            PIA = RIA*massMsun/3.1557e13*TimeUnits*dt;
            float random = float(rand())/float(RAND_MAX);
            
            if (PIA > 1.0 && UnrestrictedSN == TRUE)
            {
                int round = int(PIA);
                *nSNIA = round;
                PIA -= round;
            }
            int psn = *nSNIA;
            
            if (random < PIA)
                *nSNIA = psn+1;
        }
    }
        return SUCCESS;
    
}