/*
    Determines wind feedback parameters according to fits in Hopkins 2017:
    These fits are known to be erroneous, need to re-run and fit using SB99 sims.

    07/2019: Azton Wells
 */

#include <stdio.h>
#include <math.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "CosmologyParameters.h"
#include "StarParticleData.h"
#include "phys_constants.h"

int determineWinds(float age, float* eWinds, float* mWinds, float* zWinds,
                        float massMsun, float zZsun, float TimeUnits, float dtFixed){

    bool oldEnough = (age < 0.01)?(false):(true);
    float windE = 0,  windM = 0, windZ = 0.0;
    float wind_factor = 0.0;
    float e_factor = 0.0;
 // I dont want to deal with new particles
    // printf("Computing Winds for age = %f, Msun = %e\n", age, massMsun);
    if (StellarWinds && oldEnough && massMsun > 11){

        if (0.001 < age && age < 1.0){
            wind_factor =4.763 * min((0.01 + zZsun), 1.0) ;
        }
        if (1 <= age && age < 3.5){
            wind_factor = 4.763*min(0.01+zZsun, 1.0)* 
                pow(age, 1.45+0.08*min(log(zZsun), 1.0));
        }
        if (3.5 <= age && age < 100){
            wind_factor = 29.4*pow(age/3.5, -3.25)+0.0042;
        
        }
        if (age < 100){
            float d = powl(1+age/2.5, 1.4);
            float a50 = powl(double(age)/10.0, 5.0);
            e_factor = 5.94e4 / d + a50 +4.83;
            
        }
        if (100 <= age){
            e_factor = 4.83;
            wind_factor = 0.42*pow(age/1000, -1.1)/(19.81/log(age));
        }
        windM = massMsun * wind_factor; //Msun/Gyr
        windM = windM*dtFixed*TimeUnits/3.1557e16; //Msun
        // if (debug) printf("First winds mass = %e\nFrom wf = %f, dt=%f Z = %e\n", windM, wind_factor, dtFixed, zZsun);
        //printf("eFactor = %f age = %f\n", e_factor, age);
        if (windM > massMsun){
            printf("Winds too large Mw = %e, Mp = %e age=%f, Z = %e\n",
                windM, massMsun, age, zZsun);
            windM = 0.125*massMsun; // limit loss to huge if necessary.
        }
        windZ = max(zZsun, 0.016+0.0041*max(zZsun, 1.65)+0.0118)*windM;
        windE = e_factor * 1e12 * windM;
        //fprintf(stdout, "Age = %e Ewinds = %e Mwinds = %e Zwinds = %e  Zsun = %e\n",
          //          age, windE, windM, windZ, zZsun);
        *mWinds = windM;
        *zWinds = windZ;
        *eWinds = windE;
    }
    
    return SUCCESS;
}