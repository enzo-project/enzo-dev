

#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"


int transformComovingWithStar(float* Density, float* Metals, 
        float* Vel1, float* Vel2, float* Vel3, 
        float up, float vp, float wp,
        int sizeX, int sizeY, int sizeZ, int direction){
    /* transform velocities to momenta or back and make them comoving with the star particle */
    if (direction > 0){

        /* To comoving with star */
        for (int k = 0; k < sizeZ; k++){
            for (int j = 0; j<sizeY; j++){
                for (int i = 0; i<sizeX; i++){
                    int ind = i + j*sizeX+k*sizeX*sizeY;
                    float preV = Vel1[ind];
                    Vel1[ind] = (preV-up)*Density[ind];
                    preV = Vel2[ind];
                    Vel2[ind] = (preV-vp)*Density[ind];
                    preV = Vel3[ind];
                    Vel3[ind] = (preV-wp)*Density[ind];
                }
            }
        }
    }
    if (direction < 0){

        /* back to "lab" frame */
        for (int k = 0; k < sizeZ; k++){
            for (int j = 0; j<sizeY; j++){
                for (int i = 0; i<sizeX; i++){
                    int ind = i + j*sizeX+k*sizeX*sizeY;
                    Vel1[ind] = Vel1[ind]/Density[ind]+up;
                    Vel2[ind] = Vel2[ind]/Density[ind]+vp;
                    Vel3[ind] = Vel3[ind]/Density[ind]+wp;
                }
            }
        }
    }
    return SUCCESS;

}