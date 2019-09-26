

#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "StarParticleData.h"

int transformComovingWithStar(float* Density, float* Metals, 
        float* MetalsSNII, float* MetalsSNIA,
        float* Vel1, float* Vel2, float* Vel3, 
        float* TE, float* GE,
        float up, float vp, float wp,
        int sizeX, int sizeY, int sizeZ, int direction){
    /* transform velocities to momenta or back and make them comoving with the star particle 
        Transform metallicity fields to metal density fields*/
    int size = sizeX*sizeY*sizeZ;
    if (direction > 0){

        /* To comoving with star */
        for (int ind = 0; ind < size; ++ind){
            float mult = Density[ind];
                    TE[ind] *= mult;
                    GE[ind] *= mult;
                    float preV = Vel1[ind];
                    Vel1[ind] = (preV-up)*mult;
                    preV = Vel2[ind];
                    Vel2[ind] = (preV-vp)*mult;
                    preV = Vel3[ind];
                    Vel3[ind] = (preV-wp)*mult;
                    if(StarMakerTypeIaSNe)
                        MetalsSNIA[ind] = MetalsSNIA[ind]*mult;
                    if(StarMakerTypeIISNeMetalField)
                        MetalsSNII[ind] = Metals[ind]*mult;

        }
    }
    if (direction < 0){

        /* back to "lab" frame, metal densities back to metallicities */
        for (int ind = 0; ind < size; ++ind){
            float mult = 1./Density[ind];
                    TE[ind] *= mult;
                    GE[ind] *= mult;
                    Vel1[ind] = Vel1[ind]*mult+up;
                    Vel2[ind] = Vel2[ind]*mult+vp;
                    Vel3[ind] = Vel3[ind]*mult+wp;
                    if(StarMakerTypeIaSNe)
                        MetalsSNIA[ind] = MetalsSNIA[ind]*mult;
                    if(StarMakerTypeIISNeMetalField)
                        MetalsSNII[ind] = Metals[ind]*mult;
        }
    }
    return SUCCESS;

}