/*
    Couples the mechanical stars to the radiation machinery in ENZO by filling
    in the emissivity0 field.  
    Code must be compiled with "make emissivity-yes" and "make photon-yes".
    Usage at runtime determined by the StarMakerUseEmissivity flag.
    Unlike the CIC depositions in the rest of this module, the emissivity is set 
    solely for the cell hosting the star particle (or its kicked location).

    The radiation deposited here is time varying depending on the age of the particle (and mass).  
    A better implementation will make the individual bands of radiation  time dependent
    (more UV early, more IR late).
 */
#include <stdio.h>
#include <math.h>
#include <mpi.h>
#include <string.h>
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
//#include "gFLDProblem.h"


int MechStars_depositEmissivityField(int index, float cellwidth,
                float* emissivity0, float age, float mass)
{
    // Modeling after fFLDSplit_RadiationSource.F90 to fill in Emissivity0 field.
        //    etaconst = h_nu0*NGammaDot*specconst/dV
        // specconst ~ scaling factor for spectrum -> varies depeding on user choice of tracked radiation
        // dV ~ proper volume of cell
        // NGammaDot ~ photons per second. -> this will vary depending on age of particle!
        const float LsunToErg = 3.85e33; // erg/s
        const float evPerErg = 6.2415e11;
        const float h_nu0 = 13.6/evPerErg; // erg
        /*
            Calculate rates of photons given the age-based luminosity in Hopkins 2017.  Units are
            L_sun/M_sun.  While they are given, we dont bother with IR/optical bands here.
            We only couple the ionizing radiation, Psi_ion. The others are calculated if they happen
            to be used in the future.
         */
        float Psi_fuv = 0.0, Psi_ion = 0.0;
        if (age < 3.4){
            Psi_fuv = 271. * (1.0 + (age/3.4)*(age/3.4));
        }
        if (age < 3.5){
            Psi_ion = 500;
        }
        if (age > 3.5 && age < 25){
            Psi_ion = 60.*pow(age/3.5, -3.6)+470*pow(age/3.5, 0.045-1.82*log(age));
        }     
        if (age > 3.4){
            Psi_fuv = 572.*pow(age/3.4, -1.5);
        }
        // convert to better units
        Psi_ion *= mass; // L_sun
        Psi_ion *= LsunToErg; // erg/s
        /* 
            assuming all those photons are in the HI ionization range, the number
            of photons is 
         */
        float NGammaDot = Psi_ion / h_nu0;
        

        /*
            Select spectrum scaling based on parameters 
            (probably just HI radiation for now)
            This routine only works with HI radiation for now, as the 
            rest of the rates would require another Starburst99 sim to get
         */
        if (MechStarsRadiationSpectrum != -1){
            ENZO_FAIL("MechStars only implemented for RadHydroESpectrum = -1\n");
        }
        const float specconst = 1.0;    
        
        /*
            Apply selected to Emissivity0 in the form of etaconst.  
         */
        emissivity0[index] += pow(cellwidth, 3.0)*specconst*NGammaDot*h_nu0;

        return SUCCESS;
}

