/*
    Calculates the luminosity of a mechanical star based on Hopkins 2018
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


int MechStars_calcPhotonRates(Star* star, const float Time)
{
    // Modeling after fFLDSplit_RadiationSource.F90 to fill in Emissivity0 field.
        //    etaconst = h_nu0*NGammaDot*specconst/dV
        // specconst ~ scaling factor for spectrum -> varies depeding on user choice of tracked radiation
        // dV ~ proper volume of cell
        // NGammaDot ~ photons per second. -> this will vary depending on age of particle!
        const float LsunToErg = 3.85e33; // erg/s
        const float evPerErg = 6.2415e11;
        const float h_nu0 = 13.6/evPerErg; // erg
        float age = Time - star->ReturnBirthTime();
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
        //Psi_ion *= star->ReturnMass(); // L_sun
        return Psi_ion;

        return SUCCESS;
}

