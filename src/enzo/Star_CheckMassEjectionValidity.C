#include <stdlib.h>
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"

#include "StellarYieldsRoutines.h"


void Star::CheckMassEjectionValidity(void){

  if (abs(this->type) == PARTICLE_TYPE_INDIVIDUAL_STAR_UNRESOLVED)
    return;

  float tolerance = 1.05; // 5 percent tolerance
  int sn_fail = 0, wind_fail = 0;

  float total_wind_mass_ejection = this->InterpolateYield(1, -1);
//StellarYieldsInterpolateYield(1, this->yield_table_position[0],
//                                                                    this->yield_table_position[1],
//                                                                    this->BirthMass, this->Metallicity, -1);
  float total_sn_mass_ejection = 9999.9;

  if (this->BirthMass > IndividualStarSNIIMassCutoff &&
      this->BirthMass < IndividualStarDirectCollapseThreshold){
    total_sn_mass_ejection = this->InterpolateYield(0, -1);

//StellarYieldsInterpolateYield(0, this->yield_table_position[0],
//                                                              this->yield_table_position[1],
//                                                              this->BirthMass, this->Metallicity, -1);

  }

  if ( this->sn_mass_ejected > tolerance * total_sn_mass_ejection ){
    sn_fail = TRUE;
  }

  if (this->wind_mass_ejected > tolerance * total_wind_mass_ejection){
    wind_fail = TRUE;
  }

  if (sn_fail || wind_fail){
    this->PrintInfo();
    printf("total_wind_injection = %"ESYM"   total_sn_injection = %"ESYM"\n",total_wind_mass_ejection, total_sn_mass_ejection);
    if (sn_fail && wind_fail){
      ENZO_FAIL("Both supernova and wind mass injection exceed tolerance\n");
    } else if (sn_fail){
      ENZO_FAIL("Supernova mass injection exceeds tolerance\n");
    } else if (wind_fail){
      ENZO_FAIL("Wind mass injection exceeds tolerance\n");
    }
  }


  return;
}
