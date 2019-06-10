/***********************************************************************                                                
/                                                                                                                        
/  FIND DEAD STAR PARTICLES AND ADD THEM TO MAGNETIC SUPERNOVA LIST
/                                                                                                                       
/  written by: Iryna Butsky                                                                                                
/  date:       June 2018                                                                                          
/                                                                                                                       
/ PURPOSE: To apply magnetic supernova feedback, we need to first add
/          all stars that have recently gone supernova to the global list. 
/          
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "phys_constants.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "TopGridData.h"
#include "Grid.h"
#include "hydro_rk/SuperNova.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, double *MassUnits, FLOAT Time);

void mt_init(unsigned_int seed);
unsigned_long_int mt_random();

int grid::AddMagneticSupernovaeToList()
{

  if (ProcessorNumber != MyProcessorNumber) {
    return SUCCESS;
  }

  float random_u, random_v, random_phi, random_theta, phi_x, phi_y, phi_z;
  float sn_birthtime, sn_duration, sn_radius, sn_energy, star_birthtime, star_lifetime;

  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits, VelocityUnits, EnergyUnits;
  double  MassUnits;
  if (GetUnits(&DensityUnits, &LengthUnits,&TemperatureUnits, &TimeUnits,
               &VelocityUnits, &MassUnits, Time) == FAIL){
    fprintf(stderr, "Error in GetUnits.\n");
    return FAIL;
  }


  // if UseMagneticSupernovaFeedback > 1, then we set the magnetic feedback radius and duration 
  // below based on the resolution of the grid at the highest refinement level
  if (UseMagneticSupernovaFeedback > 1){
    sn_duration = 5.0 * dtFixed;
    MagneticSupernovaDuration = sn_duration * TimeUnits / yr_s;
    sn_radius = 1.5 * this->CellWidth[0][0];
    MagneticSupernovaRadius = sn_radius * LengthUnits / pc_cm;
  }
  else {
    // Converting time from years to seconds, then internal units
    sn_duration = MagneticSupernovaDuration * yr_s / TimeUnits;
    if(sn_duration < 5.0 * dtFixed) {
      printf("WARNING: Magnetic supernova feedback duration is less than the recommended minimum of 5 x dtFixed\n");
      printf("Current dtFixed = %e years \n", dtFixed * TimeUnits / yr_s);
    }
    // Converting radius from parsecs to cm, then internal units
    sn_radius = MagneticSupernovaRadius * pc_cm / LengthUnits;
    if(sn_radius < 1.5 * this->CellWidth[0][0]){
      printf("WARNING: Magnetic supernova feedback radius is less than the recommended minimum of 1.5 x CellWidth\n"); 
      printf("Current CellWidth = %e pc \n", this->CellWidth[0][0] * LengthUnits / pc_cm);
    }
   
  }
  // Converting energy from ergs to internal units 
  sn_energy = MagneticSupernovaEnergy / (MassUnits*VelocityUnits*VelocityUnits);

  for (int star_id = 0;  star_id < this->NumberOfParticles; star_id++){
    if (this->ParticleType[star_id] == PARTICLE_TYPE_STAR){

      star_birthtime = this->ParticleAttribute[0][star_id];
      star_lifetime = this->ParticleAttribute[1][star_id];
      sn_birthtime = star_birthtime + star_lifetime;

      // Creates a supernova with magnetic feedback set by user-defined parameters and                              
      // adds it to supernova list                                                                                  
      if((this->Time > sn_birthtime) && (this->Time < sn_birthtime + sn_duration)){
	SuperNova P = SuperNova();
        mt_init((unsigned_int) this->ParticleNumber[star_id]);

	random_u = (float)(mt_random()%32768)/32768.0; // random variable from 0 to 1                               
        random_v = (float)(mt_random()%32768)/32768.0;
        random_phi = 2*pi*random_u; // 0 to 2pi                                                                   
        random_theta = acos(2*random_v-1); // 0 to pi                          

	// Setting up randomly oriented magnetic feedback of supernova                                              
        phi_x = sin(random_theta)*cos(random_phi);
        phi_y = sin(random_theta)*sin(random_phi);
        phi_z = cos(random_theta);

        P.setValues(phi_x, phi_y, phi_z,
                    ParticlePosition[0][star_id], ParticlePosition[1][star_id], ParticlePosition[2][star_id],
                    sn_radius, sn_birthtime, sn_duration, sn_energy);
        this->MagneticSupernovaList.push_back(P);
      }
    }
  }
  return SUCCESS;
}
