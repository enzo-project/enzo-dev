#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "macros_and_parameters.h"  // for multimetals access
#include "typedefs.h"
#include "global_data.h" // not sure
#include "phys_constants.h"
// 

#define USE // not sure if needed

/* function prototypes */
int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, double *MassUnits, FLOAT Time);

float SampleIMF(void);

unsigned_long_int mt_random(void);



int individual_star_maker(// a lot of stuff in here
                         )
{
  
  double msun = 1.989e33;
  double m1   = (*d1)*POW(*x1,3)/msun; // mass of gas

  // check jeans length refinement
  if (*jlrefine > 0){
    jlsquared = ((double)((*gamma) * 3.14159621 * 1.38e-16 / 6.673e-08)/
        ((double)(*d1) * 1.673e-24)) / (POW(*x1,2) / (*mu) / POW((*jlrefine),2);
  }

  /* Loop over grid and apply star formation criteria
     Create a star particle if it matches all of the criteria */
  

//  if (*level == MaximumRefinementLevel) { // do only at max lvl

    // 3D -> 1D index adjacent cells
    xo = 1;
    yo = *nx;
    zo = (*nx) * (*ny);

    for (k = *ibuff; k < *nz - *ibuff; k++){
      for (j = *ibuff; j < *ny - *ibuff; j++){
        index = (k * (*ny) + j) * (*nx) + (*ibuff);
        for (i = *ibuff; i < *nx - *ibuff; i++){
        
        // 1) Are we in highest refinement leve
        // 2) Are we past density threshold
        // 3) Are we below the temperature threshold
        // 4) Is there enough gas to form a star? M_baryon > M_min of IMF?  

          bmass = d[index] * m1;

          if (r[index] == 0.0 && d[index] > densthresh && temp[index] <= min_temp
            && bmass > minimum_star_mass){

            // 3) only check if above true. Divergence
            if (*imethod == 2) {
              div = u[index + xo] - u[index] +
                    v[index + yo] - v[index] +
                    w[index + zo] - w[index];
            } else{
              dev = u[index + xo] - u[index - xo] +
                    v[index + yo] - v[index - yo] +
                    w[index + zo] - w[index - zo];        

            } // divergence calculation

            if (div < 0 ){ // if divergence is negative
            
              
              // -------------- create star particles --------------
              while ( ii < *nmax ||  ) {

                // here, make a call to IMF to sample the stellar mass
                starMass = SampleIMF() ;
                //
                
                 
                mp[ii] = starMass;
              
                if (mp[ii] < 2){
                  tdp[ii] = POW(mp[ii], -3.0) ; //Msun**4.0 / Lsun
                }else if (mp[ii] < 20){
                  tdp[ii] = 0.666667 * POW(mp[ii],-2.5); // Msun**3.5 / Lsun

                }else if (mp[ii] > 20){ // UNITS !!!!!!
                  tdp[ii] = 3.125e-4 ; // M_sun / L_sun
                }

                type[ii] = -(*ctype);
                tcp[ii]  = *t; // time of creation

                // give the star particle a position chosen at random over
                // the grid cell size .... random() function different 
                // than mt_random to keep repeatability of IMF draws

                // note to self... there is a bug in the .C star maker methods
                // from the fortran due to the fact that they star at 0 and 1 in 
                // index increments resepctively

                xp[ii] = (*dx) * random() + *xstart + ((float) i + 0.5)*(*dx);
                yp[ii] = (*dx) * random() + *ystart + ((float) j + 0.5)*(*dx);
                zp[ii] = (*dx) * random() + *zstart + ((float) k + 0.5)*(*dx);
                // upper.. from star_maker9 ... i think below should be i + 0.5 etc... 
                //xp[ii] = *xstart + ((float) i - 0.5)*(*dx);
                //yp[ii] = *ystart + ((float) j - 0.5)*(*dx);
                //zp[ii] = *zstart + ((float) k - 0.5)*(*dx);
               

                // Assign velocities depending on Hydro method
                // 2 = Zeus .. otherwise PPM
                // copied from pop3_maker.F
                if (*imethod == 2){
                  up[ii] = (
                       0.5 * (u[index   ] + u[index+xo])*d[index] +
                       0.5 * (u[index-xo] + u[index   ])*d[index-xo] +
                       0.5 * (u[index+xo] + u[index + xo + xo])*d[index+xo] +
                       0.5 * (u[index+yo] + u[index + xo + yo])*d[index+yo] +
                       0.5 * (u[index-yo] + u[index + xo - yo])*d[index-yo] +
                       0.5 * (u[index+zo] + u[index + xo + ko])*d[index+ko] +
                       0.5 * (u[index-zo] + u[index + xo - ko])*d[index-ko]) /
                      ( d[index] + d[index-xo] + d[index+xo] +
                        d[index-yo] + d[index+yo] +
                        d[index-ko] + d[index+ko] ); // 
                 //below copied from above... check for typos 
                  vp[ii] = (
                       0.5 * (v[index   ] + v[index+xo])*d[index] +
                       0.5 * (v[index-xo] + v[index   ])*d[index-xo] +
                       0.5 * (v[index+xo] + v[index + xo + xo])*d[index+xo] +
                       0.5 * (v[index+yo] + v[index + xo + yo])*d[index+yo] +
                       0.5 * (v[index-yo] + v[index + xo - yo])*d[index-yo] +
                       0.5 * (v[index+zo] + v[index + xo + ko])*d[index+ko] +
                       0.5 * (v[index-zo] + v[index + xo - ko])*d[index-ko]) /
                      ( d[index] + d[index-xo] + d[index+xo] +
                       d[index-yo] + d[index+yo] +
                        d[index-ko] + d[index+ko] ); // 
                  wp[ii] = (
                       0.5 * (w[index   ] + w[index+xo])*d[index] +
                       0.5 * (w[index-xo] + w[index   ])*d[index-xo] +
                       0.5 * (w[index+xo] + w[index + xo + xo])*d[index+xo] +
                       0.5 * (w[index+yo] + w[index + xo + yo])*d[index+yo] +
                       0.5 * (w[index-yo] + w[index + xo - yo])*d[index-yo] +
                       0.5 * (w[index+zo] + w[index + xo + ko])*d[index+ko] +
                       0.5 * (w[index-zo] + w[index + xo - ko])*d[index-ko]) /
                      ( d[index] + d[index-xo] + d[index+xo] +
                        d[index-yo] + d[index+yo] +
                        d[index-ko] + d[index+ko] ); // 
                }
                else{
                  
                  up[ii] = (u[index]*d[index] +
                            u[index-xo]*d[index-xo] +
                            u[index+xo]*d[index+xo] +
                            u[index-yo]*d[index-yo] +
                            u[index+yo]*d[index+yo] +
                            u[index+zo]*d[index+zo] +
                            u[index-zo]*d[index-zo] ) /
                            (d[index] + d[index-xo] + d[index+xo] +
                             d[index-yo] + d[index+yo] +
                             d[index-zo] + d[index+zo]);
                  vp[ii] = (v[index]*d[index] +
                            v[index-xo]*d[index-xo] +
                            v[index+xo]*d[index+xo] +
                            v[index-yo]*d[index-yo] +
                            v[index+yo]*d[index+yo] +
                            v[index+zo]*d[index+zo] +
                            v[index-zo]*d[index-zo] ) /
                            (d[index] + d[index-xo] + d[index+xo] +
                             d[index-yo] + d[index+yo] +
                             d[index-zo] + d[index+zo]);

                  wp[ii] = (w[index]*d[index] +
                            w[index-xo]*d[index-xo] +
                            w[index+xo]*d[index+xo] +
                            w[index-yo]*d[index-yo] +
                            w[index+yo]*d[index+yo] +
                            w[index+zo]*d[index+zo] +
                            w[index-zo]*d[index-zo] ) /
                            (d[index] + d[index-xo] + d[index+xo] +
                             d[index-yo] + d[index+yo] +
                             d[index-zo] + d[index+zo]);
                } // imethod velocity assignment
                
                // this is where code would go to assign
                // chemical tags to all of the particles 
                // depending on whether or not multimetals is ON
                if (imetal == 1){
                  metalf(ii) = metal[index];
                } else{
                  metalf(ii) = 0.0;
                }

                

                
                ii++; // increment star index
              } // end while loop for creating stars
              // ---------------------------------------------------
              
              if (ii >= *nmax) {
                fprintf(stdout, "individual_star_maker: Reached max new star particle count");
                return FAIL;
              }
                         

              

            }
              



           } // resolution and density


        } // enx x loop
      } // end y loop
    } // end z loop
//  } // end max ref level if

    if (ii > 0){
      printf("P(%"ISYM"): individual_star_maker[add]: %"ISYM" new star particles\n", *nproc, ii);
    }


    *np = ii;
    return SUCCESS;
}


float SampleIMF(void)
{
  unsigned_long_int random_int = mt_random();
  const int max_random = (1<<16);
  float x = (float) (random_int%max_random) / (float) (max_random);
  float dm = log10(IndividualStarUpperMassCutoff / IndividualStarLowerMassCutoff)/ (float) (IMF_TABLE_ENTRIES-1);
  float m;


  int width = IMF_TABLE_ENTRIES/2;
  int bin_number = IMF_TABLE_ENTRIES/2;

  while (width > 1) {
    width /= 2;
    if (x > IMFData[bin_number]){
      bin_number += width;
    } else if (x < IMFData[bin_number]){
      bin_number -= width;
    } else{
      break;
    }
  } // found the bin

  m = IndividualStarLowerMassCutoff * POW(10.0, bin_number * dm);
    
  return m;
}

// put feedback here


