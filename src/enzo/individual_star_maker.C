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

float GaussianRandomVariable(void);

unsigned_long_int mt_random(void);

float random(void);


// need to initialize random!!!


int individual_star_maker(int *nx, int *ny, int *nz, int *size,
                          float *d, float *dm, float *temp, float *u,
                          float *v, float *w,  float *dt,
                          float *dx, FLOAT *t, float *z, int *procnum,
                          float *d1, float *x1, float *v1, float *t1,
                          int *nmax, FLOAT *xstart, FLOAT *ystart,
                          FLOAT *zstart, int *ibuff, int *imethod,
                          float *mu, float *metal, int *ctype,
                          int *np, FLOAT *xp, FLOAT *yp, FLOAT *zp, float *up,
                          float *vp, float *wp, float *mp, float *tdp, float *tcp,
                          float *metalf, int *type)

 // since I would need to read in all of the chemical tracer fields in the above
 // what I can do in StarParticleHandler is check which (if any) field is defined
 // and then if it exists, pass the field
 // if it doesn't exist then just pass a zero value array of size equal to the grid size...
 // not sure how to do this completely though.


{

  const double msun = 1.989e33;
  double m1   = (*d1)*POW(*x1,3)/msun; // mass of gas

  int i, j, k, index, ii, istar, index_presf;
  int xo, yo, zo;
  float bmass, div, min_temp, star_mass, sum_mass;
  float pstar, mass_to_stars, mass_available, tdyn;
  float dtot, isosndsp2, jeansmass, star_fraction, odthreshold;
  float umean, vmean, wmean, px, py, pz, px_excess, py_excess, pz_excess;
  const double m_h = 1.673e-24;

  // check jeans length refinement
  //if (*jlrefine > 0){
  //  jlsquared = ((double)((*gamma) * 3.14159621 * 1.38e-16 / 6.673e-08)/
  //      ((double)(*d1) * 1.673e-24)) / (POW(*x1,2) / (*mu) / POW((*jlrefine),2);
  //}

  /* Loop over grid and apply star formation criteria
     Create a star particle if it matches all of the criteria */

    // 3D -> 1D index adjacent cells
    xo = 1;
    yo = *nx;
    zo = (*nx) * (*ny);

    min_temp = 1.0E4; // set conditional based on cooling and metals present

    odthreshold = StarMakerOverDensityThreshold * m_h * (*mu) / (*d1);

    for (k = *ibuff; k < *nz - *ibuff; k++){
      for (j = *ibuff; j < *ny - *ibuff; j++){
        index = (k * (*ny) + j) * (*nx) + (*ibuff);
        for (i = *ibuff; i < *nx - *ibuff; i++){


          bmass = d[index] * m1;

          // perform the following easy checks for SF before proceeding
          // 1) Are we in the highest resolution region? (yes if this is called)
          // 2) Is density greater than the density threshold?
          // 3) Is temperature < the minimum temperature?
          // 4) Do not allow star formation if minimum star particle mass is
          //    above some fraction of cell mass. This is very unlikely to occur
          //    in intended use case:
          //    (i.e. baryon mass likely always going to be > ~10 solar masses)
          if ( d[index]      > odthreshold
              && temp[index] <= min_temp
              && IndividualStarMassFracion*bmass > IndividualStarIMFLowerMassCutoff){

            // 5) only check if above true. Divergence
            //if (*imethod == 2) {
            //  div = u[index + xo] - u[index] +
            //        v[index + yo] - v[index] +
            //        w[index + zo] - w[index];
            //} else{
            //  div = u[index + xo] - u[index - xo] +
            //        v[index + yo] - v[index - yo] +
            //       w[index + zo] - w[index - zo];

            //} // divergence calculation

            // now that the above is done, compute some values
            dtot = ( d[index] + dm[index] ) * (*d1); // total density
            tdyn = sqrt(3.0 * pi / 32.0 / GravConst / dtot) / (*t1);
            isosndsp2 = sndspdC * temp[index] ; // not sure
            jeansmass = pi / (6.0 * sqrt(d[index]*(*d1)) *
                            POW(pi * isosndsp2 / GravConst),1.5) / msolar;


            if (jeansmass <= bmass){
            //if ( 1 == true){ // temporary if statement to replace div check
            //if (div < 0 ){ // if divergence is negative



/*              // probability of star formation
              if (star_fraction * bmass < IndividualStarIMFLowerMassCutoff){
                pstar = starfraction * bmass / IndividualStarIMFLowerMassCutoff;
                if (random() < pstar){
                    star_fraction = IndividualStarIMFLowerMassCutoff / bmass;
                }

              }
*/
              sum_mass = 0.0;
              index_presf = ii;
              //
              star_fraction  = min(StarMakerMassEfficienty*(*dt)/tdyn,1.0);
              mass_to_stars  = star_fraction * bmass;
              mass_available = StarMakerMassEfficiency * bmass;
              mass_to_stars  = min(mass_to_stars, mass_available);

              // if mass_to_stars greater than available mass, convert
              // all of available mass into stars
              // Frankly this is very unlikely to occur...
              if(mass_to_stars >= mass_available){
                while( ii < *nmax && mass_available > IndividualStarIMFUpperMassCutoff){
                  mp[ii] = SampleIMF();

                  sum_mass       += mp[ii]; // counter for total star mass formed this cell
                  mass_available -= mp[ii]; // reduce available mass
                  ii++;
                }
              }

              if (mass_to_stars > IndividualStarIMFUpperMassCutoff){
                while (ii < *nmax && mass_to_stars > IndividualStarIMFUpperMassCutoff){
                  mp[ii] = SampleIMF();
                  sum_mass      += mp[ii];
                  mass_to_stars -= mp[ii];
                  i++;
                }
              }

              if (mass_to_stars > IndividualStarIMFLowerMassCutoff){
                while (ii < *nmax && mass_to_stars > IndividualStarIMFLowerMassCutoff){
                    mp[ii] = SampleIMF();
                    // accept if star mass doesn't exceed available
                    if( mass_to_stars >= mp[ii]){
                      sum_mass += mp[ii];
                      mass_to_stars -= mp[ii];
                      i++;
                    }
                }
              }

              // now we are in the Goldbaum et. al. 2015 regime (star_maker_ssn.F)
              // random sample from IMF, calculate probability of that star forming
              if (mass_to_stars < IndividualStarIMFLowerMassCutoff){
                star_mass = sampleIMF();
                pstar     = mass_to_stars / star_mass;
                if (random() < pstar){
                  mp[ii] = star_mass;
                  sum_mass += mp[ii];
                  i++;
                }
              }


/*
              // -------------- create star particles --------------
              while ( ii < *nmax &&
                      // change below condition
                      sum_mass < star_fraction * bmass) {

                // here, make a call to IMF to sample the stellar mass
                mp[ii] = SampleIMF();

                if (mp[ii] < 2){
                  tdp[ii] = POW(mp[ii], -3.0) ; //Msun**4.0 / Lsun
                }else if (mp[ii] < 20){
                  tdp[ii] = 0.666667 * POW(mp[ii],-2.5); // Msun**3.5 / Lsun

                }else if (mp[ii] > 20){ // UNITS !!!!!!
                  tdp[ii] = 3.125e-4 ; // M_sun / L_sun
                }

                type[ii] = -(*ctype);
                tcp[ii]  = *t; // time of creation

                if (ii >= *nmax) {
                   fprintf(stdout, "individual_star_maker: Reached max new star particle$
                   return FAIL;
                 }

                star_mass += mp[ii];
                ii++;
              } // loop over particles done
*/

              // prepare for assigning star properties by computing the local
              // gas velocity properties
              // 2 = Zeus .. otherwise PPM
              // copied from pop3_maker.F
              if (*imethod == 2){
                umean = (
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
                vmean = (
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
                wmean = (
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
                umean = (u[index]*d[index] +
                              u[index-xo]*d[index-xo] +
                              u[index+xo]*d[index+xo] +
                              u[index-yo]*d[index-yo] +
                              u[index+yo]*d[index+yo] +
                              u[index+zo]*d[index+zo] +
                              u[index-zo]*d[index-zo] ) /
                              (d[index] + d[index-xo] + d[index+xo] +
                               d[index-yo] + d[index+yo] +
                               d[index-zo] + d[index+zo]);
                vmean = (v[index]*d[index] +
                              v[index-xo]*d[index-xo] +
                              v[index+xo]*d[index+xo] +
                              v[index-yo]*d[index-yo] +
                              v[index+yo]*d[index+yo] +
                              v[index+zo]*d[index+zo] +
                              v[index-zo]*d[index-zo] ) /
                              (d[index] + d[index-xo] + d[index+xo] +
                               d[index-yo] + d[index+yo] +
                               d[index-zo] + d[index+zo]);

                wmean = (w[index]*d[index] +
                              w[index-xo]*d[index-xo] +
                              w[index+xo]*d[index+xo] +
                              w[index-yo]*d[index-yo] +
                              w[index+yo]*d[index+yo] +
                              w[index+zo]*d[index+zo] +
                              w[index-zo]*d[index-zo] ) /
                              (d[index] + d[index-xo] + d[index+xo] +
                               d[index-yo] + d[index+yo] +
                               d[index-zo] + d[index+zo]);
              } // imethod velocity computation





              // now assign particle properties
              px = 0.0; py = 0.0; pz =0.0; // initialize momentum counters
              for (istar = index_presf; istar < ii; istar++){

                //if (mp[istar] <2){ tdp[istar] = POW(mp[istar], -3.0); }
                //else if (mp[istar] < 20){ tdp[istar] = 0.666666667 *POW(mp[istar],-2.5);} // units!!
                //else    {tdp[istar] = 3.125e-4;}
                // need to make a new particle property thats the stellar lifetime

                type[istar] = -(*ctype)
                tcp[istar]  = *t;
                tdp[istar]  = tdyn;

                // give the star particle a position chosen at random over
                // the grid cell size .... random() function different 
                // than mt_random to keep repeatability of IMF draws

                // note to self... there is a bug in the .C star maker methods
                // from the fortran due to the fact that they star at 0 and 1 in 
                // index increments resepctively

                xp[istar] = (*dx) * random() + *xstart + ((float) i + 0.5)*(*dx);
                yp[istar] = (*dx) * random() + *ystart + ((float) j + 0.5)*(*dx);
                zp[istar] = (*dx) * random() + *zstart + ((float) k + 0.5)*(*dx);
                // upper.. from star_maker9 ... i think below should be i + 0.5 etc... 
                //xp[ii] = *xstart + ((float) i - 0.5)*(*dx);
                //yp[ii] = *ystart + ((float) j - 0.5)*(*dx);
                //zp[ii] = *zstart + ((float) k - 0.5)*(*dx);

/*
                  // Assign velocities depending on Hydro method
                  // 2 = Zeus .. otherwise PPM
                  // copied from pop3_maker.F
                  if (*imethod == 2){
                    up[istar] = (
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
                    vp[istar] = (
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
                    wp[istar] = (
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
                    up[istar] = (u[index]*d[index] +
                              u[index-xo]*d[index-xo] +
                              u[index+xo]*d[index+xo] +
                              u[index-yo]*d[index-yo] +
                              u[index+yo]*d[index+yo] +
                              u[index+zo]*d[index+zo] +
                              u[index-zo]*d[index-zo] ) /
                              (d[index] + d[index-xo] + d[index+xo] +
                               d[index-yo] + d[index+yo] +
                               d[index-zo] + d[index+zo]);
                    vp[istar] = (v[index]*d[index] +
                              v[index-xo]*d[index-xo] +
                              v[index+xo]*d[index+xo] +
                              v[index-yo]*d[index-yo] +
                              v[index+yo]*d[index+yo] +
                              v[index+zo]*d[index+zo] +
                              v[index-zo]*d[index-zo] ) /
                              (d[index] + d[index-xo] + d[index+xo] +
                               d[index-yo] + d[index+yo] +
                               d[index-zo] + d[index+zo]);

                    wp[istar] = (w[index]*d[index] +
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
*/

                //
                // need to handle units... GRV returns mean of zero, std of 1
                //

                // assume velocity dispersion is isotropic in each velocity component. Multiply disp by
                // sqrt(1/3) to get disp in each component... taking above velocities as the mean
                up[istar] = umean + GaussianRandomVariable() * IndividualStarVelocityDispersion * 0.577350269;
                vp[istar] = vmean + GaussianRandomVariable() * IndividualStarVelocityDispersion * 0.577350269;
                wp[istar] = wmean + GaussianRandomVariable() * IndividualStarVelocityDispersion * 0.577350269;

                // ENSURE MOMENTUM CONSERVATION!!!!!
                // make running total of momentum in each direction
                px += up[istar]*mp[istar];
                py += vp[istar]*mp[istar];
                pz += wp[istar]*mp[istar];

                // this is where code would go to assign
                // chemical tags to all of the particles 
                // depending on whether or not multimetals is ON
                if (TestProblemData.UseMetallicityField == 1){
                  metalf[istar] = metal[index];
                  // if statements here for tagging with individual metal fields
                  // will go here
                } else{
                  metalf[istar] = 0.0;
                }

              } // end while loop for assigning particle properties
              // ---------------------------------------------------

              // ensure zero net momentum from mean velocity
              // momentum of gas converted into stars (sum_mass * umean)
              // should be equal to the total momentum of the stars
              // compute excess momentum and modify star velocity evenly (mass weighted)
              // this is not completely physical, as pre-SF and post-SF gas vel is the same

              px_excess = umean * sum_mass + px;
              py_excess = vmean * sum_mass + py;
              pz_excess = wmean * sum_mass + pz;

              // remove momentum evenly from each star, weighted by mass
              if ( abs(px_excess) > 1.0e-12) {
                for (istar = index_presf; istar < ii; istar++){
                  up[istar] += (-1.0 * px_excess) * (mp[istar] / sum_mass);
                }
              }
              if ( abs(py_excess) > 1.0e-12){
                for (istar = index_presf; istar < ii; istar++){
                  vp[istar] = (-1.0 * py_excess) * (mp[istar] / sum_mass);
                }
              }
              if ( abs(pz_excess) > 1.0E-12){
                for (istar = index_presf; istar < ii; istar++){
                  wp[istar] = (-1.0 * pz_excess) * (mp[istar] / sum_mass);
                }
              }
              // done with modifying momentum of stars

              // now remove mass from grid
              d[index] = (d[index] * (*dx)*(*dx)*(*dx) - sum_mass*msolar) /
                                    ((*dx)*(*dx)*(*dx));


            } // if jeans mass unstable


          } // resolution and density


        } // enx x loop
      } // end y loop
    } // end z loop

    // Done forming stars!!! Output and exit
    if (ii > 0){
      printf("P(%"ISYM"): individual_star_maker[add]: %"ISYM" new star particles\n", *procnum, ii);
    }
    if (ii >= *nmax){
      fprintf(stdout, "individual_star_maker: reached max new particle count!! Available: %"ISYM". Made: %"ISYM"\n", *nmax, ii);

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

// maybe need a new overloaded function of the above to
// rescale for a new maximum mass???? Not trivial if
// IMF is a broken power law


// put feedback here


float GaussianRandomVariable(void){
 // returns gaussian random variable y1

    float y1, y2, w; 
    float x1 = random(), x2 = random();

    do {
       x1 = 2.0 * random() - 1.0;
       x2 = 2.0 * random() - 1.0;
       w  = x1 * x1 + x2 * x2;
    } while ( w >= 1.0);


  w = sqrt( ( -2.0 * log( w ) ) / w );
  y1 = x1 * w;
  y2 = x2 * w;

  return y1;

}

