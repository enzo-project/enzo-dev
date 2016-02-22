#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
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
#include "phys_constants.h"
//#include "StarParticleData.h"



/* function prototypes */
int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, double *MassUnits, FLOAT Time);

float SampleIMF(void);

float GaussianRandomVariable(void);

float compute_lifetime(float *mp);

unsigned_long_int mt_random(void);




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
                          float *metalf, int *type, int *pindex)

 // since I would need to read in all of the chemical tracer fields in the above
 // what I can do in StarParticleHandler is check which (if any) field is defined
 // and then if it exists, pass the field
 // if it doesn't exist then just pass a zero value array of size equal to the grid size...
 // not sure how to do this completely though.


{

  const double msolar = 1.989e33;
  const double sndspdC = 1.3095e8;
  double m1   = (*d1)*POW(*x1,3); // mass units

  int i, j, k, index, ii=0, istar=0, index_presf=0;
  int xo, yo, zo, rsign=1;
  float bmass, div, min_temp, star_mass=0.0, sum_mass=0.0;
  float pstar, mass_to_stars, mass_available, tdyn;
  float dtot, isosndsp2, jeansmass, star_fraction, odthreshold;
  float umean, vmean, wmean, px, py, pz, px_excess, py_excess, pz_excess;
  float rnum;
  const double m_h = 1.673e-24;

  const int max_random = (1<<16);

  int form_method = -1; // tracker for debugging purposes


  if ((*dt) == 0.0){
    printf("DT EQUAL TO ZERO\N");
    return FAIL;
  }

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

    min_temp = 1.0E5; // set conditional based on cooling and metals present

    odthreshold = StarMakerOverDensityThreshold * m_h * (*mu) / (*d1);

    ii = 0; index_presf = 0;
    for (k = *ibuff; k < *nz - *ibuff; k++){
      for (j = *ibuff; j < *ny - *ibuff; j++){
//        index = (k * (*ny) + j) * (*nx) + (*ibuff);
        for (i = *ibuff; i < *nx - *ibuff; i++){
          index = i + (j + k * (*ny)) * (*nx);

          bmass = (d[index]*(*dx)*(*dx)*(*dx)) * m1 / msolar; // in solar masses

          // perform the following easy checks for SF before proceeding
          // 1) Are we in the highest resolution region? (yes if this is called)
          // 2) Is density greater than the density threshold?
          // 3) Is temperature < the minimum temperature?
          // 4) Do not allow star formation if minimum star particle mass is
          //    above some fraction of cell mass. This is very unlikely to occur
          //    in intended use case:
          //    (i.e. baryon mass likely always going to be > ~10 solar masses)

//break individual_star_maker.C:106
          sum_mass = 0.0; index_presf = ii;
          if ( d[index]      > odthreshold
              && temp[index] <= min_temp
              && IndividualStarMassFraction*bmass > IndividualStarIMFLowerMassCutoff
              && 0.5*bmass > IndividualStarIMFUpperMassCutoff){

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
            tdyn = sqrt(3.0 * pi / 32.0 / GravConst / dtot) / (*t1); // in code units
            isosndsp2 = sndspdC * temp[index] ; // not sure
            jeansmass = pi / (6.0 * sqrt(d[index]*(*d1)) *
                            POW(pi * isosndsp2 / GravConst ,1.5)) / msolar; // in solar masses

            if (jeansmass <= bmass){

              //
              star_fraction  = min(StarMakerMassEfficiency*(*dt)/tdyn, 1.0);
              mass_to_stars  = star_fraction * bmass;
              mass_available = StarMakerMassEfficiency * bmass;
              mass_to_stars  = min(mass_to_stars, mass_available);
              // if mass_to_stars greater than available mass, convert
              // all of available mass into stars
              // Frankly this is very unlikely to occur...
              // Tests as of 2/22/16 show NO SF here for at least the first 10^5 stars
              if(mass_to_stars >= mass_available){
                mass_to_stars = mass_available;
                while( ii < *nmax && mass_to_stars > IndividualStarIMFUpperMassCutoff){
                  mp[ii] = SampleIMF();

                  metalf[ii] = -1.0;

                  sum_mass       += mp[ii]; // counter for total star mass formed this cell
                  mass_to_stars  -= mp[ii]; // reduce available mass
                  ii++;
                }
              }

              // Tests (as of 2/22/16) show NO SF here for at least the first 10^5 stars
              if (mass_to_stars > IndividualStarIMFUpperMassCutoff){
                while (ii < *nmax && mass_to_stars > IndividualStarIMFUpperMassCutoff){
                  mp[ii] = SampleIMF();
                  metalf[ii] = -2.0;
                  sum_mass      += mp[ii];
                  mass_to_stars -= mp[ii];
                  ii++;
                }
              }

              if(mass_to_stars > IndividualStarIMFLowerMassCutoff){
                while( ii < *nmax && mass_to_stars > IndividualStarIMFLowerMassCutoff){

                  mp[ii]         = SampleIMF();
                  metalf[ii]     = -6.0;
                  sum_mass      += mp[ii];
                  mass_to_stars -= mp[ii];

                  ii++;

                  if (mass_to_stars < 0.0){
                    mass_to_stars = 0.0;
                  }
                }
                // maybe put a catch statement here to break if ii > *nmax??

              }

/*  Try a more stochastic method (see above) without a cutoff in allowed stars
              if (mass_to_stars > IndividualStarIMFLowerMassCutoff){
                while (ii < *nmax && mass_to_stars > IndividualStarIMFLowerMassCutoff){
                    //mp[ii] = SampleIMF();
                    star_mass = SampleIMF();
                    // accept if star mass doesn't exceed available

                    if (star_mass <= mass_to_stars){
                      mp[ii] = star_mass;
                      metalf[ii] = -3.0;
                      sum_mass      += mp[ii];
                      mass_to_stars -= mp[ii];
                      ii++;
                    } else if (abs(star_mass - mass_to_stars) < 0.1*IndividualStarIMFLowerMassCutoff){
                      mp[ii] = mass_to_stars;
                      metalf[ii] = -4.0;
                      sum_mass      += mp[ii];
                      mass_to_stars -= mp[ii];
                      ii++;
                    }
                }
              }
*/




              // now we are in the Goldbaum et. al. 2015 regime (star_maker_ssn.F)
              // random sample from IMF, calculate probability of that star forming
              if (mass_to_stars < IndividualStarIMFLowerMassCutoff && mass_to_stars > tiny_number){
                star_mass = SampleIMF();
                pstar     = mass_to_stars / star_mass;
                rnum =  (float) (random() % max_random) / (float) (max_random);
                if (rnum < pstar){
                  metalf[ii] = -5.0;
                  mp[ii] = star_mass;
                  sum_mass += mp[ii];
                  ii++;
                }
              }

/*              // -------------- create star particles --------------
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
                       0.5 * (u[index+zo] + u[index + xo + zo])*d[index+zo] +
                       0.5 * (u[index-zo] + u[index + xo - zo])*d[index-zo]) /
                      ( d[index] + d[index-xo] + d[index+xo] +
                        d[index-yo] + d[index+yo] +
                        d[index-zo] + d[index+zo] ); //
                //below copied from above... check for typos
                vmean = (
                         0.5 * (v[index   ] + v[index+xo])*d[index] +
                         0.5 * (v[index-xo] + v[index   ])*d[index-xo] +
                         0.5 * (v[index+xo] + v[index + xo + xo])*d[index+xo] +
                         0.5 * (v[index+yo] + v[index + xo + yo])*d[index+yo] +
                         0.5 * (v[index-yo] + v[index + xo - yo])*d[index-yo] +
                         0.5 * (v[index+zo] + v[index + xo + zo])*d[index+zo] +
                         0.5 * (v[index-zo] + v[index + xo - zo])*d[index-zo]) /
                        ( d[index] + d[index-xo] + d[index+xo] +
                         d[index-yo] + d[index+yo] +
                          d[index-zo] + d[index+zo] ); // 
                wmean = (
                         0.5 * (w[index   ] + w[index+xo])*d[index] +
                         0.5 * (w[index-xo] + w[index   ])*d[index-xo] +
                         0.5 * (w[index+xo] + w[index + xo + xo])*d[index+xo] +
                         0.5 * (w[index+yo] + w[index + xo + yo])*d[index+yo] +
                         0.5 * (w[index-yo] + w[index + xo - yo])*d[index-yo] +
                         0.5 * (w[index+zo] + w[index + xo + zo])*d[index+zo] +
                         0.5 * (w[index-zo] + w[index + xo - zo])*d[index-zo]) /
                        ( d[index] + d[index-xo] + d[index+xo] +
                          d[index-yo] + d[index+yo] +
                          d[index-zo] + d[index+zo] ); // 
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

                // need to make a new particle property thats the stellar lifetime

                type[istar] = (*ctype);
                tcp[istar]    = *t;
                tdp[istar]    = tdyn;
                pindex[istar] = index;
                // mass in code units:
                mp[istar] = mp[istar] * msolar / m1;

                // give the star particle a position chosen at random over
                // the grid cell size .... random() function different 
                // than mt_random to keep repeatability of IMF draws

                // note to self... there is a bug in the .C star maker methods
                // from the fortran due to the fact that they star at 0 and 1 in 
                // index increments resepctively

                rnum =  (float) (random() % max_random) / (float) (max_random);
                xp[istar] = (*dx) * rnum + *xstart + ((float) i + 0.5)*(*dx);
                rnum =  (float) (random() % max_random) / (float) (max_random);
                yp[istar] = (*dx) * rnum + *ystart + ((float) j + 0.5)*(*dx);
                rnum =  (float) (random() % max_random) / (float) (max_random);
                zp[istar] = (*dx) * rnum + *zstart + ((float) k + 0.5)*(*dx);
                // upper.. from star_maker9 ... i think below should be i + 0.5 etc... 
                //xp[ii] = *xstart + ((float) i - 0.5)*(*dx);
                //yp[ii] = *ystart + ((float) j - 0.5)*(*dx);
                //zp[ii] = *zstart + ((float) k - 0.5)*(*dx);

                //
                // need to handle units... GRV returns mean of zero, std of 1
                //

                // assume velocity dispersion is isotropic in each velocity component. Multiply disp by
                // sqrt(1/3) to get disp in each component... taking above velocities as the mean
                rnum  =  (float) (random() % max_random) / (float) (max_random);
                rsign = rnum>0.5 ? 1:-1;
                up[istar] = umean + rsign * GaussianRandomVariable() * IndividualStarVelocityDispersion * 0.577350269*1.0E5*(*t1)/(*x1);

                rnum  =  (float) (random() % max_random) / (float) (max_random);
                rsign = rnum>0.5 ? 1:-1;
                vp[istar] = vmean + rsign * GaussianRandomVariable() * IndividualStarVelocityDispersion * 0.577350269*1.0e5*(*t1)/(*x1);

                rnum  =  (float) (random() % max_random) / (float) (max_random);
                rsign = rnum>0.5 ? 1:-1;
                wp[istar] = wmean + rsign * GaussianRandomVariable() * IndividualStarVelocityDispersion * 0.577350269*1.0E5*(*t1)/(*x1);

                // ENSURE MOMENTUM CONSERVATION!!!!!
                // make running total of momentum in each direction
                px += up[istar]*mp[istar];
                py += vp[istar]*mp[istar];
                pz += wp[istar]*mp[istar];

                // this is where code would go to assign
                // chemical tags to all of the particles 
                // depending on whether or not multimetals is ON
//                if (TestProblemData.UseMetallicityField == 1){
//                  metalf[istar] = metal[index];
                  // if statements here for tagging with individual metal fields
                  // will go here
//                } else{
//                  metalf[istar] = 0.0;
//                }



              } // end while loop for assigning particle properties
              // ---------------------------------------------------

              // ensure zero net momentum from mean velocity
              // momentum of gas converted into stars (sum_mass * umean)
              // should be equal to the total momentum of the stars
              // compute excess momentum and modify star velocity evenly (mass weighted)
              // this is not completely physical, as pre-SF and post-SF gas vel is the same
              sum_mass = sum_mass * msolar / m1; // in code units


              px_excess = umean * sum_mass + px;
              py_excess = vmean * sum_mass + py;
              pz_excess = wmean * sum_mass + pz;

              // remove momentum evenly from each star
              if ( abs(px_excess) > 1.0e-12) {
                for (istar = index_presf; istar < ii; istar++){
//                  up[istar] += (-1.0 * px_excess) * (mp[istar] / sum_mass);
                  up[istar] += (-1.0 * px_excess) / (mp[istar] * (float) (ii-index_presf));
                }
              }
              if ( abs(py_excess) > 1.0e-12){
                for (istar = index_presf; istar < ii; istar++){
                  vp[istar] = (-1.0 * py_excess) / (mp[istar] * (float) (ii-index_presf));
                }
              }
              if ( abs(pz_excess) > 1.0E-12){
                for (istar = index_presf; istar < ii; istar++){
                  wp[istar] = (-1.0 * pz_excess) / (mp[istar] * (float) (ii-index_presf));
                }
              }
              // done with modifying momentum of stars

              // now remove mass from grid
              d[index] = (bmass*msolar/m1 - sum_mass) / ((*dx)*(*dx)*(*dx)) ;


//              for(istar = index_presf; istar < ii; istar++){
//                mp[istar] = mp[istar] * msolar / m1;
//              }
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

  for (int counter = 0; counter < ii; counter++){ // need to convert mass to M / dx^3
    mp[counter] = mp[counter] / ((*dx)*(*dx)*(*dx)); // code units / cell volume
  }

/*
    for (int counter = 0; counter < ii; counter++){
      if(mp[counter]*m1/msolar > IndividualStarIMFUpperMassCutoff){
        printf("individual_star_maker: mass too bigproblem lies in function m = %"FSYM"\n",mp[counter]);
      } else if (mp[counter] * m1 / msolar < 0.1){
        printf("individual_star_maker: there is a units misunderstanding\n");
      } else if (mp[counter] * m1 / msolar < IndividualStarIMFLowerMassCutoff){
        printf("individual_star_maker: toot small mass %"FSYM" %"FSYM" %"FSYM"\n",mp[counter], mp[counter]*m1, mp[counter]*m1/msolar);
      }

    }
*/


    *np = ii;
    return SUCCESS;
}


float SampleIMF(void)
{
  unsigned_long_int random_int = mt_random();
  const int max_random = (1<<16);
  float x = (float) (random_int%max_random) / (float) (max_random);
  float dm = log10(IndividualStarIMFUpperMassCutoff / IndividualStarIMFLowerMassCutoff)/ (float) (IMF_TABLE_ENTRIES-1);
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

  m = IndividualStarIMFLowerMassCutoff * POW(10.0, bin_number * dm);

  IndividualStarIMFCalls++;

  if ( m > IndividualStarIMFUpperMassCutoff || m < IndividualStarIMFLowerMassCutoff){
    printf("individual_star_maker: IMF sampling not working m %"FSYM" %"FSYM"\n",m, IndividualStarIMFUpperMassCutoff);
  }

  if (IndividualStarIMFUpperMassCutoff != 100.0 || IndividualStarIMFLowerMassCutoff != 1.0){
    printf("individual_star_maker: Bounds are wrong!!!!\n");
  }

  return m;
}

// maybe need a new overloaded function of the above to
// rescale for a new maximum mass???? Not trivial if
// IMF is a broken power law


// put feedback here


float GaussianRandomVariable(void){
 // returns gaussian random variable y1

    const int max_random = (1<<16);

    float y1, y2, w;
    float x1 = (float) (random() % max_random) / ((float) (max_random));
    float x2 = (float) (random() % max_random) / ((float) (max_random));

    do {
       x1 = 2.0 * (float) (random() % max_random) / ((float) (max_random)) - 1.0;
       x2 = 2.0 * (float) (random() % max_random) / ((float) (max_random)) - 1.0;
       w  = x1 * x1 + x2 * x2;
    } while ( w >= 1.0);


  w = sqrt( ( -2.0 * log( w ) ) / w );
  y1 = x1 * w;
  y2 = x2 * w;

  return y1;

}


/*int individual_star_feedback(int *nx, int *ny, int *nz,
                             float *d, float *u, float *v, float *w,
                             float *dx, FLOAT *current_time,
                             float *d1, float *x1, float *v1, float *t1,
                             FLOAT *xstart, FLOAT *ystart, FLOAT *zstart,
                             float *ibuff,
                             grid *tg, int *np){*/

int grid::individual_star_feedback(int *nx, int *ny, int *nz,
                                   float *dx,
                                   FLOAT *current_time, float *d1, float *x1,
                                   float *v1, float *t1, FLOAT *xstart, FLOAT *ystart,
                                   FLOAT *zstart, int *ibuff, int *np,
                                   float *ParticleMass, FLOAT *ParticlePosition[],
                                   float *ParticleVelocity[], float *ParticleAttribute[]){

// This is going to be fun.
// To Do:
//   - Compute t - t_creation for all particles
//     and compare against t_lifetime
//   - If t - t_creation < t_lifetime
//   - 
//   1) Eject mass and energy
//      -

// ParticleAttribute 0 = creation time
//                   1 = dynamical time

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, CRNum, B1Num, B2Num, B3Num;

  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;

  int CINum, NINum, OINum, MgINum, SiINum, FeINum, YINum, BaINum, LaINum, EuINum;


  this->DebugCheck("StarParticleHandler");
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }
  if (CRModel)
    if ((CRNum = FindField(CRDensity, FieldType, NumberOfBaryonFields)) < 0)
      ENZO_FAIL("Cannot Find Cosmic Rays");


  if(TestProblemData.MultiMetals >=2){

    if(IdentifyChemicalTracerSpeciesFields(CINum, NINum, OINum, MgINum, SiINum,
                                           FeINum, YINum, BaINum, LaINum, EuINum) == FAIL){
      ENZO_FAIL("Error in IdentifyChemicalTracerSpeciesFields.");
    }

  }


  float mp, particle_age, lifetime, distmass, energy, dratio;
  const double msolar = 1.989e33;                 // solar mass in cgs
  const double speed_of_light = 2.99792458e10 ;
  double m1 = (*d1)*(*x1)*(*x1)*(*x1);            // code mass units
  double e1 = m1 * (*x1) * (*x1) / (*t1) * (*t1); // code energy units

  int ip, jp, kp, index; // particle location

  printf("IF Feedback about to loop over particles %"ISYM"\n",(*np));

  // loop over all star particles
  for(int i = 0; i < (*np); i++){

    // where is the particle?
    ip = int ( (ParticlePosition[0][i] - (*xstart)) / (*dx));
    jp = int ( (ParticlePosition[1][i] - (*ystart)) / (*dx));
    kp = int ( (ParticlePosition[2][i] - (*zstart)) / (*dx));

    mp = ParticleMass[i] * (*dx) * (*dx) * (*dx); // mass in code units
    mp = mp * (m1) / msolar ; // Msun

    // warning if outside current grid
    if( ip < 0 || ip > (*nx) || jp < 0 || jp > (*ny) || kp < 0 || kp > (*nz)){
      printf("Warning: star particle is outside of grid\n");
      return FAIL;
    }

    // Modify position if too close to grid boundary
    // If too close, pretends particles are at a new position, which
    // is the to closest they can be to grid without causing
    // feedback require MPI calls. Does not actually move particle
    if (StarFeedbackDistRadius > 0){
      ip = fmax( (*ibuff) + StarFeedbackDistRadius,

               fmin( (*nx) - (*ibuff) - StarFeedbackDistRadius, ip) );
      jp = fmax( (*ibuff) + StarFeedbackDistRadius,
               fmin( (*ny) - (*ibuff) - StarFeedbackDistRadius, jp) );
      kp = fmax( (*ibuff) + StarFeedbackDistRadius,
               fmin( (*nz) - (*ibuff) - StarFeedbackDistRadius, kp) );
    }

    particle_age = abs( (*current_time) - ParticleAttribute[0][i] );
    lifetime     = compute_lifetime( &mp ) / (*t1);

    // if true, do stellar winds
//    if(abs(*current_time / lifetime) < 1.0){
      // do stellar winds and radiation here
//    }

//    printf("AJE before particle age check\n");
    // Now do end-of-life feedback, but only for massive stars
//    if(abs( *current_time / (ParticleAttribute[0][i] + lifetime) >= 1.0){
    if( particle_age > lifetime ){

      if(mp < IndividualStarTypeIIMassCutoff){

        ParticleMass[i] = 0.0; // delete star and do nothing

      } else if (IndividualStarTypeIIMassCutoff <= mp &&
                 mp < IndividualStarPSNMassCutoff){

        // mass and energy to distribute over local cells
        distmass = StarMassEjectionFraction * ParticleMass[i] / ((float) StarFeedbackDistTotalCells);
        energy   = ParticleMass[i] * StarEnergyToThermalFeedback * 
                   (speed_of_light*speed_of_light/((*v1)*(*v1))) / ((float) StarFeedbackDistTotalCells);

//        energy   = IndividualStarTypeIIEnergy / e1 / ( (float) StarFeedbackDistTotalCells);
        // add energy to surroundings
        for(int kc = kp - StarFeedbackDistRadius; kc <= kp + StarFeedbackDistRadius; kc++){
          int stepk = abs(kc - kp);
          for(int jc = jp - StarFeedbackDistRadius; jc <= jp + StarFeedbackDistRadius; jc++){
            int stepj = stepk + abs(jc - jp);
            for(int ic = ip - StarFeedbackDistRadius; ic <= ip + StarFeedbackDistRadius; ic++){
              int cellstep = stepj + abs(ic - ip);
              index = ic + (jc + kc * (*ny)) * (*nx);


              if (cellstep <= StarFeedbackDistCellStep){ // not sure
                dratio    = 1.0 / (BaryonField[DensNum][index] + distmass);
                BaryonField[TENum][index] = (BaryonField[TENum][index]*BaryonField[DensNum][index] + energy) * dratio;

                if(DualEnergyFormalism){
                  BaryonField[GENum][index] = ((BaryonField[GENum][index] * BaryonField[DensNum][index]) + energy) * dratio;
                }


                // now do metal feedback and chemical yields
                // here
                // if(TestProblemData.UseMetallicityField){
                //
                //  if(MULTIMETALS_METHOD(MULTIMETALS_ALPHA){
                //
                //    }
                //  BaryonField[CINum][index] = BaryonField[CINum][index] * 1.01;
                //  BaryonField[NINum][index] = BaryonField[NINum][index] * 1.05;
                //  BaryonField[OINum][index] = BaryonField[OINum][index] * 0.9;
                //}

                // do mass and momentum feedback
                BaryonField[Vel1Num][index] = BaryonField[Vel1Num][index] * BaryonField[DensNum][index] +
                           distmass * ParticleVelocity[0][i];
                BaryonField[Vel2Num][index] = BaryonField[Vel2Num][index] * BaryonField[DensNum][index] +
                           distmass * ParticleVelocity[1][i];
                BaryonField[Vel3Num][index] = BaryonField[Vel3Num][index] * BaryonField[DensNum][index] +
                           distmass * ParticleVelocity[2][i];

                // add mass to cells and convert vels to des again
                BaryonField[DensNum][index] += distmass;
                BaryonField[Vel1Num][index] /= BaryonField[DensNum][index];
                BaryonField[Vel2Num][index] /= BaryonField[DensNum][index];
                BaryonField[Vel3Num][index] /= BaryonField[DensNum][index];

                if (HydroMethod != 2 && DualEnergyFormalism){
                  BaryonField[TENum][index] = 0.5 * (BaryonField[Vel1Num][index]*BaryonField[Vel1Num][index] + 
                                                     BaryonField[Vel2Num][index]*BaryonField[Vel2Num][index] +
                                                     BaryonField[Vel3Num][index]*BaryonField[Vel3Num][index])
                                                   + BaryonField[GENum][index];
                }

              } // end if distribute


            } // i dist
          } // j dist
        }// k dist


        // particle is now dead
      } // end if Type II SN

      // kill particle
      ParticleMass[i] = 0.0;
    } // end end of life feedback


  } // loop over particles

  return SUCCESS;
}

//
// Function to compute star lifetime based on age
// and metallicity. Input mass in Msun.
// Age in cgs
float compute_lifetime(float *mp){

  const double msolar = 1.989e33;
  const double myr    = 3.1536e13;

  float age = 2.0 * myr;

  return age;
}
