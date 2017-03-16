#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "phys_constants.h"


/* ------------------------------------------------
 * Copy of functions in halogen to generate double power law dark
 * matter profile background
 * ------------------------------------------------
 */

double trapez(double (*function)(double), double a, double b, int n);
double DoublePower_integrandMenc(double r);
double DoublePower_integrandPot(double r);
double DoublePower_tau(double r);
double DoublePower_rho(double r);

int search_lower_bound(float *arr, float value, int low, int high,
                       int total);
float LinearInterpolationCoefficient(const int &i, const float &x1, const float *x1a);



double integrate(double (*function)(double), double a, double b){

  const int maxsteps = 28;
  const int minsteps = 5 ;
  const float tol = 1.0E-10;

  double T[3], sa, sb;

  T[0] = trapez(function,a,b,0);
  T[1] = trapez(function,a,b,1);
  T[2] = trapez(function,a,b,2);

  sa = (4.0 * T[1] - T[0])/3.0;
  sb = (4.0 * T[2] - T[1])/3.0;

  for(int i = 3; i < (maxsteps + 1); i++){
    T[0] = T[1];
    T[1] = T[2];
    T[2] = trapez(function, a, b, i);

    sa = sb;
    sb = (4.0 * T[2] - T[1])/3.0;

    if (i > minsteps){
      if ((fabs(sb-sa) < tol*fabs(sb)) || (sa == 0.0 && sb == 0.0)){
        return sb;
      }
    }

  }

  fprintf(stderr, "WARNING TOO MANY STEPS IN INTEGRATION\n");
  return sb;
}

double trapez(double (*function)(double), double a, double b, int n){

  int i,j;
  int *N;
  double deltax, sum, sumN;

  N = new int[n+1];

  if (n==0){
    delete [] N;

    return (0.5 * (b-a) *((*function)(a)+(*function)(b)));
  } else{
    for (i = 0; i < (n+1); i ++){
      N[i] = POW(2,i);
    }
    sum = 0.5 * ((*function)(a) + (*function)(b));

    for(i = 1; i < (n+1); i++){
      deltax = (b-a)/N[i];
      sumN = 0;

      for( j = 1; j < (N[i-1]+1); j++){
        sumN += (*function)(a + (2*j-1)*deltax);
      }

      sum = sum + sumN;
    }


    double value = (sum*(b-a)/N[n]);
    delete [] N;

    return value;
  }

}

int InitializeDoublePowerDarkMatter(void){

  if (DiskGravityDoublePowerMass != NULL)
    return SUCCESS;

  const float f_router = 1.0E20;
  // tabulate dark matter potential as a function of log r and
  // dark matter mass as a function of log r to the virial radius

  float *DiskGravityDoublePowerPoutr;
  DiskGravityDoublePowerMass = new float[DOUBLE_POWER_DG_POINTS];
  DiskGravityDoublePowerPot  = new float[DOUBLE_POWER_DG_POINTS];
  DiskGravityDoublePowerR    = new float[DOUBLE_POWER_DG_POINTS];
  DiskGravityDoublePowerPoutr = new float[DOUBLE_POWER_DG_POINTS];

  float rmin, rmax, dr;

  rmin = log10(0.001 * pc);
  rmax = log10(DiskGravityDarkMatterCutoffR * Mpc);

  dr   = log10((rmax - rmin) / ((float) DOUBLE_POWER_DG_POINTS - 1));

  DiskGravityDoublePowerR[0]    = POW(10.0, rmin);
  DiskGravityDoublePowerMass[0] = integrate(DoublePower_integrandMenc,
                                            DiskGravityDoublePowerR[0] - POW(10.0, dr),
                                            DiskGravityDoublePowerR[0]);

  for (int i = 1; i < DOUBLE_POWER_DG_POINTS; i++){
    DiskGravityDoublePowerR[i]     = POW(10.0, rmin + dr*i);
    DiskGravityDoublePowerMass[i] += integrate(DoublePower_integrandMenc,
                                               DiskGravityDoublePowerR[i-1],
                                               DiskGravityDoublePowerR[i]);
  }

  /* Below is a direct copy of halogen */

  float Potoutr = 0.0;

  float router = DiskGravityDarkMatterR * Mpc;

  while (DoublePower_rho(DiskGravityDarkMatterR*Mpc) / DoublePower_rho(router) < f_router){
    router = router * sqrt(10);
  }


  int index = DOUBLE_POWER_DG_POINTS - 1;
  if (DiskGravityDarkMatterBeta > 3){
    Potoutr += 4.0 * pi * GravConst * DoublePower_rho(router)*(router*router)/(2.0 - DiskGravityDarkMatterBeta);
  }
  float Potr = (-1)*GravConst*(DiskGravityDoublePowerMass[index]/DiskGravityDoublePowerR[index] + Potoutr);

  DiskGravityDoublePowerPot[index]   = Potr;
  DiskGravityDoublePowerPoutr[index] = Potoutr;

  for( int i = DOUBLE_POWER_DG_POINTS - 2; i >=0; i --){
    Potoutr = DiskGravityDoublePowerPoutr[i+1];
    Potoutr += integrate(DoublePower_integrandPot, DiskGravityDoublePowerR[i],
                        DiskGravityDoublePowerR[i+1]);

    Potr = (-1) * GravConst * (DiskGravityDoublePowerMass[i]/DiskGravityDoublePowerR[i] + Potoutr);

    DiskGravityDoublePowerPot[i] = Potr;
    DiskGravityDoublePowerPoutr[i] = Potoutr;

  }

  return SUCCESS;
}

double DoublePowerInterpolateMass(double r){

  float M;

  int index = search_lower_bound(DiskGravityDoublePowerR, r, 0, DOUBLE_POWER_DG_POINTS, DOUBLE_POWER_DG_POINTS);
  float coeff = LinearInterpolationCoefficient(index, r, DiskGravityDoublePowerR);

  M = (1.0 - coeff)*DiskGravityDoublePowerMass[index] + coeff*DiskGravityDoublePowerMass[index-1];

  return M;
}

double DoublePowerInterpolatePotential(double r){

  float pot;

  int index = search_lower_bound(DiskGravityDoublePowerR, r, 0, DOUBLE_POWER_DG_POINTS, DOUBLE_POWER_DG_POINTS);
  float coeff = LinearInterpolationCoefficient(index, r, DiskGravityDoublePowerR);

  pot = (1.0 - coeff)*DiskGravityDoublePowerPot[index] + coeff*DiskGravityDoublePowerPot[index-1];

  return pot;
}


double DoublePower_integrandMenc(double r){
  /* ----------------------------------------
   * DoublePower_integrandMenc
   * ----------------------------------------
   * Integrand to compute volume integral over
   * density profile
   * ----------------------------------------
   */
  return (4.0 * pi * r * r * DoublePower_rho(r));
}

double DoublePower_integrandPot(double r){


  return (4.0 * pi * GravConst * r * DoublePower_rho(r));
}

double DoublePower_tau(double r){
  /* -----------------------------------------
   * DoublePower_tau
   * -----------------------------------------
   * Computes the actual functional form of the
   * density profile
   * -----------------------------------------
   */

  double exp1, exp2, exp3;
  double fac1, fac2, fac3;

  exp1 = DiskGravityDarkMatterGamma;
  exp2 = DiskGravityDarkMatterAlpha;
  exp3 = (DiskGravityDarkMatterBeta - DiskGravityDarkMatterGamma)/
              DiskGravityDarkMatterAlpha;

  fac1 = POW(r / (DiskGravityDarkMatterR*Mpc), exp1);
  fac2 = 1.0 + POW(r/(DiskGravityDarkMatterR*Mpc), exp2);
  fac3 = POW(fac2,exp3);

  return fac1 * fac3;
}


double DoublePower_rho(double r){

  double fac1, fac2;

  if (DiskGravityDarkMatterBeta > 3){
    // finite mass

    return DiskGravityDarkMatterDensity / DoublePower_tau(r);
  } else {
    // Need exponential cutoff

    if (r <= DiskGravityDarkMatterCutoffR*Mpc){
      return DiskGravityDarkMatterDensity / DoublePower_tau(r);

    } else{
      ENZO_FAIL("DoublePower_DiskGravity: Cannot currently compute outside cutoff radius of halo\n");
    }
  }
}
