/***********************************************************************
/
/  GRID CLASS (FLAG CELLS TO BE REFINED BY THE JEAN'S CRITERION)
/
/  written by: Tom Abel 
/  date:       October 2010
/  modified1:  
/
/  PURPOSE: Refine on the Jeans length determined from gradients in the forces
/           rather than purely the baryons own self-gravity
/
/  RETURNS:
/    number of flagged cells, or -1 on failure
/
************************************************************************/
 
#include <stdio.h>
#include <math.h>
#include "phys_constants.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "hydro_rk/EOS.h"
 
/* function prototypes */
 
int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);
 
int FindField(int f, int farray[], int n);
 
int grid::FlagCellsToBeRefinedByTotalJeansLength()
{
  /* declarations */
 
  int i, dim;
 
  /* error check */
 
  if (FlaggingField == NULL) {
    fprintf(stderr, "Flagging Field is undefined.\n");
    return -1;
  }
 
  /* compute size */
 
  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* Compute the temperature field. */
 
  float *temperature = NULL;
  if (ProblemType != 60 && ProblemType != 61 && EOSType == 0) { //AK
    temperature = new float[size];
    if (this->ComputeTemperatureField(temperature) == FAIL) {
      fprintf(stderr, "Error in grid->ComputeTemperature.\n");
      return -1;
    }
    /* This is less efficient, but it avoids too many conditionals */
    if(JeansRefinementColdTemperature > 0.0){
      for (i = 0; i < size; i++) temperature[i] =
        JeansRefinementColdTemperature;
    }
  }
 
  /* Find fields: density, total energy, velocity1-3. */
 
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, GPotNum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }


  // set up temporary array to hold maximal relevant densities as computed from acceleration gradients
  float *MaxDensity = NULL;
  MaxDensity = new float[size];
  FLOAT CellWidthSquared = CellWidth[0][0]*CellWidth[0][0];
  int j, k, index;
  for (i=0;i<size;i++) MaxDensity[i] = tiny_number;
  float rhox, rhoy, rhoz, rhoxy, rhoxz, rhoyz,  rhomax, maxmax;

  if ((GPotNum = FindField(GravPotential, FieldType, NumberOfBaryonFields)) < 0) {
    ENZO_FAIL("Cannot find Gravitational Potential. Set WritePotential = 1 ... hack");
  };
  //  printf("GpotNum : %i\n", GPotNum);
  rhox = rhoy = rhoz = 0.;
  rhomax = tiny_number;
  maxmax = tiny_number;
  for (k = GridStartIndex[2]; k <= GridEndIndex[2]; k++)
    for (j = GridStartIndex[1]; j <= GridEndIndex[1]; j++) {
      for (i = GridStartIndex[0]; i <= GridEndIndex[0]; i++) {
	index = GRIDINDEX_NOGHOST(i,j,k);
	rhox = BaryonField[GPotNum][GRIDINDEX_NOGHOST(min(i+1,GridEndIndex[0]),j,k)] 
	  + BaryonField[GPotNum][GRIDINDEX_NOGHOST(max(i-1,GridStartIndex[0]),j,k)]
	  -2.*BaryonField[GPotNum][index];
	rhoy = BaryonField[GPotNum][GRIDINDEX_NOGHOST(i,min(j+1,GridEndIndex[1]),k)]
	  + BaryonField[GPotNum][GRIDINDEX_NOGHOST(i,max(j-1,GridStartIndex[1]),k)]
	  -2.*BaryonField[GPotNum][index];
	rhoz = BaryonField[GPotNum][GRIDINDEX_NOGHOST(i,j,min(k+1,GridEndIndex[2]))] 
	  + BaryonField[GPotNum][GRIDINDEX_NOGHOST(i,j,max(k-1,GridStartIndex[2]))]
	  -2.*BaryonField[GPotNum][index];
	rhoxy = (BaryonField[GPotNum][GRIDINDEX_NOGHOST(i+1,j+1,k)] +
		 BaryonField[GPotNum][GRIDINDEX_NOGHOST(i-1,j-1,k)] -
		 BaryonField[GPotNum][GRIDINDEX_NOGHOST(i+1,j-1,k)] -
		 BaryonField[GPotNum][GRIDINDEX_NOGHOST(i-1,j+1,k)])/4; 
	rhoyz = (BaryonField[GPotNum][GRIDINDEX_NOGHOST(i,j+1,k+1)] +
		 BaryonField[GPotNum][GRIDINDEX_NOGHOST(i,j-1,k-1)] - 
		 BaryonField[GPotNum][GRIDINDEX_NOGHOST(i,j-1,k+1)] -
		 BaryonField[GPotNum][GRIDINDEX_NOGHOST(i,j+1,k-1)])/4; 
	rhoxz = (BaryonField[GPotNum][GRIDINDEX_NOGHOST(i+1,j,k+1)] +
		 BaryonField[GPotNum][GRIDINDEX_NOGHOST(i-1,j,k-1)] +
		 BaryonField[GPotNum][GRIDINDEX_NOGHOST(i+1,j,k-1)] -
		 BaryonField[GPotNum][GRIDINDEX_NOGHOST(i-1,j,k+1)])/4; 
	
	//	  MaxDensity[index] = max(max(rhox, max(rhoy, rhoz)), tiny_number)/ GravitationalConstant/CellWidthSquared;
	MaxDensity[index] = max(rhoxz, max(rhoyz, max(rhoxy,max(max(rhox, max(rhoy, rhoz)), tiny_number))))
	  / GravitationalConstant/CellWidthSquared;
	
	// eigenvalues of tidal tensor instead of maximal component
#if 0
	double m ,c0,c1,sqrt_p,w1,w2,w0, phi,p,q,c,s;
	m = rhox+rhoy+rhoz;
	c1 = rhox*rhoy+rhox*rhoz+rhoy*rhoz - rhoxy*rhoxy - rhoyz*rhoyz - rhoxz*rhoxz;
	c0 = rhoz*rhoxy*rhoxy + rhox*rhoyz*rhoyz + rhoy*rhoxz*rhoxz -
	  rhox*rhoy*rhoz - 2.0*rhoxz*rhoxy*rhoyz;
	p = m*m - 3.0*c1;
	q = m*(p - (3.0/2.0)*c1) - (27.0/2.0)*c0;
	sqrt_p = sqrt(fabs(p));

	phi = 27.0 * ( 0.25*c1*c1*(p - c1) + c0*(q + 27.0/4.0*c0));
	phi = (1.0/3.0) * atan2(sqrt(fabs(phi)), q);
	
#define M_SQRT3    1.73205080756887729352744634151   // sqrt(3)	  	  
	c = sqrt_p*cos(phi);
	q  s = (1.0/M_SQRT3)*sqrt_p*sin(phi);
	
	w1  = (1.0/3.0)*(m - c);
	w2  = w1 + s;
	w0  = w1 + c;
	w1 -= s;
	
	printf("%g %g \n",  max(w2, max( w1, max(w0, tiny_number))),
	       max(rhoxz, max(rhoyz, max(rhoxy, max(max(rhox, max(rhoy, rhoz)), tiny_number)))));
	// using maximal eigenvalue for density estimate
	MaxDensity[index] = max(w2, max( w1, max(w0, tiny_number)))
	  / GravitationalConstant/CellWidthSquared;
#endif // use eigenvalues 
	
	rhomax = max(rhomax,BaryonField[DensNum][index]);
	maxmax = max(maxmax, MaxDensity[index]);
	
	BaryonField[NumberOfBaryonFields-1][index] = MaxDensity[index]; // for debugging
	MaxDensity[index] = max(MaxDensity[index]  , BaryonField[DensNum][index] );
	
	
	//	printf("Density:  %g  MaxDensity: %g  ratio: %g \n", BaryonField[DensNum][index],
	//	       MaxDensity[index], MaxDensity[index]/BaryonField[DensNum][index]);
	
      }
    };

  //  printf("rhomax %g , maxmax: %g ratio:: %g\n", rhomax,maxmax, maxmax/rhomax);

  /* Get density units. */
 
  float DensityUnits=1, LengthUnits=1, VelocityUnits=1, TimeUnits=1,
    TemperatureUnits=1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.\n");
  }
 
  /* Compute constant for Jean's length computation.
      l_j = sqrt((pi*k*T) / (G \rho m_p))  . */
 
  FLOAT JLSquared = (double(3.14159*1.38e-16/6.67e-8)/
                 (double(DensityUnits)*double(1.67e-24))) /
                (double(LengthUnits)*double(LengthUnits));
 
  if (ProblemType == 60 || ProblemType == 61)
    JLSquared = double(4.0*3.14159*3.14159)/GravitationalConstant; //AK

  if (EOSType > 0)
    {
      float cs,dpdrho,dpde, eint, h, rho, p;
      EOS(p, rho, eint, h, cs, dpdrho, dpde, EOSType, 1) ;
      JLSquared = cs*cs*M_PI/GravConst/DensityUnits*VelocityUnits*VelocityUnits/LengthUnits/LengthUnits; // TA
    }
  //  printf("JLSquared %g\n", JLSquared);

  /* This is the safety factor to decrease the Jean's length by. */
 
  JLSquared /= POW(RefineByJeansLengthSafetyFactor, 2);
 
/* printf("jl: JL, dx, t, d = %"GSYM" %"GSYM" %"GSYM" %"GSYM"\n", sqrt(JLSquared), CellWidth[0][3],
	 temperature[(3 + 3*GridDimension[1])*GridDimension[0]+3],
	 BaryonField[DensNum][(3 + 3*GridDimension[1])*GridDimension[0]+3]);*/
 
  /* Loop over grid. */

  for (k = GridStartIndex[2]+1; k < GridEndIndex[2]; k++)
    for (j = GridStartIndex[1]+1; j < GridEndIndex[1]; j++) {
      for (i = GridStartIndex[0]+1; i < GridEndIndex[0]; i++) {
	index = GRIDINDEX_NOGHOST(i,j,k);
	if (EOSType == 0) {
	  if (CellWidthSquared > JLSquared*temperature[index]/MaxDensity[index])
	    FlaggingField[index]++; 
	}
        // isothermal and polytropic sound speed version
	if (EOSType > 0)
	  if (CellWidthSquared > JLSquared/MaxDensity[index])
	    FlaggingField[index]++; 
      }
    }

  /* clean up */
 
  if (ProblemType != 60 && ProblemType != 61 && EOSType == 0 ) //AK
    delete temperature;
 
  delete MaxDensity;

  /* Count number of flagged Cells. */
 
  int NumberOfFlaggedCells = 0;
  for (i = 0; i < size; i++) {
    FlaggingField[i] = (FlaggingField[i] >= 1)? 1 : 0;
    NumberOfFlaggedCells += FlaggingField[i];
  }
  //  printf("NumberOfFlaggedCells %i\n", NumberOfFlaggedCells);
 
  return NumberOfFlaggedCells;
 
}
