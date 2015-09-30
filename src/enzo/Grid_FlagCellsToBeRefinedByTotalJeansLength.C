/***********************************************************************
/
/  GRID CLASS (FLAG CELLS TO BE REFINED BY THE JEAN'S CRITERION)
/
/  written by: Tom Abel 
/  date:       October 2010
/  modified1:  
/
/  PURPOSE: Refine on the Jeans length determined from the inertial tensor
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
 
extern "C" void FORTRAN_NAME(smooth2)(float *source, float *dest, int *ndim,
				      int *sdim1, int *sdim2, int *sdim3);

int grid::FlagCellsToBeRefinedByTotalJeansLength()
{
  /* declarations */
 
  int i, dim;

  /* error check */
 
  if (FlaggingField == NULL) {
    fprintf(stderr, "Flagging Field is undefined.\n");
    return -1;
  }

  if (NumberOfBaryonFields == 0) 
    return SUCCESS;
 
  /* compute size */
 
  int size = 1;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
 
  /* Compute the temperature field. */
 
  float *temperature = NULL;
  if (ProblemType != 60 && ProblemType != 61 && EOSType == 0) { //AK
    temperature = new float[size];
    if (JeansRefinementColdTemperature > 0.0) {
      for (i = 0; i < size; i++)
	temperature[i] = JeansRefinementColdTemperature;
    } else {
      if (this->ComputeTemperatureField(temperature) == FAIL)
	ENZO_FAIL("Error in grid->ComputeTemperature.");
      for (i = 0; i < size; i++) 
	temperature[i] = max(JeansRefinementColdTemperature, temperature[i]);
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
  for (i=0;i<size;i++) {
    MaxDensity[i] = BaryonField[DensNum][i];
    BaryonField[NumberOfBaryonFields-1][i]  = MaxDensity[i];
  }
  
  FLOAT CellWidthSquared = CellWidth[0][0]*CellWidth[0][0];
  int j, jj, k, index;
  float rhox, rhoy, rhoz, rhoxy, rhoxz, rhoyz,  rhomax, maxmax;

  if ((GPotNum = FindField(GravPotential, FieldType, NumberOfBaryonFields)) < 0) {
    ENZO_FAIL("Cannot find Gravitational Potential. Set WritePotential = 1 ... hack");
  };
  
  // make a copy
  float *Phi = NULL;
  Phi = new float[size];
  for (i=0;i<size;i++) 
    Phi[i] = BaryonField[GPotNum][i];
  //    Phi[i] = BaryonField[GPotNum][i];
  

  int Off[3];
  for (int dim = 0; dim < GridRank; dim++)
    Off[dim] = (GravitatingMassFieldDimension[dim] - GridDimension[dim])/2;
  
  jj = 0;
  if (PotentialField != NULL)
    for (k = 0; k < GridDimension[2]; k++)
      for (j = 0; j < GridDimension[1]; j++) {
	index = (((k+Off[2])*GravitatingMassFieldDimension[1]) + (j+Off[1]))*GravitatingMassFieldDimension[0] + Off[0] ;
	for (i = 0; i < GridDimension[0]; i++, index++)
	  Phi[jj++] = PotentialField[index];
      }


  //  FORTRAN_NAME(smooth2)(BaryonField[GPotNum] ,Phi, &GridRank, GridDimension, GridDimension+1, GridDimension+2);
      
  //  printf("GpotNum : %i\n", GPotNum);
  rhox = rhoy = rhoz = 0.;
  int ci,cj,ck,cind;
  for (k = GridStartIndex[2]+1; k < GridEndIndex[2]; k++) {
    for (j = GridStartIndex[1]+1; j < GridEndIndex[1]; j++) {
      for (i = GridStartIndex[0]+1; i < GridEndIndex[0]; i++) {

	index = GRIDINDEX_NOGHOST(i,j,k);

	ci = i; cj= j; ck = k; cind = index; //  <-  would do this if BC would be correct

	rhox = (Phi[GRIDINDEX_NOGHOST(ci+1,cj,ck)] 
		+ Phi[GRIDINDEX_NOGHOST(ci-1,cj,ck)] - 2.*Phi[cind]);
	rhoy = (Phi[GRIDINDEX_NOGHOST(ci,cj+1,ck)]
		+ Phi[GRIDINDEX_NOGHOST(ci,cj-1,ck)] - 2.*Phi[cind]);
	rhoz = (Phi[GRIDINDEX_NOGHOST(ci,cj,ck+1)] 
		+ Phi[GRIDINDEX_NOGHOST(ci,cj,ck-1)] - 2.*Phi[cind]);
	rhoxy = (Phi[GRIDINDEX_NOGHOST(ci+1,cj+1,ck)] +
		 Phi[GRIDINDEX_NOGHOST(ci-1,cj-1,ck)] -
		 Phi[GRIDINDEX_NOGHOST(ci+1,cj-1,ck)] -
		 Phi[GRIDINDEX_NOGHOST(ci-1,cj+1,ck)])/4; 
	rhoyz = (Phi[GRIDINDEX_NOGHOST(ci,cj+1,ck+1)] +
		 Phi[GRIDINDEX_NOGHOST(ci,cj-1,ck-1)] - 
		 Phi[GRIDINDEX_NOGHOST(ci,cj-1,ck+1)] -
		 Phi[GRIDINDEX_NOGHOST(ci,cj+1,ck-1)])/4; 
	rhoxz = (Phi[GRIDINDEX_NOGHOST(ci+1,cj,ck+1)] +
		 Phi[GRIDINDEX_NOGHOST(ci-1,cj,ck-1)] +
		 Phi[GRIDINDEX_NOGHOST(ci+1,cj,ck-1)] -
		 Phi[GRIDINDEX_NOGHOST(ci-1,cj,ck+1)])/4; 
	
	double det = 0.;
	det = rhox*rhoy*rhoz + rhoxy*rhoyz*rhoxz + rhoxz*rhoxy*rhoyz 
	   - rhoxz*rhoy*rhoxz - rhoxy*rhoxy*rhoz - rhox*rhoyz*rhoyz; 
	//	MaxDensity[index] = max(max(rhox, max(rhoy, rhoz)), tiny_number)/ GravitationalConstant/CellWidthSquared;
	if (det > 0.) 
	  MaxDensity[index] = pow(det, 0.33334)/ GravitationalConstant/CellWidthSquared;

	// divergence without negative values
	//	MaxDensity[index] = (max(rhox, 0)+max(rhoy, 0)+max(rhoz, 0))/ GravitationalConstant/CellWidthSquared ;
	
	// largest component
	//MaxDensity[index] = max(rhoxz, max(rhoyz, max(rhoxy,max(max(rhox, max(rhoy, rhoz)),tiny_number))))/ GravitationalConstant/CellWidthSquared ;
	
	BaryonField[NumberOfBaryonFields-1][index] = MaxDensity[index]; // for debugging copy into Debug field	
	//	BaryonField[NumberOfBaryonFields-1][index] = Phi[index]; // for debugging copy into Debug field	
	MaxDensity[index] = max(MaxDensity[index], BaryonField[DensNum][index]);
      }
    }
  };
  
  delete [] Phi;
  
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
    };
  //  printf("JLSquared %g\n", JLSquared);

  /* This is the safety factor to decrease the Jean's length by. */
 
  JLSquared /= POW(RefineByJeansLengthSafetyFactor, 2);
 
  /* printf("jl: JL, dx, t, d = %"GSYM" %"GSYM" %"GSYM" %"GSYM"\n", sqrt(JLSquared), CellWidth[0][3],
     temperature[(3 + 3*GridDimension[1])*GridDimension[0]+3],
     BaryonField[DensNum][(3 + 3*GridDimension[1])*GridDimension[0]+3]);*/
 
  /* Loop over grid. */

  for (k = GridStartIndex[2]; k < GridEndIndex[2]; k++) {
    for (j = GridStartIndex[1]; j < GridEndIndex[1]; j++) {
      for (i = GridStartIndex[0]; i < GridEndIndex[0]; i++) {
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
    };
  };

  /* clean up */
  delete [] MaxDensity ;
  
  if (ProblemType != 60 && ProblemType != 61 && EOSType == 0 ) //AK
    delete [] temperature;
  
  
  /* Count number of flagged Cells. */
  
  int NumberOfFlaggedCells = 0;
  for (i = 0; i < size; i++) {
    FlaggingField[i] = (FlaggingField[i] >= 1)? 1 : 0;
    NumberOfFlaggedCells += FlaggingField[i];
  };
  //  printf("NumberOfFlaggedCells %i\n", NumberOfFlaggedCells);
  
  return NumberOfFlaggedCells;
  
}
