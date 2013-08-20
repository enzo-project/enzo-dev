/***********************************************************************
/
/  EXTERNAL BOUNDARY CLASS (SETS OUTFLOW BOUNDARY CONDITIONS FOR GAL SIM)
/
/  written by: Greg Bryan
/  date:       May, 1995
/  modified1:  Munier Salem
/  date:       August, 2013
/
/  PURPOSE:
/
/  RETURNS: SUCCESS or FAIL
/
************************************************************************/
 
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

/* function prototypes */
 
int FindField(int f, int farray[], int n);

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, float *MassUnits, FLOAT time);
 
// Set the Left BoundaryValue of the chosen wave direction (set by
//  GalaxySimulationRPSWindSpeed) to the appropriate inflow boundary condition.
 
int ExternalBoundary::SetGalaxySimulationBoundary(FLOAT time)
{
	if( 0 == GalaxySimulationRPSWind ) return SUCCESS;

	if( MyProcessorNumber == ROOT_PROCESSOR ){
	  fprintf(stderr,"GalaxySimulationRPSWindShockSpeed = %"GSYM"\n",GalaxySimulationRPSWindShockSpeed);

	  fprintf(stderr,"GalaxySimulationRPSWindDensity = %"GSYM"\n",GalaxySimulationRPSWindDensity);
	  fprintf(stderr,"GalaxySimulationRPSWindTotalEnergy = %"GSYM"\n",GalaxySimulationRPSWindTotalEnergy);
	  fprintf(stderr,"GalaxySimulationRPSWindVelocity = %"GSYM", %"GSYM", %"GSYM"\n",
			GalaxySimulationRPSWindVelocity[0],GalaxySimulationRPSWindVelocity[1], GalaxySimulationRPSWindVelocity[2]);

	  fprintf(stderr,"GalaxySimulationPreWindDensity = %"GSYM"\n",GalaxySimulationPreWindDensity);
	  fprintf(stderr,"GalaxySimulationPreWindTotalEnergy = %"GSYM"\n",GalaxySimulationPreWindTotalEnergy);
	  fprintf(stderr,"GalaxySimulationPreWindVelocity = %"GSYM", %"GSYM", %"GSYM"\n",
			GalaxySimulationPreWindVelocity[0],GalaxySimulationPreWindVelocity[1],GalaxySimulationPreWindVelocity[2]);
	} // end if

  /* declarations */

  int i, j, dim, index;
  int NumberOfZones[MAX_DIMENSION], Offset[MAX_DIMENSION];
  float deltime, distance, pos[MAX_DIMENSION];
  const float TwoPi = 6.283185;
 
  /* Compute size of entire mesh. */
 
  int size = 1;
  for (dim = 0; dim < BoundaryRank; dim++)
    size = size*BoundaryDimension[dim];
 
  /* Find fields: density, total energy, velocity1-3. */
 
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
					 Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.\n");
  }

  /* set the appropriate BoundaryValues on the left side */
 
  for (dim = 0; dim < BoundaryRank; dim++)
    if (BoundaryDimension[dim] != 1) {
 
      /* If the BoundaryValue fields are missing, create them. */
 
      for (int field = 0; field < NumberOfBaryonFields; field++)
	if (BoundaryValue[field][dim][0] == NULL)
	  BoundaryValue[field][dim][0] =
	    new float[size/BoundaryDimension[dim]];
 
      /* Compute quantities needed for boundary face loop (below). */
 
      int dim1, dim2;
      dim1 = (dim == 0) ? 1 : 0;
      dim2 = dim1 + 1;
      dim2 = (dim2 == dim) ? dim2+1 : dim2;
      for (i = 0; i < 3; i++) {
	NumberOfZones[i] = max(BoundaryDimension[i] - 2*NumberOfGhostZones,1);
	Offset[i]        = min(NumberOfGhostZones, BoundaryDimension[i]) - 1;
      }
      pos[dim] = 0.0;
 
      /* Loop over the boundary face. */
 
      for (i = 0; i < BoundaryDimension[dim1]; i++)
	for (j = 0; j < BoundaryDimension[dim2]; j++) {
 
	  /* Compute the index into the boundary value. */
 
	  index = j*BoundaryDimension[dim1] + i;
 
	  /* Find the 3D vector from the corner to the current location. */
 
	  pos[dim1] = (float(i-Offset[dim1]))*
	    (DomainRightEdge[dim1]-DomainLeftEdge[dim1]) /
	      float(NumberOfZones[dim1]);
	  pos[dim2] = (float(j-Offset[dim2]))*
	    (DomainRightEdge[dim2]-DomainLeftEdge[dim2]) /
	      float(NumberOfZones[dim2]);
 
	  /* Compute the distance along the wave propogation vector
     *    |d| = |v_wind . x|/|v_wind| */ 		

		float vMag = sqrt( POW(GalaxySimulationRPSWindVelocity[0],2.0) +
                       POW(GalaxySimulationRPSWindVelocity[1],2.0) +
                       POW(GalaxySimulationRPSWindVelocity[2],2.0) );
		distance = 0.0;
		if( vMag > 0.0 )
			distance = fabs( GalaxySimulationRPSWindVelocity[0]*pos[0] +
			                 GalaxySimulationRPSWindVelocity[1]*pos[1] +
			                 GalaxySimulationRPSWindVelocity[2]*pos[2] )/vMag;

	  /* Find the difference between the current time and the time at
	     which the wave will reach this point. */
 
	  deltime = time - distance/GalaxySimulationRPSWindShockSpeed - GalaxySimulationRPSWindDelay;

		if( 1 == GalaxySimulationRPSWind ){

			/* Update bounds with simple shock wind */

			if (deltime > 0.0) {  // Shock has arrived, set post-shock values
				BoundaryValue[DensNum][dim][0][index] = GalaxySimulationRPSWindDensity;
				BoundaryValue[TENum][dim][0][index] = GalaxySimulationRPSWindTotalEnergy;
				BoundaryValue[Vel1Num][dim][0][index] = GalaxySimulationRPSWindVelocity[0];
				if (BoundaryRank > 1)
					BoundaryValue[Vel2Num][dim][0][index] = GalaxySimulationRPSWindVelocity[1];
				if (BoundaryRank > 2)
					BoundaryValue[Vel3Num][dim][0][index] = GalaxySimulationRPSWindVelocity[2];
			} else { // If not, set pre-shock values
				BoundaryValue[DensNum][dim][0][index] = GalaxySimulationPreWindDensity;
				BoundaryValue[TENum][dim][0]  [index] = GalaxySimulationPreWindTotalEnergy;
				BoundaryValue[Vel1Num][dim][0][index] = GalaxySimulationPreWindVelocity[0];
				if (BoundaryRank > 1)
					BoundaryValue[Vel2Num][dim][0][index] = GalaxySimulationPreWindVelocity[1];
				if (BoundaryRank > 2)
					BoundaryValue[Vel3Num][dim][0][index] = GalaxySimulationPreWindVelocity[2];
			}

		} else if( 2 == GalaxySimulationRPSWind ){

			/* Update Bounds w/ table of density and wind velocity components */
	
			static double *ICMDensityTable, *ICMTotalEnergyTable, *ICMVelocityXTable, 
				*ICMVelocityYTable, *ICMVelocityZTable, *ICMTimeTable;	
			static int loadTable = 1,ICMTableSize=0; // only when processor revs up
			if( loadTable ){
	
				/* Find units */
				float DensityUnits,LengthUnits,TemperatureUnits,TimeUnits,VelocityUnits,MassUnits,temperature;
				GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
				         &TimeUnits, &VelocityUnits, &MassUnits, time);

				char filename[] = "ICMinflow_data.out"; char line[MAX_LINE_LENGTH]; FILE *fptr;
				if ((fptr = fopen(filename, "r")) == NULL) ENZO_FAIL("Erroring opening ICMinflow_data.out");

				int index = 0;
				while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
					if (line[0] == 'N') {
						(sscanf(line, "NumberOfSteps = %d\n", &ICMTableSize));
						ICMTimeTable        = new float[ICMTableSize];
						ICMDensityTable     = new float[ICMTableSize];
						ICMTotalEnergyTable = new float[ICMTableSize];
						ICMVelocityXTable   = new float[ICMTableSize];
						ICMVelocityYTable   = new float[ICMTableSize];
						ICMVelocityZTable   = new float[ICMTableSize];
					}
					if (line[0] != '#' && line[0] != 'N') {
						
						// read values
						sscanf(line,"%"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM"",&ICMTimeTable[index],&ICMDensityTable[index], 
						       &temperature,&ICMVelocityXTable[index],&ICMVelocityYTable[index],&ICMVelocityZTable[index]);
						
						// convert to code units
						ICMTimeTable[index]        /= TimeUnits;
						ICMDensityTable[index]     /= DensityUnits;
						ICMVelocityXTable[index]   /= (LengthUnits/TimeUnits); 
						ICMVelocityYTable[index]   /= (LengthUnits/TimeUnits); 
						ICMVelocityZTable[index]   /= (LengthUnits/TimeUnits); 
						ICMTotalEnergyTable[index] = temperature/TemperatureUnits/((Gamma-1.0)*0.6);

						if (HydroMethod != 2) {
							ICMTotalEnergyTable[index] += 0.5*(   pow(ICMVelocityXTable[0],2)
							                                    + pow(ICMVelocityYTable[1],2)
							                                    + pow(ICMVelocityZTable[2],2));
						}

						index++;
					} // end data-line if
      	} // end file read while

				fclose(fptr); 

				// display table
				fprintf(stderr,"LOOKUP TABLE:\n");
				for( int i = 0 ; i < ICMTableSize ; i++ )
					fprintf(stderr,"%"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM"\n",ICMTimeTable[i],ICMDensityTable[i],
						ICMTotalEnergyTable[i],ICMVelocityXTable[i],ICMVelocityYTable[i],ICMVelocityZTable[i]);

				loadTable = 0;
			}// end load table if

/*
			if(deltime > 0.0){ // use default pre-wind values
				BoundaryValue[DensNum][dim][0][index] = GalaxySimulationPreWindDensity;
				BoundaryValue[TENum][dim][0]  [index] = GalaxySimulationPreWindTotalEnergy;
				BoundaryValue[Vel1Num][dim][0][index] = GalaxySimulationPreWindVelocity[0];
				if (BoundaryRank > 1)
					BoundaryValue[Vel2Num][dim][0][index] = GalaxySimulationPreWindVelocity[1];
				if (BoundaryRank > 2)
					BoundaryValue[Vel3Num][dim][0][index] = GalaxySimulationPreWindVelocity[2];
			} else {
				// use table and interpolate!
			}
*/

		} else {
			ENZO_FAIL("Error in ExternalBoundary_SetGalaxyBoundary: GalaxySimulationRPSWind choice invalid");
		}

	} // end loop over boundary slice
	
	} // end loop over boundary directions

  return SUCCESS;
 
}
