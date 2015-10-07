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

	/* Determine if we're using metallicity (as a color field) */
	int MetalNum = FindField(Metallicity,BoundaryFieldType,NumberOfBaryonFields);
	int UseMetallicityField = (MetalNum == -1) ? 0 : 1;
		
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

		// update bndry type (needed for restart runs)
		for( int field = 0 ; field < NumberOfBaryonFields; ++field ){
			BoundaryType[field][dim][0][index] = inflow;  // left bnd
			BoundaryType[field][dim][1][index] = outflow; // right bnd
		}
 
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
		static int hasNotArrived = 1;
		if( hasNotArrived ) 
			deltime = time - distance/GalaxySimulationRPSWindShockSpeed - GalaxySimulationRPSWindDelay;
		else
			deltime = time - distance/vMag - GalaxySimulationRPSWindDelay; // fluid travels at bulk speed
		if( hasNotArrived && deltime > 0 ) hasNotArrived = 0;

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

				char filename[] = "ICMinflow_data.in"; char line[MAX_LINE_LENGTH]; FILE *fptr;
				if ((fptr = fopen(filename, "r")) == NULL) ENZO_FAIL("Erroring opening ICMinflow_data.in");

				int f_index = 0;
				while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {
					if (line[0] == 'N') {
						(sscanf(line, "NumberOfSteps = %"ISYM"\n", &ICMTableSize));
						ICMTimeTable        = new float[ICMTableSize];
						ICMDensityTable     = new float[ICMTableSize];
						ICMTotalEnergyTable = new float[ICMTableSize];
						ICMVelocityXTable   = new float[ICMTableSize];
						ICMVelocityYTable   = new float[ICMTableSize];
						ICMVelocityZTable   = new float[ICMTableSize];
					}
					if (line[0] != '#' && line[0] != 'N') {
						
						// read values
						sscanf(line,"%"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM" %"FSYM"",&ICMTimeTable[f_index],&ICMDensityTable[f_index], 
						       &temperature,&ICMVelocityXTable[f_index],&ICMVelocityYTable[f_index],&ICMVelocityZTable[f_index]);
						
						// convert to code units
						ICMTimeTable[f_index]        /= TimeUnits;
						ICMDensityTable[f_index]     /= DensityUnits;
						ICMVelocityXTable[f_index]   /= (LengthUnits/TimeUnits); 
						ICMVelocityYTable[f_index]   /= (LengthUnits/TimeUnits); 
						ICMVelocityZTable[f_index]   /= (LengthUnits/TimeUnits); 
						ICMTotalEnergyTable[f_index] = temperature/TemperatureUnits/((Gamma-1.0)*0.6);

						if (HydroMethod != 2) {
							ICMTotalEnergyTable[f_index] += 0.5*(   pow(ICMVelocityXTable[0],2)
							                                      + pow(ICMVelocityYTable[1],2)
							                                      + pow(ICMVelocityZTable[2],2));
						}

						f_index++;
					} // end data-line if
      	} // end file read while

				fclose(fptr); 
				loadTable = 0;
			}// end load table if

			if(deltime < 0.0){

				// use default pre-wind values
				BoundaryValue[DensNum][dim][0][index] = GalaxySimulationPreWindDensity;
				BoundaryValue[TENum  ][dim][0][index] = GalaxySimulationPreWindTotalEnergy;
				BoundaryValue[Vel1Num][dim][0][index] = GalaxySimulationPreWindVelocity[0];
				if (BoundaryRank > 1)
					BoundaryValue[Vel2Num][dim][0][index] = GalaxySimulationPreWindVelocity[1];
				if (BoundaryRank > 2)
					BoundaryValue[Vel3Num][dim][0][index] = GalaxySimulationPreWindVelocity[2];

			} else {

				/* interpolate w/ lookup table */

				float t_ratio,v1,v2;
				int i1,i2=-1;
				
				// find times that bracket what we want
				while( ++i2 < ICMTableSize ) if( deltime < ICMTimeTable[i2] ) break;
				i1 = i2-1;
				
				if(i2<ICMTableSize)
					t_ratio = (deltime - ICMTimeTable[i2])/(ICMTimeTable[i1] - ICMTimeTable[i2]);
				else { // if beyond final time, just use final val
					i2--;
					t_ratio = 1.0;
				}

			
				BoundaryValue[DensNum][dim][0][index] = t_ratio*ICMDensityTable[i1]
				                                        + (1.0-t_ratio)*ICMDensityTable[i2];
				BoundaryValue[TENum  ][dim][0][index] = t_ratio*ICMTotalEnergyTable[i1] 
				                                        + (1.0-t_ratio)*ICMTotalEnergyTable[i2];
				BoundaryValue[Vel1Num][dim][0][index] = t_ratio*ICMVelocityXTable[i1]
				                                        + (1.0-t_ratio)*ICMVelocityXTable[i2];
				if (BoundaryRank > 1)
					BoundaryValue[Vel2Num][dim][0][index] = t_ratio*ICMVelocityYTable[i1]
				 	                                        + (1.0-t_ratio)*ICMVelocityYTable[i2];
				if (BoundaryRank > 2)
					BoundaryValue[Vel3Num][dim][0][index] = t_ratio*ICMVelocityZTable[i1]
						                                      + (1.0-t_ratio)*ICMVelocityZTable[i2];

				// update RPS Wind Vector for time delay calc
				if( index == 0.0 ){
					GalaxySimulationRPSWindVelocity[0] = BoundaryValue[Vel1Num][dim][0][index];
					GalaxySimulationRPSWindVelocity[1] = BoundaryValue[Vel2Num][dim][0][index];
					GalaxySimulationRPSWindVelocity[2] = BoundaryValue[Vel3Num][dim][0][index];
				}
	
			}

		} else {
			ENZO_FAIL("Error in ExternalBoundary_SetGalaxyBoundary: GalaxySimulationRPSWind choice invalid");
		}

		// update metallicity field
		if( UseMetallicityField )
			BoundaryValue[MetalNum][dim][0][index] = 1.0e-10;

		if( BoundaryValue[DensNum][dim][0][index] < 0.0 ) 
			ENZO_FAIL("Error in ExternalBoundary_SetGalaxyBoundary: Negative Density");
		if( BoundaryValue[TENum][dim][0][index] < 0.0 ) 
			ENZO_FAIL("Error in ExternalBoundary_SetGalaxyBoundary: Negative Total Energy");

    if( BoundaryValue[DensNum][dim][0][index] != BoundaryValue[DensNum][dim][0][index] )  
      ENZO_FAIL("Error in ExternalBoundary_SetGalaxyBoundary: Density NaN");
    if( BoundaryValue[TENum][dim][0][index] != BoundaryValue[TENum][dim][0][index] )  
      ENZO_FAIL("Error in ExternalBoundary_SetGalaxyBoundary: Total Energy NaN");


	} // end loop over boundary slice
	
	} // end loop over boundary directions

  return SUCCESS;
 
}
