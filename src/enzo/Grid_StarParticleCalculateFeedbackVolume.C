  #include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "list.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "CosmologyParameters.h"
#include "phys_constants.h"



int grid::StarParticleCalculateFeedbackVolume(Star *cstar, int level, float radius, float DensityUnits,
							float LengthUnits, float VelocityUnits,
							float TemperatureUnits, float TimeUnits, int &nCells, float &depositedMass,
							float &depositedMetal, FLOAT &depositedVolume)
{


	int dim, i, j, k, index;
	int sx, sy, sz;
	FLOAT delx, dely, delz, radius2, Radius, DomainWidth[MAX_DIMENSION];

	if (MyProcessorNumber != ProcessorNumber)
		return 0.0;

	/* If the radius is less than the cell width, return */

	if (radius < CellWidth[0][0])
		return 0.0;

	/* Check if sphere overlaps with this grid */

	for (dim = 0; dim < GridRank; dim++)
		if (cstar->pos[dim] - radius > GridRightEdge[dim] ||
			cstar->pos[dim] + radius < GridLeftEdge[dim])
			return 0.0;

	for (dim = 0; dim < GridRank; dim++)
		DomainWidth[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];

	int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
	if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
										 Vel3Num, TENum) == FAIL)
	{
		ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
	}

	/* Find Multi-species fields. */

	int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
		DINum, DIINum, HDINum;
	if (MultiSpecies)
		if (this->IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum,
										HeIIINum, HMNum, H2INum, H2IINum, DINum,
										DIINum, HDINum) == FAIL)
		{
			ENZO_FAIL("Error in grid->IdentifySpeciesFields.");
		}


	/* Find Metallicity or SNColour field and set flag. */

	int SNColourNum, MetalNum, Metal2Num, MBHColourNum, Galaxy1ColourNum,
		Galaxy2ColourNum, MetalIaNum, MetalIINum;
	int MetallicityField = FALSE;

	if (this->IdentifyColourFields(SNColourNum, Metal2Num, MetalIaNum,
								   MetalIINum, MBHColourNum, Galaxy1ColourNum,
								   Galaxy2ColourNum) == FAIL)
		ENZO_FAIL("Error in grid->IdentifyColourFields.\n");

	MetalNum = max(Metal2Num, SNColourNum);
	MetallicityField = (MetalNum > 0) ? TRUE : FALSE;
	if (MetalNum > 0 && SNColourNum > 0 && cstar->type == PopIII)
		MetalNum = SNColourNum;





	/***********************************************************************
                                Checking volume--The real volume the supernova
								will be put into.
  ************************************************************************/

	float outerRadius2=radius*radius;
	int GZ = NumberOfGhostZones;
		for (k = GZ; k < GridDimension[2]; k++)
		{

			delz = CellLeftEdge[2][k] + 0.5 * CellWidth[2][k] - cstar->ReturnPosition()[2];
			sz = sign(delz);
			delz = fabs(delz);
			delz = min(delz, DomainWidth[2] - delz);

			for (j = GZ; j < GridDimension[1]; j++)
			{

				dely = CellLeftEdge[1][j] + 0.5 * CellWidth[1][j] - cstar->ReturnPosition()[1];
				sy = sign(dely);
				dely = fabs(dely);
				dely = min(dely, DomainWidth[1] - dely);

				index = (k * GridDimension[1] + j) * GridDimension[0];
				for (i = GZ; i < GridDimension[0]; i++, index++)
				{

					delx = CellLeftEdge[0][i] + 0.5 * CellWidth[0][i] - cstar->ReturnPosition()[0];
					sx = sign(delx);
					delx = fabs(delx);
					delx = min(delx, DomainWidth[0] - delx);

					radius2 = delx * delx + dely * dely + delz * delz;
					if (radius2 <= outerRadius2)
					{
						depositedVolume += CellWidth[0][i]*CellWidth[1][j]*CellWidth[2][k];
                        nCells ++;
						depositedMass += BaryonField[DensNum][index]*CellWidth[0][i]*CellWidth[1][j]*CellWidth[2][k];
						depositedMetal += BaryonField[SNColourNum][index]*CellWidth[0][i]*CellWidth[1][j]*CellWidth[2][k];
					} // END if inside radius
				}	 // END i-direction
			}		  // END j-direction
		}			  // END k-direction
		return SUCCESS;
}