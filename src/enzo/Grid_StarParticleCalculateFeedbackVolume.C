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
							float TemperatureUnits, float TimeUnits)
{

	const float WhalenMaxVelocity = 35; // km/s

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







	/***********************************************************************
                                Checking volume--The real volume the supernova
								will be put into.
  ************************************************************************/

	float outerRadius2=1.2*1.2*radius*radius;
    FLOAT depositedVolume = 0.0;
    int CellsModified = 0;
		for (k = 0; k < GridDimension[2]; k++)
		{

			delz = CellLeftEdge[2][k] + 0.5 * CellWidth[2][k] - cstar->ReturnPosition()[2];
			sz = sign(delz);
			delz = fabs(delz);
			delz = min(delz, DomainWidth[2] - delz);

			for (j = 0; j < GridDimension[1]; j++)
			{

				dely = CellLeftEdge[1][j] + 0.5 * CellWidth[1][j] - cstar->ReturnPosition()[1];
				sy = sign(dely);
				dely = fabs(dely);
				dely = min(dely, DomainWidth[1] - dely);

				index = (k * GridDimension[1] + j) * GridDimension[0];
				for (i = 0; i < GridDimension[0]; i++, index++)
				{

					delx = CellLeftEdge[0][i] + 0.5 * CellWidth[0][i] - cstar->ReturnPosition()[0];
					sx = sign(delx);
					delx = fabs(delx);
					delx = min(delx, DomainWidth[0] - delx);

					radius2 = delx * delx + dely * dely + delz * delz;
					if (radius2 <= outerRadius2)
					{
						depositedVolume += CellWidth[0][i]*CellWidth[1][j]*CellWidth[2][k];
                        CellsModified ++;
					} // END if inside radius
				}	 // END i-direction
			}		  // END j-direction
		}			  // END k-direction
        return CellsModified;
}