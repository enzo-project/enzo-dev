/***********************************************************************
/
/  Set flagging field to refine around pop III particles. 
/
/  written by: Azton Wells
/  date:       Dec, 2019
/  modified1:
/
/  PURPOSE: 
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif
#include <stdlib.h>
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
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "CommunicationUtilities.h"
#include "phys_constants.h"
int GetUnits(float *DensityUnits, float *LengthUnits,
             float *TemperatureUnits, float *TimeUnits,
             float *VelocityUnits, FLOAT Time);
/* 
    New refinement criteria will refine the region around Pop III star particles.  
    This will ensure that their thermal feedback is effective, given that MaximumRefinementLevel is 
    sufficient
 */

int grid::FlagCellsToBeRefinedByPopIII(int level)
{
    /* declarations */

    float LengthUnits, TimeUnits, TemperatureUnits, VelocityUnits,
        DensityUnits;
    int i, j, k, index, dim, size = 1, NumberOfFlaggedCells = 0;
    int Start[MAX_DIMENSION], End[MAX_DIMENSION];

    FLOAT CellSize, xpos, ypos, zpos;
    /* Return if this grid is not on this processor. */

    if (MyProcessorNumber != ProcessorNumber)
        return SUCCESS;

    /* error check */
    if (FlaggingField == NULL)
    {
        fprintf(stderr, "Flagging Field is undefined.\n");
        return -1;
    }

    /* loop over dimensions - I guess this is unnecessary, 
     but it's handy to have shorter names */
    for (dim = 0; dim < MAX_DIMENSION; dim++)
    {
        Start[dim] = GridStartIndex[dim];
        End[dim] = GridEndIndex[dim];
    }

    /* compute size */
    for (int dim = 0; dim < GridRank; dim++)
        size *= GridDimension[dim];
    GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
             &TimeUnits, &VelocityUnits, Time);

    //If A: no particles or 
    // B: the resolution is such that its already as refined as PopIIISupernovaMustRefineResolution 
    int nCellsPerRadius = 1.1*float(PopIIISupernovaRadius)/(CellWidth[0][0]*LengthUnits/pc_cm)+0.5;
    if (NumberOfParticles == 0 || 2*nCellsPerRadius > PopIIISupernovaMustRefineResolution || level >= MaximumRefinementLevel){
        for (i = 0; i < size; i++)
        {
            if (FlaggingField[i] > 0){
                FlaggingField[i] = 1;
                NumberOfFlaggedCells += 1;
                }
        }
        return NumberOfFlaggedCells;
    }
    float PopIIIMustRefineLifetime = 3; //Myr
    int istar = 0;
    // printf("[%d] Checking particles %d<->%d::%f<->%f\n", level, 2*nCellsPerRadius,                        
    //                 PopIIISupernovaMustRefineResolution, 
    //                 CellWidth[0][0]*LengthUnits/pc_cm, PopIIISupernovaRadius);
    int np = 0;
    while (istar < NumberOfParticles)
    {
        if (ParticleAttribute[0][istar] > 0.0 && ParticleType[istar] == 5) // if its a star
        {
            /* factor = birthtime + lifetime + refineTime */
            float factor = (ParticleAttribute[0][istar] + ParticleAttribute[1][istar]) * TimeUnits + PopIIIMustRefineLifetime * 3.1557e13;
            float m = ParticleMass[istar]*(DensityUnits*pow(LengthUnits*CellWidth[0][0], 3))/SolarMass;
            if ((((11 < m && m < 40) || (140 < m && m < 260)) || (m < 1e-10)) && (Time * TimeUnits < factor))
            {
                // fprintf(stdout, "Flagging cells for particle m=%"FSYM" fact=%"FSYM" Time=%"FSYM"\n", m, factor/3.1557e13, Time*TimeUnits/3.1557e13);
                // refine a radius a*Pop3 SN radius
                FLOAT radius = 1.1 * PopIIISupernovaRadius * 3.086e18 / LengthUnits;
                CellSize = FLOAT(CellWidth[0][0]);

                for (k = Start[2]; k <= End[2]; k++)
                    for (j = Start[1]; j <= End[1]; j++)
                        for (i = Start[0]; i <= End[0]; i++)
                        {
                            index = i + j * GridDimension[0] +
                                    k * GridDimension[1] * GridDimension[0];

                            xpos = GridLeftEdge[0] + (FLOAT(i - Start[0]) + 0.5) * CellSize;
                            ypos = GridLeftEdge[1] + (FLOAT(j - Start[1]) + 0.5) * CellSize;
                            zpos = GridLeftEdge[2] + (FLOAT(k - Start[2]) + 0.5) * CellSize;
                            float dist = pow(xpos - ParticlePosition[0][istar], 2) 
                                        + pow(ypos - ParticlePosition[1][istar], 2) 
                                        + pow(zpos - ParticlePosition[2][istar], 2);
                            if (dist < radius * radius)
                            {
                                FlaggingField[index] += 1;
                                np += 1;
                            }
                        } // end loop over cells
            }             // end if its a star
        }                 // end if refine star
        istar++;
    } // end for particles
      //                   // sum cells flagged
    for (i = 0; i < size; i++)
    {
        if (FlaggingField[i] > 0){
            FlaggingField[i] = 1;
            NumberOfFlaggedCells += 1;
            }
    }
    // fprintf(stdout, "[ %d ] P3 flagged %d cells\n", level, NumberOfFlaggedCells);

    return NumberOfFlaggedCells;
}