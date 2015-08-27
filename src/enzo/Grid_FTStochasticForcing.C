/***********************************************************************
/
/  GRID CLASS: FTStochasticForcing
/
/  written by: Wolfram Schmidt
/  date:       October 2005
/  modified1: Jul, 2013: modified header includes for enzo 2.3 // PG
/
/  PURPOSE: computes physical force field via inverse FT of the forcing
/           spectrum onto a particular grid including ghost cells;
/           the algorithm is based on iterative phase factor
/           multiplication and works well for coarse forcing spectra
/  
/  EXPLANATION: We are not using standard FFTW here due to its overhead
/   when applied to only few modes. Since the coefficients of the inverse
/   FT are given by exp(i*(k1*x + k2*y + k3*z)) = exp(i*k1*x) * exp(i*k2*y)*...
/   we pre-calculate the local coefficients (see Grid_Phases) and iteratively
/   update the coefficients with exp(i*k1*delta_x), etc. during the 
/   inverse FT. 
/   As warned below: this is only efficient for a limited number of modes.
/
************************************************************************/


#include "preincludes.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"


/* normalisation factor for FT */
#define FT_NORM 1.4142135623730951

int FindField(int f, int farray[], int n);

int grid::FTStochasticForcing(int FieldDim)
{
    int dim, m, mi, mj, mk;
    int size = Forcing.get_NumNonZeroModes();
    int numberOfGridZones = 1;
    int numberOfPhases = 0;
    /* WARNING: for broad spectra with a large number of modes, 
       the present implementation of the FT will not work efficiently */

    if (size > MAX_FORCING_MODES) {
    if (MyProcessorNumber == ROOT_PROCESSOR) 
        printf("Number of forcing modes exceeds MAX_FORCING_MODES = %"ISYM"\n",MAX_FORCING_MODES);
    return FAIL;
    }

    int StochAccelNum;  
    if ((StochAccelNum = FindField(DrivingField1, FieldType, NumberOfBaryonFields)) < 0) {
        fprintf(stderr, "Cannot find stochastic acceleration field.\n");
        return FAIL;
    };

    StochAccelNum += FieldDim;

    if (MyProcessorNumber == ProcessorNumber) {

    float sum;
    float buf[size];
    float ModeEven[size];
    float ModeOdd[size];
    float* PhaseFctEven[GridRank];
    float* PhaseFctOdd[GridRank];

    /* copy modes from Forcing object */

    Forcing.copy_SpectrumOdd (FieldDim, ModeOdd);
    Forcing.copy_SpectrumEven(FieldDim, ModeEven);

    for (dim = 0; dim < GridRank; dim++) {
        numberOfGridZones *= GridDimension[dim];
        numberOfPhases += size*GridDimension[dim];
    }

        if (debug)
            printf("Grid patch: #zones = %"ISYM", #phases = %"ISYM"\n",
                numberOfGridZones,numberOfPhases);

    /* check if memory is allocated */

    if ((GridRank > 1) && (numberOfPhases < numberOfGridZones)) {

        for (dim = 0; dim < GridRank; dim++) {

        PhaseFctEven[dim] = new float[numberOfPhases];
        PhaseFctOdd[dim]  = new float[numberOfPhases];

        // initialize phase factors
        if (dim == GridRank-1) {
            for (m = 0; m < size; m++) {
            PhaseFctEven[dim][m] = PhaseFctInitEven[m]; 
            PhaseFctOdd [dim][m] = PhaseFctInitOdd [m];
            }
        } else {
            for (m = 0; m < size; m++) {
            PhaseFctEven[dim][m] = 1.0; 
            PhaseFctOdd [dim][m] = 0.0;
            }
        }

        // iterate phase factors
        for (int i = 1; i < GridDimension[dim]; i++) {
            for (m = 0, mi = i*size; m < size; m++, mi++) {
            buf[m] = PhaseFctEven[dim][mi-size];
            PhaseFctEven[dim][mi] = PhaseFctMultEven[dim][m] * PhaseFctEven[dim][mi-size] -
                                    PhaseFctMultOdd [dim][m] * PhaseFctOdd [dim][mi-size]; 
            PhaseFctOdd [dim][mi] = PhaseFctMultEven[dim][m] * PhaseFctOdd [dim][mi-size] +
                                    PhaseFctMultOdd [dim][m] * buf[m];          
            }           
        }
        }

        switch (GridRank) {

        /* 2D inverse FT */
            
        case 2: 
            
            // sum over all modes
            for (int j = 0; j < GridDimension[1]; j++)
            for (int i = 0; i < GridDimension[0]; i++) {
                
                for (m = 0, mi = i*size, mj = j*size, sum = 0.0; m < size; m++, mi++, mj++) {
                sum += (PhaseFctEven[0][mi] * PhaseFctEven[1][mj] - 
                    PhaseFctOdd [0][mi] * PhaseFctOdd [1][mj]) * ModeEven[m] + 
                       (PhaseFctEven[0][mi] * PhaseFctOdd [1][mj]  + 
                        PhaseFctOdd [0][mi] * PhaseFctEven[1][mj]) * ModeOdd[m];
                }

                BaryonField[StochAccelNum][i + j*GridDimension[0]] = FT_NORM * sum;
                }

            break;
            
        /* 3D inverse FT */
            
        case 3: 
            
            // sum over all modes
            for (int k = 0; k < GridDimension[2]; k++)
            for (int j = 0; j < GridDimension[1]; j++)
                for (int i = 0; i < GridDimension[0]; i++) {
                
                for (m = 0, mi = i*size, mj = j*size, mk = k*size, sum = 0.0; m < size; m++, mi++, mj++, mk++) {
                    sum += ((PhaseFctEven[0][mi] * PhaseFctEven[1][mj] - 
                         PhaseFctOdd [0][mi] * PhaseFctOdd [1][mj]) * PhaseFctEven[2][mk] -
                        (PhaseFctEven[0][mi] * PhaseFctOdd [1][mj]  + 
                         PhaseFctOdd [0][mi] * PhaseFctEven[1][mj]) * PhaseFctOdd [2][mk]) * ModeEven[m] -
                       ((PhaseFctEven[0][mi] * PhaseFctEven[1][mj] - 
                         PhaseFctOdd [0][mi] * PhaseFctOdd [1][mj]) * PhaseFctOdd [2][mk] +
                        (PhaseFctEven[0][mi] * PhaseFctOdd [1][mj] + 
                         PhaseFctOdd [0][mi] * PhaseFctEven[1][mj]) * PhaseFctEven[2][mk]) * ModeOdd[m]; 
                }

                BaryonField[StochAccelNum][i + 
                                 j*GridDimension[0] + 
                                 k*GridDimension[0]*GridDimension[1]] = FT_NORM * sum;
                }

            break;

        default:

            return FAIL;
        }           

    } else {

        for (dim = 0; dim < GridRank; dim++) {
        PhaseFctEven[dim] = new float[size];
        PhaseFctOdd[dim]  = new float[size];
        }

        switch (GridRank) {

            /* 1D inverse FT */

        case 1: 
        
            // initialize phase factors
            for (m = 0; m < size; m++) {
            PhaseFctEven[0][m] = PhaseFctInitEven[m]; 
            PhaseFctOdd [0][m] = PhaseFctInitOdd [m];
            }
            
            for (int i = 0; i < GridDimension[0]; i++) {

            // sum over all modes
            for (m = 0, sum = 0.0; m < size; m++) {
                sum += PhaseFctEven[0][m]* ModeEven[m] + PhaseFctOdd[0][m] * ModeOdd[m];
            }

            BaryonField[StochAccelNum][i] = FT_NORM * sum;
            
            // iterate phase factors
            for (m = 0; m < size; m++) {
                buf[m] = PhaseFctEven[0][m];
                PhaseFctEven[0][m] = PhaseFctMultEven[0][m] * PhaseFctEven[0][m] -
                                 PhaseFctMultOdd [0][m] * PhaseFctOdd [0][m]; 
                PhaseFctOdd [0][m] = PhaseFctMultEven[0][m] * PhaseFctOdd [0][m] +
                                 PhaseFctMultOdd [0][m] * buf[m]; 
            }           
            }
            
            break;
            
        /* 2D inverse FT */
            
        case 2: 
            
            // initialize y-phase factors (complete phase information of grid origin)
            for (m = 0; m < size; m++) {
            PhaseFctEven[1][m] = PhaseFctInitEven[m];
            PhaseFctOdd [1][m] = PhaseFctInitOdd [m];
            }
            
            for (int j = 0; j < GridDimension[1]; j++) {
            
            // initialize x-phase factors
            for (m = 0; m < size; m++) {
                PhaseFctEven[0][m] = 1.0; 
                PhaseFctOdd [0][m] = 0.0;
            }
            
            for (int i = 0; i < GridDimension[0]; i++) {
                
                // sum over all modes
                for (m = 0, sum = 0.0; m < size; m++) {
                sum += (PhaseFctEven[0][m] * PhaseFctEven[1][m] - 
                    PhaseFctOdd [0][m] * PhaseFctOdd [1][m]) * ModeEven[m] + 
                       (PhaseFctEven[0][m] * PhaseFctOdd [1][m]  + 
                        PhaseFctOdd [0][m] * PhaseFctEven[1][m]) * ModeOdd[m];
                }
                
                BaryonField[StochAccelNum][i + j*GridDimension[0]] = FT_NORM * sum;
                
                // iterate x-phase factors
                for (m = 0; m < size; m++) {
                buf[m] = PhaseFctEven[0][m];
                PhaseFctEven[0][m] = PhaseFctMultEven[0][m] * PhaseFctEven[0][m] -
                                     PhaseFctMultOdd [0][m] * PhaseFctOdd [0][m]; 
                PhaseFctOdd [0][m] = PhaseFctMultEven[0][m] * PhaseFctOdd [0][m] +
                                     PhaseFctMultOdd [0][m] * buf[m]; 
                }           
            }
            
            // iterate y-phase factors
            for (m = 0; m < size; m++) {
                buf[m] = PhaseFctEven[1][m];
                PhaseFctEven[1][m] = PhaseFctMultEven[1][m] * PhaseFctEven[1][m] -
                                 PhaseFctMultOdd [1][m] * PhaseFctOdd [1][m]; 
                PhaseFctOdd [1][m] = PhaseFctMultEven[1][m] * PhaseFctOdd [1][m] +
                                 PhaseFctMultOdd [1][m] * buf[m]; 
            }   
            }       
            
            break;
            
        /* 3D inverse FT */
            
        case 3: 
            
            // initialize z-phase factors (complete phase information of grid origin)
            for (m = 0; m < size; m++) {
            PhaseFctEven[2][m] = PhaseFctInitEven[m]; 
            PhaseFctOdd [2][m] = PhaseFctInitOdd [m];
            }

            for (int k = 0; k < GridDimension[2]; k++) {
            
            // initialize y-phase factors
            for (m = 0; m < size; m++) {
                PhaseFctEven[1][m] = 1.0;
                PhaseFctOdd [1][m] = 0.0;
            }
            
            for (int j = 0; j < GridDimension[1]; j++) {
                
                // initialize x-phase factors
                for (m = 0; m < size; m++) {
                PhaseFctEven[0][m] = 1.0; 
                PhaseFctOdd [0][m] = 0.0;
                }
                
                for (int i = 0; i < GridDimension[0]; i++) {
                
                // sum over all modes
                for (m = 0, sum = 0.0; m < size; m++) {
                    sum += ((PhaseFctEven[0][m] * PhaseFctEven[1][m] - 
                         PhaseFctOdd [0][m] * PhaseFctOdd [1][m]) * PhaseFctEven[2][m] -
                        (PhaseFctEven[0][m] * PhaseFctOdd [1][m]  + 
                         PhaseFctOdd [0][m] * PhaseFctEven[1][m]) * PhaseFctOdd [2][m]) * ModeEven[m] -
                       ((PhaseFctEven[0][m] * PhaseFctEven[1][m] - 
                         PhaseFctOdd [0][m] * PhaseFctOdd [1][m]) * PhaseFctOdd [2][m] +
                        (PhaseFctEven[0][m] * PhaseFctOdd [1][m] + 
                         PhaseFctOdd [0][m] * PhaseFctEven[1][m]) * PhaseFctEven[2][m]) * ModeOdd[m]; 
                }
                
                BaryonField[StochAccelNum][i + 
                                 j*GridDimension[0] + 
                                 k*GridDimension[0]*GridDimension[1]] = FT_NORM * sum;
                // iterate x-phase factors
                for (m = 0; m < size; m++) {
                    buf[m] = PhaseFctEven[0][m];
                    PhaseFctEven[0][m] = PhaseFctMultEven[0][m] * PhaseFctEven[0][m] -
                                     PhaseFctMultOdd [0][m] * PhaseFctOdd [0][m]; 
                    PhaseFctOdd [0][m] = PhaseFctMultEven[0][m] * PhaseFctOdd [0][m] +
                                     PhaseFctMultOdd [0][m] * buf[m]; 
                }           
                }
                
                // iterate y-phase factors
                for (m = 0; m < size; m++) {
                buf[m] = PhaseFctEven[1][m];
                PhaseFctEven[1][m] = PhaseFctMultEven[1][m] * PhaseFctEven[1][m] -
                                     PhaseFctMultOdd [1][m] * PhaseFctOdd [1][m]; 
                PhaseFctOdd [1][m] = PhaseFctMultEven[1][m] * PhaseFctOdd [1][m] +
                                     PhaseFctMultOdd [1][m] * buf[m]; 
                
                }   
            }
            
            // iterate z-phase factors
            for (m = 0; m < size; m++) {
                buf[m] = PhaseFctEven[2][m];
                PhaseFctEven[2][m] = PhaseFctMultEven[2][m] * PhaseFctEven[2][m] -
                                 PhaseFctMultOdd [2][m] * PhaseFctOdd [2][m]; 
                PhaseFctOdd [2][m] = PhaseFctMultEven[2][m] * PhaseFctOdd [2][m] +
                                 PhaseFctMultOdd [2][m] * buf[m]; 
            }   
            }
            
            break;
            
        default:

            return FAIL;
        }
    }    
    
    for (dim = 0; dim < GridRank; dim++) {
        delete [] PhaseFctEven[dim];
        delete [] PhaseFctOdd[dim];
    }
    }

    return SUCCESS;
}
