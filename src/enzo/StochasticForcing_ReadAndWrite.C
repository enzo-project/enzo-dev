/***********************************************************************
 *
 *  STOCHASTIC FORCING CLASS: ReadSpectrum, WriteSpectrum, 
 *                            WriteParameters
 *
 *  written by: Wolfram Schmidt
 *  date:       November, 2005
 *  modified1: Oct, 2014: updated to support Enzo 2.4 // P. Grete
 *
 *  PURPOSE: reads/writes the forcing spectrum from/to a file
 *
 ***********************************************************************/

#include "preincludes.h"
#include "StochasticForcing.h"
#include "global_data.h"

int StochasticForcing::ReadSpectrum(char *fname)
{
    if (MyProcessorNumber == ROOT_PROCESSOR) {

    FILE *fptr;

	if (debug) printf("ReadSpectrum: reading data\n");

    if ((fptr = fopen(fname, "r")) == NULL) {
        ENZO_VFAIL("ReadSpectrum: failed to open file  %s\n", fname)
    }


    char line[MAX_LINE_LENGTH];
    
    fgets(line, MAX_LINE_LENGTH, fptr);
    sscanf(line,"%"ISYM" %"ISYM" %"ISYM,&seed, &idum2, &iy);
    
	for (int i = 0; i < 32; i++) {
        fgets(line, MAX_LINE_LENGTH, fptr);
        sscanf(line,"%"ISYM,&iv[i]);
    }

	for (int dim = 0; dim < SpectralRank; dim++)
	    for (int m = 0; m < NumNonZeroModes; m++) {
            fgets(line, MAX_LINE_LENGTH, fptr);
            sscanf(line,"%"FSYM" %"FSYM,&SpectrumEven[dim][m],&SpectrumOdd[dim][m]);
        }
    
    fclose(fptr);
    }
    

    return SUCCESS;
}

int StochasticForcing::WriteSpectrum(char *fname)
{
    if (MyProcessorNumber == ROOT_PROCESSOR) {

    FILE *fptr;	

	if (debug) printf("WriteSpectrum: writing data to file %s\n",fname);

    if ((fptr = fopen(fname, "w")) == NULL) {
        ENZO_VFAIL("WriteSpectrum: failed to open file  %s\n", fname)
    }

    fprintf(fptr,"%"ISYM" %"ISYM" %"ISYM"\n",seed,idum2,iy);

	for (int i = 0; i < 32; i++)
        fprintf(fptr,"%"ISYM"\n",iv[i]);

	for (int dim = 0; dim < SpectralRank; dim++)
	    for (int m = 0; m < NumNonZeroModes; m++)
            fprintf(fptr,"%.16"FSYM" %.16"FSYM"\n",SpectrumEven[dim][m],SpectrumOdd[dim][m]);
    
    
    fclose(fptr);
    }

    return SUCCESS;
}

void StochasticForcing::WriteParameters(FILE *fptr)
{
    fprintf(fptr, "DrivenFlowWeight            = %"FSYM"\n", SolenoidalWeight);
    fprintf(fptr, "DrivenFlowAlpha             = %"ISYM" %"ISYM" %"ISYM"\n", alpha[0], alpha[1], alpha[2]);
    fprintf(fptr, "DrivenFlowSeed             = %"ISYM"\n", seed);
    fprintf(fptr, "DrivenFlowBandWidth         = %"FSYM" %"FSYM" %"FSYM"\n", 
	    BandWidth[0], BandWidth[1], BandWidth[2]);
    fprintf(fptr, "DrivenFlowVelocity          = %"FSYM" %"FSYM" %"FSYM"\n", 
	    IntgrVelocity[0], IntgrVelocity[1], IntgrVelocity[2]);
    fprintf(fptr, "DrivenFlowAutoCorrl         = %"FSYM" %"FSYM" %"FSYM"\n", 
	    AutoCorrlTime[0]/IntgrTime[0], AutoCorrlTime[1]/IntgrTime[1], AutoCorrlTime[2]/IntgrTime[2]);
}

