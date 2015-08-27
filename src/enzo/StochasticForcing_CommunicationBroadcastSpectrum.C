/***********************************************************************
/
/  STOCHASTIC FORCING CLASS METHOD: CommunicationBroadcastSpectrum
/
/  written by: Wolfram Schmidt
/  date:       August, 2005
/  modified1:  January, 2009
/  modified2: Oct, 2014: modified to support Enzo 2.4 // P. Grete
/
/  PURPOSE: broadcasts the forcing spectrum from root to others
/
************************************************************************/



#include "preincludes.h"
#include <stdlib.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "StochasticForcing.h"
#include "global_data.h"
#include "error.h"

/* function prototypes */

void my_exit(int status);

void StochasticForcing::CommunicationBroadcastSpectrum(void)
{
  if (NumberOfProcessors == 1) return;
    

#ifdef USE_MPI

  if (debug) printf("Broadcasting spectrum, proc # %"ISYM" \n",MyProcessorNumber);

  int dim, m, n;
  int TransferSize = 2*SpectralRank*NumNonZeroModes;
  float* buffer = NULL;

  buffer = new float[TransferSize];

  n = 0;
  for (dim = 0; dim < SpectralRank; dim++) {
      for (m = 0; m < NumNonZeroModes; m++) buffer[n++] = SpectrumEven[dim][m];
      for (m = 0; m < NumNonZeroModes; m++) buffer[n++] = SpectrumOdd [dim][m];      
  }

#ifdef MPI_INSTRUMENTATION
  starttime = MPI_Wtime();
#endif

  MPI_Datatype DataTypeFloat = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  MPI_Arg stat;

  stat = MPI_Bcast(buffer, TransferSize, DataTypeFloat, ROOT_PROCESSOR, MPI_COMM_WORLD);
    if( stat != MPI_SUCCESS ){my_exit(EXIT_FAILURE);}

#ifdef MPI_INSTRUMENTATION
  endtime = MPI_Wtime();
  timer[15]+= endtime-starttime;
  counter[15] ++;
  GlobalCommunication += endtime-starttime;
  CommunicationTime += endtime-starttime;
#endif /* MPI_INSTRUMENTATION */

  n = 0;
  for (dim = 0; dim < SpectralRank; dim++) {
      for (m = 0; m < NumNonZeroModes; m++) SpectrumEven[dim][m] = buffer[n++];
      for (m = 0; m < NumNonZeroModes; m++) SpectrumOdd [dim][m] = buffer[n++];
  }

  delete [] buffer;

#endif /* USE_MPI */
}
