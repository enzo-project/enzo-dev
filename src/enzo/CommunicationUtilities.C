/***********************************************************************
/
/  COMMUNICATION ROUTINE: FIND MINIMUM VALUE AMOUNG PROCESSORS
/
/  written by: Greg Bryan
/  date:       December, 1997
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
 
#include <stdio.h>
#include <string.h>
 
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
 
 
 
 
float CommunicationMinValue(float Value)
{
 
  if (NumberOfProcessors == 1)
    return Value;
 
  float ReturnValue = Value;
 
#ifdef USE_MPI
 
  MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  MPI_Arg Count;
 
//  printf("min: %"ISYM" sending %"FSYM"\n", MyProcessorNumber, Value);
 
#ifdef MPI_INSTRUMENTATION
  starttime = MPI_Wtime();
#endif

  Count = 1;
 
  MPI_Allreduce(&Value, &ReturnValue, Count, DataType, MPI_MIN, MPI_COMM_WORLD);
 
#ifdef MPI_INSTRUMENTATION
  endtime = MPI_Wtime();
  timer[16]+= endtime-starttime;
  counter[16] ++;
  GlobalCommunication += endtime-starttime;
  CommunicationTime += endtime-starttime;
#endif /* MPI_INSTRUMENTATION */
 
#endif /* USE_MPI */
 
  return ReturnValue;
}
 
 
float CommunicationMaxValue(float Value)
{
 
  if (NumberOfProcessors == 1)
    return Value;
 
  float ReturnValue = Value;
 
#ifdef USE_MPI
 
  MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  MPI_Arg Count;
 
//  printf("min: %"ISYM" sending %"FSYM"\n", MyProcessorNumber, Value);
 
#ifdef MPI_INSTRUMENTATION
  starttime = MPI_Wtime();
#endif

  Count = 1;
 
  MPI_Allreduce(&Value, &ReturnValue, Count, DataType, MPI_MAX, MPI_COMM_WORLD);
 
#ifdef MPI_INSTRUMENTATION
  endtime = MPI_Wtime();
  timer[16]+= endtime-starttime;
  counter[16] ++;
  GlobalCommunication += endtime-starttime;
  CommunicationTime += endtime-starttime;
#endif /* MPI_INSTRUMENTATION */
 
#endif /* USE_MPI */
 
  return ReturnValue;
}
 
 
int CommunicationSumValues(float *Values, int Number)
{
 
  if (NumberOfProcessors == 1)
    return SUCCESS;
 
#ifdef USE_MPI
 
  MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  MPI_Arg Count;
 
  int i;
  float *buffer = new float[Number];
  for (i = 0; i < Number; i++)
    buffer[i] = Values[i];
 
  Count = Number;

  MPI_Reduce(buffer, Values, Count, DataType, MPI_SUM, ROOT_PROCESSOR, MPI_COMM_WORLD);
 
  delete [] buffer;
 
#endif /* USE_MPI */
 
  return SUCCESS;
}
 
 
 
 
int CommunicationAllSumValues(float *Values, int Number)
{
 
  if (NumberOfProcessors == 1)
    return SUCCESS;
 
#ifdef USE_MPI
 
  MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  MPI_Arg Count;
 
  int i;
  float *buffer = new float[Number];
  for (i = 0; i < Number; i++)
    buffer[i] = Values[i];
 
  Count = Number;
 
  MPI_Allreduce(buffer, Values, Count, DataType, MPI_SUM, MPI_COMM_WORLD);
 
  delete [] buffer;
 
#endif /* USE_MPI */
 
  return SUCCESS;
}
 
 
 
 
int CommunicationAllSumIntegerValues(int *Values, int Number)
{
  if (NumberOfProcessors == 1)
    return SUCCESS;
 
#ifdef USE_MPI

  MPI_Datatype DataTypeInt = (sizeof(int) == 4) ? MPI_INT : MPI_LONG_LONG_INT; 
  MPI_Arg Count;
 
  int i;
  int *buffer = new int[Number];
  for (i = 0; i < Number; i++)
    buffer[i] = Values[i];
 
  Count = Number;
 
  MPI_Allreduce(buffer, Values, Count, DataTypeInt, MPI_SUM, MPI_COMM_WORLD);
 
  delete [] buffer;
 
#endif /* USE_MPI */
 
  return SUCCESS;
}

#ifdef USE_MPI
  
int CommunicationReduceValues(float *Values, int Number, 
			      MPI_Op ReduceOperation)
{
  
  if (NumberOfProcessors == 1)
    return SUCCESS;

  MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  
  int i;
  float *buffer = new float[Number];
  for (i = 0; i < Number; i++)
    buffer[i] = Values[i];

  MPI_Reduce(buffer, Values, Number, DataType, ReduceOperation, ROOT_PROCESSOR,
	     MPI_COMM_WORLD);

  delete [] buffer;

  return SUCCESS;
}

int CommunicationAllReduceValues(float *Values, int Number, 
				 MPI_Op ReduceOperation)
{
  
  if (NumberOfProcessors == 1)
    return SUCCESS;

  MPI_Datatype DataType = (sizeof(float) == 4) ? MPI_FLOAT : MPI_DOUBLE;
  
  int i;
  float *buffer = new float[Number];
  for (i = 0; i < Number; i++)
    buffer[i] = Values[i];

  MPI_Allreduce(buffer, Values, Number, DataType, ReduceOperation,
		MPI_COMM_WORLD);

  delete [] buffer;

  return SUCCESS;
}


int CommunicationAllReduceValuesFLOAT(FLOAT *Values, int Number, 
				      MPI_Op ReduceOperation)
{
  
  if (NumberOfProcessors == 1)
    return SUCCESS;

  MPI_Datatype DataType;
  switch(sizeof(FLOAT)) {
    case 4:
        DataType = MPI_FLOAT;
        break;
    case 8:
        DataType = MPI_DOUBLE;
        break;
    case 16:
        DataType = MPI_LONG_DOUBLE;
        break;
  }
   
  
  int i;
  FLOAT *buffer = new FLOAT[Number];
  for (i = 0; i < Number; i++)
    buffer[i] = Values[i];

  MPI_Allreduce(buffer, Values, Number, DataType, ReduceOperation, MPI_COMM_WORLD);

  delete [] buffer;

  return SUCCESS;
}

int CommunicationReduceValuesFLOAT(FLOAT *Values, int Number, 
			      MPI_Op ReduceOperation)
{
  
  if (NumberOfProcessors == 1)
    return SUCCESS;

  MPI_Datatype DataType;
  switch(sizeof(FLOAT)) {
    case 4:
        DataType = MPI_FLOAT;
        break;
    case 8:
        DataType = MPI_DOUBLE;
        break;
    case 16:
        DataType = MPI_LONG_DOUBLE;
        break;
  }
   
  
  int i;
  FLOAT *buffer = new FLOAT[Number];
  for (i = 0; i < Number; i++)
    buffer[i] = Values[i];

  MPI_Reduce(buffer, Values, Number, DataType, ReduceOperation, ROOT_PROCESSOR,
	     MPI_COMM_WORLD);

  delete [] buffer;

  return SUCCESS;
}

int CommunicationReduceValuesDouble(double *Values, int Number, 
				    MPI_Op ReduceOperation)
{
  
  if (NumberOfProcessors == 1)
    return SUCCESS;

  MPI_Datatype DataType;
  DataType = MPI_DOUBLE;
  
  int i;
  double *buffer = new double[Number];
  for (i = 0; i < Number; i++)
    buffer[i] = Values[i];

  MPI_Reduce(buffer, Values, Number, DataType, ReduceOperation, ROOT_PROCESSOR,
	     MPI_COMM_WORLD);

  delete [] buffer;

  return SUCCESS;
}

int CommunicationAllReduceValuesINT(int *Values, int Number, 
				    MPI_Op ReduceOperation)
{
  
  if (NumberOfProcessors == 1)
    return SUCCESS;

  int i;
  int *buffer = new int[Number];
  for (i = 0; i < Number; i++)
    buffer[i] = Values[i];

  MPI_Allreduce(buffer, Values, Number, MPI_INT, ReduceOperation,
		MPI_COMM_WORLD);

  delete [] buffer;

  return SUCCESS;
}

#endif /* USE_MPI */
  
int CommunicationSumValuesFLOAT(FLOAT *Values, int Number)
{
#ifdef USE_MPI
  return CommunicationReduceValuesFLOAT(Values, Number, MPI_SUM);
#else /* USE_MPI */
  return SUCCESS;
#endif /* USE_MPI */
}

  
int CommunicationAllSumValuesFLOAT(FLOAT *Values, int Number)
{
  
  if (NumberOfProcessors == 1)
    return SUCCESS;

#ifdef USE_MPI

  int i;
  FLOAT *buffer = new FLOAT[Number];
  for (i = 0; i < Number; i++)
    buffer[i] = Values[i];

  MPI_Datatype DataType;
  switch(sizeof(FLOAT)) {
    case 4:
        DataType = MPI_FLOAT;
        break;
    case 8:
        DataType = MPI_DOUBLE;
        break;
    case 16:
        DataType = MPI_LONG_DOUBLE;
        break;
  }
 
  MPI_Allreduce(buffer, Values, Number, DataType, MPI_SUM, MPI_COMM_WORLD);

  delete [] buffer;

#endif /* USE_MPI */

  return SUCCESS;
}

int CommunicationBarrier()
{
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif /* USE_MPI */
  return SUCCESS;
}

/* Just like grid::CommunicationMethodShouldExit (in Grid.h) but for
   any processor numbers */

int CommunicationShouldExit(int FromProc, int ToProc)
{

    /* Return if neither grid lives on this processor. */

    if (MyProcessorNumber != ToProc && 
        MyProcessorNumber != FromProc)
      return SUCCESS;

    /* If the two grids are on the same processor then return if
       either in post-receive or receive modes to avoid duplicating method
       (i.e. action is only carried out if in send mode (or send-receive)). */

    if (ToProc == FromProc)
      if (CommunicationDirection == COMMUNICATION_POST_RECEIVE ||
	  CommunicationDirection == COMMUNICATION_RECEIVE)
        return SUCCESS;

    /* If in send mode then exit if the send grid is not on this processor. */

    if (CommunicationDirection == COMMUNICATION_SEND &&
        MyProcessorNumber != FromProc)
      return SUCCESS;

    /* If in either receive phase then exit if receive grid is not on this 
       processor. */

    if ((CommunicationDirection == COMMUNICATION_RECEIVE ||
         CommunicationDirection == COMMUNICATION_POST_RECEIVE) &&
        MyProcessorNumber != ToProc)
      return SUCCESS;

    ENZO_FAIL("Error in: "__FILE__); /* i.e. method should not exit immediately. */

}
