/***********************************************************************
/
/  COMMUNICATION ROUTINES: WRAPPED AND OVERLOADED MPI ROUTINES
/
/  written by: Greg Bryan
/  date:       December, 1997
/  modified1:  June, 2009 by John H. Wise (overloaded routines)
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

#ifdef MPI_INSTRUMENTATION
#define START_TIMING starttime = MPI_Wtime();
#define END_TIMING				\
  endtime = MPI_Wtime();			\
  timer[16]+= endtime-starttime;		\
  counter[16] ++;				\
  GlobalCommunication += endtime-starttime;	\
  CommunicationTime += endtime-starttime;
#else
#define START_TIMING ;
#define END_TIMING ;
#endif

 
/***********************************************************************
                             MINIMUM VALUE
************************************************************************/
 
Eflt32 CommunicationMinValue(Eflt32 Value)
{
 
  if (NumberOfProcessors == 1)
    return Value;
 
  Eflt32 ReturnValue = Value;
#ifdef USE_MPI
  START_TIMING;
  MPI_Datatype DataType = MPI_FLOAT;
  MPI_Arg Count = 1;
  MPI_Allreduce(&Value, &ReturnValue, Count, DataType, MPI_MIN, MPI_COMM_WORLD);
  END_TIMING;
 
#endif /* USE_MPI */
 
  return ReturnValue;
}

Eflt64 CommunicationMinValue(Eflt64 Value)
{
 
  if (NumberOfProcessors == 1)
    return Value;
 
  Eflt64 ReturnValue = Value;
#ifdef USE_MPI
   MPI_Datatype DataType = MPI_DOUBLE;
  MPI_Arg Count = 1;

  START_TIMING;
  MPI_Allreduce(&Value, &ReturnValue, Count, DataType, MPI_MIN, MPI_COMM_WORLD);
  END_TIMING;
 
#endif /* USE_MPI */
 
  return ReturnValue;
}

Eflt128 CommunicationMinValue(Eflt128 Value)
{
 
  if (NumberOfProcessors == 1)
    return Value;
 
  Eflt128 ReturnValue = Value;
#ifdef USE_MPI
  MPI_Datatype DataType = MPI_LONG_DOUBLE;
  MPI_Arg Count = 1;
  
  START_TIMING;
  MPI_Allreduce(&Value, &ReturnValue, Count, DataType, MPI_MIN, MPI_COMM_WORLD);
  END_TIMING;
 
#endif /* USE_MPI */
 
  return ReturnValue;
}

Eint32 CommunicationMinValue(Eint32 Value)
{
 
  if (NumberOfProcessors == 1)
    return Value;
 
  Eint32 ReturnValue = Value;
#ifdef USE_MPI
  MPI_Datatype DataType = MPI_INT;
  MPI_Arg Count = 1;

  START_TIMING;
  MPI_Allreduce(&Value, &ReturnValue, Count, DataType, MPI_MIN, MPI_COMM_WORLD);
  END_TIMING;
 
#endif /* USE_MPI */
 
  return ReturnValue;
}
 
Eint64 CommunicationMinValue(Eint64 Value)
{
 
  if (NumberOfProcessors == 1)
    return Value;
 
  Eint64 ReturnValue = Value;
#ifdef USE_MPI
  MPI_Datatype DataType = MPI_LONG_LONG_INT;
  MPI_Arg Count = 1;

  START_TIMING;
  MPI_Allreduce(&Value, &ReturnValue, Count, DataType, MPI_MIN, MPI_COMM_WORLD);
  END_TIMING;

#endif /* USE_MPI */
 
  return ReturnValue;
}

/***********************************************************************
                             MAXIMUM VALUE
************************************************************************/
 
Eflt32 CommunicationMaxValue(Eflt32 Value)
{
 
  if (NumberOfProcessors == 1)
    return Value;
 
 
  Eflt32 ReturnValue = Value;
#ifdef USE_MPI
  MPI_Datatype DataType = MPI_FLOAT;
  MPI_Arg Count = 1;

  START_TIMING;
  MPI_Allreduce(&Value, &ReturnValue, Count, DataType, MPI_MAX, MPI_COMM_WORLD);
  END_TIMING;
 
#endif /* USE_MPI */
 
  return ReturnValue;
}

Eflt64 CommunicationMaxValue(Eflt64 Value)
{
 
  if (NumberOfProcessors == 1)
    return Value;
 
  Eflt64 ReturnValue = Value;
#ifdef USE_MPI
  MPI_Datatype DataType = MPI_DOUBLE;
  MPI_Arg Count = 1;

  START_TIMING;
  MPI_Allreduce(&Value, &ReturnValue, Count, DataType, MPI_MAX, MPI_COMM_WORLD);
  END_TIMING;
 
#endif /* USE_MPI */
 
  return ReturnValue;
}

Eflt128 CommunicationMaxValue(Eflt128 Value)
{
 
  if (NumberOfProcessors == 1)
    return Value;
 
  Eflt128 ReturnValue = Value;
#ifdef USE_MPI
  MPI_Datatype DataType = MPI_LONG_DOUBLE;
  MPI_Arg Count = 1;

  START_TIMING;
  MPI_Allreduce(&Value, &ReturnValue, Count, DataType, MPI_MAX, MPI_COMM_WORLD);
  END_TIMING;
 
#endif /* USE_MPI */
 
  return ReturnValue;
}
Eint32 CommunicationMaxValue(Eint32 Value)
{
 
  if (NumberOfProcessors == 1)
    return Value;
 
  Eint32 ReturnValue = Value;
#ifdef USE_MPI
  MPI_Datatype DataType = MPI_INT;
  MPI_Arg Count = 1;

  START_TIMING;
  MPI_Allreduce(&Value, &ReturnValue, Count, DataType, MPI_MAX, MPI_COMM_WORLD);
  END_TIMING;
 
#endif /* USE_MPI */
 
  return ReturnValue;
}

Eint64 CommunicationMaxValue(Eint64 Value)
{
 
  if (NumberOfProcessors == 1)
    return Value;
 
 
  Eint64 ReturnValue = Value;
#ifdef USE_MPI
  MPI_Datatype DataType = MPI_LONG_LONG_INT;
  MPI_Arg Count = 1;

  START_TIMING;
  MPI_Allreduce(&Value, &ReturnValue, Count, DataType, MPI_MAX, MPI_COMM_WORLD);
  END_TIMING;
 
#endif /* USE_MPI */
 
  return ReturnValue;
}

/***********************************************************************
                              SUM VALUES
************************************************************************/ 
 
int CommunicationSumValues(Eflt32 *Values, int Number)
{
 
  if (NumberOfProcessors == 1)
    return SUCCESS;
 
#ifdef USE_MPI
 
  MPI_Datatype DataType = MPI_FLOAT;
  MPI_Arg Count = Number;
 
  int i;
  Eflt32 *buffer = new Eflt32[Number];
  for (i = 0; i < Number; i++)
    buffer[i] = Values[i];
 
  START_TIMING;
  MPI_Reduce(buffer, Values, Count, DataType, MPI_SUM, ROOT_PROCESSOR, MPI_COMM_WORLD);
  END_TIMING;
 
  delete [] buffer;
 
#endif /* USE_MPI */
 
  return SUCCESS;
}

int CommunicationSumValues(Eflt64 *Values, int Number)
{
 
  if (NumberOfProcessors == 1)
    return SUCCESS;
 
#ifdef USE_MPI
 
  MPI_Datatype DataType = MPI_DOUBLE;
  MPI_Arg Count = Number;
 
  int i;
  Eflt64 *buffer = new Eflt64[Number];
  for (i = 0; i < Number; i++)
    buffer[i] = Values[i];
 
  START_TIMING;
  MPI_Reduce(buffer, Values, Count, DataType, MPI_SUM, ROOT_PROCESSOR, MPI_COMM_WORLD);
  END_TIMING;
 
  delete [] buffer;
 
#endif /* USE_MPI */
 
  return SUCCESS;
}

int CommunicationSumValues(Eflt128 *Values, int Number)
{
 
  if (NumberOfProcessors == 1)
    return SUCCESS;
 
#ifdef USE_MPI
 
  MPI_Datatype DataType = MPI_LONG_DOUBLE;
  MPI_Arg Count = Number;
 
  int i;
  Eflt128 *buffer = new Eflt128[Number];
  for (i = 0; i < Number; i++)
    buffer[i] = Values[i];
 
  START_TIMING;
  MPI_Reduce(buffer, Values, Count, DataType, MPI_SUM, ROOT_PROCESSOR, MPI_COMM_WORLD);
  END_TIMING;
 
  delete [] buffer;
 
#endif /* USE_MPI */
 
  return SUCCESS;
}

int CommunicationSumValues(Eint32 *Values, int Number)
{
 
  if (NumberOfProcessors == 1)
    return SUCCESS;
 
#ifdef USE_MPI
 
  MPI_Datatype DataType = MPI_INT;
  MPI_Arg Count = Number;
 
  int i;
  Eint32 *buffer = new Eint32[Number];
  for (i = 0; i < Number; i++)
    buffer[i] = Values[i];
 
  START_TIMING;
  MPI_Reduce(buffer, Values, Count, DataType, MPI_SUM, ROOT_PROCESSOR, MPI_COMM_WORLD);
  END_TIMING;
 
  delete [] buffer;
 
#endif /* USE_MPI */
 
  return SUCCESS;
}

int CommunicationSumValues(Eint64 *Values, int Number)
{
 
  if (NumberOfProcessors == 1)
    return SUCCESS;
 
#ifdef USE_MPI
 
  MPI_Datatype DataType = MPI_LONG_LONG_INT;
  MPI_Arg Count = Number;
 
  int i;
  Eint64 *buffer = new Eint64[Number];
  for (i = 0; i < Number; i++)
    buffer[i] = Values[i];
 
  START_TIMING;
  MPI_Reduce(buffer, Values, Count, DataType, MPI_SUM, ROOT_PROCESSOR, MPI_COMM_WORLD);
  END_TIMING;
 
  delete [] buffer;
 
#endif /* USE_MPI */
 
  return SUCCESS;
}

 
/***********************************************************************
                           SUM ALL VALUES
************************************************************************/ 
 
int CommunicationAllSumValues(Eflt32 *Values, int Number)
{
 
  if (NumberOfProcessors == 1)
    return SUCCESS;
 
#ifdef USE_MPI
 
  MPI_Datatype DataType = MPI_FLOAT;
  MPI_Arg Count = Number;
 
  int i;
  Eflt32 *buffer = new Eflt32[Number];
  for (i = 0; i < Number; i++)
    buffer[i] = Values[i];
 
  START_TIMING;
  MPI_Allreduce(buffer, Values, Count, DataType, MPI_SUM, MPI_COMM_WORLD);
  END_TIMING;
 
  delete [] buffer;
 
#endif /* USE_MPI */
 
  return SUCCESS;
}

int CommunicationAllSumValues(Eflt64 *Values, int Number)
{
 
  if (NumberOfProcessors == 1)
    return SUCCESS;
 
#ifdef USE_MPI
 
  MPI_Datatype DataType = MPI_DOUBLE;
  MPI_Arg Count = Number;
 
  int i;
  Eflt64 *buffer = new Eflt64[Number];
  for (i = 0; i < Number; i++)
    buffer[i] = Values[i];
 
  START_TIMING;
  MPI_Allreduce(buffer, Values, Count, DataType, MPI_SUM, MPI_COMM_WORLD);
  END_TIMING;
 
  delete [] buffer;
 
#endif /* USE_MPI */
 
  return SUCCESS;
}

int CommunicationAllSumValues(Eflt128 *Values, int Number)
{
 
  if (NumberOfProcessors == 1)
    return SUCCESS;
 
#ifdef USE_MPI
 
  MPI_Datatype DataType = MPI_LONG_DOUBLE;
  MPI_Arg Count = Number;
 
  int i;
  Eflt128 *buffer = new Eflt128[Number];
  for (i = 0; i < Number; i++)
    buffer[i] = Values[i];
 
  START_TIMING;
  MPI_Allreduce(buffer, Values, Count, DataType, MPI_SUM, MPI_COMM_WORLD);
  END_TIMING;
 
  delete [] buffer;
 
#endif /* USE_MPI */
 
  return SUCCESS;
}

int CommunicationAllSumValues(Eint32 *Values, int Number)
{
 
  if (NumberOfProcessors == 1)
    return SUCCESS;
 
#ifdef USE_MPI
 
  MPI_Datatype DataType = MPI_INT;
  MPI_Arg Count = Number;
 
  int i;
  Eint32 *buffer = new Eint32[Number];
  for (i = 0; i < Number; i++)
    buffer[i] = Values[i];
 
  START_TIMING;
  MPI_Allreduce(buffer, Values, Count, DataType, MPI_SUM, MPI_COMM_WORLD);
  END_TIMING;
 
  delete [] buffer;
 
#endif /* USE_MPI */
 
  return SUCCESS;
}

int CommunicationAllSumValues(Eint64 *Values, int Number)
{
 
  if (NumberOfProcessors == 1)
    return SUCCESS;
 
#ifdef USE_MPI
 
  MPI_Datatype DataType = MPI_LONG_LONG_INT;
  MPI_Arg Count = Number;
 
  int i;
  Eint64 *buffer = new Eint64[Number];
  for (i = 0; i < Number; i++)
    buffer[i] = Values[i];
 
  START_TIMING;
  MPI_Allreduce(buffer, Values, Count, DataType, MPI_SUM, MPI_COMM_WORLD);
  END_TIMING;
 
  delete [] buffer;
 
#endif /* USE_MPI */
 
  return SUCCESS;
}

#ifdef USE_MPI
/***********************************************************************
                        GENERAL REDUCE CALL
************************************************************************/ 

int CommunicationReduceValues(Eflt32 *Values, int Number, 
			      MPI_Op ReduceOperation)
{
  
  if (NumberOfProcessors == 1)
    return SUCCESS;

  MPI_Datatype DataType = MPI_FLOAT;
  MPI_Arg Count = Number;
  
  int i;
  Eflt32 *buffer = new Eflt32[Number];
  for (i = 0; i < Number; i++)
    buffer[i] = Values[i];

  START_TIMING;
  MPI_Reduce(buffer, Values, Count, DataType, ReduceOperation, ROOT_PROCESSOR,
	     MPI_COMM_WORLD);
  END_TIMING;

  delete [] buffer;

  return SUCCESS;
}


int CommunicationReduceValues(Eflt64 *Values, int Number, 
			      MPI_Op ReduceOperation)
{
  
  if (NumberOfProcessors == 1)
    return SUCCESS;

  MPI_Datatype DataType = MPI_DOUBLE;
  MPI_Arg Count = Number;
  
  int i;
  Eflt64 *buffer = new Eflt64[Number];
  for (i = 0; i < Number; i++)
    buffer[i] = Values[i];

  START_TIMING;
  MPI_Reduce(buffer, Values, Count, DataType, ReduceOperation, ROOT_PROCESSOR,
	     MPI_COMM_WORLD);
  END_TIMING;

  delete [] buffer;

  return SUCCESS;
}

int CommunicationReduceValues(Eflt128 *Values, int Number, 
			      MPI_Op ReduceOperation)
{
  
  if (NumberOfProcessors == 1)
    return SUCCESS;

  MPI_Datatype DataType = MPI_LONG_DOUBLE;
  MPI_Arg Count = Number;
  
  int i;
  Eflt128 *buffer = new Eflt128[Number];
  for (i = 0; i < Number; i++)
    buffer[i] = Values[i];

  START_TIMING;
  MPI_Reduce(buffer, Values, Count, DataType, ReduceOperation, ROOT_PROCESSOR,
	     MPI_COMM_WORLD);
  END_TIMING;

  delete [] buffer;

  return SUCCESS;
}

int CommunicationReduceValues(Eint32 *Values, int Number, 
			      MPI_Op ReduceOperation)
{
  
  if (NumberOfProcessors == 1)
    return SUCCESS;

  MPI_Datatype DataType = MPI_INT;
  MPI_Arg Count = Number;
  
  int i;
  Eint32 *buffer = new Eint32[Number];
  for (i = 0; i < Number; i++)
    buffer[i] = Values[i];

  START_TIMING;
  MPI_Reduce(buffer, Values, Count, DataType, ReduceOperation, ROOT_PROCESSOR,
	     MPI_COMM_WORLD);
  END_TIMING;

  delete [] buffer;

  return SUCCESS;
}

int CommunicationReduceValues(Eint64 *Values, int Number, 
			      MPI_Op ReduceOperation)
{
  
  if (NumberOfProcessors == 1)
    return SUCCESS;

  MPI_Datatype DataType = MPI_LONG_LONG_INT;
  MPI_Arg Count = Number;
  
  int i;
  Eint64 *buffer = new Eint64[Number];
  for (i = 0; i < Number; i++)
    buffer[i] = Values[i];

  START_TIMING;
  MPI_Reduce(buffer, Values, Count, DataType, ReduceOperation, ROOT_PROCESSOR,
	     MPI_COMM_WORLD);
  END_TIMING;

  delete [] buffer;

  return SUCCESS;
}

/***********************************************************************
                        GENERAL ALLREDUCE CALL
************************************************************************/ 

int CommunicationAllReduceValues(Eflt32 *Values, int Number, 
				 MPI_Op ReduceOperation)
{
  
  if (NumberOfProcessors == 1)
    return SUCCESS;

  MPI_Datatype DataType = MPI_FLOAT;
  MPI_Arg Count = Number;
  
  int i;
  Eflt32 *buffer = new Eflt32[Number];
  for (i = 0; i < Number; i++)
    buffer[i] = Values[i];

  START_TIMING;
  MPI_Allreduce(buffer, Values, Count, DataType, ReduceOperation,
		MPI_COMM_WORLD);
  END_TIMING;

  delete [] buffer;

  return SUCCESS;
}

int CommunicationAllReduceValues(Eflt64 *Values, int Number, 
				 MPI_Op ReduceOperation)
{
  
  if (NumberOfProcessors == 1)
    return SUCCESS;

  MPI_Datatype DataType = MPI_DOUBLE;
  MPI_Arg Count = Number;
  
  int i;
  Eflt64 *buffer = new Eflt64[Number];
  for (i = 0; i < Number; i++)
    buffer[i] = Values[i];

  START_TIMING;
  MPI_Allreduce(buffer, Values, Count, DataType, ReduceOperation,
		MPI_COMM_WORLD);
  END_TIMING;

  delete [] buffer;

  return SUCCESS;
}

int CommunicationAllReduceValues(Eflt128 *Values, int Number, 
				 MPI_Op ReduceOperation)
{
  
  if (NumberOfProcessors == 1)
    return SUCCESS;

  MPI_Datatype DataType = MPI_LONG_DOUBLE;
  MPI_Arg Count = Number;
  
  int i;
  Eflt128 *buffer = new Eflt128[Number];
  for (i = 0; i < Number; i++)
    buffer[i] = Values[i];

  START_TIMING;
  MPI_Allreduce(buffer, Values, Count, DataType, ReduceOperation,
		MPI_COMM_WORLD);
  END_TIMING;

  delete [] buffer;

  return SUCCESS;
}

int CommunicationAllReduceValues(Eint32 *Values, int Number, 
				 MPI_Op ReduceOperation)
{
  
  if (NumberOfProcessors == 1)
    return SUCCESS;

  MPI_Datatype DataType = MPI_INT;
  MPI_Arg Count = Number;
  
  int i;
  Eint32 *buffer = new Eint32[Number];
  for (i = 0; i < Number; i++)
    buffer[i] = Values[i];

  START_TIMING;
  MPI_Allreduce(buffer, Values, Count, DataType, ReduceOperation,
		MPI_COMM_WORLD);
  END_TIMING;

  delete [] buffer;

  return SUCCESS;
}

int CommunicationAllReduceValues(Eint64 *Values, int Number, 
				 MPI_Op ReduceOperation)
{
  
  if (NumberOfProcessors == 1)
    return SUCCESS;

  MPI_Datatype DataType = MPI_LONG_LONG_INT;
  MPI_Arg Count = Number;
  
  int i;
  Eint64 *buffer = new Eint64[Number];
  for (i = 0; i < Number; i++)
    buffer[i] = Values[i];

  START_TIMING;
  MPI_Allreduce(buffer, Values, Count, DataType, ReduceOperation,
		MPI_COMM_WORLD);
  END_TIMING;

  delete [] buffer;

  return SUCCESS;
}

#endif /* USE_MPI */
/***********************************************************************
                               BARRIER
************************************************************************/ 

int CommunicationBarrier()
{
#ifdef USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif /* USE_MPI */
  return SUCCESS;
}

/************************************************************************
  Just like grid::CommunicationMethodShouldExit (in Grid.h) but for
  any processor numbers
************************************************************************/

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

    return FAIL; /* i.e. method should not exit immediately. */

}
