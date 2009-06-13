#ifndef COMMUNICATION_UTILITIES_DEFINED__
#define COMMUNICATION_UTILITIES_DEFINED__

#ifdef USE_MPI
#include "mpi.h"
#endif /* USE_MPI */
#include "macros_and_parameters.h"
#include "typedefs.h"

Eflt32 CommunicationMinValue(Eflt32 Value);
Eflt64 CommunicationMinValue(Eflt64 Value);
Eflt128 CommunicationMinValue(Eflt128 Value);
Eint32 CommunicationMinValue(Eint32 Value);
Eint64 CommunicationMinValue(Eint64 Value);

Eflt32 CommunicationMaxValue(Eflt32 Value);
Eflt64 CommunicationMaxValue(Eflt64 Value);
Eflt128 CommunicationMaxValue(Eflt128 Value);
Eint32 CommunicationMaxValue(Eint32 Value);
Eint64 CommunicationMaxValue(Eint64 Value);

Eflt32 CommunicationSumValues(Eflt32 *Values, int Number);
Eflt64 CommunicationSumValues(Eflt64 *Values, int Number);
Eflt128 CommunicationSumValues(Eflt128 *Values, int Number);
Eint32 CommunicationSumValues(Eint32 *Values, int Number);
Eint64 CommunicationSumValues(Eint64 *Values, int Number);

Eflt32 CommunicationAllSumValues(Eflt32 *Values, int Number);
Eflt64 CommunicationAllSumValues(Eflt64 *Values, int Number);
Eflt128 CommunicationAllSumValues(Eflt128 *Values, int Number);
Eint32 CommunicationAllSumValues(Eint32 *Values, int Number);
Eint64 CommunicationAllSumValues(Eint64 *Values, int Number);

#ifdef USE_MPI
int CommunicationReduceValues(Eflt32 *Values, int Number, 
			      MPI_Op ReduceOperation);
int CommunicationReduceValues(Eflt64 *Values, int Number, 
			      MPI_Op ReduceOperation);
int CommunicationReduceValues(Eflt128 *Values, int Number, 
			      MPI_Op ReduceOperation);
int CommunicationReduceValues(Eint32 *Values, int Number, 
			      MPI_Op ReduceOperation);
int CommunicationReduceValues(Eint64 *Values, int Number, 
			      MPI_Op ReduceOperation);

int CommunicationAllReduceValues(Eflt32 *Values, int Number, 
				 MPI_Op ReduceOperation);
int CommunicationAllReduceValues(Eflt64 *Values, int Number, 
				 MPI_Op ReduceOperation);
int CommunicationAllReduceValues(Eflt128 *Values, int Number, 
				 MPI_Op ReduceOperation);
int CommunicationAllReduceValues(Eint32 *Values, int Number, 
				 MPI_Op ReduceOperation);
int CommunicationAllReduceValues(Eint64 *Values, int Number, 
				 MPI_Op ReduceOperation);
#endif /* USE_MPI */

int CommunicationBarrier();
int CommunicationShouldExit(int FromProc, int ToProc);

#endif
