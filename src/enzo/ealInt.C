#ifdef USE_MPI
#include <mpi.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"

#include "ealInt.h"

void my_exit(int exit_status);

ealInt::ealInt(int size){
  Size = size;
  Array = new int[size];
  Zero();
}

ealInt::ealInt( const ealInt & orig)
  :Size(orig.Size)
{
  Array = new int[Size];
  for(int i = 0; i < Size; i++)
    Array[i] = orig.Array[i];
}

ealInt::~ealInt(){
  delete[] Array;
}

void ealInt::Print(const char *fmt, FILE *stream){
  if(stream){
    for(int i = 0; i < Size; i++)
      fprintf(stream, fmt, Array[i]);
  }else{
    for(int i = 0; i < Size; i++)
      printf(fmt, Array[i]);
  }
}


void ealInt::PrintWithIndex(const char *fmt, FILE *stream){
  if(stream){
    for(int i = 0; i < Size; i++)
      fprintf(stream, fmt, i, Array[i]);
  }else{
    for(int i = 0; i < Size; i++)
      printf(fmt, i, Array[i]);
  }
}


void ealInt::Zero(){

  for(int i = 0; i < Size; i++)
    Array[i] = 0;
}

int ealInt::Sum(){
  int result = 0;

  for( int i = 0; i < Size; i++)
    result+=Array[i];

  return result;
}

int ealInt::Min(){
  int min_element = Array[0];

  for( int i=1; i < Size; i++)
    min_element = min( min_element, Array[i] );

  return min_element;
}

int ealInt::Max(){
  int max_element = Array[0];

  for( int i=1; i < Size; i++)
    max_element = max( max_element, Array[i] );

  return max_element;
}

void ealInt::ReduceSum(){

#ifdef USE_MPI
  int *RecvBuffer = new int[Size];

  int mpi_err = MPI_Allreduce( Array, RecvBuffer, Size,
			       IntDataType, MPI_SUM, MPI_COMM_WORLD );

  if( mpi_err != MPI_SUCCESS ){
    fprintf(stderr, "ealInt::ReduceSum, mpi_err = %"ISYM", exiting.\n", mpi_err);
    my_exit(EXIT_FAILURE);
  }

  for( int i = 0; i < Size; i++ )
    Array[i] = RecvBuffer[i];

  delete[] RecvBuffer;
#endif

}

void ealInt::ReduceMin(){

#ifdef USE_MPI
  int *RecvBuffer = new int[Size];
  Eint32 array_size = (Eint32) Size;
  int mpi_err = MPI_Allreduce( Array, RecvBuffer, array_size,
			       IntDataType, MPI_MIN, MPI_COMM_WORLD );

  if( mpi_err != MPI_SUCCESS ){
    fprintf(stderr, "ealInt::ReduceMax, mpi_err = %"ISYM", exiting.\n", mpi_err);
    my_exit(EXIT_FAILURE);
  }

  for( int i = 0; i < Size; i++ )
    Array[i] = RecvBuffer[i];

  delete[] RecvBuffer;
#endif

}

void ealInt::ReduceMax(){

#ifdef USE_MPI
  int *RecvBuffer = new int[Size];
  Eint32 array_size = (Eint32) Size;
  int mpi_err = MPI_Allreduce( Array, RecvBuffer, array_size,
			       IntDataType, MPI_MAX, MPI_COMM_WORLD );

  if( mpi_err != MPI_SUCCESS ){
    fprintf(stderr, "ealInt::ReduceMax, mpi_err = %"ISYM", exiting.\n", mpi_err);
    my_exit(EXIT_FAILURE);
  }

  for( int i = 0; i < Size; i++ )
    Array[i] = RecvBuffer[i];

  delete[] RecvBuffer;
#endif

}

int ealInt::GlobalSum(){
#ifdef USE_MPI
  ealInt temp( *this );
  temp.ReduceSum();
  return temp.Sum();
#else
  return this->Sum();
#endif
}

int ealInt::GlobalMin(){
#ifdef USE_MPI
  ealInt temp( *this );
  temp.ReduceMin();
  return temp.Min();
#else
  return this->Min();
#endif
}

int ealInt::GlobalMax(){
#ifdef USE_MPI
  ealInt temp( *this );
  temp.ReduceMax();
  return temp.Max();
#else
  return this->Max();
#endif
}

void ealInt::Bcast(int FromProcessor){

#ifdef USE_MPI
  Eint32 array_size = (Eint32) Size;
  int mpi_err =   MPI_Bcast( Array, array_size, 
			     IntDataType, FromProcessor,
			     MPI_COMM_WORLD);

  if( mpi_err != MPI_SUCCESS ){
    fprintf(stderr, "ealInt::Bcast, mpi_err = %"ISYM", exiting.\n", mpi_err);
    my_exit(EXIT_FAILURE);
  }
#endif

}


int &ealInt::operator[](int subscript){

  if(subscript < 0 || subscript >= Size){
    fprintf(stderr, "ealInt: subscript %"ISYM" out of range. Size = %"ISYM".\n", subscript, Size);
    my_exit(EXIT_FAILURE);
  }

  return Array[subscript];
}


int ealInt::operator[](int subscript) const
{ 

  if(subscript < 0 || subscript >=Size){
    fprintf(stderr, "ealInt: subscript %"ISYM" out of range. Size = %"ISYM".\n", subscript, Size);
    my_exit(EXIT_FAILURE);
  }

  return Array[subscript];
}

const ealInt &ealInt::operator=(const ealInt &right){

  if(Size != right.Size && right.Size != 1){
    fprintf(stderr, "ealInt::assignment Size %"ISYM" != right.Size %"ISYM".\n", Size, right.Size);
    my_exit(EXIT_FAILURE);
  }

  if( &right != this ){
    if( right.Size != 1){
      for( int i = 0; i < Size; i++ )
	Array[i] = right.Array[i];
    }else{
      for( int i = 0; i < Size; i++ )
	Array[i] = right.Array[0];
    }
  }

  return *this;
}

const ealInt &ealInt::operator=(const int &right){

  for( int i = 0; i < Size; i++ )
      Array[i] = right;

  return *this;  
}


/********** Addition **********/
const ealInt 
&ealInt::operator+=(const ealInt &right){

  if(Size != right.Size && right.Size != 1){
    fprintf(stderr, "ealInt::Math Error Size %"ISYM" != right.Size %"ISYM".\n", Size, right.Size);
    my_exit(EXIT_FAILURE);
  }

  if( &right != this ){
    if( right.Size != 1){
      for( int i = 0; i < Size; i++ )
	Array[i] += right.Array[i];
    }else{
      for( int i = 0; i < Size; i++ )
	Array[i] += right.Array[0];
    }
  }
  
  return *this;
}

const ealInt &ealInt::operator+=(const int &right){

  for( int i = 0; i < Size; i++ )
    Array[i] += right;

  return *this;
}

/********** Subtraction **********/
const ealInt 
&ealInt::operator-=(const ealInt &right){

  if(Size != right.Size && right.Size != 1){
    fprintf(stderr, "ealInt::Math Error Size %"ISYM" != right.Size %"ISYM".\n", Size, right.Size);
    my_exit(EXIT_FAILURE);
  }

  if( &right != this ){
    if( right.Size != 1){
      for( int i = 0; i < Size; i++ )
	Array[i] -= right.Array[i];
    }else{
      for( int i = 0; i < Size; i++ )
	Array[i] -= right.Array[0];
    }
  }
  
  return *this;
}

const ealInt &ealInt::operator-=(const int &right){

  for( int i = 0; i < Size; i++ )
    Array[i] -= right;

  return *this;
}

/********** Multiplication **********/
const ealInt 
&ealInt::operator*=(const ealInt &right){

  if(Size != right.Size && right.Size != 1){
    fprintf(stderr, "ealInt::Math Error Size %"ISYM" != right.Size %"ISYM".\n", Size, right.Size);
    my_exit(EXIT_FAILURE);
  }

  if( &right != this ){
    if( right.Size != 1){
      for( int i = 0; i < Size; i++ )
	Array[i] *= right.Array[i];
    }else{
      for( int i = 0; i < Size; i++ )
	Array[i] *= right.Array[0];
    }
  }
  
  return *this;
}

const ealInt &ealInt::operator*=(const int &right){

  for( int i = 0; i < Size; i++ )
    Array[i] *= right;

  return *this;
}

/********** Division **********/
const ealInt 
&ealInt::operator/=(const ealInt &right){

  if(Size != right.Size && right.Size != 1){
    fprintf(stderr, "ealInt::Math Error Size %"ISYM" != right.Size %"ISYM".\n", Size, right.Size);
    my_exit(EXIT_FAILURE);
  }

  if( &right != this ){
    if( right.Size != 1){
      for( int i = 0; i < Size; i++ )
	Array[i] /= right.Array[i];
    }else{
      for( int i = 0; i < Size; i++ )
	Array[i] /= right.Array[0];
    }
  }
  
  return *this;
}

const ealInt &ealInt::operator/=(const int &right){

  for( int i = 0; i < Size; i++ )
    Array[i] /= right;

  return *this;
}


/******** Math Ops ********/
const ealInt 
ealInt::operator+(const ealInt &right){
  ealInt result( *this );
  result += right;
  return result;
}

const ealInt ealInt::operator+(const int &right){
  ealInt result( *this );
  result += right;
  return result;
}

const ealInt 
ealInt::operator-(const ealInt &right){
  ealInt result( *this );
  result -= right;
  return result;
}

const ealInt ealInt::operator-(const int &right){
  ealInt result( *this );
  result -= right;
  return result;
}

const ealInt 
ealInt::operator*(const ealInt &right){
  ealInt result( *this );
  result *= right;
  return result;
}

const ealInt ealInt::operator*(const int &right){
  ealInt result( *this );
  result *= right;
  return result;
}

const ealInt 
ealInt::operator/(const ealInt &right){
  ealInt result( *this );
  result /= right;
  return result;
}

const ealInt ealInt::operator/(const int &right){
  ealInt result( *this );
  result /= right;
  return result;
}
