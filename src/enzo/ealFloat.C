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
#include "ealFloat.h"

void my_exit(int exit_status);

ealFloat::ealFloat(int size,
			 float start, float end,
			 const char *axis_type){

  Size = size;
  Array = new float[size];
  if(!strcmp(axis_type,"none")){
    Zero();
  }else if(!strcmp(axis_type,"linear")){
    Linear(start, end);
  }else if(!strcmp(axis_type,"log")){
    Log(start, end);
  }else{
    fprintf(stderr, "ealFloat::ealFloat, invalid axis type.\n");
    fprintf(stderr, "Given %s, valid arguments are \"log\" or \"linear\".\n", axis_type);
    my_exit(EXIT_FAILURE);
  }
}


ealFloat::ealFloat( const ealFloat & orig)
  :Size(orig.Size)
{
  Array = new float[Size];
  for(int i = 0; i < Size; i++)
    Array[i] = orig.Array[i];
}

ealFloat::~ealFloat(){
  delete[] Array;
}

void ealFloat::Print(const char *fmt, FILE *stream){
  if(stream){
    for(int i = 0; i < Size; i++)
      fprintf(stream, fmt, Array[i]);
  }else{
    for(int i = 0; i < Size; i++)
      printf(fmt, Array[i]);
  }
}

void ealFloat::PrintWithIndex(const char *fmt, FILE *stream){
  if(stream){
    for(int i = 0; i < Size; i++)
      fprintf(stream, fmt, i, Array[i]);
  }else{
    for(int i = 0; i < Size; i++)
      printf(fmt, i, Array[i]);
  }
}

void ealFloat::Zero(){

  for(int i = 0; i < Size; i++)
    Array[i] = 0.0;
}

void ealFloat::Log(float start, float end){

  int i;
  float step = 1./((float)(Size-1))*log10(end/start);

  Array[0] = start;
  Array[Size - 1] = end;

  for( i=1; i < (Size - 1); i++)
    Array[i] = Array[i-1]*POW(10., step); 
}


void ealFloat::Linear(float start, float end){

  int i;
  float step = (end - start)/((float)(Size-1));

  Array[0] = start;
  Array[Size - 1] = end;

  for( i=1; i < (Size - 1); i++)
    Array[i] = i*step + Array[0]; 
}

float ealFloat::Sum(){
  float result = 0;

  for( int i=0; i < Size; i++)
    result+=Array[i];

  return result;
}

float ealFloat::Min(){
  float min_element = Array[0];

  for( int i=1; i < Size; i++)
    min_element = min( min_element, Array[i] );

  return min_element;
}

float ealFloat::Max(){
  float max_element = Array[0];

  for( int i=1; i < Size; i++)
    max_element = max( max_element, Array[i] );

  return max_element;
}

void ealFloat::ReduceSum(){

#ifdef USE_MPI
  float *RecvBuffer = new float[Size];
  Eint32 array_size = (Eint32) Size;
  int mpi_err = MPI_Allreduce( Array, RecvBuffer, array_size,
			       FloatDataType, MPI_SUM, MPI_COMM_WORLD );

  if( mpi_err != MPI_SUCCESS ){
    fprintf(stderr, "ealFloat::ReduceSum, mpi_err = %"ISYM", exiting.\n", mpi_err);
    my_exit(EXIT_FAILURE);
  }

  for( int i = 0; i < Size; i++ )
    Array[i] = RecvBuffer[i];

  delete[] RecvBuffer;
#endif

}

void ealFloat::ReduceMin(){

#ifdef USE_MPI
  float *RecvBuffer = new float[Size];
  Eint32 array_size = (Eint32) Size;
  int mpi_err = MPI_Allreduce( Array, RecvBuffer, array_size,
			       FloatDataType, MPI_MIN, MPI_COMM_WORLD );

  if( mpi_err != MPI_SUCCESS ){
    fprintf(stderr, "ealFloat::ReduceMax, mpi_err = %"ISYM", exiting.\n", mpi_err);
    my_exit(EXIT_FAILURE);
  }

  for( int i = 0; i < Size; i++ )
    Array[i] = RecvBuffer[i];

  delete[] RecvBuffer;
#endif

}

void ealFloat::ReduceMax(){

#ifdef USE_MPI
  float *RecvBuffer = new float[Size];
  Eint32 array_size = (Eint32) Size;
  int mpi_err = MPI_Allreduce( Array, RecvBuffer, array_size,
			       FloatDataType, MPI_MAX, MPI_COMM_WORLD );

  if( mpi_err != MPI_SUCCESS ){
    fprintf(stderr, "ealFloat::ReduceMax, mpi_err = %"ISYM", exiting.\n", mpi_err);
    my_exit(EXIT_FAILURE);
  }

  for( int i = 0; i < Size; i++ )
    Array[i] = RecvBuffer[i];

  delete[] RecvBuffer;
#endif

}

float ealFloat::GlobalSum(){
#ifdef USE_MPI
  ealFloat temp( *this );
  temp.ReduceSum();
  return temp.Sum();
#else
  return this->Sum();
#endif
}

float ealFloat::GlobalMin(){
#ifdef USE_MPI
  ealFloat temp( *this );
  temp.ReduceMin();
  return temp.Min();
#else
  return this->Min();
#endif
}

float ealFloat::GlobalMax(){
#ifdef USE_MPI
  ealFloat temp( *this );
  temp.ReduceMax();
  return temp.Max();
#else
  return this->Max();
#endif
}

void ealFloat::Bcast(int FromProcessor){

#ifdef USE_MPI
  Eint32 array_size = (Eint32) Size;
  int mpi_err =   MPI_Bcast( Array, array_size, 
			     FloatDataType, FromProcessor,
			     MPI_COMM_WORLD);

  if( mpi_err != MPI_SUCCESS ){
    fprintf(stderr, "ealFloat::Bcast, mpi_err = %"ISYM", exiting.\n", mpi_err);
    my_exit(EXIT_FAILURE);
  }
#endif

}


float &ealFloat::operator[](int subscript){

  if(subscript < 0 || subscript >=Size){
    fprintf(stderr, "ealFloat: subscript %"ISYM" out of range.\n", subscript);
    my_exit(EXIT_FAILURE);
  }

  return Array[subscript];
}


float ealFloat::operator[](int subscript) const
{ 

  if(subscript < 0 || subscript >=Size){
    fprintf(stderr, "ealFloat: subscript %"ISYM" out of range.\n", subscript);
    my_exit(EXIT_FAILURE);
  }

  return Array[subscript];
}

const ealFloat &ealFloat::operator=(const ealFloat &right){

  if(Size != right.Size && right.Size != 1){
    fprintf(stderr, "ealFloat::assignment Size %"ISYM" != right.Size %"ISYM".\n", Size, right.Size);
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

const ealFloat &ealFloat::operator=(const float &right){

  for( int i = 0; i < Size; i++ )
      Array[i] = right;

  return *this;  
}


/********** Addition **********/
const ealFloat 
&ealFloat::operator+=(const ealFloat &right){

  if(Size != right.Size && right.Size != 1){
    fprintf(stderr, "ealFloat::Math Error Size %"ISYM" != right.Size %"ISYM".\n", Size, right.Size);
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

const ealFloat &ealFloat::operator+=(const float &right){

  for( int i = 0; i < Size; i++ )
    Array[i] += right;

  return *this;
}

/********** Subtraction **********/
const ealFloat 
&ealFloat::operator-=(const ealFloat &right){

  if(Size != right.Size && right.Size != 1){
    fprintf(stderr, "ealFloat::Math Error Size %"ISYM" != right.Size %"ISYM".\n", Size, right.Size);
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

const ealFloat &ealFloat::operator-=(const float &right){

  for( int i = 0; i < Size; i++ )
    Array[i] -= right;

  return *this;
}

/********** Multiplication **********/
const ealFloat 
&ealFloat::operator*=(const ealFloat &right){

  if(Size != right.Size && right.Size != 1){
    fprintf(stderr, "ealFloat::Math Error Size %"ISYM" != right.Size %"ISYM".\n", Size, right.Size);
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

const ealFloat &ealFloat::operator*=(const float &right){

  for( int i = 0; i < Size; i++ )
    Array[i] *= right;

  return *this;
}

/********** Division **********/
const ealFloat 
&ealFloat::operator/=(const ealFloat &right){

  if(Size != right.Size && right.Size != 1){
    fprintf(stderr, "ealFloat::Math Error Size %"ISYM" != right.Size %"ISYM".\n", Size, right.Size);
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

const ealFloat &ealFloat::operator/=(const float &right){

  for( int i = 0; i < Size; i++ )
    Array[i] /= right;

  return *this;
}


/******** Math Ops ********/
const ealFloat 
ealFloat::operator+(const ealFloat &right){
  ealFloat result( *this );
  result += right;
  return result;
}

const ealFloat ealFloat::operator+(const float &right){
  ealFloat result( *this );
  result += right;
  return result;
}

const ealFloat 
ealFloat::operator-(const ealFloat &right){
  ealFloat result( *this );
  result -= right;
  return result;
}

const ealFloat ealFloat::operator-(const float &right){
  ealFloat result( *this );
  result -= right;
  return result;
}

const ealFloat 
ealFloat::operator*(const ealFloat &right){
  ealFloat result( *this );
  result *= right;
  return result;
}

const ealFloat ealFloat::operator*(const float &right){
  ealFloat result( *this );
  result *= right;
  return result;
}

const ealFloat 
ealFloat::operator/(const ealFloat &right){
  ealFloat result( *this );
  result /= right;
  return result;
}

const ealFloat ealFloat::operator/(const float &right){
  ealFloat result( *this );
  result /= right;
  return result;
}
