/***********************************************************************
/
/  READ/WRITE LIST OF INTS/FLOATS
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
#include <stdio.h>
#include "ErrorExceptions.h"
 
#ifdef FAIL
#undef FAIL
#endif
#define FAIL      0
#define SUCCESS   1

 
int ReadListOfInts(FILE *fptr, int N, int nums[])
{
  for (int i = 0; i < N; i++)
    if (fscanf(fptr, "%d", nums + i) != 1) {
      fprintf(stderr, "ReadListOfInts called with %i numbers to read. %i %i %i %i %i",N,nums);
      throw(EnzoFatalException("Error in: "__FILE__));
 }
  fscanf(fptr, "\n");
  return SUCCESS;
}

int ReadListOfInts(FILE *fptr, long long int N, long long int nums[])
{
  for (long long int i = 0; i < N; i++)
    if (fscanf(fptr, "%lld", nums + i) != 1)
      throw(EnzoFatalException("Error in: "__FILE__));

  fscanf(fptr, "\n");
  return SUCCESS;
}
 
void WriteListOfInts(FILE *fptr, int N, int nums[])
{
  for (int i = 0; i < N; i++)
    fprintf(fptr, "%d ", nums[i]);
  fprintf(fptr, "\n");
}

void WriteListOfInts(FILE *fptr, long long int N, long long int nums[])
{
  for (long long int i = 0; i < N; i++)
    fprintf(fptr, "%lld ", nums[i]);
  fprintf(fptr, "\n");
}
 
int ReadListOfFloats(FILE *fptr, int N, float floats[])
{
  for (int i = 0; i < N; i++)
    if (fscanf(fptr, "%f", floats + i) != 1)
      throw(EnzoFatalException("Error in: "__FILE__));
 
  fscanf(fptr, "\n");
  return SUCCESS;
}

long long int ReadListOfFloats(FILE *fptr, long long int N, float floats[])
{
  for (long long int i = 0; i < N; i++)
    if (fscanf(fptr, "%f", floats + i) != 1)
      throw(EnzoFatalException("Error in: "__FILE__));

  fscanf(fptr, "\n");
  return SUCCESS;
}
 
void WriteListOfFloats(FILE *fptr, int N, float floats[])
{
  for (int i = 0; i < N; i++)
    fprintf(fptr, "%.7g ", floats[i]);
  fprintf(fptr, "\n");
}

void WriteListOfFloats(FILE *fptr, long long int N, float floats[])
{
  for (long long int i = 0; i < N; i++)
    fprintf(fptr, "%.7g ", floats[i]);
  fprintf(fptr, "\n");
} 


 
void WriteListOfFloats(FILE *fptr, int N, double floats[])
{
  for (int i = 0; i < N; i++)
    fprintf(fptr, "%.14g ", floats[i]);
  fprintf(fptr, "\n");
}

void WriteListOfFloats(FILE *fptr, long long int N, double floats[])
{
  for (long long int i = 0; i < N; i++)
    fprintf(fptr, "%.14g ", floats[i]);
  fprintf(fptr, "\n");
}
 
int ReadListOfFloats(FILE *fptr, int N, double floats[])
{
  for (int i = 0; i < N; i++)
    if (fscanf(fptr, "%lf", floats + i) != 1)
      throw(EnzoFatalException("Error in: "__FILE__));
 
  fscanf(fptr, "\n");
  return SUCCESS;
}
 
long long int ReadListOfFloats(FILE *fptr, long long int N, double floats[])
{
  for (long long int i = 0; i < N; i++)
    if (fscanf(fptr, "%lf", floats + i) != 1)
      throw(EnzoFatalException("Error in: "__FILE__));

  fscanf(fptr, "\n");
  return SUCCESS;
}




void WriteListOfFloats(FILE *fptr, int N, long double floats[])
{
  for (int i = 0; i < N; i++)
    fprintf(fptr, "%.21Lg ", floats[i]);
  fprintf(fptr, "\n");
}

void WriteListOfFloats(FILE *fptr, long long int N, long double floats[])
{
  for (long long int i = 0; i < N; i++)
    fprintf(fptr, "%.21Lg ", floats[i]);
  fprintf(fptr, "\n");
}


 
int ReadListOfFloats(FILE *fptr, int N, long double floats[])
{
  for (int i = 0; i < N; i++)
    if (fscanf(fptr, "%Lf", floats + i) != 1)
      throw(EnzoFatalException("Error in: "__FILE__));
 
  fscanf(fptr, "\n");
  return SUCCESS;
 
}

long long int ReadListOfFloats(FILE *fptr, long long int N, long double floats[])
{
  for (long long int i = 0; i < N; i++)
    if (fscanf(fptr, "%Lf", floats + i) != 1)
      throw(EnzoFatalException("Error in: "__FILE__));

  fscanf(fptr, "\n");
  return SUCCESS;

}

