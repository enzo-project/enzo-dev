#include <stdio.h>
 
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
 
void icol(int *x, int n, int m, FILE *log_fptr)
{
 
int nrow,mrow;
int i,j;
 
#ifdef IO_LOG
  int         io_log = 1;
#else
  int         io_log = 0;
#endif
 
//if (io_log)
//{
  nrow = n/m;
  mrow = n - nrow * m;
  if( mrow > 0 )
  {
    nrow = nrow+1;
  }
 
  fprintf(log_fptr,"\n");
 
  for(j=0;j<n;j=j+m)
  {
    for(i=j;i<min(j+m,n);i++)
    {
      fprintf(log_fptr, "%1"ISYM" ", x[i]);
    }
    fprintf(log_fptr,"\n");
  }
//}
 
}
