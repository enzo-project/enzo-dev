#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"




void allocate_memory(void)
{
  int n;

  n= Nslab[ThisTask] + Nshadow[ThisTask];

  if(n>0)
    {
      if(!(P=(struct particle_data*) malloc(n*sizeof(struct particle_data))))
	{
	  fprintf(stderr,"failed to allocate memory. (A)\n");
	  MPI_Abort(MPI_COMM_WORLD, 1); exit(1);
	}

      P-= 1;
    }
}









