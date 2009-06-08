#include<stdio.h>

#include "macros_and_parameters.h"

Eint64 FreeRealMem( char *node );


int main()
{

  char node_name[16];
  Eint64 mem;

  strcpy(node_name, "ds007");

  mem = FreeRealMem(node_name);

  fprintf(stderr, "MEM = %lld\n", mem);

}
