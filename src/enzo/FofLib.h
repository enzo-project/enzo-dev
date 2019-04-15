/* foflib.h
   Mark Krumholz, 3/20/00
   Modified by Mark Krumholz, 8/22/00
   Modfiied by Nathan Goldbaum, December 2011, to include in enzo
   Header file to accompany foflib.c */

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

int FofVar(int, FLOAT *, FLOAT *, int *, int **);
int FofVarList(int, FLOAT *, FLOAT *, int *, int **, int ***);
int Fof(int, FLOAT *, FLOAT, int *, int **);
int FofList(int, FLOAT *, FLOAT, int *, int **, int ***);
int FofPrune(int, int, int *, int **, int);
int FofListPrune(int, int, int *, int **, int ***, int);
void NearNeighbor(int, FLOAT *, int, int *);
void NearNeighborPartial(int, FLOAT *, int, int, int *, int *);
void FindNeighbor(int, FLOAT *, FLOAT, int ***, int *);
void FindNeighborPartial(int, FLOAT *, int, int *, FLOAT *, int ***, int *);
