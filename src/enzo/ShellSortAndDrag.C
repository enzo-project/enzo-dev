/***********************************************************************
/
/  (SHELL)SORT A LIST OF INTS AND DRAG A NUMBER OF FIELDS WITH IT
/
/  written by: John Wise
/  date:       May, 2010
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
 
#define SWAP(A, n1, n2, tmp) {tmp = A[n1]; A[n1] = A[n2]; A[n2] = tmp;}

// List = 32-bit int, DragList3 = 32-bit int
 
void ShellSortAndDrag(Eint32 List[], int N,
		      int NumberToDrag1, float *DragList1[],
		      int NumberToDrag2, FLOAT *DragList2[],
		      int NumberToDrag3, Eint32 *DragList3[])
{

  int a, nn, i, j, m;
  Eint32 temp1, t11, t12[NumberToDrag3];
  float temp2, t2[NumberToDrag1];
  FLOAT temp3, t3[NumberToDrag2];

  m = N;
  for (nn = 1; nn <= N; nn <<= 1) {
    m >>= 1;
    for (j = m; j < N; j++) {
      i = j - m;

      // Store j-element in temp. variables
      t11 = List[j];
      for (a = 0; a < NumberToDrag1; a++)
	t2[a] = DragList1[a][j];
      for (a = 0; a < NumberToDrag2; a++)
	t3[a] = DragList2[a][j];
      for (a = 0; a < NumberToDrag3; a++)
	t12[a] = DragList3[a][j];

      while (i >= 0 && List[i] > t11) {
//	printf("P%d: i=%d, j=%d, m=%d, Swapping %d (type=%d) and %d (type=%d)\n",
//	       MyProcessorNumber, i,j,m, i+m, List[i+m], i, List[i]);
	SWAP(List, i+m, i, temp1)
	for (a = 0; a < NumberToDrag1; a++)
	  SWAP(DragList1[a], i+m, i, temp2)
	for (a = 0; a < NumberToDrag2; a++)
	  SWAP(DragList2[a], i+m, i, temp3)
	for (a = 0; a < NumberToDrag3; a++)
	  SWAP(DragList3[a], i+m, i, temp1)
	i -= m;
      }

      // Insert stored data back in i+m element
      List[i+m] = t11;
      for (a = 0; a < NumberToDrag1; a++)
	DragList1[a][i+m] = t2[a];
      for (a = 0; a < NumberToDrag2; a++)
	DragList2[a][i+m] = t3[a];
      for (a = 0; a < NumberToDrag3; a++)
	DragList3[a][i+m] = t12[a];

    } // ENDFOR j
  } // ENDFOR nn

  return;
 
}

//----------------------------------------------------------------------
// List = 64-bit int, DragList3 = 32-bit int
 
void ShellSortAndDrag(Eint64 List[], int N,
		      int NumberToDrag1, float *DragList1[],
		      int NumberToDrag2, FLOAT *DragList2[],
		      int NumberToDrag3, Eint32 *DragList3[])
{
 
  int a, nn, i, j, m;
  Eint32 temp1_32, t12[NumberToDrag3];
  Eint64 temp1_64, t11;
  float temp2, t2[NumberToDrag1];
  FLOAT temp3, t3[NumberToDrag2];

  m = N;
  for (nn = 1; nn <= N; nn <<= 1) {
    m >>= 1;
    for (j = m; j < N; j++) {
      i = j - m;

      // Store j-element in temp. variables
      t11 = List[j];
      for (a = 0; a < NumberToDrag1; a++)
	t2[a] = DragList1[a][j];
      for (a = 0; a < NumberToDrag2; a++)
	t3[a] = DragList2[a][j];
      for (a = 0; a < NumberToDrag3; a++)
	t12[a] = DragList3[a][j];

      while (i >= 0 && List[i] > t11) {
	SWAP(List, i+m, i, temp1_64);
	for (a = 0; a < NumberToDrag1; a++)
	  SWAP(DragList1[a], i+m, i, temp2);
	for (a = 0; a < NumberToDrag2; a++)
	  SWAP(DragList2[a], i+m, i, temp3);
	for (a = 0; a < NumberToDrag3; a++)
	  SWAP(DragList3[a], i+m, i, temp1_32);
	i -= m;
      }

      // Insert stored data back in i+m element
      List[i+m] = t11;
      for (a = 0; a < NumberToDrag1; a++)
	DragList1[a][i+m] = t2[a];
      for (a = 0; a < NumberToDrag2; a++)
	DragList2[a][i+m] = t3[a];
      for (a = 0; a < NumberToDrag3; a++)
	DragList3[a][i+m] = t12[a];

    } // ENDFOR j
  } // ENDFOR nn

  return;
}

//----------------------------------------------------------------------
// List = 32-bit int, DragList3 = 64-bit int

void ShellSortAndDrag(Eint32 List[], int N,
		      int NumberToDrag1, float *DragList1[],
		      int NumberToDrag2, FLOAT *DragList2[],
		      int NumberToDrag3, Eint64 *DragList3[])
{
 
  int a, nn, i, j, m;
  Eint32 temp1_32, t11;
  Eint64 temp1_64, t12[NumberToDrag3];
  float temp2, t2[NumberToDrag1];
  FLOAT temp3, t3[NumberToDrag2];

  m = N;
  for (nn = 1; nn <= N; nn <<= 1) {
    m >>= 1;
    for (j = m; j < N; j++) {
      i = j - m;

      // Store j-element in temp. variables
      t11 = List[j];
      for (a = 0; a < NumberToDrag1; a++)
	t2[a] = DragList1[a][j];
      for (a = 0; a < NumberToDrag2; a++)
	t3[a] = DragList2[a][j];
      for (a = 0; a < NumberToDrag3; a++)
	t12[a] = DragList3[a][j];

      while (i >= 0 && List[i] > t11) {
	SWAP(List, i+m, i, temp1_32);
	for (a = 0; a < NumberToDrag1; a++)
	  SWAP(DragList1[a], i+m, i, temp2);
	for (a = 0; a < NumberToDrag2; a++)
	  SWAP(DragList2[a], i+m, i, temp3);
	for (a = 0; a < NumberToDrag3; a++)
	  SWAP(DragList3[a], i+m, i, temp1_64);
	i -= m;
      }

      // Insert stored data back in i+m element
      List[i+m] = t11;
      for (a = 0; a < NumberToDrag1; a++)
	DragList1[a][i+m] = t2[a];
      for (a = 0; a < NumberToDrag2; a++)
	DragList2[a][i+m] = t3[a];
      for (a = 0; a < NumberToDrag3; a++)
	DragList3[a][i+m] = t12[a];

    } // ENDFOR j
  } // ENDFOR nn

  return;
}

//----------------------------------------------------------------------
// List = 64-bit int, DragList3 = 64-bit int

void ShellSortAndDrag(Eint64 List[], int N,
		      int NumberToDrag1, float *DragList1[],
		      int NumberToDrag2, FLOAT *DragList2[],
		      int NumberToDrag3, Eint64 *DragList3[])
{
 
  int a, nn, i, j, m;
  Eint64 temp1, t11, t12[NumberToDrag3];
  float temp2, t2[NumberToDrag1];
  FLOAT temp3, t3[NumberToDrag2];

  m = N;
  for (nn = 1; nn <= N; nn <<= 1) {
    m >>= 1;
    for (j = m; j < N; j++) {
      i = j - m;

      // Store j-element in temp. variables
      t11 = List[j];
      for (a = 0; a < NumberToDrag1; a++)
	t2[a] = DragList1[a][j];
      for (a = 0; a < NumberToDrag2; a++)
	t3[a] = DragList2[a][j];
      for (a = 0; a < NumberToDrag3; a++)
	t12[a] = DragList3[a][j];

      while (i >= 0 && List[i] > t11) {
	SWAP(List, i+m, i, temp1);
	for (a = 0; a < NumberToDrag1; a++)
	  SWAP(DragList1[a], i+m, i, temp2);
	for (a = 0; a < NumberToDrag2; a++)
	  SWAP(DragList2[a], i+m, i, temp3);
	for (a = 0; a < NumberToDrag3; a++)
	  SWAP(DragList3[a], i+m, i, temp1);
	i -= m;
      }

      // Insert stored data back in i+m element
      List[i+m] = t11;
      for (a = 0; a < NumberToDrag1; a++)
	DragList1[a][i+m] = t2[a];
      for (a = 0; a < NumberToDrag2; a++)
	DragList2[a][i+m] = t3[a];
      for (a = 0; a < NumberToDrag3; a++)
	DragList3[a][i+m] = t12[a];

    } // ENDFOR j
  } // ENDFOR nn

  return; 
}
