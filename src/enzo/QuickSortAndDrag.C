/***********************************************************************
/
/  SORT A LIST OF INTS AND DRAG A NUMBER OF FIELDS WITH IT
/
/  written by: Greg Bryan
/  date:       July, 1995
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
 
#define SWAP(A, n1, n2, tmp) {tmp = A[n1]; A[n1] = A[n2]; A[n2] = tmp;}

// List = 32-bit int, DragList3 = 32-bit int
 
void QuickSortAndDrag(Eint32 List[], int left, int right,
		      int NumberToDrag1, float *DragList1[],
		      int NumberToDrag2, FLOAT *DragList2[],
		      int NumberToDrag3, Eint32 *DragList3[])
{
 
  int i, n, last, leftright2;
  Eint32 temp1;
  float temp2;
  FLOAT temp3;
 
  if (left >= right)
    return;
 
  leftright2 = (left + right)/2;
 
  SWAP(List, left, leftright2, temp1)
  for (n = 0; n < NumberToDrag1; n++)
    SWAP(DragList1[n], left, leftright2, temp2)
  for (n = 0; n < NumberToDrag2; n++)
    SWAP(DragList2[n], left, leftright2, temp3)
  for (n = 0; n < NumberToDrag3; n++)
    SWAP(DragList3[n], left, leftright2, temp1)
 
  last = left;
 
  for (i = left+1; i <= right; i++)
    if (List[i] < List[left]) {
      last++;
      SWAP(List, last, i, temp1)
      for (n = 0; n < NumberToDrag1; n++)
	SWAP(DragList1[n], last, i, temp2)
      for (n = 0; n < NumberToDrag2; n++)
	SWAP(DragList2[n], last, i, temp3)
      for (n = 0; n < NumberToDrag3; n++)
	SWAP(DragList3[n], last, i, temp1)
    }
 
  SWAP(List, left, last, temp1)
  for (n = 0; n < NumberToDrag1; n++)
    SWAP(DragList1[n], left, last, temp2)
  for (n = 0; n < NumberToDrag2; n++)
    SWAP(DragList2[n], left, last, temp3)
  for (n = 0; n < NumberToDrag3; n++)
    SWAP(DragList3[n], left, last, temp1)
 
  QuickSortAndDrag(List, left  , last-1,
		   NumberToDrag1, DragList1, NumberToDrag2, DragList2,
		   NumberToDrag3, DragList3);
  QuickSortAndDrag(List, last+1, right ,
		   NumberToDrag1, DragList1, NumberToDrag2, DragList2,
		   NumberToDrag3, DragList3);
 
}

//----------------------------------------------------------------------
// List = 64-bit int, DragList3 = 32-bit int
 
void QuickSortAndDrag(Eint64 List[], int left, int right,
		      int NumberToDrag1, float *DragList1[],
		      int NumberToDrag2, FLOAT *DragList2[],
		      int NumberToDrag3, Eint32 *DragList3[])
{
 
  int i, n, last, leftright2;
  Eint32 temp1_32;
  Eint64 temp1_64;
  float temp2;
  FLOAT temp3;
 
  if (left >= right)
    return;
 
  leftright2 = (left + right)/2;
 
  SWAP(List, left, leftright2, temp1_64)
  for (n = 0; n < NumberToDrag1; n++)
    SWAP(DragList1[n], left, leftright2, temp2)
  for (n = 0; n < NumberToDrag2; n++)
    SWAP(DragList2[n], left, leftright2, temp3)
  for (n = 0; n < NumberToDrag3; n++)
    SWAP(DragList3[n], left, leftright2, temp1_32)
 
  last = left;
 
  for (i = left+1; i <= right; i++)
    if (List[i] < List[left]) {
      last++;
      SWAP(List, last, i, temp1_64)
      for (n = 0; n < NumberToDrag1; n++)
	SWAP(DragList1[n], last, i, temp2)
      for (n = 0; n < NumberToDrag2; n++)
	SWAP(DragList2[n], last, i, temp3)
      for (n = 0; n < NumberToDrag3; n++)
	SWAP(DragList3[n], last, i, temp1_32)
    }
 
  SWAP(List, left, last, temp1_64)
  for (n = 0; n < NumberToDrag1; n++)
    SWAP(DragList1[n], left, last, temp2)
  for (n = 0; n < NumberToDrag2; n++)
    SWAP(DragList2[n], left, last, temp3)
  for (n = 0; n < NumberToDrag3; n++)
    SWAP(DragList3[n], left, last, temp1_32)
 
  QuickSortAndDrag(List, left  , last-1,
		   NumberToDrag1, DragList1, NumberToDrag2, DragList2,
		   NumberToDrag3, DragList3);
  QuickSortAndDrag(List, last+1, right ,
		   NumberToDrag1, DragList1, NumberToDrag2, DragList2,
		   NumberToDrag3, DragList3);
 
}

//----------------------------------------------------------------------
// List = 32-bit int, DragList3 = 64-bit int

void QuickSortAndDrag(Eint32 List[], int left, int right,
		      int NumberToDrag1, float *DragList1[],
		      int NumberToDrag2, FLOAT *DragList2[],
		      int NumberToDrag3, Eint64 *DragList3[])
{
 
  int i, n, last, leftright2;
  Eint32 temp1_32;
  Eint64 temp1_64;
  float temp2;
  FLOAT temp3;
 
  if (left >= right)
    return;
 
  leftright2 = (left + right)/2;
 
  SWAP(List, left, leftright2, temp1_32)
  for (n = 0; n < NumberToDrag1; n++)
    SWAP(DragList1[n], left, leftright2, temp2)
  for (n = 0; n < NumberToDrag2; n++)
    SWAP(DragList2[n], left, leftright2, temp3)
  for (n = 0; n < NumberToDrag3; n++)
    SWAP(DragList3[n], left, leftright2, temp1_64)
 
  last = left;
 
  for (i = left+1; i <= right; i++)
    if (List[i] < List[left]) {
      last++;
      SWAP(List, last, i, temp1_32)
      for (n = 0; n < NumberToDrag1; n++)
	SWAP(DragList1[n], last, i, temp2)
      for (n = 0; n < NumberToDrag2; n++)
	SWAP(DragList2[n], last, i, temp3)
      for (n = 0; n < NumberToDrag3; n++)
	SWAP(DragList3[n], last, i, temp1_64)
    }
 
  SWAP(List, left, last, temp1_32)
  for (n = 0; n < NumberToDrag1; n++)
    SWAP(DragList1[n], left, last, temp2)
  for (n = 0; n < NumberToDrag2; n++)
    SWAP(DragList2[n], left, last, temp3)
  for (n = 0; n < NumberToDrag3; n++)
    SWAP(DragList3[n], left, last, temp1_64)
 
  QuickSortAndDrag(List, left  , last-1,
		   NumberToDrag1, DragList1, NumberToDrag2, DragList2,
		   NumberToDrag3, DragList3);
  QuickSortAndDrag(List, last+1, right ,
		   NumberToDrag1, DragList1, NumberToDrag2, DragList2,
		   NumberToDrag3, DragList3);
 
}

//----------------------------------------------------------------------
// List = 64-bit int, DragList3 = 64-bit int

void QuickSortAndDrag(Eint64 List[], int left, int right,
		      int NumberToDrag1, float *DragList1[],
		      int NumberToDrag2, FLOAT *DragList2[],
		      int NumberToDrag3, Eint64 *DragList3[])
{
 
  int i, n, last, leftright2;
  Eint64 temp1;
  float temp2;
  FLOAT temp3;
 
  if (left >= right)
    return;
 
  leftright2 = (left + right)/2;
 
  SWAP(List, left, leftright2, temp1)
  for (n = 0; n < NumberToDrag1; n++)
    SWAP(DragList1[n], left, leftright2, temp2)
  for (n = 0; n < NumberToDrag2; n++)
    SWAP(DragList2[n], left, leftright2, temp3)
  for (n = 0; n < NumberToDrag3; n++)
    SWAP(DragList3[n], left, leftright2, temp1)
 
  last = left;
 
  for (i = left+1; i <= right; i++)
    if (List[i] < List[left]) {
      last++;
      SWAP(List, last, i, temp1)
      for (n = 0; n < NumberToDrag1; n++)
	SWAP(DragList1[n], last, i, temp2)
      for (n = 0; n < NumberToDrag2; n++)
	SWAP(DragList2[n], last, i, temp3)
      for (n = 0; n < NumberToDrag3; n++)
	SWAP(DragList3[n], last, i, temp1)
    }
 
  SWAP(List, left, last, temp1)
  for (n = 0; n < NumberToDrag1; n++)
    SWAP(DragList1[n], left, last, temp2)
  for (n = 0; n < NumberToDrag2; n++)
    SWAP(DragList2[n], left, last, temp3)
  for (n = 0; n < NumberToDrag3; n++)
    SWAP(DragList3[n], left, last, temp1)
 
  QuickSortAndDrag(List, left  , last-1,
		   NumberToDrag1, DragList1, NumberToDrag2, DragList2,
		   NumberToDrag3, DragList3);
  QuickSortAndDrag(List, last+1, right ,
		   NumberToDrag1, DragList1, NumberToDrag2, DragList2,
		   NumberToDrag3, DragList3);
 
}
 
/* Float version of the above w/o FLOAT (should just have one version). */
 
void QuickSortAndDragFloat(float List[], int left, int right, int NumberToDrag,
                           float *DragList[])
{
 
  int i, n, last, leftright2;
  float temp;
 
  if (left >= right)
    return;
 
  leftright2 = (left + right)/2;
 
  SWAP(List, left, leftright2, temp)
  for (n = 0; n < NumberToDrag; n++)
    SWAP(DragList[n], left, leftright2, temp)
 
  last = left;
 
  for (i = left+1; i <= right; i++)
    if (List[i] < List[left]) {
      last++;
      SWAP(List, last, i, temp)
      for (n = 0; n < NumberToDrag; n++)
	SWAP(DragList[n], last, i, temp)
    }
 
  SWAP(List, left, last, temp)
  for (n = 0; n < NumberToDrag; n++)
    SWAP(DragList[n], left, last, temp)
 
  QuickSortAndDragFloat(List, left  , last-1, NumberToDrag, DragList);
  QuickSortAndDragFloat(List, last+1, right , NumberToDrag, DragList);
 
}
