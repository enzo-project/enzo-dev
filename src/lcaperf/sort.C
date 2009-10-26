//----------------------------------------------------------------------
//
// File: jb.C
//
// Description: Functions used by jb-*.C utilities.  See jb.h for the list
// Description: of functions
//
//----------------------------------------------------------------------
//
// Copyright 2005 James Bordner
// Copyright 2005 Laboratory for Computational Astrophysics
// Copyright 2005 Regents of the University of California
//
//----------------------------------------------------------------------

#include "jb.h"

//------------------------------------------------------------------------

int partition_x(vector<struct reg_struct *> &y, int f, int l) {
  int up,down;
  int ipiv = int(float(rand())/RAND_MAX*(l-f-1)+f);
  struct reg_struct *piv = y[ipiv], *temp;
  swap(y[ipiv],y[f]);
  up = f;
  down = l;
  do { 
    while (y[up]->x <= piv->x && up < l) {
      up++;
    }
    while (y[down]->x > piv->x  && down > f ) {
      down--;
    }
    if (up < down ) swap(y[up],y[down]);
  } while (down > up);
  y[f] = y[down];
  y[down] = piv;
  return down;
}

void jb_sort_x(vector<struct reg_struct *> &x, int first, int last) {
  int pivIndex = 0;
  if(first < last) {
    pivIndex = partition_x(x,first, last);
    jb_sort_x(x,first,(pivIndex-1));
    jb_sort_x(x,(pivIndex+1),last);
  }
}

//------------------------------------------------------------------------

int partition_sorted(vector<struct reg_struct *> &y, int f, int l) {
  int up,down;
  int ipiv = int(float(rand())/RAND_MAX*(l-f-1)+f);
  struct reg_struct *piv = y[ipiv], *temp;
  swap(y[ipiv],y[f]);
  up = f;
  down = l;
  do { 
    while (y[up]->sorted <= piv->sorted && up < l) {
      up++;
    }
    while (y[down]->sorted > piv->sorted  && down > f ) {
      down--;
    }
    if (up < down ) swap(y[up],y[down]);
  } while (down > up);
  y[f] = y[down];
  y[down] = piv;
  return down;
}

void jb_sort_sorted(vector<struct reg_struct *> &x, int first, int last) {
  int pivIndex = 0;
  if(first < last) {
    pivIndex = partition_sorted(x,first, last);
        jb_sort_sorted(x,first,(pivIndex-1));
        jb_sort_sorted(x,(pivIndex+1),last);
    }
}
