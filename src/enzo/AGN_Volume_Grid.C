#include "AGN_Volume_Grid.h"

// Constructor that takes the grid dimensions (ncells * ncells) and the
// left and right corners of the grid [x,y,z].
AGN_Volume_Grid::AGN_Volume_Grid(int* ndim, float* gl, float* gr) {
   // Set grid parameters
   for (int i = 0; i < 3; i++) {
      dim[i] = ndim[i];

      left_corner[i] = gl[i];
      right_corner[i] = gr[i];
      }

   // Set up the cell positions and dimensions
   // cell_left gives the left corner of each cell.
   for (int d = 0; d < 3; d++) {
      dx[d] = (right_corner[d] - left_corner[d]) / ((float)dim[d]);

      cell_left[d] = new float[dim[d]];

      for (int i = 0; i < dim[d]; i++)
         cell_left[d][i] = left_corner[d] + dx[d] * (float)i;
      }
}

AGN_Volume_Grid::~AGN_Volume_Grid() {
   for (int d = 0; d < 3; d++)
      delete[] cell_left[d];
}

// Count how many corners of the cell with indicies x, y, z are enclosed within
// the cylinder cyl.
int AGN_Volume_Grid::count_enclosed_corners(int x, int y, int z, AGN_Volume_Cylinder* cyl) {
   int n = 0;

   float c[3];

   for (int k = 0; k < 2; k++) {
      for (int j = 0; j < 2; j++) {
         for (int i = 0; i < 2; i++) {
            c[0] = cell_left[0][x] + (float)i * dx[0];
            c[1] = cell_left[1][y] + (float)j * dx[1];
            c[2] = cell_left[2][z] + (float)k * dx[2];

            if (cyl -> contains_point(c) == true)
               n++;
            }
         }
      }

   return n;
}

// Return a grid showing the volume in each cell enclosed by the cylinder cyl.
// Recurses up to max_level levels, refining by refinement each time.
float*** AGN_Volume_Grid::get_intersection_volume(AGN_Volume_Cylinder* cyl, int level, int max_level, int refinement) {
   // Set up the grid to hold the volume info
   float*** vol = new float** [dim[2]];

   for (int j = 0; j < dim[2]; j++) {
      vol[j] = new float* [dim[1]];

      for (int i = 0; i < dim[1]; i++)
         vol[j][i] = new float[dim[0]];
      }
  
   for (int k = 0; k < dim[2]; k++)
      for (int j = 0; j < dim[1]; j++)
         for (int i = 0; i < dim[0]; i++)
            vol[k][j][i] = 0.0;

   // Find the volume enclosed in each cell.
   int nc;

   for (int k = 0; k < dim[2]; k++) {
      for (int j = 0; j < dim[1]; j++) {
         for (int i = 0; i < dim[0]; i++) {
            nc = count_enclosed_corners(i, j, k, cyl);

            // At the max level, estimate the enclosed area based on the number
            // of corners enclosed by the cylinder.
            if (level == max_level)
               vol[k][j][i] = dx[0] * dx[1] * dx[2] * (float)nc / 8.0;

            // If all corners are enclosed, assume the entire cell is enclosed.
            else if (nc == 8)
               vol[k][j][i] = dx[0] * dx[1] * dx[2];

            // If no corners are enclosed, assume the cell doesn't intersect the cylinder.
            else if (nc == 0)
               vol[k][j][i] = 0.0;

            else {
               float rlc[3];
               float rrc[3];

               rlc[0] = cell_left[0][i];
               rlc[1] = cell_left[1][j];
               rlc[2] = cell_left[2][k];

               rrc[0] = cell_left[0][i] + dx[0];
               rrc[1] = cell_left[1][j] + dx[1];
               rrc[2] = cell_left[2][k] + dx[2];

               int ref[3];
               ref[0] = refinement;
               ref[1] = refinement;
               ref[2] = refinement;

               AGN_Volume_Grid* refined_grid = new AGN_Volume_Grid(ref, rlc, rrc);

               float*** rvol = refined_grid -> get_intersection_volume(cyl, level+1, max_level, refinement);

               for (int kk = 0; kk < refinement; kk++)
                  for (int jj = 0; jj < refinement; jj++)
                     for (int ii = 0; ii < refinement; ii++)
                        vol[k][j][i] += rvol[kk][jj][ii];

               // Clean up memory
               delete refined_grid;

               for (int jj = 0; jj < refinement; jj++){
                  for (int ii = 0; ii < refinement; ii++) {
                     delete[] rvol[jj][ii];
                     }
                  delete[] rvol[jj];
                  }

               delete[] rvol;
               } // End if 0 < nc < 8
               
            } // End loop over i
         } // End loop over j
      } // End loop over k

   return vol;
}
