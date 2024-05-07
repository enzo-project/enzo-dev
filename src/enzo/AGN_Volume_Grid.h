#ifndef AGN_GRID_H
#define AGN_GRID_H

#include "preincludes.h"
#include "macros_and_parameters.h"
#include "AGN_Volume_Cylinder.h"

class AGN_Volume_Grid {
   public:
      AGN_Volume_Grid(int*, float*, float*);
      ~AGN_Volume_Grid();

      float*** get_intersection_volume(AGN_Volume_Cylinder*, int, int, int);
   
   private:
      int count_enclosed_corners(int, int, int, AGN_Volume_Cylinder*);

      int dim[3];
      float left_corner[3];
      float right_corner[3];

      float* cell_left[3];
      float dx[3];
};

#endif //AGN_GRID_H
