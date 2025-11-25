#ifndef AGN_CYLINDER_H
#define AGN_CYLINDER_H

#include "preincludes.h"
#include "macros_and_parameters.h"
#include "AGN_Volume_Cylinder.h"

class AGN_Volume_Cylinder {
   public:
      AGN_Volume_Cylinder(float*, float*, float, float);
      bool contains_point(float*);
      
   private:
      float base_pos[3];
      float nvec[3];
      float radius;
      float height;
};

#endif //AGN_CYLINDER_H
