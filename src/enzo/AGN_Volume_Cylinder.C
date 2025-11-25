#include "AGN_Volume_Cylinder.h"

AGN_Volume_Cylinder::AGN_Volume_Cylinder(float* b, float* n, float r, float h) {
   // Check that the parameters make sense
   if (r <= 0.0) {
      printf("Error! Cylinder radius must be positive.\n");
      exit(1);
      }

   if (h <= 0.0) {
      printf("Error! Cylinder height must be positive.\n");
      exit(1);
      }

   // Set cylinder parameters
   for (int i = 0; i < 3; i++) {
      base_pos[i] = b[i];
      nvec[i] = n[i];
      }

   radius = r;
   height = h;

   // Normalize nvec
   float ln = 0.0;

   for (int i = 0; i < 3; i++)
      ln += nvec[i] * nvec[i];

   if (ln <= 0.0) {
      printf("Error! Magnitude of nvec must be positive.\n");
      exit(1);
      }

   ln = sqrt(ln);

   for (int i = 0; i < 3; i++)
      nvec[i] /= ln;
}

// Returns true if the point x is within the cylinder.
// Returns false otherwise.
bool AGN_Volume_Cylinder::contains_point(float* x) {
   float bmx[3];
   float a[3];
   float l[3];

   // Compute bmx = b - x
   for (int i = 0; i < 3; i++)
      bmx[i] = base_pos[i] - x[i];

   // Compute a, the projection of bmx onto n
   // This is the vector from the point on the cylinder closest to x
   // to the base of the cylinder.
   float la = 0.0;

   for (int i = 0; i < 3; i++)
      la += bmx[i] * nvec[i];

   for (int i = 0; i < 3; i++)
      a[i] = nvec[i] * la;

   // Compute l, the vector from x to the closest point on th cylinder
   for (int i = 0; i < 3; i++)
      l[i] = bmx[i] - a[i];

   // Check if the closest point is above the base of the cylinder
   float adn = 0.0;

   for (int i = 0; i < 3; i++)
      adn += a[i] * nvec[i];

   if (adn > 0.0)
      return false;

   // Check if the closest point is below the top of the cylinder
   float magla = fabs(la);

   if (magla > height)
      return false;

   // Check if x is within the cylinder radius
   float ll = 0.0;

   for (int i = 0; i < 3; i++)
      ll += l[i] * l[i];

   ll = sqrt(ll);

   if (ll < radius)
      return true;

   else
      return false;
}
