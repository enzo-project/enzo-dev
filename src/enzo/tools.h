inline int array3d_index(int i, int j, int k, int dimx, int dimy, int dimz)
{
  return i+dimx*(j+k*dimy);
}

inline float my_MIN(float a, float b, float c)
{
  float r;
  r = (a < b) ? a : b;
  r = (r < c) ? r : c;
  return r;
}

inline double my_MIN(double a, double b, double c)
{
  double r;
  r = (a < b) ? a : b;
  r = (r < c) ? r : c;
  return r;
}

inline float my_MAX(float a, float b, float c)
{
  float r;
  r = (a > b) ? a : b;
  r = (r > c) ? r : c;
  return r;
}
