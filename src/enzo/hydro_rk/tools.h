inline int array3d_index(int i, int j, int k, int dimx, int dimy, int dimz)
{
  return i+dimx*(j+k*dimy);
}

inline Eflt32 my_MIN(Eflt32 a, Eflt32 b, Eflt32 c)
{
  Eflt32 r;
  r = (a < b) ? a : b;
  r = (r < c) ? r : c;
  return r;
}

inline Eflt64 my_MIN(Eflt64 a, Eflt64 b, Eflt64 c)
{
  Eflt64 r;
  r = (a < b) ? a : b;
  r = (r < c) ? r : c;
  return r;
}

inline Eflt128 my_MIN(Eflt128 a, Eflt128 b, Eflt128 c)
{
  Eflt128 r;
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
