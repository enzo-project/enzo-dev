// kernel for transpose 3d fields
template <int direction>
__global__ void transpose_cuda(float *dst, const float* __restrict src, int idim, int jdim, int kdim)
{
  int tid = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x;
  const int size = idim*jdim*kdim;
  if (tid < size) {
    int i, j, k;
    if (direction == 1) {
      j = tid % jdim;
      k = (tid % (jdim*kdim)) / jdim;
      i = tid / (jdim*kdim);
    } 
    
    if (direction == 2) {
      k = tid % kdim;
      i = (tid % (kdim*idim)) / kdim;
      j = tid / (kdim*idim);
    }
    
    int idx3 = i + (j + k*jdim)*idim;
    dst[tid] = src[idx3];
  }
}

//
// Tranpose 3d fields
// idim, jdim, kdim: x,y,z dim of the input array a
// direction = 1: x->y, y->z, z->x
// direction = 2: x->z, y->x, z->y
//
template<int direction>
void Transpose(float *dst, const float *src, int idim, int jdim, int kdim)
{
  const int size = idim*jdim*kdim;
  dim3 block;
  dim3 grid;
  block.x = 256;
  grid.x = (size + block.x - 1) / block.x;
  if (grid.x > 65535) {
    grid.y = (grid.x + 255)/256;
    grid.x = 256;
  }

  transpose_cuda<direction><<<grid, block>>>(dst, src, idim, jdim, kdim);
}
