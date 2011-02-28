/***********************************************************************
/
/  COMPUTE MHD FLUX USING GPU
/
/  written by: Peng Wang, KIPAC/Stanford
/  date:       January, 2009
/  modified1:
/
/  description: 
/      This is the heart of the MHD solver where fluxes at cell
/    interfaces are calculated using primitive variables at cell centers. 
/    The reconstruction is done using piecewise linear method (PLM) and 
/    the Riemann solver is Local Lax-Friedrichs (LLF) solver.
/      The global strcucture of the CPU code implementation and algorithms 
/    are described more in Wang, Abel & Zhang, 2008, ApJS, 176, 467:
/       http://www.iop.org/EJ/abstract/0067-0049/176/2/467/
/       or
/       http://arxiv.org/abs/astro-ph/0703742
/      The GPU implementation uses one thread per cell. The key goal
/    is to allow arbitrary grid dimension in the parallization model, 
/    which is crucial for adaptive mesh calculations where arbitrary grid 
/    size can emerge. This is achieved by two techniques: first, rearranging 
/    the input data and appending padding cells if the size does not fit. 
/    Second, with the use of method of line (MOL) spatial
/    discretization, flux calculation is made purely local (only need 4 
/    continuous cells for a flux is needed). This also makes the computational
/    model easily extend to 3D.
/
************************************************************************/

#include "../macros_and_parameters.h"
#include "../typedefs.h"
#include "../global_data.h"
#include <stdio.h>
#include <cutil.h>

// hack for making things compile
#define CUDA_BLOCK_SIZE 64
#define CUDA_GRID_SIZE 640

//#define GHOST_SIZE 4
#define PRINT_CUDA_TIMING 0
//#define Gamma 2.0
//#define Theta_Limiter 1.5
#define EOSType 3
#define EOSSoundSpeed 1.
#define NEQ_MHD 9
#define iD 0
#define iS1 1
#define iS2 2
#define iS3 3
#define iEtot 4
#define iBx 5
#define iBy 6
#define iBz 7
#define iPhi 8
#define irho 0
#define ieint 1
#define ivx 2
#define ivy 3
#define ivz 4
//#define ABS(a) ((a) >= 0 ? (a) : (-a))

// forward declarations
__global__ void ComputeInternalEnergy_kernel(float *Rho, float *Vx, float *Vy, float *Vz, float *Etot,
                                              float *Eneint, float *Bx, float *By, float *Bz, int size);
__global__ void MHDSweepX_CUDA3_kernel(float *Rho, float *Vx, float *Vy, float *Vz, float *Etot,
	    		               float *Eneint, float *Bx, float *By, float *Bz, float *Phi,
				       float *FluxD, float *FluxS1, float *FluxS2, float *FluxS3,
				       float *FluxTau, float *FluxBx, float *FluxBy, float *FluxBz,
				       float *FluxPhi, float C_h, float Gamma, int size, float cTheta_Limiter);
__global__ void MHDSweepY_CUDA3_kernel(float *Rho, float *Vx, float *Vy, float *Vz, float *Etot,
	    		               float *Eneint, float *Bx, float *By, float *Bz, float *Phi,
				       float *FluxD, float *FluxS1, float *FluxS2, float *FluxS3,
				       float *FluxTau, float *FluxBx, float *FluxBy, float *FluxBz,
				       float *FluxPhi, float C_h, float Gamma, 
				       int size, int dim0, int dim1, int dim2, float cTheta_Limiter);
__global__ void MHDSweepZ_CUDA3_kernel(float *Rho, float *Vx, float *Vy, float *Vz, float *Etot,
	    		               float *Eneint, float *Bx, float *By, float *Bz, float *Phi,
				       float *FluxD, float *FluxS1, float *FluxS2, float *FluxS3,
				       float *FluxTau, float *FluxBx, float *FluxBy, float *FluxBz,
				       float *FluxPhi, float C_h, float Gamma,
				       int size, int dim0, int dim1, int dim2, float cTheta_Limiter);
__global__ void MHDComputedUx_CUDA3_kernel(float *FluxD, float *FluxS1, float *FluxS2, float *FluxS3,
				           float *FluxTau, float *FluxBx, float *FluxBy, float *FluxBz, float *FluxPhi, 
					   float *dUD, float *dUS1, float *dUS2, float *dUS3,
				           float *dUTau, float *dUBx, float *dUBy, float *dUBz, float *dUPhi, 
					   float dtdx, int size);
__global__ void MHDComputedUy_CUDA3_kernel(float *FluxD, float *FluxS1, float *FluxS2, float *FluxS3,
				           float *FluxTau, float *FluxBx, float *FluxBy, float *FluxBz, float *FluxPhi, 
					   float *dUD, float *dUS1, float *dUS2, float *dUS3,
				           float *dUTau, float *dUBx, float *dUBy, float *dUBz, float *dUPhi, 
					   float dtdx, int size, int dim0, int dim1, int dim2);
__global__ void MHDComputedUz_CUDA3_kernel(float *FluxD, float *FluxS1, float *FluxS2, float *FluxS3,
				           float *FluxTau, float *FluxBx, float *FluxBy, float *FluxBz, float *FluxPhi, 
					   float *dUD, float *dUS1, float *dUS2, float *dUS3,
				           float *dUTau, float *dUBx, float *dUBy, float *dUBz, float *dUPhi, 
					   float dtdx, int size, int dim0, int dim1, int dim2);
__global__ void MHDUpdatePrim_CUDA3_kernel(float *Rho, float *Vx, float *Vy, float *Vz,
	  			           float *Etot, float *Bx, float *By, float *Bz, float *Phi, 
					   float *dUD, float *dUS1, float *dUS2, float *dUS3,
				           float *dUTau, float *dUBx, float *dUBy, float *dUBz, float *dUPhi, 
					   float dt, float C_h, float C_p, int size);

__device__ void LLF_PLM_MHD_CUDA3(float *Prim, float *Flux, const int &tx, 
				  const float &C_h, const float &cGamma, const float &cTheta_Limiter);
__device__ void plm_point(const float &vm1, const float &v, const float &vp1, float &vl_plm, const float &cTheta_Limiter);
__device__ void EOS(float &p, float &rho, float &e, float &cs, 
		    const float &cGamma, const int &eostype, const int &mode);
__device__ float minmod_mhd(const float &a, const float &b, const float &c);
__device__ float Max_mhd(const float &a, const float &b, const float &c);
__device__ float Min_mhd(const float &a, const float &b, const float &c);
double ReturnWallTime();

/******************************************************************************
/
/   Function MHDSweepX_CUDA3: CPU code wrap the data, send to GPU and call the kernel
/
/   Input : 
/     Prim[NEQ_MHD][GridDimension^3] : primitive variables at cell center
/     GridDimension[3] : grid dimension
/     GridStartIndex[3] : the starting index of active data
/     GridRank : dimension of the problem
/     dtdx : dt/dx
/     C_h : wave speed in Dedner's formulation
/   Output: 
/     Flux3D[NEQ_MHD][(activesize+1)^3] -- fluxes at cell interface
/
*******************************************************************************/

int MHDTimeUpdate_CUDA(float **Prim, int GridDimension[], int GridStartIndex[], int GridEndIndex[], int GridRank,
		       float dtdx, float dt, float C_h, float C_p, float cTheta_Limiter)
{

  // compute grid size
  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  cudaEvent_t start, stop;
  float elapsedTime;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);

  cudaEventRecord(start, 0);

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsedTime, start, stop);
  //  PerformanceTimers[40] += elapsedTime/1e3;
  if (PRINT_CUDA_TIMING) fprintf(stderr, "minimal measurable time:     %g \n" , elapsedTime/1e3);
  cudaEventRecord(start, 0);

#ifdef ECUDADEBUG
  printf("Theta_Limiter sent to device: %g\n", cTheta_Limiter);
  printf("Gridsize: %i\n", size);
  printf("NEQ_MHD: %i\n", NEQ_MHD);
  printf("Gridsize: %i %i %i\n", GridDimension[0], GridDimension[1], GridDimension[2]);
  for (int j=30; j < 33; j++) 
    for (int i=0; i < 9; i++) printf("Prim[%i][%i] = %g \n", i, j, Prim[i][j]);
#endif


  float *PrimDevice = 0;
  int totalsize=sizeof(float)*(NEQ_MHD+1)*size;
  if (cudaMalloc((void**)&PrimDevice, totalsize) != cudaSuccess) {
    printf("cudaMalloc for PrimDevice with size %d failed.\n", totalsize);
    return FAIL;
  }

  cudaMemcpy(PrimDevice,        Prim[0], sizeof(float)*size, cudaMemcpyHostToDevice);
  cudaMemcpy(PrimDevice+size,   Prim[1], sizeof(float)*size, cudaMemcpyHostToDevice);
  cudaMemcpy(PrimDevice+2*size, Prim[2], sizeof(float)*size, cudaMemcpyHostToDevice);
  cudaMemcpy(PrimDevice+3*size, Prim[3], sizeof(float)*size, cudaMemcpyHostToDevice);
  cudaMemcpy(PrimDevice+4*size, Prim[4], sizeof(float)*size, cudaMemcpyHostToDevice);
  cudaMemcpy(PrimDevice+5*size, Prim[5], sizeof(float)*size, cudaMemcpyHostToDevice);
  cudaMemcpy(PrimDevice+6*size, Prim[6], sizeof(float)*size, cudaMemcpyHostToDevice);
  cudaMemcpy(PrimDevice+7*size, Prim[7], sizeof(float)*size, cudaMemcpyHostToDevice);
  cudaMemcpy(PrimDevice+8*size, Prim[8], sizeof(float)*size, cudaMemcpyHostToDevice);

  // copy the tmp data to device memory
  float *Rho_Device  = PrimDevice, 
        *Vx_Device   = PrimDevice + size, 
	*Vy_Device   = PrimDevice + 2*size, 
	*Vz_Device   = PrimDevice + 3*size, 
	*Etot_Device = PrimDevice + 4*size, 
	*Bx_Device   = PrimDevice + 5*size, 
	*By_Device   = PrimDevice + 6*size, 
	*Bz_Device   = PrimDevice + 7*size, 
	*Phi_Device  = PrimDevice + 8*size,
	*Eint_Device = PrimDevice + 9*size;

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsedTime, start, stop);
  //  PerformanceTimers[34] += elapsedTime/1e3;
  if (PRINT_CUDA_TIMING) fprintf(stderr, "copying data to GPU took:    %g \n" , elapsedTime/1e3);

  cudaEventRecord(start, 0);

  float *FluxDevice = 0;
  if (cudaMalloc((void**)&FluxDevice, sizeof(float)*NEQ_MHD*size) != cudaSuccess) {
    cudaFree(PrimDevice);
    printf("cuda Malloc failed for FluxDevice with size %d\n", sizeof(float)*NEQ_MHD*size);
    return FAIL;
  }

  float *FluxD_Device   = FluxDevice, 
        *FluxS1_Device  = FluxDevice + size, 
	*FluxS2_Device  = FluxDevice + 2*size, 
	*FluxS3_Device  = FluxDevice + 3*size, 
	*FluxTau_Device = FluxDevice + 4*size, 
	*FluxBx_Device  = FluxDevice + 5*size, 
	*FluxBy_Device  = FluxDevice + 6*size, 
	*FluxBz_Device  = FluxDevice + 7*size, 
	*FluxPhi_Device = FluxDevice + 8*size;

  // allocate space for dU
  float *dUDevice = 0;
  if (cudaMalloc((void**)&dUDevice, sizeof(float)*NEQ_MHD*size) != cudaSuccess) {
    cudaFree(PrimDevice);
    cudaFree(FluxDevice);
    printf("cuda Malloc failed for dUDevice with size %d\n", sizeof(float)*NEQ_MHD*size);
    return FAIL;
  }

  float *dUD_Device   = dUDevice, 
        *dUS1_Device  = dUDevice + size, 
	*dUS2_Device  = dUDevice + 2*size, 
	*dUS3_Device  = dUDevice + 3*size, 
	*dUTau_Device = dUDevice + 4*size, 
	*dUBx_Device  = dUDevice + 5*size, 
	*dUBy_Device  = dUDevice + 6*size, 
	*dUBz_Device  = dUDevice + 7*size, 
	*dUPhi_Device = dUDevice + 8*size;

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsedTime, start, stop);
  if (PRINT_CUDA_TIMING) fprintf(stderr, "alloc space on device took:  %g \n" , elapsedTime/1e3);

  //  PerformanceTimers[35] += elapsedTime/1e3;

  // compute gridSize
  dim3 dimBlock(CUDA_BLOCK_SIZE);
  dim3 dimGrid;
  if (size <= CUDA_BLOCK_SIZE*CUDA_GRID_SIZE) {
    dimGrid.x = size/CUDA_BLOCK_SIZE+1;
    dimGrid.y = 1;
  } else {
    dimGrid.x = CUDA_GRID_SIZE;
    dimGrid.y = size/(CUDA_BLOCK_SIZE*CUDA_GRID_SIZE) + 1;
  }


  // call the kernel function to do the computation
  cudaEventRecord(start, 0);

  ComputeInternalEnergy_kernel<<<dimGrid, dimBlock>>>(Rho_Device, Vx_Device, Vy_Device, Vz_Device, Etot_Device,
                                                      Eint_Device, Bx_Device, By_Device, Bz_Device, size);
  cudaThreadSynchronize();
  MHDSweepX_CUDA3_kernel<<<dimGrid,dimBlock>>>(Rho_Device, Vx_Device, Vy_Device, Vz_Device, Etot_Device, 
					       Eint_Device, Bx_Device, By_Device, Bz_Device, Phi_Device,
					       FluxD_Device, FluxS1_Device, FluxS2_Device, FluxS3_Device,
					       FluxTau_Device, FluxBx_Device, FluxBy_Device, FluxBz_Device,
					       FluxPhi_Device, C_h, Gamma, size, cTheta_Limiter);


  // compute dU for x direction
  MHDComputedUx_CUDA3_kernel<<<dimGrid,dimBlock>>>(FluxD_Device, FluxS1_Device, FluxS2_Device, FluxS3_Device, FluxTau_Device,
  						   FluxBx_Device, FluxBy_Device, FluxBz_Device, FluxPhi_Device, 
					           dUD_Device, dUS1_Device, dUS2_Device, dUS3_Device,
				                   dUTau_Device, dUBx_Device, dUBy_Device, dUBz_Device, dUPhi_Device, 
					           dtdx, size);

  if (GridRank > 1) {
    MHDSweepY_CUDA3_kernel<<<dimGrid,dimBlock>>>(Rho_Device, Vx_Device, Vy_Device, Vz_Device, Etot_Device, 
	  				         Eint_Device, Bx_Device, By_Device, Bz_Device, Phi_Device,
		 			         FluxD_Device, FluxS1_Device, FluxS2_Device, FluxS3_Device,
					         FluxTau_Device, FluxBx_Device, FluxBy_Device, FluxBz_Device,
					         FluxPhi_Device, C_h, Gamma,
						 size, GridDimension[0], GridDimension[1], GridDimension[2], cTheta_Limiter);

    MHDComputedUy_CUDA3_kernel<<<dimGrid,dimBlock>>>(FluxD_Device, FluxS1_Device, FluxS2_Device, FluxS3_Device, FluxTau_Device,
  						     FluxBx_Device, FluxBy_Device, FluxBz_Device, FluxPhi_Device, 
					             dUD_Device, dUS1_Device, dUS2_Device, dUS3_Device,
				                     dUTau_Device, dUBx_Device, dUBy_Device, dUBz_Device, dUPhi_Device, 
					             dtdx, size, GridDimension[0], GridDimension[1], GridDimension[2]);
  }

  if (GridRank > 2) {
    MHDSweepZ_CUDA3_kernel<<<dimGrid,dimBlock>>>(Rho_Device, Vx_Device, Vy_Device, Vz_Device, Etot_Device, 
					         Eint_Device, Bx_Device, By_Device, Bz_Device, Phi_Device,
		 			         FluxD_Device, FluxS1_Device, FluxS2_Device, FluxS3_Device,
					         FluxTau_Device, FluxBx_Device, FluxBy_Device, FluxBz_Device,
					         FluxPhi_Device, C_h, Gamma, 
						 size, GridDimension[0], GridDimension[1], GridDimension[2], cTheta_Limiter);
  
    MHDComputedUz_CUDA3_kernel<<<dimGrid,dimBlock>>>(FluxD_Device, FluxS1_Device, FluxS2_Device, FluxS3_Device, FluxTau_Device,
  						     FluxBx_Device, FluxBy_Device, FluxBz_Device, FluxPhi_Device, 
					             dUD_Device, dUS1_Device, dUS2_Device, dUS3_Device,
				                     dUTau_Device, dUBx_Device, dUBy_Device, dUBz_Device, dUPhi_Device, 
					             dtdx, size, GridDimension[0], GridDimension[1], GridDimension[2]);
  }

  // update prim
  MHDUpdatePrim_CUDA3_kernel<<<dimGrid,dimBlock>>>(Rho_Device, Vx_Device, Vy_Device, Vz_Device, Etot_Device,
   					           Bx_Device, By_Device, Bz_Device, Phi_Device,	
					           dUD_Device, dUS1_Device, dUS2_Device, dUS3_Device,
				                   dUTau_Device, dUBx_Device, dUBy_Device, dUBz_Device, dUPhi_Device, 
						   dt, C_h, C_p, size);

  cudaEventRecord(stop, 0);
  cudaThreadSynchronize();
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsedTime, start, stop);
  if (PRINT_CUDA_TIMING) fprintf(stderr, "running kernel on GPU took:  %g \n" , elapsedTime/1e3);
  //  PerformanceTimers[36] += elapsedTime/1e3;

  // copy prim back to cpu
  cudaEventRecord(start, 0);

#ifdef ECUDADEBUG
  printf("SIZE: %i\n", size);
#endif

  CUDA_SAFE_CALL(cudaMemcpy(Prim[0], PrimDevice,        sizeof(float)*size, cudaMemcpyDeviceToHost));
  CUDA_SAFE_CALL(cudaMemcpy(Prim[1], PrimDevice+size,   sizeof(float)*size, cudaMemcpyDeviceToHost));
  CUDA_SAFE_CALL(cudaMemcpy(Prim[2], PrimDevice+2*size, sizeof(float)*size, cudaMemcpyDeviceToHost));
  CUDA_SAFE_CALL(cudaMemcpy(Prim[3], PrimDevice+3*size, sizeof(float)*size, cudaMemcpyDeviceToHost));
  CUDA_SAFE_CALL(cudaMemcpy(Prim[4], PrimDevice+4*size, sizeof(float)*size, cudaMemcpyDeviceToHost));
  CUDA_SAFE_CALL(cudaMemcpy(Prim[5], PrimDevice+5*size, sizeof(float)*size, cudaMemcpyDeviceToHost));
  CUDA_SAFE_CALL(cudaMemcpy(Prim[6], PrimDevice+6*size, sizeof(float)*size, cudaMemcpyDeviceToHost));
  CUDA_SAFE_CALL(cudaMemcpy(Prim[7], PrimDevice+7*size, sizeof(float)*size, cudaMemcpyDeviceToHost));
  CUDA_SAFE_CALL(cudaMemcpy(Prim[8], PrimDevice+8*size, sizeof(float)*size, cudaMemcpyDeviceToHost));


  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsedTime, start, stop);
  //  PerformanceTimers[37] += elapsedTime/1e3;    
  if (PRINT_CUDA_TIMING) fprintf(stderr, "copy data back to CPU took:  %g \n" , elapsedTime/1e3);

  cudaEventRecord(start, 0);

  CUDA_SAFE_CALL(cudaFree(PrimDevice));
  CUDA_SAFE_CALL(cudaFree(FluxDevice));
  CUDA_SAFE_CALL(cudaFree(dUDevice));

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsedTime, start, stop);
  if (PRINT_CUDA_TIMING) fprintf(stderr, "freeing memory on GPU took:  %g \n" , elapsedTime/1e3);
  //  PerformanceTimers[38] += elapsedTime/1e3;

#ifdef ECUDADEBUG
  printf("After update.\n");
  for (int j=30; j < 33; j++) 
    for (int i=0; i < 9; i++) printf("Prim[%i][%i] = %g \n", i, j, Prim[i][j]);
#endif

cudaEventDestroy( start ); 
cudaEventDestroy( stop ); 

  return SUCCESS;

}

__global__ void ComputeInternalEnergy_kernel(float *Rho, float *Vx, float *Vy, float *Vz, float *Etot,
                                             float *Eneint, float *Bx, float *By, float *Bz, int size)
{
  // get thread and block index
  const long tx = threadIdx.x;
  const long bx = blockIdx.x;
  const long by = blockIdx.y;

  int igrid = tx + bx*CUDA_BLOCK_SIZE + by*CUDA_BLOCK_SIZE*CUDA_GRID_SIZE;


  if (igrid >= size)
    return;

  // compute internal energy
  Eneint[igrid] = Etot[igrid] - 0.5*(Vx[igrid]*Vx[igrid] + Vy[igrid]*Vy[igrid] + Vz[igrid]*Vz[igrid]) -
                0.5*(Bx[igrid]*Bx[igrid] + By[igrid]*By[igrid] + Bz[igrid]*Bz[igrid])/Rho[igrid]; 

}
  
// kernel: load in data and call the computation routin
__global__ void MHDSweepX_CUDA3_kernel(float *Rho, float *Vx, float *Vy, float *Vz, float *Etot,
	    		               float *Eneint, float *Bx, float *By, float *Bz, float *Phi,
				       float *FluxD, float *FluxS1, float *FluxS2, float *FluxS3,
				       float *FluxTau, float *FluxBx, float *FluxBy, float *FluxBz,
				       float *FluxPhi, float C_h, float Gamma, int size, float cTheta_Limiter)
{
  // get thread and block index
  const long tx = threadIdx.x;
  const long bx = blockIdx.x;
  const long by = blockIdx.y;

  int igrid = tx + bx*CUDA_BLOCK_SIZE + by*CUDA_BLOCK_SIZE*CUDA_GRID_SIZE;

  // Every flux which is at cell interface is calculated 
  // by the cell to its right. This leads to the weird-looking asymmetry
  // in boundary condition loading below.
  __shared__ float PrimLine[NEQ_MHD*(CUDA_BLOCK_SIZE+3)]; // input primitive variable
  __shared__ float FluxLine[NEQ_MHD*CUDA_BLOCK_SIZE]; // output flux

  if (igrid >= 2 && igrid <= size - 2) {

  // load data from device to shared.
  int idx_prim1 = (tx+2)*NEQ_MHD;
  PrimLine[idx_prim1++] = Rho [igrid];
  PrimLine[idx_prim1++] = Eneint[igrid];
  PrimLine[idx_prim1++] = Vx  [igrid];
  PrimLine[idx_prim1++] = Vy  [igrid];
  PrimLine[idx_prim1++] = Vz  [igrid];
  PrimLine[idx_prim1++] = Bx  [igrid];
  PrimLine[idx_prim1++] = By  [igrid];
  PrimLine[idx_prim1++] = Bz  [igrid];
  PrimLine[idx_prim1  ] = Phi [igrid];

  // if the first, load in two more cells for boundary condition
  if (tx == 0 || igrid == 2) {
    for (int i = -2; i <=-1; i++) {
      const int idx_prim  = igrid + i;
      int idx_prim1 = (i+tx+2)*NEQ_MHD;
      PrimLine[idx_prim1++] = Rho [idx_prim];
      PrimLine[idx_prim1++] = Eneint[idx_prim];
      PrimLine[idx_prim1++] = Vx  [idx_prim];
      PrimLine[idx_prim1++] = Vy  [idx_prim];
      PrimLine[idx_prim1++] = Vz  [idx_prim];
      PrimLine[idx_prim1++] = Bx  [idx_prim];
      PrimLine[idx_prim1++] = By  [idx_prim];
      PrimLine[idx_prim1++] = Bz  [idx_prim];
      PrimLine[idx_prim1  ] = Phi [idx_prim];
    }
  }

  // if the last, load in one more cells for boundary condition
  if (tx == CUDA_BLOCK_SIZE - 1 || igrid == size - 2) {
    const int idx_prim  = igrid + 1;
    int idx_prim1 = (tx+3)*NEQ_MHD;
    PrimLine[idx_prim1++] = Rho [idx_prim];
    PrimLine[idx_prim1++] = Eneint[idx_prim];
    PrimLine[idx_prim1++] = Vx  [idx_prim];
    PrimLine[idx_prim1++] = Vy  [idx_prim];
    PrimLine[idx_prim1++] = Vz  [idx_prim];
    PrimLine[idx_prim1++] = Bx  [idx_prim];
    PrimLine[idx_prim1++] = By  [idx_prim];
    PrimLine[idx_prim1++] = Bz  [idx_prim];
    PrimLine[idx_prim1  ] = Phi [idx_prim];
  }
  }

  // synchronize to ensure all the data are loaded
  __syncthreads();

  if (igrid >= 2 && igrid <= size - 2) {
  // the main computation: calculating the flux at tx
    LLF_PLM_MHD_CUDA3(PrimLine, FluxLine, tx, C_h, Gamma, cTheta_Limiter);

  // copy 1D Flux back to Flux
  int idx_prim1 = tx*NEQ_MHD;
  FluxD  [igrid] = FluxLine[idx_prim1++];
  FluxS1 [igrid] = FluxLine[idx_prim1++];
  FluxS2 [igrid] = FluxLine[idx_prim1++];
  FluxS3 [igrid] = FluxLine[idx_prim1++];
  FluxTau[igrid] = FluxLine[idx_prim1++];
  FluxBx [igrid] = FluxLine[idx_prim1++];
  FluxBy [igrid] = FluxLine[idx_prim1++];
  FluxBz [igrid] = FluxLine[idx_prim1++];
  FluxPhi[igrid] = FluxLine[idx_prim1  ];
  }
}

// kernel: load in data and call the computation routin
__global__ void MHDSweepY_CUDA3_kernel(float *Rho, float *Vx, float *Vy, float *Vz, float *Etot,
	    		               float *Eneint, float *Bx, float *By, float *Bz, float *Phi,
				       float *FluxD, float *FluxS1, float *FluxS2, float *FluxS3,
				       float *FluxTau, float *FluxBx, float *FluxBy, float *FluxBz,
				       float *FluxPhi, float C_h, float Gamma, 
				       int size, int dim0, int dim1, int dim2, float cTheta_Limiter)
{
  // get thread and block index
  const long tx = threadIdx.x;
  const long bx = blockIdx.x;
  const long by = blockIdx.y;

  int igridy = tx + bx*CUDA_BLOCK_SIZE + by*CUDA_BLOCK_SIZE*CUDA_GRID_SIZE;

  int k = igridy/(dim0*dim1);
  int i = (igridy - k*dim0*dim1)/dim1;
  int j = igridy - k*dim0*dim1 - i*dim1;
  int igrid = i + (j + k*dim1) * dim0;

  // Every flux which is at cell interface is calculated 
  // by the cell to its right. This leads to the weird-looking asymmetry
  // in boundary condition loading below.
  __shared__ float PrimLine[NEQ_MHD*(CUDA_BLOCK_SIZE+3)]; // input primitive variable
  __shared__ float FluxLine[NEQ_MHD*CUDA_BLOCK_SIZE]; // output flux

  if (igridy >= 2 && igridy <= size - 2) {

  // load data from device to shared.
  int idx_prim1 = (tx+2)*NEQ_MHD;
  PrimLine[idx_prim1++] = Rho [igrid];
  PrimLine[idx_prim1++] = Eneint[igrid];
  PrimLine[idx_prim1++] = Vy  [igrid]; // vx = vy
  PrimLine[idx_prim1++] = Vz  [igrid]; // vy = vz
  PrimLine[idx_prim1++] = Vx  [igrid]; // vz = vx
  PrimLine[idx_prim1++] = By  [igrid]; // Bx = By
  PrimLine[idx_prim1++] = Bz  [igrid]; // By = Bz
  PrimLine[idx_prim1++] = Bx  [igrid]; // Bz = Bx
  PrimLine[idx_prim1  ] = Phi [igrid];

  // if the first, load in two more cells for boundary condition
  if (tx == 0 || igridy == 2) {
    for (int di = -2; di <=-1; di++) {
      int igridypi = igridy + di;
      int k = igridypi/(dim0*dim1);
      int i = (igridypi - k*dim0*dim1)/dim1;
      int j = igridypi - k*dim0*dim1 - i*dim1;
      int idx_prim = i + (j + k*dim1) * dim0;

      int idx_prim1 = (di+tx+2)*NEQ_MHD;
      PrimLine[idx_prim1++] = Rho [idx_prim];
      PrimLine[idx_prim1++] = Eneint[idx_prim];
      PrimLine[idx_prim1++] = Vy  [idx_prim];
      PrimLine[idx_prim1++] = Vz  [idx_prim];
      PrimLine[idx_prim1++] = Vx  [idx_prim];
      PrimLine[idx_prim1++] = By  [idx_prim];
      PrimLine[idx_prim1++] = Bz  [idx_prim];
      PrimLine[idx_prim1++] = Bx  [idx_prim];
      PrimLine[idx_prim1  ] = Phi [idx_prim];
    }
  }

  // if the last, load in one more cells for boundary condition
  if (tx == CUDA_BLOCK_SIZE - 1 || igridy == size - 2) {
    int igridyp1 = igridy + 1;
    int k = igridyp1/(dim0*dim1);
    int i = (igridyp1 - k*dim0*dim1)/dim1;
    int j = igridyp1 - k*dim0*dim1 - i*dim1;
    int idx_prim = i + (j + k*dim1) * dim0;

    int idx_prim1 = (tx+3)*NEQ_MHD;
    PrimLine[idx_prim1++] = Rho [idx_prim];
    PrimLine[idx_prim1++] = Eneint[idx_prim];
    PrimLine[idx_prim1++] = Vy  [idx_prim];
    PrimLine[idx_prim1++] = Vz  [idx_prim];
    PrimLine[idx_prim1++] = Vx  [idx_prim];
    PrimLine[idx_prim1++] = By  [idx_prim];
    PrimLine[idx_prim1++] = Bz  [idx_prim];
    PrimLine[idx_prim1++] = Bx  [idx_prim];
    PrimLine[idx_prim1  ] = Phi [idx_prim];
  }
  }

  // synchronize to ensure all the data are loaded
  __syncthreads();

  if (igridy >= 2 && igridy <= size - 2) {
  // the main computation: calculating the flux at tx
    LLF_PLM_MHD_CUDA3(PrimLine, FluxLine, tx, C_h, Gamma, cTheta_Limiter);

  // copy 1D Flux back to Flux
  int idx_prim1 = tx*NEQ_MHD;
  FluxD  [igrid] = FluxLine[idx_prim1++];
  FluxS2 [igrid] = FluxLine[idx_prim1++]; // vy = vx
  FluxS3 [igrid] = FluxLine[idx_prim1++]; // vz = vy
  FluxS1 [igrid] = FluxLine[idx_prim1++]; // vx = vz
  FluxTau[igrid] = FluxLine[idx_prim1++];
  FluxBy [igrid] = FluxLine[idx_prim1++]; // By = Bx
  FluxBz [igrid] = FluxLine[idx_prim1++]; // Bz = By
  FluxBx [igrid] = FluxLine[idx_prim1++]; // Bx = Bz
  FluxPhi[igrid] = FluxLine[idx_prim1  ];
  }

}

// kernel: load in data and call the computation routin
__global__ void MHDSweepZ_CUDA3_kernel(float *Rho, float *Vx, float *Vy, float *Vz, float *Etot,
	    		               float *Eneint, float *Bx, float *By, float *Bz, float *Phi,
				       float *FluxD, float *FluxS1, float *FluxS2, float *FluxS3,
				       float *FluxTau, float *FluxBx, float *FluxBy, float *FluxBz,
				       float *FluxPhi, float C_h, float Gamma, 
				       int size, int dim0, int dim1, int dim2, float cTheta_Limiter)
{
  // get thread and block index
  const long tx = threadIdx.x;
  const long bx = blockIdx.x;
  const long by = blockIdx.y;

  int igridz = tx + bx*CUDA_BLOCK_SIZE + by*CUDA_BLOCK_SIZE*CUDA_GRID_SIZE;

  int j = igridz / (dim0*dim2);
  int i = (igridz - j*dim0*dim2) / dim2;
  int k = igridz - j*dim0*dim2 - i*dim2;
  int igrid = i + (j + k*dim1) * dim0;

  // Every flux which is at cell interface is calculated 
  // by the cell to its right. This leads to the weird-looking asymmetry
  // in boundary condition loading below.
  __shared__ float PrimLine[NEQ_MHD*(CUDA_BLOCK_SIZE+3)]; // input primitive variable
  __shared__ float FluxLine[NEQ_MHD*CUDA_BLOCK_SIZE]; // output flux

  if (igridz >= 2 && igridz <= size - 2) {

  // load data from device to shared.
  int idx_prim1 = (tx+2)*NEQ_MHD;
  PrimLine[idx_prim1++] = Rho [igrid];
  PrimLine[idx_prim1++] = Eneint[igrid];
  PrimLine[idx_prim1++] = Vz  [igrid]; // vx = vz
  PrimLine[idx_prim1++] = Vx  [igrid]; // vy = vx
  PrimLine[idx_prim1++] = Vy  [igrid]; // vz = vy
  PrimLine[idx_prim1++] = Bz  [igrid]; // Bx = Bz
  PrimLine[idx_prim1++] = Bx  [igrid]; // By = Bx
  PrimLine[idx_prim1++] = By  [igrid]; // Bz = By
  PrimLine[idx_prim1  ] = Phi [igrid];
  
  // if the first, load in two more cells for boundary condition
  if (tx == 0 || igridz == 2) {
    for (int di = -2; di <=-1; di++) {

      int igridzpi = igridz + di;
      int j = igridzpi / (dim0*dim2);
      int i = (igridzpi - j*dim0*dim2) / dim2;
      int k = igridzpi - j*dim0*dim2 - i*dim2;
      int idx_prim = i + (j + k*dim1) * dim0;

      int idx_prim1 = (di+tx+2)*NEQ_MHD;
      PrimLine[idx_prim1++] = Rho [idx_prim];
      PrimLine[idx_prim1++] = Eneint[idx_prim];
      PrimLine[idx_prim1++] = Vz  [idx_prim];
      PrimLine[idx_prim1++] = Vx  [idx_prim];
      PrimLine[idx_prim1++] = Vy  [idx_prim];
      PrimLine[idx_prim1++] = Bz  [idx_prim];
      PrimLine[idx_prim1++] = Bx  [idx_prim];
      PrimLine[idx_prim1++] = By  [idx_prim];
      PrimLine[idx_prim1  ] = Phi [idx_prim];
    }
  }

  // if the last, load in one more cells for boundary condition
  if (tx == CUDA_BLOCK_SIZE - 1 || igridz == size - 2) {
    int igridzp1 = igridz + 1;
    int j = igridzp1 / (dim0*dim2);
    int i = (igridzp1 - j*dim0*dim2) / dim2;
    int k = igridzp1 - j*dim0*dim2 - i*dim2;
    int idx_prim = i + (j + k*dim1) * dim0;

    int idx_prim1 = (tx+3)*NEQ_MHD;
    PrimLine[idx_prim1++] = Rho [idx_prim];
    PrimLine[idx_prim1++] = Eneint[idx_prim];
    PrimLine[idx_prim1++] = Vz  [idx_prim];
    PrimLine[idx_prim1++] = Vx  [idx_prim];
    PrimLine[idx_prim1++] = Vy  [idx_prim];
    PrimLine[idx_prim1++] = Bz  [idx_prim];
    PrimLine[idx_prim1++] = Bx  [idx_prim];
    PrimLine[idx_prim1++] = By  [idx_prim];
    PrimLine[idx_prim1  ] = Phi [idx_prim];
  }
  } 

  // synchronize to ensure all the data are loaded
  __syncthreads();

  if (igridz >= 2 && igridz <= size - 2) {
  // the main computation: calculating the flux at tx
    LLF_PLM_MHD_CUDA3(PrimLine, FluxLine, tx, C_h, Gamma, cTheta_Limiter);

  // copy 1D Flux back to Flux
  int idx_prim1 = tx*NEQ_MHD;
  FluxD  [igrid] = FluxLine[idx_prim1++];
  FluxS3 [igrid] = FluxLine[idx_prim1++]; // vz = vx
  FluxS1 [igrid] = FluxLine[idx_prim1++]; // vx = vy
  FluxS2 [igrid] = FluxLine[idx_prim1++]; // vy = vz
  FluxTau[igrid] = FluxLine[idx_prim1++];
  FluxBz [igrid] = FluxLine[idx_prim1++]; // Bz = Bx
  FluxBx [igrid] = FluxLine[idx_prim1++]; // Bx = By
  FluxBy [igrid] = FluxLine[idx_prim1++]; // By = Bz
  FluxPhi[igrid] = FluxLine[idx_prim1  ];
  }
}

__global__ void MHDComputedUx_CUDA3_kernel(float *FluxD, float *FluxS1, float *FluxS2, float *FluxS3,
				           float *FluxTau, float *FluxBx, float *FluxBy, float *FluxBz, float *FluxPhi, 
					   float *dUD, float *dUS1, float *dUS2, float *dUS3,
				           float *dUTau, float *dUBx, float *dUBy, float *dUBz, float *dUPhi, 
					   float dtdx, int size)
{
  // get thread and block index
  const long tx = threadIdx.x;
  const long bx = blockIdx.x;
  const long by = blockIdx.y;

  int igrid = tx + bx*CUDA_BLOCK_SIZE + by*CUDA_BLOCK_SIZE*CUDA_GRID_SIZE;


  if (igrid < 2 || igrid > size - 3)
    return;

  int igridp1 = igrid + 1;
  dUD  [igrid] = (FluxD  [igrid] - FluxD  [igridp1])*dtdx;
  dUS1 [igrid] = (FluxS1 [igrid] - FluxS1 [igridp1])*dtdx;
  dUS2 [igrid] = (FluxS2 [igrid] - FluxS2 [igridp1])*dtdx;
  dUS3 [igrid] = (FluxS3 [igrid] - FluxS3 [igridp1])*dtdx;
  dUTau[igrid] = (FluxTau[igrid] - FluxTau[igridp1])*dtdx;
  dUBx [igrid] = (FluxBx [igrid] - FluxBx [igridp1])*dtdx;
  dUBy [igrid] = (FluxBy [igrid] - FluxBy [igridp1])*dtdx;
  dUBz [igrid] = (FluxBz [igrid] - FluxBz [igridp1])*dtdx;
  dUPhi[igrid] = (FluxPhi[igrid] - FluxPhi[igridp1])*dtdx;

}

__global__ void MHDComputedUy_CUDA3_kernel(float *FluxD, float *FluxS1, float *FluxS2, float *FluxS3,
				           float *FluxTau, float *FluxBx, float *FluxBy, float *FluxBz, float *FluxPhi, 
					   float *dUD, float *dUS1, float *dUS2, float *dUS3,
				           float *dUTau, float *dUBx, float *dUBy, float *dUBz, float *dUPhi, 
					   float dtdx, int size, int dim0, int dim1, int dim2)
{
  // get thread and block index
  const long tx = threadIdx.x;
  const long bx = blockIdx.x;
  const long by = blockIdx.y;

  int igridy = tx + bx*CUDA_BLOCK_SIZE + by*CUDA_BLOCK_SIZE*CUDA_GRID_SIZE;

  if (igridy < 2 || igridy > size - 3)
    return;

  int k = igridy/(dim0*dim1);
  int i = (igridy - k*dim0*dim1)/dim1;
  int j = igridy - k*dim0*dim1 - i*dim1;
  int igrid = i + (j + k*dim1) * dim0;

  int igridyp1 = igridy + 1;
  k = igridyp1/(dim0*dim1);
  i = (igridyp1 - k*dim0*dim1)/dim1;
  j = igridyp1 - k*dim0*dim1 - i*dim1;
  int igridp1 = i + (j + k*dim1) * dim0;


  dUD  [igrid] += (FluxD  [igrid] - FluxD  [igridp1])*dtdx;
  dUS1 [igrid] += (FluxS1 [igrid] - FluxS1 [igridp1])*dtdx;
  dUS2 [igrid] += (FluxS2 [igrid] - FluxS2 [igridp1])*dtdx;
  dUS3 [igrid] += (FluxS3 [igrid] - FluxS3 [igridp1])*dtdx;
  dUTau[igrid] += (FluxTau[igrid] - FluxTau[igridp1])*dtdx;
  dUBx [igrid] += (FluxBx [igrid] - FluxBx [igridp1])*dtdx;
  dUBy [igrid] += (FluxBy [igrid] - FluxBy [igridp1])*dtdx;
  dUBz [igrid] += (FluxBz [igrid] - FluxBz [igridp1])*dtdx;
  dUPhi[igrid] += (FluxPhi[igrid] - FluxPhi[igridp1])*dtdx;

}

__global__ void MHDComputedUz_CUDA3_kernel(float *FluxD, float *FluxS1, float *FluxS2, float *FluxS3,
				           float *FluxTau, float *FluxBx, float *FluxBy, float *FluxBz, float *FluxPhi, 
					   float *dUD, float *dUS1, float *dUS2, float *dUS3,
				           float *dUTau, float *dUBx, float *dUBy, float *dUBz, float *dUPhi, 
					   float dtdx, int size, int dim0, int dim1, int dim2)
{
  // get thread and block index
  const long tx = threadIdx.x;
  const long bx = blockIdx.x;
  const long by = blockIdx.y;

  int igridz = tx + bx*CUDA_BLOCK_SIZE + by*CUDA_BLOCK_SIZE*CUDA_GRID_SIZE;


  if (igridz < 2 || igridz > size - 3)
    return;

  int j = igridz / (dim0*dim2);
  int i = (igridz - j*dim0*dim2) / dim2;
  int k = igridz - j*dim0*dim2 - i*dim2;
  int igrid = i + (j + k*dim1) * dim0;

  int igridzp1 = igridz + 1;
  j = igridzp1 / (dim0*dim2);
  i = (igridzp1 - j*dim0*dim2) / dim2;
  k = igridzp1 - j*dim0*dim2 - i*dim2;
  int igridp1 = i + (j + k*dim1) * dim0;

  dUD  [igrid] += (FluxD  [igrid] - FluxD  [igridp1])*dtdx;
  dUS1 [igrid] += (FluxS1 [igrid] - FluxS1 [igridp1])*dtdx;
  dUS2 [igrid] += (FluxS2 [igrid] - FluxS2 [igridp1])*dtdx;
  dUS3 [igrid] += (FluxS3 [igrid] - FluxS3 [igridp1])*dtdx;
  dUTau[igrid] += (FluxTau[igrid] - FluxTau[igridp1])*dtdx;
  dUBx [igrid] += (FluxBx [igrid] - FluxBx [igridp1])*dtdx;
  dUBy [igrid] += (FluxBy [igrid] - FluxBy [igridp1])*dtdx;
  dUBz [igrid] += (FluxBz [igrid] - FluxBz [igridp1])*dtdx;
  dUPhi[igrid] += (FluxPhi[igrid] - FluxPhi[igridp1])*dtdx;

}

__global__ void MHDUpdatePrim_CUDA3_kernel(float *Rho, float *Vx, float *Vy, float *Vz,
	  			           float *Etot, float *Bx, float *By, float *Bz, float *Phi, 
					   float *dUD, float *dUS1, float *dUS2, float *dUS3,
				           float *dUTau, float *dUBx, float *dUBy, float *dUBz, float *dUPhi, 
					   float dt, float C_h, float C_p, int size)
{
  // get thread and block index
  const long tx = threadIdx.x;
  const long bx = blockIdx.x;
  const long by = blockIdx.y;

  int igrid = tx + bx*CUDA_BLOCK_SIZE + by*CUDA_BLOCK_SIZE*CUDA_GRID_SIZE;

  if (igrid < 2 || igrid > size - 3)
    return;

  float D, S1, S2, S3, Tau;
  D   = Rho[igrid];
  S1  = D*Vx[igrid];
  S2  = D*Vy[igrid];
  S3  = D*Vz[igrid];
  Tau = D*Etot[igrid];

  D   += dUD[igrid];
  S1  += dUS1[igrid];
  S2  += dUS2[igrid];
  S3  += dUS3[igrid];
  Tau += dUTau[igrid];

  Rho[igrid] = D;
  Vx[igrid] = S1/D;
  Vy[igrid] = S2/D;
  Vz[igrid] = S3/D;
  Etot[igrid] = Tau/D;
  
  Bx[igrid] += dUBx[igrid];
  By[igrid] += dUBy[igrid];
  Bz[igrid] += dUBz[igrid];
  Phi[igrid] += dUPhi[igrid];  
  Phi[igrid] *= expf(-dt*(C_h/C_p)*(C_h/C_p));
}

// the main computation routine: compute Flux at the left cell interface given Prim at the center
__device__ void LLF_PLM_MHD_CUDA3(float *Prim, float *Flux, 
				  const int &tx, const float &C_h, const float &cGamma, const float &cTheta_Limiter)
{
  // those crazily many fields are defined in stead of using arrays
  // to avoid letting the compiler putting arrays in local memory
  float Priml_rho, Priml_eint, Priml_vx, Priml_vy, Priml_vz, Priml_Bx, Priml_By, Priml_Bz, Priml_Phi;
  float Primr_rho, Primr_eint, Primr_vx, Primr_vy, Primr_vz, Primr_Bx, Primr_By, Primr_Bz, Primr_Phi;
  float Ul_D, Ul_S1, Ul_S2, Ul_S3, Ul_Etot, Ul_Bx, Ul_By, Ul_Bz, Ul_Phi;
  float Ur_D, Ur_S1, Ur_S2, Ur_S3, Ur_Etot, Ur_Bx, Ur_By, Ur_Bz, Ur_Phi;
  float Fl_D, Fl_S1, Fl_S2, Fl_S3, Fl_Etot, Fl_Bx, Fl_By, Fl_Bz;
  float Fr_D, Fr_S1, Fr_S2, Fr_S3, Fr_Etot, Fr_Bx, Fr_By, Fr_Bz;
  float cs, v2, p, lm_l, lp_l, lm_r, lp_r, B2; 

  // 1. do PLM reconstruction for all primitive fields

  plm_point(Prim[tx*NEQ_MHD+irho], Prim[(tx+1)*NEQ_MHD+irho], Prim[(tx+2)*NEQ_MHD+irho], Priml_rho, cTheta_Limiter);
  plm_point(Prim[(tx+3)*NEQ_MHD+irho], Prim[(tx+2)*NEQ_MHD+irho], Prim[(tx+1)*NEQ_MHD+irho], Primr_rho, cTheta_Limiter);

  plm_point(Prim[tx*NEQ_MHD+ieint], Prim[(tx+1)*NEQ_MHD+ieint], Prim[(tx+2)*NEQ_MHD+ieint], Priml_eint, cTheta_Limiter);
  plm_point(Prim[(tx+3)*NEQ_MHD+ieint], Prim[(tx+2)*NEQ_MHD+ieint], Prim[(tx+1)*NEQ_MHD+ieint], Primr_eint, cTheta_Limiter);

  plm_point(Prim[tx*NEQ_MHD+ivx], Prim[(tx+1)*NEQ_MHD+ivx], Prim[(tx+2)*NEQ_MHD+ivx], Priml_vx, cTheta_Limiter);
  plm_point(Prim[(tx+3)*NEQ_MHD+ivx], Prim[(tx+2)*NEQ_MHD+ivx], Prim[(tx+1)*NEQ_MHD+ivx], Primr_vx, cTheta_Limiter);

  plm_point(Prim[tx*NEQ_MHD+ivy], Prim[(tx+1)*NEQ_MHD+ivy], Prim[(tx+2)*NEQ_MHD+ivy], Priml_vy, cTheta_Limiter);
  plm_point(Prim[(tx+3)*NEQ_MHD+ivy], Prim[(tx+2)*NEQ_MHD+ivy], Prim[(tx+1)*NEQ_MHD+ivy], Primr_vy, cTheta_Limiter);

  plm_point(Prim[tx*NEQ_MHD+ivz], Prim[(tx+1)*NEQ_MHD+ivz], Prim[(tx+2)*NEQ_MHD+ivz], Priml_vz, cTheta_Limiter);
  plm_point(Prim[(tx+3)*NEQ_MHD+ivz], Prim[(tx+2)*NEQ_MHD+ivz], Prim[(tx+1)*NEQ_MHD+ivz], Primr_vz, cTheta_Limiter);

  plm_point(Prim[tx*NEQ_MHD+iBx], Prim[(tx+1)*NEQ_MHD+iBx], Prim[(tx+2)*NEQ_MHD+iBx], Priml_Bx, cTheta_Limiter);
  plm_point(Prim[(tx+3)*NEQ_MHD+iBx], Prim[(tx+2)*NEQ_MHD+iBx], Prim[(tx+1)*NEQ_MHD+iBx], Primr_Bx, cTheta_Limiter);

  plm_point(Prim[tx*NEQ_MHD+iBy], Prim[(tx+1)*NEQ_MHD+iBy], Prim[(tx+2)*NEQ_MHD+iBy], Priml_By, cTheta_Limiter);
  plm_point(Prim[(tx+3)*NEQ_MHD+iBy], Prim[(tx+2)*NEQ_MHD+iBy], Prim[(tx+1)*NEQ_MHD+iBy], Primr_By, cTheta_Limiter);

  plm_point(Prim[tx*NEQ_MHD+iBz], Prim[(tx+1)*NEQ_MHD+iBz], Prim[(tx+2)*NEQ_MHD+iBz], Priml_Bz, cTheta_Limiter);
  plm_point(Prim[(tx+3)*NEQ_MHD+iBz], Prim[(tx+2)*NEQ_MHD+iBz], Prim[(tx+1)*NEQ_MHD+iBz], Primr_Bz, cTheta_Limiter);

  plm_point(Prim[tx*NEQ_MHD+iPhi], Prim[(tx+1)*NEQ_MHD+iPhi], Prim[(tx+2)*NEQ_MHD+iPhi], Priml_Phi, cTheta_Limiter);
  plm_point(Prim[(tx+3)*NEQ_MHD+iPhi], Prim[(tx+2)*NEQ_MHD+iPhi], Prim[(tx+1)*NEQ_MHD+iPhi], Primr_Phi, cTheta_Limiter);

  // 2. use LLF Riemann solver to compute flux  

  // 2.1, compute Fl and Ul
  B2 = pow(Priml_Bx,2) + pow(Priml_By,2) + pow(Priml_Bz,2);
  v2 = pow(Priml_vx,2) + pow(Priml_vy,2) + pow(Priml_vz,2);

  EOS(p, Priml_rho, Priml_eint, cs, cGamma, EOSType, 2);    

  Ul_D    = Priml_rho;
  Ul_S1   = Priml_rho * Priml_vx;
  Ul_S2   = Priml_rho * Priml_vy;
  Ul_S3   = Priml_rho * Priml_vz;
  Ul_Etot = Priml_rho * (Priml_eint + 0.5*v2 + 0.5*B2/Priml_rho);
  Ul_Bx   = Priml_Bx;
  Ul_By   = Priml_By;
  Ul_Bz   = Priml_Bz;
  Ul_Phi  = Priml_Phi;

  Fl_D    = Priml_rho * Priml_vx;
  Fl_S1   = Ul_S1 * Priml_vx + p + 0.5*B2 - pow(Priml_Bx,2);
  Fl_S2   = Ul_S2 * Priml_vx - Priml_Bx*Priml_By;
  Fl_S3   = Ul_S3 * Priml_vx - Priml_Bx*Priml_Bz;
  Fl_Etot = Priml_rho*(0.5*v2 + Priml_eint + p/Priml_rho)*Priml_vx + B2*Priml_vx - 
            Priml_Bx*(Priml_Bx*Priml_vx + Priml_By*Priml_vy + Priml_Bz*Priml_vz);
  Fl_Bx  = 0.0;
  Fl_By  = Priml_vx*Priml_By - Priml_vy*Priml_Bx;
  Fl_Bz  = -Priml_vz*Priml_Bx + Priml_vx*Priml_Bz;

  // largest and smallest eigenvectors, reuse the space of cs
  cs = sqrtf(0.5 * (cs*cs+B2/Priml_rho + 
            sqrtf(fabs(pow(cs*cs+B2/Priml_rho,2) - 
	         4.0*cs*cs*pow(Priml_Bx,2)/Priml_rho))));

  lp_l = Priml_vx + cs;
  lm_l = Priml_vx - cs;

  // 2.2 Fr and Ur
  B2 = pow(Primr_Bx,2) + pow(Primr_By,2) + pow(Primr_Bz,2);
  v2 = pow(Primr_vx,2) + pow(Primr_vy,2) + pow(Primr_vz,2);

  EOS(p, Primr_rho, Primr_eint, cs, cGamma, EOSType, 2);

  Ur_D    = Primr_rho;
  Ur_S1   = Primr_rho * Primr_vx;
  Ur_S2   = Primr_rho * Primr_vy;
  Ur_S3   = Primr_rho * Primr_vz;
  Ur_Etot = Primr_rho * (Primr_eint + 0.5*v2 + 0.5*B2/Primr_rho);
  Ur_Bx   = Primr_Bx;
  Ur_By   = Primr_By;
  Ur_Bz   = Primr_Bz;
  Ur_Phi  = Primr_Phi;

  Fr_D    = Primr_rho * Primr_vx;
  Fr_S1   = Ur_S1 * Primr_vx + p + 0.5*B2 - pow(Primr_Bx,2);
  Fr_S2   = Ur_S2 * Primr_vx - Primr_Bx*Primr_By;
  Fr_S3   = Ur_S3 * Primr_vx - Primr_Bx*Primr_Bz;
  Fr_Etot = Primr_rho*(0.5*v2 + Primr_eint + p/Primr_rho)*Primr_vx + B2*Primr_vx - 
            Primr_Bx*(Primr_Bx*Primr_vx + Primr_By*Primr_vy + Primr_Bz*Primr_vz);
  Fr_Bx   = 0.0;
  Fr_By   = Primr_vx*Primr_By - Primr_vy*Primr_Bx;
  Fr_Bz   = -Primr_vz*Primr_Bx + Primr_vx*Primr_Bz;

  // largest and smallest eigenvectors, reuse the space of cs
  cs = sqrtf(0.5 * (cs*cs+B2/Primr_rho + 
            sqrtf(fabs(pow(cs*cs+B2/Primr_rho,2) - 
	          4.0*cs*cs*pow(Primr_Bx,2)/Primr_rho))));

  lp_r = Primr_vx + cs;
  lm_r = Primr_vx - cs;

  // 2.3. compute the maximum local wave speed
  lp_l = Max_mhd(0, lp_l, lp_r);
  lm_l = Max_mhd(0, -lm_l, -lm_r);
  lp_l = (lp_l > lm_l) ? lp_l : lm_l;

  // 2.4. compute the flux
  Flux[tx*NEQ_MHD+iD   ] = 0.5*(Fl_D+Fr_D-lp_l*(Ur_D-Ul_D));
  Flux[tx*NEQ_MHD+iS1  ] = 0.5*(Fl_S1+Fr_S1-lp_l*(Ur_S1-Ul_S1));
  Flux[tx*NEQ_MHD+iS2  ] = 0.5*(Fl_S2+Fr_S2-lp_l*(Ur_S2-Ul_S2));
  Flux[tx*NEQ_MHD+iS3  ] = 0.5*(Fl_S3+Fr_S3-lp_l*(Ur_S3-Ul_S3));
  Flux[tx*NEQ_MHD+iEtot] = 0.5*(Fl_Etot+Fr_Etot-lp_l*(Ur_Etot-Ul_Etot));
  Flux[tx*NEQ_MHD+iBx  ] = 0.5*(Fl_Bx+Fr_Bx-lp_l*(Ur_Bx-Ul_Bx));
  Flux[tx*NEQ_MHD+iBy  ] = 0.5*(Fl_By+Fr_By-lp_l*(Ur_By-Ul_By));
  Flux[tx*NEQ_MHD+iBz  ] = 0.5*(Fl_Bz+Fr_Bz-lp_l*(Ur_Bz-Ul_Bz));

  Flux[tx*NEQ_MHD+iBx  ] += Ul_Phi + 0.5*(Ur_Phi-Ul_Phi) - 0.5*C_h*(Ur_Bx-Ul_Bx);
  Flux[tx*NEQ_MHD+iPhi ] =  Ul_Bx  + 0.5*(Ur_Bx-Ul_Bx) - 0.5/C_h*(Ur_Phi-Ul_Phi);
  Flux[tx*NEQ_MHD+iPhi ] *= (C_h*C_h);

}

// minmod limiter function
__device__ float minmod_mhd(const float &a, const float &b, const float &c)
{
  return 0.25*(sign(a)+sign(b))*ABS((sign(a)+sign(c)))*Min_mhd(ABS(a), ABS(b), ABS(c));
}

// do PLM reconstruction for a point
__device__ void plm_point(const float &vm1, const float &v, const float &vp1, float &vl_plm, const float &cTheta_Limiter)
{
  float dv_l = (v-vm1) * cTheta_Limiter;
  float dv_r = (vp1-v) * cTheta_Limiter;
  float dv_m = 0.5*(vp1-vm1);
  
  float dv = minmod_mhd(dv_l, dv_r, dv_m);

  vl_plm = v + 0.5*dv;
}

// return the maximum of a, b, c
__device__ float Max_mhd(const float &a, const float &b, const float &c)  
{
  if (a > b) {
    if (a > c)
      return a;
    else 
      return c;
  } else {
    if (b > c)
      return b;
    else
      return c;
  }
}

// return the minimum of a, b, c
__device__ float Min_mhd(const float &a, const float &b, const float &c)
{
  if (a<b) {
    if (c<a)
      return c;
    else 
      return a;
  } else {
    if (c<b)
      return c;
    else 
      return b;
  }
}

/***********************************************************
/
/   Function EOS: compute equation of state
/
/   Input: 
/     eostype: 
/       0: ideal gas
/     mode:  
/       1: given p and rho, calculate others.
/       2: given rho and e, calculate others.
/
************************************************************/

__device__ void EOS(float &p, float &rho, float &e, float &cs, const float &cGamma, const int &eostype, const int &mode)
{

  if (eostype == 0) {
    
    if (mode == 1) {
      e = p / rho / (cGamma - 1);      
    } else if (mode == 2) {
      p = (cGamma - 1) * rho * e;
    }

    cs = sqrt(cGamma*p/rho);

  }

  if (eostype == 3) { // straight isothermal
    cs = EOSSoundSpeed;
    p = rho*cs*cs;
    e = p / ((cGamma-1.0)*rho);
  }

}
