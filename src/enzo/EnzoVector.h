/*****************************************************************************
 *                                                                           *
 * Copyright 2006 Daniel R. Reynolds                                         *
 * Copyright 2006 Laboratory for Computational Astrophysics                  *
 * Copyright 2006 Regents of the University of California                    *
 *                                                                           *
 * This software is released under the terms of the "Enzo Public License"    *
 * in the accompanying LICENSE file.                                         *
 *                                                                           *
 *****************************************************************************/
/***********************************************************************
/
/  Enzo Vector Class, used within linear and nonlinear solvers
/
/  written by: Daniel Reynolds
/  date:       May, 2006
/  modified1:  
/
/  PURPOSE: This class defines a general Enzo vector (one level only).  
/  In addition to containing general 3D local domain data, it contains 
/  a data pointer array.  For each species included in the vector, 
/  this data pointer array contains a pointer to a 3D dataset, 
/  containing the actual vector entries.  Each 3D data set is defined 
/  as generally as possible, such that in each dimension the local 
/  array length is set to be, e.g. Nx0 + Ng0l + Ng0r, where Nx0 
/  represents the number of active local entries, Nx0l represents the 
/  number of ghost cells at the left x0 face, and Nx0r represents the 
/  number of ghost cells at the right x0 face.
/
************************************************************************/

#ifndef ENZO_VECTOR_DEFINED__
#define ENZO_VECTOR_DEFINED__


class EnzoVector {
  
 private:

  // local lengths
  int Nx0;        // x0 local active length
  int Nx1;        // x1 local active length
  int Nx2;        // x2 local active length
  int Ng0l;       // x0 left face ghost cells
  int Ng0r;       // x0 right face ghost cells
  int Ng1l;       // x1 left face ghost cells
  int Ng1r;       // x1 right face ghost cells
  int Ng2l;       // x2 left face ghost cells
  int Ng2r;       // x2 right face ghost cells
  int Nspecies;   // number of species in this vector
  int Nglobal;    // total active cells (all procs)
  float **data;   // array of data arrays
  bool owndata;   // flag denoting whether vector or calling routine owns data
  
  // parallelism information (see EnzoVector_Exchange* for more info.)
  int Nbors[3][2];  // neighbor procs in each direction 
                    // (set to MPI_PROC_NULL (or -3) for none)
  int NborGhosts[3][2];  // ghost zones that each neighbor expects
                    // size of each send/recv buffer
  int x0Lsendsize, x0Rsendsize, x0Lrecvsize, x0Rrecvsize;
  int x1Lsendsize, x1Rsendsize, x1Lrecvsize, x1Rrecvsize;
  int x2Lsendsize, x2Rsendsize, x2Lrecvsize, x2Rrecvsize;
  int x0Lsendsize_comp, x0Rsendsize_comp, x0Lrecvsize_comp, x0Rrecvsize_comp;
  int x1Lsendsize_comp, x1Rsendsize_comp, x1Lrecvsize_comp, x1Rrecvsize_comp;
  int x2Lsendsize_comp, x2Rsendsize_comp, x2Lrecvsize_comp, x2Rrecvsize_comp;
                    // send/recv buffers
  float *SendBufX0l, *RecvBufX0l, *SendBufX0r, *RecvBufX0r;
  float *SendBufX1l, *RecvBufX1l, *SendBufX1r, *RecvBufX1r;
  float *SendBufX2l, *RecvBufX2l, *SendBufX2r, *RecvBufX2r;
  float *SendBufX0l_comp, *RecvBufX0l_comp, *SendBufX0r_comp, *RecvBufX0r_comp;
  float *SendBufX1l_comp, *RecvBufX1l_comp, *SendBufX1r_comp, *RecvBufX1r_comp;
  float *SendBufX2l_comp, *RecvBufX2l_comp, *SendBufX2r_comp, *RecvBufX2r_comp;
                    // MPI_Request handles (must declare as void* here since 
                    // Enzo re-maps int types, making include files messy)
  void *id_recv_x0l, *id_send_x0l, *id_recv_x0l_comp, *id_send_x0l_comp;
  void *id_recv_x0r, *id_send_x0r, *id_recv_x0r_comp, *id_send_x0r_comp;
  void *id_recv_x1l, *id_send_x1l, *id_recv_x1l_comp, *id_send_x1l_comp;
  void *id_recv_x1r, *id_send_x1r, *id_recv_x1r_comp, *id_send_x1r_comp;
  void *id_recv_x2l, *id_send_x2l, *id_recv_x2l_comp, *id_send_x2l_comp;
  void *id_recv_x2r, *id_send_x2r, *id_recv_x2r_comp, *id_send_x2r_comp;


 public:

  // Constructor
  EnzoVector(int Nx0, int Nx1, int Nx2, 
	     int Ng0l, int Ng0r, int Ng1l, 
	     int Ng1r, int Ng2l, int Ng2r, int Ns,
	     int NBx0L, int NBx0R, int NBx1L, 
	     int NBx1R, int NBx2L, int NBx2R);
  
  // Constructor (overloaded -- sets data to pre-allocated space)
  EnzoVector(int Nx0, int Nx1, int Nx2, 
	     int Ng0l, int Ng0r, int Ng1l, 
	     int Ng1r, int Ng2l, int Ng2r, int Ns,
	     int NBx0L, int NBx0R, int NBx1L, 
	     int NBx1R, int NBx2L, int NBx2R,
	     float **userdata);
  
  // Constructor (overloaded -- sets data arrays to NULL)
  EnzoVector(int Nx0, int Nx1, int Nx2, 
	     int Ng0l, int Ng0r, int Ng1l, 
	     int Ng1r, int Ng2l, int Ng2r, int Ns,
	     int NBx0L, int NBx0R, int NBx1L, 
	     int NBx1R, int NBx2L, int NBx2R, int Empty);
  
  // Destructor
  ~EnzoVector();
  

  // Vector operations

  //   Clone a vector, creating new data arrays
  EnzoVector* clone() const;

  //   Clone a vector, using provided data arrays
  EnzoVector* clone(float **userdata) const;

  //   Writes a given species to file (no ghosts)
  int write(char *outfile, int species) const;

  //   Writes a given species to file (with ghosts)
  int writeall(char *outfile, int species) const;

  //   Communicates ghost cells with neighbors
  int exchange();
  int exchange_component(int species);

  //   Initiates ghost cell communication with neighbors
  int exchange_start();
  int exchange_start_component(int species);

  //   Finishes ghost cell communication with neighbors
  int exchange_end();
  int exchange_end_component(int species);

  //   Set/Get data array for a given species
  float* GetData(int species);
  int SetData(int species, float *NewArray);

  //   Returns data values by location in 3-space and by species
  float operator()(int i, int j, int k, int s) {
    return data[s][(k*(Nx1+Ng1l+Ng1r)+j)*(Nx0+Ng0l+Ng0r)+i];
  };

  //   Returns dimensional size
  int size(int *n0, int *n1, int *n2, int *ns, int *g0l, 
	   int *g0r, int *g1l, int *g1r, int *g2l, int *g2r);

  //   Copies the values from a vector x (including ghost cells)
  int copy(EnzoVector *x);

  //   Copies the values from a vector x (including ghost cells)
  int copy_component(EnzoVector *x, int c);

  //   Vector linear sum operation, this = a*x + b*y
  int linearsum(float a, EnzoVector *x, float b, EnzoVector *y);

  //   Vector axpy operation, this += a*x
  int axpy(float a, EnzoVector *x);

  //   Vector axpy operation (single component), this += a*x
  int axpy_component(float a, EnzoVector *x, int c);

  //   Vector scale operation, this *= a
  int scale(float a);

  //   Vector component scale operation, this *= a
  int scale_component(int var, float a);

  //   Vector component log operation, this = log(this) || -12000 (if value == 0)
  int log_component(int var);

  //   Vector component exp operation, this = exp(this)
  int exp_component(int var);

  //   Vector constant operation, this(i) = a
  int constant(float a);

  //   Vector add const operation, this(i) += a
  int addconst(float a);

  //   Vector component add const operation, this(i) += a
  int addconst_component(int var, float a);

  //   Vector absolute value,  this(i) = |x(i)|
  int abs(EnzoVector *x);

  //   Vector product operation, this(i) = x(i)*y(i)
  int product(EnzoVector *x, EnzoVector *y);

  //   Vector quotient operation, this(i) = x(i)/y(i)
  //   [assumes y(i) != 0]
  int quotient(EnzoVector *x, EnzoVector *y);

  //   Vector minquotient operation, returns
  //   min(this(i)/y(i)) over all y(i)!=0
  float minquotient(EnzoVector *y);

  //   Vector inverse operation, this(i) = 1.0/x(i)
  //   [assumes x(i) != 0]
  int inverse(EnzoVector *x);

  //   Vector constraint checking operation, 
  //   this(i) = 0.0 where constraints true, 1.0 where false
  //   Test:  if c[i] =  2.0, then x[i] must be >  0.0
  //          if c[i] =  1.0, then x[i] must be >= 0.0
  //          if c[i] = -1.0, then x[i] must be <= 0.0
  //          if c[i] = -2.0, then x[i] must be <  0.0
  //   if all constraints satisfied, returns true
  bool constraintcheck(EnzoVector *c, EnzoVector *x);

  //   Vector dot-product,  dot(this,x)
  float dot(EnzoVector *x) const;

  //   Vector RMS norm,  sqrt(dot(this,this)/Nglobal)
  float rmsnorm() const;

  //   Component RMS norm,  sqrt(dot(this,this)/Nglobal)
  float rmsnorm_component(int var) const;

  //   Vector weighted RMS norm,  sqrt(dot(this*w,this*w)/Nglobal)
  float wrmsnorm(EnzoVector *w) const;

  //   Vector weighted L-2 norm, sqrt(dot(this*w,this*w))
  float wl2norm(EnzoVector *w) const;

  //   Vector L-1 norm,  sum(abs(this))
  float l1norm() const;

  //   Vector infinity (max) norm,  max(abs(this))
  float infnorm() const;
  
  //   Component infinity (max) norm,  max(abs(this[var]))
  float infnorm_component(int var) const;
  
  //   Relative pointwise difference,  max(abs(this(var)-x)/abs(this(var)))
  //   (assumes this vector is nonzero everywhere)
  float relative_difference(float *x, int var) const;
  
  //   Relative volumetric difference,  rmsnorm(abs(this(var)-x)/abs(this(var)))
  //   (assumes this vector is nonzero everywhere)
  float relative_vol_difference(float *x, int var) const;
  
  //   Vector minimum value,  min(this)
  float minval() const;
  
  //   Vector test routine
  int test();

};
  
#endif
