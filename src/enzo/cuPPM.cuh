/***********************************************************************
/
/  RECONSTRUCT EDGE STATES USING CUDA
/
/  written by: Peng Wang, NVIDIA
/  date:       May, 2012
/  modified1:
/
/  description:
/     
************************************************************************/
#include <stdio.h>
#include "fortran.def"
#include "transpose_cuda.cuh"

#define NBLOCK 128



const float ft = 4.0f/3.0f;

__forceinline__ __device__ float max3(float a1, float a2, float a3)
{
  float am = a1;
  if (am < a2) am = a2;
  if (am < a3) am = a3;
  return am;
}

__forceinline__ __device__ float min4(float a1, float a2, float a3, float a4)
{
  float am = a1;
  if (am > a2) am = a2;
  if (am > a3) am = a3;
  if (am > a4) am = a4;
  return am;
}

__forceinline__ __device__ float calc_flatten(float *pslice, float *uslice, 
                              float epsilon, float omega1, float omega2)
{
  float flattemp, qa, qb, wflag;
  qb = abs(pslice[1] - pslice[-1]) / min(pslice[1], pslice[-1]);
  if (qb > epsilon && uslice[-1] > uslice[1]) 
    wflag = 1.0f;
  else
    wflag = 0.0f;
  if (fabs(pslice[2] - pslice[-2]) /
      min(pslice[2], pslice[-2]) < epsilon)
    qa = 1.0f;
  else
    qa = (pslice[1] - pslice[-1]) /
      (pslice[2] - pslice[-2]);
  flattemp = min(1.0f, (qa-omega1)*omega2*wflag);
  flattemp = max(0, flattemp);
  return flattemp;
}


//
// Calculate flatten coefficient: i1-1:i2+1, cell-centered
//
__global__ void calc_flatten_kernel(float *uslice, float *pslice,
                                  int idim, int jdim, int kdim,
                                  int i1, int i2, int j1, int j2, int k1, int k2, 
                                  float *flatten)
{
  const float omega1 = 0.75f, omega2 = 10.0f, epsilon = 0.33f;
  
  float flatten_i, flatten_im1, flatten_ip1;
  // for (int k = k1; k <= k2; k++) {
  //   for (int j = j1; j <= j2; j++) {
  //     for (int i = i1-1; i <= i2+1; i++) {
  //       int idx3d = (k*jdim+j)*idim + i;
  int idx3d = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x + i1 - 1;
  int i = idx3d % idim;
  if (i >= i1-1 && i <= i2+1) {
    flatten_i = calc_flatten(pslice+idx3d, uslice+idx3d,
                             epsilon, omega1, omega2);
    if (i != i1 - 1)
      flatten_im1 = calc_flatten(pslice+idx3d-1, uslice+idx3d-1,
                                 epsilon, omega1, omega2);
    else 
      flatten_im1 = flatten_i;
    if (i != i2 + 1)
      flatten_ip1 = calc_flatten(pslice+idx3d+1, uslice+idx3d+1,
                                 epsilon, omega1, omega2);
    else
      flatten_ip1 = flatten_i;
    if (pslice[idx3d+1] - pslice[idx3d-1] < 0.0f)
      flatten[idx3d] = max(flatten_i, flatten_ip1);
    else
      flatten[idx3d] = max(flatten_i, flatten_im1);
    //if (fabs(flatten[idx3d]) > 1e-5) printf("%f %d\n", flatten[idx3d], i+1);
  }
  //     }
  //   }
  // }
}

//
//  Calculate diffusion coefficient: i1:i2+1, edge-centered
//
__global__ void calc_diffcoef_kernel(float *uslice, float *vslice, float *wslice, 
                                   int idim, int jdim, int kdim,
                                   int i1, int i2, int j1, int j2, int k1, int k2, 
                                   int nzz, int idir, int dimx, int dimy, int dimz,
                                   float *diffcoef)
{
  const float Kparam = 0.1f;
  for (int k = k1; k <= k2; k++) {
    int lk1 = (k1 < k && k < k2);
    for (int j = j1; j <= j2; j++) {
      int lj1 = (j1 < j && j < j2);
      for (int i = i1; i <= i2+1; i++) {
        int idx3d = (k*jdim+j)*idim + i;
        float vdiff1 = 0.0f, wdiff1 = 0.0f;
        //if (idir == 1) {
        //if (dimy > 1 && lj1)
        if (lj1)
            vdiff1 = (vslice[idx3d-dimx] + vslice[idx3d-1-dimx]) 
              - (vslice[idx3d+dimx] + vslice[idx3d-1+dimx]);
            //if (dimz > 1 && lk1)
        if (lk1)
            wdiff1 = (wslice[idx3d-dimx*dimy] + wslice[idx3d-1-dimx*dimy])
              - (wslice[idx3d+dimx*dimy] + wslice[idx3d-1+dimx*dimy]);
        float diff = uslice[idx3d-1] - uslice[idx3d];
        //if (lj1)
          diff += 0.25f*vdiff1;
          //if (lk1)
          diff += 0.25f*wdiff1;
          diffcoef[idx3d] = Kparam*max(0.0f, diff);
          //}
        // if (idir == 2) {
          
        //   if (dimz > 1 && lj1)
        //     vdiff1 = (v[idx3d-dimx*dimy] + v[idx3d-dimx-dimx*dimy])
        //       - (v[idx3d+dimx*dimy] + v[idx3d-dimx+dimx*dimy]);
        //   if (dimx > 1 && lk1)
        //     wdiff1 = (w[idx3d-1] + w[idx3d-1-dimx])
        //       - (w[idx3d+1] + w[idx3d+1-dimx]);
        // }
        // if (idir == 3) {
        //   if (dimx > 1 && lj1)
        //     vdiff1 = (v[idx3d-1
            
      }
    }
  }
}

// eqn 1.8
__forceinline__ __device__ void compute_deq(float qplus, float qmnus, float c1, float c2,
                            float &deq)
{
  if (qplus*qmnus > 0) {
    float qcent = c1*qplus + c2*qmnus;
    float qvanl = 2.0f*qplus*qmnus/(qmnus+qplus);
    float temp1 = min4(fabs(qcent), fabs(qvanl), 2.0f*fabs(qmnus), 2.0f*fabs(qplus));
    deq = temp1*sign(qcent);
  } else {
    deq = 0.0f;
  }
}

// eqn 1.10
__forceinline__ __device__ void monotonize_qlr(float q, float &ql, float &qr)
{
  float temp1 = (qr-q)*(q-ql);
  float temp2 = qr-ql;
  float temp3 = 6.0f*(q-0.5f*(qr+ql));
  if (temp1 <= 0) ql = q, qr = q;
  float temp22 = temp2*temp2;
  float temp23 = temp2*temp3;
  if (temp22 < temp23) 
    ql = 3.0f*q-2.0*qr;
  if (temp22 < -temp23)
    qr = 3.0f*q-2.0*ql;
}

// compute left and right states in edges using PPM
__device__ void intvar_cuda(float q_im3, float q_im2, float q_im1, 
                            float q_i, float q_ip1, float q_ip2,
                            int isteep, float steepen_im1, float steepen_i,
                            int iflatten, float flatten_im1, float flatten_i,
                            float c0_im1, float c0_i,
                            float c1, float c2, float c3,
                            float c4, float c5, float c6,
                            float char1, float char2,
                            float &dq_im1, float &dq_i,
                            float &ql_i, float &qr_im1,
                            float &q6_im1, float &q6_i,
                            float &qla, float &qra,
                            float &ql0, float &qr0)
{

  float ql_im1, qr_i;
  float deq_im2, deq_im1, deq_i, deq_ip1;
  float qplus, qmnus;
    
  // Compute average linear slopes (eqn 1.7 & 1.8)
  // Added van Leer slopes by JHW (June 2010, idea taken ATHENA)
  qplus=q_im1-q_im2, qmnus=q_im2-q_im3; 
  compute_deq(qplus, qmnus, c1, c2, deq_im2);
  qplus=q_i  -q_im1, qmnus=q_im1-q_im2; 
  compute_deq(qplus, qmnus, c1, c2, deq_im1);
  qplus=q_ip1-q_i  , qmnus=q_i  -q_im1; 
  compute_deq(qplus, qmnus, c1, c2, deq_i  );
  qplus=q_ip2-q_ip1, qmnus=q_ip1-q_i  ; 
  compute_deq(qplus, qmnus, c1, c2, deq_ip1);
  
  // Construct left and right cell values (eqn 1.6)
  ql_im1 = c3*q_im2 + c4*q_im1 + c5*deq_im2 + c6*deq_im1;
  ql_i   = c3*q_im1 + c4*q_i   + c5*deq_im1 + c6*deq_i  ;
  qr_i   = c3*q_i   + c4*q_ip1 + c5*deq_i   + c6*deq_ip1;
  qr_im1 = ql_i;

  // Steepen if asked for (use precomputed steepening parameter)
  if (isteep) {
    ql_im1 = (1.0f-steepen_im1)*ql_im1 + steepen_im1*(q_im2+0.5f*deq_im2);
    ql_i   = (1.0f-steepen_i  )*ql_i   + steepen_i  *(q_im1+0.5f*deq_im1);
    qr_im1 = (1.0f-steepen_im1)*qr_im1 + steepen_im1*(q_i  -0.5f*deq_i  );
    qr_i   = (1.0f-steepen_i  )*qr_i   + steepen_i  *(q_ip1-0.5f*deq_ip1);
  }

  if (iflatten) {
    ql_im1 = q_im1*flatten_im1 + ql_im1*(1.0f-flatten_im1);
    ql_i   = q_i  *flatten_i   + ql_i  *(1.0f-flatten_i  );
    qr_im1 = q_im1*flatten_im1 + qr_im1*(1.0f-flatten_im1);
    qr_i   = q_i  *flatten_i   + qr_i  *(1.0f-flatten_i  );
  }

  // Monotonize again (eqn 1.10)
  monotonize_qlr(q_im1, ql_im1, qr_im1);
  monotonize_qlr(q_i  , ql_i  , qr_i  );  

  // Ensure that L/R values lie between neighboring cell-centered values
  // Take from ATHENA, lr_states
  ql_im1 = max(min(q_im1, q_im2), ql_im1);
  ql_im1 = min(max(q_im1, q_im2), ql_im1);
  qr_im1 = max(min(q_im1, q_i  ), qr_im1);
  qr_im1 = min(max(q_im1, q_i  ), qr_im1);
  
  ql_i = max(min(q_i, q_im1), ql_i);
  ql_i = min(max(q_i, q_im1), ql_i);
  qr_i = max(min(q_i, q_ip1), qr_i);
  qr_i = min(max(q_i, q_ip1), qr_i);

  // Now construct left and right edge values (eqn 1.12 and eqn 3.3)
  q6_im1 = 6.0f*(q_im1 - 0.5f*(ql_im1 + qr_im1));
  q6_i   = 6.0f*(q_i   - 0.5f*(ql_i   + qr_i  ));

  dq_im1 = qr_im1 - ql_im1;
  dq_i   = qr_i   - ql_i  ;

  qla = qr_im1 - char1*(dq_im1 - (1.0f - ft*char1)*q6_im1);
  qra = ql_i   + char2*(dq_i   + (1.0f - ft*char2)*q6_i  );

  ql0 = qr_im1 - c0_im1*(dq_im1 - (1.0f - ft*c0_im1)*q6_im1);
  qr0 = ql_i   - c0_i  *(dq_i   + (1.0f + ft*c0_i  )*q6_i  );
}

__global__ void inteuler_kernel(float *dslice, float *pslice, int gravity,
                              float *grslice, float *geslice,
                              float *uslice, float *vslice, float *wslice, 
                              float dx, float *flatten,
                              int idim, int jdim, int kdim,
                              int i1, int i2, int j1, int j2, int k1, int k2,
                              int idual, float eta1, float eta2, 
                              int isteep, int iflatten, 
                              int iconsrec, int iposrec,
                              float dt, float gamma, int ipresfree,
                              float *dls, float *drs, 
                              float *pls, float *prs,
                              float *gels, float *gers,
                              float *uls, float *urs,
                              float *vls, float *vrs,
                              float *wls, float *wrs,
                              int ncolor, float *colslice,
                              float *colls, float *colrs)
{
  float c0_im1, c0_i, cm_im1, cm_i, cp_im1, cp_i;
  float c1, c2, c3, c4, c5, c6;
  float cs_im1, cs_i, char1, char2;
  float dd_im1, dd_i, dl_i, dr_im1, d6_im1, d6_i, dla, dra, dl0, dr0;
  float dp_im1, dp_i, pl_i, pr_im1, p6_im1, p6_i, pla, pra, pl0, pr0;
  float du_im1, du_i, ul_i, ur_im1, u6_im1, u6_i, ula, ura, ul0, ur0;
  float dv_im1, dv_i, vl_i, vr_im1, v6_im1, v6_i, vla, vra, vl0, vr0;
  float dw_im1, dw_i, wl_i, wr_im1, w6_im1, w6_i, wla, wra, wl0, wr0;
  float dge_im1, dge_i, gel_i, ger_im1, ge6_im1, ge6_i, gela, gera, gel0, ger0;
  float dcol_im1, dcol_i, coll_i, 
    colr_im1, col6_im1, col6_i, 
    colla, colra, coll0, colr0;
  float dx2i;
  float steepen_im1, steepen_i, flatten_im1, flatten_i;

  float plm, prm, plp, prp, ulm, urm, ulp, urp, cla, cra;
  float betalp, betalm, betal0, betarp, betarm, betar0;
  float f1;
  const int size3d = idim*jdim*kdim;

  // for (int k = k1; k <= k2; k++) {
  //   for (int j = j1; j <= j2; j++) {
  //     for (int i = i1; i <= i2+1; i++) {
  //int idx3d = (k*jdim+j)*idim + i;
  int idx3d = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x + i1;
  //int i = idx3d % idim;
  // int j = (idx3d % (idim*jdim)) / idim;
  // int k = idx3d / (idim*jdim);

  if (idx3d <= size3d - 3) {

    c1 = 0.5f;
    c2 = 0.5f;
    c3 = 0.5f;
    c4 = 0.5f;
    c5 = 1.0f/6.0f;
    c6 = -1.0f/6.0f;
    
    dx2i = 0.5f/dx;

    // Precompute steepening coefficients of needed (eqns 1.14-1.17, plus 3.2)
    if (isteep) {
      float d2d_im2, d2d_im1, d2d_i, d2d_ip1;
      d2d_im2 = (dslice[idx3d-1] - 2*dslice[idx3d-2] + dslice[idx3d-3])/6.0f;
      d2d_im1 = (dslice[idx3d  ] - 2*dslice[idx3d-1] + dslice[idx3d-2])/6.0f;
      d2d_i   = (dslice[idx3d+1] - 2*dslice[idx3d  ] + dslice[idx3d-1])/6.0f;
      d2d_ip1 = (dslice[idx3d+2] - 2*dslice[idx3d+1] + dslice[idx3d  ])/6.0f;
      float qc_im1, qc_i, s1_im1, s1_i, s2_im1, s2_i, qa_im1, qa_i, qb_im1, qb_i;
      qc_im1 = fabs(dslice[idx3d  ] - dslice[idx3d-2])
        - 0.01f*min(fabs(dslice[idx3d  ]), fabs(dslice[idx3d-2]));
      qc_i   = fabs(dslice[idx3d+1] - dslice[idx3d-1])
        - 0.01f*min(fabs(dslice[idx3d+1]), fabs(dslice[idx3d-1]));
      s1_im1 = (d2d_im2 - d2d_i  )/(dslice[idx3d  ] - dslice[idx3d-2] + tiny);
      s1_i   = (d2d_im1 - d2d_ip1)/(dslice[idx3d+1] - dslice[idx3d-1] + tiny);
      if (d2d_i  *d2d_im2 > 0) s1_im1 = 0.0f;
      if (d2d_ip1*d2d_im1 > 0) s1_i   = 0.0f;
      if (qc_im1 <= 0.0f) s1_im1 = 0.0f;
      if (qc_i   <= 0.0f) s1_i   = 0.0f;
      s2_im1 = max(0.0f, min(20.0f*(s1_im1-0.05f), 1.0f));
      s2_i   = max(0.0f, min(20.0f*(s1_i  -0.05f), 1.0f));
      qa_im1 = fabs(dslice[idx3d  ] - dslice[idx3d-2]) / 
        min(dslice[idx3d  ], dslice[idx3d-2]);
      qa_i   = fabs(dslice[idx3d+1] - dslice[idx3d-1]) /
        min(dslice[idx3d+1], dslice[idx3d-1]);
      qb_im1 = fabs(pslice[idx3d  ] - pslice[idx3d-2]) /
        min(pslice[idx3d  ], pslice[idx3d-2]);
      qb_i   = fabs(pslice[idx3d+1] - pslice[idx3d-1]) /
        min(pslice[idx3d+1], pslice[idx3d-1]);
      if (gamma*0.1f*qa_im1 >= qb_im1)
        steepen_im1 = s2_im1;
      else
        steepen_im1 = 0.0f;
      if (gamma*0.1f*qa_i  >= qb_i)
        steepen_i = s2_i;
      else
        steepen_i = 0.0f;
    }

    if (iflatten) {
      flatten_im1 = flatten[idx3d-1];
      flatten_i   = flatten[idx3d];
    }

    cs_im1 = sqrt(gamma*pslice[idx3d-1]/dslice[idx3d-1]);
    cs_i   = sqrt(gamma*pslice[idx3d  ]/dslice[idx3d  ]);
    char1 = max(0.0f,  dt*(uslice[idx3d-1]+cs_im1)*dx2i);
    char2 = max(0.0f, -dt*(uslice[idx3d  ]-cs_i  )*dx2i);
    cm_im1 = dt*(uslice[idx3d-1]-cs_im1)*dx2i;
    c0_im1 = dt*(uslice[idx3d-1]       )*dx2i;
    cp_im1 = dt*(uslice[idx3d-1]+cs_im1)*dx2i;
    cm_i   = dt*(uslice[idx3d  ]-cs_i  )*dx2i;
    c0_i   = dt*(uslice[idx3d  ]       )*dx2i;
    cp_i   = dt*(uslice[idx3d  ]+cs_i  )*dx2i;
        
    if (iconsrec == 1) {
    } else {
      intvar_cuda(dslice[idx3d-3], dslice[idx3d-2], dslice[idx3d-1],
                  dslice[idx3d  ], dslice[idx3d+1], dslice[idx3d+2],
                  isteep, steepen_im1, steepen_i,
                  iflatten, flatten_im1, flatten_i,
                  c0_im1, c0_i, 
                  c1, c2, c3, c4, c5, c6,
                  char1, char2,
                  dd_im1, dd_i, dl_i, dr_im1, d6_im1, d6_i,
                  dla, dra, dl0, dr0);
      intvar_cuda(pslice[idx3d-3], pslice[idx3d-2], pslice[idx3d-1],
                  pslice[idx3d  ], pslice[idx3d+1], pslice[idx3d+2],
                  isteep, steepen_im1, steepen_i,
                  iflatten, flatten_im1, flatten_i,
                  c0_im1, c0_i, 
                  c1, c2, c3, c4, c5, c6,
                  char1, char2,
                  dp_im1, dp_i, pl_i, pr_im1, p6_im1, p6_i,
                  pla, pra, pl0, pr0);
      intvar_cuda(uslice[idx3d-3], uslice[idx3d-2], uslice[idx3d-1],
                  uslice[idx3d  ], uslice[idx3d+1], uslice[idx3d+2],
                  isteep, steepen_im1, steepen_i,
                  iflatten, flatten_im1, flatten_i,
                  c0_im1, c0_i, 
                  c1, c2, c3, c4, c5, c6,
                  char1, char2,
                  du_im1, du_i, ul_i, ur_im1, u6_im1, u6_i,
                  ula, ura, ul0, ur0);
      intvar_cuda(vslice[idx3d-3], vslice[idx3d-2], vslice[idx3d-1],
                  vslice[idx3d  ], vslice[idx3d+1], vslice[idx3d+2],
                  isteep, steepen_im1, steepen_i,
                  iflatten, flatten_im1, flatten_i,
                  c0_im1, c0_i, 
                  c1, c2, c3, c4, c5, c6,
                  char1, char2,
                  dv_im1, dv_i, vl_i, vr_im1, v6_im1, v6_i,
                  vla, vra, vl0, vr0);
      intvar_cuda(wslice[idx3d-3], wslice[idx3d-2], wslice[idx3d-1],
                  wslice[idx3d  ], wslice[idx3d+1], wslice[idx3d+2],
                  isteep, steepen_im1, steepen_i,
                  iflatten, flatten_im1, flatten_i,
                  c0_im1, c0_i, 
                  c1, c2, c3, c4, c5, c6,
                  char1, char2,
                  dw_im1, dw_i, wl_i, wr_im1, w6_im1, w6_i,
                  wla, wra, wl0, wr0);
    }
    if (idual == 1)
      intvar_cuda(geslice[idx3d-3], geslice[idx3d-2], geslice[idx3d-1],
                  geslice[idx3d  ], geslice[idx3d+1], geslice[idx3d+2],
                  isteep, steepen_im1, steepen_i,
                  iflatten, flatten_im1, flatten_i,
                  c0_im1, c0_i, 
                  c1, c2, c3, c4, c5, c6,
                  char1, char2,
                  dge_im1, dge_i, gel_i, ger_im1, ge6_im1, ge6_i,
                  gela, gera, gel0, ger0);
    //
    // Correct the initial guess from the linearized gas equations
    //
    //  First, compute average over characteristic domain of dependance (3.5)
    //
    plm = pr_im1 - cm_im1*(dp_im1 - (1.0f - ft*cm_im1)*p6_im1);
    prm = pl_i   - cm_i  *(dp_i   + (1.0f + ft*cm_i  )*p6_i  );
    plp = pr_im1 - cp_im1*(dp_im1 - (1.0f - ft*cp_im1)*p6_im1);
    prp = pl_i   - cp_i  *(dp_i   + (1.0f + ft*cp_i  )*p6_i  );

    ulm = ur_im1 - cm_im1*(du_im1 - (1.0f - ft*cm_im1)*u6_im1);
    urm = ul_i   - cm_i  *(du_i   + (1.0f + ft*cm_i  )*u6_i  );
    ulp = ur_im1 - cp_im1*(du_im1 - (1.0f - ft*cp_im1)*u6_im1);
    urp = ul_i   - cp_i  *(du_i   + (1.0f + ft*cp_i  )*u6_i  );

    cla = sqrt(max(gamma*pla*dla, 0.0f));
    cra = sqrt(max(gamma*pra*dra, 0.0f));
    //
    // a) left size
    //
    f1 = 1.0f/cla;
    betalp = (ula - ulp) + (pla - plp)*f1;
    betalm = (ula - ulm) - (pla - plm)*f1;
    betal0 = (pla - pl0)*f1*f1 + 1.0f/dla - 1.0f/dl0;
    //
    // Add gravity component
    //
    if (gravity == 1) {
      betalp -= 0.25f*dt*(grslice[idx3d-1] + grslice[idx3d]);
      betalm -= 0.25f*dt*(grslice[idx3d-1] + grslice[idx3d]);
    }
    f1 = 0.5f/cla;
    betalp = -betalp*f1;
    betalm =  betalm*f1;
    if (cp_im1 <= 0) betalp = 0.0f;
    if (cm_im1 <= 0) betalm = 0.0f;
    if (c0_im1 <= 0) betal0 = 0.0f;
    //
    // b) right side
    //
    f1 = 1.0f/cra;
    betarp = (ura - urp) + (pra - prp)*f1;
    betarm = (ura - urm) - (pra - prm)*f1;
    betar0 = (pra - pr0)*f1*f1 + 1.0f/dra - 1.0f/dr0;
    //
    // Add gravity component
    //
    if (gravity == 1) {
      betarp -= 0.25f*dt*(grslice[idx3d-1] + grslice[idx3d]);
      betarm -= 0.25f*dt*(grslice[idx3d-1] + grslice[idx3d]);
    }
    f1 = 0.5f/cra;
    betarp = -betarp*f1;
    betarm =  betarm*f1;
    if (cp_i >= 0) betarp = 0.0f;
    if (cm_i >= 0) betarm = 0.0f;
    if (c0_i >= 0) betar0 = 0.0f;
    //
    // Finally, combine to create corrected left/right states (eq. 3.6)
    //
    pls[idx3d] = pla + (betalp + betalm)*cla*cla;
    // if (i == 36 && j == 5 && k == 32) 
    //   printf("pls0: %f,%f,%f,%f,%f\n", pls[idx3d],pla,betalp,betalm,cla);
    prs[idx3d] = pra + (betarp + betarm)*cra*cra;
    uls[idx3d] = ula + (betalp - betalm)*cla;
    urs[idx3d] = ura + (betarp - betarm)*cra;
    dls[idx3d] = 1.0f/(1.0f/dla - (betal0+betalp+betalm));
    drs[idx3d] = 1.0f/(1.0f/dra - (betar0+betarp+betarm));
    //
    // Take the appropriate state from the advected variables
    //
    if (uslice[idx3d-1] <= 0.0f) {
      vls [idx3d] = vla;
      wls [idx3d] = wla;
      gels[idx3d] = gela;
    } else {
      vls [idx3d] = vl0;
      wls [idx3d] = wl0;
      gels[idx3d] = gel0;
    }
    if (uslice[idx3d] >= 0.0f) {
      vrs [idx3d] = vra;
      wrs [idx3d] = wra;
      gers[idx3d] = gera;
    } else {
      vrs [idx3d] = vr0;
      wrs [idx3d] = wr0;
      gers[idx3d] = ger0;
    }
    for (int c = 0; c < ncolor; c++) {
      intvar_cuda(colslice[c*size3d+idx3d-3], colslice[c*size3d+idx3d-2], 
                  colslice[c*size3d+idx3d-1], colslice[c*size3d+idx3d  ], 
                  colslice[c*size3d+idx3d+1], colslice[c*size3d+idx3d+2],
                  isteep, steepen_im1, steepen_i,
                  iflatten, flatten_im1, flatten_i,
                  c0_im1, c0_i, 
                  c1, c2, c3, c4, c5, c6,
                  char1, char2,
                  dcol_im1, dcol_i, coll_i, colr_im1, 
                  col6_im1, col6_i, colla, colra, 
                  coll0, colr0);         

      if (uslice[idx3d-1] <= 0.0f) 
        colls[c*size3d+idx3d] = colla;
      else
        colls[c*size3d+idx3d] = coll0;
      if (uslice[idx3d] >= 0.0f) 
        colrs[c*size3d+idx3d] = colra;
      else
        colrs[c*size3d+idx3d] = colr0;
    }

    if (idual == 1) {
      if (gamma*pla/dla < eta2*ula*ula || 
          max3(fabs(cm_im1), fabs(c0_im1), fabs(cp_im1)) < 1e-3 ||
          dls[idx3d]/dla > 5.0f) {
        pls[idx3d] = pla;
        uls[idx3d] = ula;
        dls[idx3d] = dla;
      }
      if (gamma*pra/dra < eta2*ura*ura ||
          max3(fabs(cm_i), fabs(c0_i), fabs(cp_i)) < 1e-3 ||
          drs[idx3d]/dra > 5.0f) {
        prs[idx3d] = pra;
        urs[idx3d] = ura;
        drs[idx3d] = dra;
      }
    }

    // Enforce minimum values
    pls[idx3d] = max(pls[idx3d], tiny);
    prs[idx3d] = max(prs[idx3d], tiny);
    dls[idx3d] = max(dls[idx3d], tiny);
    drs[idx3d] = max(drs[idx3d], tiny);
    // 
    // If approximating pressure free conditions, then the density
    //  should be reset to the pre-corrected state
    //
    if (ipresfree == 1) {
      dls[idx3d] = dla;
      drs[idx3d] = dra;
    }
        
    // if (k == 3 && j == 7) {
    //   printf("%f ", urs[idx3d]);
    //   //printf("%f %f %f\n", betarp, betarm, betar0);
    // }
  }

  //     } // endfor i
  //   } // endfor j
  // } // endfor k

}

__global__ void twoshock_kernel(float *dls, float *drs,
                              float *pls, float *prs,
                              float *uls, float *urs,
                              int idim, int jdim, int kdim,
                              int i1, int i2, int j1, int j2, int k1, int k2,
                              float dt, float gamma, float pmin, int ipresfree,
                              float *pbar, float *ubar, 
                              int gravity, float *grslice,
                              int idual, float eta1)
{
  float cl, cr, dpdul, dpdur, ps, ubl, ubr, qa, zl, zr;
  float old_ps, delta_ps;
  int mask;
  const float tolerance = 1e-7;
  const int numiter = 8;

  const int size3d = idim*jdim*kdim;
  int idx3d = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x + i1;

  qa = (gamma + 1.0f)/(2.0f*gamma);
  if (idx3d <= size3d - 3) {
    //
    // If pressfree conditions are needed, set pbar to zero and ubar to the
    //   average of left and right velocity states
    //
    if (ipresfree == 1) {
      pbar[idx3d] = pmin;
      ubar[idx3d] = 0.5f*(uls[idx3d] + urs[idx3d]);
      pls[idx3d]  = pmin;
      prs[idx3d]  = pmin;
    } 
    // otherwise, solve Riemann problem
    else {
      //
      // First guess at pbar and left- and right-ubar
      //   (or is it a plus?: & + cr(i)*cl(i)*(uls(i,j) - urs(i,j))
      // it is a plus: see van Leer 1979, JCP 32, 101 (eq. 60 on page 109).
      //
      if (gravity == 99) {
        cl = sqrt(gamma*pls[idx3d]*dls[idx3d]);
        cr = sqrt(gamma*prs[idx3d]*drs[idx3d]);
        ps = (cr*pls[idx3d] + cl*prs[idx3d]
              + cr*cl*(uls[idx3d] - urs[idx3d]
                       + 0.5f*dt*(grslice[idx3d-1] - grslice[idx3d])
                ))/(cr+cl);
        //    ^
        //    changed to + June 2003
        if (ps < pmin) ps = pmin;
      } else {
        cl = sqrt(gamma*pls[idx3d]*dls[idx3d]);
        cr = sqrt(gamma*prs[idx3d]*drs[idx3d]);
        ps = (cr*pls[idx3d] + cl*prs[idx3d]
              + cr*cl*(uls[idx3d] - urs[idx3d]))/(cr+cl);
        //    ^
        //    changed to + June 2003
        if (ps < pmin) ps = pmin;
      }
      //
      // Newton iterations to compute successive guesses for pbar, 
      //   l and r ubar
      //   (there are some gravity changes in here but they have not really
      //    been worked out yet).
      //
      old_ps = ps;
      delta_ps = ps;
      mask = 1;
      for (int n = 2; n <= numiter; n++) {
        //  zl, zr depend on ps only
        //  ubl, ubr depend on ps only
        //  dpdul, dpdur depend on ps only
        if (mask > 0) {
          zl  = cl*sqrt((1.0f+qa*(ps/pls[idx3d]-1.0f)));
          zr  = cr*sqrt((1.0f+qa*(ps/prs[idx3d]-1.0f)));
          ubl = uls[idx3d] - (ps-pls[idx3d])/zl;
          ubr = urs[idx3d] + (ps-prs[idx3d])/zr;
        }
            
        if (mask > 0) {
          dpdul = -4.0f*zl*zl*zl/dls[idx3d]
            /(4.0f*zl*zl/dls[idx3d] - (gamma+1.0f)*(ps-pls[idx3d]));
          dpdur =  4.0f*zr*zr*zr/drs[idx3d]
            /(4.0f*zr*zr/drs[idx3d] - (gamma+1.0f)*(ps-prs[idx3d]));
          ps = ps + (ubr-ubl)*dpdur*dpdul/(dpdur-dpdul);
          if (ps < pmin) ps = pmin;
          delta_ps = ps - old_ps;
          old_ps = ps;
          if (fabs(delta_ps / ps) < tolerance) 
            mask = 0;
        }
      }
      //
      // Compute final values of resolved state
      //
      if (ps < pmin) ps = min(pls[idx3d], prs[idx3d]);
      pbar[idx3d] = ps;
      ubar[idx3d] = ubl + (ubr-ubl)*dpdur/(dpdur-dpdul);
    }
  }
  //     }
  //   }
  // }
  
}

__global__ void flux_twoshock_kernel(float *dslice, float *eslice, float *geslice, 
                                   float *uslice, float *vslice, float *wslice,
                                   float dx, float *diffcoef, 
                                   int idim, int jdim, int kdim,
                                   int i1, int i2, int j1, int j2, 
                                   int k1, int k2,
                                   float dt, float gamma, int idiff,
                                   int idual, float eta1, int ifallback, 
                                   float *dls, float *drs, 
                                   float *pls, float *prs,
                                   float *gels, float *gers,
                                   float *uls, float *urs,
                                   float *vls, float *vrs,
                                   float *wls, float *wrs,
                                   float *pbar, float *ubar,
                                   float *df, float *ef, float *uf,
                                   float *vf, float *wf, float *gef, float *ges, float *gesf,
                                   int ncolor, float *colslice, 
                                   float *colls, float *colrs, float *colf)
{
  float qc, frac;
  float sn, u0, p0, d0, c0, z0, dbar, cbar, l0, lbar, ub, vb, wb, db, pb,
    eb, geb, colb, upb, dub, duub, duvb, duwb, dueb, dugeb;
  //__shared__ float s_ub[NBLOCK+1];
  const int size3d = idim*jdim*kdim;
  int idx3d = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x + i1;
  
  float qa = (gamma + 1.0f)/(2.0f*gamma);
#ifdef RAREFACTION1
  float qb = (gamma - 1.0f)/(gamma + 1.0f);
  int cb;
#endif
  // for (int k = k1; k <= k2; k++) {
  //   for (int j = j1; j <= j2; j++) {
  //     for (int i = i1; i <= i2+1; i++) {
  //       int idx3d = (k*jdim+j)*idim + i;
  //int i = idx3d % idim;
  if (idx3d <= size3d - 3) {
    //
    // Evaluate time-averaged quantities
    //  (see Colella, Siam J Sci Stat Comput 1982, 3, 77. Appendix)
    //
    sn = sign(-ubar[idx3d]);
    //
    // Collect values of interest depending on which way fluid is flowing
    //
    if (sn < 0.0f) {
      u0 = uls[idx3d];
      p0 = pls[idx3d];
      d0 = dls[idx3d];
    } else {
      u0 = urs[idx3d];
      p0 = prs[idx3d];
      d0 = drs[idx3d];
    }
    c0 = sqrt(max(gamma*p0/d0, tiny));
    z0 = c0*d0*sqrt(max(1.0f+qa*(pbar[idx3d]/p0-1.0f), tiny));
    //
    // Compute equivalent bar (inside shock & rarefaction) values for density
    //   and sound speed
    //
    dbar = 1.0f/(1.0f/d0 - (pbar[idx3d]-p0)/max(z0*z0, tiny));
    cbar = sqrt(max(gamma*pbar[idx3d]/dbar, tiny));
    //
    // Find lambda values for the shock and rarefaction
    //
    if (pbar[idx3d] < p0) {
      l0 = u0*sn + c0;
      lbar = sn*ubar[idx3d] + cbar;
    } else {
      l0 = u0*sn + z0/d0;
      lbar = l0;
    }
    //
    // Compute values for inside a rarefaction fan
    //   (method described in Colella, 1982)
    // 
#define RAREFACTION2 //AK
#ifdef RAREFACTION2
    //
    // Compute values for inside a rarefaction fan
    //   (linear interpolation between end states, as suggested in PPM ref
    //
    frac = l0 - lbar;
    if (frac < tiny) frac = tiny;
    if (frac > 1.0f) frac = 1.0f;
    frac = -lbar/frac;
    frac = min(max(frac, 0.0f), 1.0f);
    pb = p0*frac + pbar[idx3d]*(1.0f - frac);
    db = d0*frac + dbar       *(1.0f - frac);
    ub = u0*frac + ubar[idx3d]*(1.0f - frac);
        
#endif
    //
    // Cull appropriate states depending on where eulerian position is in solution
    //   (lbar >= 0 --> inside post-shock region,
    //    l0   <  0 --> outside shock/rarefaction wave,
    //    otherwise --> inside rarefaction wave).
    //
    if (lbar >= 0.0f) { // AK
      pb = pbar[idx3d];
      db = dbar;
      ub = ubar[idx3d];
    }
    if (l0 < 0.0f) { // AK
      pb = p0;
      db = d0;
      ub = u0;
    }
    //
    // Collect values of interest depending on which way fluid is flowing
    //   (note: new placement for this - uses ub instead of ubar for consistency)
    //
    if (ub > 0.0f) { // AK
      vb = vls[idx3d];
      wb = wls[idx3d];
      geb = gels[idx3d];
    } else {
      vb = vrs[idx3d];
      wb = wrs[idx3d];
      geb = gers[idx3d];
    }
    //
    // Calculate total specific energy corresponding to this state
    //   (and the specific gas energy).
    //
    eb = pb/((gamma-1.0f)*db) + 0.5f*(ub*ub + vb*vb + wb*wb);
    //
    // Compute terms in differenced hydro equations (eq. 3.1)
    //
    if (idiff != 0) {
      //
      // ...with diffusion
      //
      upb = pb*ub;
      dub = ub*db; // AK
      duub = dub*ub + diffcoef[idx3d]*
        (dslice[idx3d-1]*uslice[idx3d-1] - dslice[idx3d]*uslice[idx3d]);
      duvb = dub*vb;
      duwb = dub*wb;
      //
      // (should we add diffusion to the cross velocities? I doubt it)
      // I don't. This doubt kills Noh test problem in 2D at high resolution.
      // Diffusion has to be added to cross velocities. !AK May 2005.
      //
      duvb = dub*vb + diffcoef[idx3d]*
        (dslice[idx3d-1]*vslice[idx3d-1] - dslice[idx3d]*vslice[idx3d]);
      duwb = dub*wb + diffcoef[idx3d]*
        (dslice[idx3d-1]*wslice[idx3d-1] - dslice[idx3d]*wslice[idx3d]);
      dueb = dub*eb + diffcoef[idx3d]*
        (dslice[idx3d-1]*eslice[idx3d-1] - dslice[idx3d]*eslice[idx3d]);
      //
      // This update must be the last !AK
      //
      dub = dub + diffcoef[idx3d]*(dslice[idx3d-1] - dslice[idx3d]);
      //
      // If using dual energy formalism, compute dugeb
      //
      if (idual == 1)
        dugeb = dub*geb + diffcoef[idx3d]*
          (dslice[idx3d-1]*geslice[idx3d-1] - dslice[idx3d]*geslice[idx3d]);
    } else { 
      //
      // ...and without 
      //
      upb  =  pb*ub;
      dub  =  ub*db;
      duub = dub*ub;
      duvb = dub*vb;
      duwb = dub*wb;
      dueb = dub*eb;

      if (idual == 1) 
        dugeb = dub*geb;
    }
    //
    // Copy into flux slices (to return to caller)
    //
    qc = dt/dx;
    df[idx3d] = qc*dub;
    ef[idx3d] = qc*(dueb + upb);
    uf[idx3d] = qc*(duub + pb);
    vf[idx3d] = qc*duvb;
    wf[idx3d] = qc*duwb;
        
    for (int c = 0; c < ncolor; c++) {
      if (ub > 0)
        colb = colls[idx3d+c*size3d];
      else
        colb = colrs[idx3d+c*size3d];
      colf[idx3d+c*size3d] = dt*ub*colb;
    }
    //
    // Do the same for the gas energy if using the dual energy formalism
    //   (note that we do not include the source term)
    // JHW (Jul 2010): Moved the source term to here.
    //
    if (idual == 1) {
      gef[idx3d] = dt*dugeb;
      gesf[idx3d] = ub;
    } 
  }

}

__global__ void flux_hllc_kernel(float *dslice, float *eslice, float *geslice, 
                               float *uslice, float *vslice, float *wslice,
                               float dx, float *diffcoef, 
                               int idim, int jdim, int kdim,
                               int i1, int i2, int j1, int j2, 
                               int k1, int k2,
                               float dt, float gamma, int idiff,
                               int idual, float eta1, int ifallback,
                               float *dls, float *drs, 
                               float *pls, float *prs,
                               float *gels, float *gers,
                               float *uls, float *urs,
                               float *vls, float *vrs,
                               float *wls, float *wrs,
                               float *df, float *ef, float *uf,
                               float *vf, float *wf, float *gef, float *ges, float *gesf,
                               int ncolor, float *colslice, 
                               float *colls, float *colrs, float *colf)
{
  // shared memory to store sl*gesl+sr*gesr for idual
  // extern __shared__ float smem[]; 
  
  const int size3d = idim*jdim*kdim;
  int idx3d = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x + i1;
  float gamma1 = gamma - 1.0f;
  float gamma1i = 1.0f / gamma1;

  // for (int k = k1; k <= k2; k++) {
  //   for (int j = j1; j <= j2; j++) {
  //     for (int i = i1; i <= i2+1; i++) {
  //       int idx3d = (k*jdim+j)*idim + i;
  // int i = idx3d % idim;
  // int j = (idx3d % (idim*jdim)) / idim;
  // int k = idx3d / (idim*jdim);
  if (idx3d <= size3d - 3) {
    //
    // Compute the Roe-averaged data from the left and right states
    //
    float sqrtdl = sqrt(dls[idx3d]);
    float sqrtdr = sqrt(drs[idx3d]);
    float isdlpdr = 1.0f / (sqrtdl + sqrtdr);
    float vroe1 = (sqrtdl*uls[idx3d] + sqrtdr*urs[idx3d])*isdlpdr;
    float vroe2 = (sqrtdl*vls[idx3d] + sqrtdr*vrs[idx3d])*isdlpdr;
    float vroe3 = (sqrtdl*wls[idx3d] + sqrtdr*wrs[idx3d])*isdlpdr;
    float v2 = vroe1*vroe1 + vroe2*vroe2 + vroe3*vroe3;
    //
    // The enthalpy H=(E+P)/d is averaged for adiabatic flows, rather
    // than E or P directly.
    // sqrtql*hl = sqrtdl*(el+pl)/dl = (el+pl)/sqrtdl
    //
    float el = gamma1i*pls[idx3d] + 0.5f*dls[idx3d]*
      (uls[idx3d]*uls[idx3d]+vls[idx3d]*vls[idx3d]+wls[idx3d]*wls[idx3d]);
    float er = gamma1i*prs[idx3d] + 0.5f*drs[idx3d]*
      (urs[idx3d]*urs[idx3d]+vrs[idx3d]*vrs[idx3d]+wrs[idx3d]*wrs[idx3d]);
    float hroe = ((el+pls[idx3d])/sqrtdl + (er+prs[idx3d])/sqrtdr)*isdlpdr;
    //
    // Compute characteristics using Roe-averaged values
    //
    float cs = sqrt(gamma1*max(hroe-0.5f*v2, tiny));
    float char1 = vroe1 - cs;
    float char2 = vroe1 + cs;
    //
    // Compute the min/max wave speeds
    //
    float csl0 = sqrt(gamma*pls[idx3d]/dls[idx3d]);
    float csr0 = sqrt(gamma*prs[idx3d]/drs[idx3d]);
    float csl = min(uls[idx3d]-csl0, char1);
    float csr = max(urs[idx3d]+csr0, char2);
    float bm = min(csl, 0.0f);
    float bp = max(csr, 0.0f);
    //
    // Compute the contact wave speed (cw) and pressure (cp). We do this
    // for all cases because we need to correct the momentum and energy
    // flux along the contact
    //
    float tl = pls[idx3d] - (csl - uls[idx3d])*dls[idx3d]*uls[idx3d];
    float tr = prs[idx3d] - (csr - urs[idx3d])*drs[idx3d]*urs[idx3d];
    float dl =  dls[idx3d]*(csl - uls[idx3d]);
    float dr = -drs[idx3d]*(csr - urs[idx3d]);
    float q1 = 1.0f / (dl+dr);
    float cw = (tr - tl)*q1;
    float cp = (dl*tr + dr*tl)*q1;
    float sl, sr, sm;
    if (cw >= 0.0f) {
      sl = cw/(cw-bm);
      sr = 0.0f;
      sm = -bm/(cw-bm);
    } else {
      sl = 0.0f;
      sr = -cw/(bp-cw);
      sm = bp/(bp-cw);
    }
    cp = max(cp, 0.0f);
    //
    // Compute diffusion terms 
    //
    if (idiff != 0) {
      // not implemented now
    }    
    // 
    //    Compute the left and right fluxes along the characteristics.
    //     Compute the color fluxes afterwards (loop n, then i) for cache
    //     efficiency.  For dual energy formalism, do not include the source
    //     term.  Treat like an advected quantity.
    //
    float dubl = dls[idx3d]*uls[idx3d];
    float dubr = drs[idx3d]*urs[idx3d];
    float dfl = dubl - bm*dls[idx3d];
    float dfr = dubr - bp*drs[idx3d];
    float ufl = dubl*(uls[idx3d] - bm) + pls[idx3d];
    float ufr = dubr*(urs[idx3d] - bp) + prs[idx3d];
    float vfl = dls[idx3d]*vls[idx3d]*(uls[idx3d] - bm);
    float vfr = drs[idx3d]*vrs[idx3d]*(urs[idx3d] - bp);
    float wfl = dls[idx3d]*wls[idx3d]*(uls[idx3d] - bm);
    float wfr = drs[idx3d]*wrs[idx3d]*(urs[idx3d] - bp);
    float efl = el*(uls[idx3d] - bm) + pls[idx3d]*uls[idx3d];
    float efr = er*(urs[idx3d] - bp) + prs[idx3d]*urs[idx3d];
    //
    // Gas energy is treated like an advected quantity
    //
    float gefl, gefr, gesl, gesr;
    if (idual == 1) {
      gefl = (uls[idx3d] - bm)*gels[idx3d]*dls[idx3d];
      gefr = (urs[idx3d] - bp)*gers[idx3d]*drs[idx3d];
      gesl = uls[idx3d] - bm;
      gesr = urs[idx3d] - bp;
    }
    //
    // Compute HLLC flux at interface
    //
    df[idx3d] = sl*dfl + sr*dfr;
    uf[idx3d] = sl*ufl + sr*ufr;
    vf[idx3d] = sl*vfl + sr*vfr;
    wf[idx3d] = sl*wfl + sr*wfr;
    ef[idx3d] = sl*efl + sr*efr;
    
    if (idual == 1) 
      gef[idx3d] = sl*gefl + sr*gefr;
    
    for (int c = 0; c < ncolor; c++) {
      float colfl = (uls[idx3d]-bm) * colls[idx3d+c*size3d];
      float colfr = (urs[idx3d]-bp) * colrs[idx3d+c*size3d];
      colf[idx3d+c*size3d] = sl*colfl + sr*colfr;
    }

    //
    //  Add the weighted contribution of the flux along the contact for
    //  the x-velocity and energy
    //
    uf[idx3d] += sm*cp;
    ef[idx3d] += sm*cp*cw;
    //
    // Calculate source term for gas energy
    //
    if (idual == 1) {
      gesf[idx3d] = sl*gesl+sr*gesr;
    }
    //
    // Finally absorb dt/dx into fluxes.  Advected quantities (including
    // the gas energy) are only multiplied by dt.
    //
    float qc = dt/dx;
    df[idx3d] *= qc;
    uf[idx3d] *= qc;
    vf[idx3d] *= qc;
    wf[idx3d] *= qc;
    ef[idx3d] *= qc;
    
    if (idual == 1) 
      gef[idx3d] = dt*gef[idx3d];

    for (int c = 0; c < ncolor; c++)
      colf[idx3d+c*size3d] *= dt;
  }
        //  }
  //  }
//  }

}

__global__ void flux_hll_kernel(float *dslice, float *eslice, float *geslice, 
                              float *uslice, float *vslice, float *wslice,
                              float dx, float *diffcoef, 
                              int idim, int jdim, int kdim,
                              int i1, int i2, int j1, int j2, 
                              int k1, int k2,
                              float dt, float gamma, int idiff,
                              int idual, float eta1, int ifallback,
                              float *dls, float *drs, 
                              float *pls, float *prs,
                              float *gels, float *gers,
                              float *uls, float *urs,
                              float *vls, float *vrs,
                              float *wls, float *wrs,
                              float *df, float *ef, float *uf,
                              float *vf, float *wf, float *gef, float *ges, float *gesf,
                              int ncolor, float *colslice, 
                              float *colls, float *colrs, float *colf)
{
  const int size3d = idim*jdim*kdim;
  int idx3d = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x + i1;
  float gamma1 = gamma - 1.0f;
  float gamma1i = 1.0f / gamma1;

  // for (int k = k1; k <= k2; k++) {
  //   for (int j = j1; j <= j2; j++) {
  //     for (int i = i1; i <= i2+1; i++) {
  //       int idx3d = (k*jdim+j)*idim + i;
  // int i = idx3d % idim;
  // int j = (idx3d % (idim*jdim)) / idim;
  // int k = idx3d / (idim*jdim);
  if (idx3d <= size3d - 3) {
    //
    // Compute the Roe-averaged data from the left and right states
    //
    float sqrtdl = sqrt(dls[idx3d]);
    float sqrtdr = sqrt(drs[idx3d]);
    float isdlpdr = 1.0f / (sqrtdl + sqrtdr);
    float vroe1 = (sqrtdl*uls[idx3d] + sqrtdr*urs[idx3d])*isdlpdr;
    float vroe2 = (sqrtdl*vls[idx3d] + sqrtdr*vrs[idx3d])*isdlpdr;
    float vroe3 = (sqrtdl*wls[idx3d] + sqrtdr*wrs[idx3d])*isdlpdr;
    float v2 = vroe1*vroe1 + vroe2*vroe2 + vroe3*vroe3;
    //
    // The enthalpy H=(E+P)/d is averaged for adiabatic flows, rather
    // than E or P directly.
    // sqrtql*hl = sqrtdl*(el+pl)/dl = (el+pl)/sqrtdl
    //
    float el = gamma1i*pls[idx3d] + 0.5f*dls[idx3d]*
      (uls[idx3d]*uls[idx3d]+vls[idx3d]*vls[idx3d]+wls[idx3d]*wls[idx3d]);
    float er = gamma1i*prs[idx3d] + 0.5f*drs[idx3d]*
      (urs[idx3d]*urs[idx3d]+vrs[idx3d]*vrs[idx3d]+wrs[idx3d]*wrs[idx3d]);
    float hroe = ((el+pls[idx3d])/sqrtdl + (er+prs[idx3d])/sqrtdr)*isdlpdr;
    //
    // Compute characteristics using Roe-averaged values
    //
    float cs = sqrt(gamma1*max(hroe-0.5f*v2, tiny));
    float char1 = vroe1 - cs;
    float char2 = vroe1 + cs;
    //
    // Compute the min/max wave speeds
    //
    float csl0 = sqrt(gamma*pls[idx3d]/dls[idx3d]);
    float csr0 = sqrt(gamma*prs[idx3d]/drs[idx3d]);
    float csl = min(uls[idx3d]-csl0, char1);
    float csr = max(urs[idx3d]+csr0, char2);
    float bm = min(csl, 0.0f);
    float bp = max(csr, 0.0f);
    float bm0 = uls[idx3d] - bm;
    float bp0 = urs[idx3d] - bp;
    //
    // Weights for different (left, right, star) regions in the cell
    //
    float q1 = (bp + bm) / (bp - bm);
    float sl = 0.5f*(1.0f + q1);
    float sr = 0.5f*(1.0f - q1);
    //
    // Compute diffusion terms
    //
    
    // Diffusion to be implemented
    
    // 
    //    Compute the left and right fluxes along the characteristics.
    //     Compute the color fluxes afterwards (loop n, then i) for cache
    //     efficiency.  For dual energy formalism, do not include the source
    //     term.  Treat like an advected quantity.
    //
    float dubl = dls[idx3d]*uls[idx3d];
    float dubr = drs[idx3d]*urs[idx3d];
    float dfl = dls[idx3d]*bm0;
    float dfr = drs[idx3d]*bp0;
    float ufl = dubl*bm0 + pls[idx3d];
    float ufr = dubr*bp0 + prs[idx3d];
    float vfl = dls[idx3d]*vls[idx3d]*bm0;
    float vfr = drs[idx3d]*vrs[idx3d]*bp0;
    float wfl = dls[idx3d]*wls[idx3d]*bm0;
    float wfr = drs[idx3d]*wrs[idx3d]*bp0;
    float efl = el*bm0 + pls[idx3d]*uls[idx3d];
    float efr = er*bp0 + prs[idx3d]*urs[idx3d];
    //
    // Gas energy is treated like an advected quantity
    //
    float gefl, gefr, gesl, gesr;
    if (idual == 1) {
      gefl = bm0*gels[idx3d]*dls[idx3d];
      gefr = bp0*gers[idx3d]*drs[idx3d];
      gesl = bm0;
      gesr = bp0;
    }
    //
    // Compute HLL flux at interface
    //
    df[idx3d] = sl*dfl + sr*dfr;
    uf[idx3d] = sl*ufl + sr*ufr;
    vf[idx3d] = sl*vfl + sr*vfr;
    wf[idx3d] = sl*wfl + sr*wfr;
    ef[idx3d] = sl*efl + sr*efr;
    
    if (idual == 1) 
      gef[idx3d] = sl*gefl + sr*gefr;
    
    for (int c = 0; c < ncolor; c++) {
      float colfl = bm0 * colls[idx3d+c*size3d];
      float colfr = bp0 * colrs[idx3d+c*size3d];
      colf[idx3d+c*size3d] = sl*colfl + sr*colfr;
    }
    //
    // Calculate source term for gas energy
    //
    if (idual == 1) {
      gesf[idx3d] = sl*gesl+sr*gesr;
    }
    //
    // Finally absorb dt/dx into fluxes.  Advected quantities (including
    // the gas energy) are only multiplied by dt.
    //
    float qc = dt/dx;
    df[idx3d] *= qc;
    uf[idx3d] *= qc;
    vf[idx3d] *= qc;
    wf[idx3d] *= qc;
    ef[idx3d] *= qc;
    
    if (idual == 1) 
      gef[idx3d] = dt*gef[idx3d];

    for (int c = 0; c < ncolor; c++)
      colf[idx3d+c*size3d] *= dt;
  }
        //  }
  //  }
//  }

}


__global__ void check_negative_states(float *dslice, float *eslice, float *df, 
                                      int ifallback, int *fallback,
                                      int idim, int jdim, int kdim,
                                      int i1, int i2, int j1, int j2, int k1, int k2)
{
  const int size3d = idim*jdim*kdim;
  int idx3d = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x + i1;
  int i = idx3d % idim; 
  if (idx3d <= size3d - 4 && i >= i1 && i <= i2) {
    if (dslice[idx3d] + (df[idx3d] - df[idx3d+1]) <= 0.0f ||
        eslice[idx3d] < 0.0f) {
      /*      if (eslice[idx3d] < 0.0f)
        printf("flux_twoshock: eslice < 0\n");
      else
        printf("flux_twoshock: dnu <= 0\n");
      */
      fallback[0] = 1;
    }
  }
}
                                      
__global__ void calc_ges(float *ges, float *gesf,
                         float *geslice, float *dslice,
                         float gamma, float dx, float dt,
                         int idim, int jdim, int kdim,
                         int i1, int i2, int j1, int j2, int k1, int k2)
{
  const int size3d = idim*jdim*kdim;
  int idx3d = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x + i1;
  int i = idx3d % idim;
  // int j = (idx3d % (idim*jdim)) / idim;
  // int k = idx3d / (idim*jdim);
  
  if (idx3d <= size3d - 4 && i >= i1 && i <= i2) {
    float pcent = max((gamma-1.0f)*geslice[idx3d]*dslice[idx3d], tiny);
    float qc = dt/dx;
    ges[idx3d] = qc*pcent*(gesf[idx3d] - gesf[idx3d+1]);
  }
}
                         
__global__ void pgas2d_dual_kernel(float *eslice, float *pslice, 
                                 float *dslice, float *geslice,
                                 float *uslice, float *vslice, float *wslice,
                                 float eta1, float eta2, 
                                 float gamma, float pmin,
                                 int idim, int jdim, int kdim,
                                 int i1, int i2, int j1, int j2, int k1, int k2)
{
  const int size3d = idim*jdim*kdim;
  int idx3d = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x + i1;
  int i = idx3d % idim;
  if (idx3d < size3d) {
    //
    // Compute the specific energy from total energy
    //
    float ke = 0.5f*(uslice[idx3d]*uslice[idx3d]+vslice[idx3d]*vslice[idx3d]+
                     wslice[idx3d]*wslice[idx3d]);
    float ge1 = eslice[idx3d] - ke;
    //
    // Find the maximum nearby total energy (not specific)
    //
    int im1 = (i-1 >= i1) ? -1 : 0;
    int ip1 = (i+1 <= i2) ? +1 : 0;
    float demax = max3(dslice[idx3d    ]*eslice[idx3d    ],
                       dslice[idx3d+im1]*eslice[idx3d+im1],
                       dslice[idx3d+ip1]*eslice[idx3d+ip1]);
    if (ge1*dslice[idx3d]/demax > eta2) geslice[idx3d] = ge1;
    //
    // If the ratio of the specific gas energy to specific total energy
    //  (in this cell) is < eta1, then use the gas energy to update the 
    //  specific energy (total energy is better).
    //
    float ge2;
    if (ge1/eslice[idx3d] > eta1) 
      ge2 = ge1;
    else
      ge2 = geslice[idx3d];
    //
    // If pressure is below the minimum, set it to the minimum
    //
    ge2 = max(ge2, pmin/((gamma-1.0f)*dslice[idx3d]));
    //
    // Update the total energy
    //
    eslice[idx3d] = eslice[idx3d] - ge1 + ge2;
    //
    // Compute the pressure with the total energy derived gas energy
    //
    //pslice[idx3d] = (gamma-1.0f)*dslice[idx3d]*ge2;
  }
}
__global__ void euler_kernel(float *dslice, float *eslice, float *grslice, 
                           float *geslice, float *uslice, float *vslice,
                           float *wslice, float dx,
                           int idim, int jdim, int kdim,
                           int i1, int i2, int j1, int j2, int k1, int k2,
                           float dt, float gamma, int gravity,
                           int idual, float eta1, float et2,
                           float *df, float *ef, float *uf, float *vf,
                           float *wf, float *gef, float *ges,
                           int ncolor, float *colslice, float *colf)
{
  float eold;
  float dnu, dnuinv, uold;
  const float min_color = 1e-5*tiny;

  const int size3d = idim*jdim*kdim;
  int idx3d = blockIdx.y*gridDim.x*blockDim.x + blockIdx.x*blockDim.x + threadIdx.x + i1;

  // for (int k = k1; k <= k2; k++) {
  //   for (int j = j1; j <= j2; j++) {
  //     for (int i = i1; i <= i2; i++) {
  //       int idx3d = (k*jdim+j)*idim + i;
  int i = idx3d % idim;
  // int j = (idx3d % (idim*jdim)) / idim;
  // int k = idx3d / (idim*jdim);
  if (idx3d <= size3d - 4 && i >= i1 && i <= i2) {
  // if (blockIdx.x == 0 && threadIdx.x == 0) {
  //    for (int k = k1; k <= k2; k++) {
  //    for (int j = j1; j <= j2; j++) {
  //      for (int i = i1; i <= i2; i++) {
  //        int idx3d = (k*jdim+j)*idim + i;
    
    dnu = dslice[idx3d] + (df[idx3d] - df[idx3d+1]);
    dnuinv = 1.0f/dnu;
    //
    // A small cheat: if idual is on, assume a cosmo sim, and institute a
    // minimum density. This should be a parameter and passed in.
    //
    if (idual == 1 && gravity == 1)
      dnu = max(dnu, 1e-3);
        
    uold  = uslice[idx3d];
    uslice[idx3d] = (uslice[idx3d]*dslice[idx3d] + 
                     (uf[idx3d] - uf[idx3d+1]))*dnuinv;
    vslice[idx3d] = (vslice[idx3d]*dslice[idx3d] +
                     (vf[idx3d] - vf[idx3d+1]))*dnuinv;
    wslice[idx3d] = (wslice[idx3d]*dslice[idx3d] + 
                     (wf[idx3d] - wf[idx3d+1]))*dnuinv;
    eold = eslice[idx3d];
    // eslice[idx3d] = max(0.1f*eslice[idx3d],
    //                     (eslice[idx3d]*dslice[idx3d] +
    //                      (ef[idx3d]-ef[idx3d+1]))*dnuinv); 
    eslice[idx3d] = max(0.1f*eold, (eold*dslice[idx3d]+(ef[idx3d]-ef[idx3d+1]))*dnuinv);

    //
    // Colour variables (note: colf already multiplied by dt)
    //
    for (int c = 0; c < ncolor; c++) {
      int colidx = idx3d + c*size3d;
      colslice[colidx] = colslice[colidx] + 
        (colf[colidx] - colf[colidx+1])/dx;
      colslice[colidx] = max(colslice[colidx], min_color);
    }
    //
    // Conservation law for gas energy, if using the dual energy formalism
    //   (this includes both the flux term and a source term -yuck).
    //   Here, we comptue the ratio of thermal energies derived the 
    //     two different ways and then use that ratio to choose how to 
    //     compute the pressure a the center of zone i. This is needed
    //     for the source term in the gas energy equation.
    //
    if (idual == 1) {
      geslice[idx3d] = max((geslice[idx3d]*dslice[idx3d] +
                            (gef[idx3d] - gef[idx3d+1])/dx + ges[idx3d])*dnuinv,
                           0.5f*geslice[idx3d]);
    }
    //
    // If there is gravity, then compute the second order correction to the
    //   acceleration due to a slope in both density and acceleration.
    //
#ifdef GRAVITY_SECOND_ORDER_CORRECTION
    //
    // Compute slopes and enforce limited monotonocity on dddx
    //
    if (gravity == 1) {
      dadx = grslice[idx3d+1] - grslice[idx3d-1];
      dddx =  dslice[idx3d+1] -  dslice[idx3d-1];
      dddx = 2.0f*( dslice[idx3d] - max(dslice[idx3d] - 0.5f*dddx,
                                        min(dslice[idx3d], dslice[idx3d-1])));
      dddx = 2.0f*(-dslice[idx3d] + max(dslice[idx3d] + 0.5f*dddx,
                                        min(dslice[idx3d], dslice[idx3d+1])));
      grslice[idx3d] = grslice[idx3d] + 
        0.5f*dadx*dddx/(12.0f*dslice[idx3d]);
    }
#endif
    //
    // If there is gravity, add the gravity terms here (eq. 3.1 or 3.8)
    //   (Note: the acceleration is already time-centered)
    //
    if (gravity == 1) {
#define GRAVITY_METHOD1
#ifdef GRAVITY_METHOD1
      uslice[idx3d] = uslice[idx3d] +
        dt*grslice[idx3d]*0.5f*(dslice[idx3d]*dnuinv + 1.0f);
      eslice[idx3d] = eslice[idx3d] +
        dt*grslice[idx3d]*0.5f*(uslice[idx3d]+uold*dslice[idx3d]*dnuinv);
      eslice[idx3d] = max(eslice[idx3d], tiny);
    }
#endif  
    //
    // Update the new density
    //
    dslice[idx3d] = dnu;
  }
  //           }
  //    }
  // }

}


void PPM_CUDA(cuPPMData &Data,
              cuPPMParameter &Para,
              int idim, int jdim, int kdim,
              int i1, int i2, int j1, int j2, int k1, int k2,
              float dx, float dt)
{
  const int size = idim*jdim*kdim;
  dim3 block;
  dim3 grid;

  block.x = NBLOCK, block.y = 1, block.z = 1;

  if (Para.PPMFlatteningParameter) {
    grid.x = ((size - 4) + block.x - 1)/block.x, grid.y = 1, grid.z = 1;
    if (grid.x > 65535)
      grid.y = (grid.x+255)/256, grid.x = 256;
    
    calc_flatten_kernel<<<grid,block>>>(Data.uslice, Data.pslice,
                                      idim, jdim, kdim,
                                      i1, i2, j1, j2, k1, k2,
                                      Data.flatten);
//    CUDA_SAFE_CALL( cudaGetLastError() );
  }

  grid.x = ((size - 5) + NBLOCK - 1)/NBLOCK, grid.y = 1, grid.z = 1;
  if (grid.x > 65535)
    grid.y = (grid.x+255)/256, grid.x = 256;

  inteuler_kernel<<<grid,block>>>(
    Data.dslice, Data.pslice, Para.GravityOn,
    Data.grslice, Data.geslice,
    Data.uslice, Data.vslice, Data.wslice,
    dx, Data.flatten, 
    idim, jdim, kdim,
    i1, i2, j1, j2, k1, k2,
    Para.DualEnergyFormalism, 
    Para.DualEnergyFormalismEta1,
    Para.DualEnergyFormalismEta2,
    Para.PPMSteepeningParameter,
    Para.PPMFlatteningParameter,
    Para.ConservativeReconstruction,
    Para.PositiveReconstruction,
    dt, Para.Gamma, Para.PressureFree,
    Data.dls, Data.drs, Data.pls, Data.prs,
    Data.gels, Data.gers, Data.uls, Data.urs,
    Data.vls, Data.vrs, Data.wls, Data.wrs,
    Para.NumberOfColours, Data.colslice, 
    Data.colls, Data.colrs);
//  CUDA_SAFE_CALL( cudaGetLastError() );

  int fallback = 0;
  switch (Para.RiemannSolver) {
  case TwoShock:
    grid.x = ((size - 5) + NBLOCK - 1)/NBLOCK, grid.y = 1, grid.z = 1;
    if (grid.x > 65535)
      grid.y = (grid.x+255)/256, grid.x = 256;
    twoshock_kernel<<<grid,block>>>(Data.dls, Data.drs,
                                  Data.pls, Data.prs,
                                  Data.uls, Data.urs,
                                  idim, jdim, kdim,
                                  i1, i2, j1, j2, k1, k2,
                                  dt, Para.Gamma, Para.MinimumPressure, 
                                  Para.PressureFree,
                                  Data.pbar, Data.ubar, 
                                  Para.GravityOn, Data.grslice,
                                  Para.DualEnergyFormalism, 
                                  Para.DualEnergyFormalismEta1);
    CUDA_SAFE_CALL( cudaGetLastError() );

    flux_twoshock_kernel<<<grid, block>>>(
      Data.dslice, Data.eslice, Data.geslice, Data.uslice, 
      Data.vslice, Data.wslice, dx, Data.diffcoef, 
      idim, jdim, kdim,
      i1, i2, j1, j2, k1, k2,
      dt, Para.Gamma, Para.PPMDiffusionParameter,
      Para.DualEnergyFormalism, Para.DualEnergyFormalismEta1,
      Para.RiemannSolverFallback,
      Data.dls, Data.drs, Data.pls, Data.prs,
      Data.gels, Data.gers, Data.uls, Data.urs,
      Data.vls, Data.vrs, Data.wls, Data.wrs,
      Data.pbar, Data.ubar,
      Data.df, Data.ef, Data.uf, Data.vf, 
      Data.wf, Data.gef, Data.ges, Data.gesf,
      Para.NumberOfColours, Data.colslice, 
      Data.colls, Data.colrs, Data.colf);
    CUDA_SAFE_CALL( cudaGetLastError() );

    grid.x = ((size - 6) + NBLOCK - 1)/NBLOCK;
    if (grid.x > 65535)
      grid.y = (grid.x+255)/256, grid.x = 256;
    check_negative_states<<<grid, block>>>
      (Data.dslice, Data.eslice, Data.df, 
       Para.RiemannSolverFallback, Data.fallback,
       idim, jdim, kdim, i1, i2, j1, j2, k1, k2);
    cudaMemcpy(&fallback, Data.fallback, sizeof(int), cudaMemcpyDeviceToHost);
    if (fallback) {
      if (Para.RiemannSolverFallback) {
        printf("falling back to HLL solver\n");
        cudaMemset(Data.fallback, 0, sizeof(int));
        goto LABEL_HLL;
      } else {
        printf("Two shock solver failed\n");
        exit(1);
      }
    }
    break;
  case HLLC:
    grid.x = ((size - 5) + NBLOCK - 1)/NBLOCK, grid.y = 1, grid.z = 1;
    if (grid.x > 65535)
      grid.y = (grid.x+255)/256, grid.x = 256;
    flux_hllc_kernel<<<grid, block>>>(
      Data.dslice, Data.eslice, Data.geslice, Data.uslice, 
      Data.vslice, Data.wslice, dx, Data.diffcoef, 
      idim, jdim, kdim,
      i1, i2, j1, j2, k1, k2,
      dt, Para.Gamma, Para.PPMDiffusionParameter,
      Para.DualEnergyFormalism, Para.DualEnergyFormalismEta1,
      Para.RiemannSolverFallback,
      Data.dls, Data.drs, Data.pls, Data.prs,
      Data.gels, Data.gers, Data.uls, Data.urs,
      Data.vls, Data.vrs, Data.wls, Data.wrs,
      Data.df, Data.ef, Data.uf, Data.vf, 
      Data.wf, Data.gef, Data.ges, Data.gesf,
      Para.NumberOfColours, Data.colslice, 
      Data.colls, Data.colrs, Data.colf);
    CUDA_SAFE_CALL( cudaGetLastError() );

    grid.x = ((size - 6) + NBLOCK - 1)/NBLOCK;
    if (grid.x > 65535)
      grid.y = (grid.x+255)/256, grid.x = 256;
    check_negative_states<<<grid, block>>>
      (Data.dslice, Data.eslice, Data.df, 
       Para.RiemannSolverFallback, Data.fallback,
       idim, jdim, kdim, i1, i2, j1, j2, k1, k2);
    cudaMemcpy(&fallback, Data.fallback, sizeof(int), cudaMemcpyDeviceToHost);
    if (fallback) {
      if (Para.RiemannSolverFallback) {
        printf("falling back to HLL solver\n");
        cudaMemset(Data.fallback, 0, sizeof(int));
        goto LABEL_HLL;
      } else {
        printf("HLLC solver failed\n");
      }
    }
    break;
  case HLL:
  LABEL_HLL:
    grid.x = ((size - 5) + NBLOCK - 1)/NBLOCK, grid.y = 1, grid.z = 1;
    if (grid.x > 65535)
      grid.y = (grid.x+255)/256, grid.x = 256;
    flux_hll_kernel<<<grid, block>>>(
      Data.dslice, Data.eslice, Data.geslice, Data.uslice, 
      Data.vslice, Data.wslice, dx, Data.diffcoef, 
      idim, jdim, kdim,
      i1, i2, j1, j2, k1, k2,
      dt, Para.Gamma, Para.PPMDiffusionParameter,
      Para.DualEnergyFormalism, Para.DualEnergyFormalismEta1,
      Para.RiemannSolverFallback,
      Data.dls, Data.drs, Data.pls, Data.prs,
      Data.gels, Data.gers, Data.uls, Data.urs,
      Data.vls, Data.vrs, Data.wls, Data.wrs,
      Data.df, Data.ef, Data.uf, Data.vf, 
      Data.wf, Data.gef, Data.ges, Data.gesf,
      Para.NumberOfColours, Data.colslice, 
      Data.colls, Data.colrs, Data.colf);
//    CUDA_SAFE_CALL( cudaGetLastError() );
    grid.x = ((size - 6) + NBLOCK - 1)/NBLOCK;
    if (grid.x > 65535)
      grid.y = (grid.x+255)/256, grid.x = 256;
    check_negative_states<<<grid, block>>>
      (Data.dslice, Data.eslice, Data.df, 
       Para.RiemannSolverFallback, Data.fallback,
       idim, jdim, kdim, i1, i2, j1, j2, k1, k2);
    cudaMemcpy(&fallback, Data.fallback, sizeof(int), cudaMemcpyDeviceToHost);
    if (fallback) {
      if (Para.RiemannSolverFallback) {
        printf("Already using HLL. No solver to fall back. Give up.\n");
        exit(1);
      } else {
        printf("HLL solver failed\n");
        exit(1);
      }
    }
    break;
  default:
    printf("Unsupported Riemann Solver on GPU: %d\n", Para.RiemannSolver);
    exit(1);
  }

  if (Para.DualEnergyFormalism) {
    grid.x = ((size - 6) + NBLOCK - 1)/NBLOCK;
    if (grid.x > 65535)
      grid.y = (grid.x+255)/256, grid.x = 256;

    calc_ges<<<grid, block>>>(Data.ges, Data.gesf, Data.geslice, Data.dslice,
                              Para.Gamma, dx, dt, idim, jdim, kdim, 
                              i1, i2, j1, j2, k1, k2);
  }

  euler_kernel<<<grid, block>>>(
    Data.dslice, Data.eslice, Data.grslice, Data.geslice,
    Data.uslice, Data.vslice, Data.wslice,
    dx, idim, jdim, kdim,
    i1, i2, j1, j2, k1, k2, dt, Para.Gamma,
    Para.GravityOn, Para.DualEnergyFormalism,
    Para.DualEnergyFormalismEta1, Para.DualEnergyFormalismEta2,
    Data.df, Data.ef, Data.uf, Data.vf, Data.wf,
    Data.gef, Data.ges, 
    Para.NumberOfColours, Data.colslice, Data.colf);
//  CUDA_SAFE_CALL( cudaGetLastError() );

  if (Para.DualEnergyFormalism) {
    grid.x = (size + NBLOCK - 1)/NBLOCK;
    if (grid.x > 65535)
      grid.y = (grid.x+255)/256, grid.x = 256;
    pgas2d_dual_kernel<<<grid, block>>>
      (Data.eslice, Data.pslice,
       Data.dslice, Data.geslice,
       Data.uslice, Data.vslice, Data.wslice,
       Para.DualEnergyFormalismEta1,
       Para.DualEnergyFormalismEta2,
       Para.Gamma, Para.MinimumPressure,
       idim, jdim, kdim,
       i1-3, i2+3, j1, j2, k1, k2);
  }

}

