#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"


// initializatiion function related to HEALPix
int mkPix2xy(long *ipix2x, long *ipix2y) 
{
  long  kpix,jpix,ix,iy,ip,id;

  for(kpix=0; kpix <= 1023; kpix++) {          
    jpix = kpix;
    ix = 0;
    iy = 0;
    ip = 1 ;            
    while (jpix != 0) { 
      id = jpix&1;
      jpix = jpix/2 ;
      ix = id*ip+ix ;
      id = jpix&1;
      jpix = jpix/2;
      iy = id*ip+iy;
      ip = 2*ip;  
    }  
    *ipix2x++ = ix;   
    *ipix2y++ = iy ;  
  }
  return SUCCESS;
}

/***********************************************************************/

int mk_xy2pix(int *x2pix, int *y2pix) {
  /* =======================================================================
   * subroutine mk_xy2pix
   * =======================================================================
   * sets the array giving the number of the pixel lying in (x,y)
   * x and y are in {1,128}
   * the pixel number is in {0,128**2-1}
   *
   * if  i-1 = sum_p=0  b_p * 2^p
   * then ix = sum_p=0  b_p * 4^p
   * iy = 2*ix
   * ix + iy in {0, 128**2 -1}
   * =======================================================================
   */
  int i, K,IP,I,J,ID;
  
  for(i = 0; i < 127; i++) x2pix[i] = 0;
  for( I=1;I<=128;I++ ) {
    J  = I-1; // pixel numbers
    K  = 0;
    IP = 1;
    while(J!=0) {
      ID = J&1;
      J  = J/2;
      K  = IP*ID+K;
      IP = IP*4;
    }
    x2pix[I-1] = K;
    y2pix[I-1] = 2*K;
  }
  return SUCCESS;
}

#include "typedefs.h"
#include "global_data.h"

/***********************************************************************/

void sub_compute_vertices(double z, double z_nv, double z_sv, double phi, double phi_nv, double phi_sv, double hdelta_phi, double (*vertex)[3]) {
    double sth = sqrt((1.0-z)*(1.0+z));

    double sth_nv = sqrt((1.0-z_nv)*(1.0+z_nv));
    vertex[0][0] = sth_nv*cos(phi_nv);
    vertex[0][1] = sth_nv*sin(phi_nv);
    vertex[0][2] = z_nv; // north vertex

    double phi_wv = phi - hdelta_phi;
    vertex[1][0] = sth*cos(phi_wv);
    vertex[1][1] = sth*sin(phi_wv);
    vertex[1][2] = z; // west vertex

    double sth_sv = sqrt((1.0-z_sv)*(1.0+z_sv));
    vertex[2][0] = sth_sv*cos(phi_sv); // south vertex
    vertex[2][1] = sth_sv*sin(phi_sv); // south vertex
    vertex[2][2] = z_sv; // south vertex

    double phi_ev = phi + hdelta_phi;
    vertex[3][0] = sth*cos(phi_ev); // east vertex
    vertex[3][1] = sth*sin(phi_ev); // east vertex
    vertex[3][2] = z; // east vertex

    return;
}


// HEALPix fuction that updates directional vector given a pixel number
#define nsMax   8192L  
#define halfpi  1.5707963
int pix2vec_nest(long nside, long ipix, FLOAT *v, double (*vertex)[3]=0)
{ 

  int npix, npface, face_num;
  int  ipf, ip_low, ip_trunc, ip_med, ip_hi;
  int     ix, iy, jrt, jr, nr, jpt, jp, kshift, nl4;
  double z, fn, fact1, fact2;
  double piover2=0.5*M_PI;
  double phi, sz, vx, vy, vz;
  int ns_max=8192;

  //static int pix2x[1024], pix2y[1024];
  //      common /pix2xy/ pix2x, pix2y

  int jrll[12]={2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4}; // in unit of nside
  int jpll[12]={1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7}; // coordinate of the lowest corner of each face in unit of nside/2

  if( nside<1 || nside>ns_max ) {
    ENZO_VFAIL("%s (%"ISYM"): nside out of range: %ld\n", __FILE__, __LINE__, nside)
  }
  npix = 12 * nside*nside;
  if( ipix<0 || ipix>npix-1 ) {
    ENZO_VFAIL("%s (%"ISYM"): ipix out of range: %ld\n", __FILE__, __LINE__, ipix)
  }

  /* initiates the array for the pixel number -> (x,y) mapping */
  //if( ipix2x[1023]<=0 ) mkPix2xy(ipix2x,ipix2y);

  fn = 1.*nside;
  fact1 = 1./(3.*fn*fn);
  fact2 = 2./(3.*fn);
  nl4   = 4*nside;

  //c     finds the face, and the number in the face
  npface = nside*nside;

  face_num = ipix/npface;//  ! face number in {0,11}
  ipf = ipix % npface;//  ! pixel number in the face {0,npface-1}

  //c     finds the x,y on the face (starting from the lowest corner)
  //c     from the pixel number
  ip_low = ipf & 0x3FF; //       ! content of the last 10 bits
  ip_trunc =   ipf>>10; //       ! truncation of the last 10 bits
  ip_med = ip_trunc & 0x3FF; //  ! content of the next 10 bits
  ip_hi  =     ip_trunc>>10;   //! content of the high weight 10 bits

  ix = 1024*pix2x[ip_hi] + 32*pix2x[ip_med] + pix2x[ip_low];
  iy = 1024*pix2y[ip_hi] + 32*pix2y[ip_med] + pix2y[ip_low];

  //c     transforms this in (horizontal, vertical) coordinates
  jrt = ix + iy;//  ! 'vertical' in {0,2*(nside-1)}
  jpt = ix - iy;//  ! 'horizontal' in {-nside+1,nside-1}

  //c     computes the z coordinate on the sphere
  //      jr =  jrll[face_num+1]*nside - jrt - 1;//   ! ring number in {1,4*nside-1}
  jr =  jrll[face_num]*nside - jrt - 1;
  //      cout << "face_num=" << face_num << endl;
  //      cout << "jr = " << jr << endl;
  //      cout << "jrll(face_num)=" << jrll[face_num] << endl;
  //      cout << "----------------------------------------------------" << endl;
  nr = nside;//                  ! equatorial region (the most frequent)
  z  = (2*nside-jr)*fact2;
  kshift = (jr - nside)&1;
  if( jr<nside ) { //then     ! north pole region
    nr = jr;
    z = 1. - nr*nr*fact1;
    kshift = 0;
  }
  else {
    if( jr>3*nside ) {// then ! south pole region
      nr = nl4 - jr;
      z = - 1. + nr*nr*fact1;
      kshift = 0;
    }
  }

  //c     computes the phi coordinate on the sphere, in [0,2Pi]
  //      jp = (jpll[face_num+1]*nr + jpt + 1 + kshift)/2;//  ! 'phi' number in the ring in {1,4*nr}
  jp = (jpll[face_num]*nr + jpt + 1 + kshift)/2;
  if( jp>nl4 ) jp = jp - nl4;
  if( jp<1 )   jp = jp + nl4;

  double delta_phi = piover2 / nr;
  phi = (jp - (kshift+1)*0.5) * delta_phi;

  if (vertex != 0) {
      double iphi_mod = (jp-1) % nr; // in {0,1,... nr-1}
      double iphi_rat = (jp-1) / nr;      // in {0,1,2,3}
      double phi_up = piover2 * (iphi_rat +  iphi_mod   /double(max(nr-1,1)));
      double phi_dn = piover2 * (iphi_rat + (iphi_mod+1)/double(nr+1));
      double z_nv, z_sv;
      iphi_rat = 0;
      iphi_mod = 0;

      if( jr<nside ) { //then     ! north pole region
          z_nv = 1.0 - (nr-1.)*(nr-1.)*fact1;
          z_sv = 1.0 - (nr+1.)*(nr+1.)*fact1;
          sub_compute_vertices(z, z_nv, z_sv, phi, phi_up, phi_dn, delta_phi/2, vertex);
      } else if( jr>3*nside ) { // then ! south pole region
          z_nv = -1.0 + (nr+1.)*(nr+1.)*fact1;
          z_sv = -1.0 + (nr-1.)*(nr-1.)*fact1;
          sub_compute_vertices(z, z_nv, z_sv, phi, phi_dn, phi_up, delta_phi/2, vertex);
      } else {
          // equatorial region
          z_nv = (2*nside-jr+1)*fact2;
          z_sv = (2*nside-jr-1)*fact2;
          if (jr == nside) { // northern transition
              z_nv =  1.0 - (nside-1.)*(nside-1.) * fact1;
              sub_compute_vertices(z, z_nv, z_sv, phi, phi_up, phi, delta_phi/2, vertex);
          } else if (jr == 3*nside) {
              z_sv = -1.0 + (nside-1.)*(nside-1.) * fact1;
              sub_compute_vertices(z, z_nv, z_sv, phi, phi, phi_up, delta_phi/2, vertex);
          } else
              sub_compute_vertices(z, z_nv, z_sv, phi, phi, phi, delta_phi/2, vertex);
      }
  }

  sz = sqrt(1.0-z*z);
  vx = sz * cos(phi);
  vy = sz * sin(phi);
  vz = z;

  // Rotate the vector by 45 degrees in the x-axis
//  float xrotation = 0;  // degrees
//  float yrotation = 0;
//  xrotation *= M_PI/180.0;
//  yrotation *= M_PI/180.0;
//
//  v[0] = vx;
//  v[1] = vy * cos(xrotation) - vz * sin(xrotation);
//  v[2] = vy * sin(xrotation) + vz * cos(xrotation);
//
//  vx = v[0] * cos(yrotation) + v[2] * sin(yrotation);
//  vy = v[1];
//  vz = -v[0] * sin(yrotation) + v[2] * cos(yrotation);

  v[0] = vx;
  v[1] = vy;
  v[2] = vz;

  return SUCCESS;
}

/***********************************************************************/

int vec2pix_nest( const long nside, FLOAT *vec, long *ipix) {

  /* =======================================================================
   * subroutine vec2pix_nest(nside, vec, ipix)
   * =======================================================================
   * gives the pixel number ipix (NESTED) corresponding to vector vec
   *
   * the computation is made to the highest resolution available (nside=8192)
   * and then degraded to that required (by integer division)
   * this doesn't cost more, and it makes sure that the treatement of round-off 
   * will be consistent for every resolution
   * =======================================================================
   */

  double z, za, z0, tt, tp, tmp, phi;
  int    face_num,jp,jm;
  long   ifp, ifm;
  int    ix, iy, ix_low, ix_hi, iy_low, iy_hi, ipf, ntt;
  double piover2 = 0.5*M_PI, twopi = 2.0*M_PI;
  int    ns_max = 8192;
  //  static int x2pix[128], y2pix[128];
  //  static char setup_done = 0;

  if( nside<1 || nside>ns_max ) {
    ENZO_VFAIL("%s (%"ISYM"): nside out of range: %ld\n", 
	    __FILE__, __LINE__, nside)
  }
  //  if( !setup_done ) {
  //    mk_xy2pix(x2pix,y2pix);
  //    setup_done = 1;
  //  }

  z   = vec[2]/sqrt(vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
  phi = 0.0;
  if (vec[0] != 0.0 || vec[1] != 0.0) {
    phi   = atan2(vec[1],vec[0]); /* in ]-pi, pi] */
    if (phi < 0.0) phi += twopi; /* in  [0, 2pi[ */
  }

  za = fabs(z);
  z0 = 2./3.;
  tt = phi / piover2; /* in [0,4[ */

  if( za<=z0 ) { /* equatorial region */

    /* (the index of edge lines increase when the longitude=phi goes up) */
    jp = (int)floor(ns_max*(0.5 + tt - z*0.75)); /* ascending edge line index */
    jm = (int)floor(ns_max*(0.5 + tt + z*0.75)); /* descending edge line index */

    /* finds the face */
    ifp = jp / ns_max; /* in {0,4} */
    ifm = jm / ns_max;

    if( ifp==ifm ) face_num = (ifp&3) + 4; /* faces 4 to 7 */
    else if( ifp<ifm ) face_num = (ifp&3); /* (half-)faces 0 to 3 */
    else face_num = (ifm&3) + 8;           /* (half-)faces 8 to 11 */

    ix = jm % ns_max;
    iy = ns_max - jp % ns_max - 1;
  }
  else { /* polar region, za > 2/3 */

    ntt = (int)floor(tt);
    if( ntt>=4 ) ntt = 3;
    tp = tt - ntt;
    tmp = sqrt( 3.*(1. - za) ); /* in ]0,1] */

    /* (the index of edge lines increase when distance from the closest pole
     * goes up)
     */
    /* line going toward the pole as phi increases */
    jp = (int)floor( ns_max * tp          * tmp ); 

    /* that one goes away of the closest pole */
    jm = (int)floor( ns_max * (1. - tp) * tmp );
    jp = (int)(jp < ns_max-1 ? jp : ns_max-1);
    jm = (int)(jm < ns_max-1 ? jm : ns_max-1);

    /* finds the face and pixel's (x,y) */
    if( z>=0 ) {
      face_num = ntt; /* in {0,3} */
      ix = ns_max - jm - 1;
      iy = ns_max - jp - 1;
    }
    else {
      face_num = ntt + 8; /* in {8,11} */
      ix =  jp;
      iy =  jm;
    }
  }

  ix_low = ix&0x7F;
  ix_hi  =     ix>>7;
  iy_low = iy&0x7F;
  iy_hi  =     iy>>7;

  ipf = (x2pix[ix_hi]+y2pix[iy_hi]) * (128 * 128)+ (x2pix[ix_low]+y2pix[iy_low]);
  ipf = (long)(ipf / pow(ns_max/nside,2));     /* in {0, nside**2 - 1} */
  *ipix =(long)( ipf + face_num*pow(nside,2)); /* in {0, 12*nside**2 - 1} */

  return SUCCESS;

}

/**********************************************************************/

int pix2coord_nest( long nside, long ipix, int xsize, int ysize, 
		    int &npts, int * &xpolygon, int * &ypolygon, 
		    int &draw_poly) {

  /*
    =======================================================================
    subroutine pix2ang_nest(nside, ipix, theta, phi)
    =======================================================================
         gives theta and phi corresponding to pixel ipix (NESTED) 
         for a parameter nside
    =======================================================================
  */
    
  int npix, npface, face_num;
  int ipf, ip_low, ip_trunc, ip_med, ip_hi;
  int ix, iy, jrt, jr, nr, jpt, jp, kshift, nl4;
  double phi, z, fn, fact1, fact2;
  double piover2=0.5*M_PI;
  int ns_max=8192;

  // ! coordinate of the lowest corner of each face
  int jrll[12]={2,2,2,2,3,3,3,3,4,4,4,4}, jpll[12]={1,3,5,7,0,2,4,6,1,3,5,7};

  if( nside<1 || nside>ns_max ) {
    ENZO_VFAIL("%s (%"ISYM"): nside out of range: %ld\n", __FILE__, __LINE__, nside)
  }
  npix = 12 * nside*nside;
  if( ipix<0 || ipix>npix-1 ) {
    ENZO_VFAIL("%s (%"ISYM"): ipix out of range: %ld\n", __FILE__, __LINE__, ipix)
  }

  /* initiates the array for the pixel number -> (x,y) mapping */
  //  if( pix2x[1023]<=0 ) mkPix2xy(pix2x,pix2y);

  fn = 1.*nside;
  fact1 = 1./(3.*fn*fn);
  fact2 = 2./(3.*fn);
  nl4   = 4*nside;

  //c     finds the face, and the number in the face
  npface = nside*nside;

  face_num = ipix/npface;//  ! face number in {0,11}
  ipf = ipix % npface;//  ! pixel number in the face {0,npface-1}

  //c     finds the x,y on the face (starting from the lowest corner)
  //c     from the pixel number
  ip_low = ipf&0x3FF;//       ! content of the last 10 bits
  ip_trunc =   ipf>>10;//       ! truncation of the last 10 bits
  ip_med = ip_trunc&0x3FF;//  ! content of the next 10 bits
  ip_hi  =     ip_trunc>>10;//! content of the high weight 10 bits

  ix = 1024*pix2x[ip_hi] + 32*pix2x[ip_med] + pix2x[ip_low];
  iy = 1024*pix2y[ip_hi] + 32*pix2y[ip_med] + pix2y[ip_low];

  //c     transforms this in (horizontal, vertical) coordinates
  jrt = ix + iy;//  ! 'vertical' in {0,2*(nside-1)}
  jpt = ix - iy;//  ! 'horizontal' in {-nside+1,nside-1}

  //c     computes the z coordinate on the sphere
  //      jr =  jrll[face_num+1]*nside - jrt - 1;//   ! ring number in {1,4*nside-1}
  jr =  jrll[face_num]*nside - jrt - 1;
  //      cout << "face_num=" << face_num << endl;
  //      cout << "jr = " << jr << endl;
  //      cout << "jrll(face_num)=" << jrll[face_num] << endl;
  //      cout << "----------------------------------------------------" << endl;
  nr = nside;//                  ! equatorial region (the most frequent)
  z  = (2*nside-jr)*fact2;
  kshift = (jr - nside)&1;
  if( jr<nside ) { //then     ! north pole region
    nr = jr;
    z = 1. - nr*nr*fact1;
    kshift = 0;
  }
  else {
    if( jr>3*nside ) {// then ! south pole region
      nr = nl4 - jr;
      z = - 1. + nr*nr*fact1;
      kshift = 0;
    }
  }

  //c     computes the phi coordinate on the sphere, in [0,2Pi]
  //      jp = (jpll[face_num+1]*nr + jpt + 1 + kshift)/2;//  ! 'phi' number in the ring in {1,4*nr}
  jp = (jpll[face_num]*nr + jpt + 1 + kshift)/2;
  if( jp>nl4 ) jp = jp - nl4;
  if( jp<1 )   jp = jp + nl4;
  phi = (jp - (kshift+1)*0.5) * (piover2 / nr);

  int i, ksum, dir, nx, ny, pole;
  int xoffset, xoffset1;
  int straight[4];
  int pole_circle[2] = {5*ysize/6, ysize/6};
  int k[2], kprime[2];
  int vert1[4] = {0,1,1,0}, vert2[4] = {0,0,1,1};
  int vert_x[4], vert_y[4];
  int *tempx = NULL, *tempy = NULL;
  int xrange[2] = {+10000000, -100000000}, 
      yrange[2] = {+10000000, -100000000};
  float PixelCenter[2];
  float xfactor = xsize/(2*M_PI), yfactor = 0.5*ysize;
  float yf_inv = 1.0/yfactor;
  float dphi2 = 0.5*xsize/nl4, dz2 = 0.5*ysize*fact2;
  float cc = M_PI / (nside*sqrt(12.0));
  double sqr_fact;

  //  if (xpolygon != NULL)
  //    delete xpolygon;
  //  if (ypolygon != NULL)
  //    delete ypolygon;

  /* Compute the vertices of the pixel */

  PixelCenter[0] = xfactor*phi;
  PixelCenter[1] = yfactor*z + ysize/2;
  draw_poly = FALSE;

  // Equatorial regions
  if (jr > nside && jr < 3*nside) {
    npts = 4;
    xpolygon = new int[npts];
    ypolygon = new int[npts];
    for (i = 0; i < 4; i++) {
      dir = (i >> 1 & 1) ? +1 : -1;
      xpolygon[i] = (int) ( PixelCenter[0] - dir*((i+1)%2)*dphi2 + 0.5);
      ypolygon[i] = (int) ( PixelCenter[1] + dir*(i%2)*dz2 + 0.5);
    } // ENDFOR vertices
    draw_poly = TRUE;
  } // ENDIF equatorial

  // Polar regions
  else {

    k[0] = (jp-1)%nr;
    k[1] = k[0] + 1;
    kprime[0] = nr - (jp-1)%nr;
    kprime[1] = kprime[0] - 1;
    xoffset = ((jp-1)/nr) * (xsize/4);
    xoffset1 = (jp/nr) * (xsize/4);

    pole = (jr <= nside) ? +1 : -1;
    for (i = 0; i < 4; i++) {
      straight[i] = FALSE;
      ksum = k[vert1[i]] + kprime[vert2[i]];
      if (ksum > 0) {
        vert_x[i] = xoffset + (int) (xfactor * (piover2 * k[vert1[i]] / ksum) + 0.5);
        if (pole > 0)
          vert_y[i] = ysize - (int) (yfactor * (ksum*ksum*fact1) + 0.5);
        else
          vert_y[i] = (int) (yfactor * (ksum*ksum*fact1) + 0.5);

        // Accounting for equatorial vertices in transition pixels
        if (pole > 0 && vert_y[i] < pole_circle[0]) {
          vert_x[i] = (int) (PixelCenter[0] + 0.5);
          vert_y[i] = (int) (PixelCenter[1] - dz2 + 0.5);
          straight[i] = TRUE;
        }
        if (pole < 0 && vert_y[i] > pole_circle[1]) {
          vert_x[i] = (int) (PixelCenter[0] + 0.5);
          vert_y[i] = (int) (PixelCenter[1] + dz2 + 0.5);
          straight[i] = TRUE;
        }

        xrange[0] = min(xrange[0], vert_x[i]);
        yrange[0] = min(yrange[0], vert_y[i]);
        xrange[1] = max(xrange[1], vert_x[i]);
        yrange[1] = max(yrange[1], vert_y[i]);

      } // ENDIF ksum>0
      else
	straight[i] = TRUE;

    } // ENDFOR vertices

    // Adjust y-range if on the top/bottom row
    if (jr == 1 || jr == nl4-1)
      yrange[pole > 0] = (pole > 0) ? ysize : 0;

    nx = xrange[1]-xrange[0]+1;
    ny = yrange[1]-yrange[0]+1;
    tempx = new int[2*ny];
    tempy = new int[2*ny];
    //    zz = new double[ny];
    npts = 0;

//    for (i = 0; i < nx; i++)
//      zz[i] = (yrange[0]+i)*yf_inv;

    /* Compute the pixel boundaries (scan down in y and compute
       boundaries for x) */

    int xcenter = (int) (PixelCenter[0] + 0.5);
    int ycenter = (int) (PixelCenter[1] + 0.5);
    int transition;
    float dxdy = pole*3.0/8.0 * (float)xsize/(float)ysize;
    float ycoord;

    // Lower half
    for (i = yrange[0]; i < ycenter; i++) {
      tempy[npts] = tempy[npts+1] = i;
      ycoord = pole*(i-ysize/2)*yf_inv;

      transition = (pole > 0 && i <= pole_circle[0]) ||
	(pole < 0 && i > pole_circle[1]);

      // Left edge (SOUTHWEST)
      if (k[0] == 0 && pole < 0) // base pixel edge
	tempx[npts] = xoffset;
      else if (transition == TRUE) // transition to equator
	tempx[npts] = (int) (xcenter - dxdy*(i-yrange[0]) + 0.5);
      else { // 1/phi^2 curve
	if (pole > 0) {
	  sqr_fact = kprime[0]*cc;
	  tempx[npts] = xoffset + 
	    (int) (xfactor * (piover2 - sqr_fact / sqrt(1-ycoord)) + 0.5);
	} else {
	  sqr_fact = k[0]*cc;
	  tempx[npts] = xoffset + 
	    (int) (xfactor * (sqr_fact / sqrt(1-ycoord)) + 0.5);
	}
      }
      npts++;

      // Right Edge (SOUTHEAST)
      if (k[0] == nr && pole < 0) // base pixel edge
	tempx[npts] = xoffset1;
      else if (transition == TRUE) // transition to equator
	tempx[npts] = (int) (xcenter + dxdy*(i-yrange[0]) + 0.5);
      else { // polar 1/phi^2 curve
	if (pole > 0) {
	  sqr_fact = k[1]*cc;
	  tempx[npts] = xoffset +
	    (int) (xfactor * (sqr_fact / sqrt(1-ycoord)) + 0.5);
	} else {
	  sqr_fact = kprime[1]*cc;
	  tempx[npts] = xoffset + 
	    (int) (xfactor * (piover2 - sqr_fact / sqrt(1-ycoord)) + 0.5);
	}
      }
      npts++;

    } // ENDFOR i (y)

    // Upper half
    for (i = ycenter; i <= yrange[1]; i++) {
      tempy[npts] = tempy[npts+1] = i;
      ycoord = pole*(i-ysize/2)*yf_inv;

      transition = (pole > 0 && i < pole_circle[0]) ||
	(pole < 0 && i >= pole_circle[1]);

      // Left edge (NORTHWEST)
      if (k[0] == 0 && pole > 0) // base pixel edge
	tempx[npts] = xoffset;
      else if (transition == TRUE) // transition to equator
	tempx[npts] = (int) (xcenter + dxdy*(yrange[1]-i) + 0.5);
      else { // 1/phi^2 curve
	if (pole > 0) {
	  sqr_fact = k[0]*cc;
	  tempx[npts] = xoffset + 
	    (int) (xfactor * (sqr_fact / sqrt(1-ycoord)) + 0.5);
	} else {
	  sqr_fact = kprime[0]*cc;
	  tempx[npts] = xoffset + 
	    (int) (xfactor * (piover2 - sqr_fact / sqrt(1-ycoord)) + 0.5);
	}
      }
      npts++;

      // Right Edge (NORTHEAST)
      if (k[0] == nr && pole > 0) // base pixel edge
	tempx[npts] = xoffset1;
      else if (transition == TRUE) // transition to equator
	tempx[npts] = (int) (xcenter - dxdy*(yrange[1]-i) + 0.5);
      else { // polar 1/phi^2 curve
	if (pole > 0) {

	  sqr_fact = kprime[1]*cc;
	  tempx[npts] = xoffset + 
	    (int) (xfactor * (piover2 - sqr_fact / sqrt(1-ycoord)) + 0.5);
	} else {
	  sqr_fact = k[1]*cc;
	  tempx[npts] = xoffset +
	    (int) (xfactor * (sqr_fact / sqrt(1-ycoord)) + 0.5);
	}
      }
      npts++;

    } // ENDFOR i (y)

    xpolygon = new int[npts];
    ypolygon = new int[npts];
    for (i = 0; i < npts; i++) {
      xpolygon[i] = tempx[i];
      ypolygon[i] = ysize-tempy[i];
    }

    delete tempx;
    delete tempy;
    //    delete zz;

  } // ENDELSE poles

  return SUCCESS;

}

