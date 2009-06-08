#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "macros_and_parameters.h"

void quickSort(int a[], int b[], int l, int r);
int partition(int a[], int b[], int l, int r);

// initializatiion function realted to HEALPix
int mkPix2xy(long *ipix2x, long *ipix2y) 
{
  long  kpix,jpix,ix,iy,ip,id;

  for(kpix=0; kpix <= 1023; kpix++) {          
    jpix = kpix;
    ix = 0;
    iy = 0;
    ip = 1 ;            
    while (jpix != 0) { 
      id = (int)fmod(jpix,2);   
      jpix = jpix/2 ;
      ix = id*ip+ix ;
      id = (int)fmod(jpix,2);
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
    J  = I-1;//            !pixel numbers
    K  = 0;//
    IP = 1;//
    truc : if( J==0 ) {
      x2pix[I-1] = K;
      y2pix[I-1] = 2*K;
    }
    else {
      ID = (int)fmod(J,2);
      J  = J/2;
      K  = IP*ID+K;
      IP = IP*4;
      goto truc;
    }
  }     

  return SUCCESS;
  
}

#include "typedefs.h"
#include "global_data.h"

/***********************************************************************/

// HEALPix fuction that updates directional vector given a pixel number
#define nsMax   8192L  
#define halfpi  1.5707963
int pix2vec_nest(long nside, long ipix, FLOAT *v)
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
      
  int jrll[12], jpll[12];// ! coordinate of the lowest corner of each face
  //      data jrll/2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4/ ! in unit of nside
  //      data jpll/1, 3, 5, 7, 0, 2, 4, 6, 1, 3, 5, 7/ ! in unit of nside/2
  jrll[0]=2;
  jrll[1]=2;
  jrll[2]=2;
  jrll[3]=2;
  jrll[4]=3;
  jrll[5]=3;
  jrll[6]=3;
  jrll[7]=3;
  jrll[8]=4;
  jrll[9]=4;
  jrll[10]=4;
  jrll[11]=4;
  jpll[0]=1;
  jpll[1]=3;
  jpll[2]=5;
  jpll[3]=7;
  jpll[4]=0;
  jpll[5]=2;
  jpll[6]=4;
  jpll[7]=6;
  jpll[8]=1;
  jpll[9]=3;
  jpll[10]=5;
  jpll[11]=7;
      
      
  if( nside<1 || nside>ns_max ) {
    fprintf(stderr, "%s (%"ISYM"): nside out of range: %ld\n", __FILE__, __LINE__, nside);
    exit(0);
  }
  npix = 12 * nside*nside;
  if( ipix<0 || ipix>npix-1 ) {
    fprintf(stderr, "%s (%"ISYM"): ipix out of range: %ld\n", __FILE__, __LINE__, ipix);
    exit(0);
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
  ipf = (int)fmod(ipix,npface);//  ! pixel number in the face {0,npface-1}

  //c     finds the x,y on the face (starting from the lowest corner)
  //c     from the pixel number
  ip_low = (int)fmod(ipf,1024);//       ! content of the last 10 bits
  ip_trunc =   ipf/1024 ;//       ! truncation of the last 10 bits
  ip_med = (int)fmod(ip_trunc,1024);//  ! content of the next 10 bits
  ip_hi  =     ip_trunc/1024   ;//! content of the high weight 10 bits

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
  kshift = (int)fmod(jr - nside, 2);
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
  double piover2 = 0.5*M_PI, pi = M_PI, twopi = 2.0*M_PI;
  int    ns_max = 8192;
//  static int x2pix[128], y2pix[128];
//  static char setup_done = 0;
  
  if( nside<1 || nside>ns_max ) {
    fprintf(stderr, "%s (%"ISYM"): nside out of range: %ld\n", __FILE__, __LINE__, nside);
    exit(0);
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
    
    if( ifp==ifm ) face_num = (int)fmod(ifp,4) + 4; /* faces 4 to 7 */
    else if( ifp<ifm ) face_num = (int)fmod(ifp,4); /* (half-)faces 0 to 3 */
    else face_num = (int)fmod(ifm,4) + 8;           /* (half-)faces 8 to 11 */
    
    ix = (int)fmod(jm, ns_max);
    iy = ns_max - (int)fmod(jp, ns_max) - 1;
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
  
  ix_low = (int)fmod(ix,128);
  ix_hi  =     ix/128;
  iy_low = (int)fmod(iy,128);
  iy_hi  =     iy/128;

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
    fprintf(stderr, "%s (%"ISYM"): nside out of range: %ld\n", __FILE__, __LINE__, nside);
    return FAIL;
  }
  npix = 12 * nside*nside;
  if( ipix<0 || ipix>npix-1 ) {
    fprintf(stderr, "%s (%"ISYM"): ipix out of range: %ld\n", __FILE__, __LINE__, ipix);
    return FAIL;
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
  ipf = (int)fmod(ipix,npface);//  ! pixel number in the face {0,npface-1}

  //c     finds the x,y on the face (starting from the lowest corner)
  //c     from the pixel number
  ip_low = (int)fmod(ipf,1024);//       ! content of the last 10 bits
  ip_trunc =   ipf/1024 ;//       ! truncation of the last 10 bits
  ip_med = (int)fmod(ip_trunc,1024);//  ! content of the next 10 bits
  ip_hi  =     ip_trunc/1024   ;//! content of the high weight 10 bits

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
  kshift = (int)fmod(jr - nside, 2);
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

  int i, ip, v, ksum, dir, nx, ny, thisx, pole, vp1, for_check;
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
  float xf_inv = 1.0/xfactor, yf_inv = 1.0/yfactor;
  float dphi2 = 0.5*xsize/nl4, dz2 = 0.5*ysize*fact2;
  float dphip = 0.125*xsize/nr;
  float cc = M_PI / (nside*sqrt(12.0));
  double *zz, sqr_fact;
  float (*RFN) (float);

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
      xpolygon[i] = (int) roundf( PixelCenter[0] - dir*((i+1)%2)*dphi2 );
      ypolygon[i] = (int) roundf( PixelCenter[1] + dir*(i%2)*dz2 );
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
	vert_x[i] = xoffset + (int) roundf(xfactor * (piover2 * k[vert1[i]] / ksum));
	if (pole > 0)
	  vert_y[i] = ysize - (int) roundf(yfactor * (ksum*ksum*fact1));
	else
	  vert_y[i] = (int) roundf(yfactor * (ksum*ksum*fact1));

	// Accounting for equatorial vertices in transition pixels
	if (pole > 0 && vert_y[i] < pole_circle[0]) {
	  vert_x[i] = (int) roundf(PixelCenter[0]);
	  vert_y[i] = (int) roundf(PixelCenter[1] - dz2);
	  straight[i] = TRUE;
	}
	if (pole < 0 && vert_y[i] > pole_circle[1]) {
	  vert_x[i] = (int) roundf(PixelCenter[0]);
	  vert_y[i] = (int) roundf(PixelCenter[1] + dz2);
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

    int xcenter = (int) roundf(PixelCenter[0]);
    int ycenter = (int) roundf(PixelCenter[1]);
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
	tempx[npts] = (int) roundf(xcenter - dxdy*(i-yrange[0]));
      else { // 1/phi^2 curve
	if (pole > 0) {
	  sqr_fact = kprime[0]*cc;
	  tempx[npts] = xoffset + 
	    (int) roundf(xfactor * (piover2 - sqr_fact / sqrt(1-ycoord)));
	} else {
	  sqr_fact = k[0]*cc;
	  tempx[npts] = xoffset + 
	    (int) roundf(xfactor * (sqr_fact / sqrt(1-ycoord)));
	}
      }
      npts++;

      // Right Edge (SOUTHEAST)
      if (k[0] == nr && pole < 0) // base pixel edge
	tempx[npts] = xoffset1;
      else if (transition == TRUE) // transition to equator
	tempx[npts] = (int) roundf(xcenter + dxdy*(i-yrange[0]));
      else { // polar 1/phi^2 curve
	if (pole > 0) {
	  sqr_fact = k[1]*cc;
	  tempx[npts] = xoffset +
	    (int) roundf(xfactor * (sqr_fact / sqrt(1-ycoord)));
	} else {
	  sqr_fact = kprime[1]*cc;
	  tempx[npts] = xoffset + 
	    (int) roundf(xfactor * (piover2 - sqr_fact / sqrt(1-ycoord)));
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
	tempx[npts] = (int) roundf(xcenter + dxdy*(yrange[1]-i));
      else { // 1/phi^2 curve
	if (pole > 0) {
	  sqr_fact = k[0]*cc;
	  tempx[npts] = xoffset + 
	    (int) roundf(xfactor * (sqr_fact / sqrt(1-ycoord)));
	} else {
	  sqr_fact = kprime[0]*cc;
	  tempx[npts] = xoffset + 
	    (int) roundf(xfactor * (piover2 - sqr_fact / sqrt(1-ycoord)));
	}
      }
      npts++;

      // Right Edge (NORTHEAST)
      if (k[0] == nr && pole > 0) // base pixel edge
	tempx[npts] = xoffset1;
      else if (transition == TRUE) // transition to equator
	tempx[npts] = (int) roundf(xcenter - dxdy*(yrange[1]-i));
      else { // polar 1/phi^2 curve
	if (pole > 0) {
	  sqr_fact = kprime[1]*cc;
	  tempx[npts] = xoffset + 
	    (int) roundf(xfactor * (piover2 - sqr_fact / sqrt(1-ycoord)));
	} else {
	  sqr_fact = k[1]*cc;
	  tempx[npts] = xoffset +
	    (int) roundf(xfactor * (sqr_fact / sqrt(1-ycoord)));
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

