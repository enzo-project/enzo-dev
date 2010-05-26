#ifndef __radtrans_healpix_h_
#define __radtrans_healpix_h_
int mkPix2xy(long *ipix2x, long *ipix2y);
int mk_xy2pix(int *x2pix, int *y2pix);
int pix2vec_nest(long nside, long ipix, FLOAT *v, double (*vertex)[3]=0);
int vec2pix_nest( const long nside, FLOAT *vec, long *ipix);
int pix2coord_nest( long nside, long ipix, int xsize, int ysize, int &npts, int * &xpolygon, int * &ypolygon, 
		    int &draw_poly);
#endif
