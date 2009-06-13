#ifdef DOUB
#include <dfftw.h>
#include <drfftw.h>
#else
#include <sfftw.h>
#include <srfftw.h>
#endif

#if defined ADD0US
void convert_to_fftw_complex(fftw_real *in, fftw_complex *out) {
#elif defined ADD2US
void convert_to_fftw_complex__(fftw_real *in, fftw_complex *out) {
#else
void convert_to_fftw_complex_(fftw_real *in, fftw_complex *out) {
#endif
  out = (fftw_complex *) in;
}

#if defined ADD0US
void convert_to_fftw_real(fftw_complex *in, fftw_real *out) {
#elif defined ADD2US
void convert_to_fftw_real__(fftw_complex *in, fftw_real *out) {
#else
void convert_to_fftw_real_(fftw_complex *in, fftw_real *out) {
#endif
  out = (fftw_real *) in;
}

#if defined ADD0US
void my_rfftwnd_mpi(rfftwnd_plan *p, int *n_fields, fftw_complex *cdata,
		    fftw_real *rdata, fftw_real *work, int *use_work,
		    int *ioutput_order) {
#elif defined ADD2US
void my_rfftwnd_mpi__(rfftwnd_plan *p, int *n_fields, fftw_complex *cdata,
		      fftw_real *rdata, fftw_real *work, int *use_work,
		      int *ioutput_order) {
#else
void my_rfftwnd_mpi_(rfftwnd_plan *p, int *n_fields, fftw_complex *cdata,
		     fftw_real *rdata, fftw_real *work, int *use_work,
		     int *ioutput_order) {
#endif
  return;
}
