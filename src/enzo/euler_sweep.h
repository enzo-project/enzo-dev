extern "C" void FORTRAN_NAME(pgas2d_dual)(
	        float *dslice, float *eslice, float *geslice, float *pslice,
		float *uslice, float *vslice, float *wslice, float *eta1, 
		float *eta2, int *idim, int *jdim, int *i1, int *i2, 
		int *j1, int *j2, float *gamma, float *pmin);

extern "C" void FORTRAN_NAME(pgas2d)(
	        float *dslice, float *eslice, float *pslice,
		float *uslice, float *vslice, float *wslice, int *idim, 
		int *jdim, int *i1, int *i2, int *j1, int *j2, float *gamma, 
		float *pmin);

extern "C" void FORTRAN_NAME(calcdiss)(
	        float *dslice, float *eslice, float *uslice, float *v, 
		float *w, float *pslice, float *dx, float *dy, float *dz, 
		int *idim, int *jdim, int *kdim, int *i1, int *i2, int *j1, 
		int *j2, int *k, int *nzz, int *idir, int *dimx, int *dimy, 
		int *dimz, float *dt, float *gamma, int *idiff, int *iflatten, 
		float *diffcoef, float *flatten);

extern "C" void FORTRAN_NAME(inteuler)(
		float *dslice, float *pslice, int *gravity, float *grslice, 
		float *geslice, float *uslice, float *vslice, float *wslice, 
		float *dxi, float *flatten, int *idim, int *jdim, int *i1, 
		int *i2, int *j1, int *j2, int *idual, float *eta1, 
		float *eta2, int *isteep, int *iflatten,
		int *iconsrec, int *iposrec,
		float *dt, float *gamma, int *ipresfree, float *dls, float *drs, 
		float *pls, float *prs, float *gels, float *gers, float *uls, 
		float *urs, float *vls, float *vrs, float *wls, float *wrs, 
		int *ncolor, float *colslice, float *colls, float *colrs);

extern "C" void FORTRAN_NAME(twoshock)(
		float *dls, float *drs, float *pls, float *prs, 
		float *uls, float *urs, int *idim, int *jdim, int *i1,
		int *i2, int *j1, int *j2, float *dt, float *gamma, 
		float *pmin, int *ipresfree, float *pbar, float *ubar,
		int *gravity, float *grslice, int *idual, float *eta1);

extern "C" void FORTRAN_NAME(flux_twoshock)(
           float *dslice, float *eslice, float *geslice, float *uslice, 
	   float *vslice, float *wslice, float *dx,
           float *diffcoef, int *idim, int *jdim, int *i1, int *i2, 
	   int *j1, int *j2, float *dt, float *gamma, int *idiff,
           int *idual, float *eta1, int *ifallback,
           float *dls, float *drs, float *pls, float *prs, 
	   float *gels, float *gers, 
           float *uls, float *urs, float *vls, float *vrs, 
	   float *wls, float *wrs, float *pbar, float *ubar,
           float *df, float *ef, float *uf, float *vf, float *wf, float *gef, 
	   float *ges,
           int *ncolor, float *colslice, float *colls, float *colrs, float *colf);

extern "C" void FORTRAN_NAME(flux_hll)(
           float *dslice, float *eslice, float *geslice, float *uslice, 
	   float *vslice, float *wslice, float *dx,
           float *diffcoef, int *idim, int *jdim, int *i1, int *i2, 
	   int *j1, int *j2, float *dt, float *gamma,
	   int *idiff, int *idual, float *eta1, int *ifallback,
           float *dls, float *drs, float *pls, float *prs, 
           float *uls, float *urs, float *vls, float *vrs, 
	   float *wls, float *wrs, float *gels, float *gers, 
           float *df, float *ef, float *uf, float *vf, float *wf, float *gef, 
	   float *ges,
           int *ncolor, float *colslice, float *colls, float *colrs, float *colf);

extern "C" void FORTRAN_NAME(flux_hllc)(
           float *dslice, float *eslice, float *geslice, float *uslice, 
	   float *vslice, float *wslice, float *dx,
           float *diffcoef, int *idim, int *jdim, int *i1, int *i2, 
	   int *j1, int *j2, float *dt, float *gamma,
	   int *idiff, int *idual, float *eta1, int *ifallback,
           float *dls, float *drs, float *pls, float *prs, 
           float *uls, float *urs, float *vls, float *vrs, 
	   float *wls, float *wrs, float *gels, float *gers, 
           float *df, float *ef, float *uf, float *vf, float *wf, float *gef, 
	   float *ges,
           int *ncolor, float *colslice, float *colls, float *colrs, float *colf,
           float *dfloor);

extern "C" void FORTRAN_NAME(euler)(
		float *dslice, float *pslice, float *grslice, 
		float *geslice, float *uslice, float *vslice, float *wslice, 
		float *dx, float *diffcoef, int *idim, int *jdim, int *i1, 
		int *i2, int *j1, int *j2, float *dt, float *gamma, int *idiff, 
		int *gravity, int *idual, float *eta1, float *eta2, float *df, 
		float *ef, float *uf, float *vf, float *wf, float *gef,
		float *ges,
		int *ncolor, float *colslice, float *colf, float *dfl, float *eceil);
