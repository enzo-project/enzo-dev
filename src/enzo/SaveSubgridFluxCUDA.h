#ifndef _SAVESUBGRIDFLUX_H_
#define _SAVESUBGRIDFLUX_H_
 
void SaveSubgridFluxCUDA(float *LeftFlux, float *RightFlux, float *Flux,
                         float *Tmp1, float *Tmp2,
                         const int dimx, const int dimy, const int dimz,
                         const int fistart, const int fiend,
                         const int fjstart, const int fjend,
                         const int lface, const int rface, const int dir);

#endif
