#ifndef _CUMHD_H_
#define _CUMHD_H_

#define MAX_SPECIES 20

typedef struct {
  // baryon fields (primitives)
  float *D, *V1, *V2, *V3, *TE, *B1, *B2, *B3, *Phi;
  // old baryon fields
  float *OldD, *OldV1, *OldV2, *OldV3, *OldTE, *OldB1, *OldB2, *OldB3, *OldPhi;
  // fluxes
  float *FluxD, *FluxS1, *FluxS2, *FluxS3, *FluxTau, *FluxB1, *FluxB2, *FluxB3, *FluxPhi;
  // dU
  float *dUD, *dUS1, *dUS2, *dUS3, *dUTau, *dUB1, *dUB2, *dUB3, *dUPhi; 
  // source terms
  float *divB, *gradPhi, *AccelerationField[3], *DrivingForce[3];
  // Color
  float *Species[MAX_SPECIES];
  float *OldSpecies[MAX_SPECIES];
  float *FluxSpecies[MAX_SPECIES];
  float *dUSpecies[MAX_SPECIES];
  float **SpeciesArray;
  float **OldSpeciesArray;
  float **FluxSpeciesArray;
  float **dUSpeciesArray;
} cuMHDData;

#endif
