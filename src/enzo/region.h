/* Region structure used for Parallel FFT */

#ifndef REGION_DEFINED__
#define REGION_DEFINED__

struct region {
  int StartIndex[MAX_DIMENSION];
  int RegionDim[MAX_DIMENSION];
  int Processor;
  float *Data;
};

/*
struct commSndRcv {
 struct region *Sends, *Receives;
 int sends, receives, SendSize, ReceiveSize;
} *cSndRcv;
*/

struct commSndRcv {
 struct region *Sends, *Receives;
 int sends, receives, SendSize, ReceiveSize;
};

/* int first_pass = 0; */

#endif
