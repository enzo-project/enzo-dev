#include "allvars.h"


 int     ThisTask, NTask;
 int     CycleNumber;

 double  Time;
 double  BoxSize;
 double  RhoCritical0;
 double  leftEdge[3], rightEdge[3];

 double  SearchRadius;

 int     ParticleTypeInFile;
 PINT    NumPart;   /* total particle number */
 int     *NpartInGrids;

 int     *Nslab, *Nshadow;
 int     Nlocal, *Noffset;
 int     Nlocal_in_file;
 int     *NtoLeft, *NtoRight;
 int     *Nslab_local, *NtoLeft_local, *NtoRight_local;

 int     Ncontrib, *ContribID, *ContribHead;

 int     NLinkAccross;

 double  GridExtension, GridCorner[3];
 int     ***GridFirst, ***GridLast, ***GridFlag;
 int     *GridNext;
 int     *Head,*Next;   /* for link-lists of groups */
 int     *Tail,*Len;
 int     Ngroups, NgroupsAll, *NgroupsList;


 struct  gr_data
 *GroupDat, *GroupDatAll;

 int    Nx,Ny,Nz;

/* Quantities for all particles */

 struct particle_data 
 *P, *Pbuf_local;

 
double  UnitLength_in_cm,
               UnitMass_in_g,
               UnitVelocity_in_cm_per_s,
               UnitTime_in_s,
               UnitTime_in_Megayears,
               UnitDensity_in_cgs,
               UnitPressure_in_cgs,
               UnitCoolingRate_in_cgs,
               UnitEnergy_in_cgs,
               G,
               Hubble,
               EnzoMassUnit,
               EnzoVelocityUnit,
               EnzoTimeUnit;

/* ----------------------------------------------------------
 *
 *  Variables for subfind
 *
 */

 int    DesDensityNgb;
 int    DesLinkNgb;
 float  Theta;
 
 double Time;
 double Omega;
 double OmegaLambda;
 double G;
 double H0;
 float  Epsilon;


 struct grouptree_data
 *GroupTree;


 int   NumInGroup;
 int   NSubGroups, NSubGroupsAll;

 int    AnzNodes;
 int    MaxNodes;


 int   *Id;
 int   *SubGroupLen, *SubGroupTag;
 int   *GroupTag, *NextFinal;
 int   *Index;
 int   *NewHead;
 int   *Node;
 int   *SortId;
 float *Density;
 float *Potential;
 float *Energy;


/* tabulated smoothing kernel */

 double  Kernel[KERNEL_TABLE+1],
               KernelDer[KERNEL_TABLE+1],
               KernelRad[KERNEL_TABLE+1];
