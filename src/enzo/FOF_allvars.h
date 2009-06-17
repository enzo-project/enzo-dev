#ifndef __FOF_ALLVARS_H
#define __FOF_ALLVARS_H

#include "FOF_nrutil.h"

/************************************************************************
  CONSTANTS
************************************************************************/

#define  KERNEL_TABLE 10000
#define  PI               3.1415927
#define  GRAVITY     6.672e-8
#define  SOLAR_MASS  1.989e33
#define  SOLAR_LUM   3.826e33
#define  RAD_CONST   7.565e-15
#define  AVOGADRO    6.0222e23
#define  BOLTZMANN   1.3806e-16
#define  GAS_CONST   8.31425e7
#define  C           2.9979e10
#define  PLANCK      6.6262e-27
#define  CM_PER_MPC  3.085678e24
#define  PROTONMASS  1.6726e-24
#define  ELECTRONMASS  9.10953e-28
#define  THOMPSON    6.65245e-25

#define  HUBBLE      3.2407789e-18   /* in h/sec */

#define  SEC_PER_MEGAYEAR   3.155e13
#define  SEC_PER_YEAR       3.155e7

#define  GAMMA         (5.0/3)
#define  GAMMA_MINUS1  (GAMMA-1)

/************************************************************************
   STRUCTURES
************************************************************************/

struct gr_data
{
  int Len, Tag;
};

struct FOF_particle_data 
{
  double  	Pos[3];
  float  	Vel[3];
  int    	Type;
  int    	ID;
  int    	MinID;
  int    	GrLen;
  float  	Mass;
  float  	Mfs, Mclouds, Sfr;
  float  	Energy;
  float  	Rho;
  int    	PartID;
  int           slab;
};

struct id_data 
{
  int    ID;
  int    index;
};
 

struct idmin_data 
{
  int    minID;
  int    index;
  int    len;
};

struct grouptree_data
{
  int lefthead, leftlen;
  int righthead, rightlen;
  int left, right;
};
 
/************************************************************************
  HALO FINDER VARIABLES
  -- Previously everything was global in the standalone version.  Put 
     most variables inside this structure.
*************************************************************************/

struct FOFData {

  double  LinkLength;
  int     GroupMinLen;
  int     MaxPlacement;
  int     Grid;

  double  Time;
  double  BoxSize;
  double  leftEdge[3], rightEdge[3];

  double  SearchRadius;

  int     NumPart;   /* total particle number */
  int     *NpartInGrids;

  int     *Nslab, *Nshadow;
  int     Nlocal, *Noffset;
  int     Nlocal_in_file;
  int     *NtoLeft, *NtoRight;
  int     *Nslab_local, *NtoLeft_local, *NtoRight_local;

  int     Ncontrib, *ContribID, *ContribHead;

  int     NLinkAcross;

  double  GridExtension, GridCorner[3];
  int     ***GridFirst, ***GridLast, ***GridFlag;
  int     *GridNext;
  int     *Head,*Next;   /* for link-lists of groups */
  int     *Tail,*Len;
  int     Ngroups, NgroupsAll, *NgroupsList;

  int Nx, Ny, Nz;

  gr_data *GroupDat, *GroupDatAll;

  FOF_particle_data *P, *Pbuf_local;

  /************************* UNITS *************************/

  double UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s, 
    UnitTime_in_s, UnitTime_in_Megayears, UnitDensity_in_cgs, 
    UnitPressure_in_cgs, UnitCoolingRate_in_cgs, UnitEnergy_in_cgs, 
    G, Hubble, EnzoMassUnit, EnzoVelocityUnit, EnzoTimeUnit;

  /************************* SUBFIND VARIABLES *************************/

  int    DesDensityNgb;
  int    DesLinkNgb;
  float  Theta;

  double Omega;
  double OmegaLambda;
  double H0;
  float  Epsilon;

  grouptree_data *GroupTree;

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

};

#include "FOF_ngbtree.h"

#endif
