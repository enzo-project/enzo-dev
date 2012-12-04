#ifndef __FOF_ALLVARS_H
#define __FOF_ALLVARS_H

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

#define  PLANCK      6.6262e-27
#define  CM_PER_MPC  3.085678e24
#define  PROTONMASS  1.6726e-24
#define  ELECTRONMASS  9.10953e-28
#define  THOMPSON    6.65245e-25

#define  HUBBLE      3.2407789e-18   /* in h/sec */

#define  SEC_PER_MEGAYEAR   3.155e13
#define  SEC_PER_YEAR       3.155e7

#ifndef __GAMMA_DEFINED_
#define __GAMMA_DEFINED_
#define  GAMMA         (5.0/3)
#define  GAMMA_MINUS1  (GAMMA-1)
#endif

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
  char    	Type;
  PINT    	ID;
  PINT    	MinID;
  Eint32    	GrLen;
  float  	Mass;
  float  	Attr[MAX_NUMBER_OF_PARTICLE_ATTRIBUTES];
  float  	Energy;
  float  	Rho;
  PINT    	PartID;
  Eint32        slab;
  //int           level;
  //int           GridID;
};

struct id_data 
{
  PINT    ID;
  PINT    index;
};
 

struct idmin_data 
{
  PINT    minID;
  PINT    index;
  int     len;
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

  double  RhoCritical0;
  double  Time;
  double  BoxSize;

  double  SearchRadius;

  PINT    NumPart;   /* total particle number */

  PINT    *Nslab, *Nshadow;
  PINT    Nlocal, *Noffset;
  PINT    Nlocal_in_file;
  PINT    *NtoLeft, *NtoRight;

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

  FOF_particle_data *P;

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

#include "FOF_nrutil.h"
#include "FOF_ngbtree.h"

#endif
