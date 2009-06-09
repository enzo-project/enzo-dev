
#include "nrsrc/nrutil.h"
#include "ngbtree.h"



#define  PI               3.1415927
#define KERNEL_TABLE 10000



extern double Omega;
extern double OmegaLambda;
extern double G;
extern double H0;
extern char path[200];
extern int  Files;




extern float Softening;
extern float SofteningMaxPhys;    /* mas physical softening lengh */
extern float Theta;
extern float Epsilon;


extern struct grouptree_data
{
  int lefthead, leftlen;
  int righthead, rightlen;
  int left, right;
} *GroupTree;


extern int    AnzNodes;
extern int    MaxNodes;





extern int Snapshot;
extern int GroupNr;

extern int DesDensityNgb;
extern int DesLinkNgb;


extern double  Time;
extern int     NumPart;
extern double  PartMass;


extern struct particle_data 
{
  float  Pos[3];
  float  Vel[3];
} *P;

extern int *Id;


extern int *GrIds;


extern int Ngroups;
extern int *GroupLen, *GroupTag, *NextFinal;
extern int *NSubGroups, *FirstEntrySubStructureList;

extern int Nsubstructures; /* cumulative number of substructures */
extern int Nsubids;  /* cumulative number of particles in substructures */


extern int *SubGroupLen, *SubGroupTag;

extern int NumInGroup;

extern int Nsubgroups_B, NsubIds_B;
extern int *SubLen_B, *SubFirst_B,*SubParent_B;
extern int *SubIds_B;
extern int *Tag_B;


extern int *Index;
extern int *NewHead;
extern int *Node;
extern int *SortId;
extern float *Density;
extern float *Potential;
extern float *Energy;



/* tabulated smoothing kernel */

extern double  Kernel[KERNEL_TABLE+1],
               KernelDer[KERNEL_TABLE+1],
               KernelRad[KERNEL_TABLE+1];











