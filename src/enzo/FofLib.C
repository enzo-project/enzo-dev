/* foflib.c
   Mark Krumholz, 3/9/00
   Modified by Mark Krumholz, 8/22/00
   Modified by Nathan Goldbaum, December 2011 (include in ENZO)

   This is a set of library routines to perform friends-of-friends
   groupings and similar operations on sets of particles. Detailed
   explanations of each routine are given below. All algorithms in the
   library run in N log N time, where N is the number of particles
   being analyzed. The basic functions in the library are: 

   - Form friends-of-friends particle groups for all particles in a
      given list, using a fixed linking length to determine which
      particles are friends. Return an array listing a group
      assignment for each particle. Routine: Fof().

   - Form friends-of-friends particle groups for all particles in a
      given list, with each particle having its own linking length to
      determine which particles are friends. Return an array listing a
      group assignment for each particle. Routine: FofVar().

   - Form friends-of-friends particle groups for all particles in a
      given list, using a fixed linking length to determine which
      particles are friends. Return an array listing a group
      assignment for each particle and an array of arrays listing the
      particles that are members of each group. Routine: FofList().

   - Form friends-of-friends particle groups for all particles in a
      given list, with each particle having its own linking length to
      determine which particles are friends. Return an array listing a
      group assignment for each particle and an array of arrays
      listing the particles that are members of each group. Routine:
      FofVarList().

   - Take the group assignments output by Fof() or FofVar() and remove
      all groups smaller than a certain specified size. Particles
      previously assigned to these groups will now be listed as
      belonging to group -1, and groups will be renumbered 0 through
      number of groups - 1. Routine: FofPrune().

   - Take the group assignments and group membership lists output by
      FofList() or FofVarList() and remove all groups smaller than a
      certain specified size. Particles previously assigned to these
      groups will now be listed as belonging to group -1, and groups
      will be renumbered 0 through number of groups -1. Group
      membership lists will be purged of groups smaller than the given
      size. Routine: FofListPrune().

   - Find the N nearest neighbors to every particle in a given
      list. Return an array listing the neighbors for each
      particle. The neighbors are ordered by distances from the
      particle to which they are a neighbor, with the closest one
      first in the list. Routine: NearNeighbor().

   - Find the N nearest neighbors to some of the particles in a given
      list. Return an array listing the neighbors for each particle
      whose neighbors were found. The neighbors are ordered by
      distances from the particle to which they are a neighbor, with
      the closest one first in the list. Routine:
      NearNeighborPartial().

   - For each particle in a given list, find all other particles
      within a specified distance. Return a list of these "neighbors"
      for each particle. The neighbors are in no particular
      order. Routine: FindNeighbor().

   - For a subset of the particles in a given list, find all other
      particles within a specified distance. Return a list of these
      "neighbors" for each particle in the subset whose neighbors are
      found. The neighbors are in no particular order. Routine:
      FindNeighborPartial().


   The syntax for each is:

   int FofVar ( int npart, Real *x, Real *link, int *group, int
     **groupsize )
     IN npart: number of particles
     IN *x: 3*npart element array of particle positions. The x, y, and
       z coordinates of particle n are elements 3*n, 3*n+1, and 3*n+2.
     IN *link: npart element array of linking lengths for particles
     OUT *group: npart element array of group number assignments for
       each particle. The pointer passed must point to valid,
       allocated memory. Groups are numbered sequentially 0 through
       ngroups-1, in no particular order.
     OUT **groupsize: a pointer to an array with number of elements
       equal to number of groups found. Each element gives the number
       of particles in that group number. When passed this pointer
       should not point to valid memory.
     RETURNS: number of groups found.

   int FofVarList ( int npart, Real *x, Real *link, int *group,
     int **groupsize, int ***grouplist )
     IN npart: number of particles
     IN *x: 3*npart element array of particle positions. The x, y, and
       z coordinates of particle n are 3*n, 3*n+1, and 3*n+2.
     IN *link: npart element array of linking lengths for particles
     OUT *group: npart element array of group number assignments for
       each particle. The pointer passed must point to valid,
       allocated memory. Groups are numbered sequentially 0 through
       ngroups-1, in no particular order.
     OUT **groupsize: a pointer to an array with number of elements
       equal to number of groups found. Each element gives the number
       of particles in that group number. When passed this routine
       should not point to valid memory.
     OUT ***grouplist: a pointer to an array with number of elements
       equal to number of groups found. Each element of the array is
       itself an array listing the indices of all particles in that
       group.
     RETURNS: number of groups found.

   int Fof ( int npart, Real *x, Real link, int *group, int
     **groupsize )
     IN npart: number of particles
     IN *x: 3*npart element array of particle positions. The x, y, and
       z coordinates of particle n are 3*n, 3*n+1, and 3*n+2.
     IN link: linking length to search
     OUT *group: npart element array of group number assignments for
       each particle. The pointer passed must point to valid,
       allocated memory. Groups are numbered sequentially 0 through
       ngroups-1, in no particular order.
     OUT **groupsize: a pointer to an array with number of elements
       equal to number of groups found. Each element gives the number
       of particles in that group number. When passed this routine
       should not point to valid memory.
     RETURNS: number of groups found.

   int FofList ( int npart, Real *x, Real link, int *group, int
     **groupsize, int ***grouplist )
     IN npart: number of particles
     IN *x: 3*npart element array of particle positions. The x, y, and
       z coordinates of particle n are 3*n, 3*n+1, and 3*n+2.
     IN link: linking length to search
     OUT *group: npart element array of group number assignments for
       each particle. The pointer passed must point to valid,
       allocated memory. Groups are numbered sequentially 0 through
       ngroups-1, in no particular order.
     OUT **groupsize: a pointer to an array with number of elements
       equal to number of groups found. Each element gives the number
       of particles in that group number. When passed this routine
       should not point to valid memory.
     OUT ***grouplist: a pointer to an array with number of elements
       equal to number of groups found. Each element of the array is
       itself an array listing the indices of all particles in that
       group.
     RETURNS: number of groups found.

   int FofPrune ( int npart, int ngroup, int *group, int **groupsize,
     int minsize)
     IN npart: number of particles.
     IN ngroup: number of groups passed in.
     INOUT *group: group affiliations of particles, as returned by
       Fof() or FofVar(). When returned, groups are renumbered from 0
       to nremaining-1.
     INOUT **groupsize: pointer to list of sizes of each group.
     IN minsize: minimum size group to keep.
     RETURNS: number of groups left.

   int FofListPrune ( int npart, int ngroup, int *group, int
     **groupsize, int ***grouplist, int minsize)
     IN npart: number of particles.
     IN ngroup: number of groups to be pruned.
     INOUT *group: group affiliations of particles, as returned by
       FofList() or FofVarList(). When returned, groups are renumbered
       from 0 to nremaining-1.
     INOUT **groupsize: pointer to list of sizes of each group.
     INOUT ***grouplist: pointer to list of particles in each group.
     IN minsize: minimum size group to keep.
     RETURNS: number of groups left.

   void NearNeighbor ( int npart, Real *x, int nneighbor, int
     *neighborlist )
     IN npart: number of particles.
     IN *x: 3*npart element array of particle positions. The x, y, and
       z coordinates of particle n are 3*n, 3*n+1, and 3*n+2.
     IN nneighbor: number of neighbors of each particle to find
     OUT *neighborlist: a pointer to a npart x nneighbor block of
       integers. The index of the jth nearnest neighbor of particle i
       is returned in element i*nneighbor+j of the array.

   void NearNeighborPartial ( int npart, Real *x, int nneighbor, int
     nsearch, int *searchlist, int *neighborlist )
     IN npart: number of particles.
     IN *x: 3*npart element array of particle positions. The x, y, and
       z coordinates of particle n are 3*n, 3*n+1, and 3*n+2.
     IN nneighbor: number of neighbors of each particle to find.
     IN nsearch: number of particles in the searchlist (see below).
     IN *searchlist: an array listing the particles whose neighbors
       are to be found.
     OUT *neighborlist: a pointer to a nsearch x nneighbor block of
       integers. The index of the jth nearnest neighbor of particle i,
       where i is the index of the particle in the search list, is
       returned in element i*nneighbor+j of the array.

   void FindNeighbor ( int npart, Real *x, Real rad, int
   ***neighborlist, int *nneighbor )
     IN npart: number of particles.
     IN *x: 3*npart element array of particle positions. The x, y, and
       z coordinates of particle n are 3*n, 3*n+1, and 3*n+2.
     IN rad: the distance within which particles must lie to be
       considered neighbors.
     OUT **neighborlist: a pointer to an array with npart
       elements. The nth element of this array is itself an array
       listing the neighbors of the nth particle. The length of this
       array is specified by nneighbor (see below). This pointer
       should not point to valid memory when passed.
     OUT *nneighbor: an array of npart elements listing the number of
       neighbors of each particle. This array must point to valid
       memory when passed.

   void FindNeighborPartial ( int npart, Real *x, int nsearch, int
     *searchlist, Real *searchrad, int ***neighborlist, int
     *nneighbor )
     IN npart: number of particles.
     IN *x: 3*npart element array of particle positions. The x, y, and
       z coordinates of particle n are 3*n, 3*n+1, and 3*n+2.
     IN nsearch: the number of particles whose neighbors are to be
       found
     IN *searchlist: an array giving the indices of the particles for
       whose neighbors the routine should search
     IN *searchrad: an nsearch element array whose nth element is the
       maximum distance from the nth particle in the search list that
       another particle can be and still be considered a neighbor
     OUT **neighborlist: a pointer to an array with npart
       elements. The nth element of this array is itself an array
       listing the neighbors of the nth particle in the search
       list. The length of this array is specified by nneighbor (see
       below). This pointer should not point to valid memory when
       passed.
     OUT *nneighbor: an array of npart elements listing the number of
       neighbors of each particle. This array must point to valid
       memory when passed.
     */

#include "preincludes.h"
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "FofLib.h"

/**************************************/
/* Macros, Utilities, and Definitions */
/**************************************/

struct treenode {
  int nidx, *idxlist;
  /* Number of particles and particle list pointed to by this
     node. nidx is non-zero only for leaves. */
  FLOAT xmin[3], xmax[3];
  /* Corners of box bounding region occupied by this node and its
     children. */
  int splitdim;
  /* Dimension along which this node splits. -1 for a leaf. */
};

#define ROOT 1
#define LEFT(i) (i<<1)
#define RIGHT(i) ((i<<1)+1)
#define PARENT(i) (i>>1)
#define SIBLING(i) ((i&1)?i-1:i+1)
#define SETNEXT(i) { while (i&1) i=i>>1; ++i; }

static FLOAT maxarg1,maxarg2;
#define FMAX(a,b) (maxarg1=(a),maxarg2=(b), (maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))

static FLOAT minarg1,minarg2;
#define FMIN(a,b) (minarg1=(a),minarg2=(b), (minarg1) < (minarg2) ? (minarg1) : (minarg2))

static int imaxarg1,imaxarg2;
#define IMAX(a,b) (imaxarg1=(a),imaxarg2=(b), (imaxarg1) > (imaxarg2) ? (imaxarg1) : (imaxarg2))

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b), (iminarg1) < (iminarg2) ? (iminarg1) : (iminarg2))

static FLOAT swaptmp;
#define SWAP(a,b) swaptmp=(a),(a)=(b),(b)=swaptmp;

static int iswaptmp;
#define ISWAP(a,b) iswaptmp=(a),(a)=(b),(b)=iswaptmp;

//static FLOAT sqrarg;
//#define SQR(a) ((sqrarg=(a)) == 0.0 ? 0.0 : sqrarg*sqrarg)
#define SQR(a) ((a)*(a))

static FLOAT *x1,*x2;
#define DISTSQR(a,b) (x1=(a),x2=(b), (x1)==(x2) ? 0.0 : SQR(x1[0]-x2[0])+SQR(x1[1]-x2[1])+SQR(x1[2]-x2[2]))

/*******************/
/* LOCAL FUNCTIONS */
/*******************/

static void ErrorHandler(char *errstr) {
  /* Error handler */
  fprintf(stderr, "Fatal error in Foflib.c: %s.\n", errstr);
  fprintf(stderr, "Exiting to system.\n");
  exit(1);
}

/* Numerical Recipes quicksort routine */
#define M 7
#define NSTACK 50
static void QuickSort(unsigned long n, FLOAT arr[])
{
  unsigned long i,ir=n,j,k,l=1,*istack;
  int jstack=0;
  FLOAT a,temp;

  if (!(istack=(unsigned long *) calloc(NSTACK,sizeof(long))))
    ErrorHandler("unable to allocate workspace in QuickSort");
  istack=istack-1;
  for (;;) {
    if (ir-l < M) {
      for (j=l+1;j<=ir;j++) {
	a=arr[j];
	for (i=j-1;i>=l;i--) {
	  if (arr[i] <= a) break;
	  arr[i+1]=arr[i];
	}
	arr[i+1]=a;
      }
      if (jstack == 0) break;
      ir=istack[jstack--];
      l=istack[jstack--];
    } else {
      k=(l+ir) >> 1;
      SWAP(arr[k],arr[l+1])
	if (arr[l] > arr[ir]) {
	  SWAP(arr[l],arr[ir])
	    }
      if (arr[l+1] > arr[ir]) {
	SWAP(arr[l+1],arr[ir])
	  }
      if (arr[l] > arr[l+1]) {
	SWAP(arr[l],arr[l+1])
	  }
      i=l+1;
      j=ir;
      a=arr[l+1];
      for (;;) {
	do i++; while (arr[i] < a);
	do j--; while (arr[j] > a);
	if (j < i) break;
	SWAP(arr[i],arr[j]);
      }
      arr[l+1]=arr[j];
      arr[j]=a;
      jstack += 2;
      if (jstack > NSTACK)
	ErrorHandler("stack size too small in QuickSort");
      if (ir-i+1 >= j-l) {
	istack[jstack]=ir;
	istack[jstack-1]=i;
	ir=j-1;
      } else {
	istack[jstack]=j-1;
	istack[jstack-1]=l;
	l=i;
      }
    }
  }
  free(istack+1);
}
#undef M
#undef NSTACK

static void ParticleSort(struct treenode node, FLOAT *x) {
  /* Routine to sort particles so half are to the left and half are to
     the right of a certain position in the specified dimension. The
     particle at exactly that position will wind up in the middle. */
  int leftend, rightend, leftptr, rightptr, middle;
  FLOAT splitval;
  int n;

  /* Initialize */
  leftend=0;
  rightend=node.nidx-1;
  middle=(node.nidx-1)>>1;

  /* Sort until separation is complete */
  while (leftend<rightend) {

    /* Initialize pointers */
    leftptr=leftend;
    rightptr=rightend;
    /* Choose partition value */
    splitval=x[3*node.idxlist[(leftptr+rightptr)/2]+node.splitdim];

    /* Put < partition value left, > partition value right */
    while (1) {
      while ((x[3*node.idxlist[leftptr]+node.splitdim]<splitval) &&
	     (leftptr<rightptr)) leftptr++;
      while ((x[3*node.idxlist[rightptr]+node.splitdim]>splitval) &&
	     (leftptr<=rightptr) && (rightptr!=0)) rightptr--;
      if (leftptr>=rightptr) break;
      ISWAP(node.idxlist[leftptr], node.idxlist[rightptr]);
      leftptr++;
      rightptr--;
    }

    /* Set rightptr to point to first element larger than partition
       value. */
    if ((x[3*node.idxlist[rightptr]+node.splitdim]>splitval) &&
	(leftptr<=rightptr)) rightptr--;
    rightptr++;
    /* Set leftend and rightend to point to enclose the part of the
       list containing the middle index. If the middle index is the
       only thing contained, then we're done. */
    if (rightptr>middle) rightend=rightptr-1;
    else leftend=rightptr;
  }
}

static struct treenode *BuildTree(FLOAT *x, int nidx, int *idxlist,
				  int leafsize) {
  /* Routine to build a binary tree of particle positions. */
  int i, n, nodes, leaves, curnode;
  int splitptr=0, splitidx;
  FLOAT width;
  struct treenode *tree;

  /* First figure out how many nodes we need and allocate */
  for (n=nidx, nodes=1, leaves=1; n>=leafsize; n=n>>1) {
    leaves=leaves<<1;
    nodes+=leaves;
  }
  if (!(tree=(struct treenode *) calloc(nodes, sizeof(struct treenode))))
    ErrorHandler("unable to allocate memory for particle tree");
  tree--; /* We want tree[1] to be first node */

  /* Initialize root node by setting bounding box */
  curnode=ROOT;
  tree[curnode].nidx=nidx; /* Record number of particles */
  tree[curnode].idxlist=idxlist; /* Record particle indices */
  for (i=0; i<3; i++) /* Set bounding box */
    tree[curnode].xmin[i]=tree[curnode].xmax[i]=x[i];
  for (n=1; n<nidx; n++) {
    for (i=0; i<3; i++) {
      tree[curnode].xmin[i]=FMIN(tree[curnode].xmin[i], x[3*n+i]);
      tree[curnode].xmax[i]=FMAX(tree[curnode].xmax[i], x[3*n+i]);
    }
  }

  while (1) {
      
    /* See if we should be a leaf or a parent */
    if (tree[curnode].nidx <= leafsize) { /* We are a leaf */
      tree[curnode].splitdim=-1; /* No split */
      SETNEXT(curnode); /* Move to next node */
    } else { /* We are a parent */
      /* Figure out which dimension to split along */
      width=tree[curnode].xmax[0]-tree[curnode].xmin[0];
      tree[curnode].splitdim=0;
      for (i=1; i<3; i++) {
	if (tree[curnode].xmax[i]-tree[curnode].xmin[i] > width) {
	  width=tree[curnode].xmax[i]-tree[curnode].xmin[i];
	  tree[curnode].splitdim=i;
	}
      }
      /* Put half of particles to the left, half to the right */
      ParticleSort(tree[curnode], x);

      /* Set bounding boxes for child nodes */
      for (i=0; i<3; i++) {
	tree[LEFT(curnode)].xmin[i]=tree[curnode].xmin[i];
	tree[LEFT(curnode)].xmax[i]=tree[curnode].xmax[i];
	tree[RIGHT(curnode)].xmin[i]=tree[curnode].xmin[i];
	tree[RIGHT(curnode)].xmax[i]=tree[curnode].xmax[i];
      }
      /* Set maximum for split dimension in left child node equal to
	 the separation value (the value of the middle node) in that
	 dimension. */
      tree[LEFT(curnode)].xmax[tree[curnode].splitdim] =
	x[3*tree[curnode].idxlist[(tree[curnode].nidx-1)/2]+
	 tree[curnode].splitdim];
      /* Set minimum for split dimension in right child node equal to
	 the separation value (the value of the middle node) in that
	 dimension. */
      tree[RIGHT(curnode)].xmin[tree[curnode].splitdim] =
	x[3*tree[curnode].idxlist[(tree[curnode].nidx-1)/2]+
	 tree[curnode].splitdim];

      /* Set particle indices for child nodes */
      tree[LEFT(curnode)].nidx=(tree[curnode].nidx+1)/2;
      tree[RIGHT(curnode)].nidx=tree[curnode].nidx/2;
      tree[LEFT(curnode)].idxlist=tree[curnode].idxlist;
      tree[RIGHT(curnode)].idxlist=tree[curnode].idxlist+
	(tree[curnode].nidx+1)/2;

      /* Pass to next node */
      curnode=LEFT(curnode);
    }

    /* Check if we're done */
    if (curnode==ROOT) break;
  }

  /* Return */
  return(tree);
}

static struct treenode *BuildVarTree(FLOAT *x, FLOAT *link, int
				     nidx, int *idxlist, int leafsize)
{
  /* Routine to build a binary tree of particle positions. */
  int i, n, nodes, leaves, curnode;
  int splitptr=0, splitidx;
  FLOAT width;
  FLOAT leftsplit, rightsplit;
  struct treenode *tree;

  /* First figure out how many nodes we need and allocate */
  for (n=nidx, nodes=1, leaves=1; n>=leafsize; n=n>>1) {
    leaves=leaves<<1;
    nodes+=leaves;
  }
  if (!(tree=(struct treenode *) calloc(nodes, sizeof(struct treenode))))
    ErrorHandler("unable to allocate memory for particle tree");
  tree--; /* We want tree[1] to be first node */

  /* Initialize root node by setting bounding box */
  curnode=ROOT;
  tree[curnode].nidx=nidx; /* Record number of particles */
  tree[curnode].idxlist=idxlist; /* Record particle indices */
  for (i=0; i<3; i++) { /* Set bounding box */
    tree[curnode].xmin[i]=x[i]-link[0];
    tree[curnode].xmax[i]=x[i]+link[0];
  }
  for (n=1; n<nidx; n++) {
    for (i=0; i<3; i++) {
      tree[curnode].xmin[i]=FMIN(tree[curnode].xmin[i],
				 x[3*n+i]-link[n]);
      tree[curnode].xmax[i]=FMAX(tree[curnode].xmax[i],
				 x[3*n+i]+link[n]);
    }
  }
  
  while (1) {
    
    /* See if we should be a leaf or a parent */
    if (tree[curnode].nidx <= leafsize) { /* We are a leaf */
      tree[curnode].splitdim=-1; /* No split */
      SETNEXT(curnode); /* Move to next node */
    }
    else { /* We are a parent */
      /* Figure out which dimension to split along */
      width=tree[curnode].xmax[0]-tree[curnode].xmin[0];
      tree[curnode].splitdim=0;
      for (i=1; i<3; i++) {
	if (tree[curnode].xmax[i]-tree[curnode].xmin[i] > width) {
	  width=tree[curnode].xmax[i]-tree[curnode].xmin[i];
	  tree[curnode].splitdim=i;
	}
      }
      /* Put half of particles to the left, half to the right */
      ParticleSort(tree[curnode], x);

      /* Set bounding boxes for child nodes */
      for (i=0; i<3; i++) {
	tree[LEFT(curnode)].xmin[i]=tree[curnode].xmin[i];
	tree[LEFT(curnode)].xmax[i]=tree[curnode].xmax[i];
	tree[RIGHT(curnode)].xmin[i]=tree[curnode].xmin[i];
	tree[RIGHT(curnode)].xmax[i]=tree[curnode].xmax[i];
      }
      /* Figure out new bounding box in split dimension by going
	 through all particles on "left" side and checking position +
	 linking length in split dimension. */
      leftsplit=x[3*tree[curnode].idxlist[0]+tree[curnode].splitdim] +
	link[tree[curnode].idxlist[0]];
      for (i=1; i<(tree[curnode].nidx-1)/2; i++)
	leftsplit=FMAX(leftsplit,
		       x[3*tree[curnode].idxlist[0]
			+ tree[curnode].splitdim] +
		       link[tree[curnode].idxlist[0]]);
      /* Same for the minimum bounding box edge of particles on the
	 "right" side. */
      rightsplit=x[3*tree[curnode].idxlist[0]+tree[curnode].splitdim] -
	link[tree[curnode].idxlist[0]];
      for (i=(tree[curnode].nidx-1)/2; i<tree[curnode].nidx; i++)
	rightsplit=FMIN(rightsplit,
		       x[3*tree[curnode].idxlist[0]
			+ tree[curnode].splitdim] -
		       link[tree[curnode].idxlist[0]]);
      /* Set maximum for split dimension in left child node equal to
	 the value for leftsplit we've just found. Ditto for
	 rightsplit. */
      tree[LEFT(curnode)].xmax[tree[curnode].splitdim] = leftsplit;
      tree[RIGHT(curnode)].xmin[tree[curnode].splitdim] = rightsplit;

      /* Set particle indices for child nodes */
      tree[LEFT(curnode)].nidx=(tree[curnode].nidx+1)/2;
      tree[RIGHT(curnode)].nidx=tree[curnode].nidx/2;
      tree[LEFT(curnode)].idxlist=tree[curnode].idxlist;
      tree[RIGHT(curnode)].idxlist=tree[curnode].idxlist+
	(tree[curnode].nidx+1)/2;

      /* Pass to next node */
      curnode=LEFT(curnode);
    }

    /* Check if we're done */
    if (curnode==ROOT) break;
  }

  /* Return */
  return(tree);
}

static void AddToGroup(struct treenode *tree, int rootnode, int *fifo,
		       int *fifotail, int *group, int groupnum, int
		       *groupsize) {
  /* Add all particles in this node and below to group */
  int n, curnode;

  /* Start with node passed to us */
  curnode=rootnode;
  while (1) {
    /* Check if this is a leaf or a parent */
    if (tree[curnode].splitdim==-1) { /* This is a leaf */
      for (n=0; n<tree[curnode].nidx; n++) {
	/* Skip already assigned particles */
	if (group[tree[curnode].idxlist[n]]!=-1) continue;
	group[tree[curnode].idxlist[n]]=groupnum; /* Add to group */
	fifo[++(*fifotail)]=tree[curnode].idxlist[n]; /* Add to fifo */
	(*groupsize)++; /* Add to groupsize */
      }
      SETNEXT(curnode);
      if ((curnode<=rootnode) || (curnode==rootnode+1)) return;
      else continue;
    } else { /* This is a parent */
      curnode=LEFT(curnode);
    }
  }
}

static int CheckContain(struct treenode node, FLOAT *pos, FLOAT
			 link) {
  /* Check to see if a node is entirely contained within a sphere
     centered on pos */
  FLOAT corner[3];

  corner[0]=node.xmin[0];
  corner[1]=node.xmin[1];
  corner[2]=node.xmin[2];
  if (DISTSQR(corner,pos)>link*link) return(0);
  corner[0]=node.xmin[0];
  corner[1]=node.xmin[1];
  corner[2]=node.xmax[2];
  if (DISTSQR(corner,pos)>link*link) return(0);
  corner[0]=node.xmin[0];
  corner[1]=node.xmax[1];
  corner[2]=node.xmin[2];
  if (DISTSQR(corner,pos)>link*link) return(0);
  corner[0]=node.xmax[0];
  corner[1]=node.xmin[1];
  corner[2]=node.xmin[2];
  if (DISTSQR(corner,pos)>link*link) return(0);
  corner[0]=node.xmin[0];
  corner[1]=node.xmax[1];
  corner[2]=node.xmax[2];
  if (DISTSQR(corner,pos)>link*link) return(0);
  corner[0]=node.xmax[0];
  corner[1]=node.xmin[1];
  corner[2]=node.xmax[2];
  if (DISTSQR(corner,pos)>link*link) return(0);
  corner[0]=node.xmax[0];
  corner[1]=node.xmax[1];
  corner[2]=node.xmin[2];
  if (DISTSQR(corner,pos)>link*link) return(0);
  corner[0]=node.xmax[0];
  corner[1]=node.xmax[1];
  corner[2]=node.xmax[2];
  if (DISTSQR(corner,pos)>link*link) return(0);
  return(1);
}

static void SearchTree(FLOAT link, FLOAT *x, FLOAT *pos, struct
		       treenode *tree, int *fifo, int *fifotail, int
		       *group, int groupnum, int *groupsize) {
  /* Search around a position pos for particles to add to group */
  int n;
  int curnode;

  /* Start with root node */
  curnode=ROOT;

  while (1) {
    /* First check if anything in this node could be in the linking
       radius. If not, move on. */
    if ((tree[curnode].xmin[0]-link>pos[0]) ||
	(tree[curnode].xmax[0]+link<pos[0]) ||
	(tree[curnode].xmin[1]-link>pos[1]) ||
	(tree[curnode].xmax[1]+link<pos[1]) ||
	(tree[curnode].xmin[2]-link>pos[2]) ||
	(tree[curnode].xmax[2]+link<pos[2])) {
      SETNEXT(curnode);
      if (curnode==ROOT) return; else continue;
    }

    /* Now check if node is entirely contained in linking radius. If
       so, add all particles in all leaves of this node to group. */
    if (CheckContain(tree[curnode],pos,link)) {
      /* This node is entirely contained, so add all */
      AddToGroup(tree, curnode, fifo, fifotail, group, groupnum,
		 groupsize);
      SETNEXT(curnode);
      if (curnode==ROOT) return; else continue;
    }

    /* If we're here, box for the node is neither contained by nor
       disjoint from the linking sphere. If this node is a leaf, we
       now check all its members. If not, we move on to this node's
       children. */
    if (tree[curnode].splitdim==-1) { /* This is a leaf */
      for (n=0; n<tree[curnode].nidx; n++) {
	/* Skip already-assigned particles */
	if (group[tree[curnode].idxlist[n]]!=-1) continue;
	/* Check if particle is in linking radius */
	if (DISTSQR(pos,x+3*tree[curnode].idxlist[n])<link*link) {
	  group[tree[curnode].idxlist[n]]=groupnum; /* Add to group */
	  fifo[++(*fifotail)]=tree[curnode].idxlist[n]; /* Add to fifo */
	  (*groupsize)++; /* Add to groupsize */
	}
      }
      SETNEXT(curnode);
      if (curnode==ROOT) return; else continue;
    } else { /* This is a parent */
      curnode=LEFT(curnode);
    }
  }
}

static void SearchVarTree(FLOAT *link, FLOAT *x, FLOAT *pos, FLOAT
			  thislink, struct treenode *tree, int *fifo,
			  int *fifotail, int *group, int groupnum, int
			  *groupsize) {
  /* Search around a position pos with linking radius thislink for
     particles to add to group */
  int n;
  int curnode;

  /* Start with root node */
  curnode=ROOT;

  while (1) {
    /* First check if anything in this node could be in the linking
       radius. If not, move on. */
    if ((tree[curnode].xmin[0]-thislink>pos[0]) ||
	(tree[curnode].xmax[0]+thislink<pos[0]) ||
	(tree[curnode].xmin[1]-thislink>pos[1]) ||
	(tree[curnode].xmax[1]+thislink<pos[1]) ||
	(tree[curnode].xmin[2]-thislink>pos[2]) ||
	(tree[curnode].xmax[2]+thislink<pos[2])) {
      SETNEXT(curnode);
      if (curnode==ROOT) return; else continue;
    }

    /* Now check if node is entirely contained in linking radius. If
       so, add all particles in all leaves of this node to group. */
    if (CheckContain(tree[curnode],pos,thislink)) {
      /* This node is entirely contained, so add all */
      AddToGroup(tree, curnode, fifo, fifotail, group, groupnum,
		 groupsize);
      SETNEXT(curnode);
      if (curnode==ROOT) return; else continue;
    }

    /* If we're here, box for the node is neither contained by nor
       disjoint from the linking sphere. If this node is a leaf, we
       now check all its members. If not, we move on to this node's
       children. */
    if (tree[curnode].splitdim==-1) { /* This is a leaf */
      for (n=0; n<tree[curnode].nidx; n++) {
	/* Skip already-assigned particles */
	if (group[tree[curnode].idxlist[n]]!=-1) continue;
	/* Check if particle is in linking radius */
	if (DISTSQR(pos,x+3*tree[curnode].idxlist[n]) <
	    SQR(FMAX(thislink, link[tree[curnode].idxlist[n]]))) {
	  group[tree[curnode].idxlist[n]]=groupnum; /* Add to group */
	  fifo[++(*fifotail)]=tree[curnode].idxlist[n]; /* Add to fifo */
	  (*groupsize)++; /* Add to groupsize */
	}
      }
      SETNEXT(curnode);
      if (curnode==ROOT) return; else continue;
    } else { /* This is a parent */
      curnode=LEFT(curnode);
    }
  }
}

FLOAT MinBoxDistSqr(FLOAT *x, FLOAT *boxmin, FLOAT *boxmax) {
  /* Return the minimum squared distance between point x and any point
     in the box whose minimum and maximum corners are given by boxmin
     and boxmax, respectively. x, boxmin, and boxmax are all 3 element
     arrays. */
  FLOAT sqrdist=0.0;
  int n;

  for (n=0; n<3; n++) {
    if (x[n]<boxmin[n]) {
      sqrdist+=SQR(x[n]-boxmin[n]);
    } else if (x[n]>boxmax[n]) {
      sqrdist+=SQR(boxmax[n]-x[n]);
    }
  }
  return(sqrdist);
}

#define HUGE 1.0e33
static void NeighborSearchTree(int nneighbor, FLOAT *x, int thisidx,
			       struct treenode *tree, int *neighbors)
{
  /* Routine to find the nneighbor particles closest to the particle
     with index thisidx. The results are returned in *neighbors. */
  int n, i;
  int curnode, tempnode;
  int nload, nload0;
  FLOAT *sqrdistances, sqrdist;

  /* The basic algorithm is to traverse the tree, from left to
     right. At each node, we check to see whether the smallest
     distance between our base particle and the bounding box of that
     node is smaller than the biggest distance currently on our list
     of neighbors. If it isn't there is no need to check that node,
     and we proceed to the next node. When we reach a leaf, we check
     all the particles in that leaf to see if they belong in the
     neighbors list. */

  /* Allocate workspace. */
  if (!(sqrdistances=(FLOAT *) calloc(nneighbor, sizeof(FLOAT))))
    ErrorHandler("unable to allocate workspace in NeighborSearch()");

  /* It's important to make good initial guesses for the locations of
     neighbors -- if not, you wind up having to traverse the entire
     tree for certain particles. For a good starting guess, traverse
     the tree to find the leaf in which the base particle resides and
     use its leafmates as initial guesses. */

  /* Traverse tree down to base particle's leaf */
  curnode=ROOT;
  while (tree[curnode].splitdim!=-1) {
    /* Determine if particle is to left or right of split, and go in
       that direction down the tree */
    if (tree[LEFT(curnode)].xmax[tree[curnode].splitdim] >=
	x[3*thisidx+tree[curnode].splitdim])
      curnode=LEFT(curnode); else curnode=RIGHT(curnode);
  }
  /* Load leafmates into neighbors and distances list */
  for (n=nload=0; (n<tree[curnode].nidx) && (nload<nneighbor); n++) {
    /* Don't count base particle as its own neighbor */
    if (tree[curnode].idxlist[n]==thisidx) continue;
    /* Load square distance lists. */
    sqrdistances[nload]= 
      DISTSQR(x+3*thisidx, x+3*tree[curnode].idxlist[n]);
    nload++;
  }
  /* If there weren't enough particles at this leaf to fill all the
     positions in the neighbors list, we need to go to another node to
     get some more. Try moving to the leaf to the left of this one. If
     we can't do that, just fill in the remaining values with
     dummies. That's ok, because we're going to traverse the left side
     of the tree first anyway and won't get screwed in this case. */
  while (nload<nneighbor) {
    /* See if we can move left */
    for (tempnode=curnode; !(tempnode & 1); tempnode=tempnode>>1);
    if (tempnode==ROOT) { /* We can't move left, we're at the leftmost
			     leaf of the tree */
      /* Fill in remaining distances with dummy values. */
      for (; nload<nneighbor; nload++) sqrdistances[nload]=HUGE;
      break;
    }
    /* We can move to the left neighbor leaf, so do so */
    curnode--;
    /* Load the particles from this leaf */
    for (n=nload, nload0=nload; n<IMIN(nload0+tree[curnode].nidx,
				       nneighbor); n++) {
      neighbors[n]=tree[curnode].idxlist[n-nload0];
      sqrdistances[n]=
	DISTSQR(x+3*thisidx, x+3*tree[curnode].idxlist[n-nload0]);
      nload++;
    }
  }
  /* We've now loaded initial values into the neighbors and distance
     to neighbors lists. Now sort these lists to put the closest
     particles at the beginning and the furthest at the end. */
  QuickSort(nneighbor, sqrdistances-1);

  /* Final step in initializations -- remove all the neighbor
     particles from the neighbor list, and put fill all distances with
     1.1 * the distance to the most distant particle in the list. That
     way we don't have to worry about checking whether we've already
     put a certain particle into the list, but we know the distances
     we have pre-loaded will not prevent us from finding particles we
     should find. */
  for (n=0; n<nneighbor; n++) {
    sqrdistances[n] = 1.1 * sqrdistances[nneighbor-1];
    neighbors[n] = -1;
  }

  /* Now we're ready to traverse the tree starting from root */
  curnode=ROOT;

  /* Start traversing tree */
  while (1) {
    /* Check if we need to investigate this node */
    if (sqrdistances[nneighbor-1] <
	MinBoxDistSqr(x+3*thisidx, tree[curnode].xmin,
		      tree[curnode].xmax)) {
      SETNEXT(curnode);
      if (curnode==ROOT) break;
      continue;
    }
    /* If we're here, we have to investigate this node. Now check if
       it's a leaf or not. */
    if (tree[curnode].splitdim!=-1) {
      /* We're not a leaf, so go to left child */
      curnode=LEFT(curnode);
      continue;
    }
    /* If we're here, we're a leaf, so check all particles in leaf */
    for (n=0; n<tree[curnode].nidx; n++) {
      /* Don't check the base particle against itself */
      if (thisidx==tree[curnode].idxlist[n]) continue;
      /* Get square of distance between particle and base particle */
      sqrdist = DISTSQR(x+3*thisidx, x+3*tree[curnode].idxlist[n]);
      if (sqrdist < sqrdistances[nneighbor-1]) {
	/* We've found a particle to be added to list */
	for (i=nneighbor-1; i>0; i--) {
	  /* Make room in list for new particle */
	  if (sqrdist>sqrdistances[i-1]) break;
	  sqrdistances[i]=sqrdistances[i-1];
	  neighbors[i]=neighbors[i-1];
	}	  
	/* Insert new particle */
	sqrdistances[i]=sqrdist;
	neighbors[i]=tree[curnode].idxlist[n];
      }
    }
    /* We're done with this leaf, so move to next node */
    SETNEXT(curnode);
    if (curnode==ROOT) break;
  }

  /* Clean up memory */
  free(sqrdistances);
}
#undef HUGE

static void AddToNeighborList(struct treenode *tree, int thisidx, int
			      rootnode, int *neighborlist, int
			      *nneighbor) {
  /* Add all particles at or beneath the current node to the neighbor
     list. */
  int n, curnode;

  /* Start with node passed to us */
  curnode=rootnode;
  while (1) {
    /* Check if this is a leaf or a parent */
    if (tree[curnode].splitdim==-1) { /* This is a leaf */
      for (n=0; n<tree[curnode].nidx; n++) {
	/* Don't count the particle as its own neighbor */
	if (tree[curnode].idxlist[n]==thisidx) continue;
	neighborlist[*nneighbor]=tree[curnode].idxlist[n];
	(*nneighbor)++;
      }
      SETNEXT(curnode);
      if ((curnode<=rootnode) || (curnode==rootnode+1)) return;
      else continue;
    } else { /* This is a parent */
      curnode=LEFT(curnode);
    }
  }
}


static int FindNeighborTree(FLOAT *x, int thisidx, FLOAT rad, struct
			    treenode *tree, int *neighborlist) {
  /* Search for particles within distance rad of particle with index
     thisidx. Add these particles to the neighborlist. */
  int n, nneighbor=0;
  int curnode;

  /* Start with root node */
  curnode=ROOT;

  while (1) {
    /* First check if anything in this node could be in the neighbor
       radius. If not, move on. */
    if ((tree[curnode].xmin[0]-rad>x[3*thisidx]) ||
	(tree[curnode].xmax[0]+rad<x[3*thisidx]) ||
	(tree[curnode].xmin[1]-rad>x[3*thisidx+1]) ||
	(tree[curnode].xmax[1]+rad<x[3*thisidx+1]) ||
	(tree[curnode].xmin[2]-rad>x[3*thisidx+2]) ||
	(tree[curnode].xmax[2]+rad<x[3*thisidx+2])) {
      SETNEXT(curnode);
      if (curnode==ROOT) break; else continue;
    }

    /* Now check if node is entirely contained in neighbor radius. If
       so, add all particles in all leaves of this node to neighbor
       list. */
    if (CheckContain(tree[curnode],x+3*thisidx,rad)) {
      /* This node is entirely contained, so add all */
      AddToNeighborList(tree,thisidx,curnode,neighborlist,&nneighbor);
      SETNEXT(curnode);
      if (curnode==ROOT) break; else continue;
    }

    /* If we're here, box for the node is neither contained by nor
       disjoint from the linking sphere. If this node is a leaf, we
       now check all its members. If not, we move on to this node's
       children. */
    if (tree[curnode].splitdim==-1) { /* This is a leaf */
      for (n=0; n<tree[curnode].nidx; n++) {
	/* Don't count the particles as it's own neighbor */
	if (tree[curnode].idxlist[n]==thisidx) continue;
	/* Check if particle is in neighbor radius */
	if (DISTSQR(x+3*thisidx,x+3*tree[curnode].idxlist[n]) <
	    rad*rad) {
	  neighborlist[nneighbor]=tree[curnode].idxlist[n];
	  nneighbor++;
	}
      }
      SETNEXT(curnode);
      if (curnode==ROOT) break; else continue;
    } else { /* This is a parent */
      curnode=LEFT(curnode);
    }
  }

  /* Return number of neighbors */
  return(nneighbor);
}
  

/********************/
/* PUBLIC FUNCTIONS */
/********************/

#define LEAFSIZE 8
int FofVar(int npart, FLOAT *x, FLOAT *link, int *group, int
	**groupsize) {
  /* Do friends of friends algorithm with variable linking length list
     link. */
  int *fifo, fifohead, fifotail, *idxlist;
  int groupnum=0;
  int n, m;
  FLOAT xmin[3], xmax[3];
  struct treenode *root;
  
  /* First allocate workspace and maximum size for output array
     groupsize */
  if (!(fifo=(int *) calloc(npart,sizeof(int))))
    ErrorHandler("unable to allocate workspace in Fof");
  if (!(idxlist=(int *) calloc(npart,sizeof(int))))
    ErrorHandler("unable to allocate workspace in Fof");
  if (!(*groupsize=(int *) calloc(npart,sizeof(int))))
    ErrorHandler("unable to allocate workspace in Fof");
  
  /* Initialize group list array and index list */
  for (n=0; n<npart; n++) {
    group[n]=-1;
    idxlist[n]=n;
  }
  
  /* Build tree of particle positions */
  root=BuildVarTree(x, link, npart, idxlist, LEAFSIZE);

  /* Now proceed through particle list */
  for (n=0; n<npart; n++) {
    if (group[n]!=-1) continue; /* Already in a group ? */
    else group[n]=groupnum; /* Assign groupnum to particle */
    (*groupsize)[groupnum]=1; /* Initial size of group */
    fifo[0]=n; /* Load first particle into fifo */
    /* Start with fifo head and tail pointers both at 0. In each step,
       find unassigned particles within linking length of particle
       whose index is fifo[fifohead]. Add them at the end of fifo and
       increment fifotail, because they need to be checked too. Once
       fifohead passes fifotail, we've checked all particles in this
       group and we're done. */
    for (fifohead=fifotail=0; fifohead<=fifotail; fifohead++)
      /* Search for friends and add to fifo */
      SearchVarTree(link, x, x+3*fifo[fifohead], link[fifo[fifohead]],
		    root, fifo, &fifotail, group, groupnum,
		    *groupsize+groupnum);
    groupnum++; /* Increment group number */
  }
  
  /* Free up unneeded memory */
  free(fifo);
  *groupsize=(int *) realloc(*groupsize, groupnum*sizeof(int));
  free(idxlist);
  free(root+1);

  return(groupnum); /* Return number of groups */
}

int FofVarList(int npart, FLOAT *x, FLOAT *link, int *group, int
	**groupsize, int ***grouplist) {
  int fifohead, fifotail, *idxlist;
  int groupnum=0;
  int n, m;
  struct treenode *root;
  
  /* First allocate workspace and maximum size for output array
     groupsize */
  if (!(idxlist=(int *) calloc(npart,sizeof(int))))
    ErrorHandler("unable to allocate workspace in FofList");
  if (!(*groupsize=(int *) calloc(npart,sizeof(int))))
    ErrorHandler("unable to allocate workspace in FofList");
  if (!(*grouplist=(int **) calloc(npart,sizeof(int *))))
    ErrorHandler("unable to allocate workspace in FofList");
  
  /* Initialize group list array */
  for (n=0; n<npart; n++) {
    group[n]=-1;
    idxlist[n]=n;
  }
  
  /* Build tree of particle positions */
  root=BuildVarTree(x, link, npart, idxlist, LEAFSIZE);

  /* Now proceed through particle list */
  for (n=0; n<npart; n++) {
    if (group[n]!=-1) continue; /* Already in a group ? */
    else group[n]=groupnum; /* Assign groupnum to particle */
    /* Allocate space for particle list */
    if (!((*grouplist)[groupnum]=(int *) calloc(npart,sizeof(int))))
      ErrorHandler("unable to allocate workspace in FofList");
    (*groupsize)[groupnum]=1; /* Initial size of group */
    (*grouplist)[groupnum][0]=n; /* Load first particle into fifo */
    /* Start with fifo head and tail pointers both at 0. In each step,
       find unassigned particles within linking length of particle
       whose index is fifo[fifohead]. Add them at the end of fifo and
       increment fifotail, because they need to be checked too. Once
       fifohead passes fifotail, we've checked all particles in this
       group and we're done. */
    for (fifohead=fifotail=0; fifohead<=fifotail; fifohead++)
      /* Search for friends and add to fifo */
      SearchVarTree(link, x, x+3*(*grouplist)[groupnum][fifohead],
		 link[(*grouplist)[groupnum][fifohead]], root,
		 (*grouplist)[groupnum], &fifotail, group, groupnum,
		 *groupsize+groupnum);
    /* Free unused parts of group list */
    (*grouplist)[groupnum]=(int *) realloc((*grouplist)[groupnum],
				   fifotail*sizeof(int));
    groupnum++; /* Increment group number */
  }
  
  /* Free up unneeded memory */
  *groupsize=(int *) realloc(*groupsize, groupnum*sizeof(int));
  *grouplist=(int **) realloc(*grouplist, groupnum*sizeof(int *));
  free(idxlist);
  free(root+1);

  return(groupnum); /* Return number of groups */
}

int Fof(int npart, FLOAT *x, FLOAT link, int *group, int **groupsize)
{
  int *fifo, fifohead, fifotail, *idxlist;
  int groupnum=0;
  int n, m;
  FLOAT xmin[3], xmax[3];
  struct treenode *root;
  
  /* First allocate workspace and maximum size for output array
     groupsize */
  if (!(fifo=(int *) calloc(npart,sizeof(int))))
    ErrorHandler("unable to allocate workspace in Fof");
  if (!(idxlist=(int *) calloc(npart,sizeof(int))))
    ErrorHandler("unable to allocate workspace in Fof");
  if (!(*groupsize=(int *) calloc(npart,sizeof(int))))
    ErrorHandler("unable to allocate workspace in Fof");
  
  /* Initialize group list array and index list */
  for (n=0; n<npart; n++) {
    group[n]=-1;
    idxlist[n]=n;
  }
  
  /* Build tree of particle positions */
  root=BuildTree(x, npart, idxlist, LEAFSIZE);

  /* Now proceed through particle list */
  for (n=0; n<npart; n++) {
    if (group[n]!=-1) continue; /* Already in a group ? */
    else group[n]=groupnum; /* Assign groupnum to particle */
    (*groupsize)[groupnum]=1; /* Initial size of group */
    fifo[0]=n; /* Load first particle into fifo */
    /* Start with fifo head and tail pointers both at 0. In each step,
       find unassigned particles within linking length of particle
       whose index is fifo[fifohead]. Add them at the end of fifo and
       increment fifotail, because they need to be checked too. Once
       fifohead passes fifotail, we've checked all particles in this
       group and we're done. */
    for (fifohead=fifotail=0; fifohead<=fifotail; fifohead++)
      /* Search for friends and add to fifo */
      SearchTree(link, x, x+3*fifo[fifohead], root, fifo, &fifotail,
		 group, groupnum, *groupsize+groupnum);
    groupnum++; /* Increment group number */
  }
  
  /* Free up unneeded memory */
  free(fifo);
  *groupsize=(int *) realloc(*groupsize, groupnum*sizeof(int));
  free(idxlist);
  free(root+1);

  return(groupnum); /* Return number of groups */
}

int FofList(int npart, FLOAT *x, FLOAT link, int *group, int
	    **groupsize, int ***grouplist)
{
  int fifohead, fifotail, *idxlist;
  int groupnum=0;
  int n, m;
  struct treenode *root;
  
  /* First allocate workspace and maximum size for output array
     groupsize */
  if (!(idxlist=(int *) calloc(npart,sizeof(int))))
    ErrorHandler("unable to allocate workspace in FofList");
  if (!(*groupsize=(int *) calloc(npart,sizeof(int))))
    ErrorHandler("unable to allocate workspace in FofList");
  if (!(*grouplist=(int **) calloc(npart,sizeof(int *))))
    ErrorHandler("unable to allocate workspace in FofList");
  
  /* Initialize group list array */
  for (n=0; n<npart; n++) {
    group[n]=-1;
    idxlist[n]=n;
  }
  
  /* Build tree of particle positions */
  root=BuildTree(x, npart, idxlist, LEAFSIZE);

  /* Now proceed through particle list */
  for (n=0; n<npart; n++) {
    if (group[n]!=-1) continue; /* Already in a group ? */
    else group[n]=groupnum; /* Assign groupnum to particle */
    /* Allocate space for particle list */
    if (!((*grouplist)[groupnum]=(int *) calloc(npart,sizeof(int))))
      ErrorHandler("unable to allocate workspace in FofList");
    (*groupsize)[groupnum]=1; /* Initial size of group */
    (*grouplist)[groupnum][0]=n; /* Load first particle into fifo */
    /* Start with fifo head and tail pointers both at 0. In each step,
       find unassigned particles within linking length of particle
       whose index is fifo[fifohead]. Add them at the end of fifo and
       increment fifotail, because they need to be checked too. Once
       fifohead passes fifotail, we've checked all particles in this
       group and we're done. */
    for (fifohead=fifotail=0; fifohead<=fifotail; fifohead++)
      /* Search for friends and add to fifo */
      SearchTree(link, x, x+3*(*grouplist)[groupnum][fifohead], root,
		 (*grouplist)[groupnum], &fifotail, group, groupnum,
		 *groupsize+groupnum);
    /* Free unused parts of group list */
    (*grouplist)[groupnum] = 
      (int *) realloc((*grouplist)[groupnum],
		      (*groupsize)[groupnum]*sizeof(int));
    groupnum++; /* Increment group number */
  }
  
  /* Free up unneeded memory */
  *groupsize=(int *) realloc(*groupsize, groupnum*sizeof(int));
  *grouplist=(int **) realloc(*grouplist, groupnum*sizeof(int *));
  free(idxlist);
  free(root+1);

  return(groupnum); /* Return number of groups */
}

int FofPrune(int npart, int ngroup, int *group, int **groupsize, int
	     minsize) {
  /* Routine to eliminate groups with fewer than minsize
     particles. Group affiliations for particles whose groups are
     eliminated is set to -1. */
  int n, newgroups;
  int *newgroup;

  /* Allocate workspace */
  if (!(newgroup=(int *) calloc(ngroup,sizeof(int))))
    ErrorHandler("unable to allocate workspace array in FofPrune");

  /* Figure out new group numbers */
  for (n=newgroups=0; n<ngroup; n++)
    if ((*groupsize)[n]<minsize) newgroup[n]=-1;
    else newgroup[n]=newgroups++;

  /* Fix particle affiliations */
  for (n=0; n<npart; n++) group[n]=newgroup[group[n]];

  /* Fix group sizes */
  for (n=0; n<ngroup; n++)
    if (newgroup[n]!=-1) (*groupsize)[newgroup[n]]=(*groupsize)[n];

  /* Free unneeded memory */
  (*groupsize)=(int *) realloc(*groupsize, newgroups*sizeof(int));
  free(newgroup);

  /* Return new number of groups */
  return(newgroups);
}

int FofListPrune(int npart, int ngroup, int *group, int **groupsize,
		 int ***grouplist, int minsize) {
  /* Routine to eliminate groups with fewer than minsize
     particles. Group affiliations for particles whose groups are
     eliminated is set to -1. */
  int n, newgroups;
  int *newgroup;

  /* Allocate workspace */
  if (!(newgroup=(int *) calloc(ngroup,sizeof(int))))
    ErrorHandler("unable to allocate workspace array in FofPrune");

  /* Figure out new group numbers */
  for (n=newgroups=0; n<ngroup; n++) {
    if ((*groupsize)[n]<minsize) newgroup[n]=-1;
    else newgroup[n]=newgroups++;
  }

  /* Fix particle affiliations */
  for (n=0; n<npart; n++) group[n]=newgroup[group[n]];

  /* Fix group sizes and lists */
  for (n=0; n<ngroup; n++)
    if (newgroup[n]!=-1) {
      (*groupsize)[newgroup[n]]=(*groupsize)[n];
      (*grouplist)[newgroup[n]]=(*grouplist)[n];
    } else free((*grouplist)[n]);

  /* Free unneeded memory */
  (*groupsize)=(int *) realloc(*groupsize,newgroups*sizeof(int));
  (*grouplist)=(int **) realloc(*grouplist,newgroups*sizeof(int *));
  free(newgroup);

  /* Return new number of groups */
  return(newgroups);
}

void NearNeighbor(int npart, FLOAT *x, int nneighbor, int
		  *neighborlist) {
  /* Routine to return the N nearest neighbors of every particle. */
  int n;
  int *idxlist;
  struct treenode *root;

  /* First allocate particle index list workspace */
  if (!(idxlist=(int *) calloc(npart+1,sizeof(int))))
    ErrorHandler("unable to allocate workspace in NearNeighbor");
  for (n=0; n<npart; n++) idxlist[n]=n;

  /* Build tree of particle positions */
  root=BuildTree(x, npart, idxlist, LEAFSIZE);

  /* Search for nearest neighbors of each particle */
  for (n=0; n<npart; n++)
    NeighborSearchTree(nneighbor, x, n, root,
		       neighborlist+nneighbor*n);

  /* Clean up memory */
  free(root+1);
  free(idxlist);
}

void NearNeighborPartial(int npart, FLOAT *x, int nneighbor, int
			 nsearch, int *searchlist, int *neighborlist) {
  /* Routine to return the N nearest neighbors of every particle in
     searchlist. */
  int n;
  int *idxlist;
  struct treenode *root;

  /* First allocate particle index list workspace */
  if (!(idxlist=(int *) calloc(npart,sizeof(int))))
    ErrorHandler("unable to allocate workspace in NearNeighborPartial");
  for (n=0; n<npart; n++) idxlist[n]=n;

  /* Build tree of particle positions */
  root=BuildTree(x, npart, idxlist, LEAFSIZE);

  /* Search for nearest neighbors of each particle in searchlist */
  for (n=0; n<nsearch; n++)
    NeighborSearchTree(nneighbor, x, searchlist[n], root,
		       neighborlist+nneighbor*n);

  /* Clean up memory */
  free(root+1);
  free(idxlist);
}

void FindNeighbor(int npart, FLOAT *x, FLOAT rad, int ***neighborlist,
		  int *nneighbor) {
  /* Routine to find all neighbors of each particle, where a neighbor
     is defined as another particle within a distance rad. */
  int n;
  int *idxlist;
  struct treenode *root;

  /* First allocate particle index list workspace */
  if (!(idxlist=(int *) calloc(npart,sizeof(int))))
    ErrorHandler("unable to allocate workspace in FindNeighbor");
  for (n=0; n<npart; n++) idxlist[n]=n;

  /* Build tree of particle positions */
  root=BuildTree(x, npart, idxlist, LEAFSIZE);

  /* Allocate array of pointers to hold number of neighbors of each
     particle. */
  if (!(*neighborlist=(int **) calloc(npart,sizeof(int *))))
    ErrorHandler("unable to allocate output array in FindNeighbor");

  /* Now go through particles finding neighbors */
  for (n=0; n<npart; n++) {
    /* Allocate enough memory to hold maximum possible number of
       neighbors */
    if (!((*neighborlist)[n]=(int *) calloc(npart,sizeof(int))))
      ErrorHandler("unable to allocate output array in FindNeighbor");
    /* Find neighbors and record how many */
    nneighbor[n]=FindNeighborTree(x, n, rad, root, (*neighborlist)[n]);
    /* Free unused parts of neighbor list */
    (*neighborlist)[n]=(int *) realloc((*neighborlist)[n],
			       nneighbor[n]*sizeof(int));
  }

  /* Clean up memory */
  free(idxlist);
  free(root+1);
}

void FindNeighborPartial(int npart, FLOAT *x, int nsearch, int
			 *searchlist, FLOAT *searchrad, int
			 ***neighborlist, int *nneighbor ) {
  /* Routine to find all neighbors of each particle in the search
     list, where a neighbor is defined as another particle within a
     distance searchrad[n] on the nth element of searchlist. */
  int n;
  int *idxlist;
  struct treenode *root;

  /* First allocate particle index list workspace */
  if (!(idxlist=(int *) calloc(npart,sizeof(int))))
    ErrorHandler("unable to allocate workspace in FindNeighbor");
  for (n=0; n<npart; n++) idxlist[n]=n;

  /* Build tree of particle positions */
  root=BuildTree(x, npart, idxlist, LEAFSIZE);

  /* Allocate array of pointers to hold number of neighbors of each
     particle. */
  if (!(*neighborlist=(int **) calloc(nsearch,sizeof(int *))))
    ErrorHandler("unable to allocate output array in FindNeighbor");

  /* Now go through particles finding neighbors */
  for (n=0; n<nsearch; n++) {
    /* Allocate enough memory to hold maximum possible number of
       neighbors */
    if (!((*neighborlist)[n]=(int *) calloc(npart,sizeof(int))))
      ErrorHandler("unable to allocate output array in FindNeighbor");
    /* Find neighbors and record how many */
    nneighbor[n]=FindNeighborTree(x, searchlist[n], searchrad[n],
				  root, (*neighborlist)[n]);
    /* Free unused parts of neighbor list */
    (*neighborlist)[n]=(int *) realloc((*neighborlist)[n],
			       nneighbor[n]*sizeof(int));
  }

  /* Clean up memory */
  free(idxlist);
  free(root+1);
}
#undef LEAFSIZE
