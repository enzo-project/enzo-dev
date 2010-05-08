#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "ngbtree.h"
#include "allvars.h"


static struct NODE 
{ 
  float center[3],len;   /* center and sidelength of treecubes */
  float xmin[3],xmax[3];
  int    count;           /* total # of particles in cell      */
  struct NODE *father,*suns[8];
  int    first;           /* index of first particle           */
} *nodes;

static int    numnodes;//, MaxNodes;



static int    N;

static int    *next;   /* Link-List for particle groups */


static int   *ngblist,numngb;
static float *r2list;

static float searchmin[3],searchmax[3];







float ngb_treefind(double xyz[3], int desngb, float hguess, int **ngblistback, float **r2listback)
{
#define  SECFACTOR  1.2
#ifndef  PI
#define  PI               3.1415927
#endif
  void   ngb_treesearch(struct NODE *this);
  float  selectb(unsigned long k, unsigned long n, float arr[],int ind[]);
  float  sr,sr2,h2max;  /* search radius */
  int    i,ind,ni,j,subnode,fak,k,rep=0;
  float  dx,dy,dz,r2;
  struct NODE *th,*nn;



  if(hguess>0)
    {
      sr=hguess;
    }
  else
    {
      /* determine estimate of local density */
      th=nodes;
      while(th->count>200)
	{
	  for(j=0,subnode=0,fak=1;j<3;j++,fak<<=1)
	    if(xyz[j]>th->center[j])
	      subnode+=fak;
	  
	  if(nn=th->suns[subnode])
	    if(nn->count>200)
	      th=nn;
	    else
	      break;
	  else
	    break;
	}

      sr=th->len * pow((3.0/(4*PI)*SECFACTOR)*desngb/((float)(th->count)),1.0/3);
    }


  do
    {
      for(k=0;k<3;k++)
	{
	  searchmin[k]=xyz[k]-sr;
	  searchmax[k]=xyz[k]+sr;
	}
      
      sr2=sr*sr;
      numngb=0;
      ngb_treesearch(nodes);  rep++;

      if(numngb<desngb)
	{
	  if(numngb>5)
	    {
	      sr*=pow((2.1*(float)desngb)/numngb,1.0/3);
	    }
	  else
	    {
	      sr*=2.0;
	    }
	  continue;
	}

      for(i=0;i<numngb;i++)
	{
	  ind=ngblist[i];
	  dx=P[ind].Pos[0]-xyz[0];
	  dy=P[ind].Pos[1]-xyz[1];
	  dz=P[ind].Pos[2]-xyz[2];
	  r2=dx*dx+dy*dy+dz*dz;

	  r2list[i]=r2;
	}
      
      h2max=selectb(desngb,numngb,r2list-1,ngblist-1);
      
      if(h2max<=sr2) 
	break;

      sr*=1.26;       /* 3th root of 2.0 */

      continue;
    }
  while(1);

  *ngblistback=ngblist;
  *r2listback=r2list;

  return h2max;
}






void ngb_treesearch(struct NODE *this)
{
  int k,p;
  struct NODE *nn;

  if(this->count==1)
    {
      for(k=0,p=this->first;k<3;k++)
	{
	  if(P[p].Pos[k]<searchmin[k])
	    return;
	  if(P[p].Pos[k]>searchmax[k])
	    return;
	}
      ngblist[numngb++]=p;
    }
  else
    {
      for(k=0;k<3;k++)
	{
	  if((this->xmax[k])<searchmin[k])
	    return;
	  if((this->xmin[k])>searchmax[k])
	    return;
	}

      for(k=0;k<3;k++)
	{
	  if((this->xmax[k])>searchmax[k])
	    break;
	  if((this->xmin[k])<searchmin[k])
	    break;
	}

      if(k>=3) 	  /* cell lies completely inside */
	{
	  p=this->first;
	  
	  for(k=0;k<this->count;k++)
	    {
	      ngblist[numngb++]=p;
	      p=next[p];
	    }
	}
      else
	{
	  for(k=0;k<8;k++) 
	    if(nn=this->suns[k])
	      {
		ngb_treesearch(nn);
	      }
	}
    }
}




void ngb_treeallocate(int npart,int maxnodes)  /* usually maxnodes=2*npart is suffiecient */
{
  int totbytes=0,bytes;

  MaxNodes=maxnodes;
  N=npart;
  

  if(!(nodes=malloc(bytes=MaxNodes*sizeof(struct NODE))))
    {
      fprintf(stdout,"Failed to allocate %d nodes (%d bytes).\n",MaxNodes,bytes);
      exit(0);
    }
  totbytes+= bytes;


  if(!(next=malloc(bytes=(N+1)*sizeof(int))))
    {
      fprintf(stdout,"Failed to allocate %d spaces for next array\n",N);
      exit(0);
    }
  totbytes+= bytes;

  if(!(ngblist=malloc(bytes=(N+1)*sizeof(int))))
    {
      fprintf(stdout,"Failed to allocate %d spaces for ngblist array\n",N);
      exit(0);
    }
  totbytes+= bytes;


  if(!(r2list=malloc(bytes=(N+1)*sizeof(float))))
    {
      fprintf(stdout,"Failed to allocate %d spaces for r2list array\n",N);
      exit(0);
    }
  totbytes+= bytes;

  /*
    printf("allocated %f Mbyte for ngbtree.\n", ((float)totbytes)/(1024.0*1024.0));
  */
}




void ngb_treefree(void)
{
  free(r2list);
  free(ngblist);
  free(next);
  free(nodes);
}










void ngb_treebuild(int Npart) 
/* packs the particles 0...Npart-1 in tree */
{
  int    i,j,k,subp,subi,p,ni,subnode,fak;
  float xmin[3],xmax[3],len,x;
  struct NODE *nfree,*th,*nn; 

  /*
  printf("Begin Ngb-tree construction.\n");
  */
  
  if(Npart<2)
    {
      fprintf(stdout,"must be at least two particles in tree.\n");
      MPI_Abort(MPI_COMM_WORLD, 123);
      exit(0);
    }


  N=Npart;


  for(j=0;j<3;j++)
    xmin[j]=xmax[j]=P[1].Pos[j];
  
  
  for(i=1;i<=Npart;i++)
    for(j=0;j<3;j++)
      {
	if(P[i].Pos[j]>xmax[j]) 
	  xmax[j]=P[i].Pos[j];
	if(P[i].Pos[j]<xmin[j]) 
	  xmin[j]=P[i].Pos[j];
      }
  for(j=1,len=xmax[0]-xmin[0];j<3;j++)
    if((xmax[j]-xmin[j])>len)
      len=xmax[j]-xmin[j];

  len*=1.0001;


  /* insert particle 1 in root node */

  nfree=nodes;

  for(j=0;j<3;j++)
    nfree->center[j]=(xmax[j]+xmin[j])/2;
  nfree->len=len;

  nfree->father=0;
  for(i=0;i<8;i++)
    nfree->suns[i]=0;
  nfree->first=1;

  nfree->count=1;

  numnodes=1;
  nfree++;


  next[1]=-1;


  for(i=2;i<=Npart;i++)  /* insert all other particles */
    {
      th=nodes;

      while(1)
	{
	  th->count++;

	  if(th->count==2)       /* cell was occupied with only one particle */
	    break;
	  
	  for(j=0,subnode=0,fak=1;j<3;j++,fak<<=1)
	    if(P[i].Pos[j] > th->center[j])
	      subnode+=fak;

	  if(nn=th->suns[subnode])
	    th=nn;
	  else
	    break;
	}


      if(th->count==2)  /* cell was occcupied with one particle */
	{
	  while(1)
	    {
	      p=th->first;

	      for(j=0,subp=0,fak=1;j<3;j++,fak<<=1)
		if(P[p].Pos[j]>th->center[j])
		  subp+=fak;

	      nfree->father=th;
	      for(j=0;j<8;j++)
		nfree->suns[j]=0;
	      
	      nfree->len=th->len/2;
    
	      for(j=0;j<3;j++)
		nfree->center[j]=th->center[j];

	      for(j=0;j<3;j++)
		if(P[p].Pos[j]>nfree->center[j])
		  nfree->center[j]+=nfree->len/2;
		else
		  nfree->center[j]-=nfree->len/2;

	      nfree->first=p;
	      nfree->count=1;

	      th->suns[subp]=nfree;

	      numnodes++;
	      nfree++;

	      if(numnodes>=MaxNodes)
		{
		  fprintf(stdout,"maximum node number %d in neighbour tree reached.\n",numnodes);
		  MPI_Abort(MPI_COMM_WORLD, 798798);
		  exit(0);
		}

	      for(j=0,subi=0,fak=1;j<3;j++,fak<<=1)
		if(P[i].Pos[j]>th->center[j])
		  subi+=fak;
	      
	      if(subi==subp)
		{
		  th=nfree-1;

		  th->count++;
		}
	      else
		break;
	    }

	}


      
      for(j=0,subi=0,fak=1;j<3;j++,fak<<=1)
	if(P[i].Pos[j]>th->center[j])
	  subi+=fak;
      
      nfree->father=th;


      p=th->first;
      for(j=0;j<(th->count-2);j++)
	{
	  p=next[p];
	}

      next[i]=next[p];
      next[p]=i;

  
      for(j=0;j<8;j++)
	nfree->suns[j]=0;

      nfree->len=th->len/2;
      for(j=0;j<3;j++)
	nfree->center[j]=th->center[j];

      for(j=0;j<3;j++)
	if(P[i].Pos[j]>nfree->center[j])
	  nfree->center[j]+=nfree->len/2;
	else
	  nfree->center[j]-=nfree->len/2;

      nfree->count=1;

      nfree->first=i;
      th->suns[subi]=nfree;
      
      numnodes++;
      nfree++;

      if(numnodes>=MaxNodes)
	{
	  fprintf(stdout,"maximum node number %d in neighbour tree reached.\n",numnodes);
	  MPI_Abort(MPI_COMM_WORLD, 137);
	  exit(0);
	}
    }

  for(ni=0,th=nodes;ni<numnodes;ni++,th++)
    {
      for(k=0;k<3;k++)
	{
	  th->xmin[k]=th->center[k]-th->len/2;
	  th->xmax[k]=th->center[k]+th->len/2;
	}
    }
 
  /*
  printf("Ngb-Tree contruction finished (%d nodes).\n",numnodes);
  */
}












