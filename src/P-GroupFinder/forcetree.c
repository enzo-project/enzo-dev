#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "forcetree.h"
#include "allvars.h"




static struct NODE 
{ 
  float center[3],len;                /* center and sidelength of treecubes */
  float mass,oc;                      /* mass and variable for opening criter*/
  float s[3];                         /* center of mass */
  struct NODE *next,*sibling,*father,*suns[8];
  int    partind;
} *nodes;


static struct NODE *last;



static int    numnodestotal;        /* total number of nodes */
//static int    MaxNodes;


static float xmin[3],xmax[3],len;


static  int    N;
static  int    First;

static float   knlrad  [KERNEL_TABLE+1],
               knlpot  [KERNEL_TABLE+1];




void force_treeallocate(int maxnodes) 
{
  int bytes;
  
  MaxNodes=maxnodes;

  if(!(nodes=malloc(bytes=MaxNodes*sizeof(struct NODE))))
    {
      printf("failed to allocate memory for %d tree-nodes (%d bytes).\n",MaxNodes,bytes);
      MPI_Abort(MPI_COMM_WORLD, 1212);
      exit(0);
    }
  /*
  printf("\nAllocated %g MByte for BH-tree.\n\n",bytes/(1024.0*1024.0));
  */
  force_setkernel();
}




void force_treefree(void)
{
  free(nodes);
}




void add_particle_props_to_node(struct NODE *no, int p)
{
  int i;

  for(i=0;i<3;i++)
    no->s[i] += P[p].Mass * (P[p].Pos[i] - no->center[i]);

  no->mass += P[p].Mass;
}



/* packs the particles of group 'gr' into into BH-trees */

int force_treebuild(int first, int len, float thetamax)
{
  int i,j,tr,n,ip;
  int subp,subi,p,ni,subnode,fak,fp;
  float x,length;
  float dx,dy,dz;
  struct NODE *nfree,*th,*nn,*ff;
  void force_setupnonrecursive(struct NODE *no);
  double drand48();

  First=first;
  N=len;

  nfree= nodes;  numnodestotal= 0;



  for(j=0;j<3;j++)                        /* find enclosing rectangle */
    xmin[j]=xmax[j]=P[First].Pos[j];

  for(i=1, p=First; i<=N; i++, p=Index[p])
    {
      for(j=0;j<3;j++)
	{
	  if(P[p].Pos[j]>xmax[j]) 
	    xmax[j]=P[p].Pos[j];
	  if(P[p].Pos[j]<xmin[j]) 
	    xmin[j]=P[p].Pos[j];
	}
    }
  
  for(j=1,length=xmax[0]-xmin[0];j<3;j++)  /* determine maxmimum externsion */
    if((xmax[j]-xmin[j])>length)
      length=xmax[j]-xmin[j];
  
  length*=1.01;
  

  /* insert first particle in root node */
 
  for(j=0;j<3;j++)
    nfree->center[j]=(xmax[j]+xmin[j])/2;
  nfree->len=length;
  
  nfree->father=0;
  for(i=0;i<8;i++)
    nfree->suns[i]=0;
  nfree->partind=First;

  nfree->mass= P[First].Mass;

  for(i=0;i<3;i++)
    nfree->s[i]= P[First].Mass*(P[First].Pos[i] - nfree->center[i]);
  
  nfree->sibling=0;

  numnodestotal++; nfree++;
  
  if(numnodestotal>=MaxNodes)
    {
      printf("111 maximum number %d of tree-nodes reached.\n",numnodestotal);
      MPI_Abort(MPI_COMM_WORLD, 121);
      exit(1);
    }



  
  /* insert all other particles */
  
  for(i=2, ip=Index[First]; i<=N; i++, ip=Index[ip])
    {
      th=nodes;
      
      while(1)
	{
	  add_particle_props_to_node(th, ip);

	  if(th->partind>=0)
	    break;
	  
	  for(j=0,subnode=0,fak=1;j<3;j++,fak<<=1)
	    if(P[ip].Pos[j]>th->center[j])
	      subnode+=fak;
	  
	  if(nn=th->suns[subnode])
	    th=nn;
	  else
	    break;
	}

      
      if(th->partind>=0)  /* cell is occcupied with one particle */
	{
	  while(1)
	    {
	      p=th->partind;

	      for(j=0,subp=0,fak=1;j<3;j++,fak<<=1)
		if(P[p].Pos[j]> th->center[j])
		  subp+=fak;

	      nfree->father=th;
	      
	      for(j=0;j<8;j++)
		nfree->suns[j]=0;
	      nfree->sibling=0;
	      
	      nfree->len=th->len/2;
    
	      for(j=0;j<3;j++)
		nfree->center[j]=th->center[j];

	      for(j=0;j<3;j++)
		if(P[p].Pos[j]>nfree->center[j])
		  nfree->center[j]+=nfree->len/2;
		else
		  nfree->center[j]-=nfree->len/2;

	      nfree->partind=p;

	      nfree->mass= P[p].Mass;
	      for(j=0;j<3;j++)
		nfree->s[j]=P[p].Mass* (P[p].Pos[j] - nfree->center[j]);
	      
  	      th->partind=-1;
	      th->suns[subp]=nfree;
      
	      numnodestotal++; nfree++;
	      if(numnodestotal>=MaxNodes)
		{
		  printf("222 maximum number %d of tree-nodes reached.\n",numnodestotal);
		  printf("i=%d ip=%d\n", i, ip);
		  MPI_Abort(MPI_COMM_WORLD, 121);
		  exit(1);
		}

	      for(j=0,subi=0,fak=1;j<3;j++,fak<<=1)
		if(P[ip].Pos[j]>th->center[j])
		  subi+=fak;

	      if(subi==subp)   /* the new particle lies in the same sub-cube */
		{
		  th=nfree-1;
		  add_particle_props_to_node(th,ip);		  
		}
	      else
		break;
	    }
	}
      
      
      for(j=0,subi=0,fak=1;j<3;j++,fak<<=1)
	if(P[ip].Pos[j]>th->center[j])
	  subi+=fak;
      
      nfree->father=th;
      
      for(j=0;j<8;j++)
	nfree->suns[j]=0;
      nfree->sibling=0;

      nfree->len=th->len/2;
      for(j=0;j<3;j++)
	nfree->center[j]=th->center[j];

      for(j=0;j<3;j++)
	if(P[ip].Pos[j]>nfree->center[j])
	  nfree->center[j]+=nfree->len/2;
	else
	  nfree->center[j]-=nfree->len/2;

      nfree->mass= P[ip].Mass;
      for(j=0;j<3;j++)
	nfree->s[j]= P[ip].Mass*(P[ip].Pos[j] - nfree->center[j]);

      nfree->partind=ip;
      th->suns[subi]=nfree;
      
      numnodestotal++; nfree++;

      if(numnodestotal>=MaxNodes)
	{
	  printf("333 maximum number %d of tree-nodes reached.\n",numnodestotal);
	  printf("i=%d ip=%d\n", i, ip);
	  MPI_Abort(MPI_COMM_WORLD, 333);
	  exit(1);
	}
    }

  
  
  /* now finish-up center-of-mass and quadrupol computation */
  
  for(i=0,th=nodes; i<numnodestotal; i++,th++)
    {
      for(j=0;j<3;j++)
	th->s[j] /= th->mass;
      
      if(th->partind<0)   /* cell contains more than one particle */
	{
	  dx=th->s[0];
	  dy=th->s[1];
	  dz=th->s[2];
	  
	  th->oc=sqrt(dx*dx+dy*dy+dz*dz);
	  th->oc += th->len/ (thetamax); 
	  th->oc *= th->oc;     /* used in cell-opening criterion */
	}

      th->s[0]+= th->center[0];
      th->s[1]+= th->center[1];
      th->s[2]+= th->center[2];
 
	

      for(j=7,nn=0;j>=0;j--)  	/* preparations for non-recursive walk */
	{
	  if(th->suns[j])
	    {
	      th->suns[j]->sibling=nn;
	      nn=th->suns[j];
	    }
	}
    }

  
  last=0;
  force_setupnonrecursive(nodes);  	/* set up non-recursive walk */
  last->next=0;

  
  for(i=0,th=nodes; i<numnodestotal; i++,th++)
    if(!(th->sibling))
      {
	ff=th;
	nn=ff->sibling;

	while(!nn)
	  {
	    ff=ff->father;
	    if(!ff)
	      break;
	    nn=ff->sibling;
	  }
	
	th->sibling=nn;
      }

  
  return numnodestotal;
}




void force_setupnonrecursive(struct NODE *no)
{
  int i;
  struct NODE *nn;
  
  if(last)
    last->next=no;

  last=no;
  
  for(i=0;i<8;i++)
    if(nn=no->suns[i])
      force_setupnonrecursive(nn);
}
 





void force_treeevaluate_potential(double *pos, float *pot, float epsilon)
{
  struct NODE *no,*nn;
  int i,k,p,ii;
  float r2,dx,dy,dz,r,fac,theta,u,h,ff;
  float wp;
  float r_inv;
  float h_inv;


  h = 2.8*epsilon;
  h_inv=1/h;

  no=nodes;

  *pot = 0;

  while(no)
    {
      dx=no->s[0]- pos[0];     /* observe the sign ! */
      dy=no->s[1]- pos[1];     /* this vector is -y in my thesis notation */
      dz=no->s[2]- pos[2];

      r2=dx*dx+dy*dy+dz*dz;

      if((p=no->partind)>=0)   /* single particle */
	{
	  r=sqrt(r2);  
	   
	  u=r*h_inv;

	  if(u>=1)
	    {
	      *pot -= no->mass/r;
	    }
	  else
	    {
	      ii = (int)(u*KERNEL_TABLE); ff=(u-knlrad[ii])*KERNEL_TABLE;
	      wp=knlpot[ii]+(knlpot[ii+1]-knlpot[ii])*ff;
	      
	      *pot += no->mass*h_inv*wp;
	    }
	  no=no->sibling;
	}
      else
	{
	  if(r2 < no->oc)
	    {
	      no=no->next;  /* open cell */
	    }
	  else
	    {
	      r=sqrt(r2);  
	  
	      u=r*h_inv;
	  
	      if(u>=1)  /* ordinary quadrupol moment */
		{
		  r_inv=1/r;
		  
		  *pot += -no->mass*r_inv;
		}
	      else    /* softened monopole moment */
		{
		  ii = (int)(u*KERNEL_TABLE); ff=(u-knlrad[ii])*KERNEL_TABLE;
		  wp=knlpot[ii]   + (knlpot[ii+1]-knlpot[ii])*ff;

		  *pot += no->mass*h_inv*wp;
		}
	      no=no->sibling;
	    }
	}
    }
}







void force_setkernel(void) 
{
  int i;
  float u;

  for(i=0;i<=KERNEL_TABLE;i++)
    {
      u=((float)i)/KERNEL_TABLE;

      knlrad[i] = u;

      if(u<=0.5)
	{
	  knlpot[i]=16.0/3*pow(u,2)-48.0/5*pow(u,4)+32.0/5*pow(u,5)-14.0/5;
	}
      else
	{
	  knlpot[i]=1.0/15/u +32.0/3*pow(u,2)-16.0*pow(u,3)+48.0/5*pow(u,4)-32.0/15*pow(u,5)-16.0/5;
	}
    }
}





















