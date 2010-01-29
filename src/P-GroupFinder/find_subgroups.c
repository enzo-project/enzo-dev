#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "forcetree.h"




void find_subgroups(void)
{
  float  *r2list;
  int    *ngblist;
  int    i,j,k,ind,signal;
  float  maxdens;
  int    maxindex;
  int    *listofdifferent;
  int    head, head_attach, head_s, head_p, ss;
  int    ndiff;
  int    deb;

  /* sort densities */
  indexx(NumInGroup, Density, Index);


  listofdifferent= ivector(0, DesLinkNgb);

  for(i=NumInGroup, AnzNodes=0, signal=0; i>=1; i--)
    {
      /*
      if(ThisTask==0)
	{
	  if((NumInGroup-i)> (signal/100.0)*NumInGroup)
	    {
	      if((signal%10)==0)
		printf("%d",signal);
	      else
		printf(".",signal);
	      fflush(stdout);
	      signal++;
	    }
	}
      */

      
      ngb_treefind(P[Index[i]].Pos, DesLinkNgb, 0, &ngblist, &r2list); 
      
      sort2_flt_int(DesLinkNgb, r2list-1, ngblist-1);

      
      /* ok, let's first see how many different groups there are */

      for(k=0, ndiff=0, deb=0; k<DesLinkNgb; k++)
	{
	  ind=ngblist[k];

	  if(ind != Index[i])
	    {
	      if(Density[ind] > Density[Index[i]])
		{
		  deb++;

		  if(Head[ind]) /* neighbor is attached to a group */
		    {
		      for(j=0;j<ndiff;j++)
			if(listofdifferent[j] == Head[ind])
			  break;
		      
		      if(j>=ndiff) /* a new group has been found */
			listofdifferent[ndiff++] = Head[ind];
		    }
		  else
		    {
		      printf("\nthis may not occur. %g %g %g\n",Density[ind],Density[Index[i]], Density[ind]-Density[Index[i]] );
		      printf("\n%d %d\n", ind, Index[i] );
		      MPI_Abort(MPI_COMM_WORLD, 8888);
		      exit(0); /* may not occur */
		    }
		}
	    }
	  
 	  if(deb>=2)
	    break;
	}
      
      
      if(ndiff==0)  /* this appears to be a lonely maximum -> new group */
	{
	  Head[Index[i]]= Tail[Index[i]]= Index[i];
	  Len[Index[i]] = 1;
	  Next[Index[i]]= 0;
	}

      if(ndiff==1) /* the particle is attached to exactly one group */
	{
	  head=listofdifferent[0];

	  Head[Index[i]]= head;
	  Next[Tail[head]]= Index[i];
	  Tail[head]= Index[i];
	  Len[head]++;
	  Next[Index[i]]= 0;
	}

      if(ndiff>1)  /* the particle merges (at least) two groups together */
	{
	  for(j=0; j<(ndiff-1); j++)
	    { 
	      head_p= listofdifferent[j];
	      head_s= listofdifferent[j+1];

	      if(Len[head_p] > Len[head_s]) /* p group is longer */
		{
		  head= head_p;
		  head_attach=head_s;
		}
	      else
		{
		  head= head_s;
		  head_attach=head_p;
		}
	      listofdifferent[j+1]=head;


	      if(Len[head_s] >= DesLinkNgb && Len[head_p] >= DesLinkNgb) /* only in this case we bother to register the merger event */
		{
		  AnzNodes++;
		  if(AnzNodes>=MaxNodes)
		    {
		      printf("MaxNodes=%d reached.\n", MaxNodes);
		      MPI_Abort(MPI_COMM_WORLD, 12345);
		      exit(0);
		    }
		  
		  GroupTree[AnzNodes].leftlen  = Len[head_s];
		  GroupTree[AnzNodes].rightlen = Len[head_p];
		  GroupTree[AnzNodes].lefthead  = Head[head_s];
		  GroupTree[AnzNodes].righthead = Head[head_p];

		  /* note: the usages of the Nodes to construct the break-up tree
		     is actually ambiguos in cases where there is more than one
		     subgroup at the end... this can actually happen! */
		  
		  GroupTree[AnzNodes].left=  Node[head_s];
		  GroupTree[AnzNodes].right= Node[head_p];
		  Node[head]= AnzNodes;
		}
	      
	      Next[ Tail[head] ] = head_attach;
	      Tail[head] = Tail[head_attach];
	      Len[head]+= Len[head_attach];
	      
	      ss=head_attach;
	      do
		{
		  Head[ss]= head;
		}
	      while(ss=Next[ss]);
	    }
	  
	  Head[Index[i]]= head;
	  Next[Tail[head]]= Index[i];
	  Tail[head]= Index[i];
	  Len[head]++;
	  Next[Index[i]]= 0;
	}
    }
  
  free_ivector(listofdifferent, 0, DesLinkNgb);
}




















