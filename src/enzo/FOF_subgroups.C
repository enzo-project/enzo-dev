#ifdef USE_MPI
#include <mpi.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

#include "FOF_allvars.h"
#include "FOF_forcetree.h"
#include "FOF_proto.h"

/************************************************************************/

void find_subgroups(FOFData &D)
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
  indexx(D.NumInGroup, D.Density, D.Index);


  listofdifferent = ivector(0, D.DesLinkNgb);

  for (i = D.NumInGroup, D.AnzNodes = 0, signal = 0; i >= 1; i--) {

    ngb_treefind(D.P, D.P[D.Index[i]].Pos, D.DesLinkNgb, 0, &ngblist, &r2list); 
      
    sort2_flt_int(D.DesLinkNgb, r2list-1, ngblist-1);


    /* ok, let's first see how many different groups there are */

    for (k = 0, ndiff = 0, deb = 0; k < D.DesLinkNgb; k++) {
      ind = ngblist[k];

      if (ind != D.Index[i]) {
	if(D.Density[ind] > D.Density[D.Index[i]]) {
	  deb++;

	  if (D.Head[ind]) { /* neighbor is attached to a group */
	    for (j = 0; j < ndiff; j++)
	      if (listofdifferent[j] == D.Head[ind])
		break;
		      
	    if (j >= ndiff) /* a new group has been found */
	      listofdifferent[ndiff++] = D.Head[ind];
	  }
	  else {
	    fprintf(stderr, "\nFOF: this may not occur. %g %g %g\n",
		   D.Density[ind], D.Density[D.Index[i]], 
		   D.Density[ind] - D.Density[D.Index[i]]);
	    fprintf(stderr, "\n%"ISYM" %"ISYM"\n", ind, D.Index[i] );
	    ENZO_FAIL("Error in FOF_subgroups!\n");
	  } // ENDELSE
	} // ENDIF larger density
      } // ENDIF
	  
      if (deb >= 2)
	break;
    } // ENDFOR maxima
      
      
    if (ndiff == 0) { /* this appears to be a lonely maximum -> new group */
      D.Head[D.Index[i]] = D.Tail[D.Index[i]] = D.Index[i];
      D.Len[D.Index[i]] = 1;
      D.Next[D.Index[i]]= 0;
    }

    if (ndiff == 1) { /* the particle is attached to exactly one group */
      head = listofdifferent[0];

      D.Head[D.Index[i]]   = head;
      D.Next[D.Tail[head]] = D.Index[i];
      D.Tail[head]	   = D.Index[i];
      D.Next[D.Index[i]]   = 0;
      D.Len[head]++;
    }

    if (ndiff > 1) { /* the particle merges (at least) two groups together */
      for (j = 0; j < ndiff-1; j++) { 
	head_p = listofdifferent[j];
	head_s = listofdifferent[j+1];

	if(D.Len[head_p] > D.Len[head_s]) { /* p group is longer */
	  head = head_p;
	  head_attach = head_s;
	}
	else {
	  head = head_s;
	  head_attach = head_p;
	}
	listofdifferent[j+1] = head;

	/* only in this case we bother to register the merger event */

	if (D.Len[head_s] >= D.DesLinkNgb && D.Len[head_p] >= D.DesLinkNgb) {
	  D.AnzNodes++;
	  if (D.AnzNodes >= D.MaxNodes) {
	    ENZO_VFAIL("MaxNodes=%"ISYM" reached.\n", D.MaxNodes)

	  }
		  
	  D.GroupTree[D.AnzNodes].leftlen   = D.Len[head_s];
	  D.GroupTree[D.AnzNodes].rightlen  = D.Len[head_p];
	  D.GroupTree[D.AnzNodes].lefthead  = D.Head[head_s];
	  D.GroupTree[D.AnzNodes].righthead = D.Head[head_p];

	  /* note: the usages of the Nodes to construct the break-up
	     tree is actually ambiguous in cases where there is more
	     than one subgroup at the end... this can actually
	     happen! */
		  
	  D.GroupTree[D.AnzNodes].left =  D.Node[head_s];
	  D.GroupTree[D.AnzNodes].right = D.Node[head_p];
	  D.Node[head] = D.AnzNodes;
	} // ENDIF merger event
	      
	D.Next[D.Tail[head]] = head_attach;
	D.Tail[head] = D.Tail[head_attach];
	D.Len[head] += D.Len[head_attach];
	      
	ss = head_attach;
	do {
	  D.Head[ss] = head;
	} while((ss = D.Next[ss]));
      }
	  
      D.Head[D.Index[i]]   = head;
      D.Next[D.Tail[head]] = D.Index[i];
      D.Tail[head]	   = D.Index[i];
      D.Next[D.Index[i]]   = 0;
      D.Len[head]++;
    }
  }
  
  free_ivector(listofdifferent, 0, (long) D.DesLinkNgb);
}




















