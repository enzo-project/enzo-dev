#ifdef USE_MPI
#include <mpi.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

#include "FOF_allvars.h"
#include "FOF_proto.h"

/************************************************************************/

void subfind(FOFData &D, int CycleNumber, FLOAT EnzoTime)
{

  int    i, k, gr, task, head, len, nsubs, offset, start=0, NSubGroupsAll=0;
  int    parent, ntot;
  char   ctype;
  float  cm[3], cmv[3], mtot, mstars;
  float  corner[3];
  FOF_particle_data *Pbuf, *partbuf;
  int    *sublen, *suboffset, *bufsublen, *bufsuboffset;

  char   *FOF_dirname = "FOF";
  char   catalogue_fname[200];
  char   particle_fname[200];
  char   halo_name[200];

  int    *TempInt;
  double *temp;

#ifdef USE_MPI
  MPI_Status status;
  MPI_Datatype IntType = (sizeof(int) == 4) ? MPI_INT : MPI_LONG_LONG_INT;
#endif

  sprintf(catalogue_fname, "%s/subgroups_%5.5d.dat", FOF_dirname, CycleNumber);
  sprintf(particle_fname, "%s/subparticles_%5.5d.dat", FOF_dirname, CycleNumber);

  if (MyProcessorNumber == ROOT_PROCESSOR) {

    if ((fd = fopen(catalogue_fname, "w")) == NULL)
      ENZO_FAIL("Unable to open FOF group file.");

    // Write header

    redshift = 1.0 / AllVars.Time - 1.0;
    fprintf(fd, "# Time     = %"PSYM"\n", EnzoTime);
    fprintf(fd, "# Redshift = %"PSYM"\n", redshift);
    //fprintf(fd, "# Number of subhalos = %"ISYM"\n", AllVars.NgroupsAll);
    fprintf(fd, "#\n");
    fprintf(fd, "# Column 1.  Subhalo number\n");
    fprintf(fd, "# Column 2.  Parent halo number\n");
    fprintf(fd, "# Column 3.  First particle in halo particle list\n");
    fprintf(fd, "#            --> All subgroup particles are consecutively listed in\n");
    fprintf(fd, "#                particle list (if written)\n");
    fprintf(fd, "# Column 4.  Number of particles\n");
    fprintf(fd, "# Column 5.  Halo mass [solar masses]\n");
    fprintf(fd, "# Column 6.  Stellar mass [solar masses]\n");
    fprintf(fd, "# Column 7.  Center of mass (x)\n");
    fprintf(fd, "# Column 8.  Center of mass (y)\n");
    fprintf(fd, "# Column 9.  Center of mass (z)\n");
    fprintf(fd, "# Column 10. Mean x-velocity [km/s]\n");
    fprintf(fd, "# Column 11. Mean y-velocity [km/s]\n");
    fprintf(fd, "# Column 12. Mean z-velocity [km/s]\n");
    fprintf(fd, "#\n");

    if (HaloFinderOutputParticleList) {

      file_id = H5Fcreate(particle_fname, H5F_ACC_TRUNC, H5P_DEFAULT);
      group_id = H5Gcreate(file_id, "/Parameters", 0);
      writeScalarAttribute(group_id, HDF5_REAL, "Redshift", &redshift);
      writeScalarAttribute(group_id, HDF5_PREC, "Time", &EnzoTime);
      H5Gclose(group_id);

    } // ENDIF output particle list

  } // ENDIF ROOT_PROCESSOR



  for (gr = D.NgroupsAll-1; gr >= 0; ) {

    for (task = 0; task < NumberOfProcessors && (gr-task) >= 0; task++) {

      if (MyProcessorNumber == task) {
	head = D.GroupDatAll[gr-task].Tag;
	len  = D.GroupDatAll[gr-task].Len;

	sublen	  = new int[len/D.DesLinkNgb];
	suboffset = new int[len/D.DesLinkNgb];
	Pbuf	  = new FOF_particle_data[len];
      } // ENDIF this task

      get_particles(task, head, len, Pbuf);

    } // ENDFOR task

    /* ok, this process got a group */

    if (gr - MyProcessorNumber >= 0) {

      len = D.GroupDatAll[gr-MyProcessorNumber].Len;

      for (k = 0; k < 3; k++)
	corner[k] = Pbuf[0].Pos[k];
	  
      for (i = 0; i < len; i++)
	for (k = 0; k < 3; k++)
	  Pbuf[i].Pos[k] = FOF_periodic(Pbuf[i].Pos[k]-corner[k], D);

      nsubs = do_subfind_in_group(D, Pbuf, len, sublen, suboffset);

      for (i = 0; i < len; i++)
	for (k = 0; k < 3; k++)
	  Pbuf[i].Pos[k] = FOF_periodic_wrap(Pbuf[i].Pos[k]+corner[k], D);

    } // ENDIF got a group

    for (task = 0; task < NumberOfProcessors && gr - task >= 0; task++) {

      parent = D.NgroupsAll-(gr-task);

      if (MyProcessorNumber == ROOT_PROCESSOR) {
	if (task == 0) {
	  for (i = 0; i < nsubs; i++) {
	    get_properties(D, Pbuf+suboffset[i], sublen[i], cm, cmv, &mtot, &mstars);

	    fprintf(fd, "%12"ISYM" %12"ISYM" %12"ISYM" %12"GOUTSYM" %12"GOUTSYM" %12"GOUTSYM" %12"GOUTSYM" %12"GOUTSYM" %12"GOUTSYM" %12"GOUTSYM" %12"GOUTSYM"\n",
		    D.NSubGroupsAll+i, PARENT, len, 
		    mtot, mstars, cm[0], cm[1], cm[2], cmv[0], cmv[1], cmv[2]);

	  } // ENDFOR subgroups

	  D.NSubGroupsAll += nsubs;
	} // ENDIF task == 0

	// task != 0
	else {
#ifdef USE_MPI
	  MPI_Recv(&nsubs, 1, MPI_INT, task, task, MPI_COMM_WORLD, &status);
	  if (nsubs > 0) {

	    bufsublen	 = new int[nsubs];
	    bufsuboffset = new int[nsubs];

	    MPI_Recv(bufsublen,    nsubs, IntType, task, task, MPI_COMM_WORLD, &status);
	    MPI_Recv(bufsuboffset, nsubs, IntType, task, task, MPI_COMM_WORLD, &status);

	    for (i = 0; i < nsubs; i++)
	      bufsuboffset[i] += start;

	  } // ENDIF subgroups

	  parent = D.NgroupsAll-(gr-task);
		  
	  partbuf = new FOF_particle_data[D.GroupDataAll[gr-task].Len];

	  MPI_Recv(partbuf, D.GroupDatAll[gr-task].Len*sizeof(FOF_particle_data), 
		   MPI_BYTE, task, task, MPI_COMM_WORLD, &status);

	  for (i = 0; i < nsubs; i++) {
	    get_properties(partbuf+bufsuboffset[i]-start, bufsublen[i], cm, 
			   cmv, &mtot, &mstars);
	    fprintf(fd, "%12"ISYM" %12"ISYM" %12"ISYM" %12"GOUTSYM" %12"GOUTSYM" %12"GOUTSYM" %12"GOUTSYM" %12"GOUTSYM" %12"GOUTSYM" %12"GOUTSYM" %12"GOUTSYM"\n",
		    D.NSubGroupsAll+i, PARENT, len, 
		    mtot, mstars, cm[0], cm[1], cm[2], cmv[0], cmv[1], cmv[2]);
	  } // ENDFOR subgroups


	  for(i=0; i<GroupDatAll[gr-task].Len; i++)
	    fwrite(&partbuf[i].Pos[0], sizeof(double), 3, fdpart);
	  for(i=0; i<GroupDatAll[gr-task].Len; i++)
	    fwrite(&partbuf[i].PartID, sizeof(int), 1, fdids);
	  for(i=0; i<GroupDatAll[gr-task].Len; i++) {
	    ctype= partbuf[i].Type;
	    fwrite(&ctype, sizeof(char), 1, fdtypes);
	  }

	  fwrite(&nsubs, sizeof(int), 1, fd);

	  start += D.GroupDatAll[gr-task].Len;
	  
	  if (nsubs > 0) {
	    delete [] bufsuboffset;
	    delete [] bufsublen;
	  }

	  delete [] partbuf;

	  D.NSubGroupsAll += nsubs; 

#endif /* USE_MPI */

	} // ENDELSE (task == 0)
      } // ENDIF ROOT_PROCESSOR

      // other processors
      else {

#ifdef USE_MPI
	if (MyProcessorNumber == task) {
	  MPI_Ssend(&nsubs, 1, MPI_INT, 0, ThisTask, MPI_COMM_WORLD);
	  if(nsubs) {
	    MPI_Ssend(sublen,    nsubs, MPI_INT, 0, ThisTask, MPI_COMM_WORLD);
	    MPI_Ssend(suboffset, nsubs, MPI_INT, 0, ThisTask, MPI_COMM_WORLD);
	  }
	  MPI_Ssend(Pbuf, GroupDatAll[gr-task].Len*sizeof(FOF_particle_data), 
		    MPI_BYTE, 0, ThisTask, MPI_COMM_WORLD);
	} // ENDIF proc == task

#endif /* USE_MPI */

      } // ENDIF other processors
    } // ENDFOR tasks

    for (task = 0; task < NumberOfProcessors && (gr-task) >= 0; task++)
      if (MyProcessorNumber == task) {
	free(Pbuf);
	free(suboffset);
	free(sublen);
      }

    gr -= task;

  } // ENDFOR groups

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fclose(fdpart);
    fclose(fdids);
    fclose(fdlen);
    fclose(fdoffset);
    fclose(fdparent);
    fclose(fdtypes);
    fclose(fdsubcenter);
    fclose(fdsubmtot);
    fclose(fdsubmgas);
    fclose(fdsubmstars);
    fclose(fdsubsfr);
    fclose(fdsubmcloud);
    fclose(fdcenter);
    fclose(fdmtot);
    fclose(fdmgas);
    fclose(fdmstars);
    fclose(fdsfr);
    fclose(fdmcloud);
    fclose(fd);

    printf("NSubGroupsAll= %"ISYM"\n", NSubGroupsAll);
      
  } // ENDIF ROOT_PROCESSOR

}





















int do_subfind_in_group(FOFData &D, FOF_particle_data *pbuf, int grlen, 
			int *sublen, int *suboffset)
{
  int *Len_bak, *Head_bak, *Tail_bak, *Next_bak;
  FOF_particle_data *P_bak;
  int i, j, oldhead, gr;
  int saved, offset, count, id;



  Len_bak  = D.Len;
  Head_bak = D.Head;
  Tail_bak = D.Len;
  Next_bak = D.Next;
  P_bak	   = D.P;


  D.P	       = pbuf-1;
  D.NumInGroup = grlen;

  D.Energy    = vector(1,  D.NumInGroup);
  D.Density   = vector(1,  D.NumInGroup);
  D.Potential = vector(1,  D.NumInGroup);
  D.Next      = ivector(1, D.NumInGroup);
  D.Head      = ivector(1, D.NumInGroup);
  D.NewHead   = ivector(1, D.NumInGroup);
  D.Tail      = ivector(1, D.NumInGroup);
  D.Len	      = ivector(1, D.NumInGroup);
  D.Index     = ivector(1, D.NumInGroup);
  D.Node      = ivector(1, D.NumInGroup);      
  D.SortId    = ivector(1, D.NumInGroup);      
  
  D.MaxNodes = D.NumInGroup / 10;
  D.GroupTree = new grouptree_data[D.MaxNodes];

  /* initially, there are no subgroups */
  
  for (i = 1; i <= D.NumInGroup; i++) {
    D.Head[i] =	0;
    D.Tail[i] =	0;
    D.Next[i] =	0;
    D.Len[i]  =	0;
    D.Node[i] =	0;
  }
	  
  if (D.NumInGroup >= 2*D.DesLinkNgb && D.NumInGroup > D.DesDensityNgb) {
    ngb_treeallocate(D, D.NumInGroup, 2*D.NumInGroup + 200);
    ngb_treebuild(D, D.NumInGroup);
      
    density(D);        /* compute local densities */
      
    find_subgroups(D);
      
    ngb_treefree();
  } // ENDIF
  else
    D.AnzNodes = 0;
   
  walk_tree_and_unbind(D);
  
  iindexx(D.NumInGroup, D.Head, D.Index);
  
  for (i = 1, D.NSubGroups = 0, oldhead = 0; i <= D.NumInGroup; i++) {
    if (D.Head[D.Index[i]] != oldhead) {
      oldhead = D.Head[D.Index[i]];
      if (oldhead != 1) /* exclude background group */
	D.NSubGroups++;
    } // ENDIF not oldhead
  } // ENDFOR

  if(debug)
    printf("---  %"ISYM" subgroups found. ----\n", D.NSubGroups);
  
  if (D.NSubGroups > 0) {

      /* determine sizes and tags of subgroups */
	  
      D.SubGroupLen = ivector(1, D.NSubGroups); 
      D.SubGroupTag = ivector(1, D.NSubGroups);
      
      for (i = 1, gr = 0, oldhead = 0; i <= D.NumInGroup; i++) {
	if (D.Head[D.Index[i]] != oldhead) {
	  oldhead = D.Head[D.Index[i]];

	  if (oldhead != 1) { /* exclude background group */
	    gr++;
	    D.SubGroupTag[gr] = i; /* first in that group */
	    D.SubGroupLen[gr] = 0;
	  } // ENDIF background group
	} // ENDIF !oldhead

	if (oldhead != 1) /* exclude background group */
	  D.SubGroupLen[gr]++;
      } // ENDFOR

      /* order the subgroups by len */
      
      sort2_int(D.NSubGroups, D.SubGroupLen, D.SubGroupTag);

      /* index will be needed in order_subgroups_potential */      
      for (i = 1; i <= D.NumInGroup; i++)
	D.Head[i] = D.Index[i];
      
      order_subgroups_by_potential(D); 

      for (i = D.NSubGroups, saved = 0, offset = 0, count = 0; i >= 1; i--) {
	for (j = 0; j < D.SubGroupLen[i]; j++) {
	  id = D.Head[D.SubGroupTag[i] + j];
	      
	  D.P[id].MinID = count++;  /* can be set newly at this point */
	} // ENDFOR j
	  
	if (D.SubGroupLen[i] >= D.DesLinkNgb) {
	  sublen[saved]	   = D.SubGroupLen[i];
	  suboffset[saved] = offset;
	      
	  offset += D.SubGroupLen[i];
	  saved++;
	} // ENDIF
      } // ENDFOR subgroups

      /* sort the particles */
      
      qsort(D.P+1, D.NumInGroup, sizeof(FOF_particle_data), comp_func_partminid);
      
      free_ivector(D.SubGroupTag, 1, D.NSubGroups);	  
      free_ivector(D.SubGroupLen, 1, D.NSubGroups); 

      D.NSubGroups = saved;
    }
  
  delete [] D.GroupTree;
  free_ivector(D.SortId, 1, D.NumInGroup);      
  free_ivector(D.Node, 1, D.NumInGroup);      
  free_ivector(D.Index, 1, D.NumInGroup);
  free_ivector(D.Len, 1, D.NumInGroup);
  free_ivector(D.Tail, 1, D.NumInGroup);
  free_ivector(D.NewHead, 1, D.NumInGroup);
  free_ivector(D.Head, 1, D.NumInGroup);
  free_ivector(D.Next, 1, D.NumInGroup);
  free_vector(D.Potential, 1,D.NumInGroup);
  free_vector(D.Density, 1, D.NumInGroup);
  free_vector(D.Energy, 1, D.NumInGroup);
  
  D.Len	 = Len_bak;
  D.Head = Head_bak;
  D.Tail = Tail_bak;
  D.Next = Next_bak;
  D.P	 = P_bak;


  return D.NSubGroups;
}













