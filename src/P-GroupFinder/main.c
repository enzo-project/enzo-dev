#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <mpi.h>

#include <sys/stat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>

#include "allvars.h"
#include "nrsrc/nrutil.h"
#include "proto.h"

double LinkLength = 0.1;   /* in terms of mean interparticle seperation */
int    GroupMinLen= 50;    /*  store only groups in the catalogue 
                               with at least this number of particles */

int    MaxPlacements= 4;
int    Grid= 256;       /* dimension of coarse grid. Note: the actual
                         size of a mesh cell will usually be set to
                         its optimal size, i.e. equal to the linking
                         distance. The coarse grid will then be
                         shifted around to cover the full volume */

void set_units(void)
{
  UnitLength_in_cm =          3.085678e21; 
  UnitMass_in_g    =          1.989e43;
  UnitVelocity_in_cm_per_s =  1.0e5;

  Theta=0.8;  /* opening angle for potential computation */

  DesDensityNgb=  32;
  DesLinkNgb=     DesDensityNgb;    /* DesLinkNgb is also the minimum
				       size of a subgroup */

  if(DesDensityNgb>GroupMinLen)
    {
      if(ThisTask==0)
	printf("One must have DesDensityNgb<=GroupMinLen\n");
      MPI_Finalize(); 
      exit(1);
    }

  UnitTime_in_s= UnitLength_in_cm / UnitVelocity_in_cm_per_s;
  UnitTime_in_Megayears= UnitTime_in_s/SEC_PER_MEGAYEAR;
  G=GRAVITY/pow(UnitLength_in_cm,3)*UnitMass_in_g*pow(UnitTime_in_s,2);
  UnitDensity_in_cgs=UnitMass_in_g/pow(UnitLength_in_cm,3);
  UnitPressure_in_cgs=UnitMass_in_g/UnitLength_in_cm/pow(UnitTime_in_s,2);
  UnitCoolingRate_in_cgs=UnitPressure_in_cgs/UnitTime_in_s;
  UnitEnergy_in_cgs=UnitMass_in_g * pow(UnitLength_in_cm,2) / pow(UnitTime_in_s,2);
  H0 = HUBBLE * UnitTime_in_s;
}


int main(int argc, char **argv)
{
  char input_fname[200];
  char cataloguetxt[200];
  char catalogue_fname[200];
  char subhalos_fname[200];
  char particles_fname[200];
  char parttypes_fname[200];
  char partids_fname[200];
  char subprop_fname[200];
  char fofprop_fname[200];
  char basename[200];
  char buf[200];
  int  Snapshot;
  int  Files;   
  int  FirstInteriorGrid;
  int  formatType;
  char path[500];
  
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);

  if(NTask<=1)
    {
      if(ThisTask==0)
	fprintf(stdout, "Number of processors MUST be a larger than 1.\n");
      MPI_Finalize(); 
      exit(1);
    }

  if(NTask&1)
    {
      if(ThisTask==0)
	fprintf(stdout, "Number of processors MUST be a multiple of 2.\n");
      MPI_Finalize(); 
      exit(1);
    }

  if(argc!=5)
    {
      if(ThisTask==0)
	{
	  fprintf(stderr,"\n\nwrong argument(s).  Specify:\n\n");
	  fprintf(stderr,"<path>      (path)\n");
	  fprintf(stderr,"<basename>  (basename of snapshot files)\n");
	  fprintf(stderr,"<num>       (number of snapshot)\n");
	  fprintf(stderr,"<format>    (0 = GADGET, 1 = Enzo/HDF4)\n");
	  fprintf(stderr,"\n\n");
	}
      MPI_Finalize(); 
      exit(1);
    }

  strcpy(path, argv[1]);
  strcpy(basename, argv[2]);
  Snapshot=atoi(argv[3]);
  formatType = atoi(argv[4]);

  if (formatType != GADGET && formatType != ENZO) {
    fprintf(stderr, "error: please specify a correct format\n"
	    "0 = GADGET, 1 = Enzo/HDF4\n");
    MPI_Finalize();
    exit(1);
  }
  
  if (formatType == GADGET) {
    sprintf(input_fname,     "%s/%s_%03d", path, basename, Snapshot);
    sprintf(particles_fname, "%s/groups/groups_%03d.pos", path, Snapshot);
    sprintf(parttypes_fname, "%s/groups/groups_%03d.types", path, Snapshot);
    sprintf(partids_fname, "%s/groups/groups_%03d.ids", path, Snapshot);
    sprintf(catalogue_fname, "%s/groups/groups_%03d.fofcat", path, Snapshot);
    sprintf(cataloguetxt, "%s/groups/groups_%03d.dat", path, Snapshot);
    sprintf(subhalos_fname,  "%s/groups/groups_%03d.subcat", path, Snapshot);
    sprintf(subprop_fname,   "%s/groups/groups_%03d.subprop", path, Snapshot);
    sprintf(fofprop_fname,   "%s/groups/groups_%03d.fofprop", path, Snapshot);
  } else if (formatType == ENZO) {
    sprintf(input_fname,     "%s/%s%4.4d", path, basename, Snapshot);
    sprintf(particles_fname, "%s/groups/groups_%4.4d.pos", path, Snapshot);
    sprintf(parttypes_fname, "%s/groups/groups_%4.4d.types", path, Snapshot);
    sprintf(partids_fname, "%s/groups/groups_%4.4d.ids", path, Snapshot);
    sprintf(catalogue_fname, "%s/groups/groups_%4.4d.fofcat", path, Snapshot);
    sprintf(cataloguetxt, "%s/groups/groups_%4.4d.dat", path, Snapshot);
    sprintf(subhalos_fname,  "%s/groups/groups_%4.4d.subcat", path, Snapshot);
    sprintf(subprop_fname,   "%s/groups/groups_%4.4d.subprop", path, Snapshot);
    sprintf(fofprop_fname,   "%s/groups/groups_%4.4d.fofprop", path, Snapshot);
  }
  
  sprintf(buf, "%s/groups", path);  mkdir(buf, 0xffff);

  set_units();

  if (formatType == GADGET)
    Files = find_files(input_fname);
  else if (formatType == ENZO) {
    Files = enzoFindFiles(input_fname);
    MarkInteriorParticles(input_fname, Files);
  }

  SearchRadius= LinkLength*BoxSize/pow(NumPart,1.0/3);
  Epsilon=      0.05*BoxSize/pow(NumPart,1.0/3);  /* softening length
						     for potential
						     computation */
  
  if(ThisTask==0)
    printf("\nComoving linking length: %g kpc/h \n\n", SearchRadius);

  if (formatType == GADGET)
    count_local_particles(input_fname, Files);
  else if (formatType == ENZO)
    enzoCountLocalParticles(input_fname, Files);

  if (formatType == GADGET)
    loadpositions(input_fname, Files);
  else if (formatType == ENZO)
    enzoLoadPositions(input_fname, Files);

  /*  adjust_sfr();*/

  marking();
 
  exchange_shadow();

  init_coarse_grid();
  
  link_local_slab(); 
    
  do
    {
      find_minids();
    }
  while(link_accross()>0);
  
  find_minids();
  
  stitch_together(); 

  compile_group_catalogue();

#ifdef FOF_ONLY
  save_groups(particles_fname, catalogue_fname, parttypes_fname, 
	      partids_fname);
#else
  subfind(particles_fname, catalogue_fname, subhalos_fname, 
	  parttypes_fname, partids_fname, subprop_fname, fofprop_fname);
  write_ascii_catalog(catalogue_fname, fofprop_fname, cataloguetxt);
#endif

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize(); /* clean up & finalize MPI */
  exit(0);
}






void save_groups(char *particles_fname, char *catalogue_fname,
		 char *parttypes_fname, char *partids_fname)
{
  FILE   *fd, *fdtypes, *fdids;
  int    i, gr, offset;
  int    ntot;
  int    head, len;
  char   ctype;
  float cm[3], mtot, mgas, mstars, sfr, mcloud;	  
  struct particle_data *Pbuf;

  if(ThisTask==0)
    {
      printf("writing group catalogue...\n");
      
      if(!(fd=fopen(catalogue_fname,"w")))
	{
	  printf("can't open file `%s`\n", catalogue_fname);
	  MPI_Abort(MPI_COMM_WORLD, 1); exit(1);
	}

      fwrite(&NgroupsAll, sizeof(int), 1, fd);
      
      for(gr=NgroupsAll-1; gr>=0; gr--)
	fwrite(&GroupDatAll[gr].Len, sizeof(int), 1, fd);
      
      for(gr=NgroupsAll-1, offset=0; gr>=0; gr--)
	{
	  fwrite(&offset, sizeof(int), 1, fd);
	  offset+= GroupDatAll[gr].Len;
	}
      fclose(fd);
    }

  if(ThisTask==0)
    {
      printf("writing group particles...\n");
  
      if(!(fd=fopen(particles_fname,"w")))
	{
	  printf("can't open file `%s`\n", particles_fname );
	  MPI_Abort(MPI_COMM_WORLD, 1); exit(1);
	}

      if(!(fdtypes=fopen(parttypes_fname,"w")))
	{
	  printf("can't open file `%s`\n", parttypes_fname );
	  MPI_Abort(MPI_COMM_WORLD, 1); exit(1);
	}

      if(!(fdids=fopen(partids_fname,"w")))
	{
	  printf("can't open file `%s`\n", partids_fname );
	  MPI_Abort(MPI_COMM_WORLD, 1); exit(1);
	}

      for(gr=NgroupsAll-1, ntot=0; gr>=0; gr--)
	ntot+= GroupDatAll[gr].Len;

      fwrite(&ntot, sizeof(int), 1, fd);
      fwrite(&ntot, sizeof(int), 1, fdtypes);
      fwrite(&ntot, sizeof(int), 1, fdids);
    }


  for(gr=NgroupsAll-1; gr>=0; gr--)
    {
      if(ThisTask==0)
	{
	  head= GroupDatAll[gr].Tag;
	  len=  GroupDatAll[gr].Len;
	  Pbuf= mymalloc(sizeof(struct particle_data)*len);
	}

      get_particles(0, head, len, Pbuf);

      if(ThisTask==0)
	{
	  /* sort the particles 
	     qsort(Pbuf, len, sizeof(struct particle_data), comp_func_partcoord);
	  */

	  for(i=0; i<len; i++)
	    fwrite(&Pbuf[i].Pos[0], sizeof(double), 3, fd);

	  for(i=0; i<len; i++)
	    fwrite(&Pbuf[i].PartID, sizeof(int), 1, fdids);

	  for(i=0; i<len; i++)
	    {
	      ctype= Pbuf[i].Type;
	      fwrite(&ctype, sizeof(char), 1, fdtypes);
	    }

	  get_properties(Pbuf, len, &cm[0], &mtot, &mgas, &mstars, &sfr, &mcloud);

	  free(Pbuf);
	}
    }

  if(ThisTask==0)
    {
      fclose(fd);
      fclose(fdids);
      fclose(fdtypes);
      printf("done.\n");
    }
}





int get_particles(int dest, int minid, int len, struct particle_data *buf) 
{
  MPI_Status status;
  int i, imin, imax, pp, nlocal, nrecv;
  struct particle_data *localbuf;

  MPI_Bcast(&minid,  1, MPI_INT, dest, MPI_COMM_WORLD);
  MPI_Bcast(&len,    1, MPI_INT, dest, MPI_COMM_WORLD);

  localbuf= mymalloc(sizeof(struct particle_data)*len);

  nlocal=0;

  if(minid>=(1+Noffset[ThisTask]) && minid<(1+Noffset[ThisTask]+Nslab[ThisTask]))
    {
      pp= Head[minid-Noffset[ThisTask]];
      do
	{
	  if(pp<=Nslab[ThisTask])
	    localbuf[nlocal++]= P[pp];
	}
      while(pp=Next[pp]);      
    }
  else
    {
      if(Ncontrib)  
	{
	  imin=0; imax=Ncontrib-1;  /* do a bisection to speed things up */

	  if(ContribID[imin]<= minid && minid<=ContribID[imax])
	    {
	      i=(imin+imax)/2;

	      while(ContribID[i] != minid)
		{
		  i=(imin+imax)/2;
	      
		  if(ContribID[i] < minid)
		    imin= i;

		  if(minid < ContribID[i])
		    imax= i;

		  if((imax-imin)<=1)
		    {
		      if(ContribID[imax] == minid)
			i=imax;
		      if(ContribID[imin] == minid)
			i=imin;
		      break;
		    }
		}

	      if(ContribID[i] == minid)
		{
		  pp= ContribHead[i];
		  do
		    {
		      if(pp<=Nslab[ThisTask])
			localbuf[nlocal++]= P[pp];
		    }
		  while(pp=Next[pp]);
		}
	    }
	}
    }

 
  if(ThisTask== dest)
    {
      for(i=0; i<nlocal; i++)
	buf[i]= localbuf[i];
      
      for(i=0; i<NTask; i++)
	{
	  if(i!=ThisTask)
	    {
	      MPI_Recv(&nrecv, 1, MPI_INT, i, i, MPI_COMM_WORLD, &status);
	      if(nrecv)
		{
		  MPI_Recv(&buf[nlocal], nrecv*sizeof(struct particle_data), 
			   MPI_BYTE, i, i, MPI_COMM_WORLD, &status);
		  nlocal+= nrecv;
		}
 	    }
	}
    }
  else
    {
      MPI_Ssend(&nlocal, 1, MPI_INT, dest, ThisTask, MPI_COMM_WORLD);
      if(nlocal)
	MPI_Ssend(localbuf, nlocal*sizeof(struct particle_data), 
		  MPI_BYTE, dest, ThisTask, MPI_COMM_WORLD);
    }

  free(localbuf);

  if(ThisTask==dest && nlocal!=len)
    {
      printf("local= %d  len=%d\n", nlocal, len);
      fflush(stdout);
      MPI_Abort(MPI_COMM_WORLD, 7777);

    }

  return len;
}





int link_accross(void)
{
  MPI_Status status;
  struct particle_data *buftoleft, *buftoright, *buffer;
  int    i, j, slab, nl, nr, nbuf, len;
  int    leftTask, rightTask;
  int    pp, newid, nlinktot;
  struct id_data *iddat;

  buftoleft=  mymalloc(NtoLeft[ThisTask]  * sizeof(struct particle_data));
  buftoright= mymalloc(NtoRight[ThisTask] * sizeof(struct particle_data));
  buffer=     mymalloc(Nshadow[ThisTask]  * sizeof(struct particle_data));

  nl= nr= nbuf= 0;
  
  for(i=1; i<=Nslab[ThisTask]; i++)
    {
      slab= (P[i].Pos[0]/BoxSize)*NTask;

      if(P[i].Pos[0] < slab*(BoxSize/NTask)+SearchRadius)
	buftoleft[nl++]= P[i];
      
      if(P[i].Pos[0] > (slab+1)*(BoxSize/NTask)-SearchRadius)
	buftoright[nr++]= P[i];
    }

  rightTask= ThisTask+1;
  if(rightTask==NTask)
    rightTask=0;

  leftTask= ThisTask-1;
  if(leftTask<0)
    leftTask=NTask-1;

  if(ThisTask&1)
    {
      MPI_Recv(&buffer[nbuf], NtoLeft[rightTask]*sizeof(struct particle_data), MPI_BYTE, 
	       rightTask, rightTask, MPI_COMM_WORLD, &status);
      nbuf+=  NtoLeft[rightTask];
      MPI_Ssend(buftoright, NtoRight[ThisTask]*sizeof(struct particle_data), MPI_BYTE, 
		rightTask, ThisTask, MPI_COMM_WORLD);
    }
  else
    {
      MPI_Ssend(buftoleft, NtoLeft[ThisTask]*sizeof(struct particle_data), MPI_BYTE, 
		leftTask,ThisTask, MPI_COMM_WORLD);
      MPI_Recv(&buffer[nbuf], NtoRight[leftTask]*sizeof(struct particle_data), MPI_BYTE, 
	       leftTask, leftTask, MPI_COMM_WORLD, &status);
      nbuf+=  NtoRight[leftTask];
    }

  if(ThisTask&1)
    {
      MPI_Recv(&buffer[nbuf], NtoRight[leftTask]*sizeof(struct particle_data), MPI_BYTE, 
	       leftTask,leftTask, MPI_COMM_WORLD, &status);
      nbuf+=  NtoRight[leftTask];
      MPI_Ssend(buftoleft, NtoLeft[ThisTask]*sizeof(struct particle_data), MPI_BYTE, 
		leftTask, ThisTask, MPI_COMM_WORLD);
    }
  else
    {
      MPI_Ssend(buftoright, NtoRight[ThisTask]*sizeof(struct particle_data), MPI_BYTE, 
		rightTask,ThisTask, MPI_COMM_WORLD);
      MPI_Recv(&buffer[nbuf], NtoLeft[rightTask]*sizeof(struct particle_data), MPI_BYTE, 
	       rightTask,rightTask, MPI_COMM_WORLD, &status);
      nbuf+=  NtoLeft[rightTask];
    }


  iddat= mymalloc(nbuf*sizeof(struct id_data));

  for(i=0; i<nbuf; i++)
    {
      iddat[i].ID=    buffer[i].MinID;
      iddat[i].index= 1+Nslab[ThisTask]+i;
    }

  qsort(iddat, nbuf, sizeof(struct id_data), comp_func);

  NLinkAccross= 0;

  for(i=0; i<(nbuf-1);)
    {
      j=i+1; 
      while(iddat[i].ID == iddat[j].ID)
	{
	  linkit(iddat[i].index, iddat[j].index);
	  j++;
	  if(j>=nbuf)
	    break;
	}
      i=j;
    }
  free(iddat);
  free(buffer);
  free(buftoright);
  free(buftoleft);

  MPI_Allreduce(&NLinkAccross, &nlinktot, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(ThisTask==0)
    printf("nlinktot= %d\n", nlinktot);

  return nlinktot;
}




void compile_group_catalogue(void)
{
  MPI_Status  status;
  int i, n, gr, tot, count;
  int nbound, Nbound;
  
  for(n=1, Ngroups=Ncontrib=nbound=0; n<=Nlocal; n++)
    {
      if(Head[n]==n)
	if(P[n].GrLen>=GroupMinLen)
	  {
	    if(P[n].MinID>=(1+Noffset[ThisTask]) && P[n].MinID<(1+Noffset[ThisTask]+Nslab[ThisTask]))
	      {
		Ngroups++;  /* task hosts the global head of the group */
		nbound+= P[n].GrLen;
	      }
	    else
	      {
		Ncontrib++;  /* task hosts a contribution to a group */
	      }
	  }
    }

  MPI_Allreduce(&Ngroups,  &NgroupsAll, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&nbound,   &Nbound, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  if(ThisTask==0)
    {
      printf("Number of groups: %d  bound: %d\n", NgroupsAll, Nbound);
    }

  
  GroupDat=  mymalloc(Ngroups*sizeof(struct gr_data));
  
  ContribID=   ivector(0, Ncontrib-1); 
  ContribHead= ivector(0, Ncontrib-1);

  for(n=1, Ngroups=Ncontrib=0; n<=Nlocal; n++)
    {
      if(Head[n]==n)
	if(P[n].GrLen>=GroupMinLen)
	  {
	    if(P[n].MinID>=(1+Noffset[ThisTask]) && P[n].MinID<(1+Noffset[ThisTask]+Nslab[ThisTask]))
	      {
		GroupDat[Ngroups].Len= P[n].GrLen;
		GroupDat[Ngroups].Tag= P[n].MinID;
		Ngroups++; 
	      }
	    else
	      {
		ContribID[Ncontrib]=   P[n].MinID; 
		ContribHead[Ncontrib]= Head[n];
		Ncontrib++;
	      }
	  }
    }

  if(Ncontrib>0)
    sort2_int(Ncontrib, ContribID-1, ContribHead-1); 

  NgroupsList= mymalloc(sizeof(int)*NTask);

  MPI_Allgather(&Ngroups, 1, MPI_INT, NgroupsList, 1, MPI_INT, MPI_COMM_WORLD);


  GroupDatAll= mymalloc(NgroupsAll*sizeof(struct gr_data)); 


  if(ThisTask==0)
    {
      for(i=0; i<Ngroups; i++)
	GroupDatAll[i]= GroupDat[i];
      
      count= Ngroups;

      for(i=1; i<NTask; i++)
	{
	  MPI_Recv(&GroupDatAll[count], NgroupsList[i]*sizeof(struct gr_data), MPI_BYTE, i, 0, MPI_COMM_WORLD, &status);
	  count+= NgroupsList[i];
	}
    }
  else
    {
      MPI_Ssend(&GroupDat[0], Ngroups*sizeof(struct gr_data), MPI_BYTE, 0, 0, MPI_COMM_WORLD);
    }
  
  if(ThisTask==0)
    {
      qsort(GroupDatAll, NgroupsAll, sizeof(struct gr_data), comp_func_gr);
      /* ok sorted by size, with secondary criterion minid */  

      printf("\n10 largest groups:\n");
      printf("--------------------------\n");
      for(i=NgroupsAll-1; i>=0; i--)
	{
	  printf("GroupNr= %2d   Size= %d\n", NgroupsAll-i, GroupDatAll[i].Len);
	  if(i<=NgroupsAll-10)
	    break;
	}
      printf("--------------------------\n");
    }      

  MPI_Bcast(GroupDatAll, NgroupsAll*sizeof(struct gr_data), MPI_BYTE, 0, MPI_COMM_WORLD); 
}




void find_minids(void)
{
  int n, pp, len, sum=0;

  for(n=1; n<=Nlocal; n++)
    {
      if(Head[n]==n)
	{
	  pp= n; len=0;
	  do
	    {
	      if(P[pp].ID>=(1+Noffset[ThisTask]) && P[pp].ID<(1+Noffset[ThisTask]+Nslab[ThisTask]))
		len++;
	    }
	  while(pp=Next[pp]);


	  P[n].MinID= P[n].ID;
	  P[n].GrLen= len; /* Len[n]; */

	  pp= n;
	  while(pp=Next[pp])
	    {
	      if(P[n].MinID> P[pp].ID)
		P[n].MinID= P[pp].ID;
	    }

	  pp= n;
	  while(pp=Next[pp])
	    {
	      P[pp].MinID= P[n].MinID;
	      P[pp].GrLen= P[n].GrLen;
	    }

	  sum+= len;
	}
    }
}





void stitch_together(void)
{
  MPI_Status status;
  struct particle_data *buftoleft, *buftoright, *buffer;
  int    i, slab, nl, nr, nbuf, len;
  int    leftTask, rightTask;
  int    pp, newid;
  struct idmin_data *iddat;


  buftoleft=  mymalloc(NtoLeft[ThisTask]  * sizeof(struct particle_data));
  buftoright= mymalloc(NtoRight[ThisTask] * sizeof(struct particle_data));
  buffer=     mymalloc(Nshadow[ThisTask]  * sizeof(struct particle_data));

  nl= nr= nbuf= 0;
  
  for(i=1; i<=Nslab[ThisTask]; i++)
    {
      slab= (P[i].Pos[0]/BoxSize)*NTask;
		  
      if(P[i].Pos[0] < slab*(BoxSize/NTask)+SearchRadius)
	buftoleft[nl++]= P[i];
      
      if(P[i].Pos[0] > (slab+1)*(BoxSize/NTask)-SearchRadius)
	buftoright[nr++]= P[i];
    }

  rightTask= ThisTask+1;
  if(rightTask==NTask)
    rightTask=0;

  leftTask= ThisTask-1;
  if(leftTask<0)
    leftTask=NTask-1;

  if(ThisTask&1)
    {
      MPI_Recv(&buffer[nbuf], NtoLeft[rightTask]*sizeof(struct particle_data), MPI_BYTE, 
	       rightTask, rightTask, MPI_COMM_WORLD, &status);
      nbuf+=  NtoLeft[rightTask];
      MPI_Ssend(buftoright, NtoRight[ThisTask]*sizeof(struct particle_data), MPI_BYTE, 
		rightTask, ThisTask, MPI_COMM_WORLD);
    }
  else
    {
      MPI_Ssend(buftoleft, NtoLeft[ThisTask]*sizeof(struct particle_data), MPI_BYTE, 
		leftTask,ThisTask, MPI_COMM_WORLD);
      MPI_Recv(&buffer[nbuf], NtoRight[leftTask]*sizeof(struct particle_data), MPI_BYTE, 
	       leftTask, leftTask, MPI_COMM_WORLD, &status);
      nbuf+=  NtoRight[leftTask];
    }

  if(ThisTask&1)
    {
      MPI_Recv(&buffer[nbuf], NtoRight[leftTask]*sizeof(struct particle_data), MPI_BYTE, 
	       leftTask,leftTask, MPI_COMM_WORLD, &status);
      nbuf+=  NtoRight[leftTask];
      MPI_Ssend(buftoleft, NtoLeft[ThisTask]*sizeof(struct particle_data), MPI_BYTE, 
		leftTask, ThisTask, MPI_COMM_WORLD);
    }
  else
    {
      MPI_Ssend(buftoright, NtoRight[ThisTask]*sizeof(struct particle_data), MPI_BYTE, 
		rightTask,ThisTask, MPI_COMM_WORLD);
      MPI_Recv(&buffer[nbuf], NtoLeft[rightTask]*sizeof(struct particle_data), MPI_BYTE, 
	       rightTask,rightTask, MPI_COMM_WORLD, &status);
      nbuf+=  NtoLeft[rightTask];
    }

  
  iddat= mymalloc(nbuf*sizeof(struct idmin_data));

  for(i=0; i<nbuf; i++)
    {
      iddat[i].minID=    buffer[i].MinID;
      iddat[i].index=    1+Nslab[ThisTask]+i;
      iddat[i].len=      buffer[i].GrLen;
    }

  qsort(iddat, nbuf, sizeof(struct idmin_data), comp_func2);

  for(i=0; i<nbuf;i++)
    {
      if(iddat[i].minID!= P[iddat[i].index].MinID || iddat[i].len>0)
	{
	  newid= P[iddat[i].index].MinID;
	  if(iddat[i].minID<newid)
	    newid= iddat[i].minID;

	  len=  P[iddat[i].index].GrLen + iddat[i].len;

	  pp= Head[iddat[i].index];
	  do
	    {
	      P[pp].MinID= newid;
	      P[pp].GrLen=   len;
	    }
	  while(pp=Next[pp]);

	  while(i<(nbuf-1))
	    {
	      if(iddat[i+1].minID == iddat[i].minID)
		i++;
	      else
		break;
	    }
	}
    }


  free(iddat);
  free(buffer);
  free(buftoright);
  free(buftoleft);
}




void exchange_shadow(void)
{
  MPI_Status status;
  struct particle_data *buftoleft, *buftoright;
  int    i, slab, nl, nr;
  int    leftTask, rightTask;

  if (ThisTask == 0) {
    fprintf(stdout, "exchanging shadows ...\n");
    fflush(stdout);
  }

  buftoleft=  mymalloc(NtoLeft[ThisTask]  * sizeof(struct particle_data));
  buftoright= mymalloc(NtoRight[ThisTask] * sizeof(struct particle_data));

  nl= nr= 0;
  
  for(i=1; i<=Nlocal; i++)
    {
      slab= (P[i].Pos[0]/BoxSize)*NTask;

      if(slab!=ThisTask)
	MPI_Abort(MPI_COMM_WORLD, 11);
		  
      if(P[i].Pos[0] < slab*(BoxSize/NTask)+SearchRadius)
	buftoleft[nl++]= P[i];
      
      if(P[i].Pos[0] > (slab+1)*(BoxSize/NTask)-SearchRadius)
	buftoright[nr++]= P[i];
    }

  if (nl!=NtoLeft[ThisTask]) {
    fprintf(stderr, "[proc %d] error: shadows don't match! "
	    "nl = %d, NtoLeft[%d] = %d\n", 
	    ThisTask, nl, ThisTask, NtoLeft[ThisTask]);
    MPI_Abort(MPI_COMM_WORLD, 11);
  }
  
  if (nr!=NtoRight[ThisTask]) {
    fprintf(stderr, "[proc %d] error: shadows don't match! "
	    "nr = %d, NtoRight[%d] = %d\n", 
	    ThisTask, nr, ThisTask, NtoRight[ThisTask]);
    MPI_Abort(MPI_COMM_WORLD, 11);
  }

  rightTask= ThisTask+1;
  if(rightTask==NTask)
    rightTask=0;

  leftTask= ThisTask-1;
  if(leftTask<0)
    leftTask=NTask-1;


  if(ThisTask&1)
    {
      MPI_Recv(&P[1+Nlocal], NtoLeft[rightTask]*sizeof(struct particle_data), MPI_BYTE, 
	       rightTask, rightTask, MPI_COMM_WORLD, &status);
      Nlocal+=  NtoLeft[rightTask];
      MPI_Ssend(buftoright, NtoRight[ThisTask]*sizeof(struct particle_data), MPI_BYTE, 
		rightTask, ThisTask, MPI_COMM_WORLD);
    }
  else
    {
      MPI_Ssend(buftoleft, NtoLeft[ThisTask]*sizeof(struct particle_data), MPI_BYTE, 
		leftTask,ThisTask, MPI_COMM_WORLD);
      MPI_Recv(&P[1+Nlocal], NtoRight[leftTask]*sizeof(struct particle_data), MPI_BYTE, 
	       leftTask, leftTask, MPI_COMM_WORLD, &status);
      Nlocal+=  NtoRight[leftTask];
    }

  if(ThisTask&1)
    {
      MPI_Recv(&P[1+Nlocal], NtoRight[leftTask]*sizeof(struct particle_data), MPI_BYTE, 
	       leftTask,leftTask, MPI_COMM_WORLD, &status);
      Nlocal+=  NtoRight[leftTask];
      MPI_Ssend(buftoleft, NtoLeft[ThisTask]*sizeof(struct particle_data), MPI_BYTE, 
		leftTask, ThisTask, MPI_COMM_WORLD);
    }
  else
    {
      MPI_Ssend(buftoright, NtoRight[ThisTask]*sizeof(struct particle_data), MPI_BYTE, 
		rightTask,ThisTask, MPI_COMM_WORLD);
      MPI_Recv(&P[1+Nlocal], NtoLeft[rightTask]*sizeof(struct particle_data), MPI_BYTE, 
	       rightTask,rightTask, MPI_COMM_WORLD, &status);
      Nlocal+=  NtoLeft[rightTask];
    }

  free(buftoright);
  free(buftoleft);
}





void link_local_slab(void)
{
  int  nx,ny,nz;
  int  iter;
  int  count;

  iter=1;

  for(GridCorner[0]=0, nx=0; nx<Nx; GridCorner[0]+= (Grid-2.0)/Grid*GridExtension, nx++)
    for(GridCorner[1]=0, ny=0; ny<Ny; GridCorner[1]+= (Grid-2.0)/Grid*GridExtension, ny++)
      for(GridCorner[2]=0, nz=0; nz<Nz; GridCorner[2]+= (Grid-2.0)/Grid*GridExtension, nz++)
	{
	  if(ThisTask==0)
	    printf("Grid placement number: %d out of %d\n", iter++, Nx*Ny*Nz);
	  
	  count= course_binning();
	  
	  if(count)
	    find_groups();
	}
  
  if(ThisTask==0)
    {
      printf("main linking done.\n"); 
      fflush(stdout);
    }
}



void init_coarse_grid(void)
{
  int i;

  GridExtension=Grid * SearchRadius;
  
  if(GridExtension> BoxSize)
    {
      Grid= BoxSize/SearchRadius;
      GridExtension=Grid * SearchRadius;
    }

  if(BoxSize/( (Grid-2.0)/Grid*GridExtension ) > MaxPlacements)
    {
      GridExtension= BoxSize/MaxPlacements *Grid/(Grid-2.001);
    }
 
  Nx= (int)(BoxSize/( (Grid-2.0)/Grid*GridExtension) + 1);
  Ny= (int)(BoxSize/( (Grid-2.0)/Grid*GridExtension) + 1);
  Nz= (int)(BoxSize/( (Grid-2.0)/Grid*GridExtension) + 1);

  if(ThisTask==0)
    printf("\n\nGrid has to be placed (%d|%d|%d) times in each dimension.\n\n", Nx, Ny, Nz);

  GridFirst=i3tensor(0, Grid-1, 0, Grid-1, 0, Grid-1);
  GridLast =i3tensor(0, Grid-1, 0, Grid-1, 0, Grid-1);
  GridFlag =i3tensor(0, Grid-1, 0, Grid-1, 0, Grid-1);

  GridNext= ivector(1, Nlocal);

  Tail=     ivector(1, Nlocal);
  Len=      ivector(1, Nlocal);
  Head=     ivector(1, Nlocal);
  Next=     ivector(1, Nlocal);

  if(ThisTask==0)
    printf("Nlocal= %d Task=%d\n", Nlocal, ThisTask);

  for(i=1; i<=Nlocal; i++)
    {
      Head[i]= Tail[i]= i;
      Next[i]= 0;
      Len[i]=  1;
    }
}





/*
void marking(void)
{
  int i;

  qsort(&P[1], Nlocal, sizeof(struct particle_data), comp_func_partcoord);
  
  for(i=1; i<=Nlocal; i++)
    P[i].ID   =  Noffset[ThisTask] + i; 
}
*/

void marking(void)
/* New code added to ensure no two particles lie on top of each other. */
{
  float posold[3];
  int   i,k,iter,idone;

  if (ThisTask == 0) {
    fprintf(stdout, "marking ...\n");
    fflush(stdout);
  }

  iter = 0;
  do 
    {
      qsort(&P[1], Nlocal, sizeof(struct particle_data), comp_func_partcoord);

      for(i=2, idone=0; i<=Nlocal; i++)
	if(fabs(P[i-1].Pos[0]-P[i].Pos[0])<1.0e-3*Epsilon &&
	   fabs(P[i-1].Pos[1]-P[i].Pos[1])<1.0e-3*Epsilon &&
	   fabs(P[i-1].Pos[2]-P[i].Pos[2])<1.0e-3*Epsilon) 
	  {
	    for(k=0; k<3; k++)
	      {
		posold[k]= P[i].Pos[k];
		P[i].Pos[k] += (0.001*Epsilon)*(2*drand48()-1);
	      }
	    /*
	      printf("Task %2d iter %d: Shifting particle %d (%e|%e|%e) -> (%e|%e|%e)\n",
	      ThisTask,iter,i,
	      posold[0], posold[1], posold[2],
	      P[i].Pos[0],P[i].Pos[1],P[i].Pos[2]);
	    */
	    
	    idone++;
	  }
      fflush(stdout);
      iter++;
    } 
  while(idone>0 && iter<10);

  qsort(&P[1], Nlocal, sizeof(struct particle_data), comp_func_partcoord);

  for(i=1; i<=Nlocal; i++)
    P[i].ID = Noffset[ThisTask] + i; 
}




double periodic(double x)
{
  if(x>0.5*BoxSize)
    x-= BoxSize;

  if(x<-0.5*BoxSize)
    x+= BoxSize;
  
  return x;
}

double  periodic_wrap(double x)
{
  while(x>BoxSize)
    x-= BoxSize;

  while(x<0)
    x+= BoxSize;
  
  return x;
}










int course_binning(void)
{
  int i,j,k,n,count;
  double fac;
  double pos[3];

  if(ThisTask==0)
    {
      printf("course binning..."); fflush(stdout);
    }

  for(i=0;i<Grid;i++)
    for(j=0;j<Grid;j++)
      for(k=0;k<Grid;k++)
	{
	  GridFirst[i][j][k]=0;
	  GridFlag[i][j][k]=0;
	}
  
  for(n=1; n<=Nlocal; n++)
    GridNext[n]=0;
 
  fac=Grid/GridExtension;

  for(n=1, count=0; n<=Nlocal; n++)
    {
      for(k=0; k<3; k++)
	{
	  pos[k]= P[n].Pos[k];
	  if(pos[k]<GridCorner[k])
	    pos[k]+= BoxSize;
	  else
	    {
	      if(pos[k]>=(GridCorner[k]+GridExtension))
		pos[k]-= BoxSize;
	    }
	}

      if(pos[0]>=GridCorner[0])
	if(pos[0]<(GridCorner[0]+GridExtension))
	  if(pos[1]>=GridCorner[1])
	    if(pos[1]<(GridCorner[1]+GridExtension))
	      if(pos[2]>=GridCorner[2])
		if(pos[2]<(GridCorner[2]+GridExtension))
		  {
		    i= (pos[0]-GridCorner[0])*fac;
		    j= (pos[1]-GridCorner[1])*fac;
		    k= (pos[2]-GridCorner[2])*fac;

		    if(GridFirst[i][j][k])
		      {
			GridNext[GridLast[i][j][k]]=n;
			GridLast[i][j][k]=n;
		      }
		    else
		      {
			GridFirst[i][j][k]=GridLast[i][j][k]=n;
		      }
		    count++;
		  }
    }

  if(ThisTask==0)
    {
      printf("done.  (count=%d)\n",count);
      fflush(stdout);
    }

  return count;
}







void find_groups(void)
{
  int i,j,k;
  int p;

  if(ThisTask==0)
    {
      printf("linking..."); fflush(stdout);
    }

  for(i=Grid-2; i>=1; i--)
    for(j=Grid-2; j>=1; j--)
      for(k=Grid-2; k>=1; k--)
	{
	  if(p=GridFirst[i][j][k])
	    {
	      do
		{
		  check_cell(p, i+1,j  ,k  );
		  check_cell(p, i+1,j+1,k  );
		  check_cell(p, i+1,j  ,k+1);
		  check_cell(p, i+1,j+1,k+1);
		  check_cell(p, i  ,j+1,k  );
		  check_cell(p, i  ,j+1,k+1);
		  check_cell(p, i  ,j  ,k+1);

		  check_cell(p, i, j ,k);    
		}
	      while(p=GridNext[p]);
	    }
	} 


  for(i=Grid-2; i>=1; i--)
    for(j=Grid-2; j>=1; j--)
      for(k=Grid-2; k>=1; k--)
	{
	  if(p=GridFirst[i][j][k])
	    {
	      do
		{
		  check_cell(p, i+1,j  ,k-1);
		  check_cell(p, i+1,j-1,k  );
		  check_cell(p, i  ,j-1,k+1);
		  check_cell(p, i-1,j+1,k+1);
		  check_cell(p, i-1,j-1,k+1);
		  check_cell(p, i+1,j-1,k+1);
		}
	      while(p=GridNext[p]);
	    }
	}

  if(ThisTask==0)
    {
      printf("done.\n");
      fflush(stdout);
    }
}




void check_cell(int p, int i, int j, int k)
{
  double r2,s2,dx,dy,dz;
  int pp,ss;
  int s, flag;

  
  if(s=GridFirst[i][j][k])
    {
      flag= Head[s];
      if(GridFlag[i][j][k])
	{
	  if(Head[s]==Head[p])
	    return; 
	}
    }
  else
    flag= 0;

  s2= SearchRadius*SearchRadius;

  while(s)
    {
      if(Head[s]!=flag)
	flag=0;

      if(Head[p]!=Head[s])  /* only if not yet linked */
	{
	  dx= P[p].Pos[0] - P[s].Pos[0];
	  dy= P[p].Pos[1] - P[s].Pos[1];
	  dz= P[p].Pos[2] - P[s].Pos[2];

	  dx= periodic(dx);
	  dy= periodic(dy);
	  dz= periodic(dz);
	  
	  r2=dx*dx + dy*dy + dz*dz;
	  
	  if(r2 < s2)  /* ok, we have a partner */
	    {
	      if(Len[Head[p]] > Len[Head[s]]) /* p group is longer */
		{
		  Next[ Tail[Head[p]] ] = Head[s];

		  Tail[Head[p]] = Tail[Head[s]];

		  Len[Head[p]]+= Len[Head[s]];

		  ss=Head[s];
		  do
		    {
		      Head[ss]=Head[p];
		    }
		  while(ss=Next[ss]);

		  flag=0;
		}
	      else
		{
		  Next[ Tail[Head[s]] ] = Head[p];

		  Tail[Head[s]] = Tail[Head[p]];

		  Len[Head[s]]+= Len[Head[p]];

		  pp=Head[p];
		  do
		    {
		      Head[pp]=Head[s];
		    }
		  while(pp=Next[pp]);

		  flag=0;
		}
	      
	      if(GridFlag[i][j][k])
		return; 
	    }
	}

      s=GridNext[s];
    }

  if(flag)
    GridFlag[i][j][k]= 1;
}




void linkit(int p, int s)
{
  int ss, pp;

  if(Head[p]!=Head[s])  /* only if not yet linked */
    {
      NLinkAccross++;

      if(Len[Head[p]] > Len[Head[s]]) /* p group is longer */
	{
	  Next[ Tail[Head[p]] ] = Head[s];
	  
	  Tail[Head[p]] = Tail[Head[s]];
	  
	  Len[Head[p]]+= Len[Head[s]];
	  
	  ss=Head[s];
	  do
	    {
	      Head[ss]=Head[p];
	    }
	  while(ss=Next[ss]);
	}
      else
	{
	  Next[ Tail[Head[s]] ] = Head[p];
	  
	  Tail[Head[s]] = Tail[Head[p]];
	  
	  Len[Head[s]]+= Len[Head[p]];
	  
	  pp=Head[p];
	  do
	    {
	      Head[pp]=Head[s];
	    }
	  while(pp=Next[pp]);
	}
    }
}


void *mymalloc(size_t size)
{
  void *p;

  p= malloc(size);
  
  if(!p)
    {
      printf("Failed to alloc %d bytes on process %d.\n", (int)size, ThisTask);
      fflush(stdout);
      MPI_Abort(MPI_COMM_WORLD, 88);
    }
  return p;
}
  

void adjust_sfr()
{
  int i;
  double tdyn, ratetdyn;

  for(i=1; i<=Nlocal; i++)
    {
      if(P[i].Type==0)
	{
	  tdyn= 1/sqrt(4*PI*G*P[i].Rho/(Time*Time*Time));

	  ratetdyn = (1-0.1) * 0.01 * P[i].Mclouds/tdyn;

	  ratetdyn *= (UnitMass_in_g/SOLAR_MASS)/(UnitTime_in_s/SEC_PER_YEAR);
	  
	  if(P[i].Sfr>1.1*ratetdyn && P[i].Rho<1.0e-5)
	    P[i].Sfr = ratetdyn;
	}
    }
  
}

void write_ascii_catalog(char *catalogue_fname, char *fofprop_fname, 
			 char *cataloguetxt)
{

  FILE *in, *out;
  int i, Ngroups, *Npart_g;
  float *CM, *Mass;
  double SolarMass = UnitMass_in_g / 1.989e33;
  float RootBoxSize[3];

  if (ThisTask == 0) {

    for (i = 0; i < 3; i++)
      RootBoxSize[i] = BoxSize / (rightEdge[i] - leftEdge[i]);

    in = fopen(fofprop_fname, "r");
    assert(in != NULL);

    fread(&Ngroups, sizeof(int), 1, in);
    CM = malloc(sizeof(float) * 3 * Ngroups);
    Mass = malloc(sizeof(float) * Ngroups);
    fread(CM, sizeof(float), 3*Ngroups, in);
    fread(Mass, sizeof(float), Ngroups, in);
    fclose(in);

    in = fopen(catalogue_fname, "r");
    assert(in != NULL);
    fread(&Ngroups, sizeof(int), 1, in);
    Npart_g = malloc(sizeof(int) * Ngroups);
    fread(Npart_g, sizeof(int), Ngroups, in);
    fclose(in);

    out = fopen(cataloguetxt, "w");
    assert(out != NULL);
    fprintf(out, "#%7s %12s %12s %12s %12s %12s\n",
	    "Halo", "Npart", "Mass [Msun]", "CM(x)", "CM(y)", "CM(z)");
    for (i = 0; i < Ngroups; i++)
      fprintf(out, "%8d %12d %12.6g %12.6f %12.6f %12.6f\n",
	      i, Npart_g[i], Mass[i]*SolarMass,
	      CM[3*i+0] / RootBoxSize[0] + leftEdge[0],
	      CM[3*i+1] / RootBoxSize[1] + leftEdge[1],
	      CM[3*i+2] / RootBoxSize[2] + leftEdge[2]);
    fclose(out);

  } // ENDIF root task

}
