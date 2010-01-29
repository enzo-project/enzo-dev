#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"


struct io_header_1
{
  int npart[6];
  double   mass[6];
  double   time;
  double   redshift;
  int flag_sfr;
  int flag_feedback;
  int npartTotal[6];
  int flag_cooling;
  int num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  int      flag_multiphase;
  int      flag_stellarage;
  int      flag_sfrhistogram;
  char     fill[84];  /* fills to 256 Bytes */
} header1;


#define BUFSIZE 10000
#define SKIP fread(&dummy, sizeof(dummy), 1, fd);

void loadpositions(char *fname, int files)
{
  FILE        *fd,*fdvel,*fdPID,*fdmass,*fdu,*fdrho,*fdmfs,*fdsfr,*fdclouds;
  MPI_Status  status;  
  char        buf[200];
  int         i,j,k,dummy;
  int         n,slab,particleid,counter;
  float       pos[3],vel[3],mass,intenergy,rho;
  float       mfs,sfr,mclouds;
  struct particle_data *Pbuf;
  int         *Nbuf;

  allocate_memory();
  Nlocal= 0;
  
  Pbuf= malloc(sizeof(struct particle_data)*BUFSIZE*NTask);
  Nbuf= malloc(sizeof(int)*NTask);
     
  if(ThisTask==0)
    {
      for(i=0; i<NTask; i++)
	Nbuf[i]=0;
      
      for(i=0; i<files; i++)
	{
	  if(files>1)
	    sprintf(buf,"%s.%d",fname,i);
	  else
	    sprintf(buf,"%s",fname);
	  
	  if(!(fd=fopen(buf,"r")))
	    {
	      printf("can't open file `%s`\n",buf);
	      MPI_Abort(MPI_COMM_WORLD, 1); exit(1);
	    }
	  
	  printf("reading `%s' ...\n",buf); fflush(stdout);


	  fdvel=fopen(buf,"r");
	  fread(&dummy, sizeof(dummy), 1, fdvel);
	  fread(&header1, sizeof(header1), 1, fdvel);
	  fread(&dummy, sizeof(dummy), 1, fdvel);
	  fread(&dummy, sizeof(dummy), 1, fdvel);
	  for(k=0; k<6; k++) /* skip pos data */
	    fseek(fdvel, header1.npart[k]*sizeof(float)*3, SEEK_CUR);
	  fread(&dummy, sizeof(dummy), 1, fdvel);
	  fread(&dummy, sizeof(dummy), 1, fdvel);

	  fdPID=fopen(buf,"r");
	  fread(&dummy, sizeof(dummy), 1, fdPID);
	  fread(&header1, sizeof(header1), 1, fdPID);
	  fread(&dummy, sizeof(dummy), 1, fdPID);
	  fread(&dummy, sizeof(dummy), 1, fdPID);
	  for(k=0; k<6; k++) /* skip pos data */
	    fseek(fdPID, header1.npart[k]*sizeof(float)*3, SEEK_CUR);
	  fread(&dummy, sizeof(dummy), 1, fdPID);
	  fread(&dummy, sizeof(dummy), 1, fdPID);
	  for(k=0; k<6; k++) /* skip vel data */
	    fseek(fdPID, header1.npart[k]*sizeof(float)*3, SEEK_CUR);
	  fread(&dummy, sizeof(dummy), 1, fdPID);
	  fread(&dummy, sizeof(dummy), 1, fdPID);

	  fdmass=fopen(buf,"r");
	  fread(&dummy, sizeof(dummy), 1, fdmass);
	  fread(&header1, sizeof(header1), 1, fdmass);
	  fread(&dummy, sizeof(dummy), 1, fdmass);
	  fread(&dummy, sizeof(dummy), 1, fdmass);
	  for(k=0; k<6; k++) /* skip pos data */
	    fseek(fdmass, header1.npart[k]*sizeof(float)*3, SEEK_CUR);
	  fread(&dummy, sizeof(dummy), 1, fdmass);
	  fread(&dummy, sizeof(dummy), 1, fdmass);
	  for(k=0; k<6; k++) /* skip vel data */
	    fseek(fdmass, header1.npart[k]*sizeof(float)*3, SEEK_CUR);
	  fread(&dummy, sizeof(dummy), 1, fdmass);
	  fread(&dummy, sizeof(dummy), 1, fdmass);
	  for(k=0; k<6; k++) /* skip id data */
	    fseek(fdmass, header1.npart[k]*sizeof(int), SEEK_CUR);
	  fread(&dummy, sizeof(dummy), 1, fdmass);
	  fread(&dummy, sizeof(dummy), 1, fdmass);

	  fdu=fopen(buf,"r");
	  fread(&dummy, sizeof(dummy), 1, fdu);
	  fread(&header1, sizeof(header1), 1, fdu);
	  fread(&dummy, sizeof(dummy), 1, fdu);
	  fread(&dummy, sizeof(dummy), 1, fdu);
	  for (k=0; k<6; k++) /* skip pos data */
	    fseek(fdu, header1.npart[k]*sizeof(float)*3, SEEK_CUR);
	  fread(&dummy, sizeof(dummy), 1, fdu);
	  fread(&dummy, sizeof(dummy), 1, fdu);
	  for (k=0; k<6; k++) /* skip vel data */
	    fseek(fdu, header1.npart[k]*sizeof(float)*3, SEEK_CUR);
	  fread(&dummy, sizeof(dummy), 1, fdu);
	  fread(&dummy, sizeof(dummy), 1, fdu);
	  for (k=0; k<6; k++) /* skip id data */
	    fseek(fdu, header1.npart[k]*sizeof(int), SEEK_CUR);
	  fread(&dummy, sizeof(dummy), 1, fdu);
	  fread(&dummy, sizeof(dummy), 1, fdu);
	  for (k=0; k<6; k++) /* skip masses */
	    if (header1.mass[k]==0)
	      fseek(fdu, header1.npart[k]*sizeof(float), SEEK_CUR); 
	  fread(&dummy, sizeof(dummy), 1, fdu);
	  fread(&dummy, sizeof(dummy), 1, fdu);

	  fdrho=fopen(buf,"r");
	  fread(&dummy, sizeof(dummy), 1, fdrho);
	  fread(&header1, sizeof(header1), 1, fdrho);
	  fread(&dummy, sizeof(dummy), 1, fdrho);
	  fread(&dummy, sizeof(dummy), 1, fdrho);
	  for (k=0; k<6; k++) /* skip pos data */
	    fseek(fdrho, header1.npart[k]*sizeof(float)*3, SEEK_CUR);
	  fread(&dummy, sizeof(dummy), 1, fdrho);
	  fread(&dummy, sizeof(dummy), 1, fdrho);
	  for (k=0; k<6; k++) /* skip vel data */
	    fseek(fdrho, header1.npart[k]*sizeof(float)*3, SEEK_CUR);
	  fread(&dummy, sizeof(dummy), 1, fdrho);
	  fread(&dummy, sizeof(dummy), 1, fdrho);
	  for (k=0; k<6; k++) /* skip id data */
	    fseek(fdrho, header1.npart[k]*sizeof(int), SEEK_CUR);
	  fread(&dummy, sizeof(dummy), 1, fdrho);
	  fread(&dummy, sizeof(dummy), 1, fdrho);
	  for (k=0; k<6; k++) /* skip masses */
	    if (header1.mass[k]==0)
	      fseek(fdrho, header1.npart[k]*sizeof(float), SEEK_CUR); 
	  fread(&dummy, sizeof(dummy), 1, fdrho);
	  fread(&dummy, sizeof(dummy), 1, fdrho);  /* skip u */
	  fseek(fdrho, header1.npart[0]*sizeof(float), SEEK_CUR);
	  fread(&dummy, sizeof(dummy), 1, fdrho);
	  fread(&dummy, sizeof(dummy), 1, fdrho);

	  if(header1.flag_sfr)
	    {
	      fdmfs=fopen(buf,"r");
	      fread(&dummy, sizeof(dummy), 1, fdmfs);
	      fread(&header1, sizeof(header1), 1, fdmfs);
	      fread(&dummy, sizeof(dummy), 1, fdmfs);
	      fread(&dummy, sizeof(dummy), 1, fdmfs);
	      for(k=0; k<6; k++) /* skip pos data */
		fseek(fdmfs, header1.npart[k]*sizeof(float)*3, SEEK_CUR);
	      fread(&dummy, sizeof(dummy), 1, fdmfs);
	      fread(&dummy, sizeof(dummy), 1, fdmfs);
	      for(k=0; k<6; k++) /* skip vel data */
		fseek(fdmfs, header1.npart[k]*sizeof(float)*3, SEEK_CUR);
	      fread(&dummy, sizeof(dummy), 1, fdmfs);
	      fread(&dummy, sizeof(dummy), 1, fdmfs);
	      for(k=0; k<6; k++) /* skip id data */
		fseek(fdmfs, header1.npart[k]*sizeof(int), SEEK_CUR);
	      fread(&dummy, sizeof(dummy), 1, fdmfs);

	      fread(&dummy, sizeof(dummy), 1, fdmfs);
	      for(k=0; k<6; k++) /* skip masses */
		if(header1.mass[k]==0)
		  fseek(fdmfs, header1.npart[k]*sizeof(float), SEEK_CUR);
	      fread(&dummy, sizeof(dummy), 1, fdmfs);

	      fread(&dummy, sizeof(dummy), 1, fdmfs);  /* skip u */
	      fseek(fdmfs, header1.npart[0]*sizeof(float), SEEK_CUR);
	      fread(&dummy, sizeof(dummy), 1, fdmfs);

	      fread(&dummy, sizeof(dummy), 1, fdmfs);  /* skip rho */
	      fseek(fdmfs, header1.npart[0]*sizeof(float), SEEK_CUR);
	      fread(&dummy, sizeof(dummy), 1, fdmfs);

	      fread(&dummy, sizeof(dummy), 1, fdmfs);  /* skip ne */
	      fseek(fdmfs, header1.npart[0]*sizeof(float), SEEK_CUR);
	      fread(&dummy, sizeof(dummy), 1, fdmfs);

	      fread(&dummy, sizeof(dummy), 1, fdmfs);  /* skip nh */
	      fseek(fdmfs, header1.npart[0]*sizeof(float), SEEK_CUR);
	      fread(&dummy, sizeof(dummy), 1, fdmfs);
      
	      fread(&dummy, sizeof(dummy), 1, fdmfs);
	    }

	  if(header1.flag_sfr)
	    {
	      fdsfr=fopen(buf,"r");
	      fread(&dummy, sizeof(dummy), 1, fdsfr);
	      fread(&header1, sizeof(header1), 1, fdsfr);
	      fread(&dummy, sizeof(dummy), 1, fdsfr);
	      fread(&dummy, sizeof(dummy), 1, fdsfr);
	      for(k=0; k<6; k++) /* skip pos data */
		fseek(fdsfr, header1.npart[k]*sizeof(float)*3, SEEK_CUR);
	      fread(&dummy, sizeof(dummy), 1, fdsfr);
	      fread(&dummy, sizeof(dummy), 1, fdsfr);
	      for(k=0; k<6; k++) /* skip vel data */
		fseek(fdsfr, header1.npart[k]*sizeof(float)*3, SEEK_CUR);
	      fread(&dummy, sizeof(dummy), 1, fdsfr);
	      fread(&dummy, sizeof(dummy), 1, fdsfr);
	      for(k=0; k<6; k++) /* skip id data */
		fseek(fdsfr, header1.npart[k]*sizeof(int), SEEK_CUR);
	      fread(&dummy, sizeof(dummy), 1, fdsfr);

	      fread(&dummy, sizeof(dummy), 1, fdsfr);
	      for(k=0; k<6; k++) /* skip masses */
		if(header1.mass[k]==0)
		  fseek(fdsfr, header1.npart[k]*sizeof(float), SEEK_CUR);
	      fread(&dummy, sizeof(dummy), 1, fdsfr);

	      fread(&dummy, sizeof(dummy), 1, fdsfr);  /* skip u */
	      fseek(fdsfr, header1.npart[0]*sizeof(float), SEEK_CUR);
	      fread(&dummy, sizeof(dummy), 1, fdsfr);

	      fread(&dummy, sizeof(dummy), 1, fdsfr);  /* skip rho */
	      fseek(fdsfr, header1.npart[0]*sizeof(float), SEEK_CUR);
	      fread(&dummy, sizeof(dummy), 1, fdsfr);

	      fread(&dummy, sizeof(dummy), 1, fdsfr);  /* skip ne */
	      fseek(fdsfr, header1.npart[0]*sizeof(float), SEEK_CUR);
	      fread(&dummy, sizeof(dummy), 1, fdsfr);

	      fread(&dummy, sizeof(dummy), 1, fdsfr);  /* skip nh */
	      fseek(fdsfr, header1.npart[0]*sizeof(float), SEEK_CUR);
	      fread(&dummy, sizeof(dummy), 1, fdsfr);

	      fread(&dummy, sizeof(dummy), 1, fdsfr);  /* skip mfs */
	      fseek(fdsfr, header1.npart[0]*sizeof(float), SEEK_CUR);
	      fread(&dummy, sizeof(dummy), 1, fdsfr);
      
	      if(header1.flag_multiphase>0)
		{
		  fread(&dummy, sizeof(dummy), 1, fdsfr);  /* skip cold cloud mass */
		  fseek(fdsfr, header1.npart[0]*sizeof(float), SEEK_CUR);
		  fread(&dummy, sizeof(dummy), 1, fdsfr);
		}
  
	      fread(&dummy, sizeof(dummy), 1, fdsfr);  /* skip hsml */
	      fseek(fdsfr, header1.npart[0]*sizeof(float), SEEK_CUR);
	      fread(&dummy, sizeof(dummy), 1, fdsfr);

	      fread(&dummy, sizeof(dummy), 1, fdsfr);
	    }


	  if(header1.flag_sfr && header1.flag_multiphase)
	    {
	      fdclouds=fopen(buf,"r");
	      fread(&dummy, sizeof(dummy), 1, fdclouds);
	      fread(&header1, sizeof(header1), 1, fdclouds);
	      fread(&dummy, sizeof(dummy), 1, fdclouds);
	      fread(&dummy, sizeof(dummy), 1, fdclouds);
	      for(k=0; k<6; k++) /* skip pos data */
		fseek(fdclouds, header1.npart[k]*sizeof(float)*3, SEEK_CUR);
	      fread(&dummy, sizeof(dummy), 1, fdclouds);
	      fread(&dummy, sizeof(dummy), 1, fdclouds);
	      for(k=0; k<6; k++) /* skip vel data */
		fseek(fdclouds, header1.npart[k]*sizeof(float)*3, SEEK_CUR);
	      fread(&dummy, sizeof(dummy), 1, fdclouds);
	      fread(&dummy, sizeof(dummy), 1, fdclouds);
	      for(k=0; k<6; k++) /* skip id data */
		fseek(fdclouds, header1.npart[k]*sizeof(int), SEEK_CUR);
	      fread(&dummy, sizeof(dummy), 1, fdclouds);

	      fread(&dummy, sizeof(dummy), 1, fdclouds);
	      for(k=0; k<6; k++) /* skip masses */
		if(header1.mass[k]==0)
		  fseek(fdclouds, header1.npart[k]*sizeof(float), SEEK_CUR);
	      fread(&dummy, sizeof(dummy), 1, fdclouds);

	      fread(&dummy, sizeof(dummy), 1, fdclouds);  /* skip u */
	      fseek(fdclouds, header1.npart[0]*sizeof(float), SEEK_CUR);
	      fread(&dummy, sizeof(dummy), 1, fdclouds);

	      fread(&dummy, sizeof(dummy), 1, fdclouds);  /* skip rho */
	      fseek(fdclouds, header1.npart[0]*sizeof(float), SEEK_CUR);
	      fread(&dummy, sizeof(dummy), 1, fdclouds);

	      fread(&dummy, sizeof(dummy), 1, fdclouds);  /* skip ne */
	      fseek(fdclouds, header1.npart[0]*sizeof(float), SEEK_CUR);
	      fread(&dummy, sizeof(dummy), 1, fdclouds);

	      fread(&dummy, sizeof(dummy), 1, fdclouds);  /* skip nh */
	      fseek(fdclouds, header1.npart[0]*sizeof(float), SEEK_CUR);
	      fread(&dummy, sizeof(dummy), 1, fdclouds);

	      fread(&dummy, sizeof(dummy), 1, fdclouds);  /* skip mfs */
	      fseek(fdclouds, header1.npart[0]*sizeof(float), SEEK_CUR);
	      fread(&dummy, sizeof(dummy), 1, fdclouds);
      
	      fread(&dummy, sizeof(dummy), 1, fdclouds);
	    }


	  fread(&dummy, sizeof(dummy), 1, fd);
	  fread(&header1, sizeof(header1), 1, fd);
	  fread(&dummy, sizeof(dummy), 1, fd);

	  fread(&dummy, sizeof(dummy), 1, fd);
	  for(k=0; k<6; k++)
	    {
	      for(n=0; n<header1.npart[k]; n++)
		{
		  fread(&pos[0],sizeof(float),3,fd);
		  fread(&vel[0],sizeof(float),3,fdvel);
		  fread(&particleid,sizeof(int),1,fdPID);

		  if(k==0) 
		    {
		      fread(&intenergy,sizeof(float),1,fdu);
		      fread(&rho      ,sizeof(float),1,fdrho);
		    }
		  else
		    intenergy = rho = 0;
		  
		  if(header1.mass[k]>0)
		    mass= header1.mass[k];
		  else
		    fread(&mass,sizeof(float),1,fdmass);

		  if(header1.flag_sfr && k==0) 
		    {
		      fread(&mfs,sizeof(float),1,fdmfs);
		      fread(&sfr,sizeof(float),1,fdsfr);
		      if(header1.flag_multiphase)
			fread(&mclouds,sizeof(float),1,fdclouds);
		      else
			mclouds= 0;
		    }
		  else
		    mfs=sfr=mclouds=0;
		  
		  slab= (pos[0]/BoxSize)*NTask;

		  for(j=0; j<3; j++) 
		    {
		      Pbuf[slab*BUFSIZE+Nbuf[slab]].Pos[j]= (double) pos[j];
		      Pbuf[slab*BUFSIZE+Nbuf[slab]].Vel[j]= vel[j];
		    }
		  Pbuf[slab*BUFSIZE+Nbuf[slab]].PartID= particleid;
		  Pbuf[slab*BUFSIZE+Nbuf[slab]].Mass= mass;
		  Pbuf[slab*BUFSIZE+Nbuf[slab]].Type= k;
		  Pbuf[slab*BUFSIZE+Nbuf[slab]].Mfs= mfs;
		  Pbuf[slab*BUFSIZE+Nbuf[slab]].Sfr= sfr;
		  Pbuf[slab*BUFSIZE+Nbuf[slab]].Energy= intenergy;
		  Pbuf[slab*BUFSIZE+Nbuf[slab]].Rho   = rho;
		  Pbuf[slab*BUFSIZE+Nbuf[slab]].Mclouds= mclouds;

		  Nbuf[slab]++;

		  if(slab==0)
		    {
		      Nbuf[slab]--;
		      P[1 + Nlocal++]= Pbuf[slab*BUFSIZE+Nbuf[slab]];
		    }
		  else
		    {
		      if(Nbuf[slab]==BUFSIZE)
			{
			  MPI_Ssend(&Nbuf[slab], 1, MPI_INT, slab, slab, MPI_COMM_WORLD);
			  MPI_Ssend(&Pbuf[slab*BUFSIZE], sizeof(struct particle_data)*Nbuf[slab],
				    MPI_BYTE, slab, slab, MPI_COMM_WORLD);
			  Nbuf[slab]= 0;
			}
		    }
		}
	    }
	  fread(&dummy, sizeof(dummy), 1, fd);

      	  fclose(fd);
      	  fclose(fdvel);
	  fclose(fdPID);
      	  fclose(fdmass);
	  fclose(fdu);
	  fclose(fdrho);
	  if(header1.flag_sfr) 
	    {
	      fclose(fdmfs);
	      fclose(fdsfr);
	      if(header1.flag_multiphase)
		fclose(fdclouds);
	    }
	}

      for(slab=1; slab<NTask; slab++)
	{
	  while(Nbuf[slab])
	    {
	      MPI_Ssend(&Nbuf[slab], 1, MPI_INT, slab, slab, MPI_COMM_WORLD);
	      MPI_Ssend(&Pbuf[slab*BUFSIZE], sizeof(struct particle_data)*Nbuf[slab],
			    MPI_BYTE, slab, slab, MPI_COMM_WORLD);
	      Nbuf[slab]= 0;
	    }
	  
	  MPI_Ssend(&Nbuf[slab], 1, MPI_INT, slab, slab, MPI_COMM_WORLD);
	}
    }
  else
    {
      do
	{
	  MPI_Recv(&n, 1, MPI_INT, 0, ThisTask, MPI_COMM_WORLD, &status);
	  if(n)
	    {
	      MPI_Recv(&P[1+ Nlocal], sizeof(struct particle_data)*n,
		       MPI_BYTE, 0, ThisTask, MPI_COMM_WORLD, &status);
	      Nlocal+= n;
	    }
	}
      while(n);
    }


  free(Nbuf);
  free(Pbuf);
}


int find_files(char *fname)
{
  FILE *fd;
  char buf[200], buf1[200];
  int j,dummy;


  sprintf(buf, "%s.%d", fname, 0);
  sprintf(buf1, "%s", fname);

  if((fd=fopen(buf,"r")))
    {
      fread(&dummy, sizeof(dummy), 1, fd);
      fread(&header1, sizeof(header1), 1, fd);
      fread(&dummy, sizeof(dummy), 1, fd);
      fclose(fd);

      for(j=0, NumPart=0; j<5; j++)
	NumPart+= header1.npartTotal[j];

      BoxSize=     header1.BoxSize;      
      Omega=       header1.Omega0;
      OmegaLambda= header1.OmegaLambda;
      Time=        header1.time;

      if(ThisTask==0)
	{
	  printf("BoxSize=%g  Omega0=%g  OmegaLambda=%g  Time=%g\n", 
		 BoxSize, Omega, OmegaLambda, Time);
	  printf("NumPart= %d\n", NumPart);
	}

      return header1.num_files;
    }

  if((fd=fopen(buf1,"r")))
    {
      fread(&dummy, sizeof(dummy), 1, fd);
      fread(&header1, sizeof(header1), 1, fd);
      fread(&dummy, sizeof(dummy), 1, fd);
      fclose(fd);

      for(j=0, NumPart=0; j<5; j++)
	NumPart+= header1.npart[j];

      BoxSize= header1.BoxSize;
      Omega=       header1.Omega0;
      OmegaLambda= header1.OmegaLambda;
      Time=        header1.time;

      if(ThisTask==0)
	{
	  printf("BoxSize= %g  Omega0=%g  OmegaLambda=%g  Time=%g\n", 
		 BoxSize, Omega, OmegaLambda, Time);
	  printf("NumPart= %d\n", NumPart);      
	}

      return 1;
    }
    

  printf("Error. Can't find snapshot!\nneither as `%s'\nnor as `%s'\n\n",
	 buf, buf1);
  
  exit(1);
  return 0;
}





void count_local_particles(char *fname, int files)
{
  FILE *fd;
  char buf[200];
  int i,j, k,dummy;
  int n, slab;
  float pos[3];

  Nslab=   malloc(sizeof(int)*NTask);
  Nshadow= malloc(sizeof(int)*NTask);
  Noffset= malloc(sizeof(int)*NTask);
  NtoLeft= malloc(sizeof(int)*NTask);
  NtoRight=malloc(sizeof(int)*NTask);

  for(i=0; i<NTask; i++)
    Nslab[i]= Nshadow[i]= Noffset[i]= NtoLeft[i]= NtoRight[i]= 0;

  if(ThisTask==0)
    {
      for(i=0; i<files; i++)
	{
	  if(files>1)
	    sprintf(buf,"%s.%d",fname,i);
	  else
	    sprintf(buf,"%s",fname);
	  
	  if(!(fd=fopen(buf,"r")))
	    {
	      printf("can't open file `%s`\n",buf);
	      MPI_Abort(MPI_COMM_WORLD, 1); exit(1);
	    }
	  
	  printf("reading `%s' ...\n", buf); fflush(stdout);

	  fread(&dummy, sizeof(dummy), 1, fd);
	  fread(&header1, sizeof(header1), 1, fd);
	  fread(&dummy, sizeof(dummy), 1, fd);
	  
	  fread(&dummy, sizeof(dummy), 1, fd);
	  for(k=0; k<6; k++)
	    {
	      for(n=0; n<header1.npart[k]; n++)
		{
		  fread(&pos[0], sizeof(float), 3, fd);
		  
		  slab= (pos[0]/BoxSize)*NTask;
		  
		  Nslab[slab]++;
		  
		  if(pos[0] < slab*(BoxSize/NTask)+SearchRadius)
		    {
		      NtoLeft[slab]++; 
		    }

		  if(pos[0] > (slab+1)*(BoxSize/NTask)-SearchRadius)
		    {
		      NtoRight[slab]++; 
		    }
		}
	    }
	  fread(&dummy, sizeof(dummy), 1, fd);
	  fclose(fd);
	}
    }

  MPI_Bcast(Nslab,    NTask, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(NtoLeft,  NTask, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(NtoRight, NTask, MPI_INT, 0, MPI_COMM_WORLD);

  for(i=0; i<NTask;i++)
    {
    if(i < NTask-1)
      Nshadow[i]+= NtoLeft[i+1];
    else
      Nshadow[i]+= NtoLeft[0];
    
    if(i > 0)
      Nshadow[i]+= NtoRight[i-1];
    else
      Nshadow[i]+= NtoRight[NTask-1];
    }

  for(i=0; i<NTask; i++)
    for(j=0, Noffset[i]=0; j<i; j++)
      Noffset[i]+= Nslab[j];
}




