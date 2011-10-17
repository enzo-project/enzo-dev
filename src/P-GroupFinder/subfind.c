#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include "h5utilities.h"


#include "allvars.h"
#include "proto.h"




void subfind(char *particles_fname, char *catalogue_fname, 
	     char *subhalo_fname, char *parttypes_fname, char *partids_fname, 
	     char *subprop_fname, char *prop_fname, char *sparticles_fname5,
	     char *scataloguetxt)
{
  MPI_Status status;
  FILE   *fd, *fdpart, *fdlen, *fdoffset, *fdparent, *fdtypes, *fdids;
  FILE   *fdsubcenter, *fdsubmtot, *fdsubmgas, *fdsubmstars, *fdsubsfr, *fdsubmcloud;
  FILE   *fdcenter, *fdmtot, *fdmgas, *fdmstars, *fdsfr, *fdmcloud;
  FILE   *fdtxt;
  char   buflen[500], bufoffset[500], bufparent[500], bufcount[500], command[2000];
  char   bufsubcenter[500], bufsubmtot[500], bufsubmgas[500], bufsubmstars[500], bufsubsfr[500], bufsubmcloud[500], commandsubprop[2000];
  char   bufcenter[500], bufmtot[500], bufmgas[500], bufmstars[500], bufsfr[500], bufmcloud[500], commandprop[2000];
  int    i, k, dim, gr, task, head, len, nsubs, offset, start=0, NSubGroupsAll=0;
  int    parent, ntot, _idx;
  char   ctype;
  float  cm[3], cmv[3], AM[3], mtot, mgas, mstars, sfr, mcloud, redshift, spin, vrms;
  float  mvir, rvir;
  float  corner[3];
  struct particle_data *Pbuf, *partbuf;
  int    *sublen, *suboffset, *bufsublen, *bufsuboffset;
  int    *fsuboffset, *fbufsuboffset;

  char   halo_name[200];
  PINT    *TempPINT;
  double *temp;
  float  *msub, *bufmsub;

  hid_t  file_id, dset_id, dspace_id, group_id;
  hsize_t hdims[2];

  sprintf(buflen,    "%s.len",    subhalo_fname);
  sprintf(bufoffset, "%s.offset", subhalo_fname);
  sprintf(bufparent, "%s.parent", subhalo_fname);
  sprintf(bufcount,  "%s.count",  subhalo_fname);
  sprintf(command, "cat %s %s %s %s > %s",
	  bufcount, buflen, bufoffset, bufparent, subhalo_fname);


  sprintf(bufsubcenter,"%s.subcenter", subhalo_fname);
  sprintf(bufsubmtot,  "%s.submtot",   subhalo_fname);
  sprintf(bufsubmgas,  "%s.submgas",   subhalo_fname);
  sprintf(bufsubmstars,"%s.submstars", subhalo_fname);
  sprintf(bufsubsfr,   "%s.subsfr",    subhalo_fname);
  sprintf(bufsubmcloud,"%s.submcloud", subhalo_fname);
  sprintf(commandsubprop, "cat %s %s %s %s %s %s %s > %s",
	  bufcount, bufsubcenter, bufsubmtot, bufsubmgas, bufsubmstars, bufsubsfr, bufsubmcloud, subprop_fname);


  sprintf(bufcenter,"%s.center", catalogue_fname);
  sprintf(bufmtot,  "%s.mtot",   catalogue_fname);
  sprintf(bufmgas,  "%s.mgas",   catalogue_fname);
  sprintf(bufmstars,"%s.mstars", catalogue_fname);
  sprintf(bufsfr,   "%s.sfr",    catalogue_fname);
  sprintf(bufmcloud,"%s.mcloud",    catalogue_fname);
  sprintf(commandprop, "cat %s %s %s %s %s %s > %s",
	  bufcenter, bufmtot, bufmgas, bufmstars, bufsfr, bufmcloud, prop_fname);

    
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

      if(!(fdcenter=fopen(bufcenter, "w")))
	{
	  printf("can't open file `%s`\n", bufcenter);
	  MPI_Abort(MPI_COMM_WORLD, 1); exit(1);
	}

      if(!(fdmtot=fopen(bufmtot, "w")))
	{
	  printf("can't open file `%s`\n", bufmtot);
	  MPI_Abort(MPI_COMM_WORLD, 1); exit(1);
	}

      if(!(fdmgas=fopen(bufmgas, "w")))
	{
	  printf("can't open file `%s`\n", bufmgas);
	  MPI_Abort(MPI_COMM_WORLD, 1); exit(1);
	}

      if(!(fdmstars=fopen(bufmstars, "w")))
	{
	  printf("can't open file `%s`\n", bufmstars);
	  MPI_Abort(MPI_COMM_WORLD, 1); exit(1);
	}

      if(!(fdsfr=fopen(bufsfr, "w")))
	{
	  printf("can't open file `%s`\n", bufsfr);
	  MPI_Abort(MPI_COMM_WORLD, 1); exit(1);
	}

      if(!(fdmcloud=fopen(bufmcloud, "w")))
	{
	  printf("can't open file `%s`\n", bufmcloud);
	  MPI_Abort(MPI_COMM_WORLD, 1); exit(1);
	}

      fwrite(&NgroupsAll, sizeof(int), 1, fdcenter);

      fprintf(stdout, "Saving subhalo list to %s\n", scataloguetxt);
      fprintf(stdout, "Saving (sub)halo particle list to %s\n", 
	      sparticles_fname5);

      if ((fdtxt = fopen(scataloguetxt, "w")) == NULL) {
	printf("can't open file `%s`\n", scataloguetxt);
	MPI_Abort(MPI_COMM_WORLD, 1); exit(1);
      }

      // Write header

      redshift = 1.0 / Time - 1.0;
      fprintf(fdtxt, "# Scale factor = %lg\n", Time);
      fprintf(fdtxt, "# Redshift = %f\n", redshift);
      //fprintf(fdtxt, "# Number of subhalos = %"ISYM"\n", AllVars.NgroupsAll);
      fprintf(fdtxt, "#\n");
      fprintf(fdtxt, "# Column 1.  Center of mass (x)\n");
      fprintf(fdtxt, "# Column 2.  Center of mass (y)\n");
      fprintf(fdtxt, "# Column 3.  Center of mass (z)\n");
      fprintf(fdtxt, "# Column 4.  Subhalo number\n");
      fprintf(fdtxt, "# Column 5.  Parent halo number\n");
      fprintf(fdtxt, "# Column 6.  First particle in halo particle list\n");
      fprintf(fdtxt, "#            --> All subgroup particles are consecutively listed in\n");
      fprintf(fdtxt, "#                particle list (if written)\n");
      fprintf(fdtxt, "# Column 7.  Number of particles\n");
      fprintf(fdtxt, "# Column 8.  Halo mass [solar masses]\n");
      fprintf(fdtxt, "# Column 9.  Stellar mass [solar masses]\n");
      fprintf(fdtxt, "# Column 10. Mean x-velocity [km/s]\n");
      fprintf(fdtxt, "# Column 11. Mean y-velocity [km/s]\n");
      fprintf(fdtxt, "# Column 12. Mean z-velocity [km/s]\n");
      fprintf(fdtxt, "# Column 13. Velocity dispersion [km/s]\n");
      fprintf(fdtxt, "# Column 14. Mean x-angular momentum [Mpc * km/s]\n");
      fprintf(fdtxt, "# Column 15. Mean y-angular momentum [Mpc * km/s]\n");
      fprintf(fdtxt, "# Column 16. Mean z-angular momentum [Mpc * km/s]\n");
      fprintf(fdtxt, "# Column 17. Spin parameter\n");
      fprintf(fdtxt, "#\n");
      fprintf(fdtxt, "# datavar lines are for partiview.  Ignore them if you're not partiview.\n");
      fprintf(fdtxt, "#\n");
      fprintf(fdtxt, "datavar 1 subhalo_number\n");
      fprintf(fdtxt, "datavar 2 parent_number\n");
      fprintf(fdtxt, "datavar 3 particle_offset\n");
      fprintf(fdtxt, "datavar 4 number_of_particles\n");
      fprintf(fdtxt, "datavar 5 halo_mass\n");
      fprintf(fdtxt, "datavar 6 stellar_mass\n");
      fprintf(fdtxt, "datavar 7 x_velocity\n");
      fprintf(fdtxt, "datavar 8 y_velocity\n");
      fprintf(fdtxt, "datavar 9 z_velocity\n");
      fprintf(fdtxt, "datavar 10 velocity_dispersion\n");
      fprintf(fdtxt, "datavar 11 x_angular_momentum\n");
      fprintf(fdtxt, "datavar 12 y_angular_momentum\n");
      fprintf(fdtxt, "datavar 13 z_angular_momentum\n");
      fprintf(fdtxt, "datavar 14 spin\n");
      fprintf(fdtxt, "\n");

      file_id = H5Fcreate(sparticles_fname5, H5F_ACC_TRUNC, H5P_DEFAULT, 
			  H5P_DEFAULT);
      group_id = H5Gcreate(file_id, "/Parameters", 0);
      writeScalarAttribute(group_id, H5T_NATIVE_FLOAT, "Redshift", &redshift);
      writeScalarAttribute(group_id, H5T_NATIVE_DOUBLE, "Scale Factor", &Time);
      writeScalarAttribute(group_id, H5T_NATIVE_INT, "Number of groups", &NgroupsAll);
      H5Gclose(group_id);

    } // ENDIF TASK 0


  if(ThisTask==0)
    {
      printf("finding subhaloes...\n");
  
      if(!(fdpart=fopen(particles_fname,"w")))
	{
	  printf("can't open file `%s`\n", particles_fname);
	  MPI_Abort(MPI_COMM_WORLD, 1); exit(1);
	}

      if(!(fdtypes=fopen(parttypes_fname,"w")))
	{
	  printf("can't open file `%s`\n", parttypes_fname);
	  MPI_Abort(MPI_COMM_WORLD, 1); exit(1);
	}

      if(!(fdids=fopen(partids_fname,"w")))
	{
	  printf("can't open file `%s`\n", partids_fname);
	  MPI_Abort(MPI_COMM_WORLD, 1); exit(1);
	}

      for(gr= NgroupsAll-1, ntot=0; gr>=0; gr--)
	ntot+= GroupDatAll[gr].Len;

      fwrite(&ntot, sizeof(int), 1, fdpart);
      fwrite(&ntot, sizeof(int), 1, fdtypes);
      fwrite(&ntot, sizeof(int), 1, fdids);


      if(!(fdlen=fopen(buflen, "w")))
	{
	  printf("can't open file `%s`\n", buflen);
	  MPI_Abort(MPI_COMM_WORLD, 1); exit(1);
	}

      if(!(fdoffset=fopen(bufoffset, "w")))
	{
	  printf("can't open file `%s`\n", bufoffset);
	  MPI_Abort(MPI_COMM_WORLD, 1); exit(1);
	}

      if(!(fdparent=fopen(bufparent, "w")))
	{
	  printf("can't open file `%s`\n", bufparent);
	  MPI_Abort(MPI_COMM_WORLD, 1); exit(1);
	}

      if(!(fdsubcenter=fopen(bufsubcenter, "w")))
	{
	  printf("can't open file `%s`\n", bufsubcenter);
	  MPI_Abort(MPI_COMM_WORLD, 1); exit(1);
	}

      if(!(fdsubmtot=fopen(bufsubmtot, "w")))
	{
	  printf("can't open file `%s`\n", bufsubmtot);
	  MPI_Abort(MPI_COMM_WORLD, 1); exit(1);
	}

      if(!(fdsubmgas=fopen(bufsubmgas, "w")))
	{
	  printf("can't open file `%s`\n", bufsubmgas);
	  MPI_Abort(MPI_COMM_WORLD, 1); exit(1);
	}

      if(!(fdsubmstars=fopen(bufsubmstars, "w")))
	{
	  printf("can't open file `%s`\n", bufsubmstars);
	  MPI_Abort(MPI_COMM_WORLD, 1); exit(1);
	}

      if(!(fdsubsfr=fopen(bufsubsfr, "w")))
	{
	  printf("can't open file `%s`\n", bufsubsfr);
	  MPI_Abort(MPI_COMM_WORLD, 1); exit(1);
	}

      if(!(fdsubmcloud=fopen(bufsubmcloud, "w")))
	{
	  printf("can't open file `%s`\n", bufsubmcloud);
	  MPI_Abort(MPI_COMM_WORLD, 1); exit(1);
	}
    }



  for(gr=NgroupsAll-1; gr>=0; )
   {
      for(task=0; task<NTask && (gr-task)>=0; task++)
	{
	  if(ThisTask==task)
	    {
	      head= GroupDatAll[gr-task].Tag;
	      len=  GroupDatAll[gr-task].Len;
	      sublen=    mymalloc(len/DesLinkNgb*sizeof(int));
	      suboffset= mymalloc(len/DesLinkNgb*sizeof(int));
	      fsuboffset= mymalloc(len/DesLinkNgb*sizeof(int));
	      msub= mymalloc(len/DesLinkNgb*sizeof(float));
	      Pbuf=      mymalloc(sizeof(struct particle_data)*len);
	    }
	  get_particles(task, head, len, Pbuf);
	}

      if((gr-ThisTask)>=0) /* ok, this process got a group */
	{
	  len=   GroupDatAll[gr-ThisTask].Len;

	  for(k=0; k<3; k++)
	    corner[k]= Pbuf[0].Pos[k];
	  
	  for(i=0; i<len; i++)
	    for(k=0; k<3; k++)
	      Pbuf[i].Pos[k]= periodic(Pbuf[i].Pos[k]-corner[k]);

	  nsubs= do_subfind_in_group(Pbuf, len, sublen, suboffset);

	  for(i=0; i<len; i++)
	    for(k=0; k<3; k++)
	      Pbuf[i].Pos[k]= periodic_wrap(Pbuf[i].Pos[k]+corner[k]);
	}

      for(task=0; task<NTask && (gr-task)>=0; task++)
	{
	  parent = NgroupsAll-(gr-task)-1;
	  if(ThisTask==0)
	    {
	      if(task==0)
		{
		  for(i=0; i<nsubs; i++)
		    {
		      get_properties(Pbuf+suboffset[i], sublen[i], cm, 
				     &mtot, &mgas, &mstars, &sfr, &mcloud,
				     1, &cmv[0], &mvir, &rvir, AM, &vrms,
				     &spin);
		      msub[i] = mtot;
		      fwrite(cm,      sizeof(float), 3, fdsubcenter);
		      fwrite(&mtot,   sizeof(float), 1, fdsubmtot);
		      fwrite(&mgas,   sizeof(float), 1, fdsubmgas);
		      fwrite(&mstars, sizeof(float), 1, fdsubmstars);
		      fwrite(&sfr,    sizeof(float), 1, fdsubsfr);
		      fwrite(&mcloud,    sizeof(float), 1, fdsubmcloud);

		      /* Write to ASCII halo catalog */

		      fprintf(fdtxt, "%12.6g %12.6g %12.6g %12d %12d %12d %12d "
			      "%12.6g %12.6g %12.6g %12.6g %12.6g %12.6g %12.6g "
			      "%12.6g %12.6g %12.6g\n",
			      cm[0], cm[1], cm[2], i, parent, suboffset[i], sublen[i], 
			      mtot, mstars, cmv[0], cmv[1], cmv[2], vrms, AM[0], 
			      AM[1], AM[2], spin);

		    }

		  for(i=0; i<nsubs; i++) {
		    fsuboffset[i] = suboffset[i];
		    suboffset[i]+= start;
		  }
		  start+= GroupDatAll[gr-task].Len;

		  len = GroupDatAll[gr-task].Len;
		  temp = (double*) malloc(3*len*sizeof(double));
		  TempPINT = (PINT*) malloc(len*sizeof(PINT));
		  _idx = 0;
		  for (dim = 0; dim < 3; dim++)
		    for (i = 0; i < len; i++, _idx++)
		      temp[_idx] = Pbuf[i].Pos[dim] / BoxSize;
		  for (i = 0; i < len; i++)
		    TempPINT[i] = Pbuf[i].PartID;

		  get_properties(Pbuf, GroupDatAll[gr-task].Len, cm, &mtot, 
				 &mgas, &mstars, &sfr, &mcloud, 1,
				 &cmv[0], &mvir, &rvir, AM, &vrms, &spin);
		  fwrite(cm,      sizeof(float), 3, fdcenter);
		  fwrite(&mtot,   sizeof(float), 1, fdmtot);
		  fwrite(&mgas,   sizeof(float), 1, fdmgas);
		  fwrite(&mstars, sizeof(float), 1, fdmstars);
		  fwrite(&sfr,    sizeof(float), 1, fdsfr);
		  fwrite(&mcloud,    sizeof(float), 1, fdmcloud);

		  /* Write to HDF5 particle list */

		  sprintf(halo_name, "Halo%8.8d", parent);
		  group_id = H5Gcreate(file_id, halo_name, 0);
		  writeScalarAttribute(group_id, H5T_NATIVE_INT, "NumberOfSubhalos", &nsubs);
		  writeScalarAttribute(group_id, H5T_NATIVE_FLOAT, "Total Mass", &mtot);
		  writeScalarAttribute(group_id, H5T_NATIVE_FLOAT, "Stellar Mass", &mstars);
		  writeScalarAttribute(group_id, H5T_NATIVE_FLOAT, "Spin parameter", &spin);
		  writeScalarAttribute(group_id, H5T_NATIVE_FLOAT, "Velocity dispersion", &vrms);
		  writeArrayAttribute(group_id, H5T_NATIVE_FLOAT, 3, "Center of mass", cm);
		  writeArrayAttribute(group_id, H5T_NATIVE_FLOAT, 3, "Mean velocity [km/s]", cmv);
		  writeArrayAttribute(group_id, H5T_NATIVE_FLOAT, 3, "Angular momentum [Mpc * km/s]", AM);

		  if (nsubs > 0) {
		    hdims[0] = (hsize_t) nsubs;
		    hdims[1] = 1;
		    dspace_id = H5Screate_simple(1, hdims, NULL);
		    dset_id = H5Dcreate(group_id, "Subhalo Mass", H5T_NATIVE_FLOAT, dspace_id,
					H5P_DEFAULT);
		    H5Dwrite(dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
			     (void*) msub);
		    H5Sclose(dspace_id);
		    H5Dclose(dset_id);

		    dspace_id = H5Screate_simple(1, hdims, NULL);
		    dset_id = H5Dcreate(group_id, "Subhalo Size", H5T_NATIVE_INT, dspace_id,
					H5P_DEFAULT);
		    H5Dwrite(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
			     (void*) sublen);
		    H5Sclose(dspace_id);
		    H5Dclose(dset_id);

		    dspace_id = H5Screate_simple(1, hdims, NULL);
		    dset_id = H5Dcreate(group_id, "Subhalo Offset", H5T_NATIVE_INT, dspace_id,
					H5P_DEFAULT);
		    H5Dwrite(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
			     (void*) fsuboffset);
		    H5Sclose(dspace_id);
		    H5Dclose(dset_id);

		  } // ENDIF nsubs>0

		  hdims[0] = 3;
		  hdims[1] = (hsize_t) len;
		  dspace_id = H5Screate_simple(2, hdims, NULL);
		  dset_id = H5Dcreate(group_id, "Particle Position", H5T_NATIVE_DOUBLE,
				      dspace_id, H5P_DEFAULT);
		  H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
			   H5P_DEFAULT, (void*) temp);
		  H5Sclose(dspace_id);
		  H5Dclose(dset_id);
	
		  hdims[0] = (hsize_t) len;
		  hdims[1] = 1;
		  dspace_id = H5Screate_simple(1, hdims, NULL);
		  dset_id = H5Dcreate(group_id, "Particle ID", HDF5_PINT, dspace_id,
				      H5P_DEFAULT);
		  H5Dwrite(dset_id, HDF5_PINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
			   (void*) TempPINT);
		  H5Sclose(dspace_id);
		  H5Dclose(dset_id);

		  H5Gclose(group_id);

		  free(temp);
		  free(TempPINT);

		  /* END HDF5 writing */

		  fwrite(sublen, sizeof(int), nsubs, fdlen);
		  fwrite(suboffset, sizeof(int), nsubs, fdoffset);
		  for(i=0, parent=NgroupsAll-(gr-task); i<nsubs; i++)
		    fwrite(&parent, sizeof(int), 1, fdparent);
		  for(i=0; i<GroupDatAll[gr-task].Len; i++)
		    fwrite(&Pbuf[i].Pos[0], sizeof(double), 3, fdpart);
		  for(i=0; i<GroupDatAll[gr-task].Len; i++)
		    fwrite(&Pbuf[i].PartID, sizeof(PINT), 1, fdids);
		  for(i=0; i<GroupDatAll[gr-task].Len; i++)
		    {
		      ctype= Pbuf[i].Type;
		      fwrite(&ctype, sizeof(char), 1, fdtypes);
		    }

		  fwrite(&nsubs, sizeof(int), 1, fd);

		  NSubGroupsAll+= nsubs;
		}
	      else
		{
		  MPI_Recv(&nsubs, 1, MPI_INT, task, task, MPI_COMM_WORLD, &status);
		  if(nsubs)
		    {
		      bufmsub=      mymalloc(nsubs*sizeof(float));
		      bufsublen=    mymalloc(nsubs*sizeof(int));
		      bufsuboffset= mymalloc(nsubs*sizeof(int));
		      fbufsuboffset= mymalloc(nsubs*sizeof(int));
		      MPI_Recv(bufsublen,    nsubs, MPI_INT, task, task, MPI_COMM_WORLD, &status);
		      MPI_Recv(bufsuboffset, nsubs, MPI_INT, task, task, MPI_COMM_WORLD, &status);
		      for(i=0; i<nsubs; i++) {
			fbufsuboffset[i] = bufsuboffset[i];
			bufsuboffset[i]+= start;
		      }
		      fwrite(bufsublen,    sizeof(int), nsubs, fdlen);
		      fwrite(bufsuboffset, sizeof(int), nsubs, fdoffset);
		    }

		  for(i=0, parent=NgroupsAll-(gr-task); i<nsubs; i++)
		    fwrite(&parent, sizeof(int), 1, fdparent);
		  
		  partbuf= mymalloc(sizeof(struct particle_data)*GroupDatAll[gr-task].Len);
		  MPI_Recv(partbuf, GroupDatAll[gr-task].Len*sizeof(struct particle_data), 
			   MPI_BYTE, task, task, MPI_COMM_WORLD, &status);
		  for(i=0; i<GroupDatAll[gr-task].Len; i++)
		    fwrite(&partbuf[i].Pos[0], sizeof(double), 3, fdpart);
		  for(i=0; i<GroupDatAll[gr-task].Len; i++)
		    fwrite(&partbuf[i].PartID, sizeof(PINT), 1, fdids);
		  for(i=0; i<GroupDatAll[gr-task].Len; i++)
		    {
		      ctype= partbuf[i].Type;
		      fwrite(&ctype, sizeof(char), 1, fdtypes);
		    }

		  fwrite(&nsubs, sizeof(int), 1, fd);
		  
		  for(i=0; i<nsubs; i++)
		    {
		      get_properties(partbuf+bufsuboffset[i]-start, bufsublen[i], cm, &mtot, 
				     &mgas, &mstars, &sfr, &mcloud,
				     1, cmv, &mvir, &rvir, AM, &vrms, &spin);
		      bufmsub[i] = mtot;
		      fwrite(cm,      sizeof(float), 3, fdsubcenter);
		      fwrite(&mtot,   sizeof(float), 1, fdsubmtot);
		      fwrite(&mgas,   sizeof(float), 1, fdsubmgas);
		      fwrite(&mstars, sizeof(float), 1, fdsubmstars);
		      fwrite(&sfr,    sizeof(float), 1, fdsubsfr);
		      fwrite(&mcloud,    sizeof(float), 1, fdsubmcloud);

		      /* Write to ASCII halo catalog */

		      fprintf(fdtxt, "%12.6g %12.6g %12.6g %12d %12d %12d %12d "
			      "%12.6g %12.6g %12.6g %12.6g %12.6g %12.6g %12.6g "
			      "%12.6g %12.6g %12.6g\n",
			      cm[0], cm[1], cm[2], i, parent, suboffset[i], sublen[i], 
			      mtot, mstars, cmv[0], cmv[1], cmv[2], vrms, AM[0], 
			      AM[1], AM[2], spin);

		    }

		  get_properties(partbuf, GroupDatAll[gr-task].Len, cm, &mtot, 
				 &mgas, &mstars, &sfr, &mcloud, 1,
				 cmv, &mvir, &rvir, AM, &vrms, &spin);
		  fwrite(cm,      sizeof(float), 3, fdcenter);
		  fwrite(&mtot,   sizeof(float), 1, fdmtot);
		  fwrite(&mgas,   sizeof(float), 1, fdmgas);
		  fwrite(&mstars, sizeof(float), 1, fdmstars);
		  fwrite(&sfr,    sizeof(float), 1, fdsfr);
		  fwrite(&mcloud,    sizeof(float), 1, fdmcloud);

		  /* Write to HDF5 particle list */

		  len = GroupDatAll[gr-task].Len;
		  parent = NgroupsAll-(gr-task)-1;
		  temp = (double*) malloc(3*len*sizeof(double));
		  TempPINT = (PINT*) malloc(len*sizeof(PINT));
		  _idx = 0;
		  for (dim = 0; dim < 3; dim++)
		    for (i = 0; i < len; i++, _idx++)
		      temp[_idx] = Pbuf[i].Pos[dim] / BoxSize;
		  for (i = 0; i < len; i++)
		    TempPINT[i] = Pbuf[i].PartID;

		  sprintf(halo_name, "Halo%8.8d", parent);
		  group_id = H5Gcreate(file_id, halo_name, 0);
		  writeScalarAttribute(group_id, H5T_NATIVE_INT, "NumberOfSubhalos", &nsubs);
		  writeScalarAttribute(group_id, H5T_NATIVE_FLOAT, "Total Mass", &mtot);
		  writeScalarAttribute(group_id, H5T_NATIVE_FLOAT, "Stellar Mass", &mstars);
		  writeScalarAttribute(group_id, H5T_NATIVE_FLOAT, "Spin parameter", &spin);
		  writeScalarAttribute(group_id, H5T_NATIVE_FLOAT, "Velocity dispersion", &vrms);
		  writeArrayAttribute(group_id, H5T_NATIVE_FLOAT, 3, "Center of mass", cm);
		  writeArrayAttribute(group_id, H5T_NATIVE_FLOAT, 3, "Mean velocity [km/s]", cmv);
		  writeArrayAttribute(group_id, H5T_NATIVE_FLOAT, 3, "Angular momentum [Mpc * km/s]", AM);

		  if (nsubs > 0) {
		    hdims[0] = (hsize_t) nsubs;
		    hdims[1] = 1;
		    dspace_id = H5Screate_simple(1, hdims, NULL);
		    dset_id = H5Dcreate(group_id, "Subhalo Mass", H5T_NATIVE_FLOAT, dspace_id,
					H5P_DEFAULT);
		    H5Dwrite(dset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
			     (void*) bufmsub);
		    H5Sclose(dspace_id);
		    H5Dclose(dset_id);

		    dspace_id = H5Screate_simple(1, hdims, NULL);
		    dset_id = H5Dcreate(group_id, "Subhalo Size", H5T_NATIVE_INT, dspace_id,
					H5P_DEFAULT);
		    H5Dwrite(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
			     (void*) sublen);
		    H5Sclose(dspace_id);
		    H5Dclose(dset_id);

		    dspace_id = H5Screate_simple(1, hdims, NULL);
		    dset_id = H5Dcreate(group_id, "Subhalo Offset", H5T_NATIVE_INT, dspace_id,
					H5P_DEFAULT);
		    H5Dwrite(dset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
			     (void*) fsuboffset);
		    H5Sclose(dspace_id);
		    H5Dclose(dset_id);

		  } // ENDIF nsubs>0

		  hdims[0] = 3;
		  hdims[1] = (hsize_t) len;
		  dspace_id = H5Screate_simple(2, hdims, NULL);
		  dset_id = H5Dcreate(group_id, "Particle Position", H5T_NATIVE_DOUBLE,
				      dspace_id, H5P_DEFAULT);
		  H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, 
			   H5P_DEFAULT, (void*) temp);
		  H5Sclose(dspace_id);
		  H5Dclose(dset_id);
	
		  hdims[0] = (hsize_t) len;
		  hdims[1] = 1;
		  dspace_id = H5Screate_simple(1, hdims, NULL);
		  dset_id = H5Dcreate(group_id, "Particle ID", HDF5_PINT, dspace_id,
				      H5P_DEFAULT);
		  H5Dwrite(dset_id, HDF5_PINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, 
			   (void*) TempPINT);
		  H5Sclose(dspace_id);
		  H5Dclose(dset_id);

		  H5Gclose(group_id);

		  free(temp);
		  free(TempPINT);

		  /* END HDF5 writing */

		  start+= GroupDatAll[gr-task].Len;

		  if(nsubs)
		    {
		      free(bufmsub);
		      free(bufsuboffset);
		      free(fbufsuboffset);
		      free(bufsublen);
		    }
		  free(partbuf);

		  NSubGroupsAll+= nsubs; 
		}
	    } // ENDIF ThisTask==0
	  else
	    {
	      if(task==ThisTask) 
		{
		  MPI_Ssend(&nsubs, 1, MPI_INT, 0, ThisTask, MPI_COMM_WORLD);
		  if(nsubs)
		    {
		      MPI_Ssend(sublen,    nsubs, MPI_INT, 0, ThisTask, MPI_COMM_WORLD);
		      MPI_Ssend(suboffset, nsubs, MPI_INT, 0, ThisTask, MPI_COMM_WORLD);
		    }
		  MPI_Ssend(Pbuf, GroupDatAll[gr-task].Len*sizeof(struct particle_data), 
			    MPI_BYTE, 0, ThisTask, MPI_COMM_WORLD);
		}
	    }
	}

      for(task=0; task<NTask && (gr-task)>=0; task++)
	{
	  if(ThisTask==task)
	    {
	      free(Pbuf);
	      free(suboffset);
	      free(fsuboffset);
	      free(sublen);
	      free(msub);
	    }
	}

      gr-=task;
    }

  if(ThisTask==0)
    {
      H5Fclose(file_id);
      fclose(fdtxt);
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

      if(!(fd=fopen(bufcount,"w")))
	{
	  printf("can't open file `%s`\n", bufcount);
	  MPI_Abort(MPI_COMM_WORLD, 1); exit(1);
	}
      fwrite(&NSubGroupsAll, sizeof(int), 1, fd);
      printf("NSubGroupsAll= %d\n", NSubGroupsAll);
      fclose(fd);
      
      system(command);
      system(commandsubprop);
      system(commandprop);

      printf("%s \n", commandsubprop);
    
      unlink(bufcount);
      unlink(buflen);
      unlink(bufoffset);
      unlink(bufparent);
      unlink(bufsubcenter);
      unlink(bufsubmtot);
      unlink(bufsubmgas);
      unlink(bufsubmstars);
      unlink(bufsubsfr);
      unlink(bufsubmcloud);
      unlink(bufcenter);
      unlink(bufmtot);
      unlink(bufmgas);
      unlink(bufmstars);
      unlink(bufsfr);
      unlink(bufmcloud);

      printf("done.\n");
    }
}





















int do_subfind_in_group(struct particle_data *pbuf, int grlen, int *sublen, int *suboffset)
{
  int *Len_bak, *Head_bak, *Tail_bak, *Next_bak;
  struct particle_data *P_bak;
  int i, j, oldhead, gr;
  int saved;
  PINT offset, count, id;



  Len_bak=   Len;
  Head_bak=  Head;
  Tail_bak=  Len;
  Next_bak=  Next;
  P_bak=     P;


  P=          pbuf-1;
  NumInGroup= grlen;

  Energy=    vector(1,  NumInGroup);
  Density=   vector(1,  NumInGroup);
  Potential= vector(1,  NumInGroup);
  Next=      ivector(1, NumInGroup);
  Head=      ivector(1, NumInGroup);
  NewHead=   ivector(1, NumInGroup);
  Tail=      ivector(1, NumInGroup);
  Len=       ivector(1, NumInGroup);
  Index=     ivector(1, NumInGroup);
  Node=      ivector(1, NumInGroup);      
  SortId=    ivector(1, NumInGroup);      
  
  GroupTree= mymalloc((MaxNodes=(NumInGroup/10))*sizeof(struct grouptree_data));
  
  for(i=1;i<=NumInGroup;i++)  /* initially, there are no subgroups */
    {
      Head[i]=Tail[i]=0;
      Next[i]=0;
      Len[i] =0;
      Node[i]=0;
    }
	  
  if(NumInGroup >= 2*DesLinkNgb && NumInGroup>DesDensityNgb)
    {
      ngb_treeallocate(NumInGroup, 2.0*NumInGroup + 200);
      ngb_treebuild(NumInGroup);
      
      density();        /* compute local densities */
      
      find_subgroups();
      
      ngb_treefree();
    }
  else
    AnzNodes=0;
   
  walk_tree_and_unbind();
  
  iindexx(NumInGroup, Head, Index);
  
  for(i=1, NSubGroups=0, oldhead=0; i<=NumInGroup; i++)
    {
      if(Head[Index[i]] != oldhead)
	{
	  oldhead= Head[Index[i]];
	  if(oldhead!=1) /* exclude background group */
	    NSubGroups++;
	}
    }

  if(ThisTask==0)
    printf("---  %d subgroups found. ----\n", NSubGroups);
  
  if(NSubGroups>0)
    {
      /* determine sizes and tags of subgroups */
	  
      SubGroupLen= ivector(1, NSubGroups); 
      SubGroupTag= ivector(1, NSubGroups);
      
      for(i=1, gr=0, oldhead=0; i<=NumInGroup; i++)
	{
	  if(Head[Index[i]] != oldhead)
	    {
	      oldhead= Head[Index[i]];
	      if(oldhead!=1) /* exclude background group */
		{
		  gr++;
		  SubGroupTag[gr]= i; /* first in that group */
		  SubGroupLen[gr]= 0;
		}
	    }
	  if(oldhead!=1) /* exclude background group */
	    SubGroupLen[gr]++;
	}
      
      sort2_int(NSubGroups, SubGroupLen, SubGroupTag);   /* order the subgroups by len */
      
      for(i=1; i<=NumInGroup; i++)
	Head[i]= Index[i];  /* index will be needed in order_subgroups_potential */
      
      order_subgroups_by_potential(); 

      for(i=NSubGroups, saved=0, offset=0, count=0; i>=1; i--) 
	{
	  for(j=0; j<SubGroupLen[i]; j++)
	    {
	      id= Head[SubGroupTag[i] + j];
	      
	      P[id].MinID= count++;  /* can be set newly at this point */
	    }
	  
	  if(SubGroupLen[i] >= DesLinkNgb)
	    {
	      sublen[saved]=    SubGroupLen[i];
	      suboffset[saved]= offset;
	      
	      offset+= SubGroupLen[i];
	      saved++;
	    }
	}

      /* sort the particles */
      
      qsort(P+1, NumInGroup, sizeof(struct particle_data), comp_func_partminid);
      

      free_ivector(SubGroupTag, 1, NSubGroups);	  
      free_ivector(SubGroupLen, 1, NSubGroups); 

      NSubGroups= saved;
    }
  
  free(GroupTree);
  free_ivector(SortId, 1, NumInGroup);      
  free_ivector(Node, 1, NumInGroup);      
  free_ivector(Index, 1, NumInGroup);
  free_ivector(Len, 1, NumInGroup);
  free_ivector(Tail, 1, NumInGroup);
  free_ivector(NewHead, 1, NumInGroup);
  free_ivector(Head, 1, NumInGroup);
  free_ivector(Next, 1, NumInGroup);
  free_vector(Potential, 1,  NumInGroup);
  free_vector(Density, 1,  NumInGroup);
  free_vector(Energy, 1,  NumInGroup);
  
  Len=  Len_bak;
  Head= Head_bak;
  Tail= Tail_bak;
  Next= Next_bak;
  P=    P_bak;


  return NSubGroups;
}













