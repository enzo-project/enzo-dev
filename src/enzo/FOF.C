/***********************************************************************
/
/  INLINE HALO FINDER (P-GROUPFINDER)
/
/  originally written by: Volker Springel
/
/  modified1: July, 2005 by John Wise (use with Enzo output)
/  modified2: June, 2009 by John Wise (put inside Enzo)
/
************************************************************************/

#ifdef USE_MPI
#include "mpi.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <hdf5.h>
#include "h5utilities.h"
#include "performance.h"

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"

#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "TopGridData.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"

#include "FOF_allvars.h"
#include "FOF_nrutil.h"
#include "FOF_proto.h"

/************************************************************************/
void FOF_Initialize(TopGridData *MetaData, 
		    LevelHierarchyEntry *LevelArray[], 
		    FOFData &D, bool SmoothData);
void FOF_Finalize(FOFData &D, LevelHierarchyEntry *LevelArray[], 
		  TopGridData *MetaData, int FOFOnly);
/************************************************************************/

int FOF(TopGridData *MetaData, LevelHierarchyEntry *LevelArray[], 
	int WroteData, int FOFOnly)
{

  if (!InlineHaloFinder)
    return SUCCESS;

  // Always run if we've just outputted and HaloFinderRunAfterOutput is on.
  if (!(WroteData && HaloFinderRunAfterOutput)) {

    if (HaloFinderCycleSkip == 0 && HaloFinderTimestep < 0)
      return SUCCESS;

    if (HaloFinderCycleSkip > 0)
      if (MetaData->CycleNumber % HaloFinderCycleSkip != 0)
	return SUCCESS;

    if (HaloFinderTimestep > 0 &&
	MetaData->Time - HaloFinderLastTime < HaloFinderTimestep)
      return SUCCESS;

    if (NumberOfProcessors & 1 && NumberOfProcessors > 1) {
      fprintf(stdout, "FOF: Number of processors (in parallel) must be "
	      "EVEN to run inline halo finder.  Turning OFF.\n");
      InlineHaloFinder = FALSE;
      return SUCCESS;
    }

  } // ENDIF force run

  LCAPERF_START("InlineHaloFinder");

  if (!ComovingCoordinates)
    fprintf(stdout, "FOF: Warning -- you're running the halo finder on a"
	    "non-cosmology run!\n");

  FOFData AllVars;

  // in terms of mean interparticle seperation
  AllVars.LinkLength = HaloFinderLinkingLength;

  // store only groups in the catalogue with at least this number of particles
  AllVars.GroupMinLen = HaloFinderMinimumSize;

  set_units(AllVars);

  /* dimension of coarse grid. Note: the actual size of a mesh cell
     will usually be set to its optimal size, i.e. equal to the
     linking distance. The coarse grid will then be shifted around to
     cover the full volume */

  AllVars.Grid = 256;
  AllVars.MaxPlacement = 4;
  
  /* Initialization :: copy particle data from enzo's grids to the
     halo finder's data structure. */

  FOF_Initialize(MetaData, LevelArray, AllVars, false);

  marking(AllVars);

  if (NumberOfProcessors > 1)
    exchange_shadow(AllVars, MetaData->TopGridDims[0], false);

  init_coarse_grid(AllVars);
  
  link_local_slab(AllVars);
    
  if (NumberOfProcessors > 1)
    do {
      find_minids(AllVars);
    } while (link_across(AllVars) > 0);
  
  find_minids(AllVars);
  
  if (NumberOfProcessors > 1)
    stitch_together(AllVars);

  compile_group_catalogue(AllVars);

  save_groups(AllVars, MetaData->CycleNumber, MetaData->Time);
  if (HaloFinderSubfind)
    subfind(AllVars, MetaData->CycleNumber, MetaData->Time);

  FOF_Finalize(AllVars, LevelArray, MetaData, FOFOnly);
  deallocate_all_memory(AllVars);

  HaloFinderLastTime = MetaData->Time;

  LCAPERF_STOP("InlineHaloFinder");

}

/************************************************************************
   FOF and I/O SUBROUTINES
 ************************************************************************/

void set_units(FOFData &AllVars) {

  AllVars.UnitLength_in_cm	   = 3.085678e21; 
  AllVars.UnitMass_in_g		   = 1.989e43;
  AllVars.UnitVelocity_in_cm_per_s = 1.0e5;

  AllVars.Theta = 0.8;  /* opening angle for potential computation */
  AllVars.DesDensityNgb = 32;

  /* DesLinkNgb is also the minimum size of a subgroup */
  AllVars.DesLinkNgb = AllVars.DesDensityNgb;

  if (AllVars.DesDensityNgb > AllVars.GroupMinLen) {
    if (MyProcessorNumber == ROOT_PROCESSOR)
    ENZO_FAIL("Must have DesDensityNgb <= GroupMinLen\n");
  }

  AllVars.UnitTime_in_s	= AllVars.UnitLength_in_cm / 
    AllVars.UnitVelocity_in_cm_per_s;

  AllVars.UnitTime_in_Megayears	 = AllVars.UnitTime_in_s / SEC_PER_MEGAYEAR;

  AllVars.G = GRAVITY / pow(AllVars.UnitLength_in_cm,3) * 
    AllVars.UnitMass_in_g * pow(AllVars.UnitTime_in_s,2);

  AllVars.UnitDensity_in_cgs = AllVars.UnitMass_in_g / 
    pow(AllVars.UnitLength_in_cm,3);

  AllVars.UnitPressure_in_cgs = AllVars.UnitMass_in_g / AllVars.UnitLength_in_cm / 
    pow(AllVars.UnitTime_in_s,2);

  AllVars.UnitCoolingRate_in_cgs = AllVars.UnitPressure_in_cgs / 
    AllVars.UnitTime_in_s;

  AllVars.UnitEnergy_in_cgs = AllVars.UnitMass_in_g * 
    pow(AllVars.UnitLength_in_cm,2) / pow(AllVars.UnitTime_in_s,2);

  AllVars.H0 = HUBBLE * AllVars.UnitTime_in_s;
}

/************************************************************************/

void save_groups(FOFData &AllVars, int CycleNumber, FLOAT EnzoTime)
{

  FILE   *fd;
  hid_t  file_id, dset_id, dspace_id, group_id;
  hsize_t hdims[2];
  herr_t status;

  int    i, gr, offset, dim, index;
  int    ntot;
  int    head, len;
  char   ctype;
  float cm[3], cmv[3], AM[3], vrms, spin, mtot, mstars, redshift;
  float mvir, rvir;
  double *temp;
  PINT   *TempPINT;

  FOF_particle_data *Pbuf;
  char   *FOF_dirname = "FOF";
  char   catalogue_fname[200];
  char   particle_fname[200];
  char   halo_name[200];

  sprintf(catalogue_fname, "%s/groups_%5.5d.dat", FOF_dirname, CycleNumber);
  sprintf(particle_fname, "%s/particles_%5.5d.h5", FOF_dirname, CycleNumber);

  if (MyProcessorNumber == ROOT_PROCESSOR) {

    if (debug)
      fprintf(stdout, "FOF: Saving halo list to %s\n", catalogue_fname);

    if ((fd = fopen(catalogue_fname, "w")) == NULL)
      ENZO_FAIL("Unable to open FOF group file.");

    // Write header

    redshift = 1.0 / AllVars.Time - 1.0;
    fprintf(fd, "# Time     = %"PSYM"\n", EnzoTime);
    fprintf(fd, "# Redshift = %"FSYM"\n", redshift);
    fprintf(fd, "# Number of halos = %"ISYM"\n", AllVars.NgroupsAll);
    fprintf(fd, "#\n");
    fprintf(fd, "# Column 1.  Center of mass (x)\n");
    fprintf(fd, "# Column 2.  Center of mass (y)\n");
    fprintf(fd, "# Column 3.  Center of mass (z)\n");
    fprintf(fd, "# Column 4.  Halo number\n");
    fprintf(fd, "# Column 5.  Number of particles\n");
    fprintf(fd, "# Column 6.  Halo mass [solar masses]\n");
    fprintf(fd, "# Column 7.  Virial mass [solar masses]\n");
    fprintf(fd, "# Column 8.  Stellar mass [solar masses]\n");
    fprintf(fd, "# Column 9.  Virial radius (r200) [kpc]\n");
    fprintf(fd, "# Column 10. Mean x-velocity [km/s]\n");
    fprintf(fd, "# Column 11. Mean y-velocity [km/s]\n");
    fprintf(fd, "# Column 12. Mean z-velocity [km/s]\n");
    fprintf(fd, "# Column 13. Velocity dispersion [km/s]\n");
    fprintf(fd, "# Column 14. Mean x-angular momentum [Mpc * km/s]\n");
    fprintf(fd, "# Column 15. Mean y-angular momentum [Mpc * km/s]\n");
    fprintf(fd, "# Column 16. Mean z-angular momentum [Mpc * km/s]\n");
    fprintf(fd, "# Column 17. Spin parameter\n");
    fprintf(fd, "#\n");
    fprintf(fd, "# datavar lines are for partiview.  Ignore them if you're not partiview.\n");
    fprintf(fd, "#\n");
    fprintf(fd, "datavar 0 halo_number\n");
    fprintf(fd, "datavar 1 number_of_particles\n");
    fprintf(fd, "datavar 2 halo_mass\n");
    fprintf(fd, "datavar 3 virial_mass\n");
    fprintf(fd, "datavar 4 stellar_mass\n");
    fprintf(fd, "datavar 5 virial_radius\n");
    fprintf(fd, "datavar 6 x_velocity\n");
    fprintf(fd, "datavar 7 y_velocity\n");
    fprintf(fd, "datavar 8 z_velocity\n");
    fprintf(fd, "datavar 9 velocity_dispersion\n");
    fprintf(fd, "datavar 10 x_angular_momentum\n");
    fprintf(fd, "datavar 11 y_angular_momentum\n");
    fprintf(fd, "datavar 12 z_angular_momentum\n");
    fprintf(fd, "datavar 13 spin\n");
    fprintf(fd, "\n");

    if (HaloFinderOutputParticleList && !HaloFinderSubfind) {

      if (debug)
	fprintf(stdout, "FOF: Saving halo particle list to %s\n", particle_fname);

      file_id = H5Fcreate(particle_fname, H5F_ACC_TRUNC, H5P_DEFAULT, 
			  H5P_DEFAULT);
      group_id = H5Gcreate(file_id, "/Parameters", 0);
      writeScalarAttribute(group_id, HDF5_REAL, "Redshift", &redshift);
      writeScalarAttribute(group_id, HDF5_PREC, "Time", &EnzoTime);
      writeScalarAttribute(group_id, HDF5_INT, "Number of groups", &AllVars.NgroupsAll);
      H5Gclose(group_id);

    } // ENDIF output particle list

  } // ENDIF ROOT_PROCESSOR

  for (gr = AllVars.NgroupsAll-1; gr >= 0; gr--) {

    if (MyProcessorNumber == ROOT_PROCESSOR) {
      head = AllVars.GroupDatAll[gr].Tag;
      len = AllVars.GroupDatAll[gr].Len;
      Pbuf = new FOF_particle_data[len];
    }

    get_particles(0, head, len, Pbuf, AllVars);

    if (MyProcessorNumber == ROOT_PROCESSOR) {

      get_properties(AllVars, Pbuf, len, false, &cm[0], &cmv[0], &mtot, &mstars, &mvir, &rvir, 
		     AM, &vrms, &spin);

      if (debug && gr == AllVars.NgroupsAll-1)
	fprintf(stdout, "FOF: Largest group has %"ISYM" particles"
		" (%"GSYM" M_sun)\n", len, mtot);


      fprintf(fd, "%12"GOUTSYM" %12"GOUTSYM" %12"GOUTSYM" %12"ISYM" %12"ISYM" %12"GOUTSYM" %12"GOUTSYM" %12"GOUTSYM" %12"GOUTSYM" %12"GOUTSYM" %12"GOUTSYM" %12"GOUTSYM" %12"GOUTSYM" %12"GOUTSYM" %12"GOUTSYM" %12"GOUTSYM" %12"GOUTSYM"\n",
	      cm[0], cm[1], cm[2], AllVars.NgroupsAll-1-gr, len, 
	      mtot, mvir, mstars, rvir, cmv[0], cmv[1], cmv[2], vrms, AM[0], AM[1], AM[2], spin);

      if (HaloFinderOutputParticleList && !HaloFinderSubfind) {

	temp = new double[3*len];
	TempPINT = new PINT[len];
	index = 0;
	for (dim = 0; dim < 3; dim++)
	  for (i = 0; i < len; i++, index++)
	    temp[index] = Pbuf[i].Pos[dim] / AllVars.BoxSize;
	for (i = 0; i < len; i++)
	  TempPINT[i] = Pbuf[i].PartID;

	sprintf(halo_name, "Halo%8.8d", AllVars.NgroupsAll-1-gr);
	group_id = H5Gcreate(file_id, halo_name, 0);
	writeScalarAttribute(group_id, HDF5_REAL, "Total Mass", &mtot);
	writeScalarAttribute(group_id, HDF5_REAL, "Stellar Mass", &mstars);
	writeScalarAttribute(group_id, HDF5_REAL, "Spin parameter", &spin);
	writeScalarAttribute(group_id, HDF5_REAL, "Velocity dispersion", &vrms);
	writeArrayAttribute(group_id, HDF5_PREC, 3, "Center of mass", cm);
	writeArrayAttribute(group_id, HDF5_REAL, 3, "Mean velocity [km/s]", cmv);
	writeArrayAttribute(group_id, HDF5_REAL, 3, "Angular momentum [Mpc * km/s]", AM);

	hdims[0] = 3;
	hdims[1] = (hsize_t) len;
	dspace_id = H5Screate_simple(2, hdims, NULL);
	dset_id = H5Dcreate(group_id, "Particle Position", H5T_NATIVE_DOUBLE, dspace_id,
			    H5P_DEFAULT);
	H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) temp);
	H5Sclose(dspace_id);
	H5Dclose(dset_id);
	
	hdims[0] = (hsize_t) len;
	hdims[1] = 1;
	dspace_id = H5Screate_simple(1, hdims, NULL);
	dset_id = H5Dcreate(group_id, "Particle ID", HDF5_PINT, dspace_id,
			    H5P_DEFAULT);
	H5Dwrite(dset_id, HDF5_PINT, H5S_ALL, H5S_ALL, H5P_DEFAULT, (VOIDP) TempPINT);
	H5Sclose(dspace_id);
	H5Dclose(dset_id);

	H5Gclose(group_id);

	delete [] temp;
	delete [] TempPINT;
	
      } // ENDIF output particle list

      delete [] Pbuf;
      
    } // ENDIF ROOT_PROCESSOR

  } // ENDFOR groups

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    fclose(fd);
    
    if (HaloFinderOutputParticleList && !HaloFinderSubfind)
      H5Fclose(file_id);

  }
  
  return;

}

/************************************************************************/

int get_particles(int dest, int minid, int len, FOF_particle_data *buf,
		  FOFData &AllVars) 
{

#ifdef USE_MPI
  MPI_Status status;
#endif
  int i, imin, imax, pp, nlocal, nrecv;
  FOF_particle_data *localbuf;

  if (NumberOfProcessors == 1) {

    // No communication required.  Just created an array of particles
    // from the linked list.

    i = 0;
    pp = AllVars.Head[minid - AllVars.Noffset[MyProcessorNumber]];
    do {
      buf[i++] = AllVars.P[pp];
    } while (pp = AllVars.Next[pp]);

    return len;

  } // ENDIF serial


#ifdef USE_MPI
  MPI_Bcast(&minid,  1, IntDataType, dest, MPI_COMM_WORLD);
  MPI_Bcast(&len,    1, IntDataType, dest, MPI_COMM_WORLD);
#endif

  localbuf = new FOF_particle_data[len];
  nlocal = 0;

  if (minid >= (1 + AllVars.Noffset[MyProcessorNumber]) && 
      minid < (1 + AllVars.Noffset[MyProcessorNumber] + 
	       AllVars.Nslab[MyProcessorNumber])) {
    pp = AllVars.Head[minid - AllVars.Noffset[MyProcessorNumber]];
    do {
      if (pp <= AllVars.Nslab[MyProcessorNumber])
	localbuf[nlocal++]= AllVars.P[pp];
    } while (pp = AllVars.Next[pp]);      
  } // ENDIF
  else if (AllVars.Ncontrib) {

    /* do a bisection to speed things up */

    imin = 0; 
    imax = AllVars.Ncontrib - 1;
    if (AllVars.ContribID[imin] <= minid && 
	minid <= AllVars.ContribID[imax]) {
      i = (imin + imax) / 2;
      while (AllVars.ContribID[i] != minid) {
	i = (imin + imax) / 2;

	if (AllVars.ContribID[i] < minid)
	  imin = i;
	if (minid < AllVars.ContribID[i])
	  imax = i;

	if ((imax - imin) <= 1) {
	  if (AllVars.ContribID[imax] == minid)
	    i = imax;
	  if (AllVars.ContribID[imin] == minid)
	    i = imin;
	  break;
	} // ENDIF
      } // ENDWHILE

      if (AllVars.ContribID[i] == minid) {
	pp = AllVars.ContribHead[i];
	do {
	  if (pp <= AllVars.Nslab[MyProcessorNumber])
	    localbuf[nlocal++] = AllVars.P[pp];
	} while (pp = AllVars.Next[pp]);
      } // ENDIF

    } // ENDIF
  } // ENDELSE IF Ncontrib

#ifdef USE_MPI
  if (MyProcessorNumber == dest) {
    for (i = 0; i < nlocal; i++)
      buf[i] = localbuf[i];
    for (i = 0; i < NumberOfProcessors; i++) {
      if (i != MyProcessorNumber) {
	MPI_Recv(&nrecv, 1, IntDataType, i, i, MPI_COMM_WORLD, &status);
	if (nrecv) {
	  MPI_Recv(&buf[nlocal], nrecv*sizeof(FOF_particle_data), 
		   MPI_BYTE, i, i, MPI_COMM_WORLD, &status);
	  nlocal += nrecv;
	} // ENDIF nrecv
      } // ENDIF other processor
    } // ENDFOR processors
  }
  else {
    MPI_Ssend(&nlocal, 1, IntDataType, dest, MyProcessorNumber, MPI_COMM_WORLD);
    if (nlocal)
      MPI_Ssend(localbuf, nlocal*sizeof(FOF_particle_data), 
		MPI_BYTE, dest, MyProcessorNumber, MPI_COMM_WORLD);
  } // ENDELSE
#endif /* USE_MPI */

  delete [] localbuf;

  if (MyProcessorNumber == dest && nlocal != len) {
    ENZO_VFAIL("local= %"ISYM"  len=%"ISYM"\n", nlocal, len)
  }

  return len;
}

/************************************************************************/

int link_across(FOFData &AllVars)
{
#ifdef USE_MPI
  MPI_Status status;
#endif
  FOF_particle_data *buftoleft, *buftoright, *buffer;
  int    i, j, slab, nl, nr, nbuf, len;
  int    leftTask, rightTask;
  int    pp, newid, nlinktot;
  id_data *iddat;

  buftoleft  = new FOF_particle_data[AllVars.NtoLeft[MyProcessorNumber]];
  buftoright = new FOF_particle_data[AllVars.NtoRight[MyProcessorNumber]];
  buffer     = new FOF_particle_data[AllVars.Nshadow[MyProcessorNumber]];

  nl = nr = nbuf= 0;
  
  for (i = 1; i <= AllVars.Nslab[MyProcessorNumber]; i++) {
    //slab = (AllVars.P[i].Pos[0] / AllVars.BoxSize) * NumberOfProcessors;
    slab = AllVars.P[i].slab;

    if (AllVars.P[i].Pos[0] < 
	slab * (AllVars.BoxSize / NumberOfProcessors) + AllVars.SearchRadius)
      buftoleft[nl++] = AllVars.P[i];
      
    if (AllVars.P[i].Pos[0] > 
	(slab+1) * (AllVars.BoxSize / NumberOfProcessors) - 
	AllVars.SearchRadius)
      buftoright[nr++]= AllVars.P[i];
  } // ENDFOR particles

  rightTask = MyProcessorNumber+1;
  if (rightTask == NumberOfProcessors)
    rightTask = 0;

  leftTask = MyProcessorNumber-1;
  if (leftTask < 0)
    leftTask = NumberOfProcessors-1;

#ifdef USE_MPI
  if (MyProcessorNumber & 1) {
    MPI_Recv(&buffer[nbuf], AllVars.NtoLeft[rightTask]*sizeof(FOF_particle_data), 
	     MPI_BYTE, rightTask, rightTask, MPI_COMM_WORLD, &status);
    nbuf += AllVars.NtoLeft[rightTask];
    MPI_Ssend(buftoright, 
	      AllVars.NtoRight[MyProcessorNumber] * sizeof(FOF_particle_data), 
	      MPI_BYTE, rightTask, MyProcessorNumber, MPI_COMM_WORLD);
  } // ENDIF odd processors
  else {
    MPI_Ssend(buftoleft, 
	      AllVars.NtoLeft[MyProcessorNumber] * sizeof(FOF_particle_data), 
	      MPI_BYTE, leftTask, MyProcessorNumber, MPI_COMM_WORLD);
    MPI_Recv(&buffer[nbuf], 
	     AllVars.NtoRight[leftTask] * sizeof(FOF_particle_data),
	     MPI_BYTE, leftTask, leftTask, MPI_COMM_WORLD, &status);
    nbuf += AllVars.NtoRight[leftTask];
  } // ENDELSE

  if (MyProcessorNumber & 1) {
    MPI_Recv(&buffer[nbuf], 
	     AllVars.NtoRight[leftTask] * sizeof(FOF_particle_data), 
	     MPI_BYTE, leftTask, leftTask, MPI_COMM_WORLD, &status);
      nbuf += AllVars.NtoRight[leftTask];
      MPI_Ssend(buftoleft, 
		AllVars.NtoLeft[MyProcessorNumber] * sizeof(FOF_particle_data), 
		MPI_BYTE, leftTask, MyProcessorNumber, MPI_COMM_WORLD);
  } // ENDIF odd processors
  else {
    MPI_Ssend(buftoright, 
	      AllVars.NtoRight[MyProcessorNumber] * sizeof(FOF_particle_data), 
	      MPI_BYTE, rightTask, MyProcessorNumber, MPI_COMM_WORLD);
    MPI_Recv(&buffer[nbuf], 
	     AllVars.NtoLeft[rightTask] * sizeof(FOF_particle_data), 
	     MPI_BYTE, rightTask, rightTask, MPI_COMM_WORLD, &status);
    nbuf += AllVars.NtoLeft[rightTask];
  } // ENDELSE
#endif /* USE_MPI */

  iddat = new id_data[nbuf];

  for (i = 0; i < nbuf; i++) {
    iddat[i].ID = buffer[i].MinID;
    iddat[i].index = 1 + AllVars.Nslab[MyProcessorNumber] + i;
  }

  qsort(iddat, nbuf, sizeof(id_data), comp_func);

  AllVars.NLinkAcross = 0;

  for (i = 0; i < nbuf-1;) {
    j = i + 1;
    while (iddat[i].ID == iddat[j].ID) {
      linkit(iddat[i].index, iddat[j].index, AllVars);
      j++;
      if (j>=nbuf)
	break;
    } // ENDWHILE
    i = j;
  } // ENDFOR

  delete [] iddat;
  delete [] buffer;
  delete [] buftoright;
  delete [] buftoleft;

#ifdef USE_MPI
  MPI_Allreduce(&AllVars.NLinkAcross, &nlinktot, 1, IntDataType, MPI_SUM, 
		MPI_COMM_WORLD);
#endif

  return nlinktot;
}




void compile_group_catalogue(FOFData &AllVars)
{
#ifdef USE_MPI
  MPI_Status  status;
#endif
  int i, n, gr, tot, count;
  int nbound, Nbound;
  
  for (n = 1, AllVars.Ngroups = AllVars.Ncontrib = nbound = 0; 
       n <= AllVars.Nlocal; n++) {
    if (AllVars.Head[n] == n)
      if (AllVars.P[n].GrLen >= AllVars.GroupMinLen) {
	if (AllVars.P[n].MinID >= (1 + AllVars.Noffset[MyProcessorNumber]) && 
	    AllVars.P[n].MinID < (1 + AllVars.Noffset[MyProcessorNumber] + 
				  AllVars.Nslab[MyProcessorNumber])) {
	  AllVars.Ngroups++;  /* task hosts the global head of the group */
	  nbound += AllVars.P[n].GrLen;
	}
	else
	  AllVars.Ncontrib++;  /* task hosts a contribution to a group */
      } // ENDIF enough particles
  } // ENDFOR

  if (NumberOfProcessors == 1) {
    AllVars.NgroupsAll = AllVars.Ngroups;
    Nbound = nbound;
  }
  else {
#ifdef USE_MPI
  MPI_Allreduce(&AllVars.Ngroups, &AllVars.NgroupsAll, 1, IntDataType, MPI_SUM, 
		MPI_COMM_WORLD);
  MPI_Allreduce(&nbound, &Nbound, 1, IntDataType, MPI_SUM, MPI_COMM_WORLD);
#endif
  } // ENDELSE serial

  AllVars.GroupDat = new gr_data[AllVars.Ngroups];
  
  AllVars.ContribID =   ivector(0, AllVars.Ncontrib-1); 
  AllVars.ContribHead = ivector(0, AllVars.Ncontrib-1);

  for (n = 1, AllVars.Ngroups = AllVars.Ncontrib = 0; 
       n <= AllVars.Nlocal; n++) {
    if (AllVars.Head[n] == n)
      if (AllVars.P[n].GrLen >= AllVars.GroupMinLen) {
	if (AllVars.P[n].MinID >= (1 + AllVars.Noffset[MyProcessorNumber]) && 
	    AllVars.P[n].MinID < (1 + AllVars.Noffset[MyProcessorNumber] + 
				  AllVars.Nslab[MyProcessorNumber])) {
	  AllVars.GroupDat[AllVars.Ngroups].Len = AllVars.P[n].GrLen;
	  AllVars.GroupDat[AllVars.Ngroups].Tag = AllVars.P[n].MinID;
	  AllVars.Ngroups++; 
	}
	else {
	  AllVars.ContribID[AllVars.Ncontrib] = AllVars.P[n].MinID; 
	  AllVars.ContribHead[AllVars.Ncontrib] = AllVars.Head[n];
	  AllVars.Ncontrib++;
	}
      } // ENDIF enough particles
  } // ENDFOR

  if (AllVars.Ncontrib > 0)
    sort2_int(AllVars.Ncontrib, AllVars.ContribID-1, AllVars.ContribHead-1); 

  AllVars.NgroupsList = new int[NumberOfProcessors];

  if (NumberOfProcessors == 1)
    AllVars.NgroupsList[0] = AllVars.Ngroups;
#ifdef USE_MPI
  else
    MPI_Allgather(&AllVars.Ngroups, 1, IntDataType, AllVars.NgroupsList, 1, 
		  IntDataType, MPI_COMM_WORLD);
#endif /* USE_MPI */

  AllVars.GroupDatAll = new gr_data[AllVars.NgroupsAll];

  if (MyProcessorNumber == ROOT_PROCESSOR) {
    for (i = 0; i < AllVars.Ngroups; i++)
      AllVars.GroupDatAll[i] = AllVars.GroupDat[i];
      
#ifdef USE_MPI
    count = AllVars.Ngroups;
    for (i = 1; i < NumberOfProcessors; i++) {
      MPI_Recv(&AllVars.GroupDatAll[count], 
	       AllVars.NgroupsList[i]*sizeof(gr_data), MPI_BYTE, 
	       i, 0, MPI_COMM_WORLD, &status);
      count += AllVars.NgroupsList[i];
    }
#endif /* USE_MPI */
  } // ENDIF root
#ifdef USE_MPI
  else
    MPI_Ssend(&AllVars.GroupDat[0], AllVars.Ngroups*sizeof(gr_data), 
	      MPI_BYTE, 0, 0, MPI_COMM_WORLD);
#endif /* USE_MPI */
  
  if (debug) {
    qsort(AllVars.GroupDatAll, AllVars.NgroupsAll, sizeof(gr_data), 
	  comp_func_gr);
    if (AllVars.NgroupsAll > 0)
      fprintf(stderr, "FOF: Found %"ISYM" groups, %"ISYM" bound particles\n",
	      AllVars.NgroupsAll, Nbound);
  } // ENDIF debug

#ifdef USE_MPI
  MPI_Bcast(AllVars.GroupDatAll, AllVars.NgroupsAll*sizeof(gr_data), 
	    MPI_BYTE, 0, MPI_COMM_WORLD); 
#endif
}

/************************************************************************/

void find_minids(FOFData &AllVars)
{
  int n, pp, len, sum = 0;

  for (n = 1; n <= AllVars.Nlocal; n++)
    if (AllVars.Head[n] == n) {
      pp = n; 
      len = 0;
      do {
	if (AllVars.P[pp].ID >= (1 + AllVars.Noffset[MyProcessorNumber]) && 
	    AllVars.P[pp].ID < (1 + AllVars.Noffset[MyProcessorNumber] + 
				AllVars.Nslab[MyProcessorNumber]))
	  len++;
      } while (pp = AllVars.Next[pp]);

      AllVars.P[n].MinID = AllVars.P[n].ID;
      AllVars.P[n].GrLen = len; /* Len[n]; */

      pp = n;
      while (pp = AllVars.Next[pp])
	if (AllVars.P[n].MinID > AllVars.P[pp].ID)
	  AllVars.P[n].MinID = AllVars.P[pp].ID;
      pp = n;
      while (pp = AllVars.Next[pp]) {
	AllVars.P[pp].MinID = AllVars.P[n].MinID;
	AllVars.P[pp].GrLen = AllVars.P[n].GrLen;
      }

      sum += len;
    } // ENDIF
}





void stitch_together(FOFData &AllVars)
{
#ifdef USE_MPI
  MPI_Status status;
#endif
  FOF_particle_data *buftoleft, *buftoright, *buffer;
  int    i, slab, nl, nr, nbuf, len;
  int    leftTask, rightTask;
  int    pp;
  PINT newid;
  idmin_data *iddat;


  buftoleft = new FOF_particle_data[AllVars.NtoLeft[MyProcessorNumber]];
  buftoright = new FOF_particle_data[AllVars.NtoRight[MyProcessorNumber]];
  buffer = new FOF_particle_data[AllVars.Nshadow[MyProcessorNumber]];

  nl = nr = nbuf = 0;
  
  for (i = 1; i <= AllVars.Nslab[MyProcessorNumber]; i++) {
    //slab = (AllVars.P[i].Pos[0] / AllVars.BoxSize) * NumberOfProcessors;
    slab = AllVars.P[i].slab;
    if (AllVars.P[i].Pos[0] < slab*(AllVars.BoxSize/NumberOfProcessors) + 
	AllVars.SearchRadius)
      buftoleft[nl++] = AllVars.P[i];
      
    if (AllVars.P[i].Pos[0] > (slab+1)*(AllVars.BoxSize/NumberOfProcessors) - 
	AllVars.SearchRadius)
      buftoright[nr++] = AllVars.P[i];
  } // ENDFOR particles

  rightTask = MyProcessorNumber + 1;
  if (rightTask == NumberOfProcessors)
    rightTask = 0;

  leftTask = MyProcessorNumber - 1;
  if (leftTask < 0)
    leftTask = NumberOfProcessors - 1;

#ifdef USE_MPI
  if (MyProcessorNumber & 1) {
    MPI_Recv(&buffer[nbuf], 
	     AllVars.NtoLeft[rightTask] * sizeof(FOF_particle_data), 
	     MPI_BYTE, rightTask, rightTask, MPI_COMM_WORLD, &status);
    nbuf +=  AllVars.NtoLeft[rightTask];
    MPI_Ssend(buftoright, 
	      AllVars.NtoRight[MyProcessorNumber] * sizeof(FOF_particle_data), 
	      MPI_BYTE, rightTask, MyProcessorNumber, MPI_COMM_WORLD);
  }
  else {
    MPI_Ssend(buftoleft, 
	      AllVars.NtoLeft[MyProcessorNumber] * sizeof(FOF_particle_data), 
	      MPI_BYTE, leftTask, MyProcessorNumber, MPI_COMM_WORLD);
    MPI_Recv(&buffer[nbuf], 
	     AllVars.NtoRight[leftTask] * sizeof(FOF_particle_data), 
	     MPI_BYTE, leftTask, leftTask, MPI_COMM_WORLD, &status);
    nbuf +=  AllVars.NtoRight[leftTask];
  } // ENDELSE

  if (MyProcessorNumber & 1) {
    MPI_Recv(&buffer[nbuf], 
	     AllVars.NtoRight[leftTask] * sizeof(FOF_particle_data), 
	     MPI_BYTE, leftTask, leftTask, MPI_COMM_WORLD, &status);
      nbuf += AllVars.NtoRight[leftTask];
      MPI_Ssend(buftoleft, 
		AllVars.NtoLeft[MyProcessorNumber] * sizeof(FOF_particle_data), 
		MPI_BYTE, leftTask, MyProcessorNumber, MPI_COMM_WORLD);
  } // ENDIF odd processor
  else {
    MPI_Ssend(buftoright, 
	      AllVars.NtoRight[MyProcessorNumber] * sizeof(FOF_particle_data), 
	      MPI_BYTE, rightTask, MyProcessorNumber, MPI_COMM_WORLD);
    MPI_Recv(&buffer[nbuf], 
	     AllVars.NtoLeft[rightTask] * sizeof(FOF_particle_data), 
	     MPI_BYTE, rightTask, rightTask, MPI_COMM_WORLD, &status);
    nbuf += AllVars.NtoLeft[rightTask];
    }
#endif /* USE_MPI */
  
  iddat = new idmin_data[nbuf];

  for (i = 0; i < nbuf; i++) {
    iddat[i].minID = buffer[i].MinID;
    iddat[i].index = 1 + AllVars.Nslab[MyProcessorNumber] + i;
    iddat[i].len   = buffer[i].GrLen;
  }

  qsort(iddat, nbuf, sizeof(idmin_data), comp_func2);

  for (i = 0; i < nbuf; i++)
    if (iddat[i].minID != AllVars.P[iddat[i].index].MinID || iddat[i].len>0) {
      newid = AllVars.P[iddat[i].index].MinID;
      if (iddat[i].minID < newid)
	newid = iddat[i].minID;

      len = AllVars.P[iddat[i].index].GrLen + iddat[i].len;
      pp = AllVars.Head[iddat[i].index];
      do {
	AllVars.P[pp].MinID = newid;
	AllVars.P[pp].GrLen = len;
      } while (pp = AllVars.Next[pp]);

      while (i < nbuf-1) {
	if (iddat[i+1].minID == iddat[i].minID)
	  i++;
	else
	  break;
      }
    } // ENDIF

  delete [] iddat;
  delete [] buffer;
  delete [] buftoright;
  delete [] buftoleft;
}

/************************************************************************/

void exchange_shadow(FOFData &AllVars, int TopGridResolution, bool SmoothData)
{
#ifdef USE_MPI
  MPI_Status status;
#endif
  FOF_particle_data *buftoleft, *buftoright;
  int    i, slab, nl, nr, region, dim;
  int    leftTask, rightTask;
  bool inside;
  double sr;
  float StaticRegionCellWidth[MAX_STATIC_REGIONS+1];

//  if (debug)
//    fprintf(stdout, "FOF: exchanging shadows ...\n");

  buftoleft =  new FOF_particle_data[AllVars.NtoLeft[MyProcessorNumber]];
  buftoright = new FOF_particle_data[AllVars.NtoRight[MyProcessorNumber]];

  // Pre-compute cell widths for each static region
  if (SmoothData) {
    StaticRegionCellWidth[0] = AllVars.BoxSize / TopGridResolution;
    for (i = 0; i < MAX_STATIC_REGIONS; i++) {
      if (StaticRefineRegionRightEdge[i][0] > 0)
	StaticRegionCellWidth[i+1] = AllVars.BoxSize / TopGridResolution /
	  pow(RefineBy, StaticRefineRegionLevel[i]+1);
      else
	StaticRegionCellWidth[i+1] = 0;
    }
  } // ENDIF SmoothData
    
  nl = nr = 0;
  
  for (i = 1; i <= AllVars.Nlocal; i++) {
    //slab = (AllVars.P[i].Pos[0] / AllVars.BoxSize) * NumberOfProcessors;
    slab = AllVars.P[i].slab;

    if (slab != MyProcessorNumber)
      ENZO_FAIL("FOF: particle not on the right processor (slab)!");

    if (SmoothData) {

      sr = AllVars.LinkLength * StaticRegionCellWidth[0];
      for (region = MAX_STATIC_REGIONS-1; region >= 0; region--) {
	if (StaticRefineRegionLevel[region] < 0)
	  continue;
	inside = true;
	for (dim = 0; dim < 3; dim++) {
	  inside &= 
	    AllVars.P[i].Pos[dim] >= StaticRefineRegionLeftEdge[region][dim] &&
	    AllVars.P[i].Pos[dim] <= StaticRefineRegionRightEdge[region][dim];
	  if (!inside)
	    break;
	} // ENDFOR dim

	if (inside) {
	  sr = AllVars.LinkLength * StaticRegionCellWidth[region];
	  break;
	}

      } // ENDFOR region

    } // ENDIF SmoothData
    else
      sr = AllVars.SearchRadius;
		  
    if (AllVars.P[i].Pos[0] < slab*(AllVars.BoxSize/NumberOfProcessors) + sr)
      buftoleft[nl++] = AllVars.P[i];
      
    if (AllVars.P[i].Pos[0] > (slab+1)*(AllVars.BoxSize/NumberOfProcessors) - sr)
      buftoright[nr++] = AllVars.P[i];

  } // ENDFOR particles

  if (nl != AllVars.NtoLeft[MyProcessorNumber]) {
    ENZO_VFAIL("[proc %"ISYM"] error: shadows don't match! "
	    "nl = %"ISYM", NtoLeft[%"ISYM"] = %"ISYM"\n", 
	    MyProcessorNumber, nl, MyProcessorNumber, 
	    AllVars.NtoLeft[MyProcessorNumber])
  }
  
  if (nr != AllVars.NtoRight[MyProcessorNumber]) {
    ENZO_VFAIL("[proc %"ISYM"] error: shadows don't match! "
	    "nr = %"ISYM", NtoRight[%"ISYM"] = %"ISYM"\n",
	    MyProcessorNumber, nr, MyProcessorNumber, 
	    AllVars.NtoRight[MyProcessorNumber])
  }

  rightTask = MyProcessorNumber + 1;
  if (rightTask == NumberOfProcessors)
    rightTask = 0;

  leftTask = MyProcessorNumber - 1;
  if (leftTask<0)
    leftTask = NumberOfProcessors - 1;

#ifdef USE_MPI
  if (MyProcessorNumber & 1) {
    MPI_Recv(&AllVars.P[1+AllVars.Nlocal], 
	     AllVars.NtoLeft[rightTask] * sizeof(FOF_particle_data), 
	     MPI_BYTE, rightTask, rightTask, MPI_COMM_WORLD, &status);
    AllVars.Nlocal += AllVars.NtoLeft[rightTask];
    MPI_Ssend(buftoright, 
	      AllVars.NtoRight[MyProcessorNumber] * sizeof(FOF_particle_data), 
	      MPI_BYTE, rightTask, MyProcessorNumber, MPI_COMM_WORLD);
  }
  else {
    MPI_Ssend(buftoleft, 
	      AllVars.NtoLeft[MyProcessorNumber] * sizeof(FOF_particle_data), 
	      MPI_BYTE, leftTask, MyProcessorNumber, MPI_COMM_WORLD);
    MPI_Recv(&AllVars.P[1+AllVars.Nlocal], 
	     AllVars.NtoRight[leftTask] * sizeof(FOF_particle_data), 
	     MPI_BYTE, leftTask, leftTask, MPI_COMM_WORLD, &status);
    AllVars.Nlocal += AllVars.NtoRight[leftTask];
  } // ENDELSE

  if (MyProcessorNumber & 1) {
    MPI_Recv(&AllVars.P[1+AllVars.Nlocal], 
	     AllVars.NtoRight[leftTask] * sizeof(FOF_particle_data), 
	     MPI_BYTE, leftTask, leftTask, MPI_COMM_WORLD, &status);
    AllVars.Nlocal += AllVars.NtoRight[leftTask];
    MPI_Ssend(buftoleft, 
	      AllVars.NtoLeft[MyProcessorNumber] * sizeof(FOF_particle_data), 
	      MPI_BYTE, leftTask, MyProcessorNumber, MPI_COMM_WORLD);
  }
  else {
    MPI_Ssend(buftoright, 
	      AllVars.NtoRight[MyProcessorNumber] * sizeof(FOF_particle_data), 
	      MPI_BYTE, rightTask, MyProcessorNumber, MPI_COMM_WORLD);
    MPI_Recv(&AllVars.P[1+AllVars.Nlocal], 
	     AllVars.NtoLeft[rightTask]*sizeof(FOF_particle_data), 
	     MPI_BYTE, rightTask, rightTask, MPI_COMM_WORLD, &status);
    AllVars.Nlocal += AllVars.NtoLeft[rightTask];
  } // ENDELSE
#endif /* USE_MPI */

  delete [] buftoright;
  delete [] buftoleft;
}

/************************************************************************/

void link_local_slab(FOFData &AllVars)
{
  int  nx,ny,nz;
  int  iter;
  int  count;

  iter = 1;

  for (AllVars.GridCorner[0] = 0, nx = 0; nx < AllVars.Nx; 
       AllVars.GridCorner[0] += (AllVars.Grid - 2.0) / AllVars.Grid * 
	 AllVars.GridExtension, nx++)
    for (AllVars.GridCorner[1] = 0, ny = 0; ny < AllVars.Ny; 
	 AllVars.GridCorner[1] += (AllVars.Grid - 2.0) / AllVars.Grid * 
	   AllVars.GridExtension, ny++)
      for (AllVars.GridCorner[2] = 0, nz = 0; nz < AllVars.Nz; 
	   AllVars.GridCorner[2] += (AllVars.Grid - 2.0) / AllVars.Grid * 
	     AllVars.GridExtension, nz++) {
//	if (debug)
//	  printf("Grid placement number: %"ISYM" out of %"ISYM"\n", iter++, 
//		 AllVars.Nx * AllVars.Ny * AllVars.Nz);
	  
	count = coarse_binning(AllVars);
	  
	if (count)
	  find_groups(AllVars);
      }
  
//  if (debug)
//    printf("main linking done.\n"); 

}

/************************************************************************/

void init_coarse_grid(FOFData &AllVars)
{
  int i;

  AllVars.GridExtension = AllVars.Grid * AllVars.SearchRadius;
  
  if (AllVars.GridExtension > AllVars.BoxSize) {
    AllVars.Grid = AllVars.BoxSize / AllVars.SearchRadius;
    AllVars.GridExtension = AllVars.Grid * AllVars.SearchRadius;
  }

  if (AllVars.BoxSize / ( (AllVars.Grid-2.0)/AllVars.Grid * 
			  AllVars.GridExtension ) > AllVars.MaxPlacement)
    AllVars.GridExtension = AllVars.BoxSize / AllVars.MaxPlacement *
      AllVars.Grid / (AllVars.Grid-2.001);
 
  AllVars.Nx = (int) (AllVars.BoxSize/( (AllVars.Grid-2.0) / AllVars.Grid * 
					AllVars.GridExtension) + 1);
  AllVars.Ny = (int) (AllVars.BoxSize/( (AllVars.Grid-2.0) / AllVars.Grid * 
					AllVars.GridExtension) + 1);
  AllVars.Nz = (int) (AllVars.BoxSize/( (AllVars.Grid-2.0) / AllVars.Grid * 
					AllVars.GridExtension) + 1);

//  if (debug)
//    printf("\nGrid has to be placed (%"ISYM"|%"ISYM"|%"ISYM") times in each dimension.\n", 
//	   AllVars.Nx, AllVars.Ny, AllVars.Nz);

  AllVars.GridFirst = i3tensor(0, AllVars.Grid-1, 0, AllVars.Grid-1, 0, 
			       AllVars.Grid-1);
  AllVars.GridLast  = i3tensor(0, AllVars.Grid-1, 0, AllVars.Grid-1, 0, 
			       AllVars.Grid-1);
  AllVars.GridFlag  = i3tensor(0, AllVars.Grid-1, 0, AllVars.Grid-1, 0, 
			       AllVars.Grid-1);

  AllVars.GridNext = ivector(1, AllVars.Nlocal);

  AllVars.Tail = ivector(1, AllVars.Nlocal);
  AllVars.Len  = ivector(1, AllVars.Nlocal);
  AllVars.Head = ivector(1, AllVars.Nlocal);
  AllVars.Next = ivector(1, AllVars.Nlocal);

//  if (debug)
//    printf("Nlocal = %"ISYM" Task = %"ISYM"\n", AllVars.Nlocal, MyProcessorNumber);

  for (i = 1; i <= AllVars.Nlocal; i++) {
    AllVars.Head[i] = i;
    AllVars.Tail[i] = i;
    AllVars.Next[i] = 0;
    AllVars.Len[i]  = 1;
  }
}

/************************************************************************/

/* New code added to ensure no two particles lie on top of each other. */

void marking(FOFData &AllVars)
{
  float posold[3];
  int   i,k,iter,idone;

//  if (debug)
//    fprintf(stdout, "marking ...\n");

  iter = 0;
  do {
    qsort(&AllVars.P[1], AllVars.Nlocal, sizeof(FOF_particle_data), 
	  comp_func_partcoord);

    for (i = 2, idone = 0; i <= AllVars.Nlocal; i++)
      if (fabs(AllVars.P[i-1].Pos[0] - AllVars.P[i].Pos[0]) < 1e-3*AllVars.Epsilon &&
	  fabs(AllVars.P[i-1].Pos[1] - AllVars.P[i].Pos[1]) < 1e-3*AllVars.Epsilon &&
	  fabs(AllVars.P[i-1].Pos[2] - AllVars.P[i].Pos[2]) < 1e-3*AllVars.Epsilon) {

	for (k = 0; k < 3; k++) {
	  posold[k]= AllVars.P[i].Pos[k];
	  AllVars.P[i].Pos[k] += (0.001*AllVars.Epsilon)*(2*drand48()-1);
	}
	idone++;
      } // ENDIF
      iter++;
  } while (idone > 0 && iter < 10);

  qsort(&AllVars.P[1], AllVars.Nlocal, sizeof(FOF_particle_data), 
	comp_func_partcoord);

  for (i = 1; i <= AllVars.Nlocal; i++)
    AllVars.P[i].ID = AllVars.Noffset[MyProcessorNumber] + i; 
} // END marking()

/************************************************************************/

double FOF_periodic(double x, double boxsize)
{
  if (x > 0.5*boxsize)
    x -= boxsize;

  if (x < -0.5*boxsize)
    x += boxsize;
  
  return x;
}

double FOF_periodic_wrap(double x, double boxsize)
{
  while (x > boxsize)
    x -= boxsize;

  while (x < 0)
    x += boxsize;
  
  return x;
}

/************************************************************************/

int coarse_binning(FOFData &AllVars)
{
  int i,j,k,n,count;
  double fac;
  double pos[3];

//  if (debug)
//    fprintf(stdout, "coarse binning...");

  for (i = 0; i < AllVars.Grid; i++)
    for (j = 0; j < AllVars.Grid; j++)
      for (k = 0; k < AllVars.Grid; k++) {
	AllVars.GridFirst[i][j][k] = 0;
	AllVars.GridFlag[i][j][k] = 0;
      }
  
  for (n = 1; n <= AllVars.Nlocal; n++)
    AllVars.GridNext[n] = 0;
 
  fac = AllVars.Grid / AllVars.GridExtension;

  for (n = 1, count = 0; n <= AllVars.Nlocal; n++) {
    for (k = 0; k < 3; k++) {
      pos[k] = AllVars.P[n].Pos[k];
      if (pos[k] < AllVars.GridCorner[k])
	pos[k] += AllVars.BoxSize;
      else if (pos[k] >= (AllVars.GridCorner[k] + AllVars.GridExtension))
	pos[k] -= AllVars.BoxSize;
    } // ENDIF k

    if (pos[0] >= AllVars.GridCorner[0])
      if (pos[0] < (AllVars.GridCorner[0] + AllVars.GridExtension))
	if (pos[1] >= AllVars.GridCorner[1])
	  if (pos[1] < (AllVars.GridCorner[1] + AllVars.GridExtension))
	    if (pos[2] >= AllVars.GridCorner[2])
	      if (pos[2] < (AllVars.GridCorner[2] + AllVars.GridExtension)) {
		i = (pos[0]-AllVars.GridCorner[0])*fac;
		j = (pos[1]-AllVars.GridCorner[1])*fac;
		k = (pos[2]-AllVars.GridCorner[2])*fac;
		if (AllVars.GridFirst[i][j][k]) {
		  AllVars.GridNext[AllVars.GridLast[i][j][k]] = n;
		  AllVars.GridLast[i][j][k] = n;
		}
		else {
		  AllVars.GridFirst[i][j][k] = n;
		  AllVars.GridLast[i][j][k] = n;
		}
		count++;
	      } // ENDIF pos[2]
  } // ENDFOR n

  if (debug1)
    printf("done.  (count=%"ISYM")\n",count);

  return count;
}

/************************************************************************/

void find_groups(FOFData &AllVars)
{
  int i,j,k;
  int p;

//  if (debug)
//    printf("linking..."); fflush(stdout);

  for (i = AllVars.Grid-2; i >= 1; i--)
    for (j = AllVars.Grid-2; j >= 1; j--)
      for (k = AllVars.Grid-2; k >= 1; k--) {
	if (p = AllVars.GridFirst[i][j][k]) {
	  do {
	    check_cell(AllVars, p, i+1,j  ,k  );
	    check_cell(AllVars, p, i+1,j+1,k  );
	    check_cell(AllVars, p, i+1,j  ,k+1);
	    check_cell(AllVars, p, i+1,j+1,k+1);
	    check_cell(AllVars, p, i  ,j+1,k  );
	    check_cell(AllVars, p, i  ,j+1,k+1);
	    check_cell(AllVars, p, i  ,j  ,k+1);
	    check_cell(AllVars, p, i  ,j  ,k);    
	  } while (p = AllVars.GridNext[p]);
	} // ENDIF
      } // ENDFOR k

  for (i = AllVars.Grid-2; i >= 1; i--)
    for (j = AllVars.Grid-2; j >= 1; j--)
      for (k = AllVars.Grid-2; k >= 1; k--) {
	if (p = AllVars.GridFirst[i][j][k]) {
	  do {
	    check_cell(AllVars, p, i+1,j  ,k-1);
	    check_cell(AllVars, p, i+1,j-1,k  );
	    check_cell(AllVars, p, i  ,j-1,k+1);
	    check_cell(AllVars, p, i-1,j+1,k+1);
	    check_cell(AllVars, p, i-1,j-1,k+1);
	    check_cell(AllVars, p, i+1,j-1,k+1);
	  } while (p = AllVars.GridNext[p]);
	} // ENDIF
      } // ENDFOR k

}

/************************************************************************/

void check_cell(FOFData &AllVars, int p, int i, int j, int k)
{
  double r2,s2,dx,dy,dz;
  int pp,ss;
  int s, flag;
  
  if (s = AllVars.GridFirst[i][j][k]) {
    flag = AllVars.Head[s];
    if (AllVars.GridFlag[i][j][k])
      if (AllVars.Head[s] == AllVars.Head[p])
	return; 
  } else
    flag = 0;

  s2 = AllVars.SearchRadius * AllVars.SearchRadius;

  while (s) {
    if (AllVars.Head[s] != flag)
      flag = 0;

    /* only if not yet linked */
    if (AllVars.Head[p] != AllVars.Head[s]) {
      
      dx = AllVars.P[p].Pos[0] - AllVars.P[s].Pos[0];
      dy = AllVars.P[p].Pos[1] - AllVars.P[s].Pos[1];
      dz = AllVars.P[p].Pos[2] - AllVars.P[s].Pos[2];
      dx = FOF_periodic(dx, AllVars.BoxSize);
      dy = FOF_periodic(dy, AllVars.BoxSize);
      dz = FOF_periodic(dz, AllVars.BoxSize);

      r2 = dx*dx + dy*dy + dz*dz;

      /* ok, we have a partner */
      if (r2 < s2) {

	/* p group is longer */

	if (AllVars.Len[AllVars.Head[p]] > AllVars.Len[AllVars.Head[s]]) {

	  AllVars.Next[ AllVars.Tail[AllVars.Head[p]] ] = AllVars.Head[s];
	  AllVars.Tail[AllVars.Head[p]] = AllVars.Tail[AllVars.Head[s]];
	  AllVars.Len[AllVars.Head[p]] += AllVars.Len[AllVars.Head[s]];
	  ss = AllVars.Head[s];
		
	  do {
	    AllVars.Head[ss]=AllVars.Head[p];
	  } while (ss = AllVars.Next[ss]);

	  flag=0;
	} // ENDIF
	else {
	  AllVars.Next[ AllVars.Tail[AllVars.Head[s]] ] = AllVars.Head[p];
	  AllVars.Tail[AllVars.Head[s]] = AllVars.Tail[AllVars.Head[p]];
	  AllVars.Len[AllVars.Head[s]] += AllVars.Len[AllVars.Head[p]];
	  pp = AllVars.Head[p];
	  do {
	    AllVars.Head[pp] = AllVars.Head[s];
	  } while (pp = AllVars.Next[pp]);

	  flag=0;
	} // ENDELSE
	      
	if (AllVars.GridFlag[i][j][k])
	  return; 
      } // ENDIF p group is longer
    } // ENDIF not linked

    s = AllVars.GridNext[s];
  } // ENDWHILE

  if (flag)
    AllVars.GridFlag[i][j][k] = 1;
}

/************************************************************************/


void linkit(int p, int s, FOFData &AllVars)
{
  int ss, pp;

  /* only if not yet linked */
  if (AllVars.Head[p]!=AllVars.Head[s]) {
    AllVars.NLinkAcross++;

    /* p group is longer */

    if (AllVars.Len[AllVars.Head[p]] > AllVars.Len[AllVars.Head[s]]) {
      AllVars.Next[ AllVars.Tail[AllVars.Head[p]] ] = AllVars.Head[s];
      AllVars.Tail[AllVars.Head[p]] = AllVars.Tail[AllVars.Head[s]];
      AllVars.Len[AllVars.Head[p]] += AllVars.Len[AllVars.Head[s]];
      ss = AllVars.Head[s];

      do {
	AllVars.Head[ss] = AllVars.Head[p];
      } while (ss = AllVars.Next[ss]);
    } // ENDIF p group is longer
    else {
      AllVars.Next[ AllVars.Tail[AllVars.Head[s]] ] = AllVars.Head[p];
      AllVars.Tail[AllVars.Head[s]] = AllVars.Tail[AllVars.Head[p]];
      AllVars.Len[AllVars.Head[s]] += AllVars.Len[AllVars.Head[p]];
      pp = AllVars.Head[p];

      do {
	AllVars.Head[pp] = AllVars.Head[s];
      } while (pp = AllVars.Next[pp]);
    } // ENDELSE
  } // ENDIF not linked
}

/************************************************************************/

void *mymalloc(size_t size)
{
  void *p;

  p = malloc(size);
  
  if (!p) {
    ENZO_VFAIL("Failed to alloc %"ISYM" bytes on process %"ISYM".\n", 
	    (int) size, MyProcessorNumber)
  }
  return p;
}

/************************************************************************/

#ifdef UNUSED
void write_ascii_catalog(char *catalogue_fname, char *fofprop_fname, 
			 char *cataloguetxt)
{

  FILE *in, *out;
  int i, Ngroups, *Npart_g;
  float *CM, *Mass;
  double SolarMass = UnitMass_in_g / 1.989e33;
  float RootBoxSize[3];

  if (MyProcessorNumber == 0) {


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
#endif /* UNUSED */
