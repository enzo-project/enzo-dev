/* This is a wrapper that runs HOP (Daniel Eisenstein and Wayne Hu astro-ph/9712200)
on enzo datasets. It reads in the particle data and hands it off to HOP.
*/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <hdf5.h>
#include <assert.h>
#include <math.h>
#define hssize_t hsize_t
//#include "extern_hdf5.h"  /*  hdf 5 prototypes */
#include "macros.h"
#include "kd.h"

/* return calls */
#define SUCCESS 1
#define FAILURE 0


int total_number_particles = 0, total_number_grids;
int total_number_dm_particles = 0;
int debug=0,nostars=1,nopacked=1;

float OmegaMatterNow = 0.3; 
float Hubble = 0.65; /* in km/s/Mpc */
float boxsize = 50.0; /* in Mpc */
float HopDensityThreshold = 80.0;

int CosmologyParameters(char *filename, float *Hubble, float *OmegaMatterNow, float *boxsize);
int NumberOfGrids(char *hierfilename);
int GetGridInfo(int numberofgrids,char *hierfilename);
int AddParticlesToArrays(int gridnum);

int InterpretCommandLine(int argc, char *argv[], char *myname, int *debug, 
    int *nostars, int *nopacked, char *inputfilename[], float *HopDensityThreshold);

float *xpositions, *ypositions, *zpositions, *particlemasses,
  *particlecrtime,*particledyntime,*particlemetalfrac;
float *dmxpositions, *dmypositions, *dmzpositions, *dmparticlemasses,
  *dmparticlecrtime,*dmparticledyntime,*dmparticlemetalfrac;
int *gridnump, array_position=0, particle_count=0;
char **gridfilenames;
long *particleID;
long *dmparticleID;

int *griddx;

float *gridlex,*gridrex;

/* Hop related stuff */

/* extern "C" void hop_main(KD kd);
extern "C" void regroup_main(float dens_outer); */

void hop_main(KD kd);
void regroup_main(float dens_outer);
void kdInit(KD &kd,int nBucket);

int main(int argc, char *argv[]){


  char *myname= argv[0];
  int i,j;
  FILE *outputfile;
  



  //char *inputfilename=argv[1];
  char *inputfilename = NULL;
  InterpretCommandLine(argc, argv, myname,
  						&debug, &nostars, &nopacked, &inputfilename, &HopDensityThreshold);
  //fprintf(stdout,"HopDensityThreshold = %f\n",HopDensityThreshold);
  
  /* get the cosmology stuff */
  int inlen = (int) strlen(inputfilename);
  char *file = new char[inlen+13];
  strcpy(file,inputfilename);
  CosmologyParameters(file, &Hubble, &OmegaMatterNow, &boxsize);
  if (debug) {
    fprintf(stdout, "Hubble %f OmegaMatterNow %f boxsize %f\n",
    Hubble, OmegaMatterNow, boxsize);
  }

  float MassConv = OmegaMatterNow / Hubble * pow(boxsize, 3) * 2.76e11;
  float DensityFactor = 1.0/(2.78e11*OmegaMatterNow*pow(Hubble, 2) *
		pow(boxsize/Hubble, 3));

  /* get hierarchy file name */
  char *hierfile = new char[inlen+13];
  strcpy(hierfile,inputfilename);
  hierfile=strcat(hierfile,".hierarchy");

  int total_number_grids = NumberOfGrids(hierfile);

  printf("there are %i total particles in %i grids\n",total_number_particles,
  	total_number_grids);

  fprintf(stderr,"**starting**\n");

  xpositions = new float[total_number_particles];
  ypositions = new float[total_number_particles];
  zpositions = new float[total_number_particles];
  particlemasses = new float[total_number_particles];
  particlecrtime = new float[total_number_particles];
  particledyntime = new float[total_number_particles];
  particlemetalfrac = new float[total_number_particles];
  particleID = new long[total_number_particles];
  /* dm only arrays */
  if (nostars) {
  fprintf(stderr,"**dm only: %d particles**\n",total_number_dm_particles);
  dmxpositions = new float[total_number_dm_particles];
  dmypositions = new float[total_number_dm_particles];
  dmzpositions = new float[total_number_dm_particles];
  dmparticlemasses = new float[total_number_dm_particles];
  dmparticlecrtime = new float[total_number_dm_particles];
  dmparticledyntime = new float[total_number_dm_particles];
  dmparticlemetalfrac = new float[total_number_dm_particles];
  dmparticleID = new long[total_number_dm_particles];
  }

  fprintf(stderr,"**ending**\n");

  for(i=0; i<total_number_particles; i++){

    /* fprintf(stderr,"%i\n",i); */

    xpositions[i] = ypositions[i] = zpositions[i] = 
      particlemasses[i] = particlecrtime[i] = 
      particledyntime[i] = particlemetalfrac[i] = -1.0;
  }

  if (nostars) {
  for(i=0; i<total_number_dm_particles; i++) {
    dmxpositions[i] = dmypositions[i] = dmzpositions[i] = 
      dmparticlemasses[i] = dmparticlecrtime[i] = 
      dmparticledyntime[i] = dmparticlemetalfrac[i] = -1.0;
      }
   }

  fprintf(stderr,"**done setting values**\n");


  /* grid file name info */
  gridfilenames = new char*[total_number_grids];

  /*  number of particles info */
  gridnump = new int[total_number_grids];
  
  /* grid info */
  griddx = new int[total_number_grids];
  gridlex = new float[total_number_grids];
  gridrex = new float[total_number_grids];


  fprintf(stderr,"**(1)**\n");


  for(i=0; i<total_number_grids; i++){
    gridnump[i]=griddx[i]=-1;
    gridlex[i]=gridrex[i]=-1.0;
    gridfilenames[i] = new char[MAX_LINE_LENGTH];
   }

  GetGridInfo(total_number_grids,hierfile);

  for(i=0; i<total_number_grids; i++)
    if(gridnump[i]>0)
      AddParticlesToArrays(i);


   for(i=0; i<total_number_particles; i++)
    if(xpositions[i] < 0.0 || xpositions[i] > 1.0 ||
       ypositions[i] < 0.0 || ypositions[i] > 1.0 ||
       zpositions[i] < 0.0 || zpositions[i] > 1.0 ||
       particleID[i] < 0 || particlemasses[i] < 0)
      fprintf(stderr,"uh-oh %f %f %f %d %f\n",xpositions[i],ypositions[i],zpositions[i],particleID[i],
      particlemasses[i]);

	/* copy values into dm-only lists */
	if (nostars) {
		j = 0;
		for (i=0; i<total_number_particles; i++) {
			if(particlemetalfrac[i]<0) {
				dmxpositions[j] = xpositions[i]; 
				dmypositions[j] = ypositions[i];
				dmzpositions[j] = zpositions[i];
				dmparticlemasses[j] = particlemasses[i];
				dmparticlecrtime[j] = particlecrtime[i];
				dmparticledyntime[j] = particledyntime[i];
				dmparticlemetalfrac[j] = particlemetalfrac[i];
				dmparticleID[j] = particleID[i];
				j+=1;
			}
		}
	}

	if (nostars) {
	  delete [] xpositions;
	  delete [] ypositions;
	  delete [] zpositions;

	  delete [] particlemasses;
	  delete [] particlecrtime;
	  delete [] particledyntime;
	  delete [] particlemetalfrac;
	  delete [] particleID;
	}
			

	/* Normalize the masses */
	
	float totalmass=0.0;
	if (nostars) {
		for(i=0; i<total_number_dm_particles; i++) {
	   		totalmass += dmparticlemasses[i];
	   	}
	} else {
		for (i=0; i<total_number_particles; i++) {
			totalmass += particlemasses[i];
		}
	}

   /* outputfile = fopen("dmposindinfo.dat","w");

  fprintf(stderr,"writing file dmposindinfo.dat\n") ;

  for(i=0; i<total_number_particles; i++)
    fprintf(outputfile,"%f\t%f\t%f\t%d\n",xpositions[i],ypositions[i],zpositions[i],particleID[i]);
    
  fclose(outputfile); */


  /* initialize the kd hop structure */

  KD kd;
  int nBucket = 16, kdcount = 0;
  kdInit(&kd, nBucket);
  if (nostars) {
    kd->nActive = total_number_dm_particles;
  } else {
    kd->nActive = total_number_particles;
  }
  kd->p = new PARTICLE[kd->nActive];
  if (kd->p == NULL) {
    fprintf(stderr, "failed allocating particles.\n");
    assert(kd->p != NULL);
  }
  
  /* printf("getparts SIZE_OF_FLOAT = %d\n",SIZE_OF_FLOAT); */
  
 	/* Copy positions into kd structure. */

   if (nostars) {
	for (i = 0; i < total_number_dm_particles; i++) {
	  kd->p[i].r[0] = dmxpositions[i];
	  kd->p[i].r[1] = dmypositions[i];
	  kd->p[i].r[2] = dmzpositions[i];
	  kd->p[i].fMass = dmparticlemasses[i]/totalmass; /* fix the mass so hop is happy*/
	  kd->p[i].iID = dmparticleID[i]; /* S. Skory */
	}
   } else {
	for (i = 0; i < total_number_particles; i++) {
	  kd->p[i].r[0] = xpositions[i];
	  kd->p[i].r[1] = ypositions[i];
	  kd->p[i].r[2] = zpositions[i];
	  kd->p[i].fMass = particlemasses[i]/totalmass; /* fix the mass so hop is happy*/
	  kd->p[i].iID = particleID[i]; /* S. Skory */
	}
   }

  /* --------------------------------------------------------------- */
  /* Call hop. */
 
  fprintf(stderr, "Calling hop...\n");
  hop_main(kd);
 
  fprintf(stderr, "Calling regroup...\n");
  regroup_main(HopDensityThreshold);
  
    /* --------------------------------------------------------------- */
  /* Read the group membership and compute group properties. */
 
  FILE *fptr;
  if ((fptr = fopen("zregroup.tag", "r")) == NULL) {
    fprintf(stderr, "Error opening regroup output zregroup.tag\n"); /* S Skory */
  }
 
  int nActive, nGroups;
  fread(&nActive, SIZE_OF_INT, 1, fptr);
  fread(&nGroups, SIZE_OF_INT, 1, fptr);
  printf("nActive = %"ISYM"(=%"ISYM")   nGroups = %"ISYM"\n", nActive, kd->nActive, nGroups);
 
  /* Allocate space and read group memberships for the particles. */
 
  int *GroupID = new int[nActive];
  if (fread(GroupID, SIZE_OF_INT, nActive, fptr) != nActive) {
    fprintf(stderr, "Error reading GroupID file zregroup.tag\n"); /* S Skory */
  }
   int *realID = new int[nActive];
  if (fread(realID, SIZE_OF_INT, nActive, fptr) != nActive) {
    fprintf(stderr, "Error reading GroupID file zregroup.tag\n"); /* S Skory */
  }
  fclose(fptr);
 
  /* Allocate space and read group memberships for the particles. */
 
  float *Density = new float[nActive];
  if ((fptr = fopen("output_hop.den", "r")) == NULL) {
    fprintf(stderr, "Error opening regroup output output_hop.den\n");
  }
  fread(&nActive, SIZE_OF_INT, 1, fptr);
  if (fread(Density, SIZE_OF_FLOAT, nActive, fptr) != nActive) {
    fprintf(stderr, "Error reading GroupID file output_hop.den\n");
  }
  fclose(fptr);
 
  /* Allocate and initialize group information. */
 
  int const NumberOfGroupProperties = 14; /* S Skory */
  float *GroupProperties[NumberOfGroupProperties];
  for (i = 0; i < NumberOfGroupProperties; i++) {
    GroupProperties[i] = new float[nGroups];
    for (j = 0; j < nGroups; j++)
      GroupProperties[i][j] = 0;
  }
 
  for (j = 0; j < nGroups; j++) { /* S Skory */
    GroupProperties[8][j] = 1.0; /* S Skory xmin -> xmax is i=9*/
    GroupProperties[10][j] = 1.0; /* S Skory ymin -> ymax is i=11*/
    GroupProperties[12][j] = 1.0; /* S Skory zmin -> zmax is i=13*/
  }

  /* Loop over particles, adding to group properties. */
  FILE *blahfp = fopen("Part-Final.txt","w"); /* S Skory */
  float Luminosity;
  for (i = 0; i < nActive; i++)
    if ((j = GroupID[i]) >= 0) {

    if (nostars) {
     fprintf(blahfp, "%"ISYM" %"ISYM" %"FSYM" %"FSYM" %"FSYM"\n", realID[i],j, dmxpositions[i],
	      dmypositions[i],dmzpositions[i]); 
	} else {
	 fprintf(blahfp, "%"ISYM" %"ISYM" %"FSYM" %"FSYM" %"FSYM"\n", realID[i],j, xpositions[i],
	      ypositions[i],zpositions[i]);
	}

      /* Total mass. */
 	if (nostars) {
      GroupProperties[0][j] += dmparticlemasses[i]*MassConv*totalmass;
    } else {
      GroupProperties[0][j] += particlemasses[i]*MassConv*totalmass;
    }
 
      /* Number of points. */
 
      GroupProperties[3][j] += 1.0;
 
      /* Max density (and position). */
	if (nostars) { 
      if (Density[i] > GroupProperties[4][j]) {
	GroupProperties[4][j] = Density[i];
	GroupProperties[5][j] = dmxpositions[i];
	GroupProperties[6][j] = dmypositions[i];
	GroupProperties[7][j] = dmzpositions[i];
      }
      /* S Skory - below */
      if (dmxpositions[i] < GroupProperties[8][j]) {
	GroupProperties[8][j] = dmxpositions[i];
      }
      if (dmypositions[i] < GroupProperties[10][j]) {
	GroupProperties[10][j] = dmypositions[i];
      }
      if (dmzpositions[i] < GroupProperties[12][j]) {
	GroupProperties[12][j] = dmzpositions[i];
      }

      if (dmxpositions[i] > GroupProperties[9][j]) {
	GroupProperties[9][j] = dmxpositions[i];
      }
      if (dmypositions[i] > GroupProperties[11][j]) {
	GroupProperties[11][j] = dmypositions[i];
      }
      if (dmzpositions[i] > GroupProperties[13][j]) {
	GroupProperties[13][j] = dmzpositions[i];
      }
    } else {
      if (Density[i] > GroupProperties[4][j]) {
	GroupProperties[4][j] = Density[i];
	GroupProperties[5][j] = xpositions[i];
	GroupProperties[6][j] = ypositions[i];
	GroupProperties[7][j] = zpositions[i];
      }
      /* S Skory - below */
      if (xpositions[i] < GroupProperties[8][j]) {
	GroupProperties[8][j] = xpositions[i];
      }
      if (ypositions[i] < GroupProperties[10][j]) {
	GroupProperties[10][j] = ypositions[i];
      }
      if (zpositions[i] < GroupProperties[12][j]) {
	GroupProperties[12][j] = zpositions[i];
      }

      if (xpositions[i] > GroupProperties[9][j]) {
	GroupProperties[9][j] = xpositions[i];
      }
      if (ypositions[i] > GroupProperties[11][j]) {
	GroupProperties[11][j] = ypositions[i];
      }
      if (zpositions[i] > GroupProperties[13][j]) {
	GroupProperties[13][j] = zpositions[i];
      }
      } // end if nostars if/else

    }

  fclose(blahfp); /* S Skory */

  /* Output group properties. */
 
  if ((fptr = fopen("HopAnalysis.out", "w")) == NULL) {
    fprintf(stderr, "Error opening regroup output HopAnalysis.out\n");
  }
 
  fprintf(fptr, "#Group     Mass      # part     max dens     x        y       z        xmin     xmax     ymin     ymax      zmin     zmax\n"); /* S Skory */
  for (j = 0; j < nGroups; j++) {
    fprintf(fptr, "%"ISYM"      ", j);
    for (i = 0; i < 14; i++) /* S Skory */
      if (i < 1 || i > 2)
	fprintf(fptr, " %.9"GSYM" ", GroupProperties[i][j]);
    fprintf(fptr, "\n");
  }
 
  fclose(fptr);
  

  delete [] dmxpositions;
  delete [] dmypositions;
  delete [] dmzpositions;

  delete [] dmparticlemasses;
  delete [] dmparticlecrtime;
  delete [] dmparticledyntime;
  delete [] dmparticlemetalfrac;
  delete [] dmparticleID;
  
}

/*------------------------- CosmologyParameters(char) --------------- 
 *
 *   Reads the hierarchy file and counts the number of grids.
 *   This will be used shortly thereafter to get all of the 
 *   grid info.
 *
 *-------------------------------------------------------------*/
int CosmologyParameters(char *filename, float *Hubble, float *OmegaMatterNow, float *boxsize){
  if(debug) fprintf(stdout,"in CosmologyParameters\n");

  FILE *file;
  char *line = new char[MAX_LINE_LENGTH];

  if(debug) fprintf(stdout,"data file: %s\n",filename);
    
  file = fopen(filename,"r");
  if (file==NULL) {
    fprintf (stderr,"Error opening file %s: exiting\n",filename);
    exit(1);
  }
 float Htemp, Otemp, Btemp;
 
  /* read through hierarchy file, counting # of grids
   lines that start with "Grid = " are the beginning of a new
   piece of grid info */
  while( fgets(line, MAX_LINE_LENGTH, file) != NULL){
    if(strncmp(line,"CosmologyHubbleConstantNow",26)==0) {
        sscanf(line,"CosmologyHubbleConstantNow = %f", Hubble);
    }

    if(strncmp(line,"CosmologyOmegaMatterNow",23)==0){
      sscanf(line,"CosmologyOmegaMatterNow    = %f",
	     OmegaMatterNow);
    }
    
    if(strncmp(line,"CosmologyComovingBoxSize",24)==0){
    	sscanf(line,"CosmologyComovingBoxSize   = %f",
    		boxsize);
    }
  }

  fclose(file);


  
  /* clean up dynamically allocated stuff */
  delete [] line;

  if(debug) fprintf(stdout,"exiting CosmologyParameters\n");

  /* return # of grids */
  return 1;
}


/*------------------------- NumberOfGrids(char) --------------- 
 *
 *   Reads the hierarchy file and counts the number of grids.
 *   This will be used shortly thereafter to get all of the 
 *   grid info.
 *
 *-------------------------------------------------------------*/
int NumberOfGrids(char *hierfilename){
  if(debug) fprintf(stderr,"in NumberOfGrids\n");

  FILE *hierfile;
  char *line = new char[MAX_LINE_LENGTH];
  int numgrids=0,npartthisgrid;
  int starparttemp=0;

  if(debug) fprintf(stderr,"hierarchy file: %s\n",hierfilename);
    
  hierfile = fopen(hierfilename,"r");
  if (hierfile==NULL) {
    fprintf (stderr,"Error opening file %s: exiting\n",hierfilename);
    exit(1);
  }

  /* read through hierarchy file, counting # of grids
   lines that start with "Grid = " are the beginning of a new
   piece of grid info */
  while( fgets(line, MAX_LINE_LENGTH, hierfile) != NULL){
    if(strncmp(line,"Grid = ",7)==0) numgrids++;

    if(strncmp(line,"NumberOfParticles",17)==0){
      sscanf(line,"NumberOfParticles   = %d",
	     &npartthisgrid);
      total_number_particles += npartthisgrid;
    }
    
    if(strncmp(line,"NumberOfStarParticles",21)==0){
    	sscanf(line,"NumberOfStarParticles      = %d",
    		&starparttemp);
    }
  }
    
  if (nostars) {
  	total_number_dm_particles = total_number_particles - starparttemp;
  }
  fclose(hierfile);

  if(debug) fprintf(stderr,"NumberOfGrids:  there are %i grids\n",numgrids);
  
  /* clean up dynamically allocated stuff */
  delete [] line;

  if(debug) fprintf(stderr,"exiting NumberOfGrids\n");

  /* return # of grids */
  return numgrids;  
}


/*----------------- GetGridInfo(int,char) ------------------------
 *
 *   Reads through grid file and retrieves all important information -
 *   grid dimensions, bounds, level, and file name, and puts them all 
 *   into arrays which were allocated in main().
 *
 *---------------------------------------------------------------- */
int GetGridInfo(int numberofgrids,char *hierfilename){
  if(debug) fprintf(stderr,"in GetGridInfo\n");
  FILE *hierfile;
  char *line = new char[MAX_LINE_LENGTH];
  char *gfname = new char[MAX_LINE_LENGTH];
  int grid=0,i;

  int griddimx,griddimy,griddimz,gridsex,gridsey,gridsez,
    grideex,grideey,grideez,nbfields, particleonly, ijunk;

  float fjunk;
  
  if(debug) fprintf(stderr,"hierarchy file: %s\n",hierfilename);

  /* open hierarchy file */
  hierfile = fopen(hierfilename,"r");

  /* read through hierarchy file, get grid info */
  while( fgets(line, MAX_LINE_LENGTH, hierfile) != NULL){

    /* lines that start with "Grid =" indicate a new set of grid info */
    if(strncmp(line,"Grid = ",7)==0){
    
      fgets(line, MAX_LINE_LENGTH, hierfile);  /* junk line */
      fgets(line, MAX_LINE_LENGTH, hierfile);  /* junk line */
      fgets(line, MAX_LINE_LENGTH, hierfile);  /* grid dim */
      sscanf(line,"GridDimension     = %d %d %d",
	     &griddimx,&griddimy,&griddimz);
      fgets(line, MAX_LINE_LENGTH, hierfile);  /* start index */
      sscanf(line,"GridStartIndex    = %d %d %d",
	     &gridsex,&gridsey,&gridsez);
      fgets(line, MAX_LINE_LENGTH, hierfile);  /* end index */
      sscanf(line,"GridEndIndex      = %d %d %d",
	     &grideex,&grideey,&grideez);
      fgets(line, MAX_LINE_LENGTH, hierfile);  /* left edge */
      sscanf(line,"GridLeftEdge      = %f %f %f",
	     &gridlex[grid],&fjunk,&fjunk);

      fgets(line, MAX_LINE_LENGTH, hierfile);  /* right edge */
      sscanf(line,"GridRightEdge     = %f %f %f",
	     &gridrex[grid],&fjunk,&fjunk);
      
      for(i=0;i<3;i++) fgets(line, MAX_LINE_LENGTH, hierfile);
      sscanf(line,"NumberOfBaryonFields = %d",&nbfields);

      if(nbfields==0){
	particleonly=1;

	fgets(line, MAX_LINE_LENGTH, hierfile);
	sscanf(line,"NumberOfParticles   = %d",
	       &gridnump[grid]);  /* number of particles on this grid */

	fgets(line, MAX_LINE_LENGTH, hierfile);
	sscanf(line,"ParticleFileName = %s",gridfilenames[grid]);

	if(debug) fprintf(stderr,"GetGridInfo:  No baryon fields!  %s %d\n",gridfilenames[grid],gridnump[grid]);

      } else {
	/* "junk" lines */
	for(i=0;i<2;i++) fgets(line, MAX_LINE_LENGTH, hierfile);
	
	/* grid file name */
	sscanf(line,"BaryonFileName = %s",gridfilenames[grid]);

	/* "junk" lines */
	for(i=0;i<5;i++) fgets(line, MAX_LINE_LENGTH, hierfile);

	sscanf(line,"NumberOfParticles   = %d",
	       &gridnump[grid]);  /* number of particles on this grid */

	if(debug) fprintf(stderr,"GetGridInfo: WITH baryon fields!  %s %d\n",gridfilenames[grid],gridnump[grid]);

      }

      /* grid dims - buffers stripped */
      griddx[grid] = 1+grideex-gridsex;

      grid++;  /* increment grid number! */
    }
  }

  fclose(hierfile);

  /* clean up dynamically allocated stuff */
  delete [] line;
  delete [] gfname;

  if(debug) fprintf(stderr,"exiting GetGridInfo\n");

  /* return code - number of grids! */
  return grid;
}

/*----------------- AddParticlesToArrays(int) ------------------------
 *
 *   Adds the particles to arrays 
 *
 *---------------------------------------------------------------- */

int AddParticlesToArrays(int gridnum){
  if(debug) fprintf(stderr,"in AddParticlesToArrays %d\n",gridnum);

  /* are there no particles on this grid?  if not, exit! */
  if( gridnump[gridnum] == 0) return SUCCESS;

  /* we only need one of all of the following */
  hid_t file_id;
  hid_t       mem_type_id;
  hid_t       mem_dsp_id, file_dsp_id;
  hid_t       dsp_id;
  hid_t       typ_id;

  hsize_t     size;
  hsize_t     dims[4];
  hsize_t     xdims[4];
  hsize_t     maxdims[4];

  herr_t      h5_status;
  herr_t      h5_error = -1;

  hid_t       pposx_dset_id,pposy_dset_id, pposz_dset_id,
  	pmass_dset_id,pcrtime_dset_id,
  	 pdyntime_dset_id,pmfrac_dset_id,pID_dset_id;

  int i, ndims;
  float *particle_mass_buffer, *particle_crtime_buffer,
  	*particle_dyntime_buffer, *particle_mfrac_buffer;;
  double *particle_posx_buffer, *particle_posy_buffer, 
    *particle_posz_buffer;
  int numberofparticles,totalnumberofstars=0;
  long *particle_ID_buffer;
  char groupname[20];
  hid_t group_id;
  
    float cellvolume, deltax;

  deltax = (gridrex[gridnum]-gridlex[gridnum])/((float)griddx[gridnum]);

  cellvolume = deltax*deltax*deltax;

  if(debug) fprintf(stderr,"GridFoo:  %f %f %d %e %e\n",
		    gridlex[gridnum],gridrex[gridnum],griddx[gridnum],
		    deltax,cellvolume);

  /* open file - only once */
  if(debug) fprintf(stderr,"AddParticlesToArrays %s %d\n",gridfilenames[gridnum],
	  gridnump[gridnum]);

  file_id = H5Fopen(gridfilenames[gridnum], H5F_ACC_RDWR, H5P_DEFAULT);
  if (file_id == h5_error){
    fprintf(stderr,"%s:%d ERROR OPENING HDF5 FILE [%s]: EXITING\n",
	    __FILE__,__LINE__,gridfilenames[gridnum]);
    exit(1);
  }
  assert( file_id != h5_error );
  
  /* open group in packed AMR file if nopacked = 0 */
  if (nopacked == 0) {
  sprintf(groupname,"/Grid%08d",gridnum+1);
  group_id = H5Gopen(file_id, groupname);
  assert( group_id != h5_error);
  } else {
  group_id = file_id;
  }

  /*------------------------------ NOTE -------------------------------
    The particle position stuff is 64 bit, whereas the creation time and
    particle masses are 32 bit.  Get particle positions and then go back and
    get the creation time and particle mass info.  */

  /* get particle positions - 64 bit! */
  pposx_dset_id = H5Dopen(group_id, "particle_position_x");
  assert( pposx_dset_id != h5_error);

  pposy_dset_id = H5Dopen(group_id, "particle_position_y");  
  assert( pposy_dset_id != h5_error);

  pposz_dset_id = H5Dopen(group_id, "particle_position_z");  
  assert( pposz_dset_id != h5_error);

  pID_dset_id = H5Dopen(group_id, "particle_index");  
  assert( pID_dset_id != h5_error);

  /* open particle pos x dataspace (to get dimensions) */
  dsp_id = H5Dget_space(pposx_dset_id);
  assert( dsp_id != h5_error );

  /* get data type (only once!) */
  typ_id = H5Dget_type(pposx_dset_id);
  assert( typ_id != h5_error );

  /* get dimensional information from dataspace (only once) */
  ndims = H5Sget_simple_extent_dims(dsp_id, xdims, maxdims);

  /* from the dimensional information, calculate the size of the buffer.
   only once! */
  size = 1;
  if(debug) fprintf(stderr,"Ndims %d\n",ndims);
  for ( i = 0; i < ndims; i++)
    {
      dims[i] = xdims[i];
      size = size * dims[i];
      if(debug) fprintf(stderr," Dim %d\n", (int) xdims[i]);
    }
  if(debug) fprintf(stderr,"Size %d\n", (int) size);

  file_dsp_id = H5Screate_simple(ndims, dims, NULL);
  assert( file_dsp_id != h5_error );

  mem_dsp_id = H5Screate_simple(1, &size, NULL);
  assert( mem_dsp_id != h5_error );


      if(debug) fprintf(stderr,"(3)\n");
      particle_posx_buffer = new double[(int) size];
      particle_posy_buffer = new double[(int) size];
      particle_posz_buffer = new double[(int) size];
      particle_ID_buffer = new long[(int) size];

      /* read particle position x field into an array */
      mem_type_id = H5Dget_type(pposx_dset_id);
      h5_status = H5Dread(pposx_dset_id, mem_type_id, 
			  mem_dsp_id, file_dsp_id, 
			  H5P_DEFAULT,particle_posx_buffer);
      if(debug) fprintf(stderr,"float read status %d for particle pos x field\n", 
			(int) h5_status);
      assert( h5_status != h5_error ); 

      /* read particle position y field into an array */
      mem_type_id = H5Dget_type(pposy_dset_id);
      h5_status = H5Dread(pposy_dset_id, mem_type_id, 
			  mem_dsp_id, file_dsp_id, 
			  H5P_DEFAULT,particle_posy_buffer);
      if(debug) fprintf(stderr,"float read status %d for particle pos y field\n", 
			(int) h5_status);
      assert( h5_status != h5_error ); 

      /* read particle position z field into an array */
      mem_type_id = H5Dget_type(pposz_dset_id);
      h5_status = H5Dread(pposz_dset_id, mem_type_id, 
			  mem_dsp_id, file_dsp_id, 
			  H5P_DEFAULT,particle_posz_buffer);
      if(debug) fprintf(stderr,"float read status %d for particle pos z field\n", 
			(int) h5_status);
      assert( h5_status != h5_error ); 

      /* read particle ID field into an array */
      mem_type_id = H5Dget_type(pID_dset_id);
      h5_status = H5Dread(pID_dset_id, mem_type_id, 
			  mem_dsp_id, file_dsp_id, 
			  H5P_DEFAULT,particle_ID_buffer);
      if(debug) fprintf(stderr,"int read status %d for particle index field\n", 
			(int) h5_status);
      assert( h5_status != h5_error ); 


  /* close hdf5 file, doing appropriate error checking */
  h5_status = H5Sclose(dsp_id);
  assert( h5_status != h5_error );
  
  h5_status = H5Tclose(typ_id);
  assert( h5_status != h5_error );

  h5_status = H5Sclose(mem_dsp_id);
  assert( h5_status != h5_error );

  h5_status = H5Sclose(file_dsp_id);
  assert( h5_status != h5_error );

  /* close all position datasets */
  h5_status = H5Dclose(pposx_dset_id);
  assert( h5_status != h5_error );

  h5_status = H5Dclose(pposy_dset_id);
  assert( h5_status != h5_error );

  h5_status = H5Dclose(pposz_dset_id);
  assert( h5_status != h5_error );

  h5_status = H5Dclose(pID_dset_id);
  assert( h5_status != h5_error );


  /* ------------------ get particle mass and creation time info */

  pmass_dset_id = H5Dopen(group_id, "particle_mass");  
  assert( pmass_dset_id != h5_error);
  
  pcrtime_dset_id = H5Dopen(group_id, "creation_time");  
  assert( pmass_dset_id != h5_error);

  pdyntime_dset_id = H5Dopen(group_id, "dynamical_time");  
  assert( pmass_dset_id != h5_error);

  pmfrac_dset_id = H5Dopen(group_id, "metallicity_fraction");  
  assert( pmass_dset_id != h5_error);

  /* open particle mass dataspace (to get dimensions) */
  /* only once! */
  dsp_id = H5Dget_space(pmass_dset_id);
  assert( dsp_id != h5_error );

  /* get data type (only once!) */
  typ_id = H5Dget_type(pmass_dset_id);
  assert( typ_id != h5_error );

  /* get dimensional information from dataspace (only once) */
  ndims = H5Sget_simple_extent_dims(dsp_id, xdims, maxdims);

  /* from the dimensional information, calculate the size of the buffer.
   only once! */
  size = 1;
  if(debug) fprintf(stderr,"Ndims %d\n",ndims);
  for ( i = 0; i < ndims; i++)
    {
      dims[i] = xdims[i];
      size = size * dims[i];
      if(debug) fprintf(stderr," Dim %d\n", (int) xdims[i]);
    }
  if(debug) fprintf(stderr,"Size %d\n", (int) size);

  file_dsp_id = H5Screate_simple(ndims, dims, NULL);
  assert( file_dsp_id != h5_error );

  mem_dsp_id = H5Screate_simple(1, &size, NULL);
  assert( mem_dsp_id != h5_error );

      particle_mass_buffer = new float[(int) size];
      particle_crtime_buffer = new float[(int) size];
      particle_dyntime_buffer = new float[(int) size];
      particle_mfrac_buffer = new float[(int) size];
      

      /* read particle mass field into an array */
      mem_type_id = H5Dget_type(pmass_dset_id);
      h5_status = H5Dread(pmass_dset_id, mem_type_id, 
			  mem_dsp_id, file_dsp_id, 
			  H5P_DEFAULT,particle_mass_buffer);
      if(debug) fprintf(stderr,"float read status %d for particle mass field\n", 
			(int) h5_status);
      assert( h5_status != h5_error );
      
      /* read particle crtime field into an array */
      mem_type_id = H5Dget_type(pcrtime_dset_id);
      h5_status = H5Dread(pcrtime_dset_id, mem_type_id, 
			  mem_dsp_id, file_dsp_id, 
			  H5P_DEFAULT,particle_crtime_buffer);
      if(debug) fprintf(stderr,"float read status %d for particle crtime field\n", 
			(int) h5_status);
      assert( h5_status != h5_error ); 


      /* read particle dyntime field into an array */
      mem_type_id = H5Dget_type(pdyntime_dset_id);
      h5_status = H5Dread(pdyntime_dset_id, mem_type_id, 
			  mem_dsp_id, file_dsp_id, 
			  H5P_DEFAULT,particle_dyntime_buffer);
      if(debug) fprintf(stderr,"float read status %d for particle dyntime field\n", 
			(int) h5_status);
      assert( h5_status != h5_error ); 

      /* read particle mfrac field into an array */
      mem_type_id = H5Dget_type(pmfrac_dset_id);
      h5_status = H5Dread(pmfrac_dset_id, mem_type_id, 
			  mem_dsp_id, file_dsp_id, 
			  H5P_DEFAULT,particle_mfrac_buffer);
      if(debug) fprintf(stderr,"float read status %d for particle mfrac field\n", 
			(int) h5_status);
      assert( h5_status != h5_error ); 


  /* close hdf5 data, doing appropriate error checking */
  h5_status = H5Sclose(dsp_id);
  assert( h5_status != h5_error );
  
  h5_status = H5Tclose(typ_id);
  assert( h5_status != h5_error );

  h5_status = H5Sclose(mem_dsp_id);
  assert( h5_status != h5_error );

  h5_status = H5Sclose(file_dsp_id);
  assert( h5_status != h5_error );
  
  h5_status = H5Dclose(pmass_dset_id);
  assert( h5_status != h5_error );

  h5_status = H5Dclose(pcrtime_dset_id);
  assert( h5_status != h5_error );

  h5_status = H5Dclose(pdyntime_dset_id);
  assert( h5_status != h5_error );

  h5_status = H5Dclose(pmfrac_dset_id);
  assert( h5_status != h5_error );

	/* close group if using packed data*/
	if (nopacked == 0) {
	h5_status = H5Gclose(group_id);
	assert(h5_status != h5_error);
	}

  /* close file */
  h5_status = H5Fclose(file_id);
  assert( h5_status != h5_error );

  numberofparticles = ((int) size);

  assert( numberofparticles == gridnump[gridnum] );

  for(i=0; i<numberofparticles; i++){
    /*  for(i=array_position; i<array_position+numberofparticles; i++){ */

    xpositions[array_position] = (double) particle_posx_buffer[i];
    ypositions[array_position] = (double) particle_posy_buffer[i];
    zpositions[array_position] = (double) particle_posz_buffer[i];
    particlemasses[array_position] = particle_mass_buffer[i]*cellvolume;
    particlecrtime[array_position] = particle_crtime_buffer[i];
    particledyntime[array_position] = particle_dyntime_buffer[i];
    particlemetalfrac[array_position] = particle_mfrac_buffer[i];
    particleID[array_position] = particle_ID_buffer[i];
    
    array_position++;
    particle_count++;

    if(array_position >= total_number_particles)
      fprintf(stderr,"HEY KIDS!  TOTAL NUM PARTICLES REACHED!\n\n");


  }

  /* clean up arrays! */
  delete [] particle_posx_buffer;
  delete [] particle_posy_buffer;
  delete [] particle_posz_buffer;
  delete [] particle_mass_buffer;
  delete [] particle_crtime_buffer;
  delete [] particle_dyntime_buffer;
  delete [] particle_mfrac_buffer;
  delete [] particle_ID_buffer;

  if(debug) fprintf(stderr,"exiting AddParticlesToArrays\n");
  return SUCCESS;


}
