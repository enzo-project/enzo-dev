#include "performance.h"
#include <math.h>
#include <stdio.h>
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "DaveTools.h"
#include "CosmologyParameters.h"

#ifdef NEW_DIVB
#include "AthenaObj.h"
#endif //NEW_DIVB
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
#ifdef HAOXU
int FindField(int f, int farray[], int n);
#endif
extern "C" void FORTRAN_NAME(curl_of_e)(float *bx, float *by, float *bz,
					float *ex, float *ey, float *ez, 
					float *dx, float *dy, float *dz, 
					int *idim, int *jdim, int *kdim,
					int *i1, int *i2, int *j1, int *j2, int *k1, int *k2,
					float *dt, MHD_Centering *method);

#ifdef HAOXU
extern "C" void FORTRAN_NAME(mhd_harten_cosmology)
  (float *d, float *e, float *vx, float *vy, float *vz,
   float *bxc, float *byc, float *bzc, 
   int * gravityon, float *gr_ax,float *gr_ay,float *gr_az,
   float *bfluxx1, float *bfluxy1, float *bfluxz1,
   float *bfluxx2, float *bfluxy2, float *bfluxz2,
   float *Fd, float *Fe, float *Fvx, float *Fvy, float *Fvz,
   int* FluxExtents, int* TotalFluxSize, int* NumberOfSubgrids,
   float *dx, float *dy, float *dz, 
   int *idim, int *jdim, int *kdim, 
   int *i1, int *i2, int *j1, int *j2, int *k1, int *k2,
   float *dt, float *gamma, int *nhy, int *rank, int *level, int * grid,
   int* MHD_Hack2d,
   FLOAT *a,   // Expansion factors  xh 
   int *idual, float *ge, int *idiff, //dual energy, flag for diffusion xh
   float *premin); //mininum pressure

// To call MHD solver of Shengtai
extern "C" void FORTRAN_NAME( mhd_li)
  (float *d, float *e, float *vx, float *vy, float *vz,
   float *bxc, float *byc, float *bzc,
   int *gravityon, float *gr_ax, float *gr_ay, float *gr_az,
   float *fx1, float *fy1, float *fz1,
   float *fx2, float *fy2, float *fz2,
   float *fd, float *fe, float *fvx, float *fvy, float *fvz, float *fge,
   int *fluxextents, int *totalfluxsize, int *nsubgrids,
   float *dx, float *dy, float *dz, int *idim, int *jdim, int *kdim,
   int *i1, int *i2, int *j1, int *j2, int *k1, int *k2, float *dt, float *gamma,
   int *nhy, int *rank, int *level, int *grid,  
   FLOAT  *a, int *idual, float *ge, int *idiff,
   float *premin, int *cosmos_equation, int *solverparameters, int *EquationOfState, float* IsothermalSoundSpeed);

//multi-species
extern "C" void FORTRAN_NAME( mhd_li_ms)
       (float *d, float *e, float *vx, float *vy, float *vz,
        float *bxc, float *byc, float *bzc,
       int *gravityon, float *gr_ax, float *gr_ay, float *gr_az,
       float *fx1, float *fy1, float *fz1,
       float *fx2, float *fy2, float *fz2,
       float *fd, float *fe, float *fvx, float *fvy, float *fvz, float *fge,
       int *fluxextents, int *totalfluxsize, int *nsubgrids,
       float *dx, float *dy, float *dz, int *idim, int *jdim, int *kdim,
       int *i1, int *i2, int *j1, int *j2, int *k1, int *k2, float *dt, float *gamma,
       int *nhy, int *rank, int *level, int *grid,
       FLOAT  *a, int *idual, float *ge, int *idiff,
       float *premin, int *cosmos_equation, int *solverparameters,
       int *numberofcolor,
       float *c0,float *c1,float *c2, float *c3, float *c4,float *c5,
        float *c6,float *c7,float *c8, float *c9, float *c10, float *c11,
        float *c12, float *c13, float *c14,
       float *fc0,float *fc1,float *fc2, float *fc3, float *fc4,float *fc5,
        float *fc6,float *fc7,float *fc8, float *fc9, float *fc10, float *fc11,
        float *fc12, float *fc13, float *fc14); 

#endif /*HAOXU */

#ifdef BIERMANN
extern "C" void FORTRAN_NAME (create_e_biermann)(float *fx1, float *fy1, float *fz1,
                                       float *fx2, float *fy2, float *fz2,
                                       float *ex, float *ey, float *ez,
                                       float *dx, float *dy, float *dz,
                                       float *density, float *ge,
                                       int *idim, int *jdim, int *kdim,
                                       int *i1, int *i2, int *j1, int *j2, int *k1, int *k2,
                                       float * dt, int *ProjectE,
                     float *gamma, float *c,float *mh,float *e,float *chi, FLOAT *a);

int MHDCosmologyGetUnits(float *DensityUnits, float *LengthUnits,
                       float *TemperatureUnits, float *TimeUnits,
                       float *VelocityUnits, FLOAT Time,
                       float *BFieldUnits);
#endif //BIERMANN

extern "C" void FORTRAN_NAME(mhd_harten)
  (float *d, float *e, float *vx, float *vy, float *vz,
   float *bxc, float *byc, float *bzc, 
   int * gravityon, float *gr_ax,float *gr_ay,float *gr_az,
   float *bfluxx1, float *bfluxy1, float *bfluxz1,
   float *bfluxx2, float *bfluxy2, float *bfluxz2,
   float *Fd, float *Fe, float *Fvx, float *Fvy, float *Fvz,
   int* FluxExtents, int* TotalFluxSize, int* NumberOfSubgrids,
   float *dx, float *dy, float *dz, 
   int *idim, int *jdim, int *kdim, 
   int *i1, int *i2, int *j1, int *j2, int *k1, int *k2,
   float *dt, float *gamma, int *nhy, int *rank, int *level, int * grid
   , int* MHD_Hack2d
#ifdef NSS
   , int * recon, int * riemann, float * eps
#endif
);

extern "C" void FORTRAN_NAME(create_e)(float *fx1, float *fy1, float *fz1,
				       float *fx2, float *fy2, float *fz2,
				       float *ex, float *ey, float *ez,
				       int *idim, int *jdim, int *kdim,
				       int *i1, int *i2, int *j1, int *j2, int *k1, int *k2,
				       float * dt, int * ProjectE);

extern "C" void FORTRAN_NAME(divb_rj)(float *fx1, float *fy1, float *fz1,
				      float *fx2, float *fy2, float *fz2,
				      float *bxb, float *byb, float *bzb, 
				      float *bxc, float *byc, float *bzc, 
				      float *dx, float *dy, float *dz, float *dt,
				      int *nx, int *ny, int *nz);

int grid::SolveMHDEquations(int CycleNumber, int NumberOfSubgrids, 
	    fluxes *SubgridFluxes[], ExternalBoundary *Exterior, int level, int grid)
  /*begin*/
{

  //MHD_Equation = 2; 

  int oot = (MyProcessorNumber == ROOT_PROCESSOR ) ? TRUE : FALSE;

  JBMEM_MESSAGE(MyProcessorNumber,"jb: SMHDEnter");  
  if( ProcessorNumber != MyProcessorNumber )
    return SUCCESS;
  
  if( TVtool("start of SMHD") == FAIL ){
    fprintf(stderr,"problem at the beginning of SMHD.\n");
    return FAIL;
  }
  if( this->CheckForNans("beginnin of SMHD") == FAIL ) {
    fprintf(stderr,"startMHD failure\n");
    return FAIL;
  }

  wall_time("Start SMHD");
  this->DebugCheck("SMHD: Before");
  /*
  fprintf(stderr,"===  SMHD n = %d L = %d g = %d proc = %d dt = %15.12e ===\n",
	  CycleNumber, level, grid, MyProcessorNumber, dtFixed);
  */
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  int i, j, k, n, dim, field, axis, face, index, z, coord, Findex, Cindex, index2, size;
  int subgrid, TotalFluxSize=0, SizeOfFace=1, SizeOfSubgrid = 0;
  int *SizeOtherSubgrids = new int [NumberOfSubgrids];
  int *TotalOffset[3], *SizeOtherDims[3], Dim[3]={1,1,1};
  int testBflux = FALSE;
  int nx = GridDimension[0] - 2*DEFAULT_GHOST_ZONES, 
    ny = GridDimension[1] - 2*DEFAULT_GHOST_ZONES, 
    nz = GridDimension[2] - 2*DEFAULT_GHOST_ZONES;
  char basename[20];
  int MyGridNumber = grid+1; //THe grid number for output.  Use an integer, level, processor, bananna, whatever.
  float dtUsed = 1.0;
  int UseDT = 1;

  //
  //Figure out which BaryonField contains which hydro variable.
  //
  
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				       Vel3Num, TENum) == FAIL) {
    fprintf(stderr, "Error in IdentifyPhysicalQuantities.\n");
    return FAIL;
  }
  
#ifdef HAOXU //multi-species
  /* If multi-species being used, then treat them as colour variables
     (note: the solver has been modified to treat these as density vars). */
  
  int NumberOfColours = 0, ColourNum, coloff[MAX_COLOR];
  float *colourpt = NULL;
  
  if (MultiSpecies > 0) {
    
    NumberOfColours = 6 + 3*(MultiSpecies-1);
    
    if ((ColourNum =
	 FindField(ElectronDensity, FieldType, NumberOfBaryonFields)) <0) {
      fprintf(stderr, "Could not find ElectronDensity.\n");
      return FAIL;
    }
    
    /* Set Offsets from the first field (here assumed to be the electron
       density) to the other multi-species fields. */
    
    colourpt = BaryonField[ColourNum];
    
    for (i = 0; i < NumberOfColours; i++)
      coloff[i] = BaryonField[ColourNum+i] - colourpt;
    
  }
  
  /* Add metallicity as a colour variable. */
  
  int MetalNum;
  
  if ((MetalNum = FindField(Metallicity, FieldType,
			    NumberOfBaryonFields)) != -1) {
    if (colourpt == NULL) {
      colourpt = BaryonField[MetalNum];
      ColourNum = MetalNum;
    }
    
    coloff[NumberOfColours++] = BaryonField[MetalNum  ] - colourpt;
    coloff[NumberOfColours++] = BaryonField[MetalNum+1] - colourpt;
    coloff[NumberOfColours++] = BaryonField[MetalNum+2] - colourpt;
    
  }
  
#endif //HAOXU
  
  // The fluxe macro is used to map the 5 dimensional, non-rectangular array to a 1 dimensional array
  // that fortran can deal with. 
  
#define fluxe(dim,coord,face,end,sub) FluxExtents[dim+3*(coord+3*(face+2*(end+2*sub)))]
  int * FluxExtents = new int[3*3*2*2*NumberOfSubgrids];
  int * FluxDims[3][3];
  float ** FluxFromSolver = new float*[NumberOfBaryonFields];

  JBMEM_MESSAGE(MyProcessorNumber,"jb: SMHDBeforeFluxAllocate");  

  for(dim=0;dim<3;dim++)
    for(coord=0;coord<3;coord++)
      for(face=0;face<2;face++)
	for(n=0;n<2;n++){
	  FluxDims[dim][coord] = new int[NumberOfSubgrids];
	  for(subgrid=0;subgrid<NumberOfSubgrids;subgrid++){
	    fluxe(dim,coord,face,n,subgrid)=0;
	    FluxDims[dim][coord][subgrid]=0;
	  }
	}

   for (subgrid = 0; subgrid < NumberOfSubgrids; subgrid++) {
    for (dim = 0; dim < GridRank; dim++)  {

      /* Compute local indicies and dimensions.  These will be used in the solver */
      /* The 1 is for fortran indicies. */
      for(coord=0;coord<3;coord++){
	fluxe(dim,coord,0,0,subgrid)=1+
	  SubgridFluxes[subgrid]->LeftFluxStartGlobalIndex[dim][coord]-
	    nlongint((CellLeftEdge[coord][0] - DomainLeftEdge[coord])/CellWidth[coord][0]);
	fluxe(dim,coord,0,1,subgrid)=1+
	  SubgridFluxes[subgrid]->LeftFluxEndGlobalIndex[dim][coord]-
	    nlongint((CellLeftEdge[coord][0] - DomainLeftEdge[coord])/CellWidth[coord][0]);
	fluxe(dim,coord,1,0,subgrid)=1+
	  SubgridFluxes[subgrid]->RightFluxStartGlobalIndex[dim][coord]-
	  nlongint((CellLeftEdge[coord][0] - DomainLeftEdge[coord])/CellWidth[coord][0]);
	fluxe(dim,coord,1,1,subgrid)=1+
	  SubgridFluxes[subgrid]->RightFluxEndGlobalIndex[dim][coord]-
	    nlongint((CellLeftEdge[coord][0] - DomainLeftEdge[coord])/CellWidth[coord][0]);
	FluxDims[dim][coord][subgrid]=
	  fluxe(dim,coord,0,1,subgrid)-fluxe(dim,coord,0,0,subgrid)+1;

      }      

      // Because we need the actual boundary... See MagneticField indexing convention.
      fluxe(dim,dim,1,0,subgrid)++;
      fluxe(dim,dim,1,1,subgrid)++;

      /* compute size (in floats) of flux storage */
      size = 1;
      for (j = 0; j < GridRank; j++)
	size *= SubgridFluxes[subgrid]->LeftFluxEndGlobalIndex[dim][j] -
	  SubgridFluxes[subgrid]->LeftFluxStartGlobalIndex[dim][j] + 1;
      
      /* set unused dims */
      
      for (j = GridRank; j < 3; j++) {
	SubgridFluxes[subgrid]->LeftFluxStartGlobalIndex[dim][j] = 0;
	SubgridFluxes[subgrid]->LeftFluxEndGlobalIndex[dim][j] = 0;
	SubgridFluxes[subgrid]->RightFluxStartGlobalIndex[dim][j] = 0;
	SubgridFluxes[subgrid]->RightFluxEndGlobalIndex[dim][j] = 0;
      }
      
      /* Allocate space (if necessary). */
      //JBMEM_MESSAGE(MyProcessorNumber,"jb: SMHDBeforeFluxAllocate loop");  
      for (field = 0; field < NumberOfBaryonFields; field++) {
	if (SubgridFluxes[subgrid]->LeftFluxes[field][dim] == NULL)
	  SubgridFluxes[subgrid]->LeftFluxes[field][dim]  = new float[size];
	if (SubgridFluxes[subgrid]->RightFluxes[field][dim] == NULL)
	  SubgridFluxes[subgrid]->RightFluxes[field][dim] = new float[size];
	for (n = 0; n < size; n++) {
	  SubgridFluxes[subgrid]->LeftFluxes[field][dim][n] = 0;
	  SubgridFluxes[subgrid]->RightFluxes[field][dim][n] = 0;
	}
      }
      //JBMEM_MESSAGE(MyProcessorNumber,"jb: SMHDBeforeFluxAllocate loop");        
      for (field = NumberOfBaryonFields; field < MAX_NUMBER_OF_BARYON_FIELDS;
	   field++) {
	SubgridFluxes[subgrid]->LeftFluxes[field][dim] = NULL;
	SubgridFluxes[subgrid]->RightFluxes[field][dim] = NULL;
      }
      
    }  // next dimension
   
    /* Clean up the remaining faces */
    
    for (dim = GridRank; dim < 3; dim++)
      for (field = 0; field < MAX_NUMBER_OF_BARYON_FIELDS; field++) {
	SubgridFluxes[subgrid]->LeftFluxes[field][dim] = NULL;
	SubgridFluxes[subgrid]->RightFluxes[field][dim] = NULL;
      }
    
   } // end of loop over subgrids
   
   //Allocate the Giant Flux Array and some Indexing Arrays
   //to keep track of positions in the Giant Flux Array

   for(dim=0;dim<3;dim++){
     TotalOffset[dim]=new int[NumberOfSubgrids];
     SizeOtherDims[dim]=new int[NumberOfSubgrids];
   }
  JBMEM_MESSAGE(MyProcessorNumber,"jb: SMHDAfterFluxCreation");  

   SizeOtherSubgrids[0]=0;

   for(subgrid=0;subgrid<NumberOfSubgrids;subgrid++){
     SizeOfSubgrid=0;
     SizeOtherDims[0][subgrid]=0;

     for(dim=0;dim<3;dim++){

       SizeOfFace=1;

       for(coord=0;coord<3;coord++)
	 SizeOfFace *= 2*FluxDims[dim][coord][subgrid];
       
       TotalFluxSize += SizeOfFace;
       SizeOfSubgrid += SizeOfFace;
       if(dim<2)
	 SizeOtherDims[dim+1][subgrid]=SizeOtherDims[dim][subgrid]+SizeOfFace;

     }//dim
     if(subgrid<NumberOfSubgrids-1)
       SizeOtherSubgrids[subgrid+1] = SizeOtherSubgrids[subgrid]+SizeOfSubgrid;

     for(dim=0;dim<3;dim++){
       TotalOffset[dim][subgrid]=SizeOtherSubgrids[subgrid]+SizeOtherDims[dim][subgrid];
     }
     
   }//subgrid

  JBMEM_MESSAGE(MyProcessorNumber,"jb: SMHDFluxFromSolver0");  
   for(field=0;field<NumberOfBaryonFields;field++){
     FluxFromSolver[field] = new float[TotalFluxSize];
     for(i=0;i<TotalFluxSize; i++)
       FluxFromSolver[field][i]=0.0;
   }

  JBMEM_MESSAGE(MyProcessorNumber,"jb: SMHDFluxFromSolver1");  
   //
   // Allocate Electric Field.  
   //

   for(field=0;field<3;field++){
     
     if(ElectricField[field] != NULL ) {
       delete [] ElectricField[field];
     }
     
     ElectricField[field] = new float[ElectricSize[field]];
     for(i=0;i<ElectricSize[field]; i++) ElectricField[field][i] = 0.0;
     
     
     if(MagneticField[field]==NULL){
       fprintf(stderr, "========== Solve MHD create Magnetic Field\nShit, that's not good.\n");
       return FAIL;
     }
     
   }


  float *MagneticFlux[3][2];

  int GridOffset[MAX_DIMENSION];
  
  //  
  //Setup arrays for Magnetic Flux
  //This is a seperate object from the Hydro flux correction because the two
  //are dealt with in different ways.  This flux is used only to calculate the Electric Field.
  //Note that MagneticFluxActual is used in order to keep the flux contiguous in memory.
  //This was done for the Athena CT machinery, which wants contiguous datasets.
  //
  
  //I'm a driver, I'm a winner.
  //Things are gonna change, I can feel it.
  //Soy un Perdador.  

  //each magnetic field has components from the flux of two 'other' axis
  JBMEM_MESSAGE(MyProcessorNumber,"jb: SMHDMagFlux0");  

#ifdef NEW_DIVB

  for(field=0;field<3;field++){
      MagneticFlux[field][0] = new float[2*MagneticSize[field]];
      MagneticFlux[field][1] =  MagneticFlux[field][0] +MagneticSize[field];
  }

#else //NEW_DIVB
  for(field=0;field<3;field++)
    for(axis=0;axis<2;axis++)
      MagneticFlux[field][axis] = new float[MagneticSize[field]];
#endif //NEW_DIVB
  for(field=0;field<3;field++)
    for( i=0;i<MagneticSize[field]; i++){
      MagneticFlux[field][0][i] = 0.0;
      MagneticFlux[field][1][i] = 0.0;
    }

  //Global memory spike here.  WHy isn't Magnetic Flux deleted?  Or is this a tool bug?
  JBMEM_MESSAGE(MyProcessorNumber,"jb: SMHDMagFlux1");  

  if( CheckForSingleGridDump(29) == TRUE) {
    sprintf(basename, "data29%02d%d.grid",CycleNumber, level);

    FILE *dummy = fopen(basename, "a");    
    if( this->WriteGrid(dummy, basename, MyGridNumber) == FAIL ){
      fprintf(stderr, "Shit.  Problem with Write Grid in SMHD.\n");
      return FAIL;
    }
    fclose(dummy);

  }  

#ifdef OLD_CENTER
  if( this->CenterMagneticField() == FAIL ) {
    fprintf(stderr," error with CenterMagneticField, first call \n");
    return FAIL;
  }
#endif //OLD_CENTER

  if( CheckForSingleGridDump(30) == TRUE) {

    sprintf(basename, "data30%02d%d.grid",CycleNumber, level);

    FILE *dummy = fopen(basename, "a");    
    if( this->WriteGrid(dummy, basename, MyGridNumber) == FAIL ){
      fprintf(stderr, "Shit.  Problem with Write Grid in SMHD.\n");
      return FAIL;
    }
    fclose(dummy);

  }  

 /* Set minimum support. */
 
    float MinimumSupportEnergyCoefficient = 0;
    if (UseMinimumPressureSupport == TRUE && level>MaximumRefinementLevel-1)
      if (this->SetMinimumSupport(MinimumSupportEnergyCoefficient) == FAIL) {
        fprintf(stderr, "Error in grid->SetMinimumSupport,\n");
        return FAIL;
      }


  int GravityOn = 0;
  if (SelfGravity || UniformGravity || PointSourceGravity)
    GravityOn = 1;

  /* if comoving coordinates, compute a and dadt at Time, Time + 0.5dt*/ 

  FLOAT a[4];
 
  if(ComovingCoordinates==1)
  {
      if (CosmologyComputeExpansionFactor(Time, &a[0], &a[1]) 
	== FAIL) {
      fprintf(stderr, "Error in CosmologyComputeExpansionFactor.\n");
      return FAIL;
      }
      if (CosmologyComputeExpansionFactor(Time+(FLOAT)0.5*dtFixed,
        &a[2], &a[3]) == FAIL) {
      fprintf(stderr, "Error in CosmologyComputeExpansionFactor.\n");
      return FAIL;
	}
  }else{
  a[0] = 1.0;
  a[1] = 0.0;
  a[2] = 1.0;
  a[3] = 0.0;
  }

    /* Create a cell width array to pass (and convert to absolute coords). */
/*
    float *CellWidthTemp[MAX_DIMENSION];
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      CellWidthTemp[dim] = new float[GridDimension[dim]];
      for (i = 0; i < GridDimension[dim]; i++)
        if (dim < GridRank)
          CellWidthTemp[dim][i] = float(a[2]*CellWidth[dim][i]);
        else
          CellWidthTemp[dim][i] = 1.0;
    }
*/

#ifdef BIERMANN
    // Biermann battery constants in cgs units, convert to enzo units later
    float speedoflight = 3.0e10;
    float hydrogenmass = 1.6733e-24;
    float electroncharge = 4.803e-10;
    float chi=1.0;

 /* Compute Units. */

  float DensityUnits = 1, LengthUnits = 1, TemperatureUnits = 1, TimeUnits = 1,
     VelocityUnits = 1, BFieldUnits = 1;


     if(ComovingCoordinates){
        if (MHDCosmologyGetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
   &TimeUnits, &VelocityUnits, Time,&BFieldUnits) == FAIL) {
     fprintf(stderr, "Error in MHD CosmologyGetUnits.\n");
      return FAIL;
    }
   

     /* Transform speed of light, hydrogen mass and electron change in units of ENZO */
    electroncharge *= TimeUnits*BFieldUnits/(speedoflight*DensityUnits*pow(LengthUnits,3));
    speedoflight /=VelocityUnits;
    hydrogenmass /= DensityUnits*pow(LengthUnits,3);
    }//ComovingCoordinates

#endif //BIERMANN


    /* Create a cell width array to pass (and convert to absolute coords). */

    float *CellWidthTemp[MAX_DIMENSION];
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      CellWidthTemp[dim] = new float[GridDimension[dim]];
      for (i = 0; i < GridDimension[dim]; i++)
        if (dim < GridRank)
          CellWidthTemp[dim][i] = float(a[2]*CellWidth[dim][i]);
        else
          CellWidthTemp[dim][i] = CellWidth[dim][i]; 
      //dcc I changed "else{CellWidthTemp[dim][i] = 1.0" to the above.
    }

#ifdef HAOXU //multispecies
   float  **ColorFlux =  new float*[15], **Color= new float*[15];
//   for(i=0;i<15;i++) {
//     ColorFlux[i] = NULL;
//     Color[i] = NULL;
//     }
   i=0;
   
  if (MultiSpecies>0) {
      Color[i] = BaryonField[ColourNum];
      ColorFlux[i++] = FluxFromSolver[ColourNum];
      Color[i] = BaryonField[ColourNum+1];
      ColorFlux[i++] = FluxFromSolver[ColourNum+1];
      Color[i] = BaryonField[ColourNum+2];
      ColorFlux[i++] = FluxFromSolver[ColourNum+2];
      Color[i] = BaryonField[ColourNum+3];
      ColorFlux[i++] = FluxFromSolver[ColourNum+3];
      Color[i] = BaryonField[ColourNum+4];
      ColorFlux[i++] = FluxFromSolver[ColourNum+4];
      Color[i] = BaryonField[ColourNum+5];
      ColorFlux[i++] = FluxFromSolver[ColourNum+5];
  if (MultiSpecies > 1) {
      Color[i] = BaryonField[ColourNum+6];
      ColorFlux[i++] = FluxFromSolver[ColourNum+6];
      Color[i] = BaryonField[ColourNum+7];
      ColorFlux[i++] = FluxFromSolver[ColourNum+7];
      Color[i] = BaryonField[ColourNum+8];
      ColorFlux[i++] = FluxFromSolver[ColourNum+8];
      }
  if (MultiSpecies > 2) {
      Color[i] = BaryonField[ColourNum+9];
      ColorFlux[i++] = FluxFromSolver[ColourNum+9];
      Color[i] = BaryonField[ColourNum+10];
      ColorFlux[i++] = FluxFromSolver[ColourNum+10];
      Color[i] = BaryonField[ColourNum+11];
      ColorFlux[i++] = FluxFromSolver[ColourNum+11];
      }
    }
    if ((MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields)) != -1) {
      Color[i] = BaryonField[MetalNum];
      ColorFlux[i++] = FluxFromSolver[MetalNum];
      Color[i] = BaryonField[MetalNum+1];
      ColorFlux[i++] = FluxFromSolver[MetalNum+1];
      Color[i] = BaryonField[MetalNum+2];
      ColorFlux[i++] = FluxFromSolver[MetalNum+2];
    }
      int NumberofColor = i;

if(MultiSpecies > 0){     
      for(i=NumberofColor;i<15;i++){
      Color[i] = new float[sizeof(BaryonField[DensNum])]; 
      ColorFlux[i] = new float [sizeof(FluxFromSolver[DensNum])];
       }
    }

#endif //HAOXU

  switch( HydroMethod ){

  case MHD_Harten:
#ifdef HAOXU
    if(DualEnergyFormalism==1){
      FORTRAN_NAME(mhd_harten_cosmology)
	(BaryonField[DensNum], BaryonField[TENum], 
	 BaryonField[Vel1Num],  BaryonField[Vel2Num], BaryonField[Vel3Num], 
	 CenteredB[0], CenteredB[1], CenteredB[2], 
	 &GravityOn,
	 AccelerationField[0],AccelerationField[1],AccelerationField[2],
	 MagneticFlux[0][0], MagneticFlux[1][0], MagneticFlux[2][0],
	 MagneticFlux[0][1], MagneticFlux[1][1], MagneticFlux[2][1],
	 FluxFromSolver[DensNum], FluxFromSolver[TENum],
	 FluxFromSolver[Vel1Num],
	 FluxFromSolver[Vel2Num],FluxFromSolver[Vel3Num],
	 FluxExtents, &TotalFluxSize, &NumberOfSubgrids,
	 CellWidthTemp[0], CellWidthTemp[1], CellWidthTemp[2],
	 GridDimension, GridDimension +1, GridDimension +2, 
	 GridStartIndex, GridEndIndex,
	 GridStartIndex+1, GridEndIndex+1,
	 GridStartIndex+2, GridEndIndex+2,
	 &dtFixed, &Gamma, &CycleNumber, &GridRank, &level, &grid,
	 &MHD_Hack2d, a,  //cosmology expandsion factor added, hx
	 &DualEnergyFormalism, BaryonField[GENum],&PPMDiffusionParameter,  //dual energy,flag for diffusion hx
	 &tiny_pressure);
    }else{
      FORTRAN_NAME(mhd_harten_cosmology)
	(BaryonField[DensNum], BaryonField[TENum], 
	 BaryonField[Vel1Num],  BaryonField[Vel2Num], BaryonField[Vel3Num], 
	 CenteredB[0], CenteredB[1], CenteredB[2], 
	 &GravityOn,
	 AccelerationField[0],AccelerationField[1],AccelerationField[2],
	 MagneticFlux[0][0], MagneticFlux[1][0], MagneticFlux[2][0],
	 MagneticFlux[0][1], MagneticFlux[1][1], MagneticFlux[2][1],
	 FluxFromSolver[DensNum], FluxFromSolver[TENum],
	 FluxFromSolver[Vel1Num],
	 FluxFromSolver[Vel2Num],FluxFromSolver[Vel3Num],
	 FluxExtents, &TotalFluxSize, &NumberOfSubgrids,
	 CellWidthTemp[0], CellWidthTemp[1], CellWidthTemp[2],
	 GridDimension, GridDimension +1, GridDimension +2, 
	 GridStartIndex, GridEndIndex,
	 GridStartIndex+1, GridEndIndex+1,
	 GridStartIndex+2, GridEndIndex+2,
	 &dtFixed, &Gamma, &CycleNumber, &GridRank, &level, &grid,
	 &MHD_Hack2d, a,  //cosmology expandsion factor added, hx
	 &DualEnergyFormalism, 0,&PPMDiffusionParameter, //dual energy,flag for diffusion hx
	 &tiny_pressure);
    }
#else /* HAOXU */
    
    FORTRAN_NAME(mhd_harten)
      (BaryonField[DensNum], BaryonField[TENum], 
       BaryonField[Vel1Num],  BaryonField[Vel2Num], BaryonField[Vel3Num], 
       CenteredB[0], CenteredB[1], CenteredB[2], 
       &GravityOn,
       AccelerationField[0],AccelerationField[1],AccelerationField[2],
       MagneticFlux[0][0], MagneticFlux[1][0], MagneticFlux[2][0],
       MagneticFlux[0][1], MagneticFlux[1][1], MagneticFlux[2][1],
       FluxFromSolver[DensNum], FluxFromSolver[TENum],
       FluxFromSolver[Vel1Num],
       FluxFromSolver[Vel2Num],FluxFromSolver[Vel3Num],
       FluxExtents, &TotalFluxSize, &NumberOfSubgrids,
       CellWidth[0], CellWidth[1], CellWidth[2],
       GridDimension, GridDimension +1, GridDimension +2, 
       GridStartIndex, GridEndIndex,
       GridStartIndex+1, GridEndIndex+1,
       GridStartIndex+2, GridEndIndex+2,
       &dtFixed, &Gamma, &CycleNumber, &GridRank, &level, &grid,
       &MHD_Hack2d
#ifdef NSS
       , MHD_Recon, MHD_Riemann, MHD_Eps
#endif
       );
    
#endif /* HAOXU */
    break;    
    
#ifdef HAOXU
  case 6:

    if(DualEnergyFormalism==1){
      if(NumberofColor==0){
	FORTRAN_NAME(mhd_li)
	  (BaryonField[DensNum], BaryonField[TENum],
	   BaryonField[Vel1Num],  BaryonField[Vel2Num], BaryonField[Vel3Num],
	   CenteredB[0], CenteredB[1], CenteredB[2],
	   &GravityOn,
	   AccelerationField[0],AccelerationField[1],AccelerationField[2],
	   MagneticFlux[0][0], MagneticFlux[1][0], MagneticFlux[2][0],
	   MagneticFlux[0][1], MagneticFlux[1][1], MagneticFlux[2][1],
	   FluxFromSolver[DensNum], FluxFromSolver[TENum],
	   FluxFromSolver[Vel1Num],
	   FluxFromSolver[Vel2Num],FluxFromSolver[Vel3Num],
	   FluxFromSolver[GENum],
	   FluxExtents, &TotalFluxSize, &NumberOfSubgrids,
	   CellWidthTemp[0], CellWidthTemp[1], CellWidthTemp[2],
	   GridDimension, GridDimension +1, GridDimension +2,
	   GridStartIndex, GridEndIndex,
	   GridStartIndex+1, GridEndIndex+1,
	   GridStartIndex+2, GridEndIndex+2,
	   &dtFixed, &Gamma, &CycleNumber, &GridRank, &level, &grid,
	   a,  //cosmology expandsion factor added, hx
	   &DualEnergyFormalism, BaryonField[GENum],&PPMDiffusionParameter, //dual energy,flag for diffusion hx
	   &tiny_pressure, &MHD_Equation, MHDLi, &EquationOfState, &IsothermalSoundSpeed);
      }else{
	FORTRAN_NAME(mhd_li_ms)
	  (BaryonField[DensNum], BaryonField[TENum],
	   BaryonField[Vel1Num],  BaryonField[Vel2Num], BaryonField[Vel3Num],
	   CenteredB[0], CenteredB[1], CenteredB[2],
	   &GravityOn,
	   AccelerationField[0],AccelerationField[1],AccelerationField[2],
	   MagneticFlux[0][0], MagneticFlux[1][0], MagneticFlux[2][0],
	   MagneticFlux[0][1], MagneticFlux[1][1], MagneticFlux[2][1],
	   FluxFromSolver[DensNum], FluxFromSolver[TENum],
	   FluxFromSolver[Vel1Num],
	   FluxFromSolver[Vel2Num],FluxFromSolver[Vel3Num],
	   FluxFromSolver[GENum],
	   FluxExtents, &TotalFluxSize, &NumberOfSubgrids,
	   CellWidthTemp[0], CellWidthTemp[1], CellWidthTemp[2],
	   GridDimension, GridDimension +1, GridDimension +2,
	   GridStartIndex, GridEndIndex,
	   GridStartIndex+1, GridEndIndex+1,
	   GridStartIndex+2, GridEndIndex+2,
	   &dtFixed, &Gamma, &CycleNumber, &GridRank, &level, &grid,
	   a,  //cosmology expandsion factor added, hx
	   &DualEnergyFormalism, BaryonField[GENum],&PPMDiffusionParameter, //dual energy,flag for diffusion hx
	   &tiny_pressure, &MHD_Equation, MHDLi, 
	   &NumberofColor, Color[0],Color[1],Color[2],Color[3],Color[4],Color[5],
	   Color[6],Color[7],Color[8],Color[9],Color[10],Color[11],Color[12],Color[13],Color[14],
	   ColorFlux[0],ColorFlux[1],ColorFlux[2],ColorFlux[3],ColorFlux[4],ColorFlux[5],
	   ColorFlux[6],ColorFlux[7],ColorFlux[8],ColorFlux[9],ColorFlux[10],ColorFlux[11],ColorFlux[12],ColorFlux[13],ColorFlux[14]);
      }    
    }else{
      float *Pointer_GE, *Flux_GE;
      Pointer_GE = new float[(*GridDimension)*(*(GridDimension+1))*(*(GridDimension+2))]; 
      Flux_GE = new float[TotalFluxSize];
      if(NumberofColor==0){
	FORTRAN_NAME(mhd_li)
	  (BaryonField[DensNum], BaryonField[TENum],
	   BaryonField[Vel1Num],  BaryonField[Vel2Num], BaryonField[Vel3Num],
	   CenteredB[0], CenteredB[1], CenteredB[2],
	   &GravityOn,
	   AccelerationField[0],AccelerationField[1],AccelerationField[2],
	   MagneticFlux[0][0], MagneticFlux[1][0], MagneticFlux[2][0],
	   MagneticFlux[0][1], MagneticFlux[1][1], MagneticFlux[2][1],
	   FluxFromSolver[DensNum], FluxFromSolver[TENum],
	   FluxFromSolver[Vel1Num],
	   FluxFromSolver[Vel2Num],FluxFromSolver[Vel3Num],
	   FluxFromSolver[GENum],
	   FluxExtents, &TotalFluxSize, &NumberOfSubgrids,
	   CellWidthTemp[0], CellWidthTemp[1], CellWidthTemp[2],
	   GridDimension, GridDimension +1, GridDimension +2,
	   GridStartIndex, GridEndIndex,
	   GridStartIndex+1, GridEndIndex+1,
	   GridStartIndex+2, GridEndIndex+2,
	   &dtFixed, &Gamma, &CycleNumber, &GridRank, &level, &grid,
	   a,  //cosmology expandsion factor added, hx
	   &DualEnergyFormalism, Pointer_GE,&PPMDiffusionParameter, //dual energy,flag for diffusion hx
	   &tiny_pressure, &MHD_Equation, MHDLi, &EquationOfState, &IsothermalSoundSpeed);
      }else{
	FORTRAN_NAME(mhd_li_ms)(BaryonField[DensNum], BaryonField[TENum], 
				BaryonField[Vel1Num],  BaryonField[Vel2Num], BaryonField[Vel3Num],
				CenteredB[0], CenteredB[1], CenteredB[2],
				&GravityOn,
				AccelerationField[0],AccelerationField[1],AccelerationField[2],
				MagneticFlux[0][0], MagneticFlux[1][0], MagneticFlux[2][0],
				MagneticFlux[0][1], MagneticFlux[1][1], MagneticFlux[2][1],
				FluxFromSolver[DensNum], FluxFromSolver[TENum],
				FluxFromSolver[Vel1Num],
				FluxFromSolver[Vel2Num],FluxFromSolver[Vel3Num],
				FluxFromSolver[GENum],
				FluxExtents, &TotalFluxSize, &NumberOfSubgrids,
				CellWidthTemp[0], CellWidthTemp[1], CellWidthTemp[2],
				GridDimension, GridDimension +1, GridDimension +2,
				GridStartIndex, GridEndIndex,
				GridStartIndex+1, GridEndIndex+1,
				GridStartIndex+2, GridEndIndex+2,
				&dtFixed, &Gamma, &CycleNumber, &GridRank, &level, &grid,
				a,  //cosmology expandsion factor added, hx
				&DualEnergyFormalism, Pointer_GE,&PPMDiffusionParameter, //dual energy,flag for diffusion hx
				&tiny_pressure, &MHD_Equation, MHDLi,   
				&NumberofColor, Color[0],Color[1],Color[2],Color[3],Color[4],Color[5],                         
				Color[6],Color[7],Color[8],Color[9],Color[10],Color[11],Color[12],Color[13],Color[14],
				ColorFlux[0],ColorFlux[1],ColorFlux[2],ColorFlux[3],ColorFlux[4],ColorFlux[5],
				ColorFlux[6],ColorFlux[7],ColorFlux[8],ColorFlux[9],ColorFlux[10],ColorFlux[11],ColorFlux[12],ColorFlux[13],ColorFlux[14]); 
      }
      delete Pointer_GE;     
      delete Flux_GE;
    }//DualEnergyFormalism
    break;
#endif //HAOXU
    
  case MHD_None:
    fprintf(stderr,"=============== NASTY KLUDGE!!! NO SOLVER!!! ==================\n");
    break;
    
  default:
    if(MyProcessorNumber == ROOT_PROCESSOR) 
      fprintf(stderr, "SolveMHDEquations:  Hydro Method is not a defined MHD Method\n");
    return FAIL;
    
  }

  if( CheckForSingleGridDump(31) == TRUE ){
    float * SavedMagneticField[3];
    for(i=0;i<3;i++) SavedMagneticField[i] = MagneticField[i];
    for(i=0;i<3;i++) MagneticField[i] = MagneticFlux[i][0];
    
    sprintf(basename, "data31%d%d.grid",CycleNumber, level);
    FILE *dummy = fopen(basename, "a");   

    if( this->WriteGrid(dummy, basename, MyGridNumber) == FAIL ){
      fprintf(stderr, "Shit.  Problem with Write Grid in SMHD.\n");
      return FAIL;
    }
    fclose(dummy);
    for(i=0;i<3;i++) MagneticField[i] = SavedMagneticField[i];

  }  
#ifdef NEW_DIVB
  AthenaObj ATH(this);
  float* Fluxes[3] = {MagneticFlux[0][0],MagneticFlux[1][0],MagneticFlux[2][0]};

#endif //NEW_DIVB  
  int CurlStart[3] = {0,0,0}, 
    CurlEnd[3] = {GridDimension[0]-1,GridDimension[1]-1,GridDimension[2]-1};
  if( HydroMethod != MHD_None )
    switch( MHD_DivB ){
      

#ifdef NEW_DIVB
    case MHD_DivB_Athena_Balsara: //4 this is the ppml naming convention.
    case MHD_DivB_Athena_LF: //2
    case MHD_DivB_Athena_Switch: //3

      ATH.ComputeElectricField(dtFixed, Fluxes);
      //<dbg> test, db21
      //MHD_Curl( GridStartIndex, GridEndIndex, 1);

      MHD_Curl( CurlStart,CurlEnd, 1);
      //</dbg>
      CenterMagneticField();
      break;
#endif //NEW_DIVB
    case MHD_DivB_Balsara:
      //ok, this is kind of a kludge.  I decided to use dt*E instead of E in the ElectricField.
      //This is due to the timestep method in enzo, and the flux correction module.  However, it wasn't originally
      //written that way, so things go sort of ugly.

      dtUsed = 1.0;
      UseDT = 1;

      //dtUsed = dtFixed
      
#ifdef BIERMANN
      
      FORTRAN_NAME(create_e_biermann)(MagneticFlux[0][0], MagneticFlux[1][0], MagneticFlux[2][0],
				      MagneticFlux[0][1], MagneticFlux[1][1], MagneticFlux[2][1],
				      ElectricField[0], ElectricField[1], ElectricField[2],
				      CellWidthTemp[0], CellWidthTemp[1], CellWidthTemp[2],
				      BaryonField[DensNum],BaryonField[GENum],
				      GridDimension, GridDimension +1, GridDimension +2,
				      GridStartIndex, GridEndIndex,
				      GridStartIndex+1, GridEndIndex+1,
				      GridStartIndex+2, GridEndIndex+2, &dtFixed, &UseDT,
				      &Gamma, &speedoflight,&hydrogenmass, &electroncharge, &chi, a);
      
      
#else //BIERMANN
      
      FORTRAN_NAME(create_e)(MagneticFlux[0][0], MagneticFlux[1][0], MagneticFlux[2][0],
			     MagneticFlux[0][1], MagneticFlux[1][1], MagneticFlux[2][1],
			     ElectricField[0], ElectricField[1], ElectricField[2],
			     GridDimension, GridDimension +1, GridDimension +2,
			     GridStartIndex, GridEndIndex,
			     GridStartIndex+1, GridEndIndex+1,
			     GridStartIndex+2, GridEndIndex+2, &dtFixed, &UseDT);
      

#endif //BIERMANN


      
      FORTRAN_NAME(curl_of_e)(MagneticField[0], MagneticField[1], MagneticField[2],
			      ElectricField[0], ElectricField[1], ElectricField[2],
			      CellWidth[0], CellWidth[1], CellWidth[2],
			      GridDimension, GridDimension +1, GridDimension +2,
			      GridStartIndex, GridEndIndex,
			      GridStartIndex+1, GridEndIndex+1,
			      GridStartIndex+2, GridEndIndex+2,
			      &dtUsed, &MHD_CenteringMethod);


#ifdef OLD_CENTER
      
      //With the new centering paradigm, it's done at the beginning of SetBoundaryConditions
      //so doing it here is redundant.
      if( this->CenterMagneticField() == FAIL ) {
	fprintf(stderr," error with CenterMagneticField\n");
	return FAIL;
      }
#endif //OLD_CENTER

	break;
	
    case MHD_DivB_RJ:
      
      if( DEFAULT_GHOST_ZONES != 2 ){
	fprintf(stderr, "This routine expects exactly 2 ghost zones.  You have %d.  Fix it.",DEFAULT_GHOST_ZONES);
	return FAIL;
      }
      
      FORTRAN_NAME(divb_rj)(MagneticFlux[0][0], MagneticFlux[1][0], MagneticFlux[2][0],
			    MagneticFlux[0][1], MagneticFlux[1][1], MagneticFlux[2][1],
			    MagneticField[0], MagneticField[1], MagneticField[2],
			    CenteredB[0], CenteredB[1],CenteredB[2],
			    CellWidth[0], CellWidth[1], CellWidth[2], &dtFixed,
			    &nx, &ny, &nz);
      
      break;
      
    case MHD_DivB_Poisson:
      fprintf(stderr, "Shit!  You need to instal the hodge projection routine.\n");
      return FAIL;
      break;

    case MHD_DivB_none:
    default:
      
      if(MyProcessorNumber == ROOT_PROCESSOR )
	fprintf(stderr, "Warning: No Div B = 0 method used\n");
      break;
      
    }//End Switch
  
  /*  
  if( CheckForSingleGridDump(39) == TRUE){
    int writetmp = MHD_WriteElectric;
    MHD_WriteElectric=TRUE;
    sprintf(basename, "data39%02d%d.grid",CycleNumber, level);
    FILE *dummy = fopen(basename, "a");    
    if( this->WriteGrid(dummy, basename, MyGridNumber) == FAIL ){
      fprintf(stderr, "Shit.  Problem with Write Grid in SMHD.\n");
      return FAIL;
    }
    fclose(dummy);
    MHD_WriteElectric=writetmp;
  }  

   */
  //

  // Fill SubgridFluxes
  //
  
  
  for( field=0;field<NumberOfBaryonFields; field++)
    for( subgrid=0;subgrid<NumberOfSubgrids; subgrid++)
      for(dim=0;dim<3;dim++)
	for(k=0;k<FluxDims[dim][2][subgrid];k++)
	  for(j=0;j<FluxDims[dim][1][subgrid];j++)
	    for(i=0;i<FluxDims[dim][0][subgrid];i++){
	      SubgridFluxes[subgrid]->LeftFluxes[field][dim]
		[i+FluxDims[dim][0][subgrid]*(j+k*FluxDims[dim][1][subgrid])]=
		
		FluxFromSolver[field][	i+
					FluxDims[dim][0][subgrid]*(j
				       +FluxDims[dim][1][subgrid]*(k
				       +FluxDims[dim][2][subgrid]*(0)))
				       +TotalOffset[dim][subgrid]];
	      
	      
	      SubgridFluxes[subgrid]->RightFluxes[field][dim]
		[i+FluxDims[dim][0][subgrid]*(j+k*FluxDims[dim][1][subgrid])]=

		FluxFromSolver[field][	i+
					FluxDims[dim][0][subgrid]*(j
				       +FluxDims[dim][1][subgrid]*(k
				       +FluxDims[dim][2][subgrid]*(1)))
					+TotalOffset[dim][subgrid]];
	      
	    }

  JBMEM_MESSAGE(MyProcessorNumber,"jb: SMHDEFlux0");  
  for(subgrid=0;subgrid<NumberOfSubgrids; subgrid++)
    for(dim=0;dim<3;dim++)
      for( field=0;field<3;field++){
	
	SubgridFluxes[subgrid]->RightElectric[field][dim]=NULL;
	SubgridFluxes[subgrid]->LeftElectric[field][dim]=NULL;
	/* this is for debugging.
	for(i=0;i<ElectricDims[field][0];i++)
	  for(j=0;j<ElectricDims[field][1];j++)
	    for(k=0;k<ElectricDims[field][2];k++){
	      if(level == -1 && field == 0 )
		ElectricField[field][i+ElectricDims[field][0]*(j+ElectricDims[field][1]*k)]
		  = 10;
	      else
		ElectricField[field][i+ElectricDims[field][0]*(j+ElectricDims[field][1]*k)]
		  = 0.0;
	      
	    }
	*/  
	//JBMEM_MESSAGE(MyProcessorNumber,"jb: SMHDEFlux1");  
	size=1;
	if( field != dim){

	  for( coord=0;coord<3;coord++){
	    Dim[coord]=FluxDims[dim][coord][subgrid] ;
	    if(coord != field && coord != dim )
	      Dim[coord]++;
	    size *= Dim[coord];
	  }
	  //fprintf(stderr, "sgfeea: level %d sub %d dim %d field %d size %d\n",
	  //level, subgrid, dim, field, size);
	  SubgridFluxes[subgrid]->LeftElectric[field][dim]=new float[size];
	  SubgridFluxes[subgrid]->RightElectric[field][dim]=new float[size];
	  
	  for(i=0;i<size;i++){
	    SubgridFluxes[subgrid]->LeftElectric[field][dim][i]=0.0;
	    SubgridFluxes[subgrid]->RightElectric[field][dim][i] = 0.0;
	  }
	  //Global memory leak here, too.  
	  //JBMEM_MESSAGE(MyProcessorNumber,"jb: SMHDEFlux2");  
	  
	  //the field along this dim isn't used.
	  
	  
	  //for your reference:
	  //fluxe(dim, coord, left or right, start or end, subgrid)
	  /*
	    fprintf(stderr, "sgfee2: level %d fluxes(%d, 2, 0, 0, %d ) %d %d\n",
	    level, dim, subgrid, fluxe(dim,2,0,0,subgrid), fluxe(dim,2,0,1,subgrid));
	    fprintf(stderr, "sgfee2: level %d fluxes(%d, 1, 0, 0, %d ) %d %d\n",
	    level, dim, subgrid, fluxe(dim,1,0,0,subgrid), fluxe(dim,1,0,1,subgrid));
	    fprintf(stderr, "sgfee2: level %d fluxes(%d, 0, 0, 0, %d ) %d %d\n",
	    level, dim, subgrid, fluxe(dim,0,0,0,subgrid), fluxe(dim,0,0,1,subgrid));
	  */
	  for(k=0;k<Dim[2];k++)
	    for(j=0;j<Dim[1];j++)
	      for(i=0;i<Dim[0];i++){
		Findex=i+Dim[0]*(j+Dim[1]*k);
		
		Cindex=(i+fluxe(dim,0,0,0,subgrid)-1)
		  +ElectricDims[field][0]*(j+fluxe(dim,1,0,0,subgrid)-1
					   +ElectricDims[field][1]*(k+fluxe(dim,2,0,0,subgrid)-1));
		
		SubgridFluxes[subgrid]->LeftElectric[field][dim][Findex]= 
		  ElectricField[field][Cindex];//*dtFixed; taken care of in the definition of ElectricField
		
		Cindex=(i+fluxe(dim,0,1,0,subgrid)-1 )
		  +ElectricDims[field][0]*(j+fluxe(dim,1,1,0,subgrid)-1
					   +ElectricDims[field][1]*(k+fluxe(dim,2,1,0,subgrid)-1));
		
		SubgridFluxes[subgrid]->RightElectric[field][dim][Findex]=
		  ElectricField[field][Cindex];// *dtFixed; taken care of in the definition of E.
		
		/*
		  fprintf(stderr,"sgfee: index %d Right, Left %f %f\n", 
		  Findex,
		  SubgridFluxes[subgrid]->LeftElectric[field][dim][Findex], 
		  SubgridFluxes[subgrid]->RightElectric[field][dim][Findex]);
		*/
	      }
	}//field==dim
	
}//E field



//if( this->MHDAnis(" SMHD: End ") == FAIL ) 
//return FAIL;
  JBMEM_MESSAGE(MyProcessorNumber,"jb: SMHDBeforeMagFluxDelete");  
#ifdef NEW_DIVB
  for(field=0;field<3;field++)
    delete [] MagneticFlux[field][0];

#else //NEW_DIVB
  for(field=0;field<3;field++)
    for(axis=0;axis<2;axis++){
      delete [] MagneticFlux[field][axis]; 
    }
#endif //NEW_DIVB 
  JBMEM_MESSAGE(MyProcessorNumber,"jb: SMHDAfterLocalFluxDelete");  
  for(field=0;field<NumberOfBaryonFields; field++){
    delete FluxFromSolver[field];
  }
  JBMEM_MESSAGE(MyProcessorNumber,"jb: SMHDAfterMagFluxDelete");  


#ifdef HAOXU

if(MultiSpecies>0){
 for(field=0;field<NumberofColor;field++) {
     ColorFlux[field] = NULL;
     Color[field] = NULL;
     }
 for(field=NumberofColor;field<15;field++){
     delete ColorFlux[field];
     delete Color[field];
     }
   }
  delete [] ColorFlux, Color;
#endif
  //IsItShit("smhd4");
  
  this->DebugCheck("SMHD: After");

  if( this->CheckForNans("end of SMHD") == FAIL ) 
	  return FAIL;

  if( TVtool("end of SMHD") == FAIL ){
    fprintf(stderr,"problem at the end of SMHD.\n");
    return FAIL;
  }
  wall_time("End SMHD");

  if( CheckForSingleGridDump(39) == TRUE){
    int writetmp = MHD_WriteElectric;
    MHD_WriteElectric=TRUE;
    sprintf(basename, "data39%02d%d.grid",CycleNumber, level);
    FILE *dummy = fopen(basename, "a");    
    if( this->WriteGrid(dummy, basename, MyGridNumber) == FAIL ){
      fprintf(stderr, "Shit.  Problem with Write Grid in SMHD.\n");
      return FAIL;
    }
    fclose(dummy);
    MHD_WriteElectric=writetmp;
  }  


  return SUCCESS;
}
