/***********************************************************************
/
/  GRID CLASS (SOLVES THE MHD EQUATIONS)
/
/  written by: David Collins
/  date:       2004-2013
/  modified1:
/
/  PURPOSE:   Fills and allocates fluxes and updates the MHD equations
/             Depricated in favor of the new C wrapper, which is now 
/             called form SolveHydroEquations.
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#include "ErrorExceptions.h"
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
#include "DebugTools.h"
#include "CosmologyParameters.h"

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

int FindField(int f, int farray[], int n);

// To call MHD solver of Shengtai
extern "C" void FORTRAN_NAME( mhd_li)
  (float *d, float *e, float *vx, float *vy, float *vz,
   float *bxc, float *byc, float *bzc,
   int *gravityon, float *gr_ax, float *gr_ay, float *gr_az,
   float *fx1, float *fy1, float *fz1,
   float *fx2, float *fy2, float *fz2,
   float *fd, float *fe, float *fvx, float *fvy, float *fvz, float *fge,
   int *fluxextents, int *totalfluxsize, int *nsubgrids,
   FLOAT *dx, FLOAT *dy, FLOAT *dz, int *idim, int *jdim, int *kdim,
   int *i1, int *i2, int *j1, int *j2, int *k1, int *k2, float *dt, float *gamma,
   int *nhy, int *rank, int *level, int *grid,  
   FLOAT  *a, int * comoving, 
   int *idual, float *ge, int *idiff,
   float *premin, 
  int *MHDCTSlopeLimiter, int *ReconstructionMethod, int *RiemannSolver,int *MHDCTDualEnergyMethod,int *MHDCTPowellSource,   
   int *EquationOfState, float* IsothermalSoundSpeed,
   int * hack);
//multi-species

extern "C" void FORTRAN_NAME( mhd_li_ms)
  (float *d, float *e, float *vx, float *vy, float *vz,
   float *bxc, float *byc, float *bzc,
   int *gravityon, float *gr_ax, float *gr_ay, float *gr_az,
   float *fx1, float *fy1, float *fz1,
   float *fx2, float *fy2, float *fz2,
   float *fd, float *fe, float *fvx, float *fvy, float *fvz, float *fge,
   int *fluxextents, int *totalfluxsize, int *nsubgrids,
   FLOAT *dx, FLOAT *dy, FLOAT *dz, int *idim, int *jdim, int *kdim,
   int *i1, int *i2, int *j1, int *j2, int *k1, int *k2, float *dt, float *gamma,
   int *nhy, int *rank, int *level, int *grid,
   FLOAT  *a, int *idual, float *ge, int *idiff,
   float *premin, 
  int *MHDCTSlopeLimiter, int *ReconstructionMethod, int *RiemannSolver,int *MHDCTDualEnergyMethod,int *MHDCTPowellSource,   
   int *numberofcolor,
   float *c0,float *c1,float *c2, float *c3, float *c4,float *c5,
   float *c6,float *c7,float *c8, float *c9, float *c10, float *c11,
   float *c12, float *c13, float *c14,
   float *fc0,float *fc1,float *fc2, float *fc3, float *fc4,float *fc5,
   float *fc6,float *fc7,float *fc8, float *fc9, float *fc10, float *fc11,
   float *fc12, float *fc13, float *fc14); 


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
int hack = FALSE;
int grid::SolveMHDEquations(int CycleNumber, int NumberOfSubgrids, 
			    fluxes *SubgridFluxes[], int level, int grid)
/*begin*/
{
  if( ProcessorNumber != MyProcessorNumber )
    return SUCCESS;
  
  fprintf(stderr,"===  SMHD n = %"ISYM" L = %"ISYM" g = %"ISYM" proc = %"ISYM" dt = %15.12e id %"ISYM" === Right (%"ESYM", %"ESYM", %"ESYM")\n",
	  CycleNumber, level, grid, MyProcessorNumber, dtFixed, GetGridID(), GridRightEdge[0], GridRightEdge[1], GridRightEdge[2]);

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  int i, j, k, n, dim, field, axis, face, index, z, coord, Findex, Cindex, index2, size;
  int subgrid, TotalFluxSize=0, SizeOfFace=1, SizeOfSubgrid = 0;
  int *SizeOtherSubgrids = new int [NumberOfSubgrids];
  int *TotalOffset[3], *SizeOtherDims[3], Dim[3]={1,1,1};
  
  int nx = GridDimension[0] - 2*NumberOfGhostZones, 
    ny = GridDimension[1] - 2*NumberOfGhostZones, 
    nz = GridDimension[2] - 2*NumberOfGhostZones;
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
  MHDCT_ConvertEnergyToConservedS();

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

    if(TestProblemData.MultiMetals || MultiMetals){
      coloff[NumberOfColours++] = BaryonField[MetalNum+1] - colourpt;
      coloff[NumberOfColours++] = BaryonField[MetalNum+2] - colourpt;
    }    

  }
  
  
  // The fluxe macro is used to map the 5 dimensional, non-rectangular array to a 1 dimensional array
  // that fortran can deal with. 
  
#define fluxe(dim,coord,face,end,sub) FluxExtents[dim+3*(coord+3*(face+2*(end+2*sub)))]
  int * FluxExtents = new int[3*3*2*2*NumberOfSubgrids];
  int * FluxDims[3][3];
  float ** FluxFromSolver = new float*[NumberOfBaryonFields];
  
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
  
  for(field=0;field<NumberOfBaryonFields;field++){
    FluxFromSolver[field] = new float[TotalFluxSize];
    for(i=0;i<TotalFluxSize; i++)
      FluxFromSolver[field][i]=0.0;
  }
  
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
      fprintf(stderr, "========== Solve MHD create Magnetic Field: Error with MagneticField\n");
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
  
  //each magnetic field has components from the flux of two 'other' axis
  
  for(field=0;field<3;field++){
    MagneticFlux[field][0] = new float[2*MagneticSize[field]];
    MagneticFlux[field][1] =  MagneticFlux[field][0] +MagneticSize[field];
  }
  
  for(field=0;field<3;field++)
    for( i=0;i<MagneticSize[field]; i++){
      MagneticFlux[field][0][i] = 0.0;
      MagneticFlux[field][1][i] = 0.0;
    }
  
  if( this->CenterMagneticField() == FAIL ) {
    fprintf(stderr," error with CenterMagneticField, first call \n");
    return FAIL;
  }
  
  
  /* Set minimum support. */
  
  float MinimumSupportEnergyCoefficient = 0;
  if (UseMinimumPressureSupport == TRUE && level>MaximumRefinementLevel-1)
    if (this->SetMinimumSupport(MinimumSupportEnergyCoefficient) == FAIL) {
      fprintf(stderr, "Error in grid->SetMinimumSupport,\n");
      return FAIL;
    }
  
  
  int GravityOn = 0;
  if (SelfGravity || UniformGravity || PointSourceGravity || ExternalGravity)
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
  
  FLOAT *CellWidthTemp[MAX_DIMENSION];
  for (dim = 0; dim < MAX_DIMENSION; dim++) {
    CellWidthTemp[dim] = new FLOAT[GridDimension[dim]];
    for (i = 0; i < GridDimension[dim]; i++)
      if (dim < GridRank)
	CellWidthTemp[dim][i] = FLOAT(a[2]*CellWidth[dim][i]);
      else
	CellWidthTemp[dim][i] = CellWidth[dim][i]; 
  }

  float  **ColorFlux =  new float*[15], **Color= new float*[15];

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

    if(TestProblemData.MultiMetals || MultiMetals){
      Color[i] = BaryonField[MetalNum+1];
      ColorFlux[i++] = FluxFromSolver[MetalNum+1];
      Color[i] = BaryonField[MetalNum+2];
      ColorFlux[i++] = FluxFromSolver[MetalNum+2];
    }
  }
  int NumberofColor = i;
  
  if(MultiSpecies > 0){     
    for(i=NumberofColor;i<15;i++){
      Color[i] = new float[sizeof(BaryonField[DensNum])]; 
      ColorFlux[i] = new float [sizeof(FluxFromSolver[DensNum])];
    }
  }

  float *Pointer_GE;
  switch( HydroMethod ){
    
  case MHD_Li:
    
    if(DualEnergyFormalism==1){
        Pointer_GE = BaryonField[GENum];
    }else /*dual energy*/{
      Pointer_GE = new float[(*GridDimension)*(*(GridDimension+1))*(*(GridDimension+2))]; 
    }
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
       &ComovingCoordinates,
	   &DualEnergyFormalism, Pointer_GE,&PPMDiffusionParameter, //dual energy,flag for diffusion hx
	   &tiny_pressure, 
       &MHDCTSlopeLimiter, &ReconstructionMethod, &RiemannSolver ,&MHDCTDualEnergyMethod, &MHDCTPowellSource,   
       &EquationOfState, &IsothermalSoundSpeed,
       &hack);
      }else /*multi species */{

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
				&tiny_pressure, 
       &MHDCTSlopeLimiter, &ReconstructionMethod, &RiemannSolver ,&MHDCTDualEnergyMethod, &MHDCTPowellSource ,
				&NumberofColor, Color[0],Color[1],Color[2],Color[3],Color[4],Color[5],                         
				Color[6],Color[7],Color[8],Color[9],Color[10],Color[11],Color[12],Color[13],Color[14],
				ColorFlux[0],ColorFlux[1],ColorFlux[2],ColorFlux[3],ColorFlux[4],ColorFlux[5],
				ColorFlux[6],ColorFlux[7],ColorFlux[8],ColorFlux[9],ColorFlux[10],ColorFlux[11],ColorFlux[12],
                ColorFlux[13],ColorFlux[14] ); 

      }
      if( DualEnergyFormalism == FALSE )
          delete Pointer_GE;     

    break;

  case NoHydro:
    fprintf(stderr,"=============== NASTY KLUDGE!!! NO SOLVER!!! ==================\n");
    break;
    
  default:
    if(MyProcessorNumber == ROOT_PROCESSOR) 
      fprintf(stderr, "SolveMHDEquations:  Hydro Method is not a defined MHD Method\n");
    return FAIL;
    
  }
  

  float* Fluxes[3] = {MagneticFlux[0][0],MagneticFlux[1][0],MagneticFlux[2][0]};

  int CurlStart[3] = {0,0,0}, 
    CurlEnd[3] = {GridDimension[0]-1,GridDimension[1]-1,GridDimension[2]-1};
  if( HydroMethod != NoHydro )
    switch( MHD_CT_Method ){

    case CT_BalsaraSpicer: //0
    case CT_Athena_LF: //1
    case CT_Athena_Switch: //2

      ComputeElectricField(dtFixed, Fluxes);

      MHD_Curl( CurlStart,CurlEnd, 1);

      CenterMagneticField();
      break;
      
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
      
      
#endif //BIERMANN


      if( this->CenterMagneticField() == FAIL ) {
	fprintf(stderr," error with CenterMagneticField\n");
	return FAIL;
      }

	break;
	
    case CT_None:
      break;
    default:
      
      if(MyProcessorNumber == ROOT_PROCESSOR )
	fprintf(stderr, "Warning: No CT method used with MHD_Li.\n");
      break;
      
    }//End Switch

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

  for(field=0;field<3;field++)
    delete [] MagneticFlux[field][0];

  for(field=0;field<NumberOfBaryonFields; field++){
    delete FluxFromSolver[field];
  }



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
  delete [] FluxExtents, FluxDims;

  MHDCT_ConvertEnergyToSpecificS();

  return SUCCESS;
}
