/***********************************************************************
/
/  GRID CLASS (SOLVE THE HYDRO EQUATIONS, SAVING SUBGRID FLUXES)
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:
/
/  PURPOSE:
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

// Solve the hydro equations with the solver, saving the subgrid fluxes
//

#include <stdio.h>
#include "ErrorExceptions.h"
#include "EnzoTiming.h"
#include "performance.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"

/* function prototypes */

int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);
int FindField(int f, int farray[], int n);

#ifdef PPM_LR
extern "C" void FORTRAN_NAME(ppm_lr)(
			  float *d, float *E, float *u, float *v, float *w,
			    float *ge,
                          int *grav, float *gr_ax, float *gr_ay, float *gr_az,
			  float *gamma, float *dt, int *cycle_number,
                            float dx[], float dy[], float dz[],
			  int *rank, int *in, int *jn, int *kn,
                            int is[], int ie[],
			  float gridvel[], int *flatten, int *ipresfree,
			  int *diff, int *steepen, int *idual, 
                            float *eta1, float *eta2,
			  int *num_subgrids, int leftface[], int rightface[],
			  int istart[], int iend[], int jstart[], int jend[],
			  float *standard, int dindex[], int Eindex[],
			  int uindex[], int vindex[], int windex[],
			    int geindex[], float *temp,
                          int *ncolour, float *colourpt, int *coloff,
                            int colindex[], int *ifloat_size);
#endif /* PPM_LR */

int grid::SolveHydroEquations(int CycleNumber, int NumberOfSubgrids, 
			      fluxes *SubgridFluxes[], int level)
{

  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber || !UseHydro)
    return SUCCESS;

  LCAPERF_START("grid_SolveHydroEquations");
  TIMER_START("SolveHydroEquations");

  this->DebugCheck("SolveHydroEquations");

  if (NumberOfBaryonFields > 0) {

    /* initialize */

    // MAX_COLOR is defined in fortran.def
    int dim, i, j, field, size, subgrid, n, colnum[MAX_COLOR];
    Elong_int GridGlobalStart[MAX_DIMENSION];
    FLOAT a = 1, dadt;

    /* Compute size (in floats) of the current grid. */

    size = 1;
    for (dim = 0; dim < GridRank; dim++)
      size *= GridDimension[dim];

    /* If multi-species being used, then treat them as colour variables
       (note: the solver has been modified to treat these as density vars). */

    int NumberOfColours = 0, ColourNum;

    // use different color fields for RadiativeTransferFLD problems
    //   first, the standard Enzo color field advection
    if (MultiSpecies > 0 && RadiativeTransferFLD != 2) {
      NumberOfColours = 6 + 3*(MultiSpecies-1);

      if ((ColourNum =
           FindField(ElectronDensity, FieldType, NumberOfBaryonFields)) < 0) {
        ENZO_FAIL("Could not find ElectronDensity.");
      }

      /* Generate an array of field numbers corresponding to the colour fields
	 (here assumed to start with ElectronDensity and continue in order). */

      for (i = 0; i < NumberOfColours; i++)
        colnum[i] = ColourNum+i;

    }
    // second, the color field advection if using RadiativeTransferFLD for 
    // a radiation propagation problem (i.e. not using ray-tracing)
    if (RadiativeTransferFLD == 2) {
      if (ImplicitProblem < 4)  {  // grey radiation problem
	
	// set the grey radiation field (required)
	if ((ColourNum =
	     FindField(RadiationFreq0, FieldType, NumberOfBaryonFields)) < 0) 
	  ENZO_FAIL("Could not find RadiationFreq0.");
	colnum[0] = ColourNum;

	// check for other chemistry fields; add if they're present
	//   ElectronDensity
	if ((ColourNum =
	     FindField(ElectronDensity, FieldType, NumberOfBaryonFields)) >= 0) 
	  colnum[++NumberOfColours] = ColourNum;
	//   HIDensity
	if ((ColourNum =
	     FindField(HIDensity, FieldType, NumberOfBaryonFields)) >= 0) 
	  colnum[++NumberOfColours] = ColourNum;
	//   HIIDensity
	if ((ColourNum =
	     FindField(HIIDensity, FieldType, NumberOfBaryonFields)) >= 0) 
	  colnum[++NumberOfColours] = ColourNum;
	//   HeIDensity
	if ((ColourNum =
	     FindField(HeIDensity, FieldType, NumberOfBaryonFields)) >= 0) 
	  colnum[++NumberOfColours] = ColourNum;
	//   HeIIDensity
	if ((ColourNum =
	     FindField(HeIIDensity, FieldType, NumberOfBaryonFields)) >= 0) 
	  colnum[++NumberOfColours] = ColourNum;
	//   HeIIIDensity
	if ((ColourNum =
	     FindField(HeIIIDensity, FieldType, NumberOfBaryonFields)) >= 0) 
	  colnum[++NumberOfColours] = ColourNum;
	//   HMDensity
	if ((ColourNum =
	     FindField(HMDensity, FieldType, NumberOfBaryonFields)) >= 0) 
	  colnum[++NumberOfColours] = ColourNum;
	//   H2IDensity
	if ((ColourNum =
	     FindField(H2IDensity, FieldType, NumberOfBaryonFields)) >= 0) 
	  colnum[++NumberOfColours] = ColourNum;
	//   H2IIDensity
	if ((ColourNum =
	     FindField(H2IIDensity, FieldType, NumberOfBaryonFields)) >= 0) 
	  colnum[++NumberOfColours] = ColourNum;
	//   DIDensity
	if ((ColourNum =
	     FindField(DIDensity, FieldType, NumberOfBaryonFields)) >= 0) 
	  colnum[++NumberOfColours] = ColourNum;
	//   DIIDensity
	if ((ColourNum =
	     FindField(DIIDensity, FieldType, NumberOfBaryonFields)) >= 0) 
	  colnum[++NumberOfColours] = ColourNum;
	//   HDIDensity
	if ((ColourNum =
	     FindField(HDIDensity, FieldType, NumberOfBaryonFields)) >= 0) 
	  colnum[++NumberOfColours] = ColourNum;
      }
    }

    /* Add "real" colour fields (metallicity, etc.) as colour variables. */

    int SNColourNum, MetalNum, MBHColourNum, Galaxy1ColourNum, Galaxy2ColourNum,
      MetalIaNum, MetalIINum; 

    if (this->IdentifyColourFields(SNColourNum, MetalNum, MetalIaNum, MetalIINum,
                MBHColourNum, Galaxy1ColourNum, Galaxy2ColourNum) == FAIL)
      ENZO_FAIL("Error in grid->IdentifyColourFields.\n");

    if (MetalNum != -1) {
      colnum[NumberOfColours++] = MetalNum;
      if (MultiMetals || TestProblemData.MultiMetals) {
	colnum[NumberOfColours++] = MetalNum+1; //ExtraType0
	colnum[NumberOfColours++] = MetalNum+2; //ExtraType1
      }
    }

    if (MetalIaNum       != -1) colnum[NumberOfColours++] = MetalIaNum;
    if (MetalIINum       != -1) colnum[NumberOfColours++] = MetalIINum;
    if (SNColourNum      != -1) colnum[NumberOfColours++] = SNColourNum;
    if (MBHColourNum     != -1) colnum[NumberOfColours++] = MBHColourNum;
    if (Galaxy1ColourNum != -1) colnum[NumberOfColours++] = Galaxy1ColourNum;
    if (Galaxy2ColourNum != -1) colnum[NumberOfColours++] = Galaxy2ColourNum;


    /* Add Simon Glover's chemistry species as color fields */

    if(TestProblemData.GloverChemistryModel){

      // Declarations for Simon Glover's cooling.
      int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
	DINum, DIINum, HDINum;

      int CINum,CIINum,OINum,OIINum,SiINum,SiIINum,SiIIINum,CHINum,CH2INum,
	CH3IINum,C2INum,COINum,HCOIINum,OHINum,H2OINum,O2INum;

      int GCM = TestProblemData.GloverChemistryModel;  // purely for convenience

      if (IdentifyGloverSpeciesFields(HIINum,HINum,H2INum,DINum,DIINum,HDINum,
				      HeINum,HeIINum,HeIIINum,CINum,CIINum,OINum,
				      OIINum,SiINum,SiIINum,SiIIINum,CHINum,CH2INum,
				      CH3IINum,C2INum,COINum,HCOIINum,OHINum,H2OINum,
				      O2INum) == FAIL) {
	ENZO_FAIL("Error in IdentifyGloverSpeciesFields.");
      }

      colnum[NumberOfColours++] = HIINum;
      colnum[NumberOfColours++] = HINum;
      colnum[NumberOfColours++] = H2INum;

      if( (GCM==1) || (GCM==2) || (GCM==3) || (GCM==7) ){
	colnum[NumberOfColours++] = DINum;
	colnum[NumberOfColours++] = DIINum;
	colnum[NumberOfColours++] = HDINum;
	colnum[NumberOfColours++] = HeINum;
	colnum[NumberOfColours++] = HeIINum;
	colnum[NumberOfColours++] = HeIIINum;
      }

      if( (GCM==3) || (GCM==5) || (GCM==7) ){
	colnum[NumberOfColours++] = COINum;
      }

      if( (GCM==2) || (GCM==3) || (GCM==7) ){
	colnum[NumberOfColours++] = CINum;
	colnum[NumberOfColours++] = CIINum;
	colnum[NumberOfColours++] = OINum;
	colnum[NumberOfColours++] = OIINum;
      }

      if( (GCM==2) || (GCM==3) ){
	colnum[NumberOfColours++] = SiINum;
	colnum[NumberOfColours++] = SiIINum;
	colnum[NumberOfColours++] = SiIIINum;
      }

      if( (GCM==3) || (GCM==7) ){
	colnum[NumberOfColours++] = CHINum;
	colnum[NumberOfColours++] = CH2INum;
	colnum[NumberOfColours++] = CH3IINum;
	colnum[NumberOfColours++] = C2INum;
	colnum[NumberOfColours++] = HCOIINum;
	colnum[NumberOfColours++] = OHINum;
	colnum[NumberOfColours++] = H2OINum;
	colnum[NumberOfColours++] = O2INum;
      }
      
    } // if(TestProblemData.GloverChemistryModel)


    /* Add shock/cosmic ray variables as a colour variable. */

    if(ShockMethod){
      int MachNum, PSTempNum,PSDenNum;
      
      if (IdentifyShockSpeciesFields(MachNum,PSTempNum,PSDenNum) == FAIL) {
	ENZO_FAIL("Error in IdentifyShockSpeciesFields.")
      }
      
      colnum[NumberOfColours++] = MachNum;
      if(StorePreShockFields){
	colnum[NumberOfColours++] = PSTempNum;
	colnum[NumberOfColours++] = PSDenNum;
      }
    }
    /* Determine if Gamma should be a scalar or a field. */
    
    int UseGammaField = FALSE;
    float *GammaField = NULL;
    if (HydroMethod == Zeus_Hydro && MultiSpecies > 1) {
      UseGammaField = TRUE;
      GammaField = new float[size];
      if (this->ComputeGammaField(GammaField) == FAIL) {
	ENZO_FAIL("Error in grid->ComputeGammaField.");
      }
    } else {
      GammaField = new float[1];
      GammaField[0] = Gamma;

    }
    
    /* Set lowest level flag (used on Zeus hydro). */

    int LowestLevel = (level > MaximumRefinementLevel-1) ? TRUE : FALSE;

    /* Set minimum support (used natively in zeus hydro). */

    float MinimumSupportEnergyCoefficient = 0;
    if (UseMinimumPressureSupport == TRUE && level > MaximumRefinementLevel-1)
      if (this->SetMinimumSupport(MinimumSupportEnergyCoefficient) == FAIL) {
	ENZO_FAIL("Error in grid->SetMinimumSupport,");
      }

    /* allocate space for fluxes */

    /* Set up our restart dump fluxes container */
    this->SubgridFluxStorage = SubgridFluxes;
    this->NumberOfSubgrids = NumberOfSubgrids;

    for (i = 0; i < NumberOfSubgrids; i++) {
      for (dim = 0; dim < GridRank; dim++)  {

	/* compute size (in floats) of flux storage */

        size = 1;
        for (j = 0; j < GridRank; j++)
          size *= SubgridFluxes[i]->LeftFluxEndGlobalIndex[dim][j] -
                  SubgridFluxes[i]->LeftFluxStartGlobalIndex[dim][j] + 1;

	/* set unused dims (for the solver, which is hardwired for 3d). */

        for (j = GridRank; j < 3; j++) {
          SubgridFluxes[i]->LeftFluxStartGlobalIndex[dim][j] = 0;
          SubgridFluxes[i]->LeftFluxEndGlobalIndex[dim][j] = 0;
          SubgridFluxes[i]->RightFluxStartGlobalIndex[dim][j] = 0;
          SubgridFluxes[i]->RightFluxEndGlobalIndex[dim][j] = 0;
        }

	/* Allocate space (if necessary). */

        for (field = 0; field < NumberOfBaryonFields; field++) {
	  //
	  if (SubgridFluxes[i]->LeftFluxes[field][dim] == NULL)
	    SubgridFluxes[i]->LeftFluxes[field][dim]  = new float[size];
	  if (SubgridFluxes[i]->RightFluxes[field][dim] == NULL)
	    SubgridFluxes[i]->RightFluxes[field][dim] = new float[size];
	  for (n = 0; n < size; n++) {
	    SubgridFluxes[i]->LeftFluxes[field][dim][n] = 0;
	    SubgridFluxes[i]->RightFluxes[field][dim][n] = 0;
	  }
        }

	for (field = NumberOfBaryonFields; field < MAX_NUMBER_OF_BARYON_FIELDS;
	     field++) {
          SubgridFluxes[i]->LeftFluxes[field][dim] = NULL;
          SubgridFluxes[i]->RightFluxes[field][dim] = NULL;
	}

      }  // next dimension

      /* make things pretty */

      for (dim = GridRank; dim < 3; dim++)
        for (field = 0; field < MAX_NUMBER_OF_BARYON_FIELDS; field++) {
          SubgridFluxes[i]->LeftFluxes[field][dim] = NULL;
          SubgridFluxes[i]->RightFluxes[field][dim] = NULL;
	}

    } // end of loop over subgrids

    /* compute global start index for left edge of entire grid 
       (including boundary zones) */

    for (dim = 0; dim < GridRank; dim++)
      GridGlobalStart[dim] = nlongint((GridLeftEdge[dim]-DomainLeftEdge[dim])/(*(CellWidth[dim]))) -
	GridStartIndex[dim];

    /* fix grid quantities so they are defined to at least 3 dims */

    for (i = GridRank; i < 3; i++) {
      GridDimension[i]   = 1;
      GridStartIndex[i]  = 0;
      GridEndIndex[i]    = 0;
      GridVelocity[i]    = 0.0;
      GridGlobalStart[i] = 0;
    }

    /* If using comoving coordinates, multiply dx by a(n+1/2).
       In one fell swoop, this recasts the equations solved by solver
       in comoving form (except for the expansion terms which are taken
       care of elsewhere). */

    if (ComovingCoordinates)
      if (CosmologyComputeExpansionFactor(Time+0.5*dtFixed, &a, &dadt) 
	  == FAIL) {
	ENZO_FAIL("Error in CosmologyComputeExpansionFactors.");
      }

    /* Create a cell width array to pass (and convert to absolute coords). */

    float *CellWidthTemp[MAX_DIMENSION];
    for (dim = 0; dim < MAX_DIMENSION; dim++) {
      CellWidthTemp[dim] = new float[GridDimension[dim]];
      for (i = 0; i < GridDimension[dim]; i++)
	if (dim < GridRank)
	  CellWidthTemp[dim][i] = float(a*CellWidth[dim][i]);
	else
	  CellWidthTemp[dim][i] = 1.0;
    }

    /* Prepare Gravity. */

    int GravityOn = 0, FloatSize = sizeof(float);
    if (SelfGravity || UniformGravity || PointSourceGravity || ExternalGravity)
      GravityOn = 1;
#ifdef TRANSFER
    if (RadiationPressure)
      GravityOn = 1;
#endif    

    //Some setup for MHDCT

    float *MagneticFlux[3][2];
    if ( UseMHDCT ) {
        MHDCT_ConvertEnergyToConservedS();  //Energy toggle.  Probably will be removed soon.
        for(field=0;field<3;field++){
            if(ElectricField[field] == NULL ) 
                ElectricField[field] = new float[ElectricSize[field]];
            for(i=0;i<ElectricSize[field]; i++) ElectricField[field][i] = 0.0;
        }
        for(field=0;field<3;field++){
          MagneticFlux[field][0] = new float[2*MagneticSize[field]];
          MagneticFlux[field][1] =  MagneticFlux[field][0] +MagneticSize[field];
          for (i=0; i< 2*MagneticSize[field]; i++) MagneticFlux[field][0][i] = 0.0;
        }
        CenterMagneticField();
#ifdef BIERMANN
        
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
          // Biermann battery constants in cgs units, convert to enzo units later
          float speedoflight = clight/VelocityUnits;
          float hydrogenmass = mh/DensityUnits*POW(LengthUnits,3);
          float electroncharge = 4.803e-10*
                TimeUnits*BFieldUnits/(speedoflight*DensityUnits*POW(LengthUnits,3));
          float chi=1.0;
        }
#endif //BIERMANN
    }//UseMHDCT



    float* Fluxes[3] = {MagneticFlux[0][0],MagneticFlux[1][0],MagneticFlux[2][0]};
    int CurlStart[3] = {0,0,0}, 
    CurlEnd[3] = {GridDimension[0]-1,GridDimension[1]-1,GridDimension[2]-1};
    if ( UseMHDCT ){
        if (HydroMethod == MHD_Li){
          this->SolveMHD_Li(CycleNumber, NumberOfSubgrids, SubgridFluxes, 
                CellWidthTemp, GridGlobalStart, GravityOn, 
                NumberOfColours, colnum, Fluxes);
        }

        if( HydroMethod != NoHydro )
            switch( MHD_CT_Method ){
                case CT_BalsaraSpicer: //1
                case CT_Athena_LF:     //2
                case CT_Athena_Switch: //3
                    ComputeElectricField(dtFixed, Fluxes);
                    break;
                case CT_Biermann:      //4
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
                    break;
#endif //BIERMANN
                case CT_None:
                default:
                    if(MyProcessorNumber == ROOT_PROCESSOR )
                        fprintf(stderr, "Warning: No CT method used with MHD_Li.\n");
                break;
            }
        MHD_Curl( CurlStart,CurlEnd, 1);
        CenterMagneticField();

        MHDCT_ConvertEnergyToSpecificS();
        for(field=0;field<3;field++){
          delete [] MagneticFlux[field][0];
        }
    }
    /* Call Solver on this grid.
       Notice that it is hard-wired for three dimensions, but it does
       the right thing for < 3 dimensions. */
    /* note: Start/EndIndex are zero based */
        
    if (HydroMethod == PPM_DirectEuler)
      this->SolvePPM_DE(CycleNumber, NumberOfSubgrids, SubgridFluxes,
                        CellWidthTemp, GridGlobalStart, GravityOn,
                        NumberOfColours, colnum,
                        MinimumSupportEnergyCoefficient);

    /* PPM LR has been withdrawn. */

    if (HydroMethod == PPM_LagrangeRemap) {
#ifdef PPM_LR
      FORTRAN_NAME(ppm_lr)(
			density, totalenergy, velocity1, velocity2, velocity3,
                          gasenergy,
			&GravityOn, AccelerationField[0],
                           AccelerationField[1],
                           AccelerationField[2],
			&Gamma, &dtFixed, &CycleNumber,
                          CellWidthTemp[0], CellWidthTemp[1], CellWidthTemp[2],
			&GridRank, &GridDimension[0], &GridDimension[1],
                           &GridDimension[2], GridStartIndex, GridEndIndex,
			GridVelocity, &PPMFlatteningParameter,
                           &PressureFree,
			&PPMDiffusionParameter, &PPMSteepeningParameter,
                           &DualEnergyFormalism, &DualEnergyFormalismEta1,
                           &DualEnergyFormalismEta2,
			&NumberOfSubgrids, leftface, rightface,
			istart, iend, jstart, jend,
			standard, dindex, Eindex, uindex, vindex, windex,
			  geindex, temp,
                        &NumberOfColours, colourpt, coloff, colindex);
#else /* PPM_LR */
      ENZO_FAIL("PPM LR is not supported.");
#endif /* PPM_LR */
    }

    if (HydroMethod == Zeus_Hydro)
      if (this->ZeusSolver(GammaField, UseGammaField, CycleNumber, 
               CellWidthTemp[0], CellWidthTemp[1], CellWidthTemp[2],
               GravityOn, NumberOfSubgrids, GridGlobalStart,
               SubgridFluxes,
               NumberOfColours, colnum, LowestLevel,
               MinimumSupportEnergyCoefficient) == FAIL)
	ENZO_FAIL("ZeusSolver() failed!\n");
	


    /* Clean up allocated fields. */

    delete [] GammaField;   

    for (dim = 0; dim < MAX_DIMENSION; dim++)
      delete [] CellWidthTemp[dim];

  /* If we're supposed to be outputting on Density, we need to update
     the current maximum value of that Density. */

    if(OutputOnDensity == 1){
      int DensNum = FindField(Density, FieldType, NumberOfBaryonFields);
      for(i = 0; i < size; i++)
        CurrentMaximumDensity =
            max(BaryonField[DensNum][i], CurrentMaximumDensity);
    }

  }  // end: if (NumberOfBaryonFields > 0)


  this->DebugCheck("SolveHydroEquations (after)");

  TIMER_STOP("SolveHydroEquations");
  LCAPERF_STOP("grid_SolveHydroEquations");
  return SUCCESS;

}
