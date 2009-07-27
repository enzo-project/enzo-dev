/***********************************************************************
/
/  GRID CLASS (INITIALIZE THE GRID FOR SEDOV BLAST WAVE TEST)
/
/  written by: Brian O'Shea
/  date:       December 2007
/  modified1:  
/
/  PURPOSE: Sets the energy in the initial explosion region, as well as
/           setting color fields, species fields, and maybe even kinetic
/           energy fields.
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/
 
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

int FindField(int field, int farray[], int numfields);

void mt_init(unsigned_int seed);

unsigned_long_int mt_random();

float cell_fraction(FLOAT cellx, FLOAT celly, FLOAT cellz, FLOAT shockx, FLOAT shocky, FLOAT shockz, FLOAT dx, FLOAT radius, int GridRank);

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, double *MassUnits, FLOAT Time);

void set_analytic_sedov(int Nbins, double *radius, double *density,
			double *pressure, double *velocity,
			double timesolve, double gamma, double beta,
			double explosion_energy, double p_ambient, 
			double rho_ambient);

double compute_sedov_v(double xi, double gamma, double alpha1, double alpha2);

double sedov_vfunc(double V, double gamma, double alpha1, double alpha2);

int grid::RadiatingShockInitializeGrid(FLOAT dr,
				       float RadiatingShockInnerDensity,
				       float RadiatingShockInnerTotalEnergy,
				       int RadiatingShockUseDensityFluctuations,
				       int RadiatingShockRandomSeed,
				       float RadiatingShockDensityFluctuationLevel,
				       int RadiatingShockInitializeWithKE,
				       int RadiatingShockUseSedovProfile,
				       FLOAT RadiatingShockSedovBlastRadius,
				       float RadiatingShockEnergy,
				       float RadiatingShockPressure,
				       float RadiatingShockKineticEnergyFraction,
				       float RadiatingShockRhoZero,
				       float RadiatingShockVelocityZero,
				       int RadiatingShockRandomSeedInitialize,
				       FLOAT RadiatingShockCenterPosition[MAX_DIMENSION])
{
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  /* declarations */
 
  int size = 1, dim;
  for (dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];

  FLOAT r,x,y,z;

  float cell_HydrogenFractionByMass, cell_DeuteriumToHydrogenRatio;

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
    DINum, DIINum, HDINum, MetalNum;
  int CINum, CIINum, OINum, OIINum, SiINum, SiIINum, SiIIINum, CHINum, CH2INum, 
    CH3IINum, C2INum, COINum, HCOIINum, OHINum, H2OINum, O2INum;

  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
        ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }

  if (MultiSpecies)
    if (IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
                      HMNum, H2INum, H2IINum, DINum, DIINum, HDINum) == FAIL) {
            ENZO_FAIL("Error in grid->IdentifySpeciesFields.");
    }

  // If using Simon Glover's non-equilibrium cooling, set up fields.
  if (GloverChemistryModel) 
    if (IdentifyGloverSpeciesFields(HIINum,HINum,H2INum,DINum,DIINum,HDINum,
				    HeINum,HeIINum,HeIIINum,CINum,CIINum,OINum,
				    OIINum,SiINum,SiIINum,SiIIINum,CHINum,CH2INum,
				    CH3IINum,C2INum,COINum,HCOIINum,OHINum,H2OINum,
				    O2INum) == FAIL) {
      fprintf(stderr,"Error in IdentifyGloverSpeciesFields.\n");
      return FAIL;
    }

  int MetallicityField = FALSE;
  if ((MetalNum = FindField(Metallicity, FieldType, NumberOfBaryonFields))
      != -1)
    MetallicityField = TRUE;
  else
    MetalNum = 0;


  FILE *fptr;
  double *sedovradius=NULL,*sedovdensity=NULL,*sedovpressure=NULL,*sedovvelocity=NULL;
  int numbins=1000;


  if(RadiatingShockUseSedovProfile){

    if(debug)
      printf("input pressure is %e\n",RadiatingShockPressure);

    float DensityUnits=1.0, LengthUnits=1.0, TemperatureUnits=1.0, TimeUnits=1.0,
      VelocityUnits=1.0;
    double MassUnits=1.0;

    if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
		 &TimeUnits, &VelocityUnits, &MassUnits, 0.0) == FAIL) {
            ENZO_FAIL("Error in GetUnits.");
    }
 
    sedovradius = new double[numbins];
    sedovdensity = new double[numbins];
    sedovpressure = new double[numbins];
    sedovvelocity = new double[numbins];

    double logminrad, logmaxrad, logdeltar;

    logmaxrad = log10( double(RadiatingShockSedovBlastRadius) * LengthUnits * 1.1);
    logminrad = log10( double(RadiatingShockSedovBlastRadius) * LengthUnits * 1.1 * 0.001);

    logdeltar = (logmaxrad - logminrad) / double(numbins);

    double rho_ambient, p_ambient, explosion_energy, shocktime;

    explosion_energy = double(RadiatingShockEnergy)*1.0e+51;
    rho_ambient = double(BaryonField[DensNum][0])*double(DensityUnits);
    p_ambient =  double(RadiatingShockPressure) * double(DensityUnits) * double(VelocityUnits) * double(VelocityUnits);

    printf("explosion_energy is %e\n",explosion_energy);
    printf("rho_ambient is %e\n",rho_ambient);
    printf("p_ambient is %e\n",p_ambient);

    shocktime = POW( explosion_energy / rho_ambient, -0.5)
      * POW( double(RadiatingShockSedovBlastRadius) * LengthUnits, 2.5) / POW( 1.15, 2.5);

    if(debug)
      printf("shock time is %e for radius %e\n",
	     shocktime, double(RadiatingShockSedovBlastRadius) * LengthUnits);
    
    for(int i=0; i<numbins; i++){
      sedovradius[i] = POW(10.0, logminrad + (double(i)+1.0)*logdeltar);
    }

    // here's where all of the magic happens
    set_analytic_sedov(numbins,sedovradius,sedovdensity,sedovpressure,sedovvelocity,
		       shocktime,double(Gamma),1.15,explosion_energy,
		       p_ambient,rho_ambient);

    // convert sedov values (CGS) to Enzo internal energy units!
    for(int i=0; i<numbins; i++){
      if(debug)
	printf("Pre-conversion (CGS) values: %e %e %e %e  (%d)\n",
	       sedovradius[i],sedovdensity[i],sedovpressure[i],sedovvelocity[i], i);

      sedovradius[i] /= double(LengthUnits);
      sedovvelocity[i] /= double(VelocityUnits);

      // convert sedovpressure to internal energy (ergs/g)
      if(sedovdensity[i] > 1.0e-30){
	sedovpressure[i] /= ( (Gamma-1.0)*sedovdensity[i]);
      } else {
	sedovpressure[i] = 1.0e-20;
      }

      sedovdensity[i] /= double(DensityUnits);
      sedovpressure[i] /= (double(VelocityUnits)*double(VelocityUnits));

      if(sedovdensity[i] < 1.0e-10) sedovdensity[i] = 1.0e-10;
      if(sedovpressure[i] < 1.0e-10) sedovpressure[i] = 1.0e-10;

      if(debug){
	printf("Post-conversion (enzo units):  %e %e %e %e (%d)\n",
	       sedovradius[i],sedovdensity[i],sedovpressure[i],sedovvelocity[i], i);
	fflush(stdout);
      }

    } // for(int i=0; i<numbins; i++)

  } // if(RadiatingShockUseSedovProfile)


  /* set fields in the initial explosion region: x^2+y^2+z^2 < dr^2. */
 
  int i, j, k;

  float therandomfraction;

  if( RadiatingShockRandomSeedInitialize == 0)
    mt_init(((unsigned_int) RadiatingShockRandomSeed));

  unsigned_long_int therandominteger;

  float outside_rho=0.0, outside_TE=0.0, outside_GE=0.0, cellfraction;

  int cellindex;

  outside_rho =  BaryonField[DensNum][0];

  if(HydroMethod==2){  // ZEUS

    outside_TE = BaryonField[TENum][0];

  } else { // PPM

    outside_TE = BaryonField[TENum][0];
    
    if(DualEnergyFormalism){
      outside_GE = BaryonField[GENum][0];
    }

  }  // if(HydroMethod==2)

  for (k = 0; k < GridDimension[2]; k++)
    for (j = 0; j < GridDimension[1]; j++)
      for (i = 0; i < GridDimension[0]; i++) {
 
	/* Compute position */
	x=y=z=0.0;

	cellindex = i + j*GridDimension[0];
	if(GridRank > 2)
	  cellindex += k*GridDimension[0]*GridDimension[1];

	if ((i+j+k) == 0) {
	  if (MultiSpecies) {
	    fprintf(stderr,"External medium: Density: %.2"ESYM, BaryonField[DensNum][cellindex]);
	    fprintf(stderr,", HI: %.2"ESYM", HII: %.2"ESYM", HeI: %.2"ESYM,BaryonField[HINum][cellindex],
		    BaryonField[HIINum][cellindex],
		    BaryonField[HeINum][cellindex]);
	    fprintf(stderr,", HeII: %.2"ESYM", HeIII: %.2"ESYM", De: %.2"ESYM,BaryonField[HeIINum][cellindex],
		    BaryonField[HeIIINum][cellindex],BaryonField[DeNum][cellindex]);
	    if (MultiSpecies > 1) {
	      fprintf(stderr,", H2I: %.2"ESYM", H2II: %.2"ESYM", HM: %.2"ESYM,BaryonField[H2INum][cellindex],
		      BaryonField[H2IINum][cellindex],BaryonField[HMNum][cellindex]);
	    }
	    if (MultiSpecies > 2) {
	      fprintf(stderr,", DI: %.2"ESYM", DII: %.2"ESYM", HDI: %.2"ESYM,BaryonField[DINum][cellindex],
		      BaryonField[DIINum][cellindex],BaryonField[HDINum][cellindex]);
	    }
	    if (TestProblemData.UseMetallicityField) {
	      fprintf(stderr,", Metal: %.2"ESYM,BaryonField[MetalNum][cellindex]);
	    }
	    fprintf(stderr,"\n");
	  }
	  else if (GloverChemistryModel) {
	    int GCM = TestProblemData.GloverChemistryModel;  // purely for convenience
	    fprintf(stderr,"      External medium: Density: %.2"ESYM, BaryonField[DensNum][cellindex]);
	    fprintf(stderr,", HI: %.2"ESYM", HII: %.2"ESYM", H2I: %.2"ESYM,BaryonField[HINum][cellindex],
		    BaryonField[HIINum][cellindex],
		    BaryonField[H2INum][cellindex]);
	    if( (GCM==1) || (GCM==2) || (GCM==3) || (GCM==7) ){
	      fprintf(stderr,", HeI: %.2"ESYM", HeII: %.2"ESYM", HeIII: %.2"ESYM,BaryonField[HeINum][cellindex],
		      BaryonField[HeIINum][cellindex],BaryonField[HeIIINum][cellindex]);
	      fprintf(stderr,", DI: %.2"ESYM", DII: %.2"ESYM", HDI: %.2"ESYM,BaryonField[DINum][cellindex],
		      BaryonField[DIINum][cellindex],BaryonField[HDINum][cellindex]);
	    }
	    fprintf(stderr,"\n");
	  }
	}
 
	x = CellLeftEdge[0][i] + 0.5*CellWidth[0][i];
	if (GridRank > 1)
	  y = CellLeftEdge[1][j] + 0.5*CellWidth[1][j];
	if (GridRank > 2)
	  z = CellLeftEdge[2][k] + 0.5*CellWidth[2][k];
	

	/* Find distance from center. */
 
	// it's REALLY r^2 right now
	r = POW(x-RadiatingShockCenterPosition[0], 2.0) +
	  POW(y-RadiatingShockCenterPosition[1], 2.0);
      
	if(GridRank > 2)
	  r+= POW(z-RadiatingShockCenterPosition[2], 2.0 );
  
	r = sqrt(r);  // ok, now it's just r
	
	//r = max(r, 0.1*CellWidth[0][0]);
	

	if (r <= dr + 1.5*CellWidth[0][i]) {

	  cellfraction = cell_fraction(x,y,z, RadiatingShockCenterPosition[0], 
				       RadiatingShockCenterPosition[1],
				       RadiatingShockCenterPosition[2], 
				       CellWidth[0][i], dr, GridRank);

	  // calculate H mass fraction and D/H ratio for this cell.
	  // It is a combination of the inner and outer values.
	  cell_HydrogenFractionByMass = ((cellfraction * TestProblemData.InnerHydrogenFractionByMass 
					  * RadiatingShockInnerDensity) +
					 ((1-cellfraction) * TestProblemData.HydrogenFractionByMass * outside_rho)) /
	    ((cellfraction * RadiatingShockInnerDensity) + ((1-cellfraction) * outside_rho));

	  cell_DeuteriumToHydrogenRatio = ((cellfraction * TestProblemData.InnerDeuteriumToHydrogenRatio * 
					    TestProblemData.InnerHydrogenFractionByMass * RadiatingShockInnerDensity) +
					   ((1-cellfraction) * TestProblemData.DeuteriumToHydrogenRatio *
					    TestProblemData.HydrogenFractionByMass * outside_rho)) /
	    ((cellfraction * TestProblemData.InnerHydrogenFractionByMass * RadiatingShockInnerDensity) +
	     ((1-cellfraction) * TestProblemData.HydrogenFractionByMass * outside_rho));
	  
	  if(RadiatingShockInitializeWithKE){

	    if(RadiatingShockUseSedovProfile){

	      int sedovindex=-1;

	      for(int nnn=0; nnn < numbins-1; nnn++)
		if(sedovradius[nnn] < r && r <= sedovradius[nnn+1]) sedovindex=nnn;

	      // if we're outside the radius of the shock, set to the outside value (last value in array)
	      // this is set above to be consistent with the desired outer density value
	      if(r >= sedovradius[numbins-1]) sedovindex = numbins-1;

	      if(sedovindex < 0 || sedovindex >= numbins){
		fprintf(stderr,"Grid:RadiatingShockInitializeGrid: Argh!  %d  %e  %e  %e\n",
			sedovindex, r, sedovradius[0], sedovradius[numbins-1]);
		return FAIL;
	      }

	      // density
	      BaryonField[DensNum][cellindex] = sedovdensity[sedovindex];

	      // x,y, and maybe z velocity.  Note that we're ADDING the supernova
	      // velocities because Grid_InitializeUniformGrid may have set the
	      // gas to have some ambient velocity 
	      BaryonField[Vel1Num][cellindex] += sedovvelocity[sedovindex]
		* (x-RadiatingShockCenterPosition[0]) / dr;

	      BaryonField[Vel2Num][cellindex] += sedovvelocity[sedovindex]
		* (y-RadiatingShockCenterPosition[1]) / dr;
	      
	      if(GridRank > 2)
		BaryonField[Vel3Num][cellindex] += sedovvelocity[sedovindex]
		  * (z-RadiatingShockCenterPosition[2]) / dr;
	      
	      if(HydroMethod == 2){
		
		// ZEUS
		BaryonField[TENum][cellindex] = sedovpressure[sedovindex];  // sedovpressure is really internal energy
		
	      } else {
		
		// PPM
		BaryonField[TENum][cellindex] = sedovpressure[sedovindex]  // sedovpressure is really internal energy  
		  + 0.5 * BaryonField[Vel1Num][cellindex] * BaryonField[Vel1Num][cellindex]
		  + 0.5 * BaryonField[Vel2Num][cellindex] * BaryonField[Vel2Num][cellindex];
		
		if(GridRank > 2)
		  BaryonField[TENum][cellindex]+= 
		    + 0.5 * BaryonField[Vel3Num][cellindex] * BaryonField[Vel3Num][cellindex];
		
	      } // if(HydroMethod == 2)
	      
	      // gas energy (PPM dual energy formalims)
	      if(DualEnergyFormalism)
		BaryonField[GENum][cellindex] = sedovpressure[sedovindex];  // sedovpressure is really internal energy

	    } else {  // if(RadiatingShockUseSedovProfile)

	      // density
	      BaryonField[DensNum][cellindex] = RadiatingShockRhoZero * r / dr;

	      // x,y, and maybe z velocity.  Note that we're ADDING the supernova
	      // velocities because Grid_InitializeUniformGrid may have set the
	      // gas to have some ambient density 
	      BaryonField[Vel1Num][cellindex] += RadiatingShockVelocityZero
		* (x-RadiatingShockCenterPosition[0]) / dr;

	      BaryonField[Vel2Num][cellindex] += RadiatingShockVelocityZero
		* (y-RadiatingShockCenterPosition[1]) / dr;
	      
	      if(GridRank > 2)
		BaryonField[Vel3Num][cellindex] += RadiatingShockVelocityZero
		  * (z-RadiatingShockCenterPosition[2]) / dr;
	      
	      if(HydroMethod == 2){
		
		// ZEUS
		BaryonField[TENum][cellindex] = RadiatingShockInnerTotalEnergy*(1.0-RadiatingShockKineticEnergyFraction);
		
	      } else {
		
		// PPM
		BaryonField[TENum][cellindex] = RadiatingShockInnerTotalEnergy*(1.0-RadiatingShockKineticEnergyFraction)
		  + 0.5 * BaryonField[Vel1Num][cellindex] * BaryonField[Vel1Num][cellindex]
		  + 0.5 * BaryonField[Vel2Num][cellindex] * BaryonField[Vel2Num][cellindex];
		
		if(GridRank > 2)
		  BaryonField[TENum][cellindex]+= 
		    + 0.5 * BaryonField[Vel3Num][cellindex] * BaryonField[Vel3Num][cellindex];
		
	      } // if(HydroMethod == 2)
	      
	      // gas energy (PPM dual energy formalims)
	      if(DualEnergyFormalism)
		BaryonField[GENum][cellindex] = RadiatingShockInnerTotalEnergy*(1.0-RadiatingShockKineticEnergyFraction);
	      
	    }  // if(RadiatingShockUseSedovProfile)

	  } // if(RadiatingShockInitializeWithKE)
	  else {  // if not, initialize just thermal energy (no kinetic energy)

	    BaryonField[DensNum][cellindex] = cellfraction*RadiatingShockInnerDensity + (1.0-cellfraction)*outside_rho;

	    if(HydroMethod==2){

	      BaryonField[TENum][cellindex] = cellfraction*RadiatingShockInnerTotalEnergy + (1.0-cellfraction)*outside_TE;

	    }
	    else{

	      // PPM - remember, we still have to take care of the kinetic energy in gas with
	      // some bulk flow
	      BaryonField[TENum][cellindex] = cellfraction*RadiatingShockInnerTotalEnergy + (1.0-cellfraction)*outside_TE
		+ 0.5 * BaryonField[Vel1Num][cellindex] * BaryonField[Vel1Num][cellindex]
		+ 0.5 * BaryonField[Vel2Num][cellindex] * BaryonField[Vel2Num][cellindex];

	      if(GridRank > 2)
		BaryonField[TENum][cellindex]+= 
		  + 0.5 * BaryonField[Vel3Num][cellindex] * BaryonField[Vel3Num][cellindex];	
	    }

	    if(DualEnergyFormalism)
	      BaryonField[GENum][cellindex] = cellfraction*RadiatingShockInnerTotalEnergy + (1.0-cellfraction)*outside_GE;

	  } // if(RadiatingShockInitializeWithKE)

	  /* If we have metals or anything turned on, set those to the user-specified fractional values */

	  // Set multispecies fields!
	  // this attempts to set them such that species conservation is maintained,
	  // using the method in CosmologySimulationInitializeGrid.C
	  if(TestProblemData.MultiSpecies) {

	    BaryonField[HIINum][cellindex] = TestProblemData.HII_Fraction_Inner * 
	      cell_HydrogenFractionByMass * BaryonField[DensNum][cellindex];
	      
	    BaryonField[HeIINum][cellindex] = TestProblemData.HeII_Fraction_Inner *
	      BaryonField[DensNum][cellindex] * (1.0-cell_HydrogenFractionByMass);
	      
	    BaryonField[HeIIINum][cellindex] = TestProblemData.HeIII_Fraction_Inner *
	      BaryonField[DensNum][cellindex] * (1.0-cell_HydrogenFractionByMass);

	    BaryonField[HeINum][cellindex] = 
	      (1.0 - cell_HydrogenFractionByMass)*BaryonField[DensNum][cellindex] -
	      BaryonField[HeIINum][cellindex] - BaryonField[HeIIINum][cellindex];
	      
	    if(TestProblemData.MultiSpecies > 1){
	      BaryonField[HMNum][cellindex] = TestProblemData.HM_Fraction_Inner *
		BaryonField[HIINum][cellindex];
		
	      BaryonField[H2INum][cellindex] = TestProblemData.H2I_Fraction_Inner *
		BaryonField[0][cellindex] * cell_HydrogenFractionByMass;
		
	      BaryonField[H2IINum][cellindex] = TestProblemData.H2II_Fraction_Inner * 2.0 *
		BaryonField[HIINum][cellindex];
	    }

	    // HI density is calculated by subtracting off the various ionized fractions
	    // from the total
	    BaryonField[HINum][cellindex] = cell_HydrogenFractionByMass*BaryonField[0][cellindex]
	      - BaryonField[HIINum][cellindex];
	    if (MultiSpecies > 1)
	      BaryonField[HINum][cellindex] -= (BaryonField[HMNum][cellindex] + BaryonField[H2IINum][cellindex]
						+ BaryonField[H2INum][cellindex]);

	    // Electron "density" (remember, this is a factor of m_p/m_e scaled from the 'normal'
	    // density for convenience) is calculated by summing up all of the ionized species.
	    // The factors of 0.25 and 0.5 in front of HeII and HeIII are to fix the fact that we're
	    // calculating mass density, not number density (because the BaryonField values are 4x as
	    // heavy for helium for a single electron)
	    BaryonField[DeNum][cellindex] = BaryonField[HIINum][cellindex] +
	      0.25*BaryonField[HeIINum][cellindex] + 0.5*BaryonField[HeIIINum][cellindex];
	    if (MultiSpecies > 1)
	      BaryonField[DeNum][cellindex] += 0.5*BaryonField[H2IINum][cellindex] -
		BaryonField[HMNum][cellindex];
	      
	    // Set deuterium species (assumed to be a negligible fraction of the total, so not
	    // counted in the conservation)
	    if(TestProblemData.MultiSpecies > 2){
	      BaryonField[DINum ][cellindex]  = cell_DeuteriumToHydrogenRatio * BaryonField[HINum][cellindex];
	      BaryonField[DIINum][cellindex] = cell_DeuteriumToHydrogenRatio * BaryonField[HIINum][cellindex];
	      BaryonField[HDINum][cellindex] = 0.75 * cell_DeuteriumToHydrogenRatio * BaryonField[H2INum][cellindex];
	    }

	  } // if(TestProblemData.MultiSpecies)

	  if(TestProblemData.UseMetallicityField){
	    BaryonField[MetalNum][cellindex] = (cellfraction * TestProblemData.MetallicityField_Fraction +
						(1-cellfraction)*tiny_number)
	      * BaryonField[DensNum][cellindex];

	    if(TestProblemData.MultiMetals){
	      BaryonField[MetalNum+1][cellindex] = (cellfraction * TestProblemData.MultiMetalsField1_Fraction +
						    (1-cellfraction)*tiny_number)
		* BaryonField[DensNum][cellindex];
	      BaryonField[MetalNum+2][cellindex] = (cellfraction * TestProblemData.MultiMetalsField2_Fraction +
						    (1-cellfraction)*tiny_number)
		* BaryonField[DensNum][cellindex];
	    }
	  } // if(TestProblemData.UseMetallicityField)

	  if(TestProblemData.GloverChemistryModel){
	    float tempHM,tempH2II;
	    int GCM = TestProblemData.GloverChemistryModel;  // purely for convenience

	    BaryonField[HIINum][cellindex] = TestProblemData.HII_Fraction_Inner * 
	      cell_HydrogenFractionByMass * BaryonField[DensNum][cellindex];
	    tempHM = TestProblemData.HM_Fraction_Inner * BaryonField[HIINum][cellindex];
	    BaryonField[H2INum][cellindex] = TestProblemData.H2I_Fraction_Inner *
		BaryonField[0][cellindex] * cell_HydrogenFractionByMass;
	    tempH2II = TestProblemData.H2II_Fraction_Inner * 2.0 * BaryonField[HIINum][cellindex];

	    // HI density is calculated by subtracting off the various ionized fractions
	    // from the total
	    BaryonField[HINum][cellindex] = cell_HydrogenFractionByMass*BaryonField[0][cellindex]
	      - BaryonField[HIINum][cellindex];
	    BaryonField[HINum][cellindex] -= (tempHM + tempH2II + BaryonField[H2INum][cellindex]);

	    if( (GCM==1) || (GCM==2) || (GCM==3) || (GCM==7) ){
	      BaryonField[DINum   ][cellindex] = TestProblemData.DI_Fraction_Inner * BaryonField[DensNum][cellindex];
	      BaryonField[DIINum  ][cellindex] = TestProblemData.DII_Fraction_Inner * BaryonField[DensNum][cellindex];
	      BaryonField[HDINum  ][cellindex] = TestProblemData.HDI_Fraction_Inner * BaryonField[DensNum][cellindex];
	      BaryonField[HeIINum][cellindex] = TestProblemData.HeII_Fraction_Inner *
		BaryonField[DensNum][cellindex] * (1.0-cell_HydrogenFractionByMass);
	      BaryonField[HeIIINum][cellindex] = TestProblemData.HeIII_Fraction_Inner *
		BaryonField[DensNum][cellindex] * (1.0-cell_HydrogenFractionByMass);
	      BaryonField[HeINum][cellindex] = 
		(1.0 - cell_HydrogenFractionByMass)*BaryonField[DensNum][cellindex] -
		BaryonField[HeIINum][cellindex] - BaryonField[HeIIINum][cellindex];
	    }
	      
	    if( (GCM==3) || (GCM==5) || (GCM==7) )
	      BaryonField[COINum  ][cellindex] = TestProblemData.COI_Fraction_Inner * BaryonField[DensNum][cellindex];
	      
	    if( (GCM==2) || (GCM==3) || (GCM==7) ){
	      BaryonField[CINum ][cellindex] = TestProblemData.CI_Fraction_Inner * BaryonField[DensNum][cellindex];
	      BaryonField[CIINum][cellindex] = TestProblemData.CII_Fraction_Inner * BaryonField[DensNum][cellindex];
	      BaryonField[OINum ][cellindex] = TestProblemData.OI_Fraction_Inner * BaryonField[DensNum][cellindex];
	      BaryonField[OIINum][cellindex] = TestProblemData.OII_Fraction_Inner * BaryonField[DensNum][cellindex];
	    }

	    if( (GCM==2) || (GCM==3) ){
	      BaryonField[SiINum  ][cellindex] = TestProblemData.SiI_Fraction_Inner * BaryonField[DensNum][cellindex];
	      BaryonField[SiIINum ][cellindex] = TestProblemData.SiII_Fraction_Inner * BaryonField[DensNum][cellindex];
	      BaryonField[SiIIINum][cellindex] = TestProblemData.SiIII_Fraction_Inner * BaryonField[DensNum][cellindex];
	    }

	    if( (GCM==3) || (GCM==7) ){
	      BaryonField[CHINum  ][cellindex] = TestProblemData.CHI_Fraction_Inner * BaryonField[DensNum][cellindex];
	      BaryonField[CH2INum ][cellindex] = TestProblemData.CH2I_Fraction_Inner * BaryonField[DensNum][cellindex];
	      BaryonField[CH3IINum][cellindex] = TestProblemData.CH3II_Fraction_Inner * BaryonField[DensNum][cellindex];
	      BaryonField[C2INum  ][cellindex] = TestProblemData.C2I_Fraction_Inner * BaryonField[DensNum][cellindex];
	      BaryonField[HCOIINum][cellindex] = TestProblemData.HCOII_Fraction_Inner * BaryonField[DensNum][cellindex];
	      BaryonField[OHINum  ][cellindex] = TestProblemData.OHI_Fraction_Inner * BaryonField[DensNum][cellindex];
	      BaryonField[H2OINum ][cellindex] = TestProblemData.H2OI_Fraction_Inner * BaryonField[DensNum][cellindex];
	      BaryonField[O2INum  ][cellindex] = TestProblemData.O2I_Fraction_Inner * BaryonField[DensNum][cellindex];
	    }

	} // if(TestProblemData.GloverChemistryModel)

	  if (MultiSpecies) {
	    fprintf(stderr,"Cell fraction: %.4"FSYM", Density: %.2"ESYM, cellfraction,BaryonField[DensNum][cellindex]);
	    fprintf(stderr,", HI: %.2"ESYM", HII: %.2"ESYM", HeI: %.2"ESYM,BaryonField[HINum][cellindex],BaryonField[HIINum][cellindex],
		    BaryonField[HeINum][cellindex]);
	    fprintf(stderr,", HeII: %.2"ESYM", HeIII: %.2"ESYM", De: %.2"ESYM,BaryonField[HeIINum][cellindex],
		    BaryonField[HeIIINum][cellindex],BaryonField[DeNum][cellindex]);
	    if (MultiSpecies > 1) {
	      fprintf(stderr,", H2I: %.2"ESYM", H2II: %.2"ESYM", HM: %.2"ESYM,BaryonField[H2INum][cellindex],
		      BaryonField[H2IINum][cellindex],BaryonField[HMNum][cellindex]);
	    }
	    if (MultiSpecies > 2) {
	      fprintf(stderr,", DI: %.2"ESYM", DII: %.2"ESYM", HDI: %.2"ESYM,BaryonField[DINum][cellindex],
		      BaryonField[DIINum][cellindex],BaryonField[HDINum][cellindex]);
	    }
	    if (TestProblemData.UseMetallicityField) {
	      fprintf(stderr,", Metal: %.2"ESYM,BaryonField[MetalNum][cellindex]);
	    }
	    fprintf(stderr,"\n");
	  }
	  else if (GloverChemistryModel) {
	    int GCM = TestProblemData.GloverChemistryModel;  // purely for convenience
	    fprintf(stderr,"Cell fraction: %.4"FSYM", Density: %.2"ESYM, cellfraction,BaryonField[DensNum][cellindex]);
	    fprintf(stderr,", HI: %.2"ESYM", HII: %.2"ESYM", H2I: %.2"ESYM,BaryonField[HINum][cellindex],BaryonField[HIINum][cellindex],
		    BaryonField[H2INum][cellindex]);
	    if( (GCM==1) || (GCM==2) || (GCM==3) || (GCM==7) ){
	      fprintf(stderr,", HeI: %.2"ESYM", HeII: %.2"ESYM", HeIII: %.2"ESYM,BaryonField[HeINum][cellindex],
		      BaryonField[HeIINum][cellindex],BaryonField[HeIIINum][cellindex]);
	      fprintf(stderr,", DI: %.2"ESYM", DII: %.2"ESYM", HDI: %.2"ESYM,BaryonField[DINum][cellindex],
		      BaryonField[DIINum][cellindex],BaryonField[HDINum][cellindex]);
	    }
	    fprintf(stderr,"\n");
	  }
	} else { // do this when r > dr

	  if(RadiatingShockUseDensityFluctuations){
	    therandominteger = mt_random();

	    therandomfraction =  float((therandominteger%32768)) /  32768.0;

	    if(therandomfraction < 0.0 || therandomfraction > 1.0){
	      fprintf(stderr,"yarr!  random number generator went bad!  %e\n",therandomfraction);
	      return FAIL;
	    }
	      
	    BaryonField[DensNum][cellindex] = BaryonField[DensNum][cellindex] * 
	      (1.0 +  RadiatingShockDensityFluctuationLevel*(therandomfraction - 0.5) );

	  } // if(RadiatingShockUseDensityFluctuations)

	    /* If we have primordial species turned on, leave their fractional abundances alone, but
	       scale them appropriately. */

	  if(TestProblemData.MultiSpecies > 0){

	    BaryonField[HINum][cellindex] *= (1.0 +  RadiatingShockDensityFluctuationLevel*(therandomfraction - 0.5) );
	    BaryonField[HIINum][cellindex] *= (1.0 +  RadiatingShockDensityFluctuationLevel*(therandomfraction - 0.5) );
	    BaryonField[HeINum][cellindex]  *= (1.0 +  RadiatingShockDensityFluctuationLevel*(therandomfraction - 0.5) );
	    BaryonField[HeIINum][cellindex] *= (1.0 +  RadiatingShockDensityFluctuationLevel*(therandomfraction - 0.5) );
	    BaryonField[HeIIINum][cellindex] *= (1.0 +  RadiatingShockDensityFluctuationLevel*(therandomfraction - 0.5) );
	    BaryonField[DeNum][cellindex] *= (1.0 +  RadiatingShockDensityFluctuationLevel*(therandomfraction - 0.5) );

	    if(TestProblemData.MultiSpecies > 1){
	      BaryonField[HMNum][cellindex] *= (1.0 +  RadiatingShockDensityFluctuationLevel*(therandomfraction - 0.5) );
	      BaryonField[H2INum][cellindex] *= (1.0 +  RadiatingShockDensityFluctuationLevel*(therandomfraction - 0.5) );
	      BaryonField[H2IINum][cellindex] *= (1.0 +  RadiatingShockDensityFluctuationLevel*(therandomfraction - 0.5) );
	    }

	    if(TestProblemData.MultiSpecies > 2){
	      BaryonField[DINum ][cellindex] *= (1.0 +  RadiatingShockDensityFluctuationLevel*(therandomfraction - 0.5) );
	      BaryonField[DIINum][cellindex] *= (1.0 +  RadiatingShockDensityFluctuationLevel*(therandomfraction - 0.5) );
	      BaryonField[HDINum][cellindex] *= (1.0 +  RadiatingShockDensityFluctuationLevel*(therandomfraction - 0.5) );
	    }

	  }  // if(TestProblemData.MultiSpecies > 0)

	  /* But if we have metals or simon glover's
	       cooling turned on, set these to tiny_number times the density field! */

	  if(TestProblemData.UseMetallicityField){
	    BaryonField[MetalNum][cellindex] = tiny_number * BaryonField[DensNum][cellindex];

	    if(TestProblemData.MultiMetals){
	      BaryonField[MetalNum+1][cellindex] = BaryonField[MetalNum+2][cellindex] = tiny_number * BaryonField[DensNum][cellindex];
		
	    }
	  } // if(TestProblemData.UseMetallicityField)

	  if(TestProblemData.GloverChemistryModel){
	    int GCM = TestProblemData.GloverChemistryModel;  // purely for convenience

	    if( (GCM==3) || (GCM==5) || (GCM==7) )
	      BaryonField[COINum  ][cellindex] = tiny_number * BaryonField[DensNum][cellindex];
	      
	    if( (GCM==2) || (GCM==3) || (GCM==7) ){
	      BaryonField[CINum ][cellindex] = BaryonField[CIINum][cellindex] = BaryonField[OINum ][cellindex] = 
		BaryonField[OIINum][cellindex] = tiny_number * BaryonField[DensNum][cellindex];
	    }

	    if( (GCM==2) || (GCM==3) ){
	      BaryonField[SiINum  ][cellindex] = BaryonField[SiIINum ][cellindex] = 
		BaryonField[SiIIINum][cellindex] = tiny_number * BaryonField[DensNum][cellindex];
	    }

	    if( (GCM==3) || (GCM==7) ){
	      BaryonField[CHINum  ][cellindex] = BaryonField[CH2INum ][cellindex] = 
		BaryonField[CH3IINum][cellindex] = BaryonField[C2INum  ][cellindex] = 
		BaryonField[HCOIINum][cellindex] = BaryonField[OHINum  ][cellindex] = 
		BaryonField[H2OINum ][cellindex] = BaryonField[O2INum  ][cellindex] = 
		tiny_number * BaryonField[DensNum][cellindex];
	    }

	  } // if(TestProblemData.GloverChemistryModel)


	}   // end of else block [if (r < dr)]

      } // for(i=0...)


  // clean up!
  if(RadiatingShockUseSedovProfile){
    delete [] sedovradius;
    delete [] sedovdensity;
    delete [] sedovpressure;
    delete [] sedovvelocity;
  }

  return SUCCESS;
}

// inputs:  
// cellx,y,z are cell center
// shocks,y,z are shock center
// dx is cell width (assume cell is always the same width and height
// radius is radius of shock
// GridRank is 1,2,3 (grid rank)
float cell_fraction(FLOAT cellx, FLOAT celly, FLOAT cellz, FLOAT shockx, FLOAT shocky, FLOAT shockz, FLOAT dx, FLOAT radius, int GridRank){

  float fraction=0.0, weighting=0.0;;
  FLOAT cellleftedge_x, cellleftedge_y, cellleftedge_z, thisx, thisy, thisz, thisradius;
  int i,j,k;

  

  cellleftedge_x = cellx - 0.5*dx;

  if(GridRank > 1)
    cellleftedge_y = celly - 0.5*dx;

  if(GridRank > 2)
    cellleftedge_z = cellz - 0.5*dx;

  for(k=0; k < 20; k++)
    for(j=0; j < 20; j++)
      for(i=0; i < 20; i++){

	thisx = cellleftedge_x + dx * (0.5 + float(i) )/20.0;  // center of subcell
	if (GridRank > 1)
	  thisy = cellleftedge_y + dx * (0.5 + float(j) )/20.0;  
	if (GridRank > 2)
	  thisz = cellleftedge_z + dx * (0.5 + float(k) )/20.0;  

	// it's REALLY r^2 right now
	thisradius = POW(thisx-shockx, 2.0) +
	  POW(thisy-shocky, 2.0);
      
	if(GridRank > 2)
	  thisradius += POW(thisz-shockz, 2.0 );
  
	thisradius = sqrt(thisradius);  // ok, now it's just r
	
	if(thisradius <= radius){
	  fraction += 1.0;
	}

	weighting += 1.0;

      }  // for(i=0; i < 20; i++)

  fraction /= weighting;

  //  fprintf(stderr,"fraction is %e\n",fraction);

  return fraction;

}


/* solution given by L.I. Sedov, 1956.  Originally implemented by
   Colin McNally, McMaster University.  Translated by BWO, July 2009.

  Nbins = # of solution points;
  radius = array of radii to solve at (cm)
  density = output array of densities (grams/cc)
  pressure = output array of pressures (ergs/vol)
  velocity = output array of velocities (cm/s)
  timesolve = time after explosion that we want the solution
  gamma = gamma value of gas
  beta = normalization for radius, velocity (depends on gamma)
  explosion_energy = energy in explosion (ergs)
  p_ambient, rho_ambient = ambient pressure, density of medium        */
void set_analytic_sedov(int Nbins, double *radius, double *density,
			  double *pressure, double *velocity,
			  double timesolve, double gamma, double beta,
			  double explosion_energy, double p_ambient, double rho_ambient)
{

  int i;
  double gamp1,gamm1,gam7,k, R2, u1, v2, p2;
  double alpha1, alpha2, alpha3, alpha4, alpha5;

  gamp1 = gamma + 1.0;
  gamm1 = gamma - 1.0;
  gam7 = 7.0 - gamma;
  k = gamp1 / gamm1;

  // shock radii, velocity, pressure
  R2 = beta * POW(explosion_energy / rho_ambient, 0.2) * POW(timesolve, 0.4);
  u1 = 0.4 * beta * POW( explosion_energy/rho_ambient, 0.2) / POW(timesolve, 0.6);
  v2 = 2.0 / gamp1 * u1;
  p2 = 2.0 * rho_ambient * u1 * u1  / gamp1;

  // exponents for self-similar solution
  alpha2 = (1.0 - gamma) / ( 2.0*(gamma - 1.0) + 3.0);
  alpha1 = 5.0 * gamma / ( 2.0 + 3.0*gamm1)*(6.0*(2.0-gamma)/(gamma*5.0*5.0)-alpha2);
  alpha3 = 3.0 / (2.0*gamm1 + 3.0);
  alpha4 = alpha1*5.0 / (2.0 - gamma);
  alpha5 = 2.0 / (gamma - 2.0);

  double xi, VV;

  for(i=1; i<Nbins;i++){
    xi = radius[i] / R2;  // dimensionless

    if(xi <= 1.0){
      
      VV = compute_sedov_v(xi,gamma,alpha1,alpha2);

      density[i] = rho_ambient * k * POW(k*(5.0*gamma/2.0*VV-1.0), alpha3)  
	* POW(k*(1.0 - 5.0*VV/2.0), alpha5)
	* POW(5.0 * gamp1 / (5.0*gamp1 - 2.0*(2.0+3.0*gamm1))*(1.0-(2.0+3.0*gamm1)/2.0*VV),  alpha4);

      pressure[i] = p2* POW( (5.0*gamp1*VV/4.0), 6.0/5.0) * POW( k*(1.0-5.0*VV/2.0), alpha5 + 1.0)
	* POW( 5.0*gamp1 / (5.0*gamp1 - 2.0*(2.0+3.0*gamm1))*1 - (2.0+3.0*gamm1)*VV/2.0, alpha4 - 2.0*alpha1);

      velocity[i] = 5.0 * gamp1 * VV * xi * v2/4.0;

    } else {
      density[i] = rho_ambient;
      pressure[i] = p_ambient;
      velocity[i] = 0.0;
    }

  } // for(i=1...)

  return;

}

double compute_sedov_v(double xi, double gamma, double alpha1, double alpha2){

  double V, tolerance = 1.0e-2, VL, VR, Vmid, xiL, xiR, ximid;
  int n_iter, n_iter_max = 100000, done;

  if(xi <= 1.0e-6){

    V = 2.0 / (5.0*gamma);

  } else {

    VL = 2.0 / (5.0*gamma);
    VR = 4.0 / (5.0*(gamma+1.0));

    xiL = sedov_vfunc(VL, gamma, alpha1, alpha2);
    xiR = sedov_vfunc(VR, gamma, alpha1, alpha2);

    n_iter = 1;

    done=0;

    while(!done){

      Vmid = 0.5 * (VL + VR);
      ximid = sedov_vfunc(Vmid, gamma, alpha1, alpha2);

      if( (fabs(ximid-xi) <= tolerance) || (n_iter >= n_iter_max) ) done = 1;

      if(!done){

	n_iter++;

	if( (xi >= xiL) && (xi <= ximid ) ){
	  VR = Vmid;
	} else if ( (xi >= ximid) && (xi <= xiR) ){
	  VL = Vmid;
	} else {
	  done = 1;
	}

	xiL = sedov_vfunc(VL, gamma, alpha1, alpha2);
	xiR = sedov_vfunc(VR, gamma, alpha1, alpha2);

      }  // if(!done)

    } // while(!done)

    if(n_iter > n_iter_max){
      printf("compute_sedov_v: exceeded max iters for xi = %e\n",xi);
      printf("xiL = %e   xiR = %e\n", xiL, xiR);
      return FAIL;
    }

    V = Vmid;
    
  } // big if/else statement

  return V;
}


double sedov_vfunc(double V, double gamma, double alpha1, double alpha2){

  double tmp, gamp1, gamm1, gam7;
  gamp1 = gamma + 1.0;
  gamm1 = gamma - 1.0;
  gam7 = 7.0 - gamma;

  tmp = POW( (5.0*gamp1*V/4.0) , -2.0/5.0 ) * POW( gamp1/gamm1*(5.0*gamma*V/2.0-1.0), -1.0*alpha2);
  tmp = tmp * POW( (5.0*gamp1 / (5.0*gamp1 - 2.0*(2.0+3.0*gamm1)) *( 1.0 - (2.0+3.0*gamm1)*V/2.0)), -1.0*alpha1);

  return tmp;

}
