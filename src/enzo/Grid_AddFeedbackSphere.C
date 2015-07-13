/***********************************************************************
/
/  GRID: ADD SPHERICAL STAR PARTICLE FEEDBACK TO CELLS
/
/  written by: John Wise
/  date:       September, 2005
/  modified1: Ji-hoon Kim to include MBH_THERMAL feedback
/             July, 2009
/  modified2: Ji-hoon Kim to include MBH_JETS feedback
/             November, 2009
/
/  PURPOSE:
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "CosmologyParameters.h"

#ifdef CONFIG_BFLOAT_4
#define TOLERANCE 1e-06
#endif
#ifdef CONFIG_BFLOAT_8
#define TOLERANCE 1e-12
#endif
#ifdef CONFIG_BFLOAT_16
#define TOLERANCE 1e-15
#endif

#define MAX_TEMPERATURE 1e8

int FindField(int field, int farray[], int numfields);

int grid::AddFeedbackSphere(Star *cstar, int level, float radius, float DensityUnits, 
			    float LengthUnits, float VelocityUnits, 
			    float TemperatureUnits, float TimeUnits, double EjectaDensity, 
			    double EjectaMetalDensity, double EjectaThermalEnergy, 
			    int &CellsModified)
{

  const float WhalenMaxVelocity = 35;		// km/s
  const double Msun = 1.989e33;
  const double c = 3.0e10;

  int dim, i, j, k, index;
  int sx, sy, sz;
  FLOAT delx, dely, delz, radius2, Radius, DomainWidth[MAX_DIMENSION];
  float coef, speed, maxVelocity;
  float OldDensity;
  float r1, norm, ramp, factor, newGE, fh;
  double increase;

  if (MyProcessorNumber != ProcessorNumber)
    return SUCCESS;

  /* Check if sphere overlaps with this grid */

  for (dim = 0; dim < GridRank; dim++)
    if (cstar->pos[dim] - radius > GridRightEdge[dim] ||
	cstar->pos[dim] + radius < GridLeftEdge[dim])
      return SUCCESS;

  for (dim = 0; dim < GridRank; dim++)
    DomainWidth[dim] = DomainRightEdge[dim] - DomainLeftEdge[dim];

  /* Find fields: density, total energy, velocity1-3. */

  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num;
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				       Vel3Num, TENum) == FAIL) {
        ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }
  
  /* Find Multi-species fields. */

  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
  if (MultiSpecies) 
    if (this->IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, 
				    HeIIINum, HMNum, H2INum, H2IINum, DINum, 
				    DIINum, HDINum) == FAIL) {
        ENZO_FAIL("Error in grid->IdentifySpeciesFields.");
    }

  /* Find Metallicity or SNColour field and set flag. */

  int SNColourNum, MetalNum, Metal2Num, MBHColourNum, Galaxy1ColourNum, 
    Galaxy2ColourNum, MetalIaNum, MetalIINum;
  int MetallicityField = FALSE;

  if (this->IdentifyColourFields(SNColourNum, Metal2Num, MetalIaNum, 
				 MetalIINum, MBHColourNum, Galaxy1ColourNum, 
				 Galaxy2ColourNum) == FAIL)
    ENZO_FAIL("Error in grid->IdentifyColourFields.\n");

  MetalNum = max(Metal2Num, SNColourNum);
  MetallicityField = (MetalNum > 0) ? TRUE : FALSE;
  if (MetalNum > 0 && SNColourNum > 0 && cstar->type == PopIII)
    MetalNum = SNColourNum;

  float BubbleVolume = (4.0 * M_PI / 3.0) * radius * radius * radius;

  /***********************************************************************
                                SUPERNOVAE
  ************************************************************************/

  // Assume that the remnant is still in the free expansion stage and
  // hasn't had any radiative losses.  In this case, the ejecta will
  // be at 3/4 the radius of the shock front (see Ostriker & McKee
  // 1988 or Tenorio-Tagle 1996).

  //const float MetalRadius = 0.75;
  const float MetalRadius = 1.0;
  float ionizedFraction = 0.999;  // Assume supernova is ionized
  float maxGE, MetalRadius2, PrimordialDensity, metallicity, fhz, fhez;
  float outerRadius2;

  if (cstar->FeedbackFlag == SUPERNOVA || 
      cstar->FeedbackFlag == CONT_SUPERNOVA) {

  // Correct for exaggerated influence radius for pair-instability supernovae
    if (cstar->FeedbackFlag == SUPERNOVA)
      radius /= 1.0;

  // Correct if the volume with 27 cells is larger than the energy bubble volume
#ifdef UNUSED
  float BoxVolume = 27 * CellWidth[0][0] * CellWidth[0][0] * CellWidth[0][0];
  float BubbleVolume = (4.0 * M_PI / 3.0) * radius * radius * radius;
  //printf("BoxVolume = %lg, BubbleVolume = %lg\n", BoxVolume, BubbleVolume);
  if (BoxVolume > BubbleVolume) {
    //printf("Reducing ejecta density by %g\n", BubbleVolume / BoxVolume);
    EjectaDensity *= BubbleVolume / BoxVolume;
    EjectaMetalDensity *= BubbleVolume / BoxVolume;
    EjectaThermalEnergy *= BubbleVolume / BoxVolume;
  }
#endif
//  if (cstar->level > level) {
//    printf("Reducing ejecta density and energy by 10%% on "
//	   "level %"ISYM" to avoid crashing.\n", level);
//    EjectaDensity *= 0.1;
//    EjectaMetalDensity *= 0.1;
//    EjectaThermalEnergy *= 0.1;
//  }

  // Correct for smaller enrichment radius
  EjectaMetalDensity *= pow(MetalRadius, -3.0);
  PrimordialDensity = EjectaDensity - EjectaMetalDensity;
  fh = CoolData.HydrogenFractionByMass;
  MetalRadius2 = radius * radius * MetalRadius * MetalRadius;
  outerRadius2 = 1.2 * 1.2 * radius * radius;

    /* Remove mass from the star that will now be added to grids. 
       Also, because EjectaDensity will be added with zero net momentum, 
       increase the particle's velocity accordingly. - Ji-hoon Kim, Sep.2009 */

//    printf("grid::AFS: before: cstar->Mass = %lf\n", cstar->Mass); 
    if (cstar->FeedbackFlag != SUPERNOVA) {
      float old_mass = (float)(cstar->Mass);
      cstar->Mass -= EjectaDensity * DensityUnits * BubbleVolume * pow(LengthUnits,3.0) / Msun;  
      float frac = old_mass / cstar->Mass;
      cstar->vel[0] *= frac;
      cstar->vel[1] *= frac;
      cstar->vel[2] *= frac;
    } // ENDIF !Supernova

//    printf("grid::AFS: after : cstar->Mass = %lf\n", cstar->Mass); 
//    printf("grid::AFS: pos = %"FSYM" %"FSYM" %"FSYM"\n", 
//	   cstar->pos[0], cstar->pos[1], cstar->pos[2]);

    maxGE = MAX_TEMPERATURE / (TemperatureUnits * (Gamma-1.0) * 0.6);

    for (k = 0; k < GridDimension[2]; k++) {

      delz = CellLeftEdge[2][k] + 0.5*CellWidth[2][k] - cstar->pos[2];
      sz = sign(delz);
      delz = fabs(delz);
      delz = min(delz, DomainWidth[2]-delz);

      for (j = 0; j < GridDimension[1]; j++) {

	dely = CellLeftEdge[1][j] + 0.5*CellWidth[1][j] - cstar->pos[1];
	sy = sign(dely);
	dely = fabs(dely);
	dely = min(dely, DomainWidth[1]-dely);

	index = (k*GridDimension[1] + j)*GridDimension[0];
	for (i = 0; i < GridDimension[0]; i++, index++) {

	  delx = CellLeftEdge[0][i] + 0.5*CellWidth[0][i] - cstar->pos[0];
	  sx = sign(delx);
	  delx = fabs(delx);
	  delx = min(delx, DomainWidth[0]-delx);

	  radius2 = delx*delx + dely*dely + delz*delz;
	  if (radius2 <= outerRadius2) {

	    r1 = sqrt(radius2) / radius;
	    norm = 0.98;
	    ramp = norm*(0.5 - 0.5 * tanh(10.0*(r1-1.0)));
//	    ramp = min(max(1.0 - (r1 - 0.8)/0.4, 0.01), 1.0);

	    /* 1/1.2^3 factor to dilute the density since we're
	       depositing a uniform ejecta in a sphere of 1.2*radius
	       without a ramp.  The ramp is only applied to the
	       energy*density factor. */
	    factor = 0.578704;

	    OldDensity = BaryonField[DensNum][index];
	    BaryonField[DensNum][index] += factor*EjectaDensity;

	    /* Add total energies of spheres together, then divide by
	       density to get specific energy */

	    if (GENum >= 0 && DualEnergyFormalism) {

	      newGE = (OldDensity * BaryonField[GENum][index] +
		       ramp * factor * EjectaDensity * EjectaThermalEnergy) /
		BaryonField[DensNum][index];
	      newGE = min(newGE, maxGE);  
//	      newGE = ramp * EjectaThermalEnergy;

	      BaryonField[GENum][index] = newGE;
	      BaryonField[TENum][index] = newGE;

	      for (dim = 0; dim < GridRank; dim++)
		BaryonField[TENum][index] += 
		  0.5 * BaryonField[Vel1Num+dim][index] * 
		  BaryonField[Vel1Num+dim][index];

	    } else {

	      newGE = (OldDensity * BaryonField[TENum][index] +
		       ramp * factor * EjectaDensity * EjectaThermalEnergy) /
		BaryonField[DensNum][index];

	      newGE = min(newGE, maxGE);  
	      BaryonField[TENum][index] = newGE;

	    } //end if(GENum >= 0 && DualEnergyFormalism)

	    /* Update species and colour fields */

	    //increase = BaryonField[DensNum][index] / OldDensity;
	    if (MetallicityField == TRUE) {
	      if (radius2 <= MetalRadius2) {
		metallicity = (BaryonField[MetalNum][index] + EjectaMetalDensity) /
		  BaryonField[DensNum][index];
	      } else {
		metallicity = BaryonField[MetalNum][index] / BaryonField[DensNum][index];
	      }
	    } else
	      metallicity = 0;

	    fhz = fh * (1-metallicity);
	    fhez = (1-fh) * (1-metallicity);

	    if (MultiSpecies) {
	      BaryonField[DeNum][index] = (1-metallicity) *
		BaryonField[DensNum][index] * ionizedFraction;
	      BaryonField[HINum][index] = 
		BaryonField[DensNum][index] * fhz * (1-ionizedFraction);
	      BaryonField[HIINum][index] =
		BaryonField[DensNum][index] * fhz * ionizedFraction;
	      BaryonField[HeINum][index] =
		0.5*BaryonField[DensNum][index] * fhez * (1-ionizedFraction);
	      BaryonField[HeIINum][index] =
		0.5*BaryonField[DensNum][index] * fhez * (1-ionizedFraction);
	      BaryonField[HeIIINum][index] =
		BaryonField[DensNum][index] * fhez * ionizedFraction;
	    }
	    if (MultiSpecies > 1) {
	      BaryonField[HMNum][index] = tiny_number;// * BaryonField[DensNum][index];
	      BaryonField[H2INum][index] = 
		tiny_number;// * BaryonField[DensNum][index];
	      BaryonField[H2IINum][index] = 
		tiny_number;// * BaryonField[DensNum][index];
	    }
	    if (MultiSpecies > 2) {
	      BaryonField[DINum][index] = BaryonField[DensNum][index] * fh *
		CoolData.DeuteriumToHydrogenRatio * (1-ionizedFraction);
	      BaryonField[DIINum][index] = BaryonField[DensNum][index] * fh *
		CoolData.DeuteriumToHydrogenRatio * ionizedFraction;
	      BaryonField[HDINum][index] = 
		tiny_number * BaryonField[DensNum][index];
	    }

	    if (MetallicityField == TRUE)
	      BaryonField[MetalNum][index] = metallicity * BaryonField[DensNum][index];

	    CellsModified++;

	  } // END if inside radius
	}  // END i-direction
      }  // END j-direction
    }  // END k-direction

  }  // END Supernova

  /***********************************************************************
                                MBH_THERMAL
  ************************************************************************/

  // Similar to Supernova, but here we assume the followings:
  // EjectaDensity = 0.0
  // EjectaMetalDensity = 0.0
  // The unit of EjectaThermalEnergy = ergs/cm^3, not ergs/g

  if (cstar->FeedbackFlag == MBH_THERMAL) {

    maxGE = MAX_TEMPERATURE / (TemperatureUnits * (Gamma-1.0) * 0.6);

    for (k = 0; k < GridDimension[2]; k++) {

      delz = CellLeftEdge[2][k] + 0.5*CellWidth[2][k] - cstar->pos[2];
      sz = sign(delz);
      delz = fabs(delz);
      delz = min(delz, DomainWidth[2]-delz);

      for (j = 0; j < GridDimension[1]; j++) {

	dely = CellLeftEdge[1][j] + 0.5*CellWidth[1][j] - cstar->pos[1];
	sy = sign(dely);
	dely = fabs(dely);
	dely = min(dely, DomainWidth[1]-dely);

	index = (k*GridDimension[1] + j)*GridDimension[0];
	for (i = 0; i < GridDimension[0]; i++, index++) {

	  delx = CellLeftEdge[0][i] + 0.5*CellWidth[0][i] - cstar->pos[0];
	  sx = sign(delx);
	  delx = fabs(delx);
	  delx = min(delx, DomainWidth[0]-delx);

	  radius2 = delx*delx + dely*dely + delz*delz;
	  if (radius2 <= outerRadius2) {

	    r1 = sqrt(radius2) / radius;
	    norm = 0.98;
	    ramp = norm*(0.5 - 0.5 * tanh(10.0*(r1-1.0)));
//          ramp = min(max(1.0 - (r1 - 0.8)/0.4, 0.01), 1.0);

	    /* 1/1.2^3 factor to dilute the density since we're
	       depositing a uniform ejecta in a sphere of 1.2*radius
	       without a ramp.  The ramp is only applied to the
	       energy*density factor. */
	    factor = 0.578704;

	    OldDensity = BaryonField[DensNum][index];
	    BaryonField[DensNum][index] += factor*EjectaDensity;

	    /* Get specific energy */

	    if (GENum >= 0 && DualEnergyFormalism) {

	      /* When injected energy is uniform throughout the volume;
		 EjectaThermalEnergy in ergs/cm3 */
	      newGE = (OldDensity * BaryonField[GENum][index] +
		       ramp * factor * EjectaThermalEnergy) /
		BaryonField[DensNum][index];

#ifdef USE_ONE_OVER_RSQUARED
	      /* When injected energy is proportional to 1/R^2;
		 EjectaThermalEnergy in ergs/cm3/(1/cm^2) */
	      newGE = (OldDensity * BaryonField[GENum][index] +
		       ramp * factor * EjectaThermalEnergy * 
		       min(1.0/radius2, 1.0/(4.0*CellWidth[0][0]*CellWidth[0][0]))) /
		BaryonField[DensNum][index];
#endif

#ifdef CONSTANT_SPECIFIC
	      /* When injected energy is proportional to the cell mass;
		 EjectaThermalEnergy in ergs/g = cm2/s2 */
	      newGE = (BaryonField[GENum][index] + EjectaThermalEnergy) *
		OldDensity / BaryonField[DensNum][index];
#endif

	      newGE = min(newGE, maxGE);  

	      BaryonField[GENum][index] = newGE;
	      BaryonField[TENum][index] = newGE;

	      for (dim = 0; dim < GridRank; dim++)
		BaryonField[TENum][index] += 
		  0.5 * BaryonField[Vel1Num+dim][index] * 
		  BaryonField[Vel1Num+dim][index];

	    } else {

	      newGE = (OldDensity * BaryonField[TENum][index] +
		       ramp * factor * EjectaThermalEnergy) /
		BaryonField[DensNum][index];

#ifdef USE_ONE_OVER_RSQUARED
	      newGE = (OldDensity * BaryonField[TENum][index] +
		       ramp * factor * EjectaThermalEnergy * 
		       min(1.0/radius2, 1.0/(4.0*CellWidth[0][0]*CellWidth[0][0]))) /
		BaryonField[DensNum][index];
#endif
	      
#ifdef CONSTANT_SPECIFIC
	      newGE = (BaryonField[TENum][index] + EjectaThermalEnergy) *
		OldDensity / BaryonField[DensNum][index];
#endif

//	      printf("grid::AddFS: rho= %"GSYM"->%"GSYM", TE= %"GSYM"->%"GSYM", drho= %"GSYM", dE= %"GSYM"\n", 
//		     OldDensity, BaryonField[DensNum][index], 
//		     BaryonField[TENum][index], newGE, EjectaDensity, EjectaThermalEnergy * 1/radius2); 

	      newGE = min(newGE, maxGE);  
	      BaryonField[TENum][index] = newGE;

	    } //end if(GENum >= 0 && DualEnergyFormalism)

	    /* Update species and colour fields */

	    //increase = BaryonField[DensNum][index] / OldDensity;
	    if (MetallicityField == TRUE) {
	      if (radius2 <= MetalRadius2) {
		metallicity = (BaryonField[MetalNum][index] + EjectaMetalDensity) /
		  BaryonField[DensNum][index];
	      } else {
		metallicity = BaryonField[MetalNum][index] / BaryonField[DensNum][index];
	      }
	    } else
	      metallicity = 0;

	    fhz = fh * (1-metallicity);
	    fhez = (1-fh) * (1-metallicity);

	    if (MultiSpecies) {
	      BaryonField[DeNum][index] = (1-metallicity) *
		BaryonField[DensNum][index] * ionizedFraction;
	      BaryonField[HINum][index] = 
		BaryonField[DensNum][index] * fhz * (1-ionizedFraction);
	      BaryonField[HIINum][index] =
		BaryonField[DensNum][index] * fhz * ionizedFraction;
	      BaryonField[HeINum][index] =
		0.5*BaryonField[DensNum][index] * fhez * (1-ionizedFraction);
	      BaryonField[HeIINum][index] =
		0.5*BaryonField[DensNum][index] * fhez * (1-ionizedFraction);
	      BaryonField[HeIIINum][index] =
		BaryonField[DensNum][index] * fhez * ionizedFraction;
	    }
	    if (MultiSpecies > 1) {
	      BaryonField[HMNum][index] = tiny_number * BaryonField[DensNum][index];
	      BaryonField[H2INum][index] = 
		tiny_number * BaryonField[DensNum][index];
	      BaryonField[H2IINum][index] = 
		tiny_number * BaryonField[DensNum][index];
	    }
	    if (MultiSpecies > 2) {
	      BaryonField[DINum][index] = BaryonField[DensNum][index] * fh *
		CoolData.DeuteriumToHydrogenRatio * (1-ionizedFraction);
	      BaryonField[DIINum][index] = BaryonField[DensNum][index] * fh *
		CoolData.DeuteriumToHydrogenRatio * ionizedFraction;
	      BaryonField[HDINum][index] = 
		tiny_number * BaryonField[DensNum][index];
	    }

	    if (MetallicityField == TRUE)
	      BaryonField[MetalNum][index] = metallicity * BaryonField[DensNum][index];

	    /* MBHColour injected */
	    if (MBHColourNum > 0)
	      BaryonField[MBHColourNum][index] += factor*EjectaDensity;

	    CellsModified++;

	  } // END if inside radius
	}  // END i-direction
      }  // END j-direction
    }  // END k-direction

  }  // END MBH_THERMAL

  /***********************************************************************
                                 MBH_JETS
  ************************************************************************/

  // Inject bipolar jets along the direction of the angular momentum 
  // vector L of the MBH particle (angular momentum accreted thus far)
  // or along the z-axis  - Ji-hoon Kim, Nov.2009

#define MAX_SUPERCELL_NUMBER 1000
  int SUPERCELL = 2; //2 for supercell of 5 cells wide = 5^3  
  int ind_cell_inside[MAX_SUPERCELL_NUMBER], ind_cell_edge[MAX_SUPERCELL_NUMBER];
  float nx_cell_edge[MAX_SUPERCELL_NUMBER], ny_cell_edge[MAX_SUPERCELL_NUMBER], 
    nz_cell_edge[MAX_SUPERCELL_NUMBER];
  int n_cell_inside = 0, n_cell_edge = 0, ibuff = NumberOfGhostZones;
  int ii, jj, kk, r_s, ic, sign;
  float m_cell_inside = 0.0, metal_cell_inside = 0.0, colour_cell_inside = 0.0, 
    metallicity_inside = 0.0, colour_inside = 0.0, rho_inside, rho_metal_inside, rho_colour_inside;
  float m_cell_edge = 0.0, metal_cell_edge = 0.0, colour_cell_edge = 0.0, 
    metallicity_edge = 0.0, colour_edge = 0.0, rho_jet, rho_metal_jet, rho_colour_jet;
  float L_x, L_y, L_z, L_s, nx_L = 0.0, ny_L = 0.0, nz_L = 0.0, costheta = cos(3.1415926/3.9);
  float EjectaMass, EjectaMetalMass, MBHJetsVelocity;

  if (cstar->FeedbackFlag == MBH_JETS) {

    /* Calculate star indicies */

    float CellWidthTemp = float(CellWidth[0][0]);
    i = (int)((cstar->pos[0] - CellLeftEdge[0][0]) / CellWidthTemp);
    j = (int)((cstar->pos[1] - CellLeftEdge[1][0]) / CellWidthTemp);
    k = (int)((cstar->pos[2] - CellLeftEdge[2][0]) / CellWidthTemp);
    
    /* Note that we need to inject feedback only for the finest grid the MBH belongs to */

    if (i < ibuff || i > GridDimension[0]-ibuff-1 ||
	j < ibuff || j > GridDimension[1]-ibuff-1 || 
	k < ibuff || k > GridDimension[2]-ibuff-1 ||
	cstar->ReturnCurrentGrid() == NULL ||
	cstar->level > level) {
//    fprintf(stdout, "grid::AddFS: MBH_JETS - MBH doesn't belong to this grid.\n"); 
      return SUCCESS;
    }

    /* Calculate mass that has accumulated thus far; since EjectaDensity is calculated 
       for isotropic MBH_THERMAL, we use it only to compute EjectaMass, not apply it directly. */

    cstar->NotEjectedMass += EjectaDensity * DensityUnits * BubbleVolume * pow(LengthUnits,3.0) / Msun;  

//    fprintf(stdout, "d1, d2, d3, i, j, k = %d %d %d / %d %d %d\n", 
//	    GridDimension[0], GridDimension[1], GridDimension[2], i,j,k);

    /* If NotEjectedMass is still smaller than the threshold, return */

    if (cstar->NotEjectedMass <= MBHFeedbackJetsThresholdMass) {
      fprintf(stdout, "grid::AddFS: MBH_JETS - accumulated mass (%g Ms) not passed threshold.\n",
	      cstar->NotEjectedMass);       
      return SUCCESS;
    }

    /* If the current grid cannot contain the whole supercell, return */

    if (i < ibuff+SUPERCELL || i > GridDimension[0]-ibuff-SUPERCELL-1 || 
	j < ibuff+SUPERCELL || j > GridDimension[1]-ibuff-SUPERCELL-1 ||
	k < ibuff+SUPERCELL || k > GridDimension[2]-ibuff-SUPERCELL-1) {
      fprintf(stdout, "grid::AddFS: MBH_JETS - supercell not contained; accumulated mass (%g Ms).\n",
	      cstar->NotEjectedMass); 
      
      // if the supercell issue hasn't allowed the jet injection for too long,
      // issue an warning signal and output the current hierarchy at CheckForOutput
      if (cstar->NotEjectedMass > 2.0 * MBHFeedbackJetsThresholdMass) { 
	fprintf(stdout, "grid::AddFS: MBH_JETS - jets haven't been ejected for too long!\n");
	OutputWhenJetsHaveNotEjected = TRUE;  
      } 
      
      // otherwise, just proceed and do it later
      return SUCCESS;
    }
    
    /* Find ejecta mass */

    EjectaMass = cstar->NotEjectedMass * Msun / DensityUnits / pow(LengthUnits,3.0); 
    EjectaMetalMass = EjectaMass * MBHFeedbackMetalYield;
    OutputWhenJetsHaveNotEjected = FALSE; 

    /* Find the directional vector n_L of angular momentum accreted thus far */

    L_x = cstar->accreted_angmom[0];
    L_y = cstar->accreted_angmom[1];
    L_z = cstar->accreted_angmom[2]; 
    L_s = sqrt(pow(L_x,2) + pow(L_y,2) + pow(L_z,2));
    nx_L = L_x/L_s;  //normalized directional vector
    ny_L = L_y/L_s;
    nz_L = L_z/L_s;

    // if MBHFeedback = 3, direct the jets always along the z-axis
    if (MBHFeedback == 3) {
      nx_L = 0.0;
      ny_L = 0.0;
      nz_L = 1.0;
    }

    // if MBHFeedback = 4 or 5, provide the jet direction with a random noise
    if (MBHFeedback == 4 || MBHFeedback == 5) {

      float theta, phi, MaximumNoiseAngle;
      if (MBHFeedback == 4) MaximumNoiseAngle = 10.0;
      if (MBHFeedback == 5) MaximumNoiseAngle = 90.0; //launching in random direction

      //find angles in spherical coordinate; n_L's are already normalized
      theta = atan(ny_L/nx_L);  
      phi   = acos(nz_L);
      printf("before: n_L = (%g, %g, %g) with theta,phi = (%g, %g)\n", nx_L, ny_L, nz_L, theta, phi);

      //add random noise to theta and phi
      srand(time(NULL));
      theta += MaximumNoiseAngle * M_PI / 180.0 * ((2.0*(float)rand()/((float)(RAND_MAX)+(float)(1))) - 1.0);  
      phi   += MaximumNoiseAngle * M_PI / 180.0 * ((2.0*(float)rand()/((float)(RAND_MAX)+(float)(1))) - 1.0);

      //get back to Cartesian coordinate; some tricks needed to preserve the signs of nx_L and ny_L
      nx_L = sign(nx_L) * fabs(cos(theta))*sin(phi);  
      ny_L = sign(ny_L) * fabs(sin(theta))*sin(phi);
      nz_L = cos(phi);
      printf("after : n_L = (%g, %g, %g) with theta,phi = (%g, %g) and random angle example = %g deg\n", 
	     nx_L, ny_L, nz_L, theta, phi,
	     MaximumNoiseAngle * ((2.0*(float)rand()/((float)(RAND_MAX)+(float)(1))) - 1.0));
    }

    /* Loop over the supercell around the MBH particle (5 * 5 * 5 = 125 cells, 
       but only the edges), and record the cells eligible for jet injection */

    for (kk = -SUPERCELL; kk <= SUPERCELL; kk++) {
      for (jj = -SUPERCELL; jj <= SUPERCELL; jj++) {
	for (ii = -SUPERCELL; ii <= SUPERCELL; ii++) {

	  if (fabs(ii) != SUPERCELL && fabs(jj) != SUPERCELL && fabs(kk) != SUPERCELL) {  //if not on edges

	    ind_cell_inside[n_cell_inside] = i+ii+(j+jj+(k+kk)*GridDimension[1])*GridDimension[0];
	    m_cell_inside += BaryonField[DensNum][ind_cell_inside[n_cell_inside]] * 
	      pow(CellWidth[0][0], 3);
	    if (MetallicityField == TRUE) 
	      metal_cell_inside += BaryonField[MetalNum][ind_cell_inside[n_cell_inside]] * 
		pow(CellWidth[0][0], 3);
	    if (MBHColourNum > 0)
	      colour_cell_inside += BaryonField[MBHColourNum][ind_cell_inside[n_cell_inside]] * 
		pow(CellWidth[0][0], 3);
	    n_cell_inside++;
	    
	  } else {  //if on edges

	    r_s = sqrt(pow(ii,2) + pow(jj,2) + pow(kk,2));	    
	    
	    if (fabs((ii*nx_L + jj*ny_L + kk*nz_L)/r_s) > costheta) { 

	      ind_cell_edge[n_cell_edge] = i+ii+(j+jj+(k+kk)*GridDimension[1])*GridDimension[0];
	      nx_cell_edge[n_cell_edge]  = ii / r_s;  //directional vector
	      ny_cell_edge[n_cell_edge]  = jj / r_s;
	      nz_cell_edge[n_cell_edge]  = kk / r_s;
	      m_cell_edge += BaryonField[DensNum][ind_cell_edge[n_cell_edge]] * 
		pow(CellWidth[0][0], 3);
	      if (MetallicityField == TRUE) 
		metal_cell_edge += BaryonField[MetalNum][ind_cell_edge[n_cell_edge]] * 
		  pow(CellWidth[0][0], 3);
	      if (MBHColourNum > 0) 
		colour_cell_edge += BaryonField[MBHColourNum][ind_cell_edge[n_cell_edge]] * 
		  pow(CellWidth[0][0], 3);
	      n_cell_edge++;
	    } 

	  }  

	}  // END ii-direction
      }  // END jj-direction
    }  // END kk-direction

//    printf("EjectaM in Msun = %g, EjectaM = %g, EjectaMetalM = %g, m_cell_edge = %g, n_cell_edge = %d\n",
//	   cstar->NotEjectedMass, EjectaMass, EjectaMetalMass, m_cell_edge, n_cell_edge); 

    /* Calculate the jet density */

    rho_jet = (m_cell_edge + EjectaMass) / 
      (n_cell_edge * pow(CellWidth[0][0], 3));

    if (MetallicityField == TRUE) {
      rho_metal_jet = (metal_cell_edge + EjectaMetalMass) / 
	(n_cell_edge * pow(CellWidth[0][0], 3));
      metallicity_edge = rho_metal_jet / rho_jet;
    } else
      metallicity_edge = 0.0;

    if (MBHColourNum > 0) {
      rho_colour_jet = (colour_cell_edge + EjectaMass) / 
	(n_cell_edge * pow(CellWidth[0][0], 3));
      colour_edge = rho_colour_jet / rho_jet;
    } else
      colour_edge = 0.0;
      
//    printf("rho_jet_prev = %g\n", m_cell_edge / (n_cell_edge * pow(CellWidth[0][0], 3)));
//    printf("rho_jet =%g, rho_metal_jet = %g, rho_colour_jet = %g, metallicity_edge = %g\n", 
//	   rho_jet, rho_metal_jet, rho_colour_jet, metallicity_edge); 

    /* Calculate MBHJetsVelocity using energy conservation below:

       MBHFeedbackEnergyCoupling * MBHFeedbackRadiativeEfficiency * Mdot * c^2 
       = 0.5 * MBHFeedbackMassEjectionFraction * Mdot * (MBHJetsVelocity)^2                
       
       Note that EjectaThermalEnergy is never used; MBHFeedbackEnergyCoupling 
       should now be calculated considering gravitational redshift (Kim et al. 2010) */

    MBHJetsVelocity = c * sqrt( 2 * MBHFeedbackEnergyCoupling * MBHFeedbackRadiativeEfficiency 
				      / MBHFeedbackMassEjectionFraction ) / VelocityUnits;

    if (MBHJetsVelocity * VelocityUnits > 0.99*c) {
      ENZO_VFAIL("grid::AddFS: MBHJetsVelocity is ultra-relativistic! (%g/ %g/ %g/ %g c)\n",
		 MBHFeedbackEnergyCoupling, MBHFeedbackRadiativeEfficiency, 
		 MBHFeedbackMassEjectionFraction, MBHJetsVelocity * VelocityUnits / c);
    }

    /* Finally, add the jet feedback at the edges (outer part of the supercell) */

    for (ic = 0; ic < n_cell_edge; ic++) {

      index = ind_cell_edge[ic];

      /* Update velocities and TE; note that we now have kinetic (jet) energy added, so 
	 for DualEnergyFormalism = 0 you don't have to update any energy field */

      sign = sign(nx_cell_edge[ic]*nx_L + ny_cell_edge[ic]*ny_L + nz_cell_edge[ic]*nz_L);

      if (GENum >= 0 && DualEnergyFormalism) 
	for (dim = 0; dim < GridRank; dim++)
	  BaryonField[TENum][index] -= 
	    0.5 * BaryonField[Vel1Num+dim][index] * 
	    BaryonField[Vel1Num+dim][index];

      /* Calculate grid velocity: the actual veloctiy injected in supercell edges.
	 This is different from MBHJetsVelocity because it is the mass-weighted average 
	 between MBHJetsVelocity and the original veloctiy in grid */  
      
      BaryonField[Vel1Num][index] = (BaryonField[DensNum][index] * BaryonField[Vel1Num][index] +
				     EjectaMass / (n_cell_edge*pow(CellWidth[0][0], 3)) *
				     sign * nx_L * MBHJetsVelocity) / 
	                            (BaryonField[DensNum][index] + 
				     EjectaMass / (n_cell_edge*pow(CellWidth[0][0], 3)));
      BaryonField[Vel2Num][index] = (BaryonField[DensNum][index] * BaryonField[Vel2Num][index] +
				     EjectaMass / (n_cell_edge*pow(CellWidth[0][0], 3)) *
				     sign * ny_L * MBHJetsVelocity) / 
	                            (BaryonField[DensNum][index] + 
				     EjectaMass / (n_cell_edge*pow(CellWidth[0][0], 3)));
      BaryonField[Vel3Num][index] = (BaryonField[DensNum][index] * BaryonField[Vel3Num][index] +
				     EjectaMass / (n_cell_edge*pow(CellWidth[0][0], 3)) *
				     sign * nz_L * MBHJetsVelocity) / 
	                            (BaryonField[DensNum][index] + 
				     EjectaMass / (n_cell_edge*pow(CellWidth[0][0], 3)));
      
      if (GENum >= 0 && DualEnergyFormalism) 
	for (dim = 0; dim < GridRank; dim++)
	  BaryonField[TENum][index] += 
	    0.5 * BaryonField[Vel1Num+dim][index] * 
	    BaryonField[Vel1Num+dim][index];

      /* Update density, species and colour fields */

      OldDensity = BaryonField[DensNum][index];
      BaryonField[DensNum][index] = rho_jet;
      increase = rho_jet / OldDensity;
//      printf("grid::AFS: increase = %lf\n", increase); 
      
      if (MultiSpecies) {
	BaryonField[DeNum][index] *= increase;
	BaryonField[HINum][index] *= increase;
	BaryonField[HIINum][index] *= increase;
	BaryonField[HeINum][index] *= increase;
	BaryonField[HeIINum][index] *= increase;
	BaryonField[HeIIINum][index] *= increase;
      }
      if (MultiSpecies > 1) {
	BaryonField[HMNum][index] *= increase;
	BaryonField[H2INum][index] *= increase;
	BaryonField[H2IINum][index] *= increase;
      }
      if (MultiSpecies > 2) {
	BaryonField[DINum][index] *= increase;
	BaryonField[DIINum][index] *= increase;
	BaryonField[HIINum][index] *= increase;
	BaryonField[HDINum][index] *= increase;
      }
      
      if (MetallicityField == TRUE)
	BaryonField[MetalNum][index] = rho_metal_jet;

      if (MBHColourNum > 0)
	BaryonField[MBHColourNum][index] = rho_colour_jet;
      
      CellsModified++;

    }  // END ic in n_cell_edge

    /* Remove ejected mass from star; now that we injected the NotEjectedMass, set it to zero */

    double old_mass = cstar->Mass;
    float old_vel1 = cstar->vel[1];
    cstar->Mass -= cstar->NotEjectedMass; 
    cstar->NotEjectedMass = 0.0;

    cstar->vel[0] *= old_mass / cstar->Mass; 
    cstar->vel[1] *= old_mass / cstar->Mass;
    cstar->vel[2] *= old_mass / cstar->Mass; 
//    fprintf(stdout, "grid::AFS:  Mass = %lf -> %lf, vel[1] = %f -> %f\n", 
//	    old_mass, cstar->Mass, old_vel1, cstar->vel[1]);  

    fprintf(stdout, "grid::AddFS: jets injected (EjectaM = %g Ms, JetsVelocity = %g, grid v = %g, rho_jet = %g) along n_L = (%g, %g, %g)\n", 
	    EjectaMass * DensityUnits * pow(LengthUnits,3.0) / Msun,
	    MBHJetsVelocity, BaryonField[Vel3Num][n_cell_edge-1], rho_jet, nx_L, ny_L, nz_L); 

  }  // END MBH_JETS

  /***********************************************************************
          STROEMGREN SPHERE FROM A POP III STAR (WHALEN ET AL. 2004)
  ************************************************************************/

  if (cstar->FeedbackFlag == STROEMGREN) {
    
    maxVelocity = 1e5*WhalenMaxVelocity / VelocityUnits;
    coef = maxVelocity / ( 0.8*0.8 + 2*0.8 );
    index = 0;

    for (k = 0; k < GridDimension[2]; k++) {

      delz = CellLeftEdge[2][k] + 0.5*CellWidth[2][k] - cstar->pos[2];
      sz = sign(delz);
      delz = fabs(delz);
      delz = min(delz, DomainWidth[2]-delz);

      for (j = 0; j < GridDimension[1]; j++) {

	dely = CellLeftEdge[1][j] + 0.5*CellWidth[1][j] - cstar->pos[1];
	sy = sign(dely);
	dely = fabs(dely);
	dely = min(dely, DomainWidth[1]-dely);

	for (i = 0; i < GridDimension[0]; i++, index++) {

	  delx = CellLeftEdge[0][i] + 0.5*CellWidth[0][i] - cstar->pos[0];
	  sx = sign(delx);
	  delx = fabs(delx);
	  delx = min(delx, DomainWidth[0]-delx);

	  radius2 = delx*delx + dely*dely + delz*delz;
	  if (radius2 <= 1.2*1.2*radius*radius) {

	    float r1 = sqrt(radius2) / radius;
	    float ramp = min(max(1.0 - (r1 - 0.8)/0.4, 0.01), 1.0);
	    
	    BaryonField[DensNum][index] = EjectaDensity * ramp;
	  
	    // Estimate of velocity profile
	    speed = coef * (r1*r1 + 2*r1) * ramp;
	    Radius = sqrt(radius2);
	    
	    if (Radius > TOLERANCE) {
	      BaryonField[Vel1Num][index] = sx*speed*delx/Radius + cstar->vel[0];
	      BaryonField[Vel2Num][index] = sy*speed*dely/Radius + cstar->vel[1];
	      BaryonField[Vel3Num][index] = sz*speed*delz/Radius + cstar->vel[2];
	    } else {
	      BaryonField[Vel1Num][index] = cstar->vel[0];
	      BaryonField[Vel2Num][index] = cstar->vel[1];
	      BaryonField[Vel3Num][index] = cstar->vel[2];
	    }

	    BaryonField[TENum][index] = EjectaThermalEnergy*ramp;
	    if (GENum >= 0 && DualEnergyFormalism)
	      BaryonField[GENum][index] = EjectaThermalEnergy*ramp;
	    
	    if (HydroMethod != Zeus_Hydro)
	      for (dim = 0; dim < GridRank; dim++)
		BaryonField[TENum][index] += 
		  0.5 * BaryonField[Vel1Num+dim][index] *
		  BaryonField[Vel1Num+dim][index];
	    
	    CellsModified++;

	  }  // END if inside radius

	}  // END i-direction
      }  // END j-direction
    }  // END k-direction
  }  // END Pop III sphere

  /***********************************************************************
          BARYON REMOVAL AFTER ACCRETION OR STAR PARTICLE CREATION
         -- Enforce a minimum temperature for cold gas accretion --
  ************************************************************************/

  float MinimumTemperature = 1.0, AdditionalEnergy, GasEnergy;

  if (cstar->FeedbackFlag == FORMATION) {

    index = 0;
    if (cstar->type == PopII)
      MinimumTemperature = (MultiSpecies > 1) ? 1e3 : 1e4;
    
    for (k = 0; k < GridDimension[2]; k++) {

      delz = CellLeftEdge[2][k] + 0.5*CellWidth[2][k] - cstar->pos[2];
      sz = sign(delz);
      delz = fabs(delz);
      delz = min(delz, DomainWidth[2]-delz);

      for (j = 0; j < GridDimension[1]; j++) {

	dely = CellLeftEdge[1][j] + 0.5*CellWidth[1][j] - cstar->pos[1];
	sy = sign(dely);
	dely = fabs(dely);
	dely = min(dely, DomainWidth[1]-dely);

	for (i = 0; i < GridDimension[0]; i++, index++) {

	  delx = CellLeftEdge[0][i] + 0.5*CellWidth[0][i] - cstar->pos[0];
	  sx = sign(delx);
	  delx = fabs(delx);
	  delx = min(delx, DomainWidth[0]-delx);

	  radius2 = delx*delx + dely*dely + delz*delz;
	  if (radius2 <= radius*radius) {

            //if (abs(cstar->type) == SimpleSource) {
	    //  factor = PopIIIStarMass/EjectaDensity;
	    //}
	    //else
	      factor = EjectaDensity / BaryonField[DensNum][index];
	    BaryonField[DensNum][index] *= factor;

	    if (MultiSpecies) {
	      BaryonField[DeNum][index] *= factor;
	      BaryonField[HINum][index] *= factor;
	      BaryonField[HIINum][index] *= factor;
	      BaryonField[HeINum][index] *= factor;
	      BaryonField[HeIINum][index] *= factor;
	      BaryonField[HeIIINum][index] *= factor;
	    }
	    if (MultiSpecies > 1) {
	      BaryonField[HMNum][index] *= factor;
	      BaryonField[H2INum][index] *= factor;
	      BaryonField[H2IINum][index] *= factor;
	    }
	    if (MultiSpecies > 2) {
	      BaryonField[DINum][index] *= factor;
	      BaryonField[DIINum][index] *= factor;
	      BaryonField[HIINum][index] *= factor;
	      BaryonField[HDINum][index] *= factor;
	    }

	    if (SNColourNum > 0)
	      BaryonField[SNColourNum][index] *= factor;
	    if (Metal2Num > 0)
	      BaryonField[Metal2Num][index] *= factor;

	    // For cold gas accretion, set a minimum temperature of
	    // 1e4 K since it has been accreted onto the star

	    if (DualEnergyFormalism) {
	      GasEnergy = BaryonField[GENum][index];
	    } else {
	      GasEnergy = BaryonField[TENum][index];
	      if (HydroMethod != Zeus_Hydro)
		for (dim = 0; dim < GridRank; dim++)
		  GasEnergy -= 0.5*BaryonField[Vel1Num+dim][index] * 
		    BaryonField[Vel1Num+dim][index];
	    }
	    AdditionalEnergy = 
	      MinimumTemperature / (TemperatureUnits * (Gamma-1.0) * 0.6) - 
	      GasEnergy;
	    AdditionalEnergy = max(AdditionalEnergy, 0.0);
	    BaryonField[TENum][index] += AdditionalEnergy;
	    if (DualEnergyFormalism)
	      BaryonField[GENum][index] += AdditionalEnergy;

	    CellsModified++;

	  }  // END if inside radius

	}  // END i-direction
      }  // END j-direction
    }  // END k-direction
  }  // END star birth

  /***********************************************************************
                             POPIII - Color field
  ************************************************************************/

  if (cstar->FeedbackFlag == COLOR_FIELD) {
    int CellsModified2 = 0;

    int ColorField = FindField(ForbiddenRefinement, FieldType, NumberOfBaryonFields); 
    if (ColorField < 0) ENZO_FAIL("Couldn't Find Color Field!");
    index = 0;

    for (k = 0; k < GridDimension[2]; k++) {

      delz = CellLeftEdge[2][k] + 0.5*CellWidth[2][k] - cstar->pos[2];
      sz = sign(delz);
      delz = fabs(delz);
      delz = min(delz, DomainWidth[2]-delz);

      for (j = 0; j < GridDimension[1]; j++) {

	dely = CellLeftEdge[1][j] + 0.5*CellWidth[1][j] - cstar->pos[1];
	sy = sign(dely);
	dely = fabs(dely);
	dely = min(dely, DomainWidth[1]-dely);

	for (i = 0; i < GridDimension[0]; i++, index++) {

	  delx = CellLeftEdge[0][i] + 0.5*CellWidth[0][i] - cstar->pos[0];
	  sx = sign(delx);
	  delx = fabs(delx);
	  delx = min(delx, DomainWidth[0]-delx);

	  radius2 = delx*delx + dely*dely + delz*delz;
	  if (radius2 <= radius*radius) {

	    BaryonField[ColorField][index] =
            BaryonField[DensNum][index];

	    CellsModified++;
	    CellsModified2++;

	  }  // END if inside radius

	}  // END i-direction
      }  // END j-direction
    }  // END k-direction
    fprintf(stderr, "CellsModified: %"ISYM"\n", CellsModified2);
  }

  /* Now it's done, unmark. */

  //cstar->FeedbackFlag = NO_FEEDBACK;

  return SUCCESS;

}
 
