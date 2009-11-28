
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

  int dim, i, j, k, index;
  int sx, sy, sz;
  FLOAT delx, dely, delz, radius2, Radius, DomainWidth[MAX_DIMENSION];
  float coef, speed, maxVelocity;
  float OldDensity, old_mass;
  float r1, norm, ramp, factor, newGE, increase, fh;

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

  int SNColourNum, MetalNum, MBHColourNum, Galaxy1ColourNum, Galaxy2ColourNum; 
  int MetallicityField = FALSE;

  if (this->IdentifyColourFields(SNColourNum, MetalNum, MBHColourNum, 
				 Galaxy1ColourNum, Galaxy2ColourNum) == FAIL) {
    fprintf(stderr, "Error in grid->IdentifyColourFields.\n");
    ENZO_FAIL("");
  }

  MetalNum = max(MetalNum, SNColourNum);
  MetallicityField = (MetalNum > 0) ? TRUE : FALSE;


  /***********************************************************************
                          SUPERNOVAE / MBH_THERMAL
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

  // Correct for exaggerated influence radius for pair-instability supernovae
  if (cstar->FeedbackFlag == SUPERNOVA)
    radius /= 8.0;

  // Correct if the volume with 27 cells is larger than the energy bubble volume
  float BoxVolume = 27 * CellWidth[0][0] * CellWidth[0][0] * CellWidth[0][0];
  float BubbleVolume = (4.0 * M_PI / 3.0) * radius * radius * radius;
  //printf("BoxVolume = %lg, BubbleVolume = %lg\n", BoxVolume, BubbleVolume);
  if (BoxVolume > BubbleVolume) {
    //printf("Reducing ejecta density by %g\n", BubbleVolume / BoxVolume);
    EjectaDensity *= BubbleVolume / BoxVolume;
    EjectaMetalDensity *= BubbleVolume / BoxVolume;
    EjectaThermalEnergy *= BubbleVolume / BoxVolume;
  }
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

  if (cstar->FeedbackFlag == SUPERNOVA || 
      cstar->FeedbackFlag == CONT_SUPERNOVA || 
      cstar->FeedbackFlag == MBH_THERMAL) {

    /* Remove mass from the star that will now be added to grids. 
       Also, because EjectaDensity will be injected with zero net momentum, 
       increase the particle's velocity accordingly. - Ji-hoon Kim, Sep.2009 */

    old_mass = cstar->Mass;
    cstar->Mass -= EjectaDensity * DensityUnits * BubbleVolume * pow(LengthUnits,3.0) / Msun;  
    cstar->vel[0] *= old_mass / cstar->Mass; 
    cstar->vel[1] *= old_mass / cstar->Mass;
    cstar->vel[2] *= old_mass / cstar->Mass; 

    //printf("SN: pos = %"FSYM" %"FSYM" %"FSYM"\n", 
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
	    //	     ramp = min(max(1.0 - (r1 - 0.8)/0.4, 0.01), 1.0);

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
	      /*
	      printf("AddSN: rho = %"GSYM"=>%"GSYM", GE = %"GSYM"=>%"GSYM", drho = %"GSYM", dE = %"GSYM"\n",
		     OldDensity, BaryonField[DensNum][index], 
		     BaryonField[GENum][index], newGE, EjectaDensity,
		     EjectaThermalEnergy);
	      */
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

#define NOT_SEDOV_TEST
#ifdef SEDOV_TEST
	      newGE = (OldDensity * BaryonField[TENum][index] +
		       ramp * factor * EjectaThermalEnergy) /
		BaryonField[DensNum][index];
#endif
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

	    /* MBHColour injected - Ji-hoon Kim, Oct.2009 */
	    if (cstar->FeedbackFlag == MBH_THERMAL && MBHColourNum > 0)
	      BaryonField[MBHColourNum][index] += factor*EjectaDensity;

	    CellsModified++;

	  } // END if inside radius
	}  // END i-direction
      }  // END j-direction
    }  // END k-direction

  }  // END Supernova & MBH_THERMAL

  /***********************************************************************
                                 MBH_JETS
  ************************************************************************/

  // Inject bipolar jets along the direction of the angular momentum 
  // vector L of the MBH particle (angular momentum accreted thus far).
  // Follow a similar approach as was adopted in star_maker8.C

#define MAX_SUPERCELL_NUMBER 1000
  int ind_cell[MAX_SUPERCELL_NUMBER];
  int n_cell, ibuff = DEFAULT_GHOST_ZONES;
  int ii, jj, kk, r_s, ic, sign;
  float nx_cell[MAX_SUPERCELL_NUMBER], ny_cell[MAX_SUPERCELL_NUMBER], nz_cell[MAX_SUPERCELL_NUMBER];
  float m_cell, m_metal_cell, v_jet, rho_jet, rho_metal_jet, costheta = cos(3.1415926/3.9);
  float L_x, L_y, L_z, L_s, nx_L, ny_L, nz_L;

  if (cstar->FeedbackFlag == MBH_JETS) {

    /* Remove mass from the star that will now be added to grids. 
       Also, because EjectaDensity will be injected with zero net momentum, 
       increase the particle's velocity accordingly.*/

    old_mass = cstar->Mass;
    cstar->Mass -= EjectaDensity * DensityUnits * BubbleVolume * pow(LengthUnits,3.0) / Msun;  
    cstar->vel[0] *= old_mass / cstar->Mass; 
    cstar->vel[1] *= old_mass / cstar->Mass;
    cstar->vel[2] *= old_mass / cstar->Mass; 

    //printf("SN: pos = %"FSYM" %"FSYM" %"FSYM"\n", 
    //	   cstar->pos[0], cstar->pos[1], cstar->pos[2]);
    maxGE = MAX_TEMPERATURE / (TemperatureUnits * (Gamma-1.0) * 0.6);


    i = (cstar->pos[0] - CellLeftEdge[0][0])/CellWidth[0][0];
    j = (cstar->pos[1] - CellLeftEdge[1][0])/CellWidth[0][0];
    k = (cstar->pos[2] - CellLeftEdge[2][0])/CellWidth[0][0];

    /* Decide whether the current grid contains the whole supercell */
    
    if (i < ibuff+2 || i > GridDimension[0]-ibuff-3 || 
	j < ibuff+2 || j > GridDimension[1]-ibuff-3 ||
	k < ibuff+2 || k > GridDimension[2]-ibuff-3) 
      return SUCCESS;
    
    /* Find the the direction of angular momentum L accreted thus far */

    L_x = 0.0;
    L_y = 0.0;
    L_z = 1.0;  //##### fix this!
    L_s = sqrt(pow(L_x,2) + pow(L_y,2) + pow(L_z,2));
    nx_L = L_x/L_s;
    ny_L = L_y/L_s;
    ny_L = L_z/L_s;

    /* Loop over the supercell around the MBH particle (5 * 5 * 5 = 125 cells, 
       but only the edges), and record the cells eligible for jet injection */

    n_cell = 0;
    m_cell = 0.0;
    m_metal_cell = 0.0;

    for (kk = -2; kk <= 2; kk++) {
      for (jj = -2; jj <= 2; jj++) {
	for (ii = -2; ii <= 2; ii++) {
	  if (fabs(ii) != 2 && fabs(jj) != 2 && fabs(kk) != 2) continue;

	  r_s = sqrt(pow(ii,2) + pow(jj,2) + pow(kk,2));	    
	  
	  if (fabs((ii*nx_L + jj*ny_L + kk*nz_L)/r_s) > costheta) { 
	    ind_cell[n_cell] = i+ii+(j+jj+(k+kk)*GridDimension[1])*GridDimension[0];
	    nx_cell[n_cell]  = ii / r_s;
	    ny_cell[n_cell]  = jj / r_s;
	    nz_cell[n_cell]  = kk / r_s;
	    m_cell += BaryonField[DensNum][ind_cell[n_cell]] * pow(CellWidth[0][0], 3);
	    if (MetallicityField == TRUE) 
	      m_metal_cell += BaryonField[MetalNum][ind_cell[n_cell]] * pow(CellWidth[0][0], 3);
	    n_cell++;
	  }
	}
      }
    }

    /* Calculate the feedback velocity assuming EjectaThermalEnergy calculated 
       in Star_CalculateFeedbackParameter is now purely kinetic */

    v_jet = sqrt(0.5 * EjectaThermalEnergy);  //##### fix this, should be lower than this!
    
    /* Calculate the jet density.  Since EjectaDensity is calculated for isotropic 
       feedback (MBH_THERMAL), we need to reestimate the jet density here; we use 
       EjectaDensity only to compute EjectaMass = EjectaDensity * BubbleVolume */

    rho_jet = (m_cell + EjectaDensity * BubbleVolume) / 
      (n_cell * pow(CellWidth[0][0], 3));

    if (MetallicityField == TRUE) 
      rho_metal_jet = (m_metal_cell + EjectaMetalDensity * BubbleVolume) / 
	(n_cell * pow(CellWidth[0][0], 3));

    fprintf(stdout, "grid::AddFeedbackSphere: Jets injected v_jet = %g, rho_jet = %g\n", 
	    v_jet, rho_jet);

    if (v_jet * VelocityUnits > 1e9) {
      fprintf(stderr, "Jet velocity in MBH_JETS feedback too high! \n");
      ENZO_FAIL("");
    }

    /* Now do the feedback */

    for (ic = 0; ic < n_cell; ic++) {

      index = ind_cell[ic];

      BaryonField[DensNum][index] = rho_jet;
      
      /* Update velocities and TE; note that we now have kinetic (jet) energy added, so 
	 for DualEnergyFormalism = 0 you don't have to update any energy field */

      sign = sign(nx_cell[ic]*nx_L + ny_cell[ic]*ny_L + nz_cell[ic]*nz_L);

      if (GENum >= 0 && DualEnergyFormalism) 
	for (dim = 0; dim < GridRank; dim++)
	  BaryonField[TENum][index] -= 
	    0.5 * BaryonField[Vel1Num+dim][index] * 
	    BaryonField[Vel1Num+dim][index];

      BaryonField[Vel1Num][index] = sign*nx_L*v_jet;
      BaryonField[Vel2Num][index] = sign*ny_L*v_jet;
      BaryonField[Vel3Num][index] = sign*nz_L*v_jet;
      
      if (GENum >= 0 && DualEnergyFormalism) 
	for (dim = 0; dim < GridRank; dim++)
	  BaryonField[TENum][index] += 
	    0.5 * BaryonField[Vel1Num+dim][index] * 
	    BaryonField[Vel1Num+dim][index];

      /* Update species and colour fields */

      if (MetallicityField == TRUE) {
	metallicity = rho_jet / rho_metal_jet;
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
	BaryonField[MetalNum][index] = metallicity * BaryonField[DensNum][index]; // =rho_metal_jet
      
      if (MBHColourNum > 0)
	BaryonField[MBHColourNum][index] += (EjectaDensity * BubbleVolume) / 
	                                     (n_cell * pow(CellWidth[0][0], 3));

      CellsModified++;

    }  // END ic in n_cell

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
	    if (GENum >= 0)
	      BaryonField[GENum][index] = EjectaThermalEnergy*ramp;
	    
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
      MinimumTemperature = 1e4;
    
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

	    if (MetallicityField == TRUE)
	      BaryonField[MetalNum][index] *= factor;

	    // For cold gas accretion, set a minimum temperature of
	    // 1e4 K since it has been accreted onto the star

	    if (DualEnergyFormalism) {
	      GasEnergy = BaryonField[GENum][index];
	    } else {
	      GasEnergy = BaryonField[TENum][index];
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
