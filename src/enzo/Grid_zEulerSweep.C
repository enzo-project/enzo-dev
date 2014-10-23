/***********************************************************************
/
/  GRID CLASS (WRAPPER FOR EULERIAN PPM SOLVER)
/
/  written by: John H. Wise
/  date:       May 2007
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
#include "euler_sweep.h"
//#include "fortran.def"

int grid::zEulerSweep(int j, int NumberOfSubgrids, fluxes *SubgridFluxes[], 
		      Elong_int GridGlobalStart[], float *CellWidthTemp[], 
		      int GravityOn, int NumberOfColours, int colnum[], float *pressure)
{

  int dim = 2, idim = 0, jdim = 1;
  int dim_p1 = dim+1;   // To match definition in calcdiss

  /* Find fields: density, total energy, velocity1-3. */
  
  int DensNum, GENum, Vel1Num, Vel2Num, Vel3Num, TENum;

  this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num, 
				   Vel3Num, TENum);

  int nxz, nyz, nzz, ixyz;

  nxz = GridEndIndex[0] - GridStartIndex[0] + 1;
  nyz = GridEndIndex[1] - GridStartIndex[1] + 1;
  nzz = GridEndIndex[2] - GridStartIndex[2] + 1;

  float MinimumPressure = tiny_number;
  
  // Copy from field to slice

  float *dslice, *eslice, *uslice, *vslice, *wslice, *grslice, *geslice, 
    *colslice, *pslice;

  int size = GridDimension[2] * GridDimension[0];
  dslice = new float[size];  
  eslice = new float[size];  
  uslice = new float[size];  
  vslice = new float[size];  
  wslice = new float[size];  
  pslice = new float[size];  
  if (GravityOn) {
    grslice = new float[size];  
  }
  if (DualEnergyFormalism) {
    geslice = new float[size];  
  }
  if (NumberOfColours > 0) {
    colslice = new float[NumberOfColours * size];  
  }

  int i, k, n, ncolour, index2, index3;

  for (i = 0; i < GridDimension[0]; i++) {
    index2 = i * GridDimension[2];
    for (k = 0; k < GridDimension[2]; k++) {
      index3 = (k*GridDimension[1] + j) * GridDimension[0] + i;
      dslice[index2+k] = BaryonField[DensNum][index3];
      eslice[index2+k] = BaryonField[TENum][index3];
      pslice[index2+k] = pressure[index3];
      vslice[index2+k] = BaryonField[Vel1Num][index3];
    } // ENDFOR i

    // Set velocities to zero if rank < 3 since hydro routines are
    // hard-coded for 3-d

    if (GridRank > 1)
      for (k = 0; k < GridDimension[2]; k++) {
	index3 = (k*GridDimension[1] + j) * GridDimension[0] + i;
	wslice[index2+k] = BaryonField[Vel2Num][index3];
      }
    else
      for (k = 0; k < GridDimension[2]; k++)
	wslice[index2+k] = 0;
  
    if (GridRank > 2)
      for (k = 0; k < GridDimension[2]; k++) {
	index3 = (k*GridDimension[1] + j) * GridDimension[0] + i;
	uslice[index2+k] = BaryonField[Vel3Num][index3];
      }
    else
      for (k = 0; k < GridDimension[2]; k++)
	uslice[index2+k] = 0;

    if (GravityOn)
      for (k = 0; k < GridDimension[2]; k++) {
	index3 = (k*GridDimension[1] + j) * GridDimension[0] + i;
	grslice[index2+k] = AccelerationField[dim][index3];
      }

    if (DualEnergyFormalism)
      for (k = 0; k < GridDimension[2]; k++) {
	index3 = (k*GridDimension[1] + j) * GridDimension[0] + i;
	geslice[index2+k] = BaryonField[GENum][index3];
      }

    for (n = 0; n < NumberOfColours; n++) {
      index2 = (n*GridDimension[0] + i) * GridDimension[2];
      for (k = 0; k < GridDimension[2]; k++) {
	index3 = (k*GridDimension[1] + j) * GridDimension[0] + i;
	colslice[index2+k] = BaryonField[colnum[n]][index3];
      }
    } // ENDFOR colours
  } // ENDFOR j

  /* Allocate memory for temporaries used in solver */

  float *dls, *drs, *flatten, *pbar, *pls, *prs, *ubar, *uls, *urs, *vls, 
    *vrs, *gels, *gers, *wls, *wrs, *diffcoef, *df, *ef, *uf, *vf, *wf, *gef,
    *ges, *colf, *colls, *colrs;

  dls = new float[size];	
  drs = new float[size];	
  flatten = new float[size];	
  pbar = new float[size];	
  pls = new float[size];	
  prs = new float[size];	
  ubar = new float[size];	
  uls = new float[size];	
  urs = new float[size];	
  vls = new float[size];	
  vrs = new float[size];	
  gels = new float[size];	
  gers = new float[size];	
  wls = new float[size];	
  wrs = new float[size];	
  diffcoef = new float[size];	
  df = new float[size];		
  ef = new float[size];		
  uf = new float[size];		
  vf = new float[size];		
  wf = new float[size];		
  gef = new float[size];	
  ges = new float[size];
  colf = new float[NumberOfColours*size];  
  colls = new float[NumberOfColours*size];  
  colrs = new float[NumberOfColours*size];  

  /* Convert start and end indexes into 1-based for FORTRAN */

  int is, ie, js, je, is_m3, ie_p3, ie_p1, k_p1;

  is = GridStartIndex[2] + 1;
  ie = GridEndIndex[2] + 1;
  js = 1;
  je = GridDimension[0];
  is_m3 = is - 3;
  ie_p1 = ie + 1;
  ie_p3 = ie + 3;
  k_p1 = j + 1;

  /* Compute the pressure on a slice */
  /*
  if (DualEnergyFormalism)
    FORTRAN_NAME(pgas2d_dual)(dslice, eslice, geslice, pslice, uslice, vslice, 
			      wslice, &DualEnergyFormalismEta1, 
			      &DualEnergyFormalismEta2, &GridDimension[2], 
			      &GridDimension[0], &is_m3, &ie_p3, &js, &je, 
			      &Gamma, &MinimumPressure);
  else
    FORTRAN_NAME(pgas2d)(dslice, eslice, pslice, uslice, vslice, 
			 wslice, &GridDimension[2], &GridDimension[0], 
			 &is_m3, &ie_p3, &js, &je, &Gamma, &MinimumPressure);
  */
  /* If requested, compute diffusion and slope flattening coefficients */

  if (PPMDiffusionParameter != 0 || PPMFlatteningParameter != 0)
    FORTRAN_NAME(calcdiss)(dslice, eslice, uslice, BaryonField[Vel1Num],
			   BaryonField[Vel2Num], pslice, CellWidthTemp[2],
			   CellWidthTemp[0], CellWidthTemp[1], 
			   &GridDimension[2], &GridDimension[0], 
			   &GridDimension[1], &is, &ie, &js, &je, &k_p1,
			   &nyz, &dim_p1, &GridDimension[0],
			   &GridDimension[1], &GridDimension[2],
			   &dtFixed, &Gamma, &PPMDiffusionParameter,
			   &PPMFlatteningParameter, diffcoef, flatten);

  /* Compute Eulerian left and right states at zone edges via interpolation */

  if (ReconstructionMethod == PPM)
    FORTRAN_NAME(inteuler)(dslice, pslice, &GravityOn, grslice, geslice, uslice,
			   vslice, wslice, CellWidthTemp[2], flatten,
			   &GridDimension[2], &GridDimension[0],
			   &is, &ie, &js, &je, &DualEnergyFormalism, 
			   &DualEnergyFormalismEta1, &DualEnergyFormalismEta2,
			   &PPMSteepeningParameter, &PPMFlatteningParameter,
			   &ConservativeReconstruction, &PositiveReconstruction,
			   &dtFixed, &Gamma, &PressureFree, 
			   dls, drs, pls, prs, gels, gers, uls, urs, vls, vrs,
			   wls, wrs, &NumberOfColours, colslice, colls, colrs);

  /* Compute (Lagrangian part of the) Riemann problem at each zone boundary */

  switch (RiemannSolver) {
  case TwoShock:
    FORTRAN_NAME(twoshock)(dls, drs, pls, prs, uls, urs,
			   &GridDimension[2], &GridDimension[0],
			   &is, &ie_p1, &js, &je,
			   &dtFixed, &Gamma, &MinimumPressure, &PressureFree,
			   pbar, ubar, &GravityOn, grslice,
			   &DualEnergyFormalism, &DualEnergyFormalismEta1);
    
    FORTRAN_NAME(flux_twoshock)(dslice, eslice, geslice, uslice, vslice, wslice,
				CellWidthTemp[2], diffcoef, 
				&GridDimension[2], &GridDimension[0],
				&is, &ie, &js, &je, &dtFixed, &Gamma,
				&PPMDiffusionParameter, &DualEnergyFormalism,
				&DualEnergyFormalismEta1,
				&RiemannSolverFallback,
				dls, drs, pls, prs, gels, gers, uls, urs,
				vls, vrs, wls, wrs, pbar, ubar,
				df, ef, uf, vf, wf, gef, ges,
				&NumberOfColours, colslice, colls, colrs, colf);
    break;

  case HLL:
    FORTRAN_NAME(flux_hll)(dslice, eslice, geslice, uslice, vslice, wslice,
			   CellWidthTemp[2], diffcoef, 
			   &GridDimension[2], &GridDimension[0],
			   &is, &ie, &js, &je, &dtFixed, &Gamma,
			   &PPMDiffusionParameter, &DualEnergyFormalism,
			   &DualEnergyFormalismEta1,
			   &RiemannSolverFallback,
			   dls, drs, pls, prs, uls, urs,
			   vls, vrs, wls, wrs, gels, gers,
			   df, uf, vf, wf, ef, gef, ges,
			   &NumberOfColours, colslice, colls, colrs, colf);
    break;

  case HLLC:
    FORTRAN_NAME(flux_hllc)(dslice, eslice, geslice, uslice, vslice, wslice,
			    CellWidthTemp[2], diffcoef, 
			    &GridDimension[2], &GridDimension[0],
			    &is, &ie, &js, &je, &dtFixed, &Gamma,
			    &PPMDiffusionParameter, &DualEnergyFormalism,
			    &DualEnergyFormalismEta1,
			    &RiemannSolverFallback,
			    dls, drs, pls, prs, uls, urs,
			    vls, vrs, wls, wrs, gels, gers,
			    df, uf, vf, wf, ef, gef, ges,
			    &NumberOfColours, colslice, colls, colrs, colf);
    break;

  } // ENDCASE

  /* Compute Eulerian fluxes and update zone-centered quantities */

  FORTRAN_NAME(euler)(dslice, eslice, grslice, geslice, uslice, vslice, wslice,
		      CellWidthTemp[2], diffcoef, 
		      &GridDimension[2], &GridDimension[0], 
		      &is, &ie, &js, &je, &dtFixed, &Gamma, 
		      &PPMDiffusionParameter, &GravityOn, &DualEnergyFormalism, 
		      &DualEnergyFormalismEta1, &DualEnergyFormalismEta2,
		      df, ef, uf, vf, wf, gef, ges,
		      &NumberOfColours, colslice, colf);

  /* If necessary, recompute the pressure to correctly set ge and e */

  if (DualEnergyFormalism)
    FORTRAN_NAME(pgas2d_dual)(dslice, eslice, geslice, pslice, uslice, vslice, 
			      wslice, &DualEnergyFormalismEta1, 
			      &DualEnergyFormalismEta2, &GridDimension[2], 
			      &GridDimension[0], &is_m3, &ie_p3, &js, &je, 
			      &Gamma, &MinimumPressure);

  /* Check this slice against the list of subgrids (all subgrid
     quantities are zero-based) */

  int jstart, jend, offset, nfi, lface, rface, lindex, rindex, 
    fistart, fiend, fjstart, fjend, clindex, crindex;
  
  for (n = 0; n < NumberOfSubgrids; n++) {

    fistart = SubgridFluxes[n]->RightFluxStartGlobalIndex[dim][idim] - 
      GridGlobalStart[idim];
    fiend = SubgridFluxes[n]->RightFluxEndGlobalIndex[dim][idim] -
      GridGlobalStart[idim];
    fjstart = SubgridFluxes[n]->RightFluxStartGlobalIndex[dim][jdim] - 
      GridGlobalStart[jdim];
    fjend = SubgridFluxes[n]->RightFluxEndGlobalIndex[dim][jdim] -
      GridGlobalStart[jdim];

    if (j >= fjstart && j <= fjend) {

      nfi = fiend - fistart + 1;
      for (i = fistart; i <= fiend; i++) {

	offset = (i-fistart) + (j-fjstart)*nfi;

	lface = SubgridFluxes[n]->LeftFluxStartGlobalIndex[dim][dim] -
	  GridGlobalStart[dim];
	lindex = i * GridDimension[dim] + lface;

	rface = SubgridFluxes[n]->RightFluxStartGlobalIndex[dim][dim] -
	  GridGlobalStart[dim] + 1;
	rindex = i * GridDimension[dim] + rface;	

	SubgridFluxes[n]->LeftFluxes [DensNum][dim][offset] = df[lindex];
	SubgridFluxes[n]->RightFluxes[DensNum][dim][offset] = df[rindex];
	SubgridFluxes[n]->LeftFluxes [TENum][dim][offset]   = ef[lindex];
	SubgridFluxes[n]->RightFluxes[TENum][dim][offset]   = ef[rindex];

	if (nxz > 1) {
	  SubgridFluxes[n]->LeftFluxes [Vel1Num][dim][offset] = vf[lindex];
	  SubgridFluxes[n]->RightFluxes[Vel1Num][dim][offset] = vf[rindex];
	} // ENDIF x-data

	if (nyz > 1) {
	  SubgridFluxes[n]->LeftFluxes [Vel2Num][dim][offset] = wf[lindex];
	  SubgridFluxes[n]->RightFluxes[Vel2Num][dim][offset] = wf[rindex];
	} // ENDIF y-data

	SubgridFluxes[n]->LeftFluxes [Vel3Num][dim][offset] = uf[lindex];
	SubgridFluxes[n]->RightFluxes[Vel3Num][dim][offset] = uf[rindex];

	if (DualEnergyFormalism) {
	  SubgridFluxes[n]->LeftFluxes [GENum][dim][offset] = gef[lindex];
	  SubgridFluxes[n]->RightFluxes[GENum][dim][offset] = gef[rindex];
	} // ENDIF DualEnergyFormalism

	for (ncolour = 0; ncolour < NumberOfColours; ncolour++) {
	  clindex = (i + ncolour * GridDimension[0]) * GridDimension[dim] +
	    lface;
	  crindex = (i + ncolour * GridDimension[0]) * GridDimension[dim] +
	    rface;

	  SubgridFluxes[n]->LeftFluxes [colnum[ncolour]][dim][offset] = 
	    colf[clindex];
	  SubgridFluxes[n]->RightFluxes[colnum[ncolour]][dim][offset] = 
	    colf[crindex];
	} // ENDFOR ncolour

      } // ENDFOR J

    } // ENDIF k inside

  } // ENDFOR n

  /* Copy from slice to field */

  for (i = 0; i < GridDimension[0]; i++) {
    index2 = i * GridDimension[2];
    for (k = 0; k < GridDimension[2]; k++) {
      index3 = (k*GridDimension[1] + j)*GridDimension[0] + i;
      BaryonField[DensNum][index3] = dslice[index2+k];
      BaryonField[TENum][index3] = eslice[index2+k];
      BaryonField[Vel1Num][index3] = vslice[index2+k];
    } // ENDFOR i

    if (GridRank > 1)
      for (k = 0; k < GridDimension[2]; k++) {
	index3 = (k*GridDimension[1] + j)*GridDimension[0] + i;
	BaryonField[Vel2Num][index3] = wslice[index2+k];
      }

    if (GridRank > 2)
      for (k = 0; k < GridDimension[2]; k++) {
	index3 = (k*GridDimension[1] + j)*GridDimension[0] + i;
	BaryonField[Vel3Num][index3] = uslice[index2+k];
      }

    if (DualEnergyFormalism)
      for (k = 0; k < GridDimension[2]; k++) {
	index3 = (k*GridDimension[1] + j)*GridDimension[0] + i;
	BaryonField[GENum][index3] = geslice[index2+k];
      }

    for (n = 0; n < NumberOfColours; n++) {
      index2 = (n*GridDimension[0] + i) * GridDimension[2];
      for (k = 0; k < GridDimension[2]; k++) {
	index3 = (k*GridDimension[1] + j) * GridDimension[0] + i;
	BaryonField[colnum[n]][index3] = colslice[index2+k];
      }
    } // ENDFOR colours

  } // ENDFOR j

  /* Delete all temporary slices */

  delete [] dslice;
  delete [] eslice;
  delete [] uslice;
  delete [] vslice;
  delete [] wslice;
  delete [] pslice;
  if (GravityOn)
    delete [] grslice;
  if (DualEnergyFormalism)
    delete [] geslice;
  if (NumberOfColours > 0)
    delete [] colslice;

  delete [] dls;
  delete [] drs;
  delete [] flatten;
  delete [] pbar;
  delete [] pls;
  delete [] prs;
  delete [] ubar;
  delete [] uls;
  delete [] urs;
  delete [] vls;
  delete [] vrs;
  delete [] gels;
  delete [] gers;
  delete [] wls;
  delete [] wrs;
  delete [] diffcoef;
  delete [] df;
  delete [] ef;
  delete [] uf;
  delete [] vf;
  delete [] wf;
  delete [] gef;
  delete [] ges;
  delete [] colf;
  delete [] colls;
  delete [] colrs;

  return SUCCESS;

}
