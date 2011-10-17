/***********************************************************************
/
/  GRID CLASS (Find Shocks)
/
/  written by: Sam Skillman
/  date:       May, 2008
/  modified1: 
/
/  PURPOSE:Finds all shock mach numbers 
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#include <math.h> 
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def"
#include "phys_constants.h"
#include "CosmologyParameters.h"

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

//Use Unsplit Temperature Jump Identification
int grid::FindShocks()
{
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  if (NumberOfBaryonFields == 0)
    return SUCCESS;
 
  this->DebugCheck((char *)"FindShocks");

  float invdx = 1./(CellWidth[0][0]);
  float inv2dx = 1./(2.*CellWidth[0][0]);
  float inv2dx2 = invdx*invdx;

  int is=GridStartIndex[0];
  int js=GridStartIndex[1];
  int ks=GridStartIndex[2];
  int ie=GridEndIndex[0];
  int je=GridEndIndex[1];
  int ke=GridEndIndex[2];

  int i, j, k, index,
    tempi, posti, prei;
  float preT, postT, tempjumpmag,
    gradtx, gradty, gradtz,
    maxdiv, temprat, tempmach;
  float num;

  int DensNum, TENum, GENum, 
    Vel1Num, Vel2Num, Vel3Num;

  int MachNum, PSTempNum, PSDenNum;

  /* Compute size (in floats) of the current grid. */
  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
  
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }
   
  // Get Shock species fields.
  if (IdentifyShockSpeciesFields(MachNum,PSTempNum,PSDenNum) == FAIL) {
      ENZO_FAIL("Error in IdentifyShockSpeciesFields.");
  }

  /* Get easy to handle pointers for each variable. */
 
  float *density     = BaryonField[DensNum];
  float *velocity1   = BaryonField[Vel1Num];
  float *velocity2   = BaryonField[Vel2Num];
  float *velocity3   = BaryonField[Vel3Num];
  float *mach        = BaryonField[MachNum];
  float *pstemp      = BaryonField[PSTempNum];
  float *psden       = BaryonField[PSDenNum];

  /* Create temperature, entropy fields */
  float *entropy = new float[size];
  float *tempgrad_dot_entropygrad = new float[size];
  double *flowdivergence = new double[size];
    /* If using cosmology, compute the expansion factor and get units. */
  
  float *temperature = new float[size]; 
  if (this->ComputeTemperatureField(temperature) == FAIL){
    ENZO_FAIL("Error in grid->ComputeTemperatureField.");
  }
  
  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
    VelocityUnits = 1, TimeUnits = 1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.");
  }
 
  // calculate cell temperature, set default values for mach
  for (i=0;i<size;i++){
    entropy[i] = tempgrad_dot_entropygrad[i] = mach[i] = 0.0;
    flowdivergence[i] = (double)(0.0);
    entropy[i] = temperature[i] / (pow(density[i],(Gamma - 1.0)));
    mach[i] = 0.0;
    if(StorePreShockFields){
      pstemp[i] = 0.0;
      psden[i] = 0.0;
    }
  }

  //Calculate temperature gradient dotted with entropy gradient
  //Calculate the flow divergence.
  //Do this for all cells except last ghost zone.
  int kstart = 1;
  int kend = GridDimension[2]-1;
  int jstart = 1;
  int jend = GridDimension[1]-1;
  if(GridRank < 3){
    kstart=0;
    kend=1;
  }
  if(GridRank < 2){
    jstart=0;
    jend=1;
  }

  for(k=kstart;k<kend;k++){
    for(j=jstart; j<jend;j++){
      for(i=1; i<GridDimension[0]-1;i++){
	
	index = i + GridDimension[0]*(j + GridDimension[1]*k);
	
	tempgrad_dot_entropygrad[index] = inv2dx2*
	  ((max(ShockTemperatureFloor,temperature[index+1])-max(ShockTemperatureFloor,temperature[index-1]))*
	   (entropy[index+1]-entropy[index-1]));
	if (GridRank > 1)
	  tempgrad_dot_entropygrad[index] += inv2dx2*
	    ((max(ShockTemperatureFloor,temperature[index+GridDimension[0]])-
	      max(ShockTemperatureFloor,temperature[index-GridDimension[0]]))*
	    (entropy[index+GridDimension[0]]-
	     entropy[index-GridDimension[0]]));
	if (GridRank > 2)
	  tempgrad_dot_entropygrad[index] += inv2dx2*
	    ((max(ShockTemperatureFloor,temperature[index+GridDimension[0]*GridDimension[1]])-
	      max(ShockTemperatureFloor,temperature[index-GridDimension[0]*GridDimension[1]]))*
	     (entropy[index+GridDimension[0]*GridDimension[1]]-
	      entropy[index-GridDimension[0]*GridDimension[1]]));
	
	flowdivergence[index] = (double)(inv2dx)*
	  ((double)(velocity1[index+1]) - (double)(velocity1[index-1]));
	if (GridRank > 1)
	  flowdivergence[index] += (double)(inv2dx)*
	    ((double)(velocity2[index+GridDimension[0]])- 
	     (double)(velocity2[index-GridDimension[0]]));   
	if (GridRank > 2)
	  flowdivergence[index] += (double)(inv2dx)*
	    ((double)(velocity3[index+GridDimension[0]*GridDimension[1]])- 
	     (double)(velocity3[index-GridDimension[0]*GridDimension[1]]));



      }
    }
  }
  
//Pseudo-Code
/* -----------------Pseudo Code------------------- /
   Check a cell if it is is a shock
   If not, continue
   If yes, then determine temperature gradient
   Search for pre/postshock cell
   If we come to a higher divergence, break
   If we come to a local minimum in divergence, make it the pre/post-shock
      cell, and stop looking
   If we come to a non-negative divergence, make it the pre/post-shock
      cell, and stop looking
   Make sure that the temperature/density keep increasing or decreasing
   After the ends are found, compute the mach number.

   In order to get statistics in terms of pre-shock quantities later, put 
     the shock in the pre-shock cell(i.e. mach number, cr...)
   Done!
/ ----------------------------------------------- */
  if(GridRank < 3){
    ks=0;
    ke=0;
  }
  if(GridRank < 2){
    js=0;
    je=0;
  }
  for(k=ks; k<=ke;k++){
    for(j=js; j<=je;j++){
      for(i=is; i<=ie;i++){
	
 	index = i + GridDimension[0]*(j + GridDimension[1]*k);

	if(tempgrad_dot_entropygrad[index] <= 0.0 || 
	   flowdivergence[index] >= 0.0)
	  continue;

	gradtx = gradty = gradtz = 0.0;

	preT = temperature[index];
	postT = temperature[index];	

	tempjumpmag = 
	  (max(ShockTemperatureFloor,temperature[index+1])-max(ShockTemperatureFloor,temperature[index-1]))*
	  (max(ShockTemperatureFloor,temperature[index+1])-max(ShockTemperatureFloor,temperature[index-1]));
	if (GridRank > 1)
	  tempjumpmag += 
	    (max(ShockTemperatureFloor,temperature[index+GridDimension[0]])-
	     max(ShockTemperatureFloor,temperature[index-GridDimension[0]]))*
	    (max(ShockTemperatureFloor,temperature[index+GridDimension[0]])-
	     max(ShockTemperatureFloor,temperature[index-GridDimension[0]]));
	if (GridRank > 2)
	  tempjumpmag += 
	    (max(ShockTemperatureFloor,temperature[index+GridDimension[0]*GridDimension[1]])-
	     max(ShockTemperatureFloor,temperature[index-GridDimension[0]*GridDimension[1]]))*
	    (max(ShockTemperatureFloor,temperature[index+GridDimension[0]*GridDimension[1]])-
	     max(ShockTemperatureFloor,temperature[index-GridDimension[0]*GridDimension[1]]));

	tempjumpmag = sqrt(tempjumpmag);

	gradtx = 
	  (max(ShockTemperatureFloor,temperature[index+1])-max(ShockTemperatureFloor,temperature[index-1]))/
	  tempjumpmag;
	if (GridRank > 1)
	  gradty = 
	    (max(ShockTemperatureFloor,temperature[index+GridDimension[0]])-
	     max(ShockTemperatureFloor,temperature[index-GridDimension[0]]))/tempjumpmag;
	if (GridRank > 2)
	  gradtz = 
	    (max(ShockTemperatureFloor,temperature[index+GridDimension[0]*GridDimension[1]])-
	     max(ShockTemperatureFloor,temperature[index-GridDimension[0]*GridDimension[1]]))/
	    tempjumpmag;

	
	num=0.0;
	maxdiv = flowdivergence[index];
	tempi = index;
	while(true){
	  //Find next post-cell along temperature gradient
	  //Make sure you are still in the grid
	  if( ((i+(int)(num*gradtx)) > (GridDimension[0]-1)) ||
	      ((i+(int)(num*gradtx)) < 0) )
	    break;
	  posti = index + (int)(num*gradtx);
	  
	  if (GridRank > 1){
	    if( ((j+(int)(num*gradty)) > (GridDimension[1]-1)) ||
		((j+(int)(num*gradty)) < 0) )
	      break;
	    posti += ((int)(num*gradty))*GridDimension[0];
	  }
	  if (GridRank > 2){
	    if( ((k+(int)(num*gradtz)) > (GridDimension[2]-1))  ||
		((k+(int)(num*gradtz)) < 0) )
	      break;
	    posti += ((int)(num*gradtz))*GridDimension[0]*GridDimension[1];
	  }

	  //If we haven't gone anywhere, increment num.
	  if(posti == tempi){
	    num++;
	    continue;
	  }
	  //Make sure temperature keeps increasing
	  if(temperature[posti] < postT){
	    posti = tempi;
	    break;
	  }
	  //Check for a shock in the current cell.  If not, set postT
	  //and break out.
	  if(tempgrad_dot_entropygrad[posti] <= 0.0 || 
	     flowdivergence[posti] >= 0.0){
	    postT = temperature[posti];	
	    break;
	  }
	  //Check for better center of the shock.  If so, get out.
	  if(flowdivergence[posti] < flowdivergence[index]){
	    num=-1;
	    break;
	  }
	  //Check for local maximum in divergence.  If so, set 
	  //postT and break out
	  if(flowdivergence[posti] < maxdiv){
	    //  postT = temperature[tempi];  //Debatable 
	    postT = temperature[posti];
	    break;
	  }
	  //Update temporary i, maximum divergence, and increment num.
	  tempi=posti;
	  maxdiv = flowdivergence[posti];
	  postT = temperature[posti];
	  num++;
	}
	//If a center was found, continue to next cell.
	if(num == -1)
	  continue;
	

	//Now find pre-shock cell
	num=0.0;
	maxdiv = flowdivergence[index];
	tempi = index;
	while(true){
	  //Find next pre-cell along max(ShockTemperatureFloor,temperature gradient
	  //Make sure you are still in the grid
	  if( ((i-(int)(num*gradtx)) > (GridDimension[0]-1)) ||
	      ((i-(int)(num*gradtx)) < 0) )
 	    break;
	  prei = index - (int)(num*gradtx);
	  
	  if (GridRank > 1){
	    if( ((j-(int)(num*gradty)) > (GridDimension[1]-1)) ||
		((j-(int)(num*gradty)) < 0) )
	      break;
	    prei -= ((int)(num*gradty))*GridDimension[0];
	  }
	  if (GridRank > 2){
	    if( ((k-(int)(num*gradtz)) > (GridDimension[2]-1))  ||
		((k-(int)(num*gradtz)) < 0) )
	      break;
	    prei -= ((int)(num*gradtz))*GridDimension[0]*GridDimension[1];
	  }

	  //If we haven't gone anywhere, increment num.
	  if(prei == tempi){
	    num++;
	    continue;
	  }
	  //Make sure temperature keeps decreasing
	  if(temperature[prei] > preT){
	    prei = tempi;
	    break;
	  }
	  //Check for a shock in the current cell.  If not, set preT
	  //and break out.
	  if(tempgrad_dot_entropygrad[prei] <= 0.0 || 
	     flowdivergence[prei] >= 0.0){
	    preT = temperature[prei];	
	    break;
	  }
	  //Check for better center of the shock.  If so, get out.
	  if(flowdivergence[prei] < flowdivergence[index]){
	    num=-1;
	    break;
	  }
	  //Check for local maximum in divergence.  If so, set 
	  //preT and break out
	  if(flowdivergence[prei] < maxdiv){
	    // preT = temperature[tempi];  //Debatable 
	    preT = temperature[prei];
	    break;
	  }
	  //Update temporary i, maximum divergence, and increment num.
	  tempi=prei;
	  maxdiv = flowdivergence[prei];
	  preT = temperature[prei];
	  num++;
	}
	//If a center was found, continue to next cell.
	if(num == -1)
	  continue;

	temprat = max(ShockTemperatureFloor,postT)/(max(ShockTemperatureFloor,preT));
	//temprat = max(postT,ShockTemperatureFloor)/(max(preT,ShockTemperatureFloor));
	
	if(temprat < 1.0)
	  continue;

	if(density[posti] < density[prei])
	  continue;

	tempmach = 
	  sqrt(( 8.0*temprat - 7.0e0 + 
		 sqrt( (7.0e0 - 8.0e0*temprat)*(7.0e0 - 8.0e0*temprat)
		       + 15.0e0) )/5.0e0); 
	if(tempmach <= 1.0)
	  continue;

	mach[index] = tempmach;

	if(StorePreShockFields){
	  pstemp[index] = max(temperature[prei],ShockTemperatureFloor);
	  psden[index] = density[prei];
	}
      }
    }
  }
  
  /* deallocate temporary space */
  
  delete [] temperature;
  delete [] flowdivergence;
  delete [] tempgrad_dot_entropygrad;
  delete [] entropy;
  
  return SUCCESS;
}

//Use Unsplit Velocity Jump Identification
int grid::FindVelShocks()
{
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  if (NumberOfBaryonFields == 0)
    return SUCCESS;
 
  this->DebugCheck((char *)"FindShocks");

  float inv2dx = 1./(2.*CellWidth[0][0]);

  int is=GridStartIndex[0];
  int js=GridStartIndex[1];
  int ks=GridStartIndex[2];
  int ie=GridEndIndex[0];
  int je=GridEndIndex[1];
  int ke=GridEndIndex[2];

  int i, j, k, index,
    tempi, posti, prei;
  float Csound, centervelx,centervely,centervelz,
    preV, postV,veljumpmag,
    v1jump,v2jump,v3jump,
    gradvx, gradvy, gradvz,
    maxdiv, thisjump, oldjump, velmach;
  float num;

  int DensNum, TENum, GENum, 
    Vel1Num, Vel2Num, Vel3Num;

  int MachNum, PSTempNum, PSDenNum;

  /* Compute size (in floats) of the current grid. */
  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
  
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }
   
  // Get CR species fields.

  if (IdentifyShockSpeciesFields(MachNum,PSTempNum,PSDenNum) == FAIL) {
    ENZO_FAIL("Error in IdentifyCRSpeciesFields.");
  }

  /* Get easy to handle pointers for each variable. */
 
  float *density     = BaryonField[DensNum];
  float *velocity1   = BaryonField[Vel1Num];
  float *velocity2   = BaryonField[Vel2Num];
  float *velocity3   = BaryonField[Vel3Num];
  float *mach        = BaryonField[MachNum];
  float *pstemp      = BaryonField[PSTempNum];
  float *psden       = BaryonField[PSDenNum];

  /* Create temperature, entropy fields */
  float *entropy = new float[size];
  float *tempgrad_dot_entropygrad = new float[size];
  double *flowdivergence = new double[size];
    /* If using cosmology, compute the expansion factor and get units. */
  
  float *temperature = new float[size];
  
  if (this->ComputeTemperatureField(temperature) == FAIL){
    ENZO_FAIL("Error in grid->ComputeTemperatureField.");
  }
  
  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
    VelocityUnits = 1, TimeUnits = 1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.");
  }
 
  // calculate cell temperature, set default values for mach
  for (i=0;i<size;i++){
    entropy[i] = tempgrad_dot_entropygrad[i] = mach[i] = 0.0;
    flowdivergence[i] = (double)(0.0);
    entropy[i] = temperature[i] / (pow(density[i],(Gamma - 1.0)));
    mach[i] = 0.0;
    if(StorePreShockFields){
      pstemp[i] = 0.0;
      psden[i] = 0.0;
    }
  }

  //  fprintf(stdout, "ShockTemperatureFloor= %e TemperatureUnits= %e\n",ShockTemperatureFloor,
  //  TemperatureUnits);


  //Calculate temperature gradient dotted with entropy gradient
  //Calculate the flow divergence.
  //Do this for all cells except last ghost zone.
  int kstart = 1;
  int kend = GridDimension[2]-1;
  int jstart = 1;
  int jend = GridDimension[1]-1;
  if(GridRank < 3){
    kstart=0;
    kend=1;
  }
  if(GridRank < 2){
    jstart=0;
    jend=1;
  }

  for(k=kstart;k<kend;k++){
    for(j=jstart; j<jend;j++){
      for(i=1; i<GridDimension[0]-1;i++){
	
	index = i + GridDimension[0]*(j + GridDimension[1]*k);
	
	flowdivergence[index] = (double)(inv2dx)*
	  ((double)(velocity1[index+1]) - (double)(velocity1[index-1]));
	if (GridRank > 1)
	  flowdivergence[index] += (double)(inv2dx)*
	    ((double)(velocity2[index+GridDimension[0]])- 
	     (double)(velocity2[index-GridDimension[0]]));   
	if (GridRank > 2)
	  flowdivergence[index] += (double)(inv2dx)*
	    ((double)(velocity3[index+GridDimension[0]*GridDimension[1]])- 
	     (double)(velocity3[index-GridDimension[0]*GridDimension[1]]));
      }
    }
  }
  
//Pseudo-Code
/* ----------------------------Pseudo Code---------------------------- /
 * Check a cell if it is is a shock If not, continue If yes, then
 * determine temperature * gradient Search for pre/postshock cell If we
 * come to a higher * divergence, break If we come to a local minimum in
 * divergence, make * it the pre/post-shock cell, and stop looking If we
 * come to a * non-negative divergence, make it the pre/post-shock cell,
 * and stop * looking Make sure that the temperature/density keep
 * increasing or * decreasing After the ends are found, compute the mach
 * number.  Then * compute the CR energy injected, multiplied by the
 * timestep!  In order to get statistics in terms of pre-shock
 * quantities later, put the shock in the pre-shock cell(i.e. mach
 * number, cr...)  Done!  /
 * ------------------------------------------------------------------ */
  if(GridRank < 3){
    ks=0;
    ke=0;
  }
  if(GridRank < 2){
    js=0;
    je=0;
  }
  for(k=ks; k<=ke;k++){
    for(j=js; j<=je;j++){
      for(i=is; i<=ie;i++){
	
 	index = i + GridDimension[0]*(j + GridDimension[1]*k);

	if(flowdivergence[index] >= 0.0)
	  continue;

	gradvx = gradvy = gradvz = 0.0;

	preV = index;
	postV = index;	

	veljumpmag = 			
	  (velocity1[index+1]-velocity1[index-1])*
	  (velocity1[index+1]-velocity1[index-1]);	
	if (GridRank > 1)
	  veljumpmag += 
	    (velocity2[index+GridDimension[0]]-
	     velocity2[index-GridDimension[0]])*
	    (velocity2[index+GridDimension[0]]-
	     velocity2[index-GridDimension[0]]);
	if (GridRank > 2)
	  veljumpmag += 
	    (velocity3[index+GridDimension[0]*GridDimension[1]]-
	     velocity3[index-GridDimension[0]*GridDimension[1]])*
	    (velocity3[index+GridDimension[0]*GridDimension[1]]-
	     velocity3[index-GridDimension[0]*GridDimension[1]]);

	veljumpmag = sqrt(veljumpmag);
	oldjump = 0.0;
	thisjump = 0.0;

	if(veljumpmag == 0.0)
	  continue;

	gradvx = 
	  (velocity1[index+1]-velocity1[index-1])/veljumpmag;
	if (GridRank > 1)
	  gradvy = 
	    (velocity2[index+GridDimension[0]]-
	     velocity2[index-GridDimension[0]])/veljumpmag;
	if (GridRank > 2)
	  gradvz = 
	    (velocity3[index+GridDimension[0]*GridDimension[1]]-
	     velocity3[index-GridDimension[0]*GridDimension[1]])/veljumpmag;

	//Correct gradvx to point in the same direction as temperature gradient
	
	gradvx = ( (temperature[index+1] > temperature[index-1]) 
		   ? (sqrt(gradvx*gradvx)) : (-sqrt(gradvx*gradvx)) );
	if(GridRank > 1)
	  gradvy = ( (temperature[index+GridDimension[0]] > 
		      temperature[index-GridDimension[0]]) 
		     ? (sqrt(gradvy*gradvy)) : (-sqrt(gradvy*gradvy)) );
	if(GridRank > 2)
	  gradvz = ( (temperature[index+GridDimension[0]*GridDimension[1]] > 
		      temperature[index-GridDimension[0]*GridDimension[1]]) 
		     ? (sqrt(gradvz*gradvz)) : (-sqrt(gradvz*gradvz)) );


	num=0.0;
	centervelx = velocity1[index];
	if(GridRank > 1)
	  centervely = velocity2[index];
	if(GridRank > 2)
	  centervelz = velocity3[index];

	maxdiv = flowdivergence[index];
	tempi = index;
	while(true){
	  //Find next post-cell along velocity gradient
	  //Make sure you are still in the grid
	  if( ((i+(int)(num*gradvx)) > (GridDimension[0]-1)) ||
	      ((i+(int)(num*gradvx)) < 0) )
	    break;
	  posti = index + (int)(num*gradvx);
	  
	  if (GridRank > 1){
	    if( ((j+(int)(num*gradvy)) > (GridDimension[1]-1)) ||
		((j+(int)(num*gradvy)) < 0) )
	      break;
	    posti += ((int)(num*gradvy))*GridDimension[0];
	  }
	  if (GridRank > 2){
	    if( ((k+(int)(num*gradvz)) > (GridDimension[2]-1))  ||
		((k+(int)(num*gradvz)) < 0) )
	      break;
	    posti += ((int)(num*gradvz))*GridDimension[0]*GridDimension[1];
	  }

	  //If we haven't gone anywhere, increment num.
	  if(posti == tempi){
	    num++;
	    continue;
	  }
	  v1jump = gradvx*(velocity1[posti] - centervelx);
	  if(GridRank > 1)
	    v2jump = gradvy*(velocity2[posti] - centervely);
	  if(GridRank > 2)	    
	    v3jump = gradvz*(velocity3[posti] - centervelz);
	  
	  thisjump = v1jump*v1jump;
	  if(GridRank > 1)
	    thisjump += v2jump*v2jump;
	  if(GridRank > 2)
	    thisjump += v3jump*v3jump;
	  thisjump = sqrt(thisjump);

	  //Make sure velocity (dot) velgrad keeps increasing
	  if(thisjump < oldjump){
	    posti = tempi;
	    break;
	  }
	  //Make sure density/temperature keeps increasing
	  if(temperature[posti] < temperature[tempi] ||
	     density[posti] < density[tempi]){
	    posti = tempi;
	    break;
	  }
	  oldjump = thisjump;
	  //Check for a shock in the current cell.  If not, set postV
	  //and break out. Use last actual shocked cell as cell, not
	  //the next one.
	  if(flowdivergence[posti] >= 0.0){
	    posti = tempi;	
	    break;
	  }
	  //Check for better center of the shock.  If so, get out.
	  if(flowdivergence[posti] < flowdivergence[index]){
	    num=-1;
	    break;
	  }
	  //Check for local maximum in divergence.  If so, set 
	  //postV and break out
	  if(flowdivergence[posti] < maxdiv){
	    //  postV = temperature[tempi];  //Debatable 
	    postV = posti;
	    break;
	  }
	  //Update temporary i, maximum divergence, and increment num.
	  tempi=posti;
	  maxdiv = flowdivergence[posti];
	  num++;
	}
	//If a center was found, continue to next cell.
	if(num == -1)
	  continue;

	thisjump = 0.0;

	//Now find pre-shock cell
	num=0.0;
	maxdiv = flowdivergence[index];
	tempi = index;
	while(true){
	  //Find next pre-cell along max(ShockTemperatureFloor,temperature gradient
	  //Make sure you are still in the grid
	  if( ((i-(int)(num*gradvx)) > (GridDimension[0]-1)) ||
	      ((i-(int)(num*gradvx)) < 0) )
 	    break;
	  prei = index - (int)(num*gradvx);
	  
	  if (GridRank > 1){
	    if( ((j-(int)(num*gradvy)) > (GridDimension[1]-1)) ||
		((j-(int)(num*gradvy)) < 0) )
	      break;
	    prei -= ((int)(num*gradvy))*GridDimension[0];
	  }
	  if (GridRank > 2){
	    if( ((k-(int)(num*gradvz)) > (GridDimension[2]-1))  ||
		((k-(int)(num*gradvz)) < 0) )
	      break;
	    prei -= ((int)(num*gradvz))*GridDimension[0]*GridDimension[1];
	  }

	  //If we haven't gone anywhere, increment num.
	  if(prei == tempi){
	    num++;
	    continue;
	  }
	  v1jump = gradvx*(velocity1[posti] - velocity1[prei]);
	  if(GridRank > 1)
	    v2jump = gradvy*(velocity2[posti] - velocity2[prei]);
	  if(GridRank > 2)	    
	    v3jump = gradvz*(velocity3[posti] - velocity3[prei]);
	  
	  thisjump = v1jump*v1jump;
	  if(GridRank > 1)
	    thisjump += v2jump*v2jump;
	  if(GridRank > 2)
	    thisjump += v3jump*v3jump;
	  thisjump = sqrt(thisjump);

	  //Make sure velocity jump keeps increasing
	  if(thisjump < oldjump){
	    prei = tempi;
	    break;
	  }
	  //Make sure density/temperature keeps decreasing
	  if(temperature[prei] > temperature[tempi] ||
	     density[prei] > density[tempi]){
	    posti = tempi;
	    break;
	  }
	  //Check for a shock in the current cell.  If not, set preV
	  //and break out.
	  if(flowdivergence[prei] >= 0.0){
	    preV = prei;
	    break;
	  }
	  //Check for better center of the shock.  If so, get out.
	  if(flowdivergence[prei] < flowdivergence[index]){
	    num=-1;
	    break;
	  }
	  //Check for local maximum in divergence.  If so, set 
	  //preV and break out
	  if(flowdivergence[prei] < maxdiv){
	    // preV = temperature[tempi];  //Debatable 
	    preV = prei;
	    break;
	  }
	  //Update temporary i, maximum divergence, and increment num.
	  tempi=prei;
	  maxdiv = flowdivergence[prei];
	  oldjump=thisjump;
	  num++;
	}
	//If a center was found, continue to next cell.
	if(num == -1)
	  continue;

	//temprat = max(ShockTemperatureFloor,postV)/(max(ShockTemperatureFloor,preV));
	//temprat = max(postV,ShockTemperatureFloor)/(max(preV,ShockTemperatureFloor));
	
 	if(max(temperature[posti],ShockTemperatureFloor) <= max(temperature[prei],ShockTemperatureFloor))
 	  continue;

 	if(density[posti] <= density[prei])
 	  continue;

	v1jump = gradvx*(velocity1[posti] - velocity1[prei]);
	if(GridRank > 1)
	  v2jump = gradvy*(velocity2[posti] - velocity2[prei]);
	if(GridRank > 2)	    
	  v3jump = gradvz*(velocity3[posti] - velocity3[prei]);
	
	thisjump = v1jump*v1jump;
	if(GridRank > 1)
	  thisjump += v2jump*v2jump;
	if(GridRank > 2)
	  thisjump += v3jump*v3jump;
	thisjump = sqrt(thisjump);

	//Speed of sound in code units of velocity

	//Figure out which one is the pre-shock cell:
	prei = ( (temperature[prei] < temperature[posti]) ? prei : posti);

	Csound = sqrt(Gamma*kboltz*temperature[prei]/(Mu*mh))/VelocityUnits;

	velmach = (1.0/3.0)*(2.0*thisjump/Csound + 
			     sqrt(9.0+4.0*thisjump*thisjump/Csound/Csound));

	if(velmach < 1.0)
	  continue;

	mach[index] = velmach; 

	if(StorePreShockFields){
	  pstemp[index] = max(temperature[prei],ShockTemperatureFloor);
	  psden[index] = density[prei];
	}
      }
    }
  }
  
  /* deallocate temporary space */
  
  delete [] temperature;
  delete [] flowdivergence;
  delete [] tempgrad_dot_entropygrad;
  delete [] entropy;
  
  return SUCCESS;
}

//Use Split Velocity Jump Identification
int grid::FindVelSplitShocks()
{
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  if (NumberOfBaryonFields == 0)
    return SUCCESS;
 
  this->DebugCheck((char *)"FindShocks");

  float inv2dx = 1./(2.*CellWidth[0][0]);

  int is=GridStartIndex[0];
  int js=GridStartIndex[1];
  int ks=GridStartIndex[2];
  int ie=GridEndIndex[0];
  int je=GridEndIndex[1];
  int ke=GridEndIndex[2];

  int i, j, k, index,
    prei;
  float v1jump, v2jump, v3jump, Csound;

  int DensNum, TENum, GENum, 
    Vel1Num, Vel2Num, Vel3Num;

  int MachNum, PSTempNum, PSDenNum;

  /* Compute size (in floats) of the current grid. */
  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
  
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }
   
  // Get CR species fields.

  if (IdentifyShockSpeciesFields(MachNum,PSTempNum,PSDenNum) == FAIL) {
    ENZO_FAIL("Error in IdentifyShockSpeciesFields.");
  }

  /* Get easy to handle pointers for each variable. */
 
  float *density     = BaryonField[DensNum];
  float *velocity1   = BaryonField[Vel1Num];
  float *velocity2   = BaryonField[Vel2Num];
  float *velocity3   = BaryonField[Vel3Num];
  float *mach        = BaryonField[MachNum];
  float *pstemp      = BaryonField[PSTempNum];
  float *psden       = BaryonField[PSDenNum];

  /* Create temperature, entropy fields */
  float *entropy = new float[size];
  float *tempgrad_dot_entropygrad = new float[size];
  double *flowdivergence = new double[size];
    /* If using cosmology, compute the expansion factor and get units. */
  
  float *temperature = new float[size];
  
  if (this->ComputeTemperatureField(temperature) == FAIL){
    ENZO_FAIL("Error in grid->ComputeTemperatureField.");
  }
  
  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
    VelocityUnits = 1, TimeUnits = 1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.");
  }
 
  // calculate cell temperature, set default values for mach
  for (i=0;i<size;i++){
    entropy[i] = tempgrad_dot_entropygrad[i] = mach[i] = 0.0;
    flowdivergence[i] = (double)(0.0);
    entropy[i] = temperature[i] / (pow(density[i],(Gamma - 1.0)));
    mach[i] = 0.0;
    if(StorePreShockFields){
      pstemp[i] = 0.0;
      psden[i] = 0.0;
    }
  }

  //  fprintf(stdout, "ShockTemperatureFloor= %e TemperatureUnits= %e\n",ShockTemperatureFloor, TemperatureUnits);


  //Calculate temperature gradient dotted with entropy gradient
  //Calculate the flow divergence.
  //Do this for all cells except last ghost zone.
  int kstart = 1;
  int kend = GridDimension[2]-1;
  int jstart = 1;
  int jend = GridDimension[1]-1;
  if(GridRank < 3){
    kstart=0;
    kend=1;
  }
  if(GridRank < 2){
    jstart=0;
    jend=1;
  }

  for(k=kstart;k<kend;k++){
    for(j=jstart; j<jend;j++){
      for(i=1; i<GridDimension[0]-1;i++){
	
	index = i + GridDimension[0]*(j + GridDimension[1]*k);
	
	flowdivergence[index] = (double)(inv2dx)*
	  ((double)(velocity1[index+1]) - (double)(velocity1[index-1]));
	if (GridRank > 1)
	  flowdivergence[index] += (double)(inv2dx)*
	    ((double)(velocity2[index+GridDimension[0]])- 
	     (double)(velocity2[index-GridDimension[0]]));   
	if (GridRank > 2)
	  flowdivergence[index] += (double)(inv2dx)*
	    ((double)(velocity3[index+GridDimension[0]*GridDimension[1]])- 
	     (double)(velocity3[index-GridDimension[0]*GridDimension[1]]));
      }
    }
  }
  
//Pseudo-Code
/* -----------------Pseudo Code------------------- /
   Check a cell if it is is a shock
   If not, continue
   If yes, then determine temperature gradient
   Search for pre/postshock cell
   If we come to a higher divergence, break
   If we come to a local minimum in divergence, make it the pre/post-shock
      cell, and stop looking
   If we come to a non-negative divergence, make it the pre/post-shock
      cell, and stop looking
   Make sure that the temperature/density keep increasing or decreasing
   After the ends are found, compute the mach number.
   Then compute the CR energy injected, multiplied by the timestep!

   In order to get statistics in terms of pre-shock quantities later, put 
     the shock in the pre-shock cell(i.e. mach number, cr...)
   Done!
/ ----------------------------------------------- */

//Different starting/ending to make sure we examine all pairs.  For 2 or 1D calculations, we only want to analyze k=j=0
  if(GridRank < 3){
    ks=1;
    ke=-1;
  }
  if(GridRank < 2){
    js=1;
    je=-1;
  }
  for(k=ks-1; k<=ke+1;k++){
    for(j=js-1; j<=je+1;j++){
      for(i=is-1; i<=ie+1;i++){
	
 	index = i + GridDimension[0]*(j + GridDimension[1]*k);

	if(flowdivergence[index] >= 0.0)
	  continue;

	//x-direction
	v1jump = fabs(velocity1[index+1]-velocity1[index-1]);
	if(velocity1[index+1] > velocity1[index-1])
	  prei = index-1;
	else
	  prei = index+1;
	
	Csound = sqrt(Gamma*kboltz*temperature[prei]/
		      (Mu*mh))/VelocityUnits;
	mach[prei] = (1.0/3.0)*(2.0*v1jump/Csound + 
				sqrt(9.0+4.0*v1jump*v1jump/Csound/Csound));
	
	//y-direction

	if (GridRank > 1){
	  v2jump = fabs(velocity2[index+GridDimension[0]]-velocity2[index-GridDimension[0]]);
	  if(velocity2[index+GridDimension[0]] > velocity2[index-GridDimension[0]])
	    prei = index-GridDimension[0];
	  else
	    prei = index+GridDimension[0];
	  Csound = sqrt(Gamma*kboltz*temperature[prei]/
			(Mu*mh))/VelocityUnits;
	  mach[prei] = sqrt(mach[prei]*mach[prei] + 
			    (1.0/3.0)*(2.0*v2jump/Csound + 
				       sqrt(9.0+4.0*v2jump*v2jump/Csound/Csound))*
			    (1.0/3.0)*(2.0*v2jump/Csound + 
				       sqrt(9.0+4.0*v2jump*v2jump/Csound/Csound)));
	}	
	
	//z-direction

	if (GridRank > 2){
	  v3jump = fabs(velocity3[index+GridDimension[0]*GridDimension[1]]-
			velocity3[index-GridDimension[0]*GridDimension[1]]);
	  if(velocity3[index+GridDimension[0]*GridDimension[1]] > 
	     velocity3[index-GridDimension[0]*GridDimension[1]])
	    prei = index-GridDimension[0]*GridDimension[1];
	  else
	    prei = index+GridDimension[0]*GridDimension[1];
	  Csound = sqrt(Gamma*kboltz*temperature[prei]/
			(Mu*mh))/VelocityUnits;
	  mach[prei] = sqrt(mach[prei]*mach[prei] + 
			    (1.0/3.0)*(2.0*v3jump/Csound + 
				       sqrt(9.0+4.0*v3jump*v3jump/Csound/Csound))*
			    (1.0/3.0)*(2.0*v3jump/Csound + 
				       sqrt(9.0+4.0*v3jump*v3jump/Csound/Csound)));
	}	

	if(StorePreShockFields){
	  pstemp[index] = max(temperature[prei],ShockTemperatureFloor);
	  psden[index] = density[prei];
	}	
	
      }
    }
  }
  
  /* deallocate temporary space */
  
  delete [] temperature;
  delete [] flowdivergence;
  delete [] tempgrad_dot_entropygrad;
  delete [] entropy;
  
  return SUCCESS;
}

//Use Split Temperature Jump Identification
int grid::FindTempSplitShocks()
{
  /* Return if this doesn't concern us. */
 
  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;
 
  if (NumberOfBaryonFields == 0)
    return SUCCESS;
 
  this->DebugCheck((char *)"FindShocks");

  float invdx = 1./(CellWidth[0][0]);
  float inv2dx = 1./(2.*CellWidth[0][0]);
  float inv2dx2 = invdx*invdx;

  int is=GridStartIndex[0];
  int js=GridStartIndex[1];
  int ks=GridStartIndex[2];
  int ie=GridEndIndex[0];
  int je=GridEndIndex[1];
  int ke=GridEndIndex[2];

  int i, j, k, index,posti, prei;
  float maxdiv;

  //Specific for Split Temperature:
  float tempden, temptemp, mach1,
    postT, preT,
    temprat;
  int centerfound;

  int DensNum, TENum, GENum, 
    Vel1Num, Vel2Num, Vel3Num;

  int MachNum, PSTempNum, PSDenNum;

  /* Compute size (in floats) of the current grid. */
  int size = 1;
  for (int dim = 0; dim < GridRank; dim++)
    size *= GridDimension[dim];
  
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
				       Vel3Num, TENum) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }
   
  // Get CR species fields.

  if (IdentifyShockSpeciesFields(MachNum,PSTempNum,PSDenNum) == FAIL) {
    ENZO_FAIL("Error in IdentifyShockSpeciesFields.");
  }

  /* Get easy to handle pointers for each variable. */
 
  float *density     = BaryonField[DensNum];
  float *velocity1   = BaryonField[Vel1Num];
  float *velocity2   = BaryonField[Vel2Num];
  float *velocity3   = BaryonField[Vel3Num];
  float *mach        = BaryonField[MachNum];
  float *pstemp      = BaryonField[PSTempNum];
  float *psden       = BaryonField[PSDenNum];

  /* Create temperature, entropy fields */
  float *entropy = new float[size];
  float *tempgrad_dot_entropygrad = new float[size];
  double *flowdivergence = new double[size];
    /* If using cosmology, compute the expansion factor and get units. */
  
  float *temperature = new float[size];
  
  if (this->ComputeTemperatureField(temperature) == FAIL){
    ENZO_FAIL("Error in grid->ComputeTemperatureField.");
  }
  
  float TemperatureUnits = 1, DensityUnits = 1, LengthUnits = 1,
    VelocityUnits = 1, TimeUnits = 1;

  if (GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	       &TimeUnits, &VelocityUnits, Time) == FAIL) {
    ENZO_FAIL("Error in GetUnits.");
  }
 
  // calculate cell temperature, set default values for mach
  for (i=0;i<size;i++){
    entropy[i] = tempgrad_dot_entropygrad[i] = mach[i] = 0.0;
    flowdivergence[i] = (double)(0.0);
    entropy[i] = temperature[i] / (pow(density[i],(Gamma - 1.0)));
    mach[i] = 0.0;
    if(StorePreShockFields){
      pstemp[i] = 0.0;
      psden[i] = 0.0;
    }
  }
  
  //Calculate temperature gradient dotted with entropy gradient
  //Calculate the flow divergence.
  //Do this for all cells except last ghost zone.
  int kstart = 1;
  int kend = GridDimension[2]-1;
  int jstart = 1;
  int jend = GridDimension[1]-1;
  if(GridRank < 3){
    kstart=0;
    kend=1;
  }
  if(GridRank < 2){
    jstart=0;
    jend=1;
  }

  for(k=kstart;k<kend;k++){
    for(j=jstart; j<jend;j++){
      for(i=1; i<GridDimension[0]-1;i++){
	
	index = i + GridDimension[0]*(j + GridDimension[1]*k);

	tempgrad_dot_entropygrad[index] = inv2dx2*
	  ((max(ShockTemperatureFloor,temperature[index+1])-max(ShockTemperatureFloor,temperature[index-1]))*
	   (entropy[index+1]-entropy[index-1]));
	if (GridRank > 1)
	  tempgrad_dot_entropygrad[index] += inv2dx2*
	    ((max(ShockTemperatureFloor,temperature[index+GridDimension[0]])-
	      max(ShockTemperatureFloor,temperature[index-GridDimension[0]]))*
	     (entropy[index+GridDimension[0]]-
	      entropy[index-GridDimension[0]]));
	if (GridRank > 2)
	  tempgrad_dot_entropygrad[index] += inv2dx2*
	    ((max(ShockTemperatureFloor,temperature[index+GridDimension[0]*GridDimension[1]])-
	      max(ShockTemperatureFloor,temperature[index-GridDimension[0]*GridDimension[1]]))*
	     (entropy[index+GridDimension[0]*GridDimension[1]]-
	      entropy[index-GridDimension[0]*GridDimension[1]]));
	
	flowdivergence[index] = (double)(inv2dx)*
	  ((double)(velocity1[index+1]) - (double)(velocity1[index-1]));
	if (GridRank > 1)
	  flowdivergence[index] += (double)(inv2dx)*
	    ((double)(velocity2[index+GridDimension[0]])- 
	     (double)(velocity2[index-GridDimension[0]]));   
	if (GridRank > 2)
	  flowdivergence[index] += (double)(inv2dx)*
	    ((double)(velocity3[index+GridDimension[0]*GridDimension[1]])- 
	     (double)(velocity3[index-GridDimension[0]*GridDimension[1]]));
      }
    }
  }
  
//Pseudo-Code
/* -----------------Pseudo Code------------------- /
   Check a cell if it is is a shock
   If not, continue
   If yes, then determine temperature gradient
   Search for pre/postshock cell
   If we come to a higher divergence, break
   If we come to a local minimum in divergence, make it the pre/post-shock
      cell, and stop looking
   If we come to a non-negative divergence, make it the pre/post-shock
      cell, and stop looking
   Make sure that the temperature/density keep increasing or decreasing
   After the ends are found, compute the mach number.
   Then compute the CR energy injected, multiplied by the timestep!

   In order to get statistics in terms of pre-shock quantities later, put 
     the shock in the pre-shock cell(i.e. mach number, cr...)
   Done!
/ ----------------------------------------------- */

//Different starting/ending to make sure we examine all pairs.  For 2 or 1D calculations, we only want to analyze k=j=0
  if(GridRank < 3){
    ks=1;
    ke=-1;
  }
  if(GridRank < 2){
    js=1;
    je=-1;
  }
  for(k=ks-1; k<=ke+1;k++){
    for(j=js-1; j<=je+1;j++){
      for(i=is-1; i<=ie+1;i++){
	
 	index = i + GridDimension[0]*(j + GridDimension[1]*k);

	if(tempgrad_dot_entropygrad[index] <= 0.0 ||
	   flowdivergence[index] >= 0.0)
	  continue;

	//x-direction

	posti = index;
	tempden = density[index];
	maxdiv = flowdivergence[index];
	temptemp = temperature[index];
	while(true){
	  if(temperature[index+1] >= temperature[index-1]){
	    posti++;
	  }else{
	    posti--; 
	  }
	  //Make sure temperature/density keeps increasing
	  if(temperature[posti] < temptemp ||
	     density[posti] < tempden){
	    break; 
	  }else{
	    temptemp = temperature[posti];
	    tempden = density[posti];
	  }
	    
	  if(flowdivergence[posti] < maxdiv || 
	     flowdivergence[posti] >= 0.0){
	    postT = temperature[posti];
	    break;
	  }
	  if(flowdivergence[posti] < flowdivergence[index]){
	    centerfound = 1;
	    break;
	  }
	  maxdiv = flowdivergence[posti];
	}
	prei = index;
	tempden = density[index];
	temptemp = temperature[index];
	while(true){
	  if(temperature[index+1] < temperature[index-1]){
	    prei++;
	  }else{
	    prei--; 
	  }
	  //Make sure temperature/density keeps decreasing
	  if(temperature[prei] > temptemp ||
	     density[prei] > tempden){
	    break; 
	  }else{
	    temptemp = temperature[prei];
	    tempden = density[prei];
	  }
	  if(flowdivergence[prei] < maxdiv || 
	     flowdivergence[prei] >= 0.0){
	    preT = temperature[prei];
	    break;
	  }
	  if(flowdivergence[prei] < flowdivergence[index]){
	    centerfound = 1;
	    break;
	  }
	  maxdiv = flowdivergence[prei];
	}
	if(centerfound !=1 && 
	   density[posti]>density[prei] &&
	   temperature[posti]>temperature[prei]){
	  temprat = max(postT,ShockTemperatureFloor)/(max(preT,ShockTemperatureFloor));
	  
	  mach[index] = 
	    sqrt(( 8.0*temprat - 7.0e0 + 
		   sqrt( (7.0e0 - 8.0e0*temprat)*(7.0e0 - 8.0e0*temprat) 
			 + 15.0e0) )/5.0e0);  
	  
	}

	
	//y-direction
	if(GridRank > 1){
	  posti = index;
	  tempden = density[index];
	  temptemp = temperature[index];
	  while(true){
	    maxdiv = flowdivergence[index];
	    if(temperature[index+GridDimension[0]] >= 
	       temperature[index-GridDimension[0]]){
	      posti+=GridDimension[0];
	    }else{
	      posti-=GridDimension[0]; 
	    }
	    //Make sure temperature/density keeps increasing
	    if(temperature[posti] < temptemp ||
	       density[posti] < tempden){
	      break; 
	    }else{
	      temptemp = temperature[posti];
	      tempden = density[posti];
	    }

	    if(flowdivergence[posti] < maxdiv || 
	       flowdivergence[posti] >= 0.0){
	      postT = temperature[posti];
	      break;
	    }
	    if(flowdivergence[posti] < flowdivergence[index]){
	      centerfound = 1;
	      break;
	    }
	    maxdiv = flowdivergence[posti];
	  }
	  prei = index;
	  tempden = density[index];
	  temptemp = temperature[index];
	  while(true){
	    if(temperature[index+GridDimension[0]] < 
	       temperature[index-GridDimension[0]]){
	      prei+=GridDimension[0];
	    }else{
	      prei-=GridDimension[0]; 
	    }
	    //Make sure temperature/density keeps decreasing
	    if(temperature[prei] > temptemp ||
	       density[prei] > tempden){
	      break; 
	    }else{
	      temptemp = temperature[prei];
	      tempden = density[prei];
	    }

	    if(flowdivergence[prei] < maxdiv || 
	       flowdivergence[prei] >= 0.0){
	      preT = temperature[prei];
	      break;
	    }
	    if(flowdivergence[prei] < flowdivergence[index]){
	      centerfound = 1;
	      break;
	    }
	    maxdiv = flowdivergence[prei];
	  }
	  if(centerfound !=1 && 
	   density[posti]>density[prei] &&
	   temperature[posti]>temperature[prei]){
	    temprat = max(postT,ShockTemperatureFloor)/(max(preT,ShockTemperatureFloor));
	  
	    mach1 = 
	      sqrt(( 8.0*temprat - 7.0e0 + 
		     sqrt( (7.0e0 - 8.0e0*temprat)*(7.0e0 - 8.0e0*temprat) 
			   + 15.0e0) )/5.0e0);  
	    if(mach1 > mach[index])
	      mach[index]=mach1;
	  }
	}
	
	//z-direction
	if(GridRank > 2){
	  posti = index;
	  tempden = density[index];
	  temptemp = temperature[index];
	  while(true){
	    maxdiv = flowdivergence[index];
	    if(temperature[index+GridDimension[0]*GridDimension[1]] >= 
	       temperature[index-GridDimension[0]*GridDimension[1]]){
	      posti+=GridDimension[0]*GridDimension[1];
	    }else{
	      posti-=GridDimension[0]*GridDimension[1]; 
	    }
	    //Make sure temperature/density keeps increasing
	    if(temperature[posti] < temptemp ||
	       density[posti] < tempden){
	      break; 
	    }else{
	      temptemp = temperature[posti];
	      tempden = density[posti];
	    }

	    if(flowdivergence[posti] < maxdiv || 
	       flowdivergence[posti] >= 0.0){
	      postT = temperature[posti];
	      break;
	    }
	    if(flowdivergence[posti] < flowdivergence[index]){
	      centerfound = 1;
	      break;
	    }
	    maxdiv = flowdivergence[posti];
	  }
	  prei = index;
	  tempden = density[index];
	  temptemp = temperature[index];
	  while(true){
	    if(temperature[index+GridDimension[0]*GridDimension[1]] < 
	       temperature[index-GridDimension[0]*GridDimension[1]]){
	      prei+=GridDimension[0]*GridDimension[1];
	    }else{
	      prei-=GridDimension[0]*GridDimension[1]; 
	    }
	    //Make sure temperature/density keeps decreasing
	    if(temperature[prei] > temptemp ||
	       density[prei] > tempden){
	      break; 
	    }else{
	      temptemp = temperature[prei];
	      tempden = density[prei];
	    }

	    if(flowdivergence[prei] < maxdiv || 
	       flowdivergence[prei] >= 0.0){
	      preT = temperature[prei];
	      break;
	    }
	    if(flowdivergence[prei] < flowdivergence[index]){
	      centerfound = 1;
	      break;
	    }
	    maxdiv = flowdivergence[prei];
	  }
	  if(centerfound !=1 && 
	   density[posti]>density[prei] &&
	   temperature[posti]>temperature[prei]){
	    temprat = max(postT,ShockTemperatureFloor)/(max(preT,ShockTemperatureFloor));
	  
	    mach1 = 
	      sqrt(( 8.0*temprat - 7.0e0 + 
		     sqrt( (7.0e0 - 8.0e0*temprat)*(7.0e0 - 8.0e0*temprat) 
			   + 15.0e0) )/5.0e0);  
	    if(mach1 > mach[index])
	      mach[index]=mach1;
	  }
	}  // GridRank > 2
	if(StorePreShockFields){
	  pstemp[index] = max(temperature[prei],ShockTemperatureFloor);
	  psden[index] = density[prei];
	}
      } // for i
    } // for j
  } // for k
  
  /* deallocate temporary space */
  
  delete [] temperature;
  delete [] flowdivergence;
  delete [] tempgrad_dot_entropygrad;
  delete [] entropy;

  return SUCCESS;
}
