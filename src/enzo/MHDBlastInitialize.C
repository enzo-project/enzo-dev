/***********************************************************************
/
/  Initialize the MHD Blast Test.
/
/  written by: David Collins
/  date:       2004-2013
/  modified1:
/
/  PURPOSE:  MHDBlast Test is a general purpose 2 state problem initializer
/
/  PARAMETERS:   
/ 
/       MHDBlastInitStyle  Shape of discontinuity
/                          0 = sphere: 
/                          1,2,3=rectangular slice along x, y,z 
/                          40,41,42=cylander along x,y,z:  
/                          5,6,7,8 = Index Tests: 10*i + 100*j + 1000*k, i,j,k 
/
/
/       MHDBlastCenter     Center in spatial units. 
/       MHDBlastRadius     in space units OF THE LONGEST AXIS, 
/                          MHDBlastInitStyle = 1,2,3 the width of the slab 
/                          MHDBlastInitStyle = 40,41,42 the radius of the infinite cylander
/                          MHDBlastInitStyle = 0, radius of sphere
/
/       Density, Pressure, Magnetic Field, Velocity can be set with the following.  
/       For all fields, one side of the discontinuity is denoted A, one is B.
/       MHDBlastD[A,B] 
/       MHDBlastVelocity[A,B] 
/       MHDBlastB[A,B]
/       MHDBlastGasEnergy[A,B]
/       MHDBlastP[A,B]   (Gas energy is checked first, then pressure)
/       
/       Fields may be perturbed to seed instabilities or linear characteristic advection using
/       the following:
/       MHDBlastPerturbAmplitude  
/       MHDBlastPerturbWavelength     
/       MHDBlastPerturbMethod    
/                               1:  white noise in the velocity
/                               2:  plane symmetric noise in the velocity
/                               7:  plane symmetric, energy preserving
/                               4,5,6: sinusoidal perturbation in vy
/                               8:  Perturbation for RT instability, as in Stone & Gardiner 2007
/                               80: square wave in the left moving fast characteristic
/                               81: square wave in the left moving alfven characteristic
/                               82: square wave in the left moving slow characteristic
/                               83: square wave in the contact discontinuity
/                               84: square wave in the right moving fast characteristic
/                               85: square wave in the right moving alfven characteristic
/                               86: square wave in the right moving slow characteristic
/                               70: sine wave in the left moving fast characteristic
/                               71: sine wave in the left moving alfven characteristic
/                               72: sine wave in the left moving slow characteristic
/                               73: sine wave in the contact discontinuity
/                               74: sine wave in the right moving fast characteristic
/                               75: sine wave in the right moving alfven characteristic
/                               76: sine wave in the right moving slow characteristic
/
/
/
/  RETURNS:
/    SUCCESS or FAIL
/
************************************************************************/

#include <math.h>
#include <string.h>
#include <stdio.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"

/* This initializes genearlized discontinuities.  Good for blast waves, Kelving Helmholtz, shock tubes.  Probably others.

Code flow:
1.) Declare/Define Defaults for  parameters.
2.) Read parameters from file.
3.) Calculate TotalEnergy from quantity given (options are Total Energy, Gas Energy, Pressure.)
4.) Set up data labels, units.  
5.) Declare Hierarchy Object.
6.) Define linked list
7.) Call initializer on each level.
8.) Project to parent.

*/

//in MHD_ObliqueRoutines.C
void RotateVector( float * Vector, float * Normal);
int SetupNormal(float Normal[], float MHDBlastCenter[3], TopGridData & MetaData);

int MHDBlastInitialize(FILE *fptr, FILE *Outfptr, HierarchyEntry &TopGrid,
		       TopGridData &MetaData, ExternalBoundary &Exterior)
{

  //
  //
  // Labels and Units.  (For IO.)
  // 
  
  char *DensName = "Density";
  char *TEName = "TotalEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  char *GPotName = "GravPotential";
  char *BxName = "Bx";
  char *ByName = "By";
  char *BzName = "Bz";
  char *PhiName = "Phi";
  
  int i=0, j=0;

    DataLabel[i++] = DensName;
    DataUnits[j++] = NULL;
  if( EquationOfState == 0 ){
    DataLabel[i++] = TEName;
    DataUnits[j++] = NULL;
  }
    DataLabel[i++] = Vel1Name;
    DataUnits[j++] = NULL;
    DataLabel[i++] = Vel2Name;
    DataUnits[j++] = NULL;
    DataLabel[i++] = Vel3Name;
    DataUnits[j++] = NULL;
  if( UseMHD ){
    DataUnits[i] = NULL;
    DataLabel[i++] = BxName;
    DataUnits[i] = NULL;
    DataLabel[i++] = ByName;
    DataUnits[i] = NULL;
    DataLabel[i++] = BzName;
    DataUnits[i] = NULL;
    DataLabel[i++] = PhiName;
  }


  if(DualEnergyFormalism ){
    char *GEName = "GasEnergy";
    DataLabel[i++] = GEName;
    DataUnits[j++] = NULL;   
  }

  if (WritePotential){
    DataLabel[i++] = GPotName;
    DataUnits[j++] = NULL;
  }
  
  if ( UseMHDCT ){
  MHDLabel[0] = "BxF";
  MHDLabel[1] = "ByF";
  MHDLabel[2] = "BzF";
  
  MHDcLabel[0] = "Bx";
  MHDcLabel[1] = "By";
  MHDcLabel[2] = "Bz";
  
  MHDeLabel[0] = "Ex";
  MHDeLabel[1] = "Ey";
  MHDeLabel[2] = "Ez";
  
  MHDUnits[0] = "None";
  MHDUnits[1] = "None";
  MHDUnits[2] = "None";
  
  MHDeUnits[0] = "None";
  MHDeUnits[1] = "None";
  MHDeUnits[2] = "None";
  }
  

  // General control variable
  int dim;

  // Parameters and their defaults.
  char line[MAX_LINE_LENGTH];
  int ret = 0, GasFlag = 0, Pflag=0, TotalFlag=0;
  int ObsFlag = 0;
  int RefineOnStartup = FALSE;

  float DensityA = 1.0666,
    DensityB = 1.0,
    GasEnergyA = 3.666,
    GasEnergyB = 1000.666,
    TotalEnergyA = 1.0,
    TotalEnergyB = 1.0;

  float PressureA, PressureB;
  float VelocityA[3] = {0.666, 0.666, 0.666};
  float VelocityB[3] = {0.3666, 0.3666, 0.3666};
  float BA[3]  = {0.5, 0.0, 0.0};
  float BB[3]  = {0.5, 0.0, 0.0};
  float EnergyA, EnergyB;
  float Radius = 4.0;
  int InitStyle = 0, PerturbMethod = -1;
  float MHDBlastCenter[3] = {0.5,0.5,0.5};
  float PerturbAmplitude = 0.0;
  float PerturbWavelength[3] = {0.0, 0.0, 0.0};

  FLOAT MHDBlastSubgridLeft[3]  = {DomainLeftEdge[0] ,DomainLeftEdge[1] ,DomainLeftEdge[2]};
  FLOAT MHDBlastSubgridRight[3] = {DomainRightEdge[0],DomainRightEdge[1],DomainRightEdge[2]};
  
  //Obsolete variable names.
  float Pressure0, Pressure1;
  float B0[3],B1[3],Energy0, Energy1;
  float Density0,Density1, GasEnergy0, GasEnergy1, TotalEnergy0,TotalEnergy1;
    
  //
  // Read Parameter File.
  //

  while (fgets(line, MAX_LINE_LENGTH, fptr) != NULL) {

    ret = 0;

    //I changed some nomenclature, to make things easier on myself
    //This checks for old nominclature.
    ObsFlag = 0;
    
    ret += sscanf(line, "MHDBlastDA = %"PSYM, &DensityA);
    ret += sscanf(line, "MHDBlastDB = %"PSYM, &DensityB);

    ret += sscanf(line, "MHDBlastBA = %"PSYM" %"PSYM" %"PSYM, BA, BA+1, BA+2);
    ret += sscanf(line, "MHDBlastBB = %"PSYM" %"PSYM" %"PSYM, BB, BB+1, BB+2);

    ret += sscanf(line, "MHDBlastVelocityA = %"PSYM" %"PSYM" %"PSYM,
		  VelocityA, VelocityA+1,VelocityA+2);
    ret += sscanf(line, "MHDBlastVelocityB = %"PSYM" %"PSYM" %"PSYM,
		  VelocityB, VelocityB+1,VelocityB+2);
    
    Pflag += sscanf(line, "MHDBlastPA = %"PSYM, &Pressure0);
    Pflag += sscanf(line, "MHDBlastPB = %"PSYM, &Pressure1);

    GasFlag += sscanf(line, "MHDBlastGasEnergyA = %"PSYM, &GasEnergyA);
    GasFlag += sscanf(line, "MHDBlastGasEnergyB = %"PSYM, &GasEnergyB);

    TotalFlag += sscanf(line, "MHDBlastTotalEnergyA = %"PSYM, &TotalEnergyA);
    TotalFlag += sscanf(line, "MHDBlastTotalEnergyB = %"PSYM, &TotalEnergyB);
    
    ////

    ret += sscanf(line, "MHDBlastRadius = %"PSYM, &Radius);

    ret += sscanf(line, "MHDBlastInitStyle = %"ISYM"", &InitStyle);

    ret += sscanf(line, "MHDBlastCenter = %"PSYM" %"PSYM" %"PSYM,
		  MHDBlastCenter, MHDBlastCenter+1,MHDBlastCenter+2);

    ret += sscanf(line, "MHDBlastSubgridLeft  = %"PSYM" %"PSYM" %"PSYM,
		  MHDBlastSubgridLeft, MHDBlastSubgridLeft +1 , MHDBlastSubgridLeft +2);
    ret += sscanf(line, "MHDBlastSubgridRight = %"PSYM" %"PSYM" %"PSYM,
		  MHDBlastSubgridRight, MHDBlastSubgridRight +1 , MHDBlastSubgridRight +2);

    ret += sscanf(line, "MHDBlastPerturbAmplitude      = %"PSYM, &PerturbAmplitude);
    ret += sscanf(line, "MHDBlastPerturbMethod         = %"ISYM"", &PerturbMethod);
    ret += sscanf(line, "MHDBlastPerturbWavelength      = %"PSYM" %"PSYM" %"PSYM,
                  PerturbWavelength,PerturbWavelength+1,PerturbWavelength+2);

    ret += sscanf(line, "MHDBlastRefineOnStartup  = %"ISYM"", &RefineOnStartup);

  }//line loop

  //Re scale the subgrid edges to line up with the parent grid.
  // nCellsL and nCellsR are the number of cells from the domain left edge.

  int nCellsL[3],nCellsR[3];
  int nCells[3] = {0,0,0};  
  for( dim = 0; dim < 3; dim++){
    nCellsL[dim]= nint(( MHDBlastSubgridLeft[dim] - DomainLeftEdge[dim] )/
		     (DomainRightEdge[dim]-DomainLeftEdge[dim])*MetaData.TopGridDims[dim]);

    MHDBlastSubgridLeft[dim]=max( nCellsL[dim]*(DomainRightEdge[dim]-DomainLeftEdge[dim])/MetaData.TopGridDims[dim],
				  DomainLeftEdge[dim]);
    
    nCellsR[dim] = nint(( MHDBlastSubgridRight[dim] - DomainLeftEdge[dim] )/
			(DomainRightEdge[dim]-DomainLeftEdge[dim])*MetaData.TopGridDims[dim]);
    
    MHDBlastSubgridRight[dim] = min( nCellsR[dim]*(DomainRightEdge[dim]-DomainLeftEdge[dim])/MetaData.TopGridDims[dim],
				     DomainRightEdge[dim]);
    nCells[dim] =  nint( (MHDBlastSubgridRight[dim]-MHDBlastSubgridLeft[dim])/
      (DomainRightEdge[dim]-DomainLeftEdge[dim])*MetaData.TopGridDims[dim] );
    
  }

  if( RefineOnStartup == 1 ){
    fprintf(stderr,"Subgrid Left %"GSYM" %"GSYM" %"GSYM"\n", MHDBlastSubgridLeft[0], MHDBlastSubgridLeft[1], MHDBlastSubgridLeft[2]);
    fprintf(stderr,"Subgrid Right %"GSYM" %"GSYM" %"GSYM"\n",MHDBlastSubgridRight[0],MHDBlastSubgridRight[1],MHDBlastSubgridRight[2]);
    fprintf(stderr,"nCells %"ISYM" %"ISYM" %"ISYM"\n", nCells[0], nCells[1], nCells[2]);
  }

  // Long Dimension is used to convert the radius from Physical units to Grid Units;
  // We want the axis, though, so figure out which is the longest edge (in Grid Units) 
  // then figure out which one it is.  A more elegant solution would be welcomed.

  int LongDimension = 0;
  LongDimension = (nCells[0] > nCells[1] ) ? nCells[0] : nCells[1];
  LongDimension = (LongDimension > nCells[2] ) ? LongDimension : nCells[2];
  for( dim=0; dim<3; dim++)
    if( LongDimension == nCells[dim] ){
      LongDimension = dim;
      break;
    }


  //
  // Calculate Total Energy.
  //


  if( Pflag > 0 ) {
      GasEnergyA = Pressure0/((Gamma-1)*DensityA); 
      GasEnergyB = Pressure1/((Gamma-1)*DensityB); 
  }
     

  //The variable stored is Gas+Kinetic+Magnetic Energy.
  if( GasFlag > 0 || Pflag > 0){
      Energy0 = GasEnergyA + 
	0.5*(VelocityA[0]*VelocityA[0] + VelocityA[1]*VelocityA[1] + VelocityA[2]*VelocityA[2])
	+0.5*(BA[0]*BA[0]+BA[1]*BA[1]+BA[2]*BA[2])/DensityA;
      Energy1 = GasEnergyB + 
	0.5*(VelocityB[0]*VelocityB[0] + VelocityB[1]*VelocityB[1] + VelocityB[2]*VelocityB[2])
	+0.5*(BB[0]*BB[0]+BB[1]*BB[1]+BB[2]*BB[2])/DensityB;
  }

  if( TotalFlag > 0){
    Energy0=TotalEnergyA;
    Energy1=TotalEnergyB;
  }


  //
  // Initialize the top grid.  Cant' decide if I want a uniform grid here or MHDBlastInitialize.
  //

  if( TopGrid.GridData->MHDBlastInitializeGrid(DensityA, DensityB,
					       Energy0,  Energy1,
					       VelocityA, VelocityB,
					       BA, BB, 
					       Radius, MHDBlastCenter, LongDimension,
					       PerturbAmplitude, PerturbMethod,PerturbWavelength,
					       InitStyle) == FAIL )
    ENZO_FAIL("MHDBlastInitialize:  Error in MHDBlastInitializeGrid.");

  //
  // Generate Hierarchy.
  //
  if( RefineOnStartup == 1 ){
    //Create as many subgrids as there are refinement levels 
    //needed to resolve the initial explosion region upon the start-up. 
    
    HierarchyEntry ** Subgrid;
    if (MaximumRefinementLevel > 0) 
      Subgrid   = new HierarchyEntry*[MaximumRefinementLevel];
    
    //
    //Create new HierarchyEntries.  Note that 'lev' loops are only for the SUBGRIDS.
    //
    
    int lev;
    int NumberOfSubgridZones[3], SubgridDims[3];
    
    for (lev = 0; lev < MaximumRefinementLevel; lev++) 
      Subgrid[lev] = new HierarchyEntry;
    
    for (lev = 0; lev < MaximumRefinementLevel; lev++) {
      
      //Calculate number of cells on this level.
      
      for (dim = 0; dim < MetaData.TopGridRank; dim++)
	NumberOfSubgridZones[dim] = nCells[dim]*POW(RefineBy, lev + 1);
      
      fprintf(stderr,"uncle MHDBlast:: Level[%"ISYM"]: NumberOfSubgridZones[0] = %"ISYM"\n", lev+1, 
	      NumberOfSubgridZones[0]);
      
      if (NumberOfSubgridZones[0] > 0) {
	
	// fill them out 
	
	if (lev == 0)
	  TopGrid.NextGridNextLevel  = Subgrid[0];
	Subgrid[lev]->NextGridThisLevel = NULL;
	if (lev == MaximumRefinementLevel-1)
	  Subgrid[lev]->NextGridNextLevel = NULL;
	else
	  Subgrid[lev]->NextGridNextLevel = Subgrid[lev+1];
	if (lev == 0)
	  Subgrid[lev]->ParentGrid        = &TopGrid;
	else
	  Subgrid[lev]->ParentGrid        = Subgrid[lev-1];
	
	//  compute the dimensions and left/right edges for the subgrid 
	
	for (dim = 0; dim < MetaData.TopGridRank; dim++) {
	  SubgridDims[dim] = NumberOfSubgridZones[dim] + 2*NumberOfGhostZones;
	}
	
	// create a new subgrid and initialize it 
	
	Subgrid[lev]->GridData = new grid;
	Subgrid[lev]->GridData->InheritProperties(TopGrid.GridData);
	Subgrid[lev]->GridData->PrepareGrid(MetaData.TopGridRank, SubgridDims,
					    MHDBlastSubgridLeft,MHDBlastSubgridRight, 0);
	

	if( Subgrid[lev]->GridData->MHDBlastInitializeGrid(DensityA, DensityB,
							   Energy0,  Energy1,
							   VelocityA, VelocityB,
							   BA, BB, 
							   Radius, MHDBlastCenter, LongDimension,
							   PerturbAmplitude, PerturbMethod,PerturbWavelength,
							   InitStyle) == FAIL )
	  ENZO_FAIL("MHDBlastInitialize: Error in MHDBlastInitializeGrid.");
	
      }//NumberOfSubgridZones > 0
      else{
	printf("SedovBlast: single grid start-up.\n");
      }
      
    }//level
    
    // Make sure each grid has the best data with respect to the finer grids.
    // This projection juggle is to ensure that, regardless of how the hierarchy is evolved, the field gets projected
    // properly here.
    
    int MHD_ProjectEtmp = MHD_ProjectE;
    int MHD_ProjectBtmp = MHD_ProjectB;
    MHD_ProjectE=FALSE;
    MHD_ProjectB=TRUE;
    
    for (lev = MaximumRefinementLevel - 1; lev > 0; lev--)
      if (Subgrid[lev]->GridData->ProjectSolutionToParentGrid(
							      *(Subgrid[lev-1]->GridData))
	  == FAIL) 
	ENZO_FAIL("Error in ProjectSolutionToParentGrid.");
    
    // set up the root grid 
    
    if (MaximumRefinementLevel > 0) {
      if (Subgrid[0]->GridData->ProjectSolutionToParentGrid(*(TopGrid.GridData))
	  == FAIL) 
	ENZO_FAIL("Error in ProjectSolutionToParentGrid.");
    }
    
    // Put the projection options back to the inital.
    MHD_ProjectE = MHD_ProjectEtmp;
    MHD_ProjectB = MHD_ProjectBtmp;
    
  }//RefineOnStartup

  return SUCCESS;
}








