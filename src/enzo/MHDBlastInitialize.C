#include <math.h>
#include <string.h>
#include <stdio.h>
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

  fprintf(stderr, "====================================== \n");
  fprintf(stderr, "====================================== \n");
  fprintf(stderr, "========= MHDBlastInitialize ========= \n");  
  fprintf(stderr, "====================================== \n");
  fprintf(stderr, "======== MyProcessorNumber %"ISYM"  ========= \n",
	  MyProcessorNumber);
  fprintf(stderr, "====================================== \n");
  fprintf(stderr, "====================================== \n");


  //
  //
  // Labels and Units.  (For IO.)
  // 
  
  char *DensName = "Density";
  
  char *TEName = "TotalEnergy";
  char *Vel1Name = "x-velocity";
  char *Vel2Name = "y-velocity";
  char *Vel3Name = "z-velocity";
  
  if(DualEnergyFormalism ){
    char *GEName = "GasEnergy";
    DataLabel[5] = GEName;
    DataUnits[5] = NULL;   
  }



  if( EquationOfState == 0 ){

    DataLabel[0] = DensName;
    DataLabel[1] = TEName;
    DataLabel[2] = Vel1Name;
    DataLabel[3] = Vel2Name;
    DataLabel[4] = Vel3Name;

  }else if (EquationOfState == 1){
    DataLabel[0] = DensName;
    DataLabel[1] = Vel1Name;
    DataLabel[2] = Vel2Name;
    DataLabel[3] = Vel3Name;
  }

  DataUnits[0] = NULL;
  DataUnits[1] = NULL;
  DataUnits[2] = NULL;
  DataUnits[3] = NULL;
  DataUnits[4] = NULL;
  
  
  MHDLabel[0] = "MagneticField_F_1";
  MHDLabel[1] = "MagneticField_F_2";
  MHDLabel[2] = "MagneticField_F_3";
  
  MHDcLabel[0] = "MagneticField_C_1";
  MHDcLabel[1] = "MagneticField_C_2";
  MHDcLabel[2] = "MagneticField_C_3";
  
  MHDeLabel[0] = "ElectricField_1";
  MHDeLabel[1] = "ElectricField_2";
  MHDeLabel[2] = "ElectricField_3";
  
  MHDUnits[0] = "FourPiGauss";
  MHDUnits[1] = "FourPiGauss";
  MHDUnits[2] = "FourPiGauss";
  
  MHDeUnits[0] = "FourPiGauss";
  MHDeUnits[1] = "FourPiGauss";
  MHDeUnits[2] = "FourPiGauss";
  

  // General controll variable
  //

  int dim;

  //
  // Parameters and their defaults.
  // 


  char line[MAX_LINE_LENGTH];
  int ret = 0, GasFlag = 0, Pflag=0, TotalFlag=0;
  int ObsFlag = 0;
  int RefineOnStartup = FALSE;
  // Or 1.0
  float fpi = 1.0;//1/(sqrt(4.0*3.14159265)*4.0*3.14159265);
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
  //now this is a global value.  Becuase I'm a jerk.

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

    //I changed some nominclature, to make things easier on myself
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
    
    /////
    ObsFlag  += sscanf(line, "MHDBlastD0 = %"PSYM, &Density0);
    ObsFlag += sscanf(line, "MHDBlastD1 = %"PSYM, &Density1);
    
    ObsFlag += sscanf(line, "MHDBlastB0 = %"PSYM" %"PSYM" %"PSYM, B0, B0+1, B0+2);
    ObsFlag+= sscanf(line, "MHDBlastB1 = %"PSYM" %"PSYM" %"PSYM, B1, B1+1, B1+2);
    
    ObsFlag += sscanf(line, "MHDBlastP0 = %"PSYM, &Pressure0);
    ObsFlag += sscanf(line, "MHDBlastP1 = %"PSYM, &Pressure1);

    ObsFlag += sscanf(line, "MHDBlastGasEnergy0 = %"PSYM, &GasEnergy0);
    ObsFlag += sscanf(line, "MHDBlastGasEnergy1 = %"PSYM, &GasEnergy1);

    ObsFlag += sscanf(line, "MHDBlastTotalEnergy0 = %"PSYM, &TotalEnergy0);
    ObsFlag += sscanf(line, "MHDBlastTotalEnergy1 = %"PSYM, &TotalEnergy1);

    if( ObsFlag != 0 ){
      fprintf(stderr," Obsolete nomenclature in blast init file..  Fix.  \n");
      fprintf(stderr," I changed all the '0' to 'A', all the '1' to 'B'\n");
      fprintf(stderr," so MHDBlastD0 is now MHDBlastDA\n");
      return FAIL;
    }

    if( sscanf(line, "PerturbAmplitude      = %"PSYM, &PerturbAmplitude) != 0 ||
	sscanf(line, "PerturbMethod     = %"ISYM"", &PerturbMethod) != 0           ||
	sscanf(line, "PerturbWavelength = %"PSYM,PerturbWavelength)      != 0 ){
      fprintf(stderr," Parameter renamed: PerturbAmplitude(,Method,Wavelength) -> MHDBlastPerturbAmplitude(etc).  Fix it.\n");
      return FAIL;
    }

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
    ret += sscanf(line, "MHDBlastPerturbMethod     = %"ISYM"", &PerturbMethod);
    ret += sscanf(line, "MHDBlastPerturbWavelength = %"PSYM" %"PSYM" %"PSYM,
                  PerturbWavelength,PerturbWavelength+1,PerturbWavelength+2);

    ret += sscanf(line, "MHDBlastRefineOnStartup  = %"ISYM"", &RefineOnStartup);

#ifdef Unsupported
    //see comment below
    ret += sscanf(line, "MHDBlastNormal = %"PSYM" %"PSYM" %"PSYM,
		  MHDBlastNormal, MHDBlastNormal+1,MHDBlastNormal+2);
#endif
  }//line loop

  //
  //Re scale the subgrid edges to line up with the parent grid.
  // nCellsL and nCellsR are the number of cells from the domain left edge.
  //

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
  // Long Dimension is used to conver the radius from Physical units to Grid Units;
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

    GasEnergyA = Pressure0/(Gamma-1);
    GasEnergyB = Pressure1/(Gamma-1);
    if( useMHDCT != 1 ){
      GasEnergyA = Pressure0/((Gamma-1)*DensityA);
      GasEnergyB = Pressure1/((Gamma-1)*DensityB);
    }

  }

  //The variable stored is Gas+Kinetic+Magnetic Energy.
  if( GasFlag > 0 || Pflag > 0){
    if( useMHDCT == 1 ){
      Energy0 = GasEnergyA + 
	0.5*DensityA*(VelocityA[0]*VelocityA[0] + VelocityA[1]*VelocityA[1] + VelocityA[2]*VelocityA[2])
	+0.5*(BA[0]*BA[0]+BA[1]*BA[1]+BA[2]*BA[2]);
      Energy1 = GasEnergyB + 
	0.5*DensityB*(VelocityB[0]*VelocityB[0] + VelocityB[1]*VelocityB[1] + VelocityB[2]*VelocityB[2])
	+0.5*(BB[0]*BB[0]+BB[1]*BB[1]+BB[2]*BB[2]);
    }else{
      Energy0 = GasEnergyA + 
	0.5*(VelocityA[0]*VelocityA[0] + VelocityA[1]*VelocityA[1] + VelocityA[2]*VelocityA[2]);
      Energy1 = GasEnergyB + 
	0.5*(VelocityB[0]*VelocityB[0] + VelocityB[1]*VelocityB[1] + VelocityB[2]*VelocityB[2]);
    }
  }
  //<dbg>
  //fprintf(stderr,"klown: pa %"GSYM" pb %"GSYM" ea %"GSYM" eb %"GSYM"\n",
  //Pressure0, Pressure1, Energy0, Energy1);
  //</dbg>
  if( TotalFlag > 0){
    Energy0=TotalEnergyA;
    Energy1=TotalEnergyB;
    
  }

#ifdef Unsupported
  //Oblique shocks need special treatment at the boundary.  Initially this was done in a very
  //unclean way.  
  if( InitStyle == 10 ){
    //
    // Set up normal vector.  
    // Normal[0] * x + Normal[1]*y + Normal[2]*z + Normal[3] = 0
    // x,y,z in Code Units from Domain Wall.
    // The plane is defined to go throught the MHDBlastCenter, so

    //<dbg>
    fprintf(stderr,"normal before %"GSYM" %"GSYM" %"GSYM"\n", MHDBlastNormal[0], MHDBlastNormal[1], MHDBlastNormal[2]);
    //</dbg>


    //in MHD_ObliqueRoutines
    if( SetupNormal(MHDBlastNormal, MHDBlastCenter, MetaData) == FAIL ){
      fprintf(stderr,"Failure in MHDBlastNormalSetup.\n");
      return FAIL;
    }

    //<dbg>
    fprintf(stderr,"normal after %"GSYM" %"GSYM" %"GSYM" %"GSYM"\n", 
	    MHDBlastNormal[0], MHDBlastNormal[1], MHDBlastNormal[2], MHDBlastNormal[3]);
    //</dbg>

    //The initial conditions get rotated, assuming they were initially alligned 
    //with the X axis, in the cyclic order (so VelocityA = Vx, Vy, Vz)
    //The Vx axis is now aligned with the normal component.
    //No rotation about the normal is performed.
    
    RotateVector( VelocityA, MHDBlastNormal);
    RotateVector( VelocityB, MHDBlastNormal);
    RotateVector( BA, MHDBlastNormal);
    RotateVector( BB, MHDBlastNormal);
    
  }//rotated vector.
#else
  float MHDBlastNormal[4] = {0.0, 0.0, 0.0, 0.0}; //3 cmpts of normal for rotation + offset.
#endif //Unsupported
  //
  // Initialize the top grid.  Cant' decide if I want a uniform grid here or MHDBlastInitialize.
  //

  if( TopGrid.GridData->MHDBlastInitializeGrid(DensityA, DensityB,
					       Energy0,  Energy1,
					       VelocityA, VelocityB,
					       BA, BB, 
					       Radius, MHDBlastCenter, LongDimension,
					       PerturbAmplitude, PerturbMethod,PerturbWavelength,
					       InitStyle, MHDBlastNormal) == FAIL )
    {
      fprintf(stderr, "MHDBlastInitialize: Shit.  Error in MHDBlastInitializeGrid.\n");
      return FAIL;
    }


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
	  SubgridDims[dim] = NumberOfSubgridZones[dim] + 2*DEFAULT_GHOST_ZONES;
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
							   InitStyle, MHDBlastNormal) == FAIL )
	  {
	    fprintf(stderr, "MHDBlastInitialize: Shit.  Error in MHDBlastInitializeGrid.\n");
	    return FAIL;
	  }


	
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
	  == FAIL) {
	fprintf(stderr, "Error in ProjectSolutionToParentGrid.\n");
	return FAIL;
      }
    
    // set up the root grid 
    
    if (MaximumRefinementLevel > 0) {
      if (Subgrid[0]->GridData->ProjectSolutionToParentGrid(*(TopGrid.GridData))
	  == FAIL) {
	fprintf(stderr, "Error in ProjectSolutionToParentGrid.\n");
	return FAIL;
      }
    }
    
    
    // Put the projection options back to the inital.
    MHD_ProjectE = MHD_ProjectEtmp;
    MHD_ProjectB = MHD_ProjectBtmp;
    
    
  }//RefineOnStartup

  return SUCCESS;
}








