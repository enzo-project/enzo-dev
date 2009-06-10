#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include <map>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "LevelHierarchy.h"
#include "TopGridData.h"
#include "CosmologyParameters.h"
#include "AnalysisBaseClass.h"

std::map<HierarchyEntry *, int> OriginalGridID;
#ifdef USE_HDF5_GROUPS
std::map<HierarchyEntry *, int> OriginalTaskID;
#endif

void my_exit(int exit_status);
int CosmologyComputeExpansionFactor(FLOAT time, FLOAT *a, FLOAT *dadt);

void FlagLevel( HierarchyEntry *Grid, HierarchyEntry ***GridArray, int *dx,
		float cell_width );

AnalysisBaseClass::AnalysisBaseClass( TopGridData *metadata,
			    HierarchyEntry *topgrid,
			    int level,
			    FLOAT *Left,
			    FLOAT *Right):
#ifdef USE_BAD_VALUES
  G(6.67e-8),
  MPC_CM(3.0856e+24),
  MSOLAR_G(1.989e+33)
#else
  G(6.6742e-8),
  MPC_CM(3.0856776e+24),
  MSOLAR_G(1.989e+33)
#endif
{

  int i;
  MetaData = metadata;
  TopGrid = topgrid;

  /* Create some calculated data to save flops. */
  for( i=0; i < MAX_DIMENSION; i++)
    DomainDelta[i] = DomainRightEdge[i] - DomainLeftEdge[i];  

  if( Left ){
    for( i = 0; i < MAX_DIMENSION; i++ )
      AnalysisRegionLeftEdge[i] = Left[i];
  }else{
    for( i = 0; i < MAX_DIMENSION; i++ )
      AnalysisRegionLeftEdge[i] = DomainLeftEdge[i];
  }

  if( Right ){
    for( i = 0; i < MAX_DIMENSION; i++ )
      AnalysisRegionRightEdge[i] = Right[i];
  }else{
    for( i = 0; i < MAX_DIMENSION; i++ )
      AnalysisRegionRightEdge[i] = DomainRightEdge[i];
  }

  if( level < 0 ){
    MaximumAnalysisLevel = MaximumRefinementLevel;
  }else{
    MaximumAnalysisLevel = level;
  }

  CurrentRedshift = 0.0;
  if (ComovingCoordinates) {
    FLOAT a =1.0, dadt, t = topgrid->GridData->ReturnTime();
    
    if (CosmologyComputeExpansionFactor(t, &a, &dadt) == FAIL) {
      fprintf(stderr, "Error in ComputeExpansionFactor.\n");
      my_exit(EXIT_FAILURE);
    }

    CurrentRedshift = (1.0+InitialRedshift)/a - 1.0;

    /* Thought this was a good idea, but it isn't 
       if(CurrentRedshift < 0)
       CurrentRedshift = 0.0;
    */
  }
  
  FastFindGridArray = NULL;
  FastFindDeltaN = NULL;
  FastFindCellWidth = NULL;

  MyGridArray = NULL;

  GridsPerProcessor = NULL;

}

AnalysisBaseClass::~AnalysisBaseClass(){

  delete[] FastFindGridArray;
  delete[] FastFindDeltaN;
  delete[] FastFindCellWidth;
  delete[] GridsPerProcessor;
  delete[] MyGridArray;
}


void AnalysisBaseClass::PrintParameters(){

  //printf("Enzo Analysis class parameters:\n");
  printf("MaximumAnalysisLevel                 = %"ISYM"\n", MaximumAnalysisLevel);
  printf("CurrentRedshift                      = %"FSYM"\n", CurrentRedshift); 
}

void AnalysisBaseClass::PrintGridInfo(grid *Grid){

  int i;
  int grid_rank;
  int grid_dims[MAX_DIMENSION];
  FLOAT grid_left[MAX_DIMENSION];
  FLOAT grid_right[MAX_DIMENSION];

  Grid->ReturnGridInfo(&grid_rank, grid_dims, grid_left, grid_right);

  printf("Grid information:\n");
  printf("  Rank      = %"ISYM"\n", grid_rank);
  printf("  Dims      = " );
  for( i = 0; i < MAX_DIMENSION; i++ )
    printf(" %4"ISYM, grid_dims[i]);
  printf("\n" );

  printf("  LeftEdge  = " );
  for( i = 0; i < MAX_DIMENSION; i++ )
    printf(" %8.6"FSYM, grid_left[i]);
  printf("\n" );

  printf("  RightEdge = " );
  for( i = 0; i < MAX_DIMENSION; i++ )
    printf(" %8.6"FSYM, grid_right[i]);
  printf("\n" );
}

void AnalysisBaseClass::PrintTrueFalse(FILE *buffer, bool testvalue,
				  const char *fmt){

  if(testvalue){
    fprintf(buffer, fmt, "TRUE");
  }else{
    fprintf(buffer, fmt, "FALSE");
  }
}
 

void AnalysisBaseClass::LoopNumberOfGrids( HierarchyEntry *Grid,
					   int current_level){
  
  while( Grid ){
    if ( Grid->GridData->IsInVolume(AnalysisRegionLeftEdge,
				    AnalysisRegionRightEdge)){
      if(current_level < MaximumAnalysisLevel){
	LoopNumberOfGrids(Grid->NextGridNextLevel, current_level + 1);
      }    
      NumberOfGrids++;
    }
    
    Grid = Grid->NextGridThisLevel;
  }
}

int AnalysisBaseClass::NumberOfGridsInVolume(  ){
  NumberOfGrids = 0;
  LoopNumberOfGrids( TopGrid, 0 );
  return NumberOfGrids;
}

void AnalysisBaseClass::LoopCountGridsPerProcessor( HierarchyEntry *Grid,
						    int current_level){
  
  int proc_index;

  while( Grid ){
    if ( Grid->GridData->IsInVolume(AnalysisRegionLeftEdge,
				    AnalysisRegionRightEdge)){      
      if(current_level < MaximumAnalysisLevel){
	LoopCountGridsPerProcessor(Grid->NextGridNextLevel, current_level + 1);
      }
      proc_index = Grid->GridData->ReturnProcessorNumber();
      GridsPerProcessor[proc_index]++;      
    }

    Grid = Grid->NextGridThisLevel;
  }
}

void AnalysisBaseClass::CountGridsPerProcessor(){
  if(!GridsPerProcessor)
    GridsPerProcessor = new Eint32[NumberOfProcessors];
  
  for(int i = 0; i < NumberOfProcessors; i++)
    GridsPerProcessor[i] = 0;

  LoopCountGridsPerProcessor( TopGrid, 0 );
}

void AnalysisBaseClass::LoopBuildMyGridArray(HierarchyEntry *Grid,
					     int current_level,
					     int *current_index){

  while( Grid ){
    if ( Grid->GridData->IsInVolume(AnalysisRegionLeftEdge,
				    AnalysisRegionRightEdge)){      
      if(current_level < MaximumAnalysisLevel){
	LoopBuildMyGridArray(Grid->NextGridNextLevel, current_level + 1, 
			     current_index );
      }    
      if(Grid->GridData->ReturnProcessorNumber() == MyProcessorNumber){
	MyGridArray[(*current_index)++] = Grid;	
      }
    }    
    Grid = Grid->NextGridThisLevel;
  }
}


void AnalysisBaseClass::BuildMyGridArray(  ){
  int index = 0;

  if(GridsPerProcessor == NULL)
    CountGridsPerProcessor();

  delete[] MyGridArray;
  MyGridArray = new HierarchyEntry*[GridsPerProcessor[MyProcessorNumber]];
  
  LoopBuildMyGridArray( TopGrid, 0, &index);

}


void AnalysisBaseClass::FlagGridCells(HierarchyEntry *Grid){
  
  Grid->GridData->ClearFlaggingField();

  HierarchyEntry *SubGrid = Grid->NextGridNextLevel;
  
  while(SubGrid != NULL){
    if(Grid->GridData->FlagRefinedCells(SubGrid->GridData)==FAIL)
      my_exit( EXIT_FAILURE );
    SubGrid = SubGrid->NextGridThisLevel;
  }
}

void AnalysisBaseClass::SetAnalysisRegion( FLOAT Left[MAX_DIMENSION],
			  FLOAT Right[MAX_DIMENSION] ){
  
  for( int i = 0; i < MAX_DIMENSION; i++ ){
    AnalysisRegionLeftEdge[i] = Left[i];
    AnalysisRegionRightEdge[i] = Right[i];
  }
}

void AnalysisBaseClass::InitializeFastFind(){

  // If we're re-initializing, clean up previous
  if(FastFindGridArray){
    delete [] FastFindGridArray;
    delete [] FastFindDeltaN;
    delete [] FastFindCellWidth;    
  }

  int i, j, num_cells = 1;

  HierarchyEntry *this_grid, *next_level;
  
  FastFindDeltaN = new int[MAX_DIMENSION];
  FastFindCellWidth = new float[MAX_DIMENSION];
  
  // hard coded to 16 for unigrid case
  
  for( i=0; i < MAX_DIMENSION; i++){
    if(MaximumAnalysisLevel == 0){
      num_cells*= 16;
      FastFindCellWidth[i] = (DomainRightEdge[i] - DomainLeftEdge[i])/16.0;      
      FastFindDeltaN[i] = 16;
    }else{
      num_cells*= MetaData->TopGridDims[i];
      FastFindCellWidth[i] = (DomainRightEdge[i] - DomainLeftEdge[i])/((float)(MetaData->TopGridDims[i]));
      FastFindDeltaN[i] = MetaData->TopGridDims[i];
    }
  }

  FastFindGridArray = new HierarchyEntry*[num_cells];

  this_grid = TopGrid;

  // hard coded one level, since that's all we can resolve with topgriddims

  while( this_grid ){
    this_grid->GridData->FlagGridArray( &FastFindGridArray, 
					FastFindDeltaN, FastFindCellWidth, 
					this_grid );
    if(MaximumAnalysisLevel > 0){
      next_level = this_grid->NextGridNextLevel;
      while( next_level ){
	next_level->GridData->FlagGridArray( &FastFindGridArray, 
					     FastFindDeltaN, FastFindCellWidth,
					     next_level );
	next_level = next_level->NextGridThisLevel;
      }
    }
    this_grid = this_grid->NextGridThisLevel;
  }

  return;
}


HierarchyEntry *AnalysisBaseClass::ContainingGrid( HierarchyEntry *Grid, float *point){

  while( Grid != NULL ){
    if( Grid->GridData->PointInGrid( point ))
      break;
    Grid = Grid->NextGridThisLevel;
  }

  return Grid;
}

HierarchyEntry *AnalysisBaseClass::FindGrid(float *point){

  HierarchyEntry *NextGrid;
  // Find level 0 parent or exit
  HierarchyEntry *ParentGrid = ContainingGrid( TopGrid, point );

  if( ParentGrid )
    while( ParentGrid->NextGridNextLevel ){
      NextGrid = ContainingGrid( ParentGrid->NextGridNextLevel, point );
      if( NextGrid )
	ParentGrid = NextGrid;
      else
	break;
    }
  else
    my_exit(EXIT_FAILURE);
  
  return ParentGrid;
}


HierarchyEntry *AnalysisBaseClass::FastFindGrid(float *point){

  int i, i_l[MAX_DIMENSION];
  HierarchyEntry *NextGrid;
  HierarchyEntry *Grid;

  if(!FastFindGridArray)
    InitializeFastFind();

  for(i=0; i < MAX_DIMENSION; i++){
    i_l[i] = int((point[i] - DomainLeftEdge[i])/((float)(FastFindCellWidth[i])));
  }

  Grid = FastFindGridArray[i_l[0] + i_l[1]*FastFindDeltaN[0] +
				   i_l[2]*FastFindDeltaN[0]*FastFindDeltaN[1]];

  if(MaximumAnalysisLevel > 1){
    
    if( Grid ){
      while( Grid->NextGridNextLevel ){
	NextGrid = ContainingGrid( Grid->NextGridNextLevel, point );
	if( NextGrid ){
	  Grid = NextGrid;
	}else{
	  break;
	}
      }
    }else{
      my_exit(EXIT_FAILURE);
    }
  }

  return Grid;
}
