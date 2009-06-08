#ifndef ANALYSIS_BASE_CLASS_H
#define ANALYSIS_BASE_CLASS_H

/***********************************************************************
/  
/ ANALYSIS SUPERCLASS
/
************************************************************************/

class AnalysisBaseClass{
  
 public:
  AnalysisBaseClass( TopGridData *metadata,
		HierarchyEntry *topgrid,
		int level = -1,
		FLOAT *Left = NULL,
		FLOAT *Right = NULL);

  ~AnalysisBaseClass();

  void PrintParameters();
  void PrintGridInfo( grid *Grid);
  void SetAnalysisRegion( FLOAT Left[MAX_DIMENSION],
			  FLOAT Right[MAX_DIMENSION] );
  int NumberOfGridsInVolume( );
  void InitializeFastFind();
  int IndexGrids();
  void IndexGridsByProcessor();
  void BuildMyGridArray();
  void CountGridsPerProcessor();
  HierarchyEntry *FindGrid(float *point);
  HierarchyEntry *FastFindGrid(float *point);

 protected:

  const float G;
  const float MPC_CM;
  const float MSOLAR_G;

  TopGridData *MetaData;
  HierarchyEntry *TopGrid;
  int MaximumAnalysisLevel;
  FLOAT AnalysisRegionLeftEdge[MAX_DIMENSION];
  FLOAT AnalysisRegionRightEdge[MAX_DIMENSION];
  FLOAT DomainDelta[MAX_DIMENSION];

  void PrintTrueFalse(FILE *buffer, bool testvalue,
		      const char *fmt="%s\n");

  void LoopNumberOfGrids( HierarchyEntry *Grid,
				     int current_level);
  void LoopIndexGrids( HierarchyEntry *Grid,
		       int current_level,
		       int *current_index );
  void LoopIndexGridsByProcessor(HierarchyEntry *Grid,
				 int current_level,
				 int *current_indexes);
  void LoopCountGridsPerProcessor(HierarchyEntry *Grid,
				  int current_level);
  void LoopBuildMyGridArray(HierarchyEntry *Grid,
			    int current_level,
			    int *current_index);

  void FlagGridCells(HierarchyEntry *Grid);

  float CurrentRedshift;

  void HDF5CreateFile( char *name );
  hid_t HDF5OpenFile( char *name );
  void HDF5CloseFile( hid_t file_id );
  hid_t HDF5CreateGroup( hid_t file_id, char *name );
  hid_t HDF5OpenGroup( hid_t file_id, char *name );
  void HDF5CloseGroup( hid_t group_id );
  void HDF5MakeDataset( hid_t loc_id, 
			char *dset_name,
			int rank, hsize_t dims[], 
			float *data,
			char *units = NULL);
  void HDF5MakeDataset( hid_t group_id, 
			char *dset_name,
			int rank, hsize_t dims[], 
			int *data,
			char *units = NULL);

  void HDF5ReadDataset( hid_t group_id,
			char *dset_name,
			float *data );
 
  void HDF5WriteStringAttr(hid_t dset_id, char *Alabel, char *String);

  HierarchyEntry *ContainingGrid( HierarchyEntry *Grid, float *point);

  HierarchyEntry **FastFindGridArray;
  int *FastFindDeltaN;
  float *FastFindCellWidth;

  HierarchyEntry **MyGridArray;

  int NumberOfGrids;
  Eint32 *GridsPerProcessor;

 public:
  inline float ReturnCurrentRedshift(void){ return CurrentRedshift;};
  inline int ReturnMaximumAnalysisLevel(void){ return MaximumAnalysisLevel;};
  inline void SetMaximumAnalysisLevel(int max_level)
  { MaximumAnalysisLevel = max_level;};
  inline Eint32 ReturnGridsPerProcessor(int proc_id)
  { return GridsPerProcessor[proc_id];};

};

#endif
