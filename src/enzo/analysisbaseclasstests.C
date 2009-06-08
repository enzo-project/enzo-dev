#include "enzo_unit_tests.h"

#include "AnalysisBaseClass.h"

// Start the overall test suite
START_TESTS()
  CommunicationInitialize(&argc, &argv);
  // Initialize Communications   

  TopGridData MetaData;
  HierarchyEntry TopGrid;
  ExternalBoundary Exterior;

  SetDefaultGlobalValues(MetaData);

  LoadGridDataAtStart = FALSE;

#ifndef USE_HDF5_GROUPS
  char ParameterFile[MAX_LINE_LENGTH]          = "collapse_unpacked";
  ReadAllData(ParameterFile, &TopGrid, MetaData, &Exterior);
#else
  char ParameterFile[MAX_LINE_LENGTH]          = "collapse_packed";
Group_ReadAllData(ParameterFile, &TopGrid, MetaData, &Exterior);
#endif

  AnalysisBaseClass *analysis_obj = NULL;



START_TEST("init analysis")
  analysis_obj = new AnalysisBaseClass(&MetaData, &TopGrid);
ASSERT_EQ(analysis_obj->ReturnMaximumAnalysisLevel(), MaximumRefinementLevel);

END_TEST()

// A new group of tests, with an identifier
START_TEST("good point")
  float good_point[] = {0.25, 0.25, 0.25};
  ASSERT(analysis_obj->FindGrid(good_point));
END_TEST()

START_TEST("good point - subgrid")
  float good_point[] = {0.6, 0.6, 0.6};
  ASSERT_EQUALS(analysis_obj->FindGrid(good_point), TopGrid.NextGridNextLevel);
END_TEST()

START_TEST("fast find grid array")
  float good_point[] = {0.6, 0.6, 0.6};
   HierarchyEntry *found_grid;
   analysis_obj->InitializeFastFind();
found_grid = analysis_obj->FastFindGrid(good_point);
   ASSERT_EQ( TopGrid.NextGridNextLevel, found_grid);
END_TEST()

START_TEST("level 0 fast find grid array")
  float good_point[] = {0.6, 0.6, 0.6};
   HierarchyEntry *found_grid;
   analysis_obj->SetMaximumAnalysisLevel(0);
   analysis_obj->InitializeFastFind();
   found_grid = analysis_obj->FastFindGrid(good_point);
   analysis_obj->SetMaximumAnalysisLevel(MaximumRefinementLevel);
   ASSERT_EQ( &TopGrid, found_grid);
END_TEST()

START_TEST("count grids per processor")
   //   analysis_obj->PrintParameters();
   analysis_obj->CountGridsPerProcessor();
   ASSERT_GT( analysis_obj->ReturnGridsPerProcessor(0), 0);
END_TEST()

  CommunicationFinalize();
// End the overall test suite
END_TESTS()
