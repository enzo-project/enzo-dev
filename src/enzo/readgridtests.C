#include "enzo_unit_tests.h"

extern std::map<HierarchyEntry *, int> OriginalGridID;
#ifdef USE_HDF5_GROUPS
extern std::map<HierarchyEntry *, int> OriginalTaskID;
#endif

// Start the overall test suite
START_TESTS()
  CommunicationInitialize(&argc, &argv);
  // Initialize Communications   

#ifdef TASKMAP
  START_TEST("no testing - taskmap on")
  ASSERT(1);
  END_TEST()
#else

  TopGridData MetaData;
  HierarchyEntry TopGrid;
  ExternalBoundary Exterior;

// the data (parameter files, etc) for this tests are
// in the lcatests tests/Unit/ReadGridTests directory
#ifndef USE_HDF5_GROUPS
  char ParameterFile[MAX_LINE_LENGTH]          = "collapse_unpacked";
#else
  char ParameterFile[MAX_LINE_LENGTH]          = "collapse_packed";
#endif

  SetDefaultGlobalValues(MetaData);

  LoadGridDataAtStart = FALSE;

// A new group of tests, with an identifier
START_TEST("read data")
#ifndef USE_HDF5_GROUPS
  ASSERT(ReadAllData(ParameterFile, &TopGrid, MetaData, &Exterior));
#else
  ASSERT(Group_ReadAllData(ParameterFile, &TopGrid, MetaData, &Exterior));
#endif
  ASSERT_EQUALS_INT(MaximumRefinementLevel, 9);
END_TEST()

START_TEST("grid data not read in")
  EnzoArrayFloat *dens = TopGrid.GridData->CreateFieldArrayFloat(Density);
  ASSERT_EQ(dens->Array, NULL);
  delete dens;
END_TEST()

START_TEST("original grid id")
  ASSERT_EQ(OriginalGridID[&TopGrid], 1);
  ASSERT_EQ(OriginalGridID[TopGrid.NextGridNextLevel], 2);
END_TEST()

#ifdef USE_HDF5_GROUPS
START_TEST("original task id")
  ASSERT_EQ(OriginalTaskID[&TopGrid], 0);
  ASSERT_EQ(OriginalTaskID[TopGrid.NextGridNextLevel], 1);
END_TEST()
#endif

START_TEST("read in grid data")
  // this should be grid 5
  HierarchyEntry *grid_to_read = TopGrid.NextGridNextLevel->
  NextGridNextLevel->NextGridNextLevel;
#ifndef USE_HDF5_GROUPS
  int read_result =  grid_to_read->GridData->ReadGrid(NULL, 5, FALSE, TRUE);
#else
  hid_t file_id = H5Fopen("collapse_packed.cpu0001", H5F_ACC_RDWR, H5P_DEFAULT);
int read_result =  grid_to_read->GridData->Group_ReadGrid(NULL, 5, file_id, FALSE, TRUE);
hid_t h5_status = H5Fclose(file_id);
#endif

  EnzoArrayFloat *dens = grid_to_read->GridData->CreateFieldArrayFloat(Density);
if(MyProcessorNumber == grid_to_read->GridData->ReturnProcessorNumber()){
  ASSERT_NOT_EQUALS(dens->Array, NULL);
 }else{
  ASSERT_EQ(NULL, NULL);
 }
delete dens;
END_TEST()

#endif // TASMAP

  CommunicationFinalize();
// End the overall test suite
END_TESTS()

