#include "enzo_unit_tests.h"

#define GRID_DIMS 10
grid *Uniform3DGrid();
grid *Linear3DGrid();

// Start the overall test suite
START_TESTS(){

  CommunicationInitialize(&argc, &argv);
  if(NumberOfProcessors != 2){
    if(MyProcessorNumber == ROOT_PROCESSOR){
      fprintf(stderr, "This test is intended for 2 MPI tasks only!\n");
    }
    my_exit(EXIT_FAILURE);
  }

  TopGridData MetaData;
  int set_result = SetDefaultGlobalValues(MetaData);

  {
    START_TEST("get uniform grid"){
      
      grid *agrid = Uniform3DGrid();
      EnzoArrayFloat *dens = agrid->CreateFieldArrayFloat(Density);

      if(MyProcessorNumber == 0){
	// Ensure that the density field is there
	ASSERT_EQUALS_FLOAT(dens->Array[0], 1.0);
      }else{
	// Ensure that the density field is not there
	ASSERT_EQUALS(dens->Array, NULL);
      }

      delete dens;    
      delete agrid;

    }END_TEST()

  {
    START_TEST("communication move grid"){
      
      grid *agrid = Uniform3DGrid();
      EnzoArrayFloat *dens = agrid->CreateFieldArrayFloat(Density);

      if(MyProcessorNumber == 0){
	// density field size is 10**3
	int size = 1000;
	int bad_count = 0;
	for(int i = 0; i < size; i++){
	  if(dens->Array[i] != 1.0){
	    bad_count++;
	    // fprintf(stderr, "Bad value %"FSYM" at %"ISYM"\n", dens->Array[i], i);
	  }
	}
	ASSERT_EQUALS_INT(bad_count, 0);
      }else{
	// Ensure that the density field is not there
	ASSERT_EQUALS(dens->Array, NULL);
      }

      delete dens;

      // Move the grid
      agrid->CommunicationMoveGrid(1);
      ASSERT_EQUALS_INT(agrid->ReturnProcessorNumber(), 1);

      dens = agrid->CreateFieldArrayFloat(Density);
      EnzoArrayFloat *vel3 = agrid->CreateFieldArrayFloat(Velocity3);
      if(MyProcessorNumber == 1){
	// Ensure that the density field is there
	ASSERT_EQUALS_FLOAT(dens->Array[0], 1.0);
	// density field size is (3 + 4 + 3)**3
	int size = 1000;
	int dens_bad_count = 0;
	for(int i = 0; i < size; i++){
	  if(dens->Array[i] != 1.0){
	    dens_bad_count++;
	    // fprintf(stderr, "Bad value %"FSYM" at %"ISYM"\n", dens->Array[i], i);
	  }
	}
	ASSERT_EQUALS_INT(dens_bad_count, 0);

	int vel3_bad_count = 0;
	for(int i = 0; i < size; i++){
	  if(vel3->Array[i] != 6.0){
	    vel3_bad_count++;
	    // fprintf(stderr, "Bad value %"FSYM" at %"ISYM"\n", vel3->Array[i], i);
	  }
	}
	ASSERT_EQUALS_INT(vel3_bad_count, 0);

      }else{
	// Ensure that the density field is not there
	ASSERT_EQUALS(dens->Array, NULL);
      }
      
      delete vel3;
      delete dens;    
      delete agrid;

    }END_TEST()
       }

  {
    START_TEST("communication move linear grid"){
      
      grid *agrid = Linear3DGrid();

      EnzoArrayFloat *dens = agrid->CreateFieldArrayFloat(Density);

      if(MyProcessorNumber == 0){
	ASSERT_EQUALS_FLOAT(dens->Array[0], 100.0);
      }else{
	// Ensure that the density field is not there
	ASSERT_EQUALS(dens->Array, NULL);
      }

      delete dens;

      agrid->CommunicationMoveGrid(1);
      ASSERT_EQUALS_INT(agrid->ReturnProcessorNumber(), 1);

      dens = agrid->CreateFieldArrayFloat(Density);

      if(MyProcessorNumber == 1){
	int size = 1000;
	int dens_bad_count = 0;
	for(int i = 0; i < size; i++){
	  if(dens->Array[i] != i + 100.0){
	    dens_bad_count++;
	    // fprintf(stderr, "Bad value %"FSYM" at %"ISYM"\n", dens->Array[i], i);
	  }
	}
	ASSERT_EQUALS_INT(dens_bad_count, 0);

      }else{
	// Ensure that the density field is not there
	ASSERT_EQUALS(dens->Array, NULL);
      }

      delete dens;    
      delete agrid;

    }END_TEST()
       }

  }

  // End the overall test suite
  CommunicationFinalize();

 }END_TESTS()


grid *Uniform3DGrid(){
  // Create a new 3D grid
  float dens = 1.0, total_energy = 0.5, internal_energy = 0.0;
  float vel[3];
  int dims[3];
  FLOAT left[3], right[3];

  grid *uniformgrid = new grid;
  uniformgrid->SetProcessorNumber(0);

  int i, j, k, rank = 3;
  int index;
   
  for (i = 0; i < rank; i++) {
    dims[i] = GRID_DIMS;
    left[i] = 0.0;
    right[i] = 1.0;
    vel[i] = (i+1) * 2.0;
  }
    
  NumberOfParticleAttributes = 0;
  uniformgrid->PrepareGrid(3, dims, 
			  left, right, 0);

  int result = uniformgrid->InitializeUniformGrid(dens, total_energy, internal_energy, vel);
  assert(result != FAIL);
  
  return uniformgrid;

 }

grid *Linear3DGrid(){
  // Create a new 3D grid
  float dens = 1.0, total_energy = 0.5, internal_energy = 0.0;
  float vel[3];
  int dims[3];
  FLOAT left[3], right[3];

  grid *lineargrid = new grid;
  lineargrid->SetProcessorNumber(0);

  int i, j, k, rank = 3;
  int index;
   
  for (i = 0; i < rank; i++) {
    dims[i] = GRID_DIMS;
    left[i] = 0.0;
    right[i] = 1.0;
    vel[i] = (i+1) * 2.0;
  }
    
  NumberOfParticleAttributes = 0;
  lineargrid->PrepareGrid(3, dims, 
			  left, right, 0);

  int result = lineargrid->InitializeUniformGrid(dens, total_energy, internal_energy, vel);
  assert(result != FAIL);

  if(MyProcessorNumber == 0){
    EnzoArrayFloat *dens_array = lineargrid->CreateFieldArrayFloat(Density);
    
    for( i = 0; i < GRID_DIMS*GRID_DIMS*GRID_DIMS; i++ )
      dens_array->Array[i] = 100.0 + i;
  }
  
  return lineargrid;

 }

