#include "enzo_unit_tests.h"

grid *Linear3DGrid();
grid *Linear3DSubGrid();

// Start the overall test suite
START_TESTS(){

  CommunicationInitialize(&argc, &argv);

  TopGridData MetaData;
  int set_result = SetDefaultGlobalValues(MetaData);

  {
    START_TEST("get density float array"){
      
      grid *agrid = Linear3DGrid();
      
      EnzoArrayFloat *dens = agrid->CreateFieldArrayFloat(Density);
    
      Eint32 index = 7 + 8*134 + 9*134*134;
    
      ASSERT_EQUALS_INT(dens->Rank, 3);
      ASSERT_EQUALS_INT(dens->Dimension[0], 134);
      ASSERT_EQUALS_INT(dens->Dimension[1], 134);
      ASSERT_EQUALS_INT(dens->Dimension[2], 134);
    
      ASSERT_EQUALS_INT(dens->StartIndex[0], 3);
      ASSERT_EQUALS_INT(dens->EndIndex[2], 130);
    
      ASSERT_EQUALS_FLOAT(dens->Array[index], 7 + 8*1000 + 9*1000000);
    
      delete agrid;
      delete dens;
    }END_TEST()
       }

  {
    START_TEST("get density float array by name"){
  
      grid *agrid = Linear3DGrid();

      EnzoArrayFloat *dens = agrid->CreateFieldArrayFloat("Density");
    
      Eint32 index = 7 + 8*134 + 9*134*134;
     
      ASSERT_EQUALS_INT(dens->Rank, 3);
      ASSERT_EQUALS_INT(dens->Dimension[0], 134);
      ASSERT_EQUALS_INT(dens->Dimension[1], 134);
      ASSERT_EQUALS_INT(dens->Dimension[2], 134);

      ASSERT_EQUALS_INT(dens->StartIndex[0], 3);
      ASSERT_EQUALS_INT(dens->EndIndex[2], 130);

      ASSERT_EQUALS_FLOAT(dens->Array[index], 7 + 8*1000 + 9*1000000);

      delete agrid;
      delete dens;
    }END_TEST()
       }
  {
    START_TEST("get particle type array"){
    
      grid *agrid = Linear3DGrid();

      EnzoArrayInt *type = agrid->CreateFieldArrayInt(gParticleType);

      ASSERT_EQUALS_INT(type->Rank, 1);
      ASSERT_EQUALS_INT(type->Dimension[0], 2);

      ASSERT_EQUALS_INT(type->StartIndex[0], 0);
      ASSERT_EQUALS_INT(type->EndIndex[0], 1);

      ASSERT_EQUALS_INT(type->Array[0], 1);
      ASSERT_EQUALS_INT(type->Array[1], 3);

      delete agrid;
    }END_TEST()
       }
  {
    START_TEST("get particle type array by name"){
  
      grid *agrid = Linear3DGrid();

      EnzoArrayInt *type = agrid->CreateFieldArrayInt("ParticleType");

      ASSERT_EQUALS_INT(type->Rank, 1);
      ASSERT_EQUALS_INT(type->Dimension[0], 2);

      ASSERT_EQUALS_INT(type->StartIndex[0], 0);
      ASSERT_EQUALS_INT(type->EndIndex[0], 1);

      ASSERT_EQUALS_INT(type->Array[0], 1);
      ASSERT_EQUALS_INT(type->Array[1], 3);

      delete agrid;
    }END_TEST()
       }
  {
    START_TEST("check for NULL - float, int"){
  
      grid *agrid = Linear3DGrid();

      EnzoArrayFloat *ptr = agrid->CreateFieldArrayFloat(gParticleType);

      ASSERT_EQUALS(ptr, NULL);

      delete agrid;
    }END_TEST()
       {
       }
    START_TEST("check for NULL - int, float"){
  
      grid *agrid = Linear3DGrid();

      EnzoArrayInt *ptr = agrid->CreateFieldArrayInt(Density);

      ASSERT_EQUALS(ptr, NULL);

      delete agrid;
      delete ptr;
    }END_TEST()
       {
       }
    START_TEST("get velocity array"){
  
      grid *agrid = Linear3DGrid();

      EnzoArrayFloat *vel = agrid->CreateFieldArrayFloat(gVelocity);
    
      Eint32 index = 7 + 8*134 + 9*134*134;
     
      ASSERT_EQUALS_INT(vel->Rank, 3);
      ASSERT_EQUALS_INT(vel->Dimension[0], 134);
      ASSERT_EQUALS_INT(vel->Dimension[1], 134);
      ASSERT_EQUALS_INT(vel->Dimension[2], 134);

      ASSERT_EQUALS_INT(vel->StartIndex[0], 3);
      ASSERT_EQUALS_INT(vel->EndIndex[2], 130);

      ASSERT_EQUALS(vel->Array, NULL);
      ASSERT_EQUALS_FLOAT(vel->Vector[0][index], 0.125);
      ASSERT_EQUALS_FLOAT(vel->Vector[1][index+3], 0.25);
      ASSERT_EQUALS_FLOAT(vel->Vector[2][index-30], 0.375);

      delete agrid;
      delete vel;
    }END_TEST()
       }
  {
    START_TEST("get velocity array by name"){
  
      grid *agrid = Linear3DGrid();

      EnzoArrayFloat *vel = agrid->CreateFieldArrayFloat("Velocity");
    
      if(vel == NULL)
	printf("uh, oh\n");

      Eint32 index = 7 + 8*134 + 9*134*134;
     
      ASSERT_EQUALS_INT(vel->Rank, 3);
      ASSERT_EQUALS_INT(vel->Dimension[0], 134);
      ASSERT_EQUALS_INT(vel->Dimension[1], 134);
      ASSERT_EQUALS_INT(vel->Dimension[2], 134);

      ASSERT_EQUALS_INT(vel->StartIndex[0], 3);
      ASSERT_EQUALS_INT(vel->EndIndex[2], 130);

      ASSERT_EQUALS(vel->Array, NULL);
      ASSERT_EQUALS_FLOAT(vel->Vector[0][index], 0.125);
      ASSERT_EQUALS_FLOAT(vel->Vector[1][index+3], 0.25);
      ASSERT_EQUALS_FLOAT(vel->Vector[2][index-30], 0.375);

      delete agrid;
      delete vel;
    }END_TEST()
       }

  {
    START_TEST("flag refined grid cells"){
  
      grid *parent_grid = Linear3DGrid();
      grid *sub_grid = Linear3DSubGrid();

      parent_grid->ClearFlaggingField();
      parent_grid->FlagRefinedCells(sub_grid);

      EnzoArrayInt *flagging_field = parent_grid->CreateFieldArrayInt("FlaggingField");

      ASSERT_NOT_EQUALS(NULL, flagging_field->Array);

      int i_x, i_y, i_z;
      i_x = 6;
      i_y = 12;
      i_z = 14;
      int index = (i_z + 3)*134*134 + (i_y + 3)*134 + (i_x + 3);

      ASSERT_EQUALS_INT(flagging_field->Array[index], 1);
      ASSERT_EQUALS_INT(flagging_field->Array[index - 1], 0);

      i_x = 69;
      i_y = 75;
      i_z = 77;
      index = (i_z + 3)*134*134 + (i_y + 3)*134 + (i_x + 3);

      ASSERT_EQUALS_INT(flagging_field->Array[index], 1);
      ASSERT_EQUALS_INT(flagging_field->Array[index + 1], 0);

      delete parent_grid;
      delete sub_grid;
      delete flagging_field;
    }END_TEST()
       }

  {
    START_TEST("point in grid"){
  
      grid *agrid = Linear3DGrid();
      float point[] = {0., 0.125, 0.};

      ASSERT_EQUALS(agrid->PointInGrid(point), TRUE);

      point[0] = -0.125;
      ASSERT_EQUALS(agrid->PointInGrid(point), FALSE);

      FLOAT pointF[] = {0., 0.125, 0.};

      ASSERT_EQUALS(agrid->PointInGrid(pointF), TRUE);

      pointF[0] = -0.125;
      ASSERT_EQUALS(agrid->PointInGrid(pointF), FALSE);

      delete agrid;
    }END_TEST()
       }

  {
    START_TEST("grid is in volume"){
  
      grid *agrid = Linear3DGrid();
      float left[] = {0., 0.125, 0.};
      float right[] = {0.2, 0.25, 0.1};

      ASSERT_EQUALS(agrid->IsInVolume(left, right), TRUE);

      left[0] = -0.125;
      left[1] = -0.25;
      left[2] = -0.5;

      right[0] = -0.025;
      right[1] = -0.15;
      right[2] = -0.4;

      ASSERT_EQUALS(agrid->IsInVolume(left, right), FALSE);

      FLOAT leftF[] = {0., 0.125, 0.};
      FLOAT rightF[] = {0.2, 0.125, 0.1};

      ASSERT_EQUALS(agrid->IsInVolume(leftF, rightF), TRUE);

      leftF[0] = -0.125;
      leftF[1] = -0.25;
      leftF[2] = -0.5;

      rightF[0] = -0.025;
      rightF[1] = -0.15;
      rightF[2] = -0.4;

      ASSERT_EQUALS(agrid->IsInVolume(leftF, rightF), FALSE);

      delete agrid;
    }END_TEST()
       }

  {
    START_TEST("flag grid array"){
  
      grid *subgrid = Linear3DSubGrid();

      int i, num_cells = 128*128*128;
      HierarchyEntry** ff_array = new HierarchyEntry*[num_cells];

      for(i = 0; i < num_cells; i++)
	ff_array[i] = NULL;

      int ff_deltan[] = {128, 128, 128};
      float ff_cellwidth[] = {1./128,  1./128,  1./128};
      
      HierarchyEntry he;

      subgrid->FlagGridArray( &ff_array, ff_deltan, ff_cellwidth,
			      &he );

      ASSERT_EQUALS(ff_array[0], NULL);
      ASSERT_EQUALS(ff_array[num_cells - 1], NULL);

      // no ghost zones in the flagging thingy
      // check the corners to make sure the subgrid was flagged
      int i_x, i_y, i_z;
      i_x = 6;
      i_y = 12;
      i_z = 14;
      int index = (i_z)*128*128 + (i_y)*128 + (i_x);

      ASSERT_EQUALS_INT(ff_array[index], &he);
      ASSERT_EQUALS_INT(ff_array[index - 1], NULL);

      i_x = 69;
      i_y = 75;
      i_z = 77;
      index = (i_z)*128*128 + (i_y)*128 + (i_x);

      ASSERT_EQUALS_INT(ff_array[index], &he);
      ASSERT_EQUALS_INT(ff_array[index + 1], NULL);

      delete [] ff_array;
      delete subgrid;
    }END_TEST()
       }

  // End the overall test suite
  CommunicationFinalize();

 }END_TESTS()


grid *Linear3DGrid(){
  // Create a new 3D grid
  float dens = M_PI, total_energy = 0.5, internal_energy = 0.0;
  float vel[3];
  float BField[3];
  int dims[3];
  FLOAT left[3], right[3];

  grid *lineargrid = new grid;
  int i, j, k, rank = 3;
  int index;
   
  for (i = 0; i < rank; i++) {
    dims[i] = 134;
    left[i] = 0.0;
    right[i] = 1.0;
    vel[i] = (i+1) * 0.125;
    BField[i] = 0.1;
  }
    
  NumberOfParticleAttributes = 0;
  lineargrid->PrepareGrid(3, dims, 
			  left, right, 2);

  int result = lineargrid->InitializeUniformGrid
    (dens, total_energy, internal_energy, vel, BField);
  assert(result != FAIL);
  //  float *dens_field = lineargrid->AccessDensity();

  EnzoArrayFloat *dens_field = lineargrid->CreateFieldArrayFloat("Density");

  for (k = 3; k <= 130; k++) {
    for (j = 3; j <= 130; j++) {
      index =  k*(134)*(134) +
	j*(134) + 3;
      for (i = 3; i <= 130; i++, index++) {
	dens_field->Array[index] = (float)(i + 1000*j + 1000000*k);
      }
    }
  }

  // add particles
  float *mass = new float[2];
  int *number = new int[2];
  number[0] = 0;
  number[0] = 1;

  int *type = new int[2];
  type[0] = 1;
  type[1] = 3;

  float **null_ptr = new float*[3];
  null_ptr[0] = null_ptr[1] = null_ptr[2] = NULL;

  FLOAT **null_ptrF = new FLOAT*[3];
  null_ptrF[0] = null_ptrF[1] = null_ptrF[2] = NULL;

  lineargrid->SetParticlePointers(mass, number, type, null_ptrF, null_ptr, null_ptr);
  
  delete dens_field;
  
  return lineargrid;

 }


grid *Linear3DSubGrid(){
  // Create a new 3D grid
  float dens = M_PI, total_energy = 0.5, internal_energy = 0.0;
  float vel[3];
  int dims[3];
  FLOAT left[3], right[3];

  grid *lineargrid = new grid;
  int i, j, k, rank = 3;
  int index;
  float BField[3];
   
  for (i = 0; i < rank; i++) {
    dims[i] = 134;
    vel[i] = (i+1) * 0.125;
    BField[i] = 0.1;
  }

  // 6 cells over
  left[0] = 0.046875;
  right[0] = 0.546875;

  // 12 cells over
  left[1] = 0.09375;
  right[1] = 0.59375;

  // 14 cells up
  left[2] = 0.109375;
  right[2] = 0.609375;

  NumberOfParticleAttributes = 0;
  lineargrid->PrepareGrid(3, dims, 
			  left, right, 2);

  int result = lineargrid->InitializeUniformGrid
    (dens, total_energy, internal_energy, vel, BField);
  assert(result != FAIL);
  
  return lineargrid;

 }
