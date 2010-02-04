/***********************************************************************
/
/  CREATE (GET) A FIELD ARRAY
/
/  written by: Rick Wagner
/  date:       February, 2008
/  modified1:
/
/  PURPOSE: To simplify access to grid data.
/
/  grid::CreateFieldArray"Type"(field_type field), et. al.
/
/  This returns an array object defined in EnzoArray.h, whose
/  Array member is either a float *, int *, or FLOAT *, depending
/  on the field being accessed. There is also a Vector member,
/  which is used for Velocity, ParticlePosition, Accel, etc.
/
/  The argument field maps to the fields defined in typedefs.h, 
/  which includes some additional pseudo-fields to get to other
/  data (like the particles). This is also spelled out in 
/  the field_map_list below.
/
/  To extend the types of data returned, do the following:
/
/  For a baryon field, make sure it's added to the list in
/  typedefs.h, and the FieldUndefined is incremented.
/
/  For other fields, or pseudo-fields:
/
/   1. Add the field or pseudo-field in typedefs.h
/
/   2. Increment NUM_FIELDS
/
/   3. Add the field, its name, and its type to field_map_list.
/      0 = float, 1 = FLOAT, 2 = int
/      There's no order dependence here.
/
/   4. Add a new case to the switch statement in
/      the approiate CreateFieldArray"Type" function.
/      I.e., int fields go in CreateFieldArrayInt,
/      float fields go in CreateFieldArrayFloat.
/      Look at the existing cases for how to do this.
/
/      Again, this doesn't need to be done for baryon fields,
/      since they're easy to find. 
/
/      The CreateFieldArray"Type"(field_type field) functions
/      are the workers. All the others are just for
/      convenience.
/
************************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"

#include "Grid.h"

int FindField(int field, int farray[], int numfields);

struct field_map{
  field_type field;
  char *field_name;
  // 0 for float
  // 1 for FLOAT
  // 2 for int
  int data_type;
};

// This looks psychotic, but it's probably the 
// best list of possible grid fields and their
// data types. It may even be better than
// Grid.h
const field_map field_map_list[] = {
  {Density, "Density", 0},
  {TotalEnergy,"TotalEnergy", 0},
  {InternalEnergy,"InternalEnergy", 0},
  {Pressure,"Pressure", 0},
  {Velocity1,"Velocity1", 0},
  {Velocity2,"Velocity2", 0},
  {Velocity3,"Velocity3", 0},
  {ElectronDensity,"ElectronDensity", 0},
  {HIDensity,"HIDensity", 0},
  {HIIDensity,"HIIDensity", 0},
  {HeIDensity,"HeIDensity", 0},
  {HeIIDensity,"HeIIDensity", 0},
  {HeIIIDensity,"HeIIIDensity", 0},
  {HMDensity,"HMDensity", 0},
  {H2IDensity,"H2IDensity", 0},
  {H2IIDensity,"H2IIDensity", 0},
  {DIDensity,"DIDensity", 0},
  {DIIDensity,"DIIDensity", 0},
  {HDIDensity,"HDIDensity", 0},
  {Metallicity,"Metallicity", 0},
  {ExtraType0,"ExtraType0", 0},
  {ExtraType1,"ExtraType1", 0},
  {GravPotential,"GravPotential", 0},
  {Acceleration0,"Acceleration0", 0},
  {Acceleration1,"Acceleration1", 0},
  {Acceleration2,"Acceleration2", 0},
  {gParticlePosition,"ParticlePosition", 1},
  {gParticleVelocity,"ParticleVelocity", 0},
  {gParticleMass,"ParticleMass", 0},
  {gParticleAcceleration,"ParticleAcceleration", 0},
  {gParticleNumber,"ParticleNumber", 2},
  {gParticleType,"ParticleType", 2},
  {gParticleAttribute,"ParticleAttribute", 0},
  {gPotentialField,"PotentialField", 0},
  {gAccelerationField,"AccelerationField", 0},
  {gGravitatingMassField,"GravitatingMassField", 0},
  {gFlaggingField,"FlaggingField", 3},
  {gVelocity, "Velocity", 0}};

field_type get_field_id(char *field_name){
  
  // loops over list above
  // if the field isn't found, returns -1
  int i;
  field_type field_id = -1;
  
  for(i = 0; i < FieldUndefined; i++){
    if(strcmp(field_map_list[i].field_name, field_name) == 0){
      field_id = field_map_list[i].field;
      break;
    }
  }

  // not found gets a null pointer
  return field_id;
}


EnzoArray<FLOAT> *grid::CreateFieldArrayFLOAT(char *field_name){

  field_type field_id = get_field_id(field_name);

  return this->CreateFieldArrayFLOAT(field_id);
}

EnzoArray<float> *grid::CreateFieldArrayFloat(char *field_name){

  field_type field_id = get_field_id(field_name);

  return this->CreateFieldArrayFloat(field_id);
}

EnzoArray<int> *grid::CreateFieldArrayInt(char *field_name){

  field_type field_id = get_field_id(field_name);

  return this->CreateFieldArrayInt(field_id);
}

EnzoArray<bool> *grid::CreateFieldArrayBool(char *field_name){

  field_type field_id = get_field_id(field_name);

  return this->CreateFieldArrayBool(field_id);
}

EnzoArray<float> *grid::CreateFieldArrayFloat(field_type field){

  int i, dims[MAX_DIMENSION], sindex[MAX_DIMENSION], eindex[MAX_DIMENSION];
  EnzoArray<float> *array = NULL;
  FLOAT cell_width[] = {0, 0, 0};
  field_type field_index;
  
  for(i = 0; i < this->GridRank; i++){
    cell_width[i] = this->CellWidth[i][0];
  }

  // default case, accessing a baryon field
  // Acceleration2 is the highest of the baryon fields
  if( (0 <= field) && (field <= Acceleration2) ){
    field_index = FindField(field, 
			    this->FieldType, 
			    this->NumberOfBaryonFields);  
    
    if(field_index != -1){
      
      array = new EnzoArray<float>(this->GridRank,
			     this->GridDimension,
			     this->GridStartIndex,
			     this->GridEndIndex,
			     cell_width);
      
      array->Array = BaryonField[field_index];
    }
  }else{
    // accessing one of the pseudo-fields
    switch (field){      
    case gParticleVelocity:
      if(this->NumberOfParticles > 0){
  	dims[0] = this->GridRank;
  	dims[1] = this->NumberOfParticles;
  	sindex[0] = 0;
  	sindex[1] = 0;
  	eindex[0] = this->GridRank - 1;
  	eindex[1] = this->NumberOfParticles - 1;
	
  	array = new EnzoArray<float>(2, dims, sindex, eindex);
  	for(i = 0; i < this->GridRank; i++){
  	  array->Vector[i] = this->ParticleVelocity[i];
  	}
      }
      break;

    case gParticleMass:
      if(this->NumberOfParticles > 0){
  	dims[0] = this->NumberOfParticles;
  	sindex[0] = 0;
  	eindex[0] = this->NumberOfParticles - 1;
	
  	array = new EnzoArray<float>(1, dims, sindex, eindex);
  	array->Array = this->ParticleMass;

      }
      break;

    case gParticleAcceleration:
      if(this->NumberOfParticles > 0){
  	dims[0] = this->GridRank;
  	dims[1] = this->NumberOfParticles;
  	sindex[0] = 0;
  	sindex[1] = 0;
  	eindex[0] = this->GridRank - 1;
  	eindex[1] = this->NumberOfParticles - 1;

  	array = new EnzoArray<float>(2, dims, sindex, eindex);
  	for(i = 0; i < this->GridRank; i++){
  	  array->Vector[i] = this->ParticleAcceleration[i];
  	}
      }
      break;
      
    case gParticleAttribute:
      if(this->NumberOfParticles > 0){
  	dims[0] = NumberOfParticleAttributes;
  	dims[1] = this->NumberOfParticles;
  	sindex[0] = 0;
  	eindex[0] = NumberOfParticleAttributes - 1;
  	sindex[0] = 0;
  	eindex[0] = this->NumberOfParticles - 1;
	
  	array = new EnzoArray<float>(1, dims, sindex, eindex);
  	for(i = 0; i < NumberOfParticleAttributes; i++){
  	  array->Vector[i] = this->ParticleAttribute[i];
  	}

      }
      break;
      
    case gPotentialField:
      if(this->PotentialField){
  	for(i = 0; i < this->GridRank; i++){
  	  cell_width[i] = GravitatingMassFieldCellSize;
   	  sindex[i] = this->GridStartIndex[i] + GRAVITY_BUFFER_SIZE;
   	  eindex[i] = this->GridEndIndex[i] + GRAVITY_BUFFER_SIZE;
  	}
  	array = new EnzoArray<float>(this->GridRank,
			       this->GravitatingMassFieldDimension,
			       sindex, eindex,
			       cell_width);
	
  	array->Array = this->PotentialField;
      }
      break;
      
    case gAccelerationField:
      if(this->AccelerationField[0]){
	array = new EnzoArray<float>(this->GridRank,
			       this->GridDimension,
			       this->GridStartIndex,
			       this->GridEndIndex,
			       cell_width);
	
  	for(i = 0; i < this->GridRank; i++){
  	  array->Vector[i] = this->AccelerationField[i];
  	}
      }
      break;
      
    case gGravitatingMassField:
      if(this->GravitatingMassField){
        for(i = 0; i < this->GridRank; i++){
	  cell_width[i] = GravitatingMassFieldCellSize;
	  sindex[i] = this->GridStartIndex[i] + GRAVITY_BUFFER_SIZE;
	  eindex[i] = this->GridEndIndex[i] + GRAVITY_BUFFER_SIZE;
        }
  	array = new EnzoArray<float>(this->GridRank,
			       this->GravitatingMassFieldDimension,
			       sindex, eindex,
			       cell_width);
	
  	array->Array = this->GravitatingMassField;
      }
      break;
      
    case gVelocity:
      // This breaks the pattern. Here, the rank, dims,
      // and indexes refer to the individual arrays of the
      // Vector
      array = new EnzoArray<float>(this->GridRank,
			     this->GridDimension,
			     this->GridStartIndex,
			     this->GridEndIndex,
			     cell_width);
      
      field_index = FindField(Velocity1, 
  			      this->FieldType, 
  			      this->NumberOfBaryonFields);  
      
      array->Vector[0] = this->BaryonField[field_index];
      if(this->GridRank > 1){
  	field_index = FindField(Velocity2, 
  				this->FieldType, 
  				this->NumberOfBaryonFields);  
	
  	array->Vector[1] = this->BaryonField[field_index];
  	if(this->GridRank == 3){
  	  field_index = FindField(Velocity3, 
  				  this->FieldType, 
  				  this->NumberOfBaryonFields);  
	  
  	  array->Vector[2] = this->BaryonField[field_index];	    
  	}
      }
      break;
    }
  }
  
  return array;
}

EnzoArray<FLOAT> *grid::CreateFieldArrayFLOAT(field_type field){
  
  int i, dims[MAX_DIMENSION], sindex[MAX_DIMENSION], eindex[MAX_DIMENSION];
  EnzoArray<FLOAT> *array = NULL;
  //  float cell_width[] = {0, 0, 0};
  field_type field_index;
  
  //  for(i = 0; i < this->GridRank; i++){
  //    cell_width[i] = this->CellWidth[i][0];
  //  }
  
  switch (field){
  case gParticlePosition:
    if(this->NumberOfParticles > 0){
      dims[0] = this->GridRank;
      dims[1] = this->NumberOfParticles;
      sindex[0] = 0;
      sindex[1] = 0;
      eindex[0] = this->GridRank - 1;
      eindex[1] = this->NumberOfParticles - 1;
      
      array = new EnzoArray<FLOAT>(2, dims, sindex, eindex);
      for(i = 0; i < this->GridRank; i++){
	array->Vector[i] = this->ParticlePosition[i];
      }
    }
    break;
  }
  
  return array;
}

EnzoArray<int> *grid::CreateFieldArrayInt(field_type field){

  int i, dims[MAX_DIMENSION], sindex[MAX_DIMENSION], eindex[MAX_DIMENSION];
  EnzoArray<int> *array = NULL;
  FLOAT cell_width[] = {0, 0, 0};
  field_type field_index;
  
  for(i = 0; i < this->GridRank; i++){
    cell_width[i] = this->CellWidth[i][0];
  }

  switch (field){
    
  case gParticleType:
    if(this->NumberOfParticles > 0){
      dims[0] = this->NumberOfParticles;
      sindex[0] = 0;
      eindex[0] = this->NumberOfParticles - 1;
      
      array = new EnzoArray<int>(1, dims, sindex, eindex);
      array->Array = this->ParticleType;      
    }
    break;

  case gFlaggingField:
    if(this->FlaggingField){
      array = new EnzoArray<int>(this->GridRank,
			   this->GridDimension,
			   this->GridStartIndex,
			   this->GridEndIndex,
			   cell_width);

      array->Array = this->FlaggingField;
    }
    break;

    
  }
  
  return array;
}

EnzoArray<PINT> *grid::CreateFieldArrayPINT(field_type field){

  int i, dims[MAX_DIMENSION], sindex[MAX_DIMENSION], eindex[MAX_DIMENSION];
  EnzoArray<PINT> *array = NULL;
  FLOAT cell_width[] = {0, 0, 0};
  field_type field_index;
  
  for(i = 0; i < this->GridRank; i++){
    cell_width[i] = this->CellWidth[i][0];
  }

  switch (field){
    
  case gParticleNumber:
    if(this->NumberOfParticles > 0){
      dims[0] = this->NumberOfParticles;
      sindex[0] = 0;
      eindex[0] = this->NumberOfParticles - 1;
      
      array = new EnzoArray<PINT>(1, dims, sindex, eindex);
      array->Array = this->ParticleNumber;
    }
    break;
    
  }
  
  return array;
}

EnzoArray<bool> *grid::CreateFieldArrayBool(field_type field){

  int i, dims[MAX_DIMENSION], sindex[MAX_DIMENSION], eindex[MAX_DIMENSION];
  EnzoArray<bool> *array = NULL;
  FLOAT cell_width[] = {0, 0, 0};
  field_type field_index;
  
  for(i = 0; i < this->GridRank; i++){
    cell_width[i] = this->CellWidth[i][0];
  }

  // No booleans yet.
  
  return array;
}
