/***********************************************************************
/
/  EXTERNAL BOUNDARY CLASS
/
/  written by: Greg Bryan
/  date:       November, 1994
/  modified1:  
/
/  PURPOSE:
/
************************************************************************/

#ifndef EXTERNAL_BOUNDARY_DEFINED__
#define EXTERNAL_BOUNDARY_DEFINED__


class ExternalBoundary
{
 private:
  int  BoundaryRank;                      // This is the rank and dimension 
  int  BoundaryDimension[MAX_DIMENSION];  //  of the grid to which the boundary
                                          //  values apply
  int  NumberOfBaryonFields;              // Number of boundary fields
  int  BoundaryFieldType[MAX_NUMBER_OF_BARYON_FIELDS];           
                                          // Field types

  boundary_type *BoundaryType[MAX_NUMBER_OF_BARYON_FIELDS][MAX_DIMENSION][2];
                                          /* Type of boundary value:
                                              1 - reflecting
					      2 - outflow
					      3 - inflow
					      4 - periodic   */

  boundary_type ParticleBoundaryType;

  int  MagneticBoundaryDims[3][MAX_DIMENSION];  //  of the grid to which the boundary
                                          //  values apply 
  boundary_type MagneticBoundaryType[3][MAX_DIMENSION][2];
  float *MagneticBoundaryValue[3][3][2];

  float *BoundaryValue[MAX_NUMBER_OF_BARYON_FIELDS][MAX_DIMENSION][2];  
					  // boundary values for inflow (3)

  friend class grid;						    

 public:
//
// Constructor (set rank and dims from TopGrid)
//
  ExternalBoundary();
//
// Prepare External Boundary (set dims, etc.) based on TopGrid in argument.
//
  int Prepare(class grid *TopGrid);
//
// Checks if the External Boundary has been prepared
//
  int AmIPrepared() {return (BoundaryRank > 0) ? TRUE : FALSE;};
//
// Set one face of external boundaries to a constant value 
//  (Note: this is not suitable for setting inflow conditions as
//         BoundaryValue is not set).
//  Returns SUCCESS or FAIL.
//
  int InitializeExternalBoundaryFace(int Dimension, 
			      boundary_type LeftBoundaryType,
			      boundary_type RightBoundaryType, 
			      float LeftBoundaryValue[],
			      float RightBoundaryValue[]);
//
//  Initialize particle boundary conditions
//
  void InitializeExternalBoundaryParticles(boundary_type ParticleBoundary)
    {ParticleBoundaryType = ParticleBoundary;};
//
// Read an external boundary
//
  int ReadExternalBoundary(FILE *fptr, int ReadText = TRUE, int ReadData = TRUE);
//
// Read an external boundary for hdf4
//
  int ReadExternalBoundaryHDF4(FILE *fptr);
//
// Write an external boundary
//
  int WriteExternalBoundary(FILE *fptr, char *hdfname);
//
// Given a pointer to a field and its field type, find the equivalent
//   field type in the list of boundary's and apply that boundary value/type.
//   Returns: 0 on failure
//
  int SetExternalBoundary(int FieldRank, int GridDims[], int GridOffset[],
                          int StartIndex[], int EndIndex[],
                          float *Field, int FieldType);
  int SetMagneticBoundary(int FieldRank, int GridDims[], int GridOffset[],
                          int StartIndex[], int EndIndex[],
                          float *Field, int FieldType);
//
// This routine handle the boundary conditions for particles.  The conditions
//   are assumed to be the same as the mass field.
//
  int SetExternalBoundaryParticles(int FieldRank, int NumberOfParticles,
                                   FLOAT *Position[], float *Velocity[]);
//
// Finds and returns the indexes to commonly used physical quantities.
//
  int IdentifyPhysicalQuantities(int &DensNum, int &GENum, int &Vel1Num, 
                                 int &Vel2Num, int &Vel3Num, int &TENum);
  int IdentifyPhysicalQuantities(int &DensNum, int &GENum, int &Vel1Num, 
                                 int &Vel2Num, int &Vel3Num, int &TENum,
				 int &B1Num, int&B2Num, int &B3Num, int &PhiNum);

//
/************************************************************************/
//
// WavePool test problem:
//  This routine sets up the inflow boundary conditions to model an inflowing
//   linear wave (from the left boundary).  See also WavePoolGlobalData.h.
//
  int SetWavePoolBoundary(FLOAT time);
//
// ShockPool test problem:
//  This routine sets up the inflow boundary conditions to model an inflowing
//   shock wave (from the left boundary).  See also ShockPoolGlobalData.h.
//
  int SetShockPoolBoundary(FLOAT time);
//
// DoubleMach problem:
//  This routine sets up the necessary inflow boundary conditions.
//
  int SetDoubleMachBoundary(FLOAT time, FLOAT CellLeftEdge[], 
                            FLOAT CellWidth[]);

// Wengen colliding flow with pseudo cooling has  a particular shape.
  int SetWengenCollidingFlowBoundary(FLOAT time, FLOAT CellLeftEdge[], 
                            FLOAT CellWidth[]);

/*  RandomForcing tricks. */

  int AppendForcingToBaryonFields();
  int DetachForcingFromBaryonFields();

  int AddField(int FieldType);
  int DeleteObsoleteFields(int *ObsoleteFields, 
			   int NumberOfObsoleteFields);

};


#endif
