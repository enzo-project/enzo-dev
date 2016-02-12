/***********************************************************************
/
/ GRID CLASS (IDENTIFY THE CHEMICAL TRACER FIELDS DEPENDING ON MODEL)
/
/ written by: Andrew Emerick
/ date      : February, 2016
/ modified1 :
/
/ PURPOSE:
/
/ NOTE:
/
/***********************************************************************/

#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "fortran.def" // really?


int FindField(int f, int farray[], int n);

int grid::IdentifyChemicalTracerSpeciesFields(int  &CINum, int  &NINum, int  &OINum,
                                              int &MgINum, int &SiINum, int &FeINum,
                                              int  &YINum, int &BaINum, int &LaINum,
                                              int &EuINum)
{

  CINum = NINum = OINum = MgINum = SiINum = FeINum = YINum = BaINum = LaINum = EuINum = 0;

  //
  if(MULTIMETALS_METHOD(MULTIMETALS_ALPHA)){
     CINum = FindField( CIDensity, FieldType, NumberOfBaryonFields);
     NINum = FindField( NIDensity, FieldType, NumberOfBaryonFields);
     OINum = FindField( OIDensity, FieldType, NumberOfBaryonFields);
    MgINum = FindField(MgIDensity, FieldType, NumberOfBaryonFields);
    SiINum = FindField(SiIDensity, FieldType, NumberOfBaryonFields);
    FeINum = FindField(FeIDensity, FieldType, NumberOfBaryonFields);
  }

  if(MULTIMETALS_METHOD(MULTIMETALS_SPROCESS)){
     YINum = FindField( YIDensity, FieldType, NumberOfBaryonFields);
    BaINum = FindField(BaIDensity, FieldType, NumberOfBaryonFields);
    LaINum = FindField(LaIDensity, FieldType, NumberOfBaryonFields);
  }

  if(MULTIMETALS_METHOD(MULTIMETALS_RPROCESS)){
    EuINum = FindField(EuIDensity, FieldType, NumberOfBaryonFields);
  }

  if(( CINum < 0) || ( NINum < 0) || ( OINum < 0) || (MgINum < 0) || (SiINum < 0) ||
     (FeINum < 0) || ( YINum < 0) || (BaINum < 0) || (LaINum < 0) || (EuINum < 0)){

    ENZO_VFAIL("Error identifying species for Chemical Tracers = %"ISYM".\n",
                TestProblemData.MultiMetals)
  }


  return SUCCESS;
}
