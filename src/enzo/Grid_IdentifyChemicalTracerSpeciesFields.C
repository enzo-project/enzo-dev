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

int grid::IdentifyChemicalTracerSpeciesFieldsByNumber(int &field_num, const int &atomic_number){

  this->IdentifyChemicalTracerSpeciesFieldsByNumber(field_num, atomic_number, 0);

}

int grid::IdentifyChemicalTracerSpeciesFieldsByNumber(int &field_num,
                                                      const int &atomic_number,
                                                      const int &ion_level){
 /* -------------------------------------------------------------------------------------
  * IdentifyChemicalTracerSpeciesByNumber
  * ------------------------------------------------------------------------------------
  * Author: A. Emerick - 04/21/16
  *
  * Lookup table function to grab the baryon field index number corresponding to
  * the desired chemical species tracer fields (inlcuding H and He). Elements are 
  * found by atomic number. If H or He are requested, last parameter is used (optional)
  * ---------------------------------------------------------------------------------- */

  field_num = -1;


  switch(atomic_number){

    /* Probably will never be used, but allow atomic number of 0 to refer to metallicity */
    case  0: field_num = FindField( Metallicity, FieldType, NumberOfBaryonFields); break;

    case  1: // Hydrogen a little complicated since we can track ionization
      switch( ion_level ){
        case 1 : field_num = FindField( HIDensity, FieldType, NumberOfBaryonFields); break;
        case 0 : // if ion_level is default value, assume ionized
        case 2 : field_num = FindField(HIIDensity, FieldType, NumberOfBaryonFields); break;
      }
      break; // end H case

    case  2: // Helium a little complicated since we can track ionization
      switch( ion_level ){
        case 1 : field_num = FindField(  HeIDensity, FieldType, NumberOfBaryonFields); break;
        case 2 : field_num = FindField( HeIIDensity, FieldType, NumberOfBaryonFields); break;
        case 0 : // if default, assume ionized
        case 3 : field_num = FindField(HeIIIDensity, FieldType, NumberOfBaryonFields); break;
      }
      break; // end He case

    /* only one case for all other species */

    case  6: field_num = FindField( CIDensity, FieldType, NumberOfBaryonFields); break;
    case  7: field_num = FindField( NIDensity, FieldType, NumberOfBaryonFields); break;
    case  8: field_num = FindField( OIDensity, FieldType, NumberOfBaryonFields); break;

    case 12: field_num = FindField(MgIDensity, FieldType, NumberOfBaryonFields); break;

    case 14: field_num = FindField(SiIDensity, FieldType, NumberOfBaryonFields); break;

    case 26: field_num = FindField(FeIDensity, FieldType, NumberOfBaryonFields); break;

    case 39: field_num = FindField(YIDensity, FieldType, NumberOfBaryonFields); break;

    case 56: field_num = FindField(BaIDensity, FieldType, NumberOfBaryonFields); break;
    case 57: field_num = FindField(LaIDensity, FieldType, NumberOfBaryonFields); break;

    case 63: field_num = FindField(EuIDensity, FieldType, NumberOfBaryonFields); break;
  }

  if(field_num < 0){
    ENZO_FAIL("Error in IdentifyChemicalTracerSpeciesByNumber. No species found");
  }

  return SUCCESS;
}

/* AJE overhaul 04 21 16
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
*/
