/***********************************************************************
/
/ GRID CLASS (IDENTIFY THE CHEMICAL TRACER FIELDS DEPENDING ON MODEL)
/
/ written by: Andrew Emerick
/ date      : February, 2016
/ modified1 : April, 2016
/
/ PURPOSE:
/
/ NOTE: A few helper functions exist here to help with identifying
/       the chemical tracer species fields and to minimize the amount of
/       changes needed to define new fields. These functions can be
/       used to set the FieldType and DataLabel arrays at initialization
/       following examples in ChemicalEvolutionTestInitialize.C and
/       Grid_ChemicalEvolutionTestInitializeGrid.C
/
/       If one wants to add new species field, this requires only modifying
/       typedefs.h to include the new field, and adding in the proper
/       case (coded by atomic number) in the switch statements below.
/       This is a centralized way of handling many different labels for
/       both the baryon fields and paticles in a consistent manner.
/
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

int ChemicalSpeciesBaryonFieldNumber(const int &atomic_number){

  int field_num;

  /* Look up and return field index defined in typedefs.h */
  switch(atomic_number){
    case 1 :
    case 2 :
      if (TestProblemData.MultiSpecies == 0){
        ENZO_FAIL("ChemicalSpeciesBaryonFieldNumber: Multispecies must be ON to track H and He yields");
      }
      break;

    case  6 : field_num = CIDensity; break;
    case  7 : field_num = NDensity; break;
    case  8 : field_num = OIDensity; break;

    case 12 : field_num = MgDensity; break;

    case 14 : field_num = SiIDensity; break;

    case 26 : field_num = FeDensity; break;

    case 39 : field_num = YDensity; break;

    case 56 : field_num = BaDensity; break;
    case 57 : field_num = LaDensity; break;

    case 63 : field_num = EuDensity; break;

    default:
      ENZO_FAIL("ChemicalSpeciesBaryonFieldNumber: Cannot find field");

  }

  return field_num;

}

char* ChemicalSpeciesBaryonFieldLabel(const int &atomic_number){

  char *label = {};

  /* For a given atomic number, return the name we should assign to the baryon field
     DataLabel */
  switch(atomic_number){
    case 1 : // Handled via multispecies - this can be re-worked in the case when MultiSpecies
    case 2 : // is off, making passive H and He tracers, but why would one do this?
      if (TestProblemData.MultiSpecies == 0){
        ENZO_FAIL("ChemicalSpeciesBaryonFieldLabel: MultiSpecies must be ON to track yields for H or He");
      }
      break;

    case  6 : label = "C_Density"; break;
    case  7 : label = "N_Density"; break;
    case  8 : label = "O_Density"; break;

    case 12 : label = "Mg_Density"; break;

    case 14 : label = "Si_Density"; break;

    case 26 : label = "Fe_Density"; break;

    case 39 : label = "Y_Density"; break;

    case 56 : label = "Ba_Density"; break;
    case 57 : label = "La_Density"; break;

    case 63 : label = "Eu_Density"; break;

    default:
      ENZO_FAIL("Error in ChemicalSpeciesBaryonFieldLabel - Label not found\n");
  }

  return label;
}

char* ChemicalSpeciesParticleLabel(const int &atomic_number){

  char *label = {};

  /* For a given atomic number, return the name we should assign
     to the particle attribute for chemical tagging purposes */

  switch(atomic_number){
    case  1 : label = "H_Fraction" ; break;
    case  2 : label = "He_Fraction"; break;

    case  6 : label = "C_Fraction" ; break;
    case  7 : label = "N_Fraction" ; break;
    case  8 : label = "O_Fraction" ; break;

    case 12 : label = "Mg_Fraction"; break;

    case 14 : label = "Si_Fraction"; break;

    case 26 : label = "Fe_Fraction"; break;

    case 39 : label = "Y_Fraction" ; break;

    case 56 : label = "Ba_Fraction"; break;
    case 57 : label = "La_Fraction"; break;

    case 63 : label = "Eu_Fraction"; break;

    default:
      ENZO_FAIL("Error in ChemicalSpeciesParticleLabel - Label not found\n");
  }

  return label;
}


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

    /* for all other cases can do field_num = FindField(ChemicalSpeciesBaryonFieldNumber(atomic_number), FieldType, NumberO */

    case  6: field_num = FindField( CIDensity, FieldType, NumberOfBaryonFields); break;
    case  7: field_num = FindField( NDensity, FieldType, NumberOfBaryonFields); break;
    case  8: field_num = FindField( OIDensity, FieldType, NumberOfBaryonFields); break;

    case 12: field_num = FindField(MgDensity, FieldType, NumberOfBaryonFields); break;

    case 14: field_num = FindField(SiIDensity, FieldType, NumberOfBaryonFields); break;

    case 26: field_num = FindField(FeDensity, FieldType, NumberOfBaryonFields); break;

    case 39: field_num = FindField(YDensity, FieldType, NumberOfBaryonFields); break;

    case 56: field_num = FindField(BaDensity, FieldType, NumberOfBaryonFields); break;
    case 57: field_num = FindField(LaDensity, FieldType, NumberOfBaryonFields); break;

    case 63: field_num = FindField(EuDensity, FieldType, NumberOfBaryonFields); break;
  }

  if(field_num < 0){
    ENZO_FAIL("Error in IdentifyChemicalTracerSpeciesByNumber. No species found");
  }

  return SUCCESS;
}
