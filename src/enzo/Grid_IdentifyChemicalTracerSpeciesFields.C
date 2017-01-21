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

    case  3 : field_num = LiDensity; break;
    case  4 : field_num = BeDensity; break;
    case  5 : field_num = BDensity; break;
    case  6 : field_num = CDensity; break;
    case  7 : field_num = NDensity; break;
    case  8 : field_num = ODensity; break;
    case  9 : field_num = FDensity; break;
    case 10 : field_num = NeDensity; break;
    case 11 : field_num = NaDensity; break;
    case 12 : field_num = MgDensity; break;
    case 13 : field_num = AlDensity; break;
    case 14 : field_num = SiDensity; break;
    case 15 : field_num = PDensity; break;
    case 16 : field_num = SDensity; break;
    case 17 : field_num = ClDensity; break;
    case 18 : field_num = ArDensity; break;
    case 19 : field_num = KDensity; break;
    case 20 : field_num = CaDensity; break;
    case 21 : field_num = ScDensity; break;
    case 22 : field_num = TiDensity; break;
    case 23 : field_num = VDensity; break;
    case 24 : field_num = CrDensity; break;
    case 25 : field_num = MnDensity; break;
    case 26 : field_num = FeDensity; break;
    case 27 : field_num = CoDensity; break;
    case 28 : field_num = NiDensity; break;
    case 29 : field_num = CuDensity; break;
    case 30 : field_num = ZnDensity; break;
    case 31 : field_num = GaDensity; break;
    case 32 : field_num = GeDensity; break;
    case 33 : field_num = AsDensity; break;
    case 34 : field_num = SeDensity; break;
    case 35 : field_num = BrDensity; break;
    case 36 : field_num = KrDensity; break;
    case 37 : field_num = RbDensity; break;
    case 38 : field_num = SrDensity; break;
    case 39 : field_num = YDensity; break;
    case 40 : field_num = ZrDensity; break;
    case 41 : field_num = NbDensity; break;
    case 42 : field_num = MoDensity; break;
    case 43 : field_num = TcDensity; break;
    case 44 : field_num = RuDensity; break;
    case 45 : field_num = RhDensity; break;
    case 46 : field_num = PdDensity; break;
    case 47 : field_num = AgDensity; break;
    case 48 : field_num = CdDensity; break;
    case 49 : field_num = InDensity; break;
    case 50 : field_num = SnDensity; break;
    case 51 : field_num = SbDensity; break;
    case 52 : field_num = TeDensity; break;
    case 53 : field_num =  IDensity; break;
    case 54 : field_num = XeDensity; break;
    case 55 : field_num = CsDensity; break;
    case 56 : field_num = BaDensity; break;
    case 57 : field_num = LaDensity; break;
    case 58 : field_num = CeDensity; break;
    case 59 : field_num = PrDensity; break;
    case 60 : field_num = NdDensity; break;
    case 61 : field_num = PmDensity; break;
    case 62 : field_num = SmDensity; break;
    case 63 : field_num = EuDensity; break;
    case 64 : field_num = GdDensity; break;
    case 65 : field_num = TbDensity; break;
    case 66 : field_num = DyDensity; break;
    case 67 : field_num = HoDensity; break;
    case 68 : field_num = ErDensity; break;
    case 69 : field_num = TmDensity; break;
    case 70 : field_num = YbDensity; break;
    case 71 : field_num = LuDensity; break;
    case 72 : field_num = HfDensity; break;
    case 73 : field_num = TaDensity; break;
    case 74 : field_num =  WDensity; break;
    case 75 : field_num = ReDensity; break;
    case 76 : field_num = OsDensity; break;
    case 77 : field_num = IrDensity; break;
    case 78 : field_num = PtDensity; break;
    case 79 : field_num = AuDensity; break;
    case 80 : field_num = HgDensity; break;
    case 81 : field_num = TlDensity; break;
    case 82 : field_num = PbDensity; break;
    case 83 : field_num = BiDensity; break;

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

    case  3 : label = "Li_Density"; break;
    case  4 : label = "Be_Density"; break;
    case  5 : label =  "B_Density"; break;
    case  6 : label =  "C_Density"; break;
    case  7 : label =  "N_Density"; break;
    case  8 : label =  "O_Density"; break;
    case  9 : label =  "F_Density"; break;
    case 10 : label = "Ne_Density"; break;
    case 11 : label = "Na_Density"; break;
    case 12 : label = "Mg_Density"; break;
    case 13 : label = "Al_Density"; break;
    case 14 : label = "Si_Density"; break;
    case 15 : label =  "P_Density"; break;
    case 16 : label =  "S_Density"; break;
    case 17 : label = "Cl_Density"; break;
    case 18 : label = "Ar_Density"; break;
    case 19 : label =  "K_Density"; break;
    case 20 : label = "Ca_Density"; break;
    case 21 : label = "Sc_Density"; break;
    case 22 : label = "Ti_Density"; break;
    case 23 : label =  "V_Density"; break;
    case 24 : label = "Cr_Density"; break;
    case 25 : label = "Mn_Density"; break;
    case 26 : label = "Fe_Density"; break;
    case 27 : label = "Co_Density"; break;
    case 28 : label = "Ni_Density"; break;
    case 29 : label = "Cu_Density"; break;
    case 30 : label = "Zn_Density"; break;
    case 31 : label = "Ga_Density"; break;
    case 32 : label = "Ge_Density"; break;
    case 33 : label = "As_Density"; break;
    case 34 : label = "Se_Density"; break;
    case 35 : label = "Br_Density"; break;
    case 36 : label = "Kr_Density"; break;
    case 37 : label = "Rb_Density"; break;
    case 38 : label = "Sr_Density"; break;
    case 39 : label =  "Y_Density"; break;
    case 40 : label = "Zr_Density"; break;
    case 41 : label = "Nb_Density"; break;
    case 42 : label = "Mo_Density"; break;
    case 43 : label = "Tc_Density"; break;
    case 44 : label = "Ru_Density"; break;
    case 45 : label = "Rh_Density"; break;
    case 46 : label = "Pd_Density"; break;
    case 47 : label = "Ag_Density"; break;
    case 48 : label = "Cd_Density"; break;
    case 49 : label = "In_Density"; break;
    case 50 : label = "Sn_Density"; break;
    case 51 : label = "Sb_Density"; break;
    case 52 : label = "Te_Density"; break;
    case 53 : label =  "I_Density"; break;
    case 54 : label = "Xe_Density"; break;
    case 55 : label = "Cs_Density"; break;
    case 56 : label = "Ba_Density"; break;
    case 57 : label = "La_Density"; break;
    case 58 : label = "Ce_Density"; break;
    case 59 : label = "Pr_Density"; break;
    case 60 : label = "Nd_Density"; break;
    case 61 : label = "Pm_Density"; break;
    case 62 : label = "Sm_Density"; break;
    case 63 : label = "Eu_Density"; break;
    case 64 : label = "Gd_Density"; break;
    case 65 : label = "Tb_Density"; break;
    case 66 : label = "Dy_Density"; break;
    case 67 : label = "Ho_Density"; break;
    case 68 : label = "Er_Density"; break;
    case 69 : label = "Tm_Density"; break;
    case 70 : label = "Yb_Density"; break;
    case 71 : label = "Lu_Density"; break;
    case 72 : label = "Hf_Density"; break;
    case 73 : label = "Ta_Density"; break;
    case 74 : label =  "W_Density"; break;
    case 75 : label = "Re_Density"; break;
    case 76 : label = "Os_Density"; break;
    case 77 : label = "Ir_Density"; break;
    case 78 : label = "Pt_Density"; break;
    case 79 : label = "Au_Density"; break;
    case 80 : label = "Hg_Density"; break;
    case 81 : label = "Tl_Density"; break;
    case 82 : label = "Pb_Density"; break;
    case 83 : label = "Bi_Density"; break;

    default:
      ENZO_FAIL("Error in ChemicalSpeciesBaryonFieldLabel - Label not found\n");
  }

  return label;
}

char* IndividualStarTableIDLabel(const int &num){

  char *label = {};
  switch(num){
    case 0 : label = "se_table_M_pos"; break;
    case 1 : label = "se_table_Z_pos"; break;
    case 2 : label = "rad_table_T_pos"; break;
    case 3 : label = "rad_table_g_pos"; break;
    case 4 : label = "rad_table_Z_pos"; break;
    case 5 : label = "yield_table_M_pos"; break;
    case 6 : label = "yield_table_Z_pos"; break;
    default:
      ENZO_FAIL("Looking for too large of a number for table position labels");
  }

  return label;
}

char* ChemicalSpeciesParticleLabel(const int &atomic_number){

  char *label = {};

  /* For a given atomic number, return the name we should assign
     to the particle attribute for chemical tagging purposes */

  switch(atomic_number){
    case  1 : label = "particle_H_fraction" ; break;
    case  2 : label = "particle_He_fraction"; break;
    case  3 : label = "particle_Li_fraction"; break;
    case  4 : label = "particle_Be_fraction"; break;
    case  5 : label =  "particle_B_fraction"; break;
    case  6 : label =  "particle_C_fraction"; break;
    case  7 : label =  "particle_N_fraction"; break;
    case  8 : label =  "particle_O_fraction"; break;
    case  9 : label =  "particle_F_fraction"; break;
    case 10 : label = "particle_Ne_fraction"; break;
    case 11 : label = "particle_Na_fraction"; break;
    case 12 : label = "particle_Mg_fraction"; break;
    case 13 : label = "particle_Al_fraction"; break;
    case 14 : label = "particle_Si_fraction"; break;
    case 15 : label =  "particle_P_fraction"; break;
    case 16 : label =  "particle_S_fraction"; break;
    case 17 : label = "particle_Cl_fraction"; break;
    case 18 : label = "particle_Ar_fraction"; break;
    case 19 : label =  "particle_K_fraction"; break;
    case 20 : label = "particle_Ca_fraction"; break;
    case 21 : label = "particle_Sc_fraction"; break;
    case 22 : label = "particle_Ti_fraction"; break;
    case 23 : label =  "particle_V_fraction"; break;
    case 24 : label = "particle_Cr_fraction"; break;
    case 25 : label = "particle_Mn_fraction"; break;
    case 26 : label = "particle_Fe_fraction"; break;
    case 27 : label = "particle_Co_fraction"; break;
    case 28 : label = "particle_Ni_fraction"; break;
    case 29 : label = "particle_Cu_fraction"; break;
    case 30 : label = "particle_Zn_fraction"; break;
    case 31 : label = "particle_Ga_fraction"; break;
    case 32 : label = "particle_Ge_fraction"; break;
    case 33 : label = "particle_As_fraction"; break;
    case 34 : label = "particle_Se_fraction"; break;
    case 35 : label = "particle_Br_fraction"; break;
    case 36 : label = "particle_Kr_fraction"; break;
    case 37 : label = "particle_Rb_fraction"; break;
    case 38 : label = "particle_Sr_fraction"; break;
    case 39 : label =  "particle_Y_fraction"; break;
    case 40 : label = "particle_Zr_fraction"; break;
    case 41 : label = "particle_Nb_fraction"; break;
    case 42 : label = "particle_Mo_fraction"; break;
    case 43 : label = "particle_Tc_fraction"; break;
    case 44 : label = "particle_Ru_fraction"; break;
    case 45 : label = "particle_Rh_fraction"; break;
    case 46 : label = "particle_Pd_fraction"; break;
    case 47 : label = "particle_Ag_fraction"; break;
    case 48 : label = "particle_Cd_fraction"; break;
    case 49 : label = "particle_In_fraction"; break;
    case 50 : label = "particle_Sn_fraction"; break;
    case 51 : label = "particle_Sb_fraction"; break;
    case 52 : label = "particle_Te_fraction"; break;
    case 53 : label =  "particle_I_fraction"; break;
    case 54 : label = "particle_Xe_fraction"; break;
    case 55 : label = "particle_Cs_fraction"; break;
    case 56 : label = "particle_Ba_fraction"; break;
    case 57 : label = "particle_La_fraction"; break;
    case 58 : label = "particle_Ce_fraction"; break;
    case 59 : label = "particle_Pr_fraction"; break;
    case 60 : label = "particle_Nd_fraction"; break;
    case 61 : label = "particle_Pm_fraction"; break;
    case 62 : label = "particle_Sm_fraction"; break;
    case 63 : label = "particle_Eu_fraction"; break;
    case 64 : label = "particle_Gd_fraction"; break;
    case 65 : label = "particle_Tb_fraction"; break;
    case 66 : label = "particle_Dy_fraction"; break;
    case 67 : label = "particle_Ho_fraction"; break;
    case 68 : label = "particle_Er_fraction"; break;
    case 69 : label = "particle_Tm_fraction"; break;
    case 70 : label = "particle_Yb_fraction"; break;
    case 71 : label = "particle_Lu_fraction"; break;
    case 72 : label = "particle_Hf_fraction"; break;
    case 73 : label = "particle_Ta_fraction"; break;
    case 74 : label =  "particle_W_fraction"; break;
    case 75 : label = "particle_Re_fraction"; break;
    case 76 : label = "particle_Os_fraction"; break;
    case 77 : label = "particle_Ir_fraction"; break;
    case 78 : label = "particle_Pt_fraction"; break;
    case 79 : label = "particle_Au_fraction"; break;
    case 80 : label = "particle_Hg_fraction"; break;
    case 81 : label = "particle_Tl_fraction"; break;
    case 82 : label = "particle_Pb_fraction"; break;
    case 83 : label = "particle_Bi_fraction"; break;


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
    case  3 : field_num = FindField( LiDensity, FieldType, NumberOfBaryonFields); break;
    case  4 : field_num = FindField( BeDensity, FieldType, NumberOfBaryonFields); break;
    case  5 : field_num = FindField( BDensity, FieldType, NumberOfBaryonFields); break;
    case  6 : field_num = FindField( CDensity, FieldType, NumberOfBaryonFields); break;
    case  7 : field_num = FindField( NDensity, FieldType, NumberOfBaryonFields); break;
    case  8 : field_num = FindField( ODensity, FieldType, NumberOfBaryonFields); break;
    case  9 : field_num = FindField( FDensity, FieldType, NumberOfBaryonFields); break;
    case 10 : field_num = FindField( NeDensity, FieldType, NumberOfBaryonFields); break;
    case 11 : field_num = FindField( NaDensity, FieldType, NumberOfBaryonFields); break;
    case 12 : field_num = FindField( MgDensity, FieldType, NumberOfBaryonFields); break;
    case 13 : field_num = FindField( AlDensity, FieldType, NumberOfBaryonFields); break;
    case 14 : field_num = FindField( SiDensity, FieldType, NumberOfBaryonFields); break;
    case 15 : field_num = FindField( PDensity, FieldType, NumberOfBaryonFields); break;
    case 16 : field_num = FindField( SDensity, FieldType, NumberOfBaryonFields); break;
    case 17 : field_num = FindField( ClDensity, FieldType, NumberOfBaryonFields); break;
    case 18 : field_num = FindField( ArDensity, FieldType, NumberOfBaryonFields); break;
    case 19 : field_num = FindField( KDensity, FieldType, NumberOfBaryonFields); break;
    case 20 : field_num = FindField( CaDensity, FieldType, NumberOfBaryonFields); break;
    case 21 : field_num = FindField( ScDensity, FieldType, NumberOfBaryonFields); break;
    case 22 : field_num = FindField( TiDensity, FieldType, NumberOfBaryonFields); break;
    case 23 : field_num = FindField( VDensity, FieldType, NumberOfBaryonFields); break;
    case 24 : field_num = FindField( CrDensity, FieldType, NumberOfBaryonFields); break;
    case 25 : field_num = FindField( MnDensity, FieldType, NumberOfBaryonFields); break;
    case 26 : field_num = FindField( FeDensity, FieldType, NumberOfBaryonFields); break;
    case 27 : field_num = FindField( CoDensity, FieldType, NumberOfBaryonFields); break;
    case 28 : field_num = FindField( NiDensity, FieldType, NumberOfBaryonFields); break;
    case 29 : field_num = FindField( CuDensity, FieldType, NumberOfBaryonFields); break;
    case 30 : field_num = FindField( ZnDensity, FieldType, NumberOfBaryonFields); break;
    case 31 : field_num = FindField( GaDensity, FieldType, NumberOfBaryonFields); break;
    case 32 : field_num = FindField( GeDensity, FieldType, NumberOfBaryonFields); break;
    case 33 : field_num = FindField( AsDensity, FieldType, NumberOfBaryonFields); break;
    case 34 : field_num = FindField( SeDensity, FieldType, NumberOfBaryonFields); break;
    case 35 : field_num = FindField( BrDensity, FieldType, NumberOfBaryonFields); break;
    case 36 : field_num = FindField( KrDensity, FieldType, NumberOfBaryonFields); break;
    case 37 : field_num = FindField( RbDensity, FieldType, NumberOfBaryonFields); break;
    case 38 : field_num = FindField( SrDensity, FieldType, NumberOfBaryonFields); break;
    case 39 : field_num = FindField( YDensity, FieldType, NumberOfBaryonFields); break;
    case 40 : field_num = FindField( ZrDensity, FieldType, NumberOfBaryonFields); break;
    case 41 : field_num = FindField( NbDensity, FieldType, NumberOfBaryonFields); break;
    case 42 : field_num = FindField( MoDensity, FieldType, NumberOfBaryonFields); break;
    case 43 : field_num = FindField( TcDensity, FieldType, NumberOfBaryonFields); break;
    case 44 : field_num = FindField( RuDensity, FieldType, NumberOfBaryonFields); break;
    case 45 : field_num = FindField( RhDensity, FieldType, NumberOfBaryonFields); break;
    case 46 : field_num = FindField( PdDensity, FieldType, NumberOfBaryonFields); break;
    case 47 : field_num = FindField( AgDensity, FieldType, NumberOfBaryonFields); break;
    case 48 : field_num = FindField( CdDensity, FieldType, NumberOfBaryonFields); break;
    case 49 : field_num = FindField( InDensity, FieldType, NumberOfBaryonFields); break;
    case 50 : field_num = FindField( SnDensity, FieldType, NumberOfBaryonFields); break;
    case 51 : field_num = FindField( SbDensity, FieldType, NumberOfBaryonFields); break;
    case 52 : field_num = FindField( TeDensity, FieldType, NumberOfBaryonFields); break;
    case 53 : field_num = FindField( IDensity, FieldType, NumberOfBaryonFields); break;
    case 54 : field_num = FindField( XeDensity, FieldType, NumberOfBaryonFields); break;
    case 55 : field_num = FindField( CsDensity, FieldType, NumberOfBaryonFields); break;
    case 56 : field_num = FindField( BaDensity, FieldType, NumberOfBaryonFields); break;
    case 57 : field_num = FindField( LaDensity, FieldType, NumberOfBaryonFields); break;
    case 58 : field_num = FindField( CeDensity, FieldType, NumberOfBaryonFields); break;
    case 59 : field_num = FindField( PrDensity, FieldType, NumberOfBaryonFields); break;
    case 60 : field_num = FindField( NdDensity, FieldType, NumberOfBaryonFields); break;
    case 61 : field_num = FindField( PmDensity, FieldType, NumberOfBaryonFields); break;
    case 62 : field_num = FindField( SmDensity, FieldType, NumberOfBaryonFields); break;
    case 63 : field_num = FindField( EuDensity, FieldType, NumberOfBaryonFields); break;
    case 64 : field_num = FindField( GdDensity, FieldType, NumberOfBaryonFields); break;
    case 65 : field_num = FindField( TbDensity, FieldType, NumberOfBaryonFields); break;
    case 66 : field_num = FindField( DyDensity, FieldType, NumberOfBaryonFields); break;
    case 67 : field_num = FindField( HoDensity, FieldType, NumberOfBaryonFields); break;
    case 68 : field_num = FindField( ErDensity, FieldType, NumberOfBaryonFields); break;
    case 69 : field_num = FindField( TmDensity, FieldType, NumberOfBaryonFields); break;
    case 70 : field_num = FindField( YbDensity, FieldType, NumberOfBaryonFields); break;
    case 71 : field_num = FindField( LuDensity, FieldType, NumberOfBaryonFields); break;
    case 72 : field_num = FindField( HfDensity, FieldType, NumberOfBaryonFields); break;
    case 73 : field_num = FindField( TaDensity, FieldType, NumberOfBaryonFields); break;
    case 74 : field_num = FindField( WDensity, FieldType, NumberOfBaryonFields); break;
    case 75 : field_num = FindField( ReDensity, FieldType, NumberOfBaryonFields); break;
    case 76 : field_num = FindField( OsDensity, FieldType, NumberOfBaryonFields); break;
    case 77 : field_num = FindField( IrDensity, FieldType, NumberOfBaryonFields); break;
    case 78 : field_num = FindField( PtDensity, FieldType, NumberOfBaryonFields); break;
    case 79 : field_num = FindField( AuDensity, FieldType, NumberOfBaryonFields); break;
    case 80 : field_num = FindField( HgDensity, FieldType, NumberOfBaryonFields); break;
    case 81 : field_num = FindField( TlDensity, FieldType, NumberOfBaryonFields); break;
    case 82 : field_num = FindField( PbDensity, FieldType, NumberOfBaryonFields); break;
    case 83 : field_num = FindField( BiDensity, FieldType, NumberOfBaryonFields); break;

  }

  if(field_num < 0){
    ENZO_FAIL("Error in IdentifyChemicalTracerSpeciesByNumber. No species found");
  }

  return SUCCESS;
}
