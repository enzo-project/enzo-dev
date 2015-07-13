/* 
 * MakeFieldConservative
 *
 * Returns True if field should be multiplied by Density for various operations
 * that typically expect conservative fields.  This includes interpolation and flux correction.
 *
 * David Collins. 2011-05-25
 */
#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"


int MakeFieldConservative(field_type field){
    int MultiplyField = TRUE, counter = 0;
    //types_to_skip lists types not to be multiplied.
    //Please leave FieldUndefined as the last element, it terminates the loop.
    //
    field_type  types_to_skip[] = {Bfield1, Bfield2, Bfield3, PhiField,
        DrivingField1, DrivingField2, DrivingField3, GravPotential, DebugField,
    FieldUndefined};

    if( FieldTypeIsDensity(field) ) MultiplyField = FALSE;
    if( FieldTypeIsRadiation(field) ) MultiplyField = FALSE;
    if( FieldTypeNoInterpolate( field ) ) MultiplyField = FALSE;
    if( (field >= Velocity1 && field <= Velocity3) && HydroMethod == Zeus_Hydro ) {
        MultiplyField = FALSE;
    }
    while( types_to_skip[counter] != FieldUndefined ){
        if( field == types_to_skip[counter++] ){
            MultiplyField = FALSE;
            break;
        }
    }

    if( field == TotalEnergy && HydroMethod == MHD_Li && MHDCTUseSpecificEnergy == FALSE ){
        MultiplyField = FALSE;
    }

    return MultiplyField;
}
