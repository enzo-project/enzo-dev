/***********************************************************************
/
/  PROBLEM TYPE CLASS
/
/  written by: Matthew Turk, Oliver Hahn
/  date:       July, 2010
/
/  PURPOSE:
/
************************************************************************/

#ifdef NEW_PROBLEM_TYPES
#include <string>
#include <map>
#include <iostream>
#include <stdexcept>
#include <stdio.h>
#include <math.h>

#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"

#include "ProblemType.h"

/* This will return a map associating strings with plugin_creator structs,
   and the main plugin map will be declared statically.  */
std::map< std::string, EnzoProblemType_creator *>& 
get_problem_types()
{
    static std::map< std::string, EnzoProblemType_creator* > problem_type_map;
    return problem_type_map;
}

/* This takes a string, grabs the (static) plugin map defined above, and
   returns the plugin creator for that. */
EnzoProblemType *select_problem_type( std::string problem_type_name)
{
    EnzoProblemType_creator *ept_creator = get_problem_types()
            [ problem_type_name ];

    /* Simply throw an error if no such plugin exists... */
    if( !ept_creator )
    {   
        ENZO_FAIL("Unknown output plug-in.");
    }

    EnzoProblemType *ptype = ept_creator->create();

    return ptype;

}

#endif
