/***********************************************************************
/
/  EVENT HOOK HANDLERS
/
/  written by: Matthew Turk
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
#include "EventHooks.h"

void RunEventHooks(std::string event_name, HierarchyEntry *Grids[],
                    TopGridData &MetaData)
{
    if (event_hooks.empty()) {
        return;
    }
    std::pair<std::multimap<std::string, std::string>::iterator,
              std::multimap<std::string, std::string>::iterator> it;

    it = event_hooks.equal_range(event_name);

    std::multimap<std::string, std::string>::iterator itr = it.first;

    if (itr == event_hooks.end())
    {
        /* This is commented out until further notice */
        /*std::cout << "Event plugin for hook " << event_name;
        std::cout << " not found." << std::endl;*/
        return;
    }
    for (; itr != it.second ; itr++)
    {
        std::string plugin_name = (*itr).second;
        /* Debugging statement is next */
        /*std::cout << "Loading event plugin for hook " << event_name;
        std::cout << " with name " << plugin_name << std::endl;*/
        plugin_function the_plugin = plugins[plugin_name];
        the_plugin(Grids, MetaData);
    }

}

void RegisterEventHook(std::string event_name, std::string plugin_name)
{
    std::cout << "Registering " << plugin_name << " for " << event_name << std::endl;

    event_hooks.insert(std::pair<std::string, std::string>
                    (event_name, plugin_name));
                    
    if (event_hooks.empty()) {
        std::cout << "EVENT HOOKS IS EMPTY1" << std::endl;
    }
    if (event_hooks.empty()) {
        std::cout << "EVENT HOOKS IS EMPTY2" << std::endl;
    }
}

void RegisterEventPlugin(std::string plugin_name, plugin_function the_plugin)
{
    std::cout << "Registering plugin " << plugin_name << std::endl;

    plugins.insert(std::pair<std::string, plugin_function >
                    (plugin_name, the_plugin));
}

#endif
