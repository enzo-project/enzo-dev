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
#ifndef __EVENT_HOOKS__
#define __EVENT_HOOKS__

/*
   Plugins maybe should instead be classes, so that they can store state and be
   instantiated by things like EnzoProblemType instances and be registered in
   that fashion.
*/
typedef void(*plugin_function)(HierarchyEntry *Grids[], TopGridData &MetaData);

std::multimap<std::string, std::string> &get_event_hooks();
std::map<std::string, plugin_function> &get_plugins();
void RunEventHooks(std::string event_name, HierarchyEntry *Grid[],
                    TopGridData &MetaData);
void RegisterEventPlugin(std::string plugin_name, plugin_function the_plugin);
void RegisterEventHook(std::string event_name, std::string plugin_name);

/*
   This should really be stored in a static variable returned by a function,
   but I was having a hard time getting that to be immutable and passed through
   pointers ...
*/

EXTERN std::multimap<std::string, std::string> event_hooks;
EXTERN std::map<std::string, plugin_function> plugins;
#endif
#endif
